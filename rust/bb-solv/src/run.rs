//! Running AMSOL7.1.
//!
//! [`AmsolRunner::solvate`] writes the two input files into a caller-supplied directory, runs the
//! configured executable once per solvent with that directory as the child's working directory, and
//! captures each run's stdout. The two runs are independent and execute concurrently.
//!
//! Concurrency is bounded by a permit count, not by the calling thread pool: `max_concurrent`
//! AMSOL processes may be resident at once across all callers of one runner. Callers can therefore
//! drive molecules through their own parallelism without oversubscribing the machine, and a thread
//! blocked here waits on a permit rather than occupying a compute worker.

use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::{Condvar, Mutex};
use std::time::Duration;

use wait_timeout::ChildExt;

use crate::config::{ConfigError, SolvConfig};
use crate::input::{amsol_input, Solvent};
use crate::zmatrix::ZLine;

#[derive(Debug)]
pub enum RunError {
    Config(ConfigError),
    /// An input or output file could not be written, or the executable could not be spawned.
    Io {
        what: String,
        source: io::Error,
    },
    /// A run exceeded the configured timeout and was killed.
    Timeout {
        solvent: Solvent,
        secs: u64,
    },
    /// A run exited with a non-success status.
    Failed {
        solvent: Solvent,
        code: Option<i32>,
    },
}

impl std::fmt::Display for RunError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RunError::Config(e) => write!(f, "{e}"),
            RunError::Io { what, source } => write!(f, "{what}: {source}"),
            RunError::Timeout { solvent, secs } => {
                write!(f, "AMSOL {solvent:?} run exceeded {secs}s and was killed")
            }
            RunError::Failed { solvent, code } => match code {
                Some(c) => write!(f, "AMSOL {solvent:?} run exited with status {c}"),
                None => write!(f, "AMSOL {solvent:?} run terminated by signal"),
            },
        }
    }
}

impl std::error::Error for RunError {}

impl From<ConfigError> for RunError {
    fn from(e: ConfigError) -> Self {
        RunError::Config(e)
    }
}

/// Paths to the captured AMSOL output of one molecule.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SolvationRun {
    pub water: PathBuf,
    pub hexadecane: PathBuf,
}

/// A counting semaphore bounding resident AMSOL processes.
struct Permits {
    free: Mutex<usize>,
    available: Condvar,
}

impl Permits {
    fn new(n: usize) -> Self {
        Permits {
            free: Mutex::new(n),
            available: Condvar::new(),
        }
    }

    fn acquire(&self) {
        let mut free = self.free.lock().expect("permit mutex");
        while *free == 0 {
            free = self.available.wait(free).expect("permit mutex");
        }
        *free -= 1;
    }

    fn release(&self) {
        *self.free.lock().expect("permit mutex") += 1;
        self.available.notify_one();
    }
}

/// Runs AMSOL under a fixed configuration, bounding how many processes are resident at once.
pub struct AmsolRunner {
    cfg: SolvConfig,
    permits: Permits,
}

impl AmsolRunner {
    pub fn new(cfg: SolvConfig) -> Self {
        let n = cfg.max_concurrent.max(1);
        AmsolRunner {
            cfg,
            permits: Permits::new(n),
        }
    }

    pub fn config(&self) -> &SolvConfig {
        &self.cfg
    }

    /// Write the AMSOL input files for one molecule into `dir`.
    pub fn write_inputs(
        &self,
        dir: &Path,
        name: &str,
        charge: i32,
        zlines: &[ZLine],
    ) -> Result<(), RunError> {
        for solvent in [Solvent::Water, Solvent::Hexadecane] {
            let path = dir.join(solvent.input_filename());
            let text = amsol_input(&self.cfg, solvent, name, charge, zlines);
            std::fs::write(&path, text).map_err(|source| RunError::Io {
                what: format!("writing {}", path.display()),
                source,
            })?;
        }
        Ok(())
    }

    /// Run one solvent's calculation: `<exe> < temp.in-* > temp.o-*`, with `dir` as the child's
    /// working directory. Holds a permit for the duration.
    fn run_one(&self, dir: &Path, solvent: Solvent) -> Result<PathBuf, RunError> {
        let exe = self.cfg.amsol_exe()?.to_path_buf();
        let in_path = dir.join(solvent.input_filename());
        let out_path = dir.join(solvent.output_filename());

        let stdin = File::open(&in_path).map_err(|source| RunError::Io {
            what: format!("opening {}", in_path.display()),
            source,
        })?;
        let stdout = File::create(&out_path).map_err(|source| RunError::Io {
            what: format!("creating {}", out_path.display()),
            source,
        })?;

        let mut cmd = Command::new(&exe);
        cmd.current_dir(dir)
            .stdin(Stdio::from(stdin))
            .stdout(Stdio::from(stdout))
            .stderr(Stdio::null());
        if let Some(v) = self.ld_library_path() {
            cmd.env("LD_LIBRARY_PATH", v);
        }

        self.permits.acquire();
        let result = self.spawn_and_wait(&mut cmd, solvent, &exe);
        self.permits.release();
        result.map(|()| out_path)
    }

    fn spawn_and_wait(
        &self,
        cmd: &mut Command,
        solvent: Solvent,
        exe: &Path,
    ) -> Result<(), RunError> {
        let mut child = cmd.spawn().map_err(|source| RunError::Io {
            what: format!("spawning {}", exe.display()),
            source,
        })?;
        let timeout = Duration::from_secs(self.cfg.timeout_secs);
        match child.wait_timeout(timeout) {
            Ok(Some(status)) if status.success() => Ok(()),
            Ok(Some(status)) => Err(RunError::Failed {
                solvent,
                code: status.code(),
            }),
            Ok(None) => {
                let _ = child.kill();
                let _ = child.wait();
                Err(RunError::Timeout {
                    solvent,
                    secs: self.cfg.timeout_secs,
                })
            }
            Err(source) => Err(RunError::Io {
                what: format!("waiting on {}", exe.display()),
                source,
            }),
        }
    }

    /// `LD_LIBRARY_PATH` for the child: the configured entries ahead of any inherited value.
    fn ld_library_path(&self) -> Option<std::ffi::OsString> {
        if self.cfg.ld_library_path.is_empty() {
            return None;
        }
        let mut parts = self.cfg.ld_library_path.clone();
        if let Some(existing) = std::env::var_os("LD_LIBRARY_PATH") {
            parts.extend(std::env::split_paths(&existing));
        }
        std::env::join_paths(parts).ok()
    }

    /// Write the inputs for one molecule and run both solvents concurrently.
    ///
    /// `dir` receives `temp.in-wat`, `temp.in-hex`, `temp.o-wat` and `temp.o-hex`, and is the
    /// working directory of both child processes.
    pub fn solvate(
        &self,
        dir: &Path,
        name: &str,
        charge: i32,
        zlines: &[ZLine],
    ) -> Result<SolvationRun, RunError> {
        self.write_inputs(dir, name, charge, zlines)?;
        let (water, hexadecane) = std::thread::scope(|s| {
            let w = s.spawn(|| self.run_one(dir, Solvent::Water));
            let h = s.spawn(|| self.run_one(dir, Solvent::Hexadecane));
            (
                w.join().expect("water run thread"),
                h.join().expect("hexadecane run thread"),
            )
        });
        Ok(SolvationRun {
            water: water?,
            hexadecane: hexadecane?,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{AmsolFile, SolvConfig};
    use crate::zmatrix::to_zmatrix;
    use std::io::Write;
    use std::os::unix::fs::PermissionsExt;

    /// A stand-in for AMSOL: a shell script with the same stdin/stdout contract.
    fn stub(dir: &Path, body: &str) -> PathBuf {
        let p = dir.join("stub.sh");
        let mut f = File::create(&p).expect("create stub");
        writeln!(f, "#!/bin/sh\n{body}").expect("write stub");
        let mut perm = f.metadata().expect("metadata").permissions();
        perm.set_mode(0o755);
        std::fs::set_permissions(&p, perm).expect("chmod stub");
        p
    }

    fn cfg_with(exe: PathBuf, timeout_secs: u64) -> SolvConfig {
        SolvConfig::resolve(
            &AmsolFile {
                exe: Some(exe),
                timeout_secs: Some(timeout_secs),
                ..Default::default()
            },
            None,
        )
        .expect("resolved")
    }

    fn zlines() -> Vec<ZLine> {
        to_zmatrix(
            &[6, 1, 1],
            &[0.0, 0.0, 0.0, 1.09, 0.0, 0.0, -0.36, 1.03, 0.0],
            &[[0, 1], [0, 2]],
        )
        .expect("zmatrix")
    }

    #[test]
    fn captures_output_per_solvent() {
        let d = tempfile::tempdir().expect("tempdir");
        let exe = stub(d.path(), "cat");
        let r = AmsolRunner::new(cfg_with(exe, 30));
        let run = r.solvate(d.path(), "m", 0, &zlines()).expect("solvated");

        let wat = std::fs::read_to_string(&run.water).expect("read wat");
        let hex = std::fs::read_to_string(&run.hexadecane).expect("read hex");
        // the stub echoes its stdin, so each output is that solvent's input
        assert!(wat.contains("& SOLVNT=WATER"), "water output: {wat:.80}");
        assert!(hex.contains("& SOLVNT=GENORG"), "hex output: {hex:.80}");
        assert!(run.water.ends_with("temp.o-wat"));
        assert!(run.hexadecane.ends_with("temp.o-hex"));
    }

    #[test]
    fn timeout_kills_the_run() {
        let d = tempfile::tempdir().expect("tempdir");
        let exe = stub(d.path(), "sleep 30");
        let r = AmsolRunner::new(cfg_with(exe, 1));
        let err = r
            .solvate(d.path(), "m", 0, &zlines())
            .expect_err("must time out");
        assert!(
            matches!(err, RunError::Timeout { secs: 1, .. }),
            "expected timeout, got {err}"
        );
    }

    #[test]
    fn nonzero_exit_is_an_error() {
        let d = tempfile::tempdir().expect("tempdir");
        let exe = stub(d.path(), "exit 3");
        let r = AmsolRunner::new(cfg_with(exe, 30));
        let err = r
            .solvate(d.path(), "m", 0, &zlines())
            .expect_err("must fail");
        assert!(
            matches!(err, RunError::Failed { code: Some(3), .. }),
            "expected exit 3, got {err}"
        );
    }

    #[test]
    fn child_runs_in_the_supplied_directory_not_the_process_cwd() {
        let d = tempfile::tempdir().expect("tempdir");
        // the stub writes a relative path, which lands in the child's working directory
        let exe = stub(d.path(), "cat > /dev/null; echo here > witness.txt");
        let r = AmsolRunner::new(cfg_with(exe, 30));
        r.solvate(d.path(), "m", 0, &zlines()).expect("solvated");
        assert!(
            d.path().join("witness.txt").is_file(),
            "child cwd was not the run directory"
        );
        assert!(
            !Path::new("witness.txt").exists(),
            "child leaked into the parent's cwd"
        );
    }

    #[test]
    fn missing_executable_is_reported() {
        let d = tempfile::tempdir().expect("tempdir");
        let cfg = SolvConfig::resolve(&AmsolFile::default(), None).expect("resolved");
        let r = AmsolRunner::new(cfg);
        let err = r
            .solvate(d.path(), "m", 0, &zlines())
            .expect_err("must fail");
        assert!(
            matches!(err, RunError::Config(ConfigError::NoAmsolExe)),
            "expected NoAmsolExe, got {err}"
        );
    }

    #[test]
    fn permits_bound_resident_processes() {
        let d = tempfile::tempdir().expect("tempdir");
        // each run records its start, sleeps, then records its end; with one permit the two
        // solvent runs cannot overlap
        let exe = stub(
            d.path(),
            "cat > /dev/null; echo start >> log.txt; sleep 1; echo end >> log.txt",
        );
        let cfg = SolvConfig::resolve(
            &AmsolFile {
                exe: Some(exe),
                timeout_secs: Some(30),
                max_concurrent: Some(1),
                ..Default::default()
            },
            None,
        )
        .expect("resolved");
        AmsolRunner::new(cfg)
            .solvate(d.path(), "m", 0, &zlines())
            .expect("solvated");
        let log = std::fs::read_to_string(d.path().join("log.txt")).expect("read log");
        let seq: Vec<&str> = log.split_whitespace().collect();
        assert_eq!(
            seq,
            vec!["start", "end", "start", "end"],
            "runs overlapped despite max_concurrent = 1"
        );
    }
}
