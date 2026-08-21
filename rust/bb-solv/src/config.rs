//! Solvation configuration and its resolution order.
//!
//! A setting is taken from the first source that supplies it: an explicit field set on
//! [`SolvConfig`], then the environment, then a `betterbuilder.toml` file, then the defaults in
//! this module. The defaults reproduce the values the AMSOL input files carry, so an absent config
//! file changes nothing.
//!
//! Environment variables: `BB_AMSOL_EXE`, `BB_AMSOL_TIMEOUT` (seconds),
//! `BB_AMSOL_MAX_CONCURRENT`, `BB_CONFIG` (path to the config file).

use std::path::{Path, PathBuf};

use serde::Deserialize;

/// Keyword text following `CHARGE=<q>` on the AMSOL keyword line.
pub const DEFAULT_METHOD: &str = "AM1 1SCF TLIMIT=15 GEO-OK SM5.42R";

/// Continuation lines (written with a leading `& `) selecting the water solvent.
pub const DEFAULT_WATER_SOLVENT: &[&str] = &["SOLVNT=WATER"];

/// Continuation lines (written with a leading `& `) selecting the hexadecane solvent.
pub const DEFAULT_HEXADECANE_SOLVENT: &[&str] = &[
    "SOLVNT=GENORG IOFR=1.4345 ALPHA=0.00 BETA=0.00 GAMMA=38.93",
    "DIELEC=2.06 FACARB=0.00 FEHALO=0.00 DEV",
];

/// Seconds allowed for one AMSOL run before it is killed. Bounds a hung process; a leg that runs
/// long is not one that has hung.
pub const DEFAULT_TIMEOUT_SECS: u64 = 600;

/// AMSOL processes allowed to run at once.
pub const DEFAULT_MAX_CONCURRENT: usize = 4;

/// The config file's `[amsol]` table. Every field is optional; an absent field falls through to
/// the next source in the resolution order.
#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AmsolFile {
    pub exe: Option<PathBuf>,
    /// Prepended to `LD_LIBRARY_PATH` for the AMSOL child process. The AMSOL7.1 binary is linked
    /// against the g77 runtime `libg2c.so.0`, which is not present on current systems and ships in
    /// the pipeline's `extralibs` directory.
    pub ld_library_path: Option<Vec<PathBuf>>,
    pub timeout_secs: Option<u64>,
    pub max_concurrent: Option<usize>,
    pub method: Option<String>,
    pub water_solvent: Option<Vec<String>>,
    pub hexadecane_solvent: Option<Vec<String>>,
}

/// A parsed config file.
#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ConfigFile {
    #[serde(default)]
    pub amsol: AmsolFile,
}

#[derive(Debug)]
pub enum ConfigError {
    /// The config file could not be read.
    Read {
        path: PathBuf,
        source: std::io::Error,
    },
    /// The config file is not valid TOML, or carries an unknown key.
    Parse { path: PathBuf, message: String },
    /// An environment variable holds a value that does not parse as the expected type.
    Env { var: &'static str, value: String },
    /// No AMSOL executable was supplied by any source.
    NoAmsolExe,
}

impl std::fmt::Display for ConfigError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ConfigError::Read { path, source } => {
                write!(f, "reading {}: {source}", path.display())
            }
            ConfigError::Parse { path, message } => {
                write!(f, "parsing {}: {message}", path.display())
            }
            ConfigError::Env { var, value } => {
                write!(f, "{var} is set to {value:?}, which does not parse")
            }
            ConfigError::NoAmsolExe => write!(
                f,
                "no AMSOL executable configured: pass one explicitly, set BB_AMSOL_EXE, \
                 or set amsol.exe in betterbuilder.toml"
            ),
        }
    }
}

impl std::error::Error for ConfigError {}

/// Resolved solvation settings.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SolvConfig {
    /// Path to the AMSOL executable. [`SolvConfig::resolve`] leaves this `None` when no source
    /// supplied one; [`SolvConfig::amsol_exe`] turns that into [`ConfigError::NoAmsolExe`].
    pub amsol_exe: Option<PathBuf>,
    /// Prepended to the AMSOL child's `LD_LIBRARY_PATH`; empty to leave it untouched.
    pub ld_library_path: Vec<PathBuf>,
    pub timeout_secs: u64,
    pub max_concurrent: usize,
    pub method: String,
    pub water_solvent: Vec<String>,
    pub hexadecane_solvent: Vec<String>,
}

impl Default for SolvConfig {
    fn default() -> Self {
        SolvConfig {
            amsol_exe: None,
            ld_library_path: Vec::new(),
            timeout_secs: DEFAULT_TIMEOUT_SECS,
            max_concurrent: DEFAULT_MAX_CONCURRENT,
            method: DEFAULT_METHOD.to_string(),
            water_solvent: DEFAULT_WATER_SOLVENT
                .iter()
                .map(|s| s.to_string())
                .collect(),
            hexadecane_solvent: DEFAULT_HEXADECANE_SOLVENT
                .iter()
                .map(|s| s.to_string())
                .collect(),
        }
    }
}

/// Read and parse a config file.
pub fn load_file(path: &Path) -> Result<ConfigFile, ConfigError> {
    let text = std::fs::read_to_string(path).map_err(|source| ConfigError::Read {
        path: path.to_path_buf(),
        source,
    })?;
    toml::from_str(&text).map_err(|e| ConfigError::Parse {
        path: path.to_path_buf(),
        message: e.to_string(),
    })
}

/// The config file to read: `BB_CONFIG` if set, otherwise `betterbuilder.toml` in the current
/// directory when it exists.
pub fn default_config_path() -> Option<PathBuf> {
    if let Some(p) = std::env::var_os("BB_CONFIG") {
        return Some(PathBuf::from(p));
    }
    let cwd = PathBuf::from("betterbuilder.toml");
    cwd.is_file().then_some(cwd)
}

impl SolvConfig {
    /// Resolve settings from `explicit`, the environment, `file`, and the defaults, in that order.
    ///
    /// Fields set on `explicit` win outright. `file` is the parsed config file, if there is one.
    pub fn resolve(explicit: &AmsolFile, file: Option<&ConfigFile>) -> Result<Self, ConfigError> {
        let f = file.map(|c| &c.amsol);
        let d = SolvConfig::default();

        let env_path = std::env::var_os("BB_AMSOL_EXE").map(PathBuf::from);
        let amsol_exe = explicit
            .exe
            .clone()
            .or(env_path)
            .or_else(|| f.and_then(|f| f.exe.clone()));

        let timeout_secs = match explicit.timeout_secs {
            Some(v) => v,
            None => match env_u64("BB_AMSOL_TIMEOUT")? {
                Some(v) => v,
                None => f.and_then(|f| f.timeout_secs).unwrap_or(d.timeout_secs),
            },
        };

        let max_concurrent = match explicit.max_concurrent {
            Some(v) => v,
            None => match env_usize("BB_AMSOL_MAX_CONCURRENT")? {
                Some(v) => v,
                None => f.and_then(|f| f.max_concurrent).unwrap_or(d.max_concurrent),
            },
        };

        Ok(SolvConfig {
            amsol_exe,
            ld_library_path: explicit
                .ld_library_path
                .clone()
                .or_else(|| env_path_list("BB_AMSOL_LD_LIBRARY_PATH"))
                .or_else(|| f.and_then(|f| f.ld_library_path.clone()))
                .unwrap_or(d.ld_library_path),
            timeout_secs,
            max_concurrent: max_concurrent.max(1),
            method: explicit
                .method
                .clone()
                .or_else(|| f.and_then(|f| f.method.clone()))
                .unwrap_or(d.method),
            water_solvent: explicit
                .water_solvent
                .clone()
                .or_else(|| f.and_then(|f| f.water_solvent.clone()))
                .unwrap_or(d.water_solvent),
            hexadecane_solvent: explicit
                .hexadecane_solvent
                .clone()
                .or_else(|| f.and_then(|f| f.hexadecane_solvent.clone()))
                .unwrap_or(d.hexadecane_solvent),
        })
    }

    /// Resolve using the config file named by [`default_config_path`].
    pub fn from_env(explicit: &AmsolFile) -> Result<Self, ConfigError> {
        let file = match default_config_path() {
            Some(p) => Some(load_file(&p)?),
            None => None,
        };
        SolvConfig::resolve(explicit, file.as_ref())
    }

    /// The configured AMSOL executable, or [`ConfigError::NoAmsolExe`].
    pub fn amsol_exe(&self) -> Result<&Path, ConfigError> {
        self.amsol_exe.as_deref().ok_or(ConfigError::NoAmsolExe)
    }
}

/// A colon-separated path list from the environment.
fn env_path_list(var: &str) -> Option<Vec<PathBuf>> {
    let v = std::env::var_os(var)?;
    let list: Vec<PathBuf> = std::env::split_paths(&v).collect();
    (!list.is_empty()).then_some(list)
}

fn env_u64(var: &'static str) -> Result<Option<u64>, ConfigError> {
    parse_env(var)
}
fn env_usize(var: &'static str) -> Result<Option<usize>, ConfigError> {
    parse_env(var)
}
fn parse_env<T: std::str::FromStr>(var: &'static str) -> Result<Option<T>, ConfigError> {
    match std::env::var(var) {
        Ok(s) => s
            .parse()
            .map(Some)
            .map_err(|_| ConfigError::Env { var, value: s }),
        Err(_) => Ok(None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_reproduce_the_amsol_input_values() {
        let c = SolvConfig::resolve(&AmsolFile::default(), None).expect("resolved");
        assert_eq!(c.method, "AM1 1SCF TLIMIT=15 GEO-OK SM5.42R");
        assert_eq!(c.water_solvent, vec!["SOLVNT=WATER"]);
        assert_eq!(
            c.hexadecane_solvent,
            vec![
                "SOLVNT=GENORG IOFR=1.4345 ALPHA=0.00 BETA=0.00 GAMMA=38.93",
                "DIELEC=2.06 FACARB=0.00 FEHALO=0.00 DEV",
            ]
        );
        assert_eq!(c.timeout_secs, 600);
    }

    #[test]
    fn file_overrides_defaults_and_explicit_overrides_file() {
        let file: ConfigFile = toml::from_str(
            r#"
            [amsol]
            exe = "/from/file/amsol7.1"
            timeout_secs = 90
            method = "AM1 1SCF"
            "#,
        )
        .expect("parsed");

        let c = SolvConfig::resolve(&AmsolFile::default(), Some(&file)).expect("resolved");
        assert_eq!(
            c.amsol_exe.as_deref(),
            Some(Path::new("/from/file/amsol7.1"))
        );
        assert_eq!(c.timeout_secs, 90);
        assert_eq!(c.method, "AM1 1SCF");
        // untouched by the file, so still the default
        assert_eq!(c.water_solvent, vec!["SOLVNT=WATER"]);

        let explicit = AmsolFile {
            exe: Some(PathBuf::from("/explicit/amsol7.1")),
            timeout_secs: Some(5),
            ..Default::default()
        };
        let c = SolvConfig::resolve(&explicit, Some(&file)).expect("resolved");
        assert_eq!(
            c.amsol_exe.as_deref(),
            Some(Path::new("/explicit/amsol7.1"))
        );
        assert_eq!(c.timeout_secs, 5);
        // the file still supplies what the explicit config left unset
        assert_eq!(c.method, "AM1 1SCF");
    }

    #[test]
    fn missing_exe_is_reported_not_defaulted() {
        let c = SolvConfig::resolve(&AmsolFile::default(), None).expect("resolved");
        assert!(matches!(c.amsol_exe(), Err(ConfigError::NoAmsolExe)));
    }

    #[test]
    fn unknown_keys_are_rejected() {
        let e = toml::from_str::<ConfigFile>("[amsol]\nexee = \"typo\"\n");
        assert!(e.is_err(), "an unknown key must not be silently ignored");
    }

    #[test]
    fn max_concurrent_is_never_zero() {
        let explicit = AmsolFile {
            max_concurrent: Some(0),
            ..Default::default()
        };
        let c = SolvConfig::resolve(&explicit, None).expect("resolved");
        assert_eq!(c.max_concurrent, 1);
    }
}
