"""Single source for the oracle container + AMSOL toolchain locations (validation/benchmark only).

The container path lives here under ONE env-var name (`BB_ORACLE_CONTAINER`); the AMSOL toolchain under
`BB_AMSOL_EXE` / `BB_AMSOL_LD_LIBRARY_PATH`. Every Python gate imports these instead of re-hardcoding the
`amurray2` absolute paths. Shell gates (`rust/gates.sh`, `validation/parity/dump.sh`) use the same
`BB_ORACLE_CONTAINER` name so there is exactly one knob.
"""
import os

APPTAINER = os.environ.get("APPTAINER", "/usr/bin/apptainer")
IMG = os.environ.get("BB_ORACLE_CONTAINER",
                     "/nfs/home/amurray2/toolchain/images/build_macrocycle_final.sif")
AMSOL_EXE = os.environ.get("BB_AMSOL_EXE", "/nfs/home/amurray2/toolchain/amsol/amsol7.1")
AMSOL_LIB = os.environ.get("BB_AMSOL_LD_LIBRARY_PATH", "/nfs/home/amurray2/toolchain/amsol/lib")


def amsol_env(base=None):
    """`base` env dict (default os.environ) overlaid with the AMSOL exe/lib the RDKit-free pipeline needs."""
    return dict(base if base is not None else os.environ,
                BB_AMSOL_EXE=AMSOL_EXE, BB_AMSOL_LD_LIBRARY_PATH=AMSOL_LIB)


def apptainer_exec(cmd, binds=()):
    """`[APPTAINER, exec, --bind b..., IMG, *cmd]` — the one shape for running inside the oracle container."""
    args = [APPTAINER, "exec"]
    for b in binds:
        args += ["--bind", b]
    return args + [IMG] + list(cmd)
