# Toolchain — the C++ dependencies `bb-rdkit` builds against

`bb-embed`, `bb-core` and `bb-solv` are pure Rust and need nothing here. Only `bb-rdkit` links C++,
and these scripts build what it links: a C++20 compiler, Boost, Eigen, and the patched RDKit.

Everything installs under `$BB_TOOLCHAIN_ROOT` (default `$HOME/toolchain`), one prefix per package.
Nothing is written outside it and no step needs root.

## Build

```sh
./build-all.sh              # gcc → eigen → boost → rdkit, in dependency order
. ../toolchain/env.sh       # point cargo at the result
cd ../rust && cargo test
```

Each step is idempotent and exits early if its prefix is already populated, so a re-run after a
failure resumes rather than restarting. Build trees go in `$BB_BUILD_SCRATCH` (default under `/tmp`)
and are not preserved. `BB_BUILD_JOBS` (default 32) sets parallelism — the default leaves headroom
on a shared machine rather than taking every core.

| File | Role |
|---|---|
| `versions.env` | Versions, the RDKit pin, and install prefixes. Every other script sources it. |
| `build-gcc.sh` | GCC, full 3-stage bootstrap |
| `build-eigen.sh` | Eigen headers |
| `build-boost.sh` | Boost, built with the GCC above |
| `build-rdkit.sh` | RDKit at the pin, with the diffs from `../rdkit-patch` |
| `build-all.sh` | All of the above in order |
| `env.sh` | Exports the prefixes `bb-rdkit/build.rs` reads |

## Why these versions

`versions.env` holds them; the constraints behind them:

- **Boost ≥ 1.81** is a hard floor — RDKit's `CMakeLists.txt` sets `RDK_BOOST_VERSION "1.81.0"`.
  Distribution packages older than that (Rocky 8 ships 1.66, EPEL 1.78) fail at configure time.
- **GCC ≥ 10** is required by RDKit's headers; GCC 8 and 9 fail on its type traits regardless of
  `-std=`.
- GCC 13.4.0 and Boost 1.82.0 match the environment the reference RDKit build used, so the
  toolchain is not a variable when comparing conformer output against the oracle.

## Notes that cost time to discover

- **Boost is found in CMake config mode**, so RDKit needs `CMAKE_PREFIX_PATH` pointing at the Boost
  prefix. `BOOST_ROOT` alone drives the older module mode and is ignored.
- **The option is `RDK_BUILD_PYTHON_WRAPPERS`**, plural. CMake ignores unknown `-D` variables
  silently, so the singular spelling leaves the Python wrapper enabled and pulls in a
  `boost_python` component that is not built here.
- **`CMAKE_INSTALL_RPATH` must include RDKit's own `lib`.** Linkers emit `RUNPATH` by default, which
  is not consulted for transitive dependencies, so without it the RDKit libraries cannot find each
  other at load time even when the consumer's rpath lists them.
- **The C++ runtime must come from this GCC.** Binaries built here reference symbols
  (`std::__throw_bad_array_new_length`, GCC 11+) absent from older system `libstdc++`. `env.sh` sets
  `BB_CXX_PREFIX` and `build.rs` links and rpaths it.

## Rebuilding RDKit only

```sh
./build-rdkit.sh    # re-applies the patch over a clean checkout
```
