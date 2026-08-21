//! Build the cxx bridge and link the patched RDKit (2026.09.1pre + amide torsions), resolved from
//! the environment:
//!   BB_RDKIT_ROOT   — patched RDKit tree (contains `Code/`, `lib/`, `External/…`). Falls back to
//!                     RDKit's standard `RDBASE`; panics if neither is set.
//!   BB_DEPS_PREFIX  — prefix providing boost + eigen headers/libs. Falls back to probing
//!                     /opt/homebrew, /usr/local, /usr; panics if none provides eigen3.
//!
//! Both resolved paths are re-exported as `cargo:root` / `cargo:deps` via the `links` key, reaching
//! dependents as `DEP_RDKITPATCHED_ROOT` / `DEP_RDKITPATCHED_DEPS`.

use std::env;
use std::path::Path;

fn rdkit_root() -> String {
    env::var("BB_RDKIT_ROOT")
        .or_else(|_| env::var("RDBASE"))
        .unwrap_or_else(|_| {
            panic!(
                "bb-rdkit: patched RDKit not found. Set BB_RDKIT_ROOT (or RDBASE) to the RDKit build \
                 root — the directory containing Code/ and lib/. RDKit is an external dependency and \
                 is intentionally not bundled or hardcoded."
            )
        })
}

/// A prefix shared by boost and eigen, used when neither has its own. Probes the standard system
/// locations for `<prefix>/include/eigen3`.
fn shared_deps_prefix() -> Option<String> {
    if let Ok(p) = env::var("BB_DEPS_PREFIX") {
        return Some(p);
    }
    ["/opt/homebrew", "/usr/local", "/usr"]
        .into_iter()
        .find(|c| Path::new(&format!("{c}/include/eigen3")).exists())
        .map(str::to_string)
}

/// Prefix for one dependency: its own variable, else the shared prefix.
fn dep_prefix(var: &str, what: &str) -> String {
    if let Ok(p) = env::var(var) {
        return p;
    }
    shared_deps_prefix().unwrap_or_else(|| {
        panic!(
            "bb-rdkit: could not locate {what}. Set {var} to its install prefix, or \
             BB_DEPS_PREFIX to a prefix providing both boost and eigen. \
             `. toolchain/env.sh` sets these from the toolchain built by toolchain/build-all.sh."
        )
    })
}

fn main() {
    let rdkit = rdkit_root();
    let boost = dep_prefix("BB_BOOST_PREFIX", "boost");
    let eigen = dep_prefix("BB_EIGEN_PREFIX", "eigen");

    cxx_build::bridge("src/lib.rs")
        .file("src/bridge.cc")
        .std("c++20") // RDKit 2026.09 headers require C++20
        .include("..") // resolves include!("bb-rdkit/src/bridge.h")
        .include(format!("{rdkit}/Code"))
        .include(format!(
            "{rdkit}/External/RingFamilies/RingDecomposerLib/src/RingDecomposerLib"
        ))
        .include(format!("{boost}/include"))
        .include(format!("{eigen}/include/eigen3"))
        .compile("bb_rdkit_bridge");

    // The C++ runtime RDKit and the bridge were built against. Linked and rpath'd explicitly so
    // the system libstdc++ is not used when it is older than the compiler that built RDKit.
    if let Ok(cxx) = env::var("BB_CXX_PREFIX") {
        for dir in ["lib64", "lib"] {
            let p = format!("{cxx}/{dir}");
            if Path::new(&p).is_dir() {
                println!("cargo:rustc-link-search=native={p}");
                println!("cargo:rustc-link-arg=-Wl,-rpath,{p}");
            }
        }
        println!("cargo:rustc-link-lib=dylib=stdc++");
    }

    // Link search: the patched RDKit dylibs + boost.
    println!("cargo:rustc-link-search=native={rdkit}/lib");
    println!("cargo:rustc-link-search=native={boost}/lib");
    // The RDKit components the setup extraction needs (SMILES → bounds/torsions/chiral).
    for lib in [
        "RDKitSmilesParse",
        "RDKitGraphMol",
        "RDKitRDGeneral",
        "RDKitRDGeometryLib",
        "RDKitDistGeomHelpers",
        "RDKitDistGeometry",
        "RDKitForceFieldHelpers",
        "RDKitForceField",
        "RDKitSubstructMatch",
    ] {
        println!("cargo:rustc-link-lib=dylib={lib}");
    }
    // rpath so the dylibs are found at runtime (applies to this crate's bin/tests).
    println!("cargo:rustc-link-arg=-Wl,-rpath,{rdkit}/lib");
    println!("cargo:rustc-link-arg=-Wl,-rpath,{boost}/lib");

    // Re-export resolved paths to dependents via `links = "rdkitpatched"` → DEP_RDKITPATCHED_ROOT/DEPS.
    println!("cargo:root={rdkit}");
    println!("cargo:deps={boost}");

    for var in [
        "BB_RDKIT_ROOT",
        "RDBASE",
        "BB_DEPS_PREFIX",
        "BB_BOOST_PREFIX",
        "BB_EIGEN_PREFIX",
        "BB_CXX_PREFIX",
    ] {
        println!("cargo:rerun-if-env-changed={var}");
    }
    println!("cargo:rerun-if-changed=src/lib.rs");
    println!("cargo:rerun-if-changed=src/bridge.cc");
    println!("cargo:rerun-if-changed=src/bridge.h");
}
