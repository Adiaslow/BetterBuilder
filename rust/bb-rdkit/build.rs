//! Build the cxx bridge and link the from-source patched RDKit (2026.09.1pre + amide torsions).
//!
//! No paths are baked into the repo — the patched RDKit is an external dependency the builder
//! provides, so this resolves it from the environment and fails with a clear message if absent:
//!   BB_RDKIT_ROOT   — patched RDKit tree (contains `Code/`, `lib/`, `External/…`). Falls back to
//!                     RDKit's standard `RDBASE`.
//!   BB_DEPS_PREFIX  — prefix providing boost + eigen headers/libs. Falls back to probing the
//!                     standard system prefixes (/opt/homebrew, /usr/local, /usr).
//!
//! The resolved paths are re-exported (`cargo:root` / `cargo:deps` via the `links` key) so the
//! Python-extension crate can add the runtime rpath without duplicating this logic.

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

fn deps_prefix() -> String {
    if let Ok(p) = env::var("BB_DEPS_PREFIX") {
        return p;
    }
    for cand in ["/opt/homebrew", "/usr/local", "/usr"] {
        if Path::new(&format!("{cand}/include/eigen3")).exists() {
            return cand.to_string();
        }
    }
    panic!(
        "bb-rdkit: could not locate boost + eigen. Set BB_DEPS_PREFIX to a prefix providing \
         <prefix>/include/eigen3 and boost headers."
    )
}

fn main() {
    let rdkit = rdkit_root();
    let deps = deps_prefix();

    cxx_build::bridge("src/lib.rs")
        .file("src/bridge.cc")
        .std("c++20") // RDKit 2026.09 requires C++20 (constexpr virtual in Geometry/point.h)
        .include("..") // resolves include!("bb-rdkit/src/bridge.h")
        .include(format!("{rdkit}/Code"))
        .include(format!(
            "{rdkit}/External/RingFamilies/RingDecomposerLib/src/RingDecomposerLib"
        ))
        .include(format!("{deps}/include")) // boost
        .include(format!("{deps}/include/eigen3"))
        .compile("bb_rdkit_bridge");

    // Link search: the patched RDKit dylibs + the dep prefix (boost).
    println!("cargo:rustc-link-search=native={rdkit}/lib");
    println!("cargo:rustc-link-search=native={deps}/lib");
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
    println!("cargo:rustc-link-arg=-Wl,-rpath,{deps}/lib");

    // Re-export resolved paths to dependents via `links = "rdkitpatched"` → DEP_RDKITPATCHED_ROOT/DEPS.
    println!("cargo:root={rdkit}");
    println!("cargo:deps={deps}");

    println!("cargo:rerun-if-env-changed=BB_RDKIT_ROOT");
    println!("cargo:rerun-if-env-changed=RDBASE");
    println!("cargo:rerun-if-env-changed=BB_DEPS_PREFIX");
    println!("cargo:rerun-if-changed=src/lib.rs");
    println!("cargo:rerun-if-changed=src/bridge.cc");
    println!("cargo:rerun-if-changed=src/bridge.h");
}
