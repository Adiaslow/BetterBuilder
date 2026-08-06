//! Add the RDKit dylib rpath to the Python extension module.
//!
//! `bb-rdkit`'s `rustc-link-lib`/`-search` directives propagate to this cdylib (so it links), but
//! `rustc-link-arg` (the runtime rpath) does not — so the final artifact that loads RDKit at import
//! time must add the rpath itself. The paths come from `bb-rdkit`'s build script via the `links`
//! mechanism (`DEP_RDKITPATCHED_ROOT` / `DEP_RDKITPATCHED_DEPS`), so there is no duplicated resolution
//! and nothing machine-specific in this file.

use std::env;

fn main() {
    // Set by bb-rdkit's build script (links = "rdkitpatched"); always present when bb-rdkit built.
    if let Ok(rdkit) = env::var("DEP_RDKITPATCHED_ROOT") {
        println!("cargo:rustc-link-arg=-Wl,-rpath,{rdkit}/lib");
    }
    if let Ok(deps) = env::var("DEP_RDKITPATCHED_DEPS") {
        println!("cargo:rustc-link-arg=-Wl,-rpath,{deps}/lib");
    }
}
