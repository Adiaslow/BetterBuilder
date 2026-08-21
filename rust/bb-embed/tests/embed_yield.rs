//! Embed-yield sanity gate: the reject/retry filters (init eigenvalue reject + per-atom energy
//! reject, wired for conformance target B) must not *collapse* conformer generation. This is not a
//! distributional gate (that is Phase 3 vs Divya's ensembles); it is a failable floor that catches a
//! filter/cap regression turning the embed into all-empty output.
//!
//! For a handful of real macrocycles it embeds several conformers each and requires every molecule to
//! yield at least one accepted (non-empty, n*3) conformer. If a filter over-rejects or the retry cap
//! is too tight, molecules go to zero yield and this fails. The actual per-molecule yield is printed.

use bb_spec::build_native;

const MOLS: &[&str] = &[
    "CC1CC=CCCC2(CCN(C(=O)c3ncc(F)cc3F)CC2)C(=O)NC[C@H]2CC[C@@H]2NC1=O",
    "COc1ccc2nc(C(=O)N3CC4(CCC=CCN5C(=O)C(=O)N(CC(=O)NCCC6CN(CCCO6)C4=O)C5=O)C3)cn2n1",
    "CC(=O)c1cnc2sc(C(=O)N3CC=CCN4C(=O)C(=O)N(CC(=O)NC[C@@H]5CCCN(C5)C(=O)C34)CC3CC3)cc2n1",
    "O=C1CCCCCCCCCCNC(=O)CCCCCCCCCN1",
    "CC1(C)CC=CCOCC2(CCCN(C(=O)CCOc3ccc(F)c(Cl)c3)C2)C(=O)NCC2CC2CNC1=O",
    "N#Cc1ccccc1",
];

const N_CONF: usize = 8;

#[test]
fn embed_does_not_collapse() {
    let mut total = 0usize;
    let mut nonempty = 0usize;
    let mut zero_yield: Vec<String> = Vec::new();

    for (i, smi) in MOLS.iter().enumerate() {
        let spec = match build_native(smi) {
            Ok(s) => s,
            Err(_) => continue,
        };
        let need = spec.n_atoms * 3;
        let confs = bb_embed::embed(&spec, N_CONF, 0x5EED_0000 ^ i as u64);
        let ok = confs.iter().filter(|c| c.coords.len() == need).count();
        total += confs.len();
        nonempty += ok;
        println!("  {ok}/{} non-empty (n={}) {}", confs.len(), spec.n_atoms, &smi[..smi.len().min(48)]);
        if ok == 0 {
            zero_yield.push((*smi).to_string());
        }
    }

    println!("aggregate embed yield: {nonempty}/{total} non-empty conformers");
    assert!(
        zero_yield.is_empty(),
        "embed collapsed to zero yield on {} molecule(s): {:?}",
        zero_yield.len(),
        zero_yield.iter().map(|s| &s[..s.len().min(40)]).collect::<Vec<_>>()
    );
}
