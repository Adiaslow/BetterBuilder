//! Assemble a [`TypedMol`] from BetterBuilder's own pipeline — the wiring that makes the whole chain
//! run RDKit-free end to end: perception ([`bb_perceive`]) → SYBYL atom/bond typing + names
//! ([`crate::sybyl`]) → conformers ([`bb_embed`], via [`bb_spec`]). The result feeds both output
//! serializers ([`crate::write_mol2`] and [`crate::build`] → db2) with no derivation of one from the
//! other. Charges come separately from solvation ([`bb_solv::solv::SolvFile`], from `bb-solv`).

use crate::mol::TypedMol;
use crate::sybyl;
use bb_perceive::smarts_match::Perceived;
use bb_solv::solv::SolvFile;

/// Reshape an embed result into `confs[c][i] = [x,y,z]` (per conformer, per atom) — the form the
/// serializers consume. The single source for the `Conformer`→`Vec<Vec<[f64;3]>>` extraction.
pub fn confs_from_embed(raw: &[bb_embed::Conformer], n_atoms: usize) -> Vec<Vec<[f64; 3]>> {
    raw.iter().map(|c| (0..n_atoms).map(|i| c.atom(i)).collect()).collect()
}

/// Build a typed molecule from a parsed perception and a set of conformations. `name` is the molecule
/// identity; `confs[c][i]` is atom `i`'s coordinate in conformation `c`, in perception (atom) order.
pub fn assemble(name: &str, smiles: &str, p: &Perceived, confs: Vec<Vec<[f64; 3]>>) -> TypedMol {
    let n = p.atomic_numbers.len();
    let (types, amide) = sybyl::assign_sybyl_full(p);
    let names = sybyl::atom_names(p);
    let btypes = sybyl::bond_types(p, &types, &amide);

    // per-conformer UCSF torsion strain (matched once, evaluated per conformation) — the db2 S-records
    let hits = bb_strain::match_library(p);
    let (total_strain, max_strain): (Vec<f64>, Vec<f64>) =
        confs.iter().map(|c| bb_strain::strain(&hits, c)).unzip();

    let nconf = confs.len();
    let mut mol = TypedMol {
        name: name.to_string(),
        prot_name: "fake".into(),
        smiles: smiles.to_string(),
        longname: name.to_string(),
        atom_num: (1..=n).collect(),
        atom_name: names,
        atom_type: types,
        atom_bonds: Vec::new(),
        bond_num: (1..=p.bonds.len()).collect(),
        bond_start: p.bonds.iter().map(|&(a, _)| a + 1).collect(),
        bond_end: p.bonds.iter().map(|&(_, b)| b + 1).collect(),
        bond_type: btypes,
        atom_xyz: confs,
        input_total_strain: total_strain,
        input_max_strain: max_strain,
        input_hydrogens: vec![0; nconf],
    };
    mol.build_adjacency();
    mol
}

/// SMILES → typed molecule with BetterBuilder's own conformers: perceive, embed via the ETKDG
/// core-pin recipe ([`bb_spec::build_native`] + [`bb_embed::embed_recipe`]), and assemble. `seed`
/// drives the deterministic embedding.
pub fn from_smiles(name: &str, smiles: &str, seed: u64) -> Result<TypedMol, String> {
    // One perception, shared: the spec builder returns the `Perceived` it built from, so typing and
    // strain reuse it rather than re-perceiving the same SMILES.
    let (spec, p) = bb_spec::build_native_perceived(smiles)?;
    let confs = confs_from_embed(&bb_embed::embed_recipe(&spec, seed), spec.n_atoms);
    Ok(assemble(name, smiles, &p, confs))
}

/// The per-molecule pipeline, RDKit-free: SMILES → one [`Ligand`] (its mol2 + db2). Perceive, embed our
/// own conformers, solvate one of them for AMSOL charges ([`crate::solvation`], `workdir` holds the
/// AMSOL scratch), assemble the [`TypedMol`] (typing + strain), and serialize both members from it —
/// with no derivation between them. `name`'s base becomes the member prefix; `seed` drives embedding.
pub fn ligand_from_smiles(name: &str, prot_id: i32, smiles: &str, seed: u64, workdir: &std::path::Path) -> Result<crate::Ligand, String> {
    // One perception, shared across embedding, typing and strain (see `from_smiles`).
    let (spec, p) = bb_spec::build_native_perceived(smiles)?;
    let raw = bb_embed::embed_recipe(&spec, seed);
    if raw.is_empty() {
        return Err("no conformers embedded".into());
    }
    let confs = confs_from_embed(&raw, spec.n_atoms);

    // AMSOL charges from the first conformer (per-atom, shared across the ensemble)
    let solv = crate::solvation::solvate(&spec, &raw[0].coords, workdir)?;
    if solv.atoms.len() != spec.n_atoms {
        return Err(format!("solv atoms {} != {}", solv.atoms.len(), spec.n_atoms));
    }

    // The molecule identity carries the protomer id (`name.prot_id`), as Divya's driver names it; the
    // tarball member is `name.prot_id.C`, so `Ligand` keeps the bare name + `prot_id` separately.
    let full = format!("{name}.{prot_id}");
    let mol = assemble(&full, smiles, &p, confs);
    let db2 = build_blocks_db2(&mol, &solv, spec.sidechain_confs as usize)?;
    let charges: Vec<f64> = solv.atoms.iter().map(|a| a.charge).collect();
    let mol2 = crate::write_mol2(&mol, &charges, 0);
    Ok(crate::Ligand { name: name.to_string(), prot_id, formal_charge: solv.charge.round() as i32, mol2, db2 })
}

// File-render precisions — the number of decimals in the SIBLING writer format strings that the db2
// build must read *through* (her mol2db2 parses the rendered mol2/.solv, so native has to build from
// the same rounded values). Each const is tied to one format string; change the string, change the
// const here. (docs/2026-08-21-code-audit.md finding 10 — removes the bare 100.0/10000.0/10.0/1000.0
// factors that re-encoded these by hand.)
pub const COORD_DP: i32 = 4; //          mol2/SDF coords                     (mol2_out.rs %.4f)
pub const SOLV_CHARGE_DP: i32 = 4; //    .solv per-atom charge               (solv.rs %8.4f)
pub const SOLV_ATOM_DP: i32 = 2; //      .solv per-atom pol/area/apol/atomic (solv.rs %8.2f)
pub const SOLV_TOTAL_CHARGE_DP: i32 = 1; // .solv header total charge        (solv.rs %4.1f)
pub const SOLV_TOTAL_DP: i32 = 2; //     .solv header area / totals          (solv.rs %8.2f)
pub const STRAIN_DP: i32 = 3; //         db2 S-record strain                 (write.rs %+11.3f)

/// Round `v` to `dp` decimal places (`10^dp` is exact for these small `dp`, so this is bit-identical
/// to the previous inline `(v*100.0).round()/100.0` etc.).
#[inline]
pub fn round_dp(v: f64, dp: i32) -> f64 {
    let f = 10f64.powi(dp);
    (v * f).round() / f
}

/// Round a solvation file in place to the precisions its `.solv` render uses, so a db2 built from it is
/// byte-identical to one her `mol2db2` builds by parsing the rendered `.solv`.
pub fn round_solv_to_file_precision(solv: &mut SolvFile) {
    for a in solv.atoms.iter_mut() {
        a.charge = round_dp(a.charge, SOLV_CHARGE_DP);
        a.diff_pol = round_dp(a.diff_pol, SOLV_ATOM_DP);
        a.area = round_dp(a.area, SOLV_ATOM_DP);
        a.diff_apol = round_dp(a.diff_apol, SOLV_ATOM_DP);
        a.diff_atomic_solv = round_dp(a.diff_atomic_solv, SOLV_ATOM_DP);
    }
    // header totals — the db2 M2 line reads these from the file (else totals differ ~0.001 at M2 %+10.3f)
    solv.charge = round_dp(solv.charge, SOLV_TOTAL_CHARGE_DP);
    solv.area = round_dp(solv.area, SOLV_TOTAL_DP);
    solv.tot_diff_pol = round_dp(solv.tot_diff_pol, SOLV_TOTAL_DP);
    solv.tot_diff_apol = round_dp(solv.tot_diff_apol, SOLV_TOTAL_DP);
    solv.tot_diff_pol_plus_apol = round_dp(solv.tot_diff_pol_plus_apol, SOLV_TOTAL_DP);
}

/// Round every conformer's coords in place to the mol2/SDF file precision.
pub fn round_coords_to_file_precision(mol: &mut TypedMol) {
    for conf in mol.atom_xyz.iter_mut() {
        for a in conf.iter_mut() {
            *a = [round_dp(a[0], COORD_DP), round_dp(a[1], COORD_DP), round_dp(a[2], COORD_DP)];
        }
    }
}

/// Iterate `mol`'s core-seed blocks (`build_ligands.py` lines 534-620: her `mol2db2` runs once per
/// block of `sidechain_confs` conformers — they share that block's frozen pin_atoms, a rigid
/// component). `embed_recipe` emits confs block-ordered, so chunking `atom_xyz` by `sidechain_confs`
/// recovers her blocks. `f(block_index, &TypedMol)` sees a molecule whose coords/strain are that
/// block's window; the topology is cloned ONCE (not per block — finding 5), and strain is rounded to
/// the S-record precision (idempotent under the `%.3f` S-record, so this doesn't change the db2 while
/// letting the mol2 Strain comment's `%.6f` round-trip land on the same value — finding 2). `mol`
/// should already be coord-rounded ([`round_coords_to_file_precision`]).
pub fn for_each_block<E>(
    mol: &TypedMol,
    sidechain_confs: usize,
    mut f: impl FnMut(usize, &TypedMol) -> Result<(), E>,
) -> Result<(), E> {
    let nside = sidechain_confs.max(1);
    let nconf = mol.atom_xyz.len();
    let has_strain = mol.input_total_strain.len() == nconf;
    let has_hydrogens = mol.input_hydrogens.len() == nconf;
    let mut bmol = mol.clone(); // topology cloned once; per-conf fields sliced to the block below
    let (mut lo, mut b) = (0usize, 0usize);
    while lo < nconf {
        let hi = (lo + nside).min(nconf);
        bmol.atom_xyz = mol.atom_xyz[lo..hi].to_vec();
        // Every per-conformer field is windowed to THIS block so `build`'s block-local set index reads
        // the block's own values, not block 0's. (input_hydrogens is all-zero today, but slicing it
        // keeps that from being a load-bearing invariant.)
        if has_strain {
            bmol.input_total_strain =
                mol.input_total_strain[lo..hi].iter().map(|&v| round_dp(v, STRAIN_DP)).collect();
            bmol.input_max_strain =
                mol.input_max_strain[lo..hi].iter().map(|&v| round_dp(v, STRAIN_DP)).collect();
        }
        if has_hydrogens {
            bmol.input_hydrogens = mol.input_hydrogens[lo..hi].to_vec();
        }
        f(b, &bmol)?;
        lo = hi;
        b += 1;
    }
    Ok(())
}

/// Build the molecule's db2 the way `build_ligands.py` does: per core-seed block, from FILE-precision
/// inputs, concatenating the per-block db2 text. See [`for_each_block`] / [`round_solv_to_file_precision`].
pub fn build_blocks_db2(mol: &TypedMol, solv: &SolvFile, sidechain_confs: usize) -> Result<String, String> {
    let mut solv = solv.clone();
    round_solv_to_file_precision(&mut solv);
    let mut mol = mol.clone();
    round_coords_to_file_precision(&mut mol);
    let mut out = String::new();
    for_each_block(&mol, sidechain_confs, |_b, bmol| {
        out.push_str(&bb_db2::write_entry(&crate::build(bmol, &solv)?));
        Ok::<(), String>(())
    })?;
    Ok(out)
}

/// The protomer id for the next occurrence of `name`, matching `build_ligands.py`: the first line
/// for a name is `.0`, a second (another protonation state of the same molecule) `.1`, and so on.
/// `occ` carries the running per-name counts across a batch.
pub fn next_prot_id(occ: &mut std::collections::HashMap<String, i32>, name: &str) -> i32 {
    let c = occ.entry(name.to_string()).or_insert(-1);
    *c += 1;
    *c
}

/// The whole pipeline, RDKit-free: SMILES → the DOCK output tarball bytes (one ligand).
pub fn build_tarball_from_smiles(name: &str, smiles: &str, seed: u64, workdir: &std::path::Path) -> Result<Vec<u8>, String> {
    let lig = ligand_from_smiles(name, 0, smiles, seed, workdir)?;
    crate::write_tarball(std::slice::from_ref(&lig)).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::next_prot_id;
    use std::collections::HashMap;

    #[test]
    fn prot_ids_count_by_occurrence() {
        let mut occ = HashMap::new();
        // repeated names get 0,1,2…; distinct names each restart at 0 — no member-name collisions.
        assert_eq!(next_prot_id(&mut occ, "mc0001"), 0);
        assert_eq!(next_prot_id(&mut occ, "mc0001"), 1);
        assert_eq!(next_prot_id(&mut occ, "mc0002"), 0);
        assert_eq!(next_prot_id(&mut occ, "mc0001"), 2);
        assert_eq!(next_prot_id(&mut occ, "mc0002"), 1);
    }
}
