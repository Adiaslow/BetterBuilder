//! Assemble a [`TypedMol`] from BetterBuilder's own pipeline — the wiring that makes the whole chain
//! run RDKit-free end to end: perception ([`bb_perceive`]) → SYBYL atom/bond typing + names
//! ([`crate::sybyl`]) → conformers ([`bb_embed`], via [`bb_spec`]). The result feeds both output
//! serializers ([`crate::write_mol2`] and [`crate::build`] → db2) with no derivation of one from the
//! other. Charges come separately from solvation ([`bb_solv::solv::SolvFile`], from `bb-solv`).

use crate::mol::TypedMol;
use crate::name::{MoleculeName, NameError};
use crate::sybyl;
use crate::tarball::ChargeLetter;
use bb_perceive::smarts_match::Perceived;
use bb_solv::solv::SolvFile;

/// Reshape an embedded recipe into `confs[c][i] = [x,y,z]` (per conformer, per atom), block after
/// block — the form the serializers consume. The single source for the `Block`→`Vec<Vec<[f64;3]>>`
/// extraction; [`block_lens`] gives the partition of these conformers into blocks.
pub fn confs_from_blocks(blocks: &[bb_embed::Block], n_atoms: usize) -> Vec<Vec<[f64; 3]>> {
    blocks
        .iter()
        .flat_map(|b| &b.conformers)
        .map(|c| (0..n_atoms).map(|i| c.atom(i)).collect())
        .collect()
}

/// The number of conformers in each block, in block order: the partition of [`confs_from_blocks`]'s
/// conformers that [`for_each_block`] and [`build_blocks_db2`] build one db2 hierarchy per part of.
pub fn block_lens(blocks: &[bb_embed::Block]) -> Vec<usize> {
    blocks.iter().map(|b| b.conformers.len()).collect()
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
    let blocks = bb_embed::embed_recipe(&spec, seed).map_err(|e| e.to_string())?;
    Ok(assemble(name, smiles, &p, confs_from_blocks(&blocks, spec.n_atoms)))
}

/// The per-molecule pipeline, RDKit-free: SMILES → one [`crate::Ligand`] (its mol2 + db2). Perceive, embed our
/// own conformers, solvate one of them for AMSOL charges ([`crate::solvation`], `workdir` holds the
/// AMSOL scratch), assemble the [`TypedMol`] (typing + strain), and serialize both members from it —
/// with no derivation between them. `name` becomes the member prefix; `seed` drives embedding. A net
/// formal charge with no member-name letter ([`crate::tarball::ChargeLetter`]) is an error for this
/// molecule, so the archive never has to name it.
pub fn ligand_from_smiles(name: &MoleculeName, prot_id: i32, smiles: &str, seed: u64, workdir: &std::path::Path) -> Result<crate::Ligand, String> {
    // One perception, shared across embedding, typing and strain (see `from_smiles`).
    let (spec, p) = bb_spec::build_native_perceived(smiles)?;
    let blocks = bb_embed::embed_recipe(&spec, seed).map_err(|e| e.to_string())?;
    let confs = confs_from_blocks(&blocks, spec.n_atoms);

    // AMSOL charges from the first conformer (per-atom, shared across the ensemble)
    let solv = crate::solvation::solvate(&spec, &blocks[0].conformers[0].coords, workdir)?;
    if solv.atoms.len() != spec.n_atoms {
        return Err(format!("solv atoms {} != {}", solv.atoms.len(), spec.n_atoms));
    }

    // The molecule identity carries the protomer id (`name.prot_id`), as Divya's driver names it; the
    // tarball member is `name.prot_id.C`, so `Ligand` keeps the bare name + `prot_id` separately.
    let full = format!("{name}.{prot_id}");
    let mol = assemble(&full, smiles, &p, confs);
    let db2 = build_blocks_db2(&mol, &solv, &block_lens(&blocks))?;
    let charges: Vec<f64> = solv.atoms.iter().map(|a| a.charge).collect();
    let mol2 = crate::write_mol2(&mol, &charges, 0);
    let charge = ChargeLetter::new(solv.charge.round() as i32).map_err(|e| e.to_string())?;
    Ok(crate::Ligand { name: name.clone(), prot_id, charge, mol2, db2 })
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

/// Iterate `mol`'s core blocks (`build_ligands.py` lines 534-620: her `mol2db2` runs once per block
/// of conformers that share one core's frozen pin_atoms, a rigid component). `block_lens` partitions
/// `mol`'s conformers, in order, into those blocks ([`block_lens`] of the embedded recipe); a
/// partition that does not cover the conformers exactly, or has an empty block, is an error.
/// `f(block_index, &TypedMol)` sees a molecule whose coords/strain are that block's window; the
/// topology is cloned ONCE (not per block — finding 5), and strain is rounded to the S-record precision
/// (idempotent under the `%.3f` S-record, so this doesn't change the db2 while letting the mol2 Strain
/// comment's `%.6f` round-trip land on the same value — finding 2). `mol` should already be
/// coord-rounded ([`round_coords_to_file_precision`]).
pub fn for_each_block<E: From<String>>(
    mol: &TypedMol,
    block_lens: &[usize],
    mut f: impl FnMut(usize, &TypedMol) -> Result<(), E>,
) -> Result<(), E> {
    let nconf = mol.atom_xyz.len();
    if block_lens.contains(&0) || block_lens.iter().sum::<usize>() != nconf {
        return Err(E::from(format!("block lengths {block_lens:?} do not partition {nconf} conformers")));
    }
    let has_strain = mol.input_total_strain.len() == nconf;
    let has_hydrogens = mol.input_hydrogens.len() == nconf;
    let mut bmol = mol.clone(); // topology cloned once; per-conf fields sliced to the block below
    let mut lo = 0usize;
    for (b, &len) in block_lens.iter().enumerate() {
        let hi = lo + len;
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
    }
    Ok(())
}

/// Build the molecule's db2 the way `build_ligands.py` does: per core block (`block_lens`, as in
/// [`for_each_block`]), from FILE-precision inputs, concatenating the per-block db2 text. See
/// [`round_solv_to_file_precision`].
pub fn build_blocks_db2(mol: &TypedMol, solv: &SolvFile, block_lens: &[usize]) -> Result<String, String> {
    let mut solv = solv.clone();
    round_solv_to_file_precision(&mut solv);
    let mut mol = mol.clone();
    round_coords_to_file_precision(&mut mol);
    let mut out = String::new();
    for_each_block(&mol, block_lens, |_b, bmol| {
        out.push_str(&bb_db2::write_entry(&crate::build(bmol, &solv)?));
        Ok::<(), String>(())
    })?;
    Ok(out)
}

/// One molecule to build: an input line's SMILES, its checked name, and its protomer id.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InputMolecule {
    pub smiles: String,
    pub name: MoleculeName,
    pub prot_id: i32,
}

/// Why an input line yields no molecule to build.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum InputError {
    #[error("line {line}: no name after the SMILES")]
    NoName { line: usize },
    #[error("line {line}: {source}")]
    Name { line: usize, source: NameError },
    #[error("line {line}: protomer id `{field}` is not a whole number from 0 to {}", i32::MAX)]
    ProtId { line: usize, field: String },
    #[error("line {line}: protomer id {from} and every id above it are taken for this name")]
    NoFreeProtId { line: usize, from: i32 },
}

/// The protomer id a line gets: `from` if no earlier line of the same name has it, else the next id
/// above it that is free. `None` only when every id from `from` up is taken.
fn claim_prot_id(taken: &mut std::collections::HashSet<i32>, from: i32) -> Option<i32> {
    let mut id = from;
    while taken.contains(&id) {
        id = id.checked_add(1)?;
    }
    taken.insert(id);
    Some(id)
}

/// Read a `bb-build` input: `smiles name [prot_id]` per line (further fields ignored, blank lines
/// skipped). Every other line yields its molecule, or why it cannot be built — so the caller reports
/// and skips just that line. Protomer ids are unique per name: a line keeps its own id (the third
/// field, which `build_ligands.py` reads) unless an earlier line of that name has it, in which case it
/// takes the next free id above it (protomer expansion upstream repeats ids); a line without one takes
/// the lowest free id, so a file of bare names numbers each name `.0`, `.1`, ….
pub fn read_input(text: &str) -> Vec<Result<InputMolecule, InputError>> {
    let mut taken: std::collections::HashMap<MoleculeName, std::collections::HashSet<i32>> = std::collections::HashMap::new();
    text.lines()
        .enumerate()
        .filter_map(|(i, text)| {
            let mut fields = text.split_whitespace();
            let smiles = fields.next()?; // blank line
            let line = i + 1;
            Some((|| {
                let name = fields.next().ok_or(InputError::NoName { line })?;
                let name = MoleculeName::new(name).map_err(|source| InputError::Name { line, source })?;
                let from = match fields.next() {
                    None => 0,
                    Some(f) => f.parse::<i32>().ok().filter(|&id| id >= 0).ok_or_else(|| InputError::ProtId { line, field: f.to_string() })?,
                };
                let prot_id = claim_prot_id(taken.entry(name.clone()).or_default(), from).ok_or(InputError::NoFreeProtId { line, from })?;
                Ok(InputMolecule { smiles: smiles.to_string(), name, prot_id })
            })())
        })
        .collect()
}

/// The whole pipeline, RDKit-free: SMILES → the DOCK output tarball bytes (one ligand).
pub fn build_tarball_from_smiles(name: &MoleculeName, smiles: &str, seed: u64, workdir: &std::path::Path) -> Result<Vec<u8>, String> {
    let lig = ligand_from_smiles(name, 0, smiles, seed, workdir)?;
    crate::write_tarball(std::slice::from_ref(&lig)).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::{block_lens, build_blocks_db2, confs_from_blocks, for_each_block, read_input, InputError, InputMolecule};
    use crate::mol::TypedMol;
    use crate::name::NameError;

    /// One atom, `n` conformers: conformer `c` has the atom at x = c, total strain c + 0.12345 and max
    /// strain c + 0.54321.
    fn one_atom_confs(n: usize) -> TypedMol {
        TypedMol {
            atom_xyz: (0..n).map(|c| vec![[c as f64, 0.0, 0.0]]).collect(),
            input_total_strain: (0..n).map(|c| c as f64 + 0.12345).collect(),
            input_max_strain: (0..n).map(|c| c as f64 + 0.54321).collect(),
            input_hydrogens: vec![0; n],
            ..Default::default()
        }
    }

    /// Each block sees exactly its own conformers and their strain, the strain rounded to the db2
    /// S-record's 3 decimals.
    #[test]
    fn for_each_block_gives_each_block_its_own_window() {
        let mut seen = Vec::new();
        for_each_block(&one_atom_confs(5), &[2, 3], |b, bmol| {
            let x: Vec<f64> = bmol.atom_xyz.iter().map(|c| c[0][0]).collect();
            seen.push((b, x, bmol.input_total_strain.clone(), bmol.input_max_strain.clone()));
            Ok::<(), String>(())
        })
        .expect("[2, 3] partitions 5 conformers");
        assert_eq!(
            seen,
            [
                (0, vec![0.0, 1.0], vec![0.123, 1.123], vec![0.543, 1.543]),
                (1, vec![2.0, 3.0, 4.0], vec![2.123, 3.123, 4.123], vec![2.543, 3.543, 4.543]),
            ]
        );
    }

    /// Block lengths that do not cover the conformers exactly, or include an empty block, are rejected
    /// before any block is visited.
    #[test]
    fn for_each_block_rejects_lengths_that_do_not_partition_the_conformers() {
        let mol = one_atom_confs(5);
        for lens in [&[2, 2][..], &[2, 4], &[], &[5, 0], &[0, 5]] {
            let mut visited = false;
            let r = for_each_block(&mol, lens, |_, _| {
                visited = true;
                Ok::<(), String>(())
            });
            assert!(r.is_err(), "{lens:?} accepted for 5 conformers");
            assert!(!visited, "{lens:?}: a block was visited before the lengths were rejected");
        }
    }

    /// `build_ligands.py` writes one db2 entry per core block, so a recipe of `core_seeds` blocks gives
    /// `core_seeds` entries, each closed by its `E` line.
    #[test]
    fn build_blocks_db2_writes_one_entry_per_block() {
        let smiles = "C1CCOCC1";
        let (spec, p) = bb_spec::build_native_perceived(smiles).expect("spec");
        let blocks = bb_embed::embed_recipe(&spec, 210185).expect("the recipe completes");
        let mol = super::assemble("thp", smiles, &p, confs_from_blocks(&blocks, spec.n_atoms));
        let solv = bb_solv::solv::SolvFile {
            name: "thp".into(),
            charge: 0.0,
            area: 0.0,
            tot_diff_pol: 0.0,
            tot_diff_apol: 0.0,
            tot_diff_pol_plus_apol: 0.0,
            atoms: vec![
                bb_solv::solv::SolvAtom { charge: 0.0, diff_pol: 0.0, area: 0.0, diff_apol: 0.0, diff_atomic_solv: 0.0 };
                spec.n_atoms
            ],
        };
        let db2 = build_blocks_db2(&mol, &solv, &block_lens(&blocks)).expect("db2");
        assert_eq!(db2.lines().filter(|l| *l == "E").count(), spec.core_seeds as usize);
    }

    fn ok(r: &Result<InputMolecule, InputError>) -> (&str, &str, i32) {
        let m = r.as_ref().expect("a molecule");
        (m.smiles.as_str(), m.name.as_str(), m.prot_id)
    }

    #[test]
    fn prot_ids_count_by_occurrence() {
        // repeated names get 0,1,2…; distinct names each restart at 0 — no member-name collisions.
        let got = read_input("C mc0001\nCC mc0001\nN mc0002\nO mc0001\nS mc0002\n");
        let ids: Vec<(&str, i32)> = got.iter().map(|r| (ok(r).1, ok(r).2)).collect();
        assert_eq!(ids, [("mc0001", 0), ("mc0001", 1), ("mc0002", 0), ("mc0001", 2), ("mc0002", 1)]);
    }

    /// The agreed rule: a line keeps its own protomer id unless that name already has it, then takes
    /// the next free id above; a line without one takes the lowest free id. (Expected values derive
    /// from that rule as stated, so this covers carrying it out, not the rule itself.)
    #[test]
    fn given_prot_ids_are_kept_and_repeats_take_the_next_free_id() {
        let got = read_input("C x 1\nC x 1\nC x 2\nC x\nC x 1\nN y 0\nN y\nO z 7\n");
        let ids: Vec<(&str, i32)> = got.iter().map(|r| (ok(r).1, ok(r).2)).collect();
        assert_eq!(
            ids,
            [
                ("x", 1), // its own id
                ("x", 2), // repeats 1: next free above
                ("x", 3), // its own 2 is taken now: next free above
                ("x", 0), // no id: lowest free
                ("x", 4), // repeats 1 again: 2 and 3 are taken too (three repeats never collide)
                ("y", 0),
                ("y", 1),
                ("z", 7), // ids need not start at 0
            ]
        );
    }

    #[test]
    fn a_prot_id_that_is_not_a_whole_number_is_reported() {
        let max = i32::MAX.to_string();
        let got = read_input(&format!("C x -1\nC x 1a\nC x 2147483648\nC x 1.0\nC w {max}\nC w {max}\nC x 3\n"));
        for (r, (line, field)) in got.iter().zip([(1, "-1"), (2, "1a"), (3, "2147483648"), (4, "1.0")]) {
            assert_eq!(r, &Err(InputError::ProtId { line, field: field.into() }));
        }
        assert_eq!(ok(&got[4]), ("C", "w", i32::MAX));
        assert_eq!(got[5], Err(InputError::NoFreeProtId { line: 6, from: i32::MAX }));
        assert_eq!(ok(&got[6]), ("C", "x", 3), "rejected lines claim no id");
    }

    #[test]
    fn every_input_line_is_a_molecule_or_a_reported_reason() {
        let got = read_input("C ok 0 extra fields\n\n   \nCC\nN a/b\nO ..\nS ok\n");
        assert_eq!(got.len(), 5, "blank lines yield nothing; every other line yields one result");
        assert_eq!(ok(&got[0]), ("C", "ok", 0));
        assert_eq!(got[1], Err(InputError::NoName { line: 4 }));
        assert_eq!(got[2], Err(InputError::Name { line: 5, source: NameError::Separator("a/b".into()) }));
        assert_eq!(got[3], Err(InputError::Name { line: 6, source: NameError::DotComponent("..".into()) }));
        assert_eq!(ok(&got[4]), ("S", "ok", 1), "a rejected line takes no protomer id");
    }
}
