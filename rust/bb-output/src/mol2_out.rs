//! The charge-bearing mol2 output artifact — a faithful port of `mol2amsol.write_mol2` (the writer
//! `process_amsol_mol2.py3.py` uses to emit `output.mol2`, the tarball's mol2 member). It serializes a
//! [`TypedMol`] (a single conformation) plus per-atom charges; in the pipeline those charges are the
//! AMSOL values (our `.solv`), and the residue name is the molecule name. One residue is assumed, as
//! the writer hard-codes. Atom/bond numbering is renumbered sequentially exactly as the source does.

use crate::mol::TypedMol;
use std::fmt::Write as _;

/// Serialize `mol` (conformation `conf`) with `charges` to `output.mol2` text.
pub fn write_mol2(mol: &TypedMol, charges: &[f64], conf: usize) -> String {
    let n = mol.atom_num.len();
    let resname = &mol.name;
    let mut s = String::new();

    // MOLECULE: name, counts (nres hard-coded via a single residue), type, charge type
    s.push_str("@<TRIPOS>MOLECULE\n");
    let _ = writeln!(s, "{}", mol.name);
    let _ = writeln!(s, "{:<5} {:<5} {:<5} 0     0", n, mol.bond_num.len(), 1);
    s.push_str("SMALL\n");
    s.push_str("USER_CHARGES\n");

    // ATOM: "%-6d %-4s %9.4f %9.4f %9.4f %-5s %4s %6s %9.4f" — id name x y z sybyl resid resname charge
    s.push_str("@<TRIPOS>ATOM\n");
    #[allow(clippy::needless_range_loop)] // `i` indexes several parallel arrays (name/type/xyz/charge)
    for i in 0..n {
        let xyz = mol.atom_xyz[conf][i];
        let _ = writeln!(
            s,
            "{:<6} {:<4} {:9.4} {:9.4} {:9.4} {:<5} {:>4} {:>6} {:9.4}",
            i + 1,
            mol.atom_name[i],
            xyz[0],
            xyz[1],
            xyz[2],
            mol.atom_type[i],
            1, // resid (single residue)
            resname,
            charges[i],
        );
    }

    // BOND: "%-5d %-5d %-5d %s" — num a1 a2 type (sequential numbering)
    s.push_str("@<TRIPOS>BOND\n");
    for i in 0..mol.bond_num.len() {
        let _ = writeln!(s, "{:<5} {:<5} {:<5} {}", i + 1, mol.bond_start[i], mol.bond_end[i], mol.bond_type[i]);
    }

    // SUBSTRUCTURE: one residue, first atom = 1
    s.push_str("@<TRIPOS>SUBSTRUCTURE\n");
    let res3: String = resname.chars().take(3).collect();
    let _ = writeln!(s, "{:<3} {:<5} {:<5} RESIDUE    1   A     {:<5} 1", 1, resname, 1, res3);
    s
}

/// Serialize EVERY conformation of `mol` as concatenated TRIPOS blocks — a multimol2, the input
/// `mol2db2.py` consumes (one MOLECULE block per conformer, identical topology, differing coords).
/// Used by the corpus-scale writer gate to feed native's own ensemble to the container's writer.
///
/// Each block carries a per-conformer `Strain_Energy:` header line (total, max) — her `mol2.py`
/// phase-1 parser reads the mol_name line FIRST, then a `Strain`-prefixed line sets
/// `inputTotalStrain=float(tokens[2])`, `inputMaxStrain=float(tokens[3])` and ends the header; absent,
/// it defaults to the 9999.99 sentinel. Written at `%.6f`; the strain reaching here has already been
/// rounded to the S-record display precision (`STRAIN_DP`) by [`crate::assemble::for_each_block`], so the
/// `%.6f` mol2 value and the `%.3f` db2 S-record round-trip to the same number (no half-way flip). `mol`'s
/// `input_*_strain` must already be sliced to this block's conformers.
pub fn write_mol2_all(mol: &TypedMol, charges: &[f64]) -> String {
    (0..mol.atom_xyz.len())
        .map(|c| {
            let block = write_mol2(mol, charges, c);
            if c < mol.input_total_strain.len() {
                let mut it = block.splitn(3, '\n');
                let l0 = it.next().unwrap_or(""); // @<TRIPOS>MOLECULE
                let l1 = it.next().unwrap_or(""); // mol name (her parser sets self.name from this)
                let rest = it.next().unwrap_or("");
                format!(
                    "{l0}\n{l1}\nStrain_Energy: 0 {:.6} {:.6}\n{rest}",
                    mol.input_total_strain[c], mol.input_max_strain[c]
                )
            } else {
                block
            }
        })
        .collect()
}
