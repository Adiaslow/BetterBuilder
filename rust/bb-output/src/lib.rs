//! bb-output — produce the DOCK build-pipeline **output artifacts** from a typed in-memory molecule.
//!
//! The deliverable is the tarball DOCK consumes: a charge-bearing mol2 and a `.db2`. Both are
//! serializations of the *same* [`mol::TypedMol`] (atoms + SYBYL types, bonds, conformations) plus its
//! `.solv` — we never derive one from the other. [`build`] assembles a [`bb_db2::Db2Entry`] (which
//! [`bb_db2`] serializes to the exact db2 bytes) by typing every atom ([`sybyl2dock`], [`color`]) and
//! constructing the conformer [`hierarchy`]; [`write_mol2`] writes the mol2; [`assemble`] runs the
//! pipeline from SMILES to both; [`tarball`] packs them under each molecule's [`name`]. The molecule
//! comes from BetterBuilder's own pipeline in production; the gate feeds an identical one from the
//! container's authoritative parse, so the output is checked byte-for-byte against the pipeline's own
//! db2 builder (and per-stage against its `hierarchy.Hierarchy` internals).

pub mod assemble;
pub mod buckets;
pub mod clash;
pub mod color;
pub mod hierarchy;
pub mod mol;
pub mod mol2_out;
pub mod name;
pub mod solvation;
pub mod sybyl;
pub mod sybyl2dock;
pub mod tarball;
pub mod unionfind;
pub mod zinc22;

pub use mol2_out::{write_mol2, write_mol2_all};
pub use name::MoleculeName;
pub use tarball::{write_tarball, Ligand};

use bb_db2::model::{Atom, Bond, Conf, Coord, Db2Entry, Rigid, Set};
use bb_solv::solv::SolvFile;
use mol::TypedMol;

/// Distances below which two positions are "the same" (`mol2db2.py` `--disttol` default).
pub const DISTTOL: f64 = 0.001;

/// Build a complete [`Db2Entry`] from a typed molecule and its solvation ([`bb_solv::solv::SolvFile`],
/// the single `.solv` model shared with the solvation stage). `mol2.atom_bonds` must be populated
/// ([`TypedMol::build_adjacency`]).
pub fn build(mol2: &TypedMol, solv: &SolvFile) -> Result<Db2Entry, String> {
    let n = mol2.atom_num.len();
    if solv.atoms.len() != n {
        return Err(format!("solv atoms {} != mol2 atoms {}", solv.atoms.len(), n));
    }

    // per-atom DOCK vdw type + pharmacophore color
    let mut dock_num = Vec::with_capacity(n);
    let mut color_num = Vec::with_capacity(n);
    for &num in &mol2.atom_num {
        dock_num.push(sybyl2dock::convert_atom(mol2, num)?);
        color_num.push(color::convert_atom(mol2, num));
    }

    let built = hierarchy::build(mol2, DISTTOL);

    // A records
    let atoms = (0..n)
        .map(|i| Atom {
            num: mol2.atom_num[i] as i32,
            name: mol2.atom_name[i].clone(),
            sybyl: mol2.atom_type[i].clone(),
            vdwtype: dock_num[i],
            color: color_num[i],
            charge: solv.atoms[i].charge,
            polar_solv: solv.atoms[i].diff_pol,
            apolar_solv: solv.atoms[i].diff_apol,
            total_solv: solv.atoms[i].diff_atomic_solv,
            surface: solv.atoms[i].area,
        })
        .collect();
    // B records
    let bonds = (0..mol2.bond_num.len())
        .map(|i| Bond {
            num: mol2.bond_num[i] as i32,
            a1: mol2.bond_start[i] as i32,
            a2: mol2.bond_end[i] as i32,
            btype: mol2.bond_type[i].clone(),
        })
        .collect();
    // X records (flexible coordinates)
    let coords = (0..built.out_atoms)
        .map(|k| Coord {
            num: (k + 1) as i32,
            atom: (built.out_atom_orig_atom[k] + 1) as i32,
            conf: built.out_atom_conf_num[k] as i32,
            xyz: built.out_atom_xyz[k],
        })
        .collect();
    // R records (rigid core, heavy atoms)
    let rigid = built
        .heavy_rigid_atom_nums
        .iter()
        .enumerate()
        .map(|(i, &atom)| Rigid {
            num: (i + 1) as i32,
            color: color_num[atom],
            xyz: mol2.atom_xyz[0][atom],
        })
        .collect();
    // C records (conformations)
    let confs = (0..built.num_confs)
        .map(|c| Conf {
            num: (c + 1) as i32,
            start: (built.conf_num_atom_list[c][0] + 1) as i32,
            end: (built.conf_num_atom_list[c][1] + 1) as i32,
        })
        .collect();
    // S records (sets = input conformers), sorted by set index
    let mut sets = Vec::new();
    for set_idx in 0..built.set_to_confs.len() {
        let confs_1based: Vec<i32> = built.set_to_confs[set_idx].iter().map(|&c| (c + 1) as i32).collect();
        if confs_1based.is_empty() {
            continue;
        }
        let member_lines: Vec<Vec<i32>> = confs_1based.chunks(8).map(<[i32]>::to_vec).collect();
        sets.push(Set {
            num: (set_idx + 1) as i32,
            broken: i32::from(built.broken_sets.contains(&set_idx)),
            hydrogens: mol2.input_hydrogens[set_idx],
            total_strain: mol2.input_total_strain[set_idx],
            max_strain: mol2.input_max_strain[set_idx],
            member_lines,
        });
    }

    Ok(Db2Entry {
        color_table: Vec::new(),
        name: mol2.name.clone(),
        protname: mol2.prot_name.clone(),
        total_charge: solv.charge,
        total_polar_solv: solv.tot_diff_pol,
        total_apolar_solv: solv.tot_diff_apol,
        total_solv: solv.tot_diff_pol_plus_apol,
        surface_area: solv.area,
        smiles: mol2.smiles.clone(),
        longname: mol2.longname.clone(),
        m5: 999.999,
        extra_m: Vec::new(),
        atoms,
        bonds,
        coords,
        rigid,
        confs,
        sets,
        cluster_lines: Vec::new(),
    })
}
