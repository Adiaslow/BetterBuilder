//! Assembling a [`MoleculeSpec`] entirely from the pure-Rust ports — no RDKit.
//!
//! Every field is sourced from [`bb_perceive`] (perception, experimental torsions, 1-3 angles, the
//! two-stage recipe) and [`bb_bounds`] (the full `setTopolBounds` matrix, chirality, stereo, UFF
//! inversion coefficients), each gated against RDKit. The bounds are triangle-smoothed the way the
//! pipeline expects (`triangleSmoothBounds`, tol 0), reproduced by
//! [`bb_core::smooth::triangle_smooth_f64`]. This is the RDKit-free replacement for what the
//! `bb-rdkit` bridge used to build.

use bb_core::{Angle, ChiralSet, ExpTorsion, Improper, MoleculeSpec, StereoDoubleBond};
use bb_perceive::smarts_match::Perceived;

/// Build a [`MoleculeSpec`] entirely from the pure-Rust ports — no RDKit bridge. Every field is
/// sourced from `bb-perceive` (perception, torsions, angles, recipe) and `bb-bounds` (bounds
/// matrix, chirality, stereo, UFF inversion coefficients), each gated against RDKit.
pub fn build_native(smiles: &str) -> Result<MoleculeSpec, String> {
    build_native_perceived(smiles).map(|(spec, _)| spec)
}

/// Like [`build_native`], but also returns the [`Perceived`] the spec was built from. The downstream
/// pipeline (typing, strain) needs the same perception, so this hands back the one already computed
/// here rather than making the caller re-perceive — one perception is the single source of truth for
/// the molecule's chemistry, shared by the embedding spec and the output serializers.
pub fn build_native_perceived(smiles: &str) -> Result<(MoleculeSpec, Perceived), String> {
    let g = bb_perceive::parse(smiles).map_err(|e| e.to_string())?;
    let counts = bb_perceive::hydrogens::counts(&g);
    let ex = bb_perceive::addhs::add_hs(&g, &counts);
    let mol = bb_perceive::smarts_match::perceive(smiles).map_err(|e| e.to_string())?;
    let rings = bb_perceive::sssr::symmetrized_sssr_masked(g.atoms.len(), &g.bonds, &g.dative_donor);
    let n = mol.atomic_numbers.len();

    // bounds: full setTopolBounds, triangle-smoothed in f64 then cast once — matching the bridge,
    // which smooths the `double` matrix before casting to f32 (byte-identical spec bounds).
    let mut raw = bb_bounds::bounds::bounds_full_checked(&mol)?;
    // Keep the pre-smoothing matrix: the core-pin recipe tightens it (not the smoothed `bounds`) to
    // match RDKit's coordMap path (see MoleculeSpec::raw_bounds). f64 for the embed, f32 for symmetry.
    let raw_bounds_f64: Vec<f64> = raw.clone();
    let raw_bounds: Vec<f32> = raw.iter().map(|&x| x as f32).collect();
    bb_core::smooth::triangle_smooth_f64(&mut raw, n, 0.0);
    let bounds: Vec<f32> = raw.iter().map(|&x| x as f32).collect();
    // The embed reads these f64 bounds (via ub64/lb64) to match RDKit's force-field precision.
    let bounds_f64: Vec<f64> = raw;

    // experimental torsions + basic-knowledge terms; impropers expanded to three permutation contribs
    let lib = bb_perceive::torsion_lib::patterns();
    let asts = bb_perceive::torsions::compile(&lib);
    let (mut ts, mut done) = bb_perceive::torsions::assign(&mol, &lib, &asts);
    let (bk, imps) = bb_perceive::torsions::basic_knowledge(&mol, &mut done);
    ts.extend(bk);
    let exp_torsions = ts
        .iter()
        .map(|t| ExpTorsion {
            atoms: [
                t.atoms[0] as u32,
                t.atoms[1] as u32,
                t.atoms[2] as u32,
                t.atoms[3] as u32,
            ],
            v: t.v,
            signs: t.signs,
        })
        .collect();
    const PERM: [[usize; 3]; 3] = [[0, 2, 3], [0, 3, 2], [2, 3, 0]];
    let mut impropers = Vec::new();
    for imp in &imps {
        let (res, c0, c1, c2) =
            bb_bounds::uff::calc_inversion_coefficients(imp.atomic_num, imp.is_c_bound_to_sp2_o);
        let fc = res * 10.0;
        for p in &PERM {
            impropers.push(Improper {
                atoms: [
                    imp.atoms[p[0]] as u32,
                    imp.atoms[1] as u32,
                    imp.atoms[p[1]] as u32,
                    imp.atoms[p[2]] as u32,
                ],
                c0,
                c1,
                c2,
                fc,
            });
        }
    }

    // 1-3 angles (collectBondsAndAngles), over the post-AddHs molecule
    let mut border: Vec<u8> = g
        .bond_kinds
        .iter()
        .enumerate()
        .map(|(bi, k)| {
            let (x, y) = g.bonds[bi];
            let ar = g.atoms[x].aromatic_as_written && g.atoms[y].aromatic_as_written;
            let c = bb_perceive::hydrogens::contribution(*k, ar);
            if (c - 2.0).abs() < 1e-9 {
                2
            } else if (c - 3.0).abs() < 1e-9 {
                3
            } else {
                1
            }
        })
        .collect();
    border.resize(ex.bonds.len(), 1);
    let mut degs = vec![0usize; ex.atomic_numbers.len()];
    for &(x, y) in &ex.bonds {
        degs[x] += 1;
        degs[y] += 1;
    }
    let angles = bb_perceive::angles::collect(&ex, &border, &degs)
        .iter()
        .map(|a| Angle {
            atoms: [a.atoms[0] as u32, a.atoms[1] as u32, a.atoms[2] as u32],
            triple: a.triple,
        })
        .collect();

    // two-stage recipe: seed/sidechain counts and the pinned core
    let exo = bb_perceive::recipe::count_exo_rotatable(&g, &rings);
    let (core_seeds, sidechain_confs) = bb_perceive::recipe::recipe(exo);
    let pin_atoms: Vec<u32> = bb_perceive::recipe::pin_atoms(&g, &rings, &mol.bond_aromatic)
        .iter()
        .map(|&a| a as u32)
        .collect();

    // chirality and stereo constraints
    let conv = |c: &bb_bounds::chirality::ChiralSet| ChiralSet {
        center: c.center,
        atoms: c.atoms,
        // exact: the chiral volume bounds are integer constants (2/5/100/0), so the widen is lossless
        vol_lo: c.vol_lo as f64,
        vol_hi: c.vol_hi as f64,
        fused_small_rings: c.fused_small_rings,
    };
    let (chiral, tetra) = bb_bounds::chirality::find_chiral_sets(&mol);
    let (dbe, sdb) = bb_bounds::chirality::find_double_bonds(&mol);

    let spec = MoleculeSpec {
        n_atoms: n,
        dim: 4,
        atomic_numbers: mol.atomic_numbers.clone(),
        formal_charge: mol.charges.iter().map(|&c| i32::from(c)).sum(),
        bounds,
        bounds_f64,
        raw_bounds,
        raw_bounds_f64,
        chiral_sets: chiral.iter().map(conv).collect(),
        tetrahedral_centers: tetra.iter().map(conv).collect(),
        exp_torsions,
        impropers,
        // Emit each bond as (begin, end). RDKit's `metalBondCleanup` orients a dative bond
        // donor→metal (setBeginAtom(donor)/setEndAtom(metal)), so a dative bond leads with its
        // donor; all others keep their parsed order.
        bonds: mol
            .bonds
            .iter()
            .enumerate()
            .map(|(bi, &(a, b))| match g.dative_donor.get(bi).copied().flatten() {
                Some(d) if d == b => [b as u32, a as u32],
                _ => [a as u32, b as u32],
            })
            .collect(),
        angles,
        bounds_force_scaling: 1.0,
        pin_atoms,
        core_seeds,
        sidechain_confs,
        double_bond_ends: dbe,
        stereo_double_bonds: sdb
            .iter()
            .map(|s| StereoDoubleBond {
                atoms: s.atoms,
                sign: s.sign,
            })
            .collect(),
    };
    Ok((spec, mol))
}
