//! Native SYBYL atom typing — a faithful port of `rdkit_mol2_writer.py`'s `_assign_sybyl_types_to_mol`
//! (the pipeline's RDKit-based typer), operating on [`bb_perceive`]'s perception (already gated against
//! RDKit). Two stages, exactly as the source: functional-group special cases assigned first
//! (`_pre_assign_atom_types`), then a per-atom guesser for the rest (`_guess_sybyl_atom_type`). Because
//! our atom order matches RDKit's, the result is gated per-atom against `rdkit_mol2_writer`'s own
//! `3d/name.mol2` output.

use bb_perceive::hybrid::Hybridization;
use bb_perceive::smarts_match::Perceived;

fn elem(z: u8) -> &'static str {
    match z {
        1 => "H",
        6 => "C",
        7 => "N",
        8 => "O",
        9 => "F",
        15 => "P",
        16 => "S",
        17 => "Cl",
        35 => "Br",
        53 => "I",
        _ => "*",
    }
}

/// Per-atom neighbor list: `(neighbor, bond_order)`, using the given (Kekulé) bond orders.
fn adjacency(p: &Perceived, orders: &[u8]) -> Vec<Vec<(usize, u8)>> {
    let mut adj = vec![Vec::new(); p.atomic_numbers.len()];
    for (bi, &(a, b)) in p.bonds.iter().enumerate() {
        adj[a].push((b, orders[bi]));
        adj[b].push((a, orders[bi]));
    }
    adj
}

/// SYBYL type for each atom, in perception (= RDKit) order. Convenience over [`assign_sybyl_full`].
pub fn assign_sybyl(p: &Perceived) -> Vec<String> {
    assign_sybyl_full(p).0
}

/// mol2 atom name: element symbol + 1-based atom index (`O1`, `C2`, … `H10`), as `rdkit_mol2_writer`
/// names them.
pub fn atom_names(p: &Perceived) -> Vec<String> {
    p.atomic_numbers.iter().enumerate().map(|(i, &z)| format!("{}{}", elem(z), i + 1)).collect()
}

/// SYBYL bond type for each bond, reproducing `rdkit_mol2_writer._bond_section`: `ar` when both ends
/// are `.ar` (or the carboxylate C.2/O.co2 pair), `am` for amide C–N bonds, else the integer order.
pub fn bond_types(p: &Perceived, types: &[String], amide_bonds: &[(usize, usize)]) -> Vec<String> {
    let kek = bb_perceive::kekulize::kekule_bond_orders(p);
    p.bonds
        .iter()
        .enumerate()
        .map(|(bi, &(a, b))| {
            let (ta, tb) = (&types[a], &types[b]);
            let both_ar = ta.ends_with(".ar") && tb.ends_with(".ar");
            let carboxylate = matches!((ta.as_str(), tb.as_str()), ("C.2", "O.co2") | ("O.co2", "C.2"));
            if both_ar || carboxylate {
                "ar".into()
            } else if amide_bonds.contains(&(a.min(b), a.max(b))) {
                "am".into()
            } else {
                kek[bi].to_string()
            }
        })
        .collect()
}

/// SYBYL type for each atom plus the amide C–N bonds recorded while typing. Aromatic bonds are
/// Kekulé-resolved first (`bb_perceive::kekulize`), reproducing RDKit's
/// `Chem.Kekulize(clearAromaticFlags=True)` before typing — so `aro6` and the guesser see the same
/// single/double structure RDKit's typer does.
pub fn assign_sybyl_full(p: &Perceived) -> (Vec<String>, Vec<(usize, usize)>) {
    let mut amide_bonds: Vec<(usize, usize)> = Vec::new();
    let n = p.atomic_numbers.len();
    let z = &p.atomic_numbers;
    let kek = bb_perceive::kekulize::kekule_bond_orders(p);
    let adj = adjacency(p, &kek);
    // Kekulé order between two atoms, for the aro6 ring-alternation test
    let mut bond_of: std::collections::HashMap<(usize, usize), u8> = std::collections::HashMap::new();
    for (bi, &(a, b)) in p.bonds.iter().enumerate() {
        bond_of.insert((a, b), kek[bi]);
        bond_of.insert((b, a), kek[bi]);
    }
    let deg = |i: usize| adj[i].len();
    // number of H neighbours (post-AddHs, so explicit)
    let h_neighbors = |i: usize| adj[i].iter().filter(|&&(j, _)| z[j] == 1).count();

    let mut t: Vec<Option<String>> = vec![None; n];
    let set = |t: &mut Vec<Option<String>>, i: usize, v: &str| t[i] = Some(v.to_string());

    // --- pre-assign: functional groups, in source order (later overwrites earlier) ---
    for c in 0..n {
        if z[c] != 6 || deg(c) != 3 {
            continue;
        }
        let o_dbl = adj[c].iter().find(|&&(j, o)| z[j] == 8 && o == 2).map(|&(j, _)| j);
        let Some(o) = o_dbl else { continue };
        // amide [C;X3](=O)N: any N neighbour
        for &(nn, _) in &adj[c] {
            if z[nn] == 7 {
                set(&mut t, c, "C.2");
                set(&mut t, o, "O.2");
                set(&mut t, nn, "N.am");
                amide_bonds.push((c.min(nn), c.max(nn)));
            }
        }
        // carboxylate [C;X3](=O)[O-]: single-bonded O with -1 charge
        for &(oo, o_ord) in &adj[c] {
            if z[oo] == 8 && oo != o && o_ord == 1 && p.charges[oo] == -1 {
                set(&mut t, c, "C.2");
                set(&mut t, o, "O.co2");
                set(&mut t, oo, "O.co2");
            }
        }
        // carboxylic acid [C;X3](=O)[OH1]: single-bonded O bearing exactly one H
        for &(oo, o_ord) in &adj[c] {
            if z[oo] == 8 && oo != o && o_ord == 1 && h_neighbors(oo) == 1 {
                set(&mut t, c, "C.2");
                set(&mut t, o, "O.2");
                set(&mut t, oo, "O.2");
            }
        }
    }
    // phosphate [P](=O)(O)(O)(O): P with a =O and three more O
    for pp in 0..n {
        if z[pp] != 15 {
            continue;
        }
        let os: Vec<usize> = adj[pp].iter().filter(|&&(j, _)| z[j] == 8).map(|&(j, _)| j).collect();
        let o_dbl = adj[pp].iter().find(|&&(j, o)| z[j] == 8 && o == 2).map(|&(j, _)| j);
        if let Some(od) = o_dbl {
            if os.len() >= 4 {
                set(&mut t, pp, "P.3");
                set(&mut t, od, "O.2");
                for o in os {
                    if o != od {
                        set(&mut t, o, "O.3");
                    }
                }
            }
        }
    }
    // sulfone [S;X4](=O)(=O)
    for s in 0..n {
        if z[s] != 16 || deg(s) != 4 {
            continue;
        }
        let dbl_os: Vec<usize> = adj[s].iter().filter(|&&(j, o)| z[j] == 8 && o == 2).map(|&(j, _)| j).collect();
        if dbl_os.len() >= 2 {
            set(&mut t, s, "S.o2");
            for o in dbl_os {
                set(&mut t, o, "O.2");
            }
        }
    }
    // sulfoxide [S;X3;!$(S(=O)=O)]=O — only when both still unassigned
    for s in 0..n {
        if z[s] != 16 || deg(s) != 3 {
            continue;
        }
        let dbl_os: Vec<usize> = adj[s].iter().filter(|&&(j, o)| z[j] == 8 && o == 2).map(|&(j, _)| j).collect();
        if dbl_os.len() == 1 && t[s].is_none() && t[dbl_os[0]].is_none() {
            set(&mut t, s, "S.o");
            set(&mut t, dbl_os[0], "O.2");
        }
    }
    // aro6: a 6-ring whose bonds ALTERNATE single/double in the Kekulé form (RDKit's
    // `A1=A-A=A-A=A1`). O -> O.3 (pyrylium), else symbol + ".ar".
    for ring in &p.rings {
        if ring.len() != 6 {
            continue;
        }
        let ords: Option<Vec<u8>> = (0..6).map(|k| bond_of.get(&(ring[k], ring[(k + 1) % 6])).copied()).collect();
        let Some(ords) = ords else { continue }; // ring atoms not consecutively bonded
        let alternates = (0..6).all(|k| {
            let (a, b) = (ords[k], ords[(k + 1) % 6]);
            (a == 1 || a == 2) && (b == 1 || b == 2) && a != b
        });
        if alternates {
            for &a in ring {
                if z[a] == 8 {
                    set(&mut t, a, "O.3");
                } else {
                    set(&mut t, a, &format!("{}.ar", elem(z[a])));
                }
            }
        }
    }

    // --- guesser for anything still unassigned ---
    for i in 0..n {
        if t[i].is_some() {
            continue;
        }
        let sym = elem(z[i]);
        let ch = p.charges[i];
        let hyb = p.hybridization[i];
        let guess = match sym {
            "P" => "P.3".to_string(),
            "F" | "Cl" | "Br" | "I" | "H" => sym.to_string(),
            "C" => {
                if ch == 1 {
                    "C.cat".into()
                } else {
                    match hyb {
                        Hybridization::Sp3 => "C.3",
                        Hybridization::Sp2 => "C.2",
                        Hybridization::Sp => "C.1",
                        _ => "C.3",
                    }
                    .into()
                }
            }
            "N" => {
                if ch == 1 && deg(i) == 4 {
                    "N.4".into()
                } else {
                    match hyb {
                        Hybridization::Sp => "N.1".into(),
                        // source: N.2 if any bond order >= 2, else N.pl3
                        Hybridization::Sp2 => {
                            if adj[i].iter().any(|&(_, o)| o >= 2) {
                                "N.2".into()
                            } else {
                                "N.pl3".into()
                            }
                        }
                        _ => "N.3".into(),
                    }
                }
            }
            "O" => {
                if hyb == Hybridization::Sp2 {
                    if ch == 1 {
                        "O.3".into()
                    } else {
                        "O.2".into()
                    }
                } else {
                    "O.3".into()
                }
            }
            "S" => {
                if hyb == Hybridization::Sp2 {
                    "S.2".into()
                } else {
                    "S.3".into()
                }
            }
            other => other.to_string(),
        };
        t[i] = Some(guess);
    }

    (t.into_iter().map(|x| x.unwrap()).collect(), amide_bonds)
}
