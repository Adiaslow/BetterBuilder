//! Pharmacophore color assignment — a faithful port of `mol2db2_py3_strain/atom_color_table.py`.
//!
//! Seven color classes; the default is `neutral`. A fixed rule list is applied in order and the
//! **last** matching rule wins. A rule is either `(sybyl_prefix → color)` (the atom's SYBYL type
//! starts with the prefix), or a bond rule `(sybyl_prefix, dist, other, color)`: `dist == -1` fires
//! when the atom is *not* bonded (1 away) to `other`, and `dist > 0` fires when it *is* bonded
//! exactly `dist` away. Bond checks go through [`TypedMol::bonded_to`].

use crate::mol::TypedMol;

fn color_int(name: &str) -> i32 {
    match name {
        "positive" => 1,
        "negative" => 2,
        "acceptor" => 3,
        "donor" => 4,
        "ester_o" => 5,
        "amide_o" => 6,
        _ => 7, // neutral
    }
}

enum Rule {
    /// prefix → color
    Prefix(&'static str, &'static str),
    /// prefix, distance, other-type, color (distance -1 = "not bonded 1 away")
    Bond(&'static str, i32, &'static str, &'static str),
}

/// `rulesTableDefault`, in order.
const RULES: &[Rule] = &[
    Rule::Prefix("N.4", "positive"),
    Rule::Prefix("O.co2", "negative"),
    Rule::Prefix("O.2", "acceptor"),
    Rule::Prefix("O.3", "acceptor"),
    Rule::Prefix("S.2", "acceptor"),
    Rule::Prefix("N.ar", "acceptor"),
    Rule::Bond("P.3", 1, "O.co2", "negative"),
    Rule::Bond("S.o2", 1, "O.co2", "negative"),
    Rule::Bond("N.2", 1, "H", "donor"),
    Rule::Bond("N.am", 1, "H", "donor"),
    Rule::Bond("N.pl3", 1, "H", "donor"),
    Rule::Bond("O.3", 1, "H", "donor"),
    Rule::Bond("N.ar", -1, "H", "acceptor"),
    Rule::Bond("N.ar", -1, "C.3", "acceptor"),
    Rule::Bond("N.ar", 1, "H", "donor"),
    Rule::Bond("O.3", 1, "H", "donor"),
    Rule::Bond("O.2", 2, "O.3", "ester_o"),
    Rule::Bond("O.2", 2, "N.pl3", "amide_o"),
    Rule::Bond("O.2", 2, "N.am", "amide_o"),
    Rule::Bond("O.2", 2, "N.3", "amide_o"),
];

/// Color integer for atom `atom_num` (1-based), reproducing `convertMol2color`.
pub fn convert_atom(mol2: &TypedMol, atom_num: usize) -> i32 {
    let actual = &mol2.atom_type[atom_num - 1];
    let mut last = "neutral";
    for rule in RULES {
        match rule {
            Rule::Prefix(prefix, color) => {
                if actual.starts_with(prefix) {
                    last = color;
                }
            }
            Rule::Bond(prefix, dist, other, color) => {
                if actual.starts_with(prefix) {
                    if *dist == -1 {
                        if !mol2.bonded_to(atom_num, other, 1, None) {
                            last = color;
                        }
                    } else if mol2.bonded_to(atom_num, other, *dist as usize, None) {
                        last = color;
                    }
                }
            }
        }
    }
    color_int(last)
}
