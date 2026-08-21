//! SYBYL → DOCK vdw-type conversion — a faithful port of `mol2db2_py3_strain/sybyl2dock.py`.
//!
//! A fixed table maps each SYBYL atom type to a DOCK integer vdw type. The one context-dependent
//! case is hydrogen: the table has `H`→6 and `H-C`→7, and a hydrogen counts as `H-C` when it is
//! bonded to a carbon (any type containing `C`). This is the only "special key" (the only default
//! entry containing `-`), so the bond check runs exactly when the atom's SYBYL type is `H`.

use crate::mol::TypedMol;

/// The `convertTypesDefault` table, verbatim — SYBYL type string → DOCK vdw integer.
pub fn convert_type(sybyl: &str) -> Option<i32> {
    Some(match sybyl {
        "C.3" => 5,
        "C.2" | "C.ar" | "C.1" | "C.cat" => 1,
        "N.3" => 10,
        "N.2" | "N.1" | "N.ar" | "N.pl3" | "N.am" => 8,
        "O.3" => 12,
        "O.2" | "O.co2" | "O.spc" | "O.t3p" => 11,
        "S.3" | "S.2" | "S.o" | "S.o2" => 14,
        "P.3" => 13,
        "H" | "H.spc" | "H.t3p" => 6,
        "Br" => 17,
        "Cl" => 16,
        "F" => 15,
        "I" => 18,
        "N.4" => 9,
        "LP" | "Du" | "Du.C" | "ANY" | "HEV" | "HET" | "HAL" | "Cr.oh" | "Cr.th" | "Se" | "Fe"
        | "Cu" | "Sn" | "Mo" | "Mn" | "Co.oh" => 25,
        "Na" | "K" => 19,
        "Ca" => 21,
        "Li" | "Al" | "Mg" => 20,
        "Si" => 24,
        "Zn" => 26,
        _ => return None,
    })
}

/// DOCK vdw type for atom `atom_num` (1-based), reproducing `convertMol2atomNum`.
pub fn convert_atom(mol2: &TypedMol, atom_num: usize) -> Result<i32, String> {
    let sybyl = &mol2.atom_type[atom_num - 1];
    // The only special key is `H-C` (justFirstPart `H`), so the bond branch fires only for type `H`.
    if sybyl == "H" {
        // `H-C`: bonded (1 away) to an atom whose type contains `C` → type 7, else `H` → 6.
        if mol2.bonded_to(atom_num, "C", 1, None) {
            return Ok(7);
        }
    }
    convert_type(sybyl).ok_or_else(|| format!("unknown SYBYL type {sybyl:?}"))
}
