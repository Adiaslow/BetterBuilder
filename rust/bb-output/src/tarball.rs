//! Assemble the pipeline's output tarball — a faithful reproduction of `build_ligands.py`'s
//! `output.tar.gz`. Each built ligand contributes two members, `<dir>/<name>.<prot_id>.<C>.mol2` (the
//! charge-bearing mol2) and `<dir>/<name>.<prot_id>.<C>.db2.gz` (the db2 text — the driver stores it
//! *uncompressed* under a `.gz` name, so we match that), where `C` is the [`ChargeLetter`] and `<dir>` is
//! the ligand's [`archive_dir`]. The archive is a gzipped tar. The driver's own tar is left truncated by
//! a close-ordering bug; ours is finalized properly — an intentional deviation, a fix not a divergence.
//!
//! A [`Ligand`] holds only checked parts (its [`MoleculeName`] and [`ChargeLetter`]), so writing it
//! cannot fail for one molecule; the only failures left are the archive's own I/O.

use crate::name::MoleculeName;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::Write;

/// One ligand's output: identity + the two already-built artifact texts.
pub struct Ligand {
    pub name: MoleculeName,
    pub prot_id: i32,
    pub charge: ChargeLetter,
    pub mol2: String,
    pub db2: String,
}

/// The letter a member name gives the ligand's net formal charge: `chr(78 + charge)` (N neutral,
/// M −1, O +1, …), as the ZINC-22 3D build names its members (docking-org/zinc22-3d
/// `generate/build_ligands.py`). It is defined for the charges that map to a capital letter, −13 (A)
/// to +12 (Z); any other charge would put punctuation, or a `/`, into a file name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ChargeLetter(char);

/// A net formal charge with no [`ChargeLetter`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
#[error("net formal charge {0} has no member-name letter (defined for -13 to +12)")]
pub struct ChargeOutOfRange(pub i32);

impl ChargeLetter {
    pub fn new(formal_charge: i32) -> Result<Self, ChargeOutOfRange> {
        if (-13..=12).contains(&formal_charge) {
            // in range, so 78 + charge is 65..=90, `A`..=`Z`
            Ok(ChargeLetter(char::from((78 + formal_charge) as u8)))
        } else {
            Err(ChargeOutOfRange(formal_charge))
        }
    }

    pub fn as_char(self) -> char {
        self.0
    }
}

/// `mol_fullname` = `name.prot_id.C`.
pub fn full_name(name: &MoleculeName, prot_id: i32, charge: ChargeLetter) -> String {
    format!("{name}.{prot_id}.{}", charge.as_char())
}

/// Where a ligand's members sit in the archive. An identifier that itself encodes a storage location
/// is filed there — today that is a ZINC-22 id with a tranche ([`crate::zinc22::tranche_dir`]); every
/// other name goes to the top level, `.`, as `build_ligands.py` files it.
pub fn archive_dir(name: &MoleculeName) -> String {
    crate::zinc22::tranche_dir(name.as_str()).unwrap_or_else(|| ".".to_string())
}

fn add(builder: &mut tar::Builder<impl Write>, name: &str, data: &[u8]) -> std::io::Result<()> {
    let mut header = tar::Header::new_gnu();
    header.set_size(data.len() as u64);
    header.set_mode(0o644);
    header.set_mtime(0); // deterministic: the tar carries no wall-clock time
    header.set_cksum();
    builder.append_data(&mut header, name, data)
}

/// Build the gzipped output tarball for a set of ligands, returning its bytes. Members are emitted in
/// the given order; for each ligand a `.mol2` then a `.db2.gz` (matching the driver's per-ligand order).
pub fn write_tarball(ligands: &[Ligand]) -> std::io::Result<Vec<u8>> {
    let gz = GzEncoder::new(Vec::new(), Compression::default());
    let mut builder = tar::Builder::new(gz);
    for lig in ligands {
        let dir = archive_dir(&lig.name);
        let fullname = full_name(&lig.name, lig.prot_id, lig.charge);
        add(&mut builder, &format!("{dir}/{fullname}.mol2"), lig.mol2.as_bytes())?;
        add(&mut builder, &format!("{dir}/{fullname}.db2.gz"), lig.db2.as_bytes())?;
    }
    let gz = builder.into_inner()?;
    gz.finish()
}

#[cfg(test)]
mod tests {
    use super::{archive_dir, full_name, write_tarball, ChargeLetter, ChargeOutOfRange, Ligand};
    use crate::name::MoleculeName;

    fn name(s: &str) -> MoleculeName {
        MoleculeName::new(s).expect("valid name")
    }

    /// `chr(78 + charge)`: N neutral, M -1, O +1 (the ZINC-22 3D build's member names), A and Z at
    /// the ends of the range; outside it there is no letter.
    #[test]
    fn charge_letters_cover_exactly_the_capital_letters() {
        for (charge, letter) in [(-13, 'A'), (-1, 'M'), (0, 'N'), (1, 'O'), (12, 'Z')] {
            assert_eq!(ChargeLetter::new(charge).map(ChargeLetter::as_char), Ok(letter), "charge {charge}");
        }
        for charge in [-14, 13, -31, i32::MIN, i32::MAX] {
            assert_eq!(ChargeLetter::new(charge), Err(ChargeOutOfRange(charge)));
        }
    }

    #[test]
    fn a_name_without_an_encoded_location_goes_to_the_top_level() {
        for n in ["mc0001_n5_t2_ha41", "CSLB0000477NCT", "ZINC000012345678", "ZINC1Z0000000001", "ZINCnot-an-id"] {
            assert_eq!(archive_dir(&name(n)), ".", "{n}");
        }
    }

    #[test]
    fn a_zinc22_id_goes_to_its_tranche_directory() {
        // fixture record: tranche H01M000, sub-id 1 (validation/fixtures/zinc22/tranche_dirs.tsv)
        assert_eq!(archive_dir(&name("ZINC150000000001")), "H01/H01M000/00/01");
    }

    #[test]
    fn member_names_follow_the_full_name() {
        assert_eq!(full_name(&name("benzoic"), 0, ChargeLetter::new(0).unwrap()), "benzoic.0.N");
        assert_eq!(full_name(&name("x"), 2, ChargeLetter::new(-1).unwrap()), "x.2.M");
    }

    fn members(ligands: &[Ligand]) -> Vec<String> {
        let bytes = write_tarball(ligands).expect("tarball");
        let mut ar = tar::Archive::new(flate2::read::GzDecoder::new(&bytes[..]));
        ar.entries().unwrap().map(|e| e.unwrap().path().unwrap().to_string_lossy().into_owned()).collect()
    }

    fn ligand(n: &str) -> Ligand {
        Ligand { name: name(n), prot_id: 0, charge: ChargeLetter::new(0).unwrap(), mol2: "m".into(), db2: "d".into() }
    }

    #[test]
    fn the_archive_files_each_ligand_by_its_name() {
        let got = members(&[ligand("mc0001"), ligand("ZINC150000000001")]);
        assert_eq!(
            got,
            [
                "mc0001.0.N.mol2", // `./` normalizes away (POSIX), as the driver's `./mc0001...` does
                "mc0001.0.N.db2.gz",
                "H01/H01M000/00/01/ZINC150000000001.0.N.mol2",
                "H01/H01M000/00/01/ZINC150000000001.0.N.db2.gz",
            ]
        );
    }

    proptest::proptest! {
        /// Any string either is rejected as a name or yields an archive whose members are relative,
        /// stay inside the archive (no `..`, no root), and are named after it — never a panic or a
        /// refused member. (A `./` prefix is inside the archive: tar keeps it on names past 100 bytes,
        /// and the driver writes every top-level member that way.)
        #[test]
        fn any_accepted_name_packs_safely(s in "\\PC{0,40}|ZINC[0-9a-zA-Z]{12}|[./]{1,4}") {
            use std::path::Component;
            if let Ok(n) = MoleculeName::new(&s) {
                for m in members(&[ligand(n.as_str())]) {
                    let p = std::path::Path::new(&m);
                    proptest::prop_assert!(p.is_relative(), "{m}");
                    proptest::prop_assert!(p.components().all(|c| matches!(c, Component::Normal(_) | Component::CurDir)), "{m}");
                    proptest::prop_assert!(m.ends_with(&format!("{s}.0.N.mol2")) || m.ends_with(&format!("{s}.0.N.db2.gz")), "{m}");
                }
            }
        }
    }
}
