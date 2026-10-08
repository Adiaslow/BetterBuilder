//! A molecule's name from the input, checked once where it is read. The name becomes a file name:
//! the molecule's work directory (`<name>.<prot_id>`) and its members in the output archive
//! (`<name>.<prot_id>.<C>.<ext>`). A name that cannot be a single file-name component would put files
//! outside the work directory or make the archive refuse the member, so it is rejected here and that
//! molecule is reported and skipped; everything downstream only ever holds a [`MoleculeName`].

use std::fmt;

/// Why an input name cannot name the molecule's files.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum NameError {
    #[error("the name is empty")]
    Empty,
    #[error("name `{0}` is `.` or `..`, which name directories, not files")]
    DotComponent(String),
    #[error("name `{0}` contains `/`, which would make it a path")]
    Separator(String),
    #[error("name {0:?} contains a NUL byte, which no file name can hold")]
    Nul(String),
}

/// A molecule name that is usable as a single file-name component: not empty, not `.` or `..`, and
/// free of `/` and NUL (the only bytes a POSIX file name cannot contain).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct MoleculeName(String);

impl MoleculeName {
    pub fn new(name: &str) -> Result<Self, NameError> {
        if name.is_empty() {
            Err(NameError::Empty)
        } else if name == "." || name == ".." {
            Err(NameError::DotComponent(name.to_string()))
        } else if name.contains('/') {
            Err(NameError::Separator(name.to_string()))
        } else if name.contains('\0') {
            Err(NameError::Nul(name.to_string()))
        } else {
            Ok(MoleculeName(name.to_string()))
        }
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for MoleculeName {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::{MoleculeName, NameError};

    #[test]
    fn names_that_are_file_names_are_accepted_unchanged() {
        for name in ["mc0001_n5_t2_ha41", "CSLB0000477NCT", "ZINChg0000000001", "x.y", "a b", "αβγ", "-", "..."] {
            assert_eq!(MoleculeName::new(name).map(|n| n.as_str().to_string()), Ok(name.to_string()));
        }
    }

    #[test]
    fn names_that_cannot_name_a_file_are_rejected() {
        assert_eq!(MoleculeName::new(""), Err(NameError::Empty));
        assert_eq!(MoleculeName::new("."), Err(NameError::DotComponent(".".into())));
        assert_eq!(MoleculeName::new(".."), Err(NameError::DotComponent("..".into())));
        for name in ["a/b", "/abs", "../up", "trailing/"] {
            assert_eq!(MoleculeName::new(name), Err(NameError::Separator(name.into())));
        }
        assert_eq!(MoleculeName::new("a\0b"), Err(NameError::Nul("a\0b".into())));
    }
}
