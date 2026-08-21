//! The experimental-torsion library, read from RDKit's own data files.
//!
//! `rdkit-patch/torsionPreferences_v2.in` and `_macrocycles.in` are RDKit's files verbatim — C++
//! string-literal arrays — including Divya's amide reweighting (force constants 80.0 and 130.9 in
//! place of 8.0 and 13.9). Each entry is a SMARTS whose atom maps 1..4 mark the torsion quadruple,
//! followed by six (sign, V) pairs.
//!
//! Read rather than transcribed, so the patch stays visible in the file it was made to and a
//! regenerated RDKit checkout can be diffed against it.

/// One library entry.
#[derive(Debug, Clone, PartialEq)]
pub struct TorsionPattern {
    pub smarts: String,
    pub signs: [i8; 6],
    pub v: [f64; 6],
    /// Whether the pattern constrains ring size (`r{N-}`), which the SMARTS engine cannot express.
    pub ring_min: Option<u8>,
}

const V2: &str = include_str!("../../../rdkit-patch/torsionPreferences_v2.in");
const MACRO: &str = include_str!("../../../rdkit-patch/torsionPreferences_macrocycles.in");

/// Extract the quoted payloads from a C++ string-literal array, ignoring comments.
fn literals(src: &str) -> Vec<String> {
    let mut out = Vec::new();
    for line in src.lines() {
        let trimmed = line.trim_start();
        if !trimmed.starts_with('"') {
            continue;
        }
        // take the first quoted run only; trailing // comments may contain quotes of their own
        let body = &trimmed[1..];
        let Some(end) = body.find('"') else { continue };
        out.push(body[..end].trim_end_matches("\\n").to_string());
    }
    out
}

/// `r{9-}` and friends: the minimum ring size a pattern demands, if any.
fn ring_min_of(smarts: &str) -> Option<u8> {
    let idx = smarts.find("r{")?;
    let rest = &smarts[idx + 2..];
    let end = rest.find('}')?;
    let spec = &rest[..end];
    spec.trim_end_matches('-').parse().ok()
}

fn parse_entry(line: &str) -> Option<TorsionPattern> {
    let mut it = line.split_whitespace();
    let smarts = it.next()?.to_string();
    let mut signs = [0i8; 6];
    let mut v = [0f64; 6];
    for k in 0..6 {
        signs[k] = it.next()?.parse().ok()?;
        v[k] = it.next()?.parse().ok()?;
    }
    let ring_min = ring_min_of(&smarts);
    Some(TorsionPattern {
        smarts,
        signs,
        v,
        ring_min,
    })
}

/// Every pattern, in file order: the general library then the macrocycle library, matching the
/// order RDKit concatenates them.
pub fn patterns() -> Vec<TorsionPattern> {
    let mut out = Vec::new();
    for src in [V2, MACRO] {
        for line in literals(src) {
            if let Some(p) = parse_entry(&line) {
                out.push(p);
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn library_reads() {
        let p = patterns();
        assert!(p.len() > 700, "expected the full library, got {}", p.len());
        // every entry carries six (sign, V) pairs
        for e in &p {
            assert!(e.signs.iter().all(|s| (-1..=1).contains(s)), "{}", e.smarts);
        }
    }

    #[test]
    fn divyas_patch_is_present() {
        let p = patterns();
        let patched: Vec<&TorsionPattern> = p
            .iter()
            .filter(|e| e.v.iter().any(|&x| (x - 80.0).abs() < 1e-9 || (x - 130.9).abs() < 1e-9))
            .collect();
        assert!(
            !patched.is_empty(),
            "the amide reweighting (80.0 / 130.9) is missing from the vendored library"
        );
        // and it is on macrocyclic amide patterns
        assert!(patched.iter().all(|e| e.ring_min == Some(9)), "patched entries should be r{{9-}}");
    }

    #[test]
    fn ring_ranges_are_recognised() {
        assert_eq!(ring_min_of("[C;r{9-}:2]"), Some(9));
        assert_eq!(ring_min_of("[CX3:1]"), None);
        let p = patterns();
        let with_ring = p.iter().filter(|e| e.ring_min.is_some()).count();
        assert!(with_ring > 300, "expected the macrocycle library, got {with_ring}");
    }
}
