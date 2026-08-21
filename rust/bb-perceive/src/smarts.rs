//! SMARTS parsing for the experimental-torsion library.
//!
//! The grammar lives in `src/grammar/smarts.pest` and is the specification we claim parity with:
//! it covers exactly the primitives RDKit's torsion library uses, so an unsupported construct is a
//! parse error rather than a silent misreading.
//!
//! Predicates are evaluated later against our own perception — rings, hybridization, hydrogen
//! counts — all of which reproduce RDKit exactly, rather than against a second implementation's
//! interpretation of them.

use pest::Parser;
use pest_derive::Parser;

#[derive(Parser)]
#[grammar = "grammar/smarts.pest"]
struct SmartsParser;

#[derive(Debug)]
pub struct ParseError(pub String);

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "SMARTS parse error: {}", self.0)
    }
}

impl std::error::Error for ParseError {}

use crate::smarts_ast::{AtomPred, BondPred, Pattern};
use pest::iterators::Pair;
use std::collections::HashMap;

/// Parse a SMARTS into a [`Pattern`].
pub fn parse(smarts: &str) -> Result<Pattern, ParseError> {
    let mut top = SmartsParser::parse(Rule::pattern, smarts)
        .map_err(|e| ParseError(e.to_string()))?;
    let pattern = top.next().ok_or_else(|| ParseError("empty".into()))?;
    let mut b = Builder::default();
    b.walk_top(pattern)?;
    b.close_rings();
    Ok(Pattern { atoms: b.atoms, maps: b.maps, bonds: b.bonds })
}

#[derive(Default)]
struct Builder {
    atoms: Vec<AtomPred>,
    maps: Vec<Option<u8>>,
    bonds: Vec<(usize, usize, BondPred)>,
    open_rings: HashMap<String, (usize, BondPred)>,
    /// closures held back until the whole pattern is parsed, as `(digit, opener, closer, bond)`
    closures: Vec<(u32, usize, usize, BondPred)>,
}

impl Builder {
    /// Append the ring-closure bonds, ascending by closure digit.
    ///
    /// `SmilesParseOps::CloseMolRings` runs after the whole string is parsed and walks the atom
    /// bookmarks, which are keyed by the closure digit, so every closure bond follows every chain
    /// bond and closures are ordered by digit rather than by position. Bond order decides which
    /// neighbour a match traversal reaches first, so this sequence is load-bearing.
    fn close_rings(&mut self) {
        self.closures.sort_by_key(|&(digit, ..)| digit);
        for (_, opener, closer, bond) in std::mem::take(&mut self.closures) {
            self.bonds.push((opener, closer, bond));
        }
    }

    fn walk_top(&mut self, pattern: Pair<Rule>) -> Result<(), ParseError> {
        let mut prev: Option<usize> = None;
        for p in pattern.into_inner() {
            match p.as_rule() {
                Rule::atom => prev = Some(self.push_atom(p)?),
                Rule::chain => self.walk_chain(p, &mut prev)?,
                Rule::EOI => {}
                r => return Err(ParseError(format!("unexpected {r:?}"))),
            }
        }
        Ok(())
    }

    fn walk_chain(&mut self, chain: Pair<Rule>, prev: &mut Option<usize>) -> Result<(), ParseError> {
        let mut pending = BondPred::Default;
        let mut saw_bond = false;
        for p in chain.into_inner() {
            match p.as_rule() {
                Rule::bond => {
                    pending = lower_bond(p)?;
                    saw_bond = true;
                }
                Rule::atom => {
                    let idx = self.push_atom(p)?;
                    if let Some(a) = *prev {
                        self.bonds.push((a, idx, pending.clone()));
                    }
                    *prev = Some(idx);
                    pending = BondPred::Default;
                    saw_bond = false;
                }
                Rule::branch => {
                    let anchor = *prev;
                    let mut inner_prev = anchor;
                    let mut first_bond = BondPred::Default;
                    let mut seen_first = false;
                    for q in p.into_inner() {
                        match q.as_rule() {
                            Rule::bond if !seen_first => first_bond = lower_bond(q)?,
                            Rule::atom if !seen_first => {
                                let idx = self.push_atom(q)?;
                                if let Some(a) = anchor {
                                    self.bonds.push((a, idx, first_bond.clone()));
                                }
                                inner_prev = Some(idx);
                                seen_first = true;
                            }
                            Rule::chain => self.walk_chain(q, &mut inner_prev)?,
                            _ => {}
                        }
                    }
                    // a branch does not advance the main chain
                }
                Rule::ring_closure => {
                    let text = p.as_str();
                    let digits: String = text.chars().filter(|c| c.is_ascii_digit()).collect();
                    let bond = if saw_bond { pending.clone() } else { BondPred::Default };
                    if let Some((other, ob)) = self.open_rings.remove(&digits) {
                        let use_bond = if matches!(bond, BondPred::Default) { ob } else { bond };
                        if let Some(a) = *prev {
                            let digit = digits.parse::<u32>().unwrap_or(u32::MAX);
                            self.closures.push((digit, other, a, use_bond));
                        }
                    } else if let Some(a) = *prev {
                        self.open_rings.insert(digits, (a, bond));
                    }
                    pending = BondPred::Default;
                    saw_bond = false;
                }
                _ => {}
            }
        }
        Ok(())
    }

    fn push_atom(&mut self, atom: Pair<Rule>) -> Result<usize, ParseError> {
        let mut pred = AtomPred::Any;
        let mut map = None;
        for p in atom.into_inner() {
            match p.as_rule() {
                Rule::organic_atom => pred = lower_symbol(p.as_str()),
                Rule::bracket_atom => {
                    for q in p.into_inner() {
                        match q.as_rule() {
                            // `[H]` on its own is a hydrogen atom, the same as `[#1]`; `H`
                            // alongside other terms (`[CH]`, `[C;H]`) is "exactly one hydrogen".
                            // Verified against RDKit: [H] and [#1] select identical atom sets.
                            Rule::atom_expr if q.as_str().trim() == "H" => {
                                pred = AtomPred::AtomicNum(1);
                            }
                            Rule::atom_expr => pred = lower_atom_expr(q)?,
                            Rule::atom_map => {
                                map = q.as_str().trim_start_matches(':').parse().ok();
                            }
                            _ => {}
                        }
                    }
                }
                _ => {}
            }
        }
        self.atoms.push(pred);
        self.maps.push(map);
        Ok(self.atoms.len() - 1)
    }
}

fn lower_symbol(s: &str) -> AtomPred {
    match s {
        "*" => AtomPred::Any,
        "a" => AtomPred::AnyAromatic,
        "A" => AtomPred::AnyAliphatic,
        _ => {
            let aromatic = s.chars().next().is_some_and(|c| c.is_ascii_lowercase());
            let cap = {
                let mut c = s.chars();
                match c.next() {
                    Some(f) => f.to_uppercase().collect::<String>() + c.as_str(),
                    None => String::new(),
                }
            };
            let z = mendeleev::ALL_ELEMENTS
                .iter()
                .find(|e| e.symbol() == cap)
                .map_or(0, |e| e.atomic_number() as u8);
            AtomPred::Symbol { z, aromatic: Some(aromatic) }
        }
    }
}

fn lower_atom_expr(p: Pair<Rule>) -> Result<AtomPred, ParseError> {
    let mut terms = Vec::new();
    for q in p.into_inner() {
        if q.as_rule() == Rule::atom_or {
            terms.push(lower_atom_or(q)?);
        }
    }
    Ok(if terms.len() == 1 { terms.pop().expect("one") } else { AtomPred::And(terms) })
}

fn lower_atom_or(p: Pair<Rule>) -> Result<AtomPred, ParseError> {
    let mut terms = Vec::new();
    for q in p.into_inner() {
        if q.as_rule() == Rule::atom_and {
            terms.push(lower_atom_and(q)?);
        }
    }
    Ok(if terms.len() == 1 { terms.pop().expect("one") } else { AtomPred::Or(terms) })
}

fn lower_atom_and(p: Pair<Rule>) -> Result<AtomPred, ParseError> {
    let mut terms = Vec::new();
    for q in p.into_inner() {
        if q.as_rule() == Rule::atom_not {
            terms.push(lower_atom_not(q)?);
        }
    }
    Ok(if terms.len() == 1 { terms.pop().expect("one") } else { AtomPred::And(terms) })
}

fn lower_atom_not(p: Pair<Rule>) -> Result<AtomPred, ParseError> {
    let negated = p.as_str().starts_with('!');
    let inner = p
        .into_inner()
        .find(|q| q.as_rule() == Rule::atom_prim)
        .ok_or_else(|| ParseError("missing primitive".into()))?;
    let prim = lower_atom_prim(inner)?;
    Ok(if negated { AtomPred::Not(Box::new(prim)) } else { prim })
}

fn digits_of(s: &str) -> Option<u8> {
    let d: String = s.chars().filter(|c| c.is_ascii_digit()).collect();
    if d.is_empty() { None } else { d.parse().ok() }
}

fn lower_atom_prim(p: Pair<Rule>) -> Result<AtomPred, ParseError> {
    let text = p.as_str().to_string();
    let inner = p.into_inner().next();
    let Some(q) = inner else { return Ok(lower_symbol(&text)) };
    Ok(match q.as_rule() {
        Rule::recursive => {
            // strip exactly one "$(" and its matching ")" — trimming all trailing parens would
            // unbalance a sub-pattern that legitimately ends in one, e.g. $([CX3](=[!O]))
            let raw = q.as_str();
            let inner = raw
                .strip_prefix("$(")
                .and_then(|r| r.strip_suffix(')'))
                .ok_or_else(|| ParseError(format!("malformed recursive query {raw}")))?;
            AtomPred::Recursive(Box::new(parse(inner)?))
        }
        Rule::atomic_num => AtomPred::AtomicNum(digits_of(q.as_str()).unwrap_or(0)),
        Rule::total_h => AtomPred::TotalH(digits_of(q.as_str())),
        Rule::connectivity => AtomPred::Connectivity(digits_of(q.as_str()).unwrap_or(0)),
        Rule::degree => AtomPred::Degree(digits_of(q.as_str()).unwrap_or(0)),
        Rule::valence => AtomPred::Valence(digits_of(q.as_str()).unwrap_or(0)),
        Rule::ring_count => AtomPred::RingCount(digits_of(q.as_str())),
        // bare `r` is ring membership, the same question `R` asks
        Rule::ring_size => match digits_of(q.as_str()) {
            Some(n) => AtomPred::RingSize(n),
            None => AtomPred::RingCount(None),
        },
        Rule::ring_size_range => {
            let d: String = q.as_str().chars().skip(2).take_while(|c| c.is_ascii_digit()).collect();
            AtomPred::MinRingSize(d.parse().unwrap_or(0))
        }
        Rule::ring_bond_count => AtomPred::RingBondCount(digits_of(q.as_str()).unwrap_or(0)),
        Rule::hybridization => AtomPred::Hybridization(digits_of(q.as_str()).unwrap_or(0)),
        Rule::charge => {
            let n = digits_of(q.as_str()).unwrap_or(1) as i8;
            AtomPred::Charge(if q.as_str().starts_with('-') { -n } else { n })
        }
        Rule::symbol => lower_symbol(q.as_str()),
        r => return Err(ParseError(format!("unhandled primitive {r:?}"))),
    })
}

fn lower_bond(p: Pair<Rule>) -> Result<BondPred, ParseError> {
    let mut terms = Vec::new();
    for q in p.into_inner() {
        if q.as_rule() == Rule::bond_or {
            terms.push(lower_bond_or(q)?);
        }
    }
    Ok(if terms.len() == 1 { terms.pop().expect("one") } else { BondPred::And(terms) })
}

fn lower_bond_or(p: Pair<Rule>) -> Result<BondPred, ParseError> {
    let mut terms = Vec::new();
    for q in p.into_inner() {
        if q.as_rule() == Rule::bond_and {
            terms.push(lower_bond_and(q)?);
        }
    }
    Ok(if terms.len() == 1 { terms.pop().expect("one") } else { BondPred::Or(terms) })
}

fn lower_bond_and(p: Pair<Rule>) -> Result<BondPred, ParseError> {
    let mut terms = Vec::new();
    for q in p.into_inner() {
        if q.as_rule() == Rule::bond_not {
            terms.push(lower_bond_not(q)?);
        }
    }
    Ok(if terms.len() == 1 { terms.pop().expect("one") } else { BondPred::And(terms) })
}

fn lower_bond_not(p: Pair<Rule>) -> Result<BondPred, ParseError> {
    let negated = p.as_str().starts_with('!');
    let prim = p
        .into_inner()
        .find(|q| q.as_rule() == Rule::bond_prim)
        .map(|q| match q.as_str() {
            "-" => BondPred::Single,
            "=" => BondPred::Double,
            "#" => BondPred::Triple,
            ":" => BondPred::Aromatic,
            "~" => BondPred::Any,
            "@" => BondPred::Ring,
            "/" => BondPred::Up,
            "\\" => BondPred::Down,
            _ => BondPred::Default,
        })
        .unwrap_or(BondPred::Default);
    Ok(if negated { BondPred::Not(Box::new(prim)) } else { prim })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn whole_torsion_library_parses() {
        let pats = crate::torsion_lib::patterns();
        let mut ok = 0;
        let mut failed: Vec<(String, String)> = Vec::new();
        for p in &pats {
            match parse(&p.smarts) {
                Ok(_) => ok += 1,
                Err(e) => failed.push((p.smarts.clone(), e.0.lines().next().unwrap_or("").into())),
            }
        }
        println!("  {ok}/{} torsion patterns parse", pats.len());
        for (s, e) in failed.iter().take(6) {
            println!("    FAIL {s}\n      {e}");
        }
        assert_eq!(failed.len(), 0, "{} patterns fail to parse", failed.len());
    }
}

#[cfg(test)]
mod lowering {
    use super::*;
    use crate::smarts_ast::{AtomPred, BondPred};

    #[test]
    fn quadruple_and_predicates() {
        let p = parse("[O:1]=[CX3;r{9-}:2](a)@;-[NX3H0;r:3][!#1:4]").expect("parses");
        // four mapped atoms plus the aromatic branch neighbour
        assert_eq!(p.atoms.len(), 5, "{:?}", p.atoms);
        assert_eq!(p.torsion_quad(), Some([0, 1, 3, 4]));
        // the ring-size range survives as a predicate rather than being lifted out
        let has_min_ring = format!("{:?}", p.atoms[1]).contains("MinRingSize(9)");
        assert!(has_min_ring, "{:?}", p.atoms[1]);
        // the bond conjunction is a real AND node
        let (_, _, b) = p.bonds.iter().find(|(a, b, _)| *a == 1 && *b == 3).expect("2-3 bond");
        assert!(
            matches!(b, BondPred::And(_)),
            "expected a conjunction, got {b:?}"
        );
    }

    #[test]
    fn recursive_subpattern_is_nested() {
        let p = parse("[$(C=O):1][O:2]").expect("parses");
        assert!(matches!(p.atoms[0], AtomPred::Recursive(_)), "{:?}", p.atoms[0]);
    }

    #[test]
    fn whole_library_lowers() {
        let pats = crate::torsion_lib::patterns();
        let mut quads = 0;
        let mut no_quad = Vec::new();
        let mut failed = Vec::new();
        for p in &pats {
            match parse(&p.smarts) {
                Ok(ast) => match ast.torsion_quad() {
                    Some(_) => quads += 1,
                    None => no_quad.push(p.smarts.clone()),
                },
                Err(e) => failed.push((p.smarts.clone(), e.0.lines().next().unwrap_or("").to_string())),
            }
        }
        for (s, e) in failed.iter().take(4) {
            println!("    LOWER FAIL {s}  -> {e}");
        }
        assert!(failed.is_empty(), "{} patterns fail to lower", failed.len());
        println!("  {quads}/{} patterns yield a 1-2-3-4 quadruple", pats.len());
        for s in no_quad.iter().take(4) {
            println!("    no quad: {s}");
        }
    }
}
