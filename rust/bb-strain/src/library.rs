//! The UCSF Torsion Library (`TL_2.1_VERSION_6.xml`), parsed into ordered torsion rules.
//!
//! The file is RDKit's/Shoichet's torsion-strain library, vendored verbatim under `rdkit-patch/`. Rule
//! order encodes specificity, which the strain lookup uses to keep the most-specific match per torsion:
//! `TL_Functions.TL_Lookup.__init__` walks every non-`GG` hierarchy class in document order, then the
//! general `GG` class last, assigning each rule an increasing index. We reproduce that order exactly.
//!
//! Each rule is `exact` (a 36-bin energy histogram, interpolated) or `approximate` (a list of angle
//! peaks with a `β₁·δ² + β₂·δ⁴` penalty), matched by a SMARTS with `:1`–`:4` torsion atom maps.

use std::sync::OnceLock;

const LIBRARY_XML: &str = include_str!("../../../rdkit-patch/torsion_library_TL2.1_v6.xml");

/// An angle peak for an `approximate` rule.
#[derive(Debug, Clone)]
pub struct Angle {
    pub theta_0: f64,
    pub tolerance2: f64,
    pub beta_1: f64,
    pub beta_2: f64,
}

/// A rule's energy model.
#[derive(Debug, Clone)]
pub enum Method {
    /// 36 bins (10° each from −180°), each an (energy, lower, upper) triple.
    Exact { energy: Vec<f64>, lower: Vec<f64>, upper: Vec<f64> },
    /// angle peaks for the not-as-small-angle approximation.
    Approximate { angles: Vec<Angle> },
}

/// One torsion rule: its SMARTS, energy model, specificity index (lower = more specific), and whether
/// it comes from the general `GG` class (vs a specific class).
#[derive(Debug, Clone)]
pub struct Rule {
    pub smarts: String,
    pub method: Method,
    pub index: usize,
    pub is_general: bool,
}

fn parse_rule(node: roxmltree::Node, index: usize, is_general: bool) -> Option<Rule> {
    // some library SMARTS carry trailing whitespace; RDKit's MolFromSmarts ignores it, so trim to
    // match (an untrimmed pattern is rejected by the parser and the rule would be silently dropped)
    let smarts = node.attribute("smarts")?.trim().to_string();
    let method = match node.attribute("method")? {
        "exact" => {
            let hist = node.children().find(|c| c.has_tag_name("histogram_converted"))?;
            let (mut energy, mut lower, mut upper) = (Vec::new(), Vec::new(), Vec::new());
            for bin in hist.children().filter(|c| c.has_tag_name("bin")) {
                energy.push(bin.attribute("energy")?.parse().ok()?);
                lower.push(bin.attribute("lower")?.parse().ok()?);
                upper.push(bin.attribute("upper")?.parse().ok()?);
            }
            Method::Exact { energy, lower, upper }
        }
        _ => {
            let list = node.children().find(|c| c.has_tag_name("angleList"))?;
            let mut angles = Vec::new();
            for a in list.children().filter(|c| c.has_tag_name("angle")) {
                angles.push(Angle {
                    theta_0: a.attribute("theta_0")?.parse().ok()?,
                    tolerance2: a.attribute("tolerance2")?.parse().ok()?,
                    beta_1: a.attribute("beta_1")?.parse().ok()?,
                    beta_2: a.attribute("beta_2")?.parse().ok()?,
                });
            }
            Method::Approximate { angles }
        }
    };
    Some(Rule { smarts, method, index, is_general })
}

/// The parsed library, in the lookup's specificity order (non-`GG` classes first, `GG` last).
pub fn library() -> &'static Vec<Rule> {
    static LIB: OnceLock<Vec<Rule>> = OnceLock::new();
    LIB.get_or_init(|| {
        let doc = roxmltree::Document::parse(LIBRARY_XML).expect("parse torsion library");
        let root = doc.root_element();
        let classes: Vec<roxmltree::Node> = root.children().filter(|c| c.has_tag_name("hierarchyClass")).collect();
        let mut rules = Vec::new();
        let mut idx = 0;
        // specific classes (document order), then the general GG class
        for pass_general in [false, true] {
            for class in &classes {
                let is_gg = class.attribute("name") == Some("GG");
                if is_gg != pass_general {
                    continue;
                }
                for tr in class.descendants().filter(|c| c.has_tag_name("torsionRule")) {
                    if let Some(rule) = parse_rule(tr, idx, pass_general) {
                        rules.push(rule);
                    }
                    idx += 1; // index advances per rule slot, matching the Python counter
                }
            }
        }
        rules
    })
}
