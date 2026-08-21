//! Minimal 3-vector primitives shared across the embed's force field and acceptance checks. Hoisted
//! here so `cross`/`dot`/`sub`/`scale`/`norm` have one definition (they were re-implemented per module,
//! with a dim-aware vs dim-unaware `sub` that was easy to confuse). These operate on plain `[f64; 3]`
//! arrays; extracting a point from a strided coordinate slice stays with each caller (the stride is
//! theirs to know), but those extractors build on [`sub`] here.

/// Cross product `a × b`.
#[inline]
pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Dot product `a · b`.
#[inline]
pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Component-wise difference `a − b`.
#[inline]
pub fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scalar multiple `a · s`.
#[inline]
pub fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Euclidean norm `‖a‖`.
#[inline]
pub fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
