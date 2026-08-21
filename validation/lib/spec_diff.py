"""Single source for the bridge-vs-native MoleculeSpec field canonicalisers.

Two gates diff the RDKit-bridge spec against the native spec: `validation/hunt/differential.py`
(the wired divergence hunt, hard/stereo/numeric split) and `validation/parity/gate_spec.py` (the
per-field per-molecule model). They apply DIFFERENT diff policies, but the tricky, drift-prone parts
— the stereo-double-bond canonicalisation ([[stereo-representation-equivalence]]) and the
torsion/improper set canonicalisation — must be IDENTICAL between them, or the equivalence proof
could differ by gate. Those canonicalisers live here; each gate keeps its own `compare()` policy on
top of them.
"""
import collections


def tup(x):
    """Atom-tuple key for a constraint record (dict with 'atoms') or a raw list."""
    return tuple(x["atoms"]) if isinstance(x, dict) else tuple(x)


def canon_set(items):
    """Order-insensitive canonical form of a list of lists/tuples."""
    return sorted(tuple(x) if isinstance(x, list) else x for x in items)


def num_diff(a, b, tol):
    """True if two float vectors differ by more than `tol` (absolute) + `tol`·relative."""
    if a is None or b is None or len(a) != len(b):
        return True
    return any(abs(x - y) > tol + 1e-4 * max(abs(x), abs(y)) for x, y in zip(a, b))


def canon_stereo(sdb, bonds):
    """Representation-invariant form of the stereo double bonds. A bond (a,b,c,d) has central pair
    (b,c) and reference substituents a (on b) and d (on c); RDKit and native may pick DIFFERENT
    same-end substituents, but switching a reference to the other substituent on that end flips the
    E/Z sign. Canonicalise to the lowest-indexed neighbour on each end with the sign adjusted → two
    encodings of the SAME geometry map to the SAME tuple. Returns {(b,c): (ref_b, ref_c, canon_sign)}.

    Neighbours come from `bonds` (always emitted by both sides). AddHs appends hydrogens after every
    heavy atom, so their indices are the highest — `min()` therefore never picks an H over a heavy
    substituent, matching the heavy-neighbour reference either gate would otherwise choose.
    """
    adj = collections.defaultdict(set)
    for i, j in bonds:
        adj[i].add(j)
        adj[j].add(i)
    out = {}
    for e in sdb:
        a, b, c, d = e["atoms"]
        s = e["sign"]
        rb = min(adj[b] - {c}) if adj[b] - {c} else a   # canonical ref on b's end
        rc = min(adj[c] - {b}) if adj[c] - {b} else d
        flips = (1 if a != rb else 0) + (1 if d != rc else 0)
        out[(b, c)] = (rb, rc, s * (-1) ** flips)
    return out


def canon_chiral(cs):
    """Canonical set form of chiral_sets / tetrahedral_centers (center + atoms + rounded volumes)."""
    return sorted(
        (c["center"], tuple(c["atoms"]), round(c["vol_lo"], 6), round(c["vol_hi"], 6),
         c["fused_small_rings"])
        for c in cs
    )


def canon_tors(ts):
    """Canonical set form of the experimental torsions.

    A torsion (i,j,k,l) and its reversal (l,k,j,i) are the same physical term: the signed dihedral is
    invariant under reversal (b1,b2,b3 -> -b3,-b2,-b1 gives n1,n2 -> -n2,-n1, leaving
    atan2(dot(cross(n1,n2),b2), dot(n1,n2)) unchanged), and V/signs index the multiplicity m, not the
    atoms, so E = Σ V[m](1+s[m]cos(mφ)) and its gradient are identical. RDKit and the native ETKDG
    enumeration pick opposite ring traversal directions, so canonicalise each quadruple to the smaller
    of (atoms, reversed atoms).
    """
    def key_atoms(a):
        a = tuple(a)
        return min(a, a[::-1])
    return sorted(
        (key_atoms(t["atoms"]), tuple(round(v, 4) for v in t["v"]), tuple(t["signs"]))
        for t in ts
    )


def canon_imp(imps):
    """Canonical set form of the impropers (atoms + rounded coefficients + force constant)."""
    return sorted(
        (tuple(i["atoms"]), round(i["c0"], 4), round(i["c1"], 4), round(i["c2"], 4),
         round(i["fc"], 4))
        for i in imps
    )
