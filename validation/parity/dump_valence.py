"""Dump RDKit's per-element valence data from its own accessors."""
import json, sys
from rdkit.Chem import GetPeriodicTable
pt = GetPeriodicTable()
out = {}
for z in range(1, 119):
    try:
        out[z] = {
            "symbol": pt.GetElementSymbol(z),
            "default_valence": pt.GetDefaultValence(z),
            "valence_list": list(pt.GetValenceList(z)),
            "n_outer_elecs": pt.GetNOuterElecs(z),
            "rvdw": pt.GetRvdw(z),
        }
    except Exception as e:
        out[z] = {"error": str(e)}
json.dump(out, open(sys.argv[1], "w"), sort_keys=True)
print(f"  dumped {len(out)} elements")
for z in (1,5,6,7,8,9,14,15,16,17,35,53):
    r = out[z]
    print(f"    Z={z:<3} {r['symbol']:<3} default={r['default_valence']:<3} "
          f"list={r['valence_list']} outer={r['n_outer_elecs']}")
