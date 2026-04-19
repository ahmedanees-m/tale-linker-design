import sys, numpy as np
sys.path.insert(0, 'src')
from Bio.PDB import PDBParser

p = PDBParser(QUIET=True)
s = p.get_structure('3V6T', 'data/pdb_cleaned/3V6T.pdb')
model = list(s.get_models())[0]
c_coords = np.array([29.606, -40.355, 26.639])
_dna_set = {"DA","DT","DC","DG","A","T","C","G"}

# Chain H in full
chain = model['H']
pts = []
for i, r in enumerate(chain.get_residues()):
    if r.get_resname().strip() in _dna_set and "P" in r:
        xyz = np.array(r["P"].get_vector().get_array())
        pts.append((i, r.get_id()[1], xyz))

print(f"Chain H: {len(pts)} phosphates with P")
for raw_idx, seqid, xyz in pts:
    d = np.linalg.norm(xyz - c_coords)
    print(f"  raw_idx={raw_idx} seqid={seqid} dist={d:.1f}A")

anchor = min(range(len(pts)), key=lambda i: np.linalg.norm(pts[i][2]-c_coords))
anchor_raw = pts[anchor][0]
print(f"\nAnchor: raw_idx={anchor_raw} (bp_offset=0)")
for raw_idx, seqid, xyz in pts:
    offset = raw_idx - anchor_raw
    print(f"  offset={offset}  seqid={seqid}")
