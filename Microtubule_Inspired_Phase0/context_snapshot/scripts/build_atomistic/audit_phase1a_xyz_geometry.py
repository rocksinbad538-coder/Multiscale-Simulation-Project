from pathlib import Path
import math, csv, sys
from collections import Counter

if len(sys.argv) != 3:
    print("Usage: python audit_phase1a_xyz_geometry.py input.xyz output.csv")
    sys.exit(1)

inp = Path(sys.argv[1])
out = Path(sys.argv[2])

lines = inp.read_text().splitlines()
natoms = int(lines[0].strip())

atoms = []
for i, line in enumerate(lines[2:2+natoms], start=1):
    p = line.split()
    atoms.append((i, p[0], float(p[1]), float(p[2]), float(p[3])))

def d(a,b):
    return math.sqrt((a[2]-b[2])**2 + (a[3]-b[3])**2 + (a[4]-b[4])**2)

counts = Counter(a[1] for a in atoms)

bn = []
bh = []
nh = []
overlaps = []

for i,a in enumerate(atoms):
    for b in atoms[i+1:]:
        dist = d(a,b)
        pair = ''.join(sorted([a[1], b[1]]))
        if dist < 0.50:
            overlaps.append((a[0],a[1],b[0],b[1],dist))
        if set([a[1],b[1]]) == set(["B","N"]) and dist < 2.0:
            bn.append(dist)
        if set([a[1],b[1]]) == set(["B","H"]) and dist < 1.6:
            bh.append(dist)
        if set([a[1],b[1]]) == set(["N","H"]) and dist < 1.4:
            nh.append(dist)

def stats(vals):
    if not vals:
        return ("NA","NA","NA","0")
    return (f"{min(vals):.5f}", f"{sum(vals)/len(vals):.5f}", f"{max(vals):.5f}", str(len(vals)))

rows = [
    ("source", str(inp)),
    ("natoms", str(natoms)),
    ("element_counts", dict(counts)),
    ("B-N_min_mean_max_count", stats(bn)),
    ("B-H_min_mean_max_count", stats(bh)),
    ("N-H_min_mean_max_count", stats(nh)),
    ("overlaps_lt_0p50A_count", str(len(overlaps))),
]

out.parent.mkdir(parents=True, exist_ok=True)
with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric","value"])
    w.writerows(rows)

print(out.read_text())
if overlaps[:20]:
    print("First overlaps <0.50 A:")
    for row in overlaps[:20]:
        print(row)
