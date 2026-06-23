from pathlib import Path
import argparse
import math

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--pair", nargs=2, type=int, required=True, help="PYR residue ids, e.g. 3 4")
    p.add_argument("--cutoff-nm", type=float, default=0.50)
    p.add_argument("--outdir", required=True)
    p.add_argument("--tag", required=True)
    p.add_argument("--gro", default="runs/phase1A/accepted/hybrid_hbnBonded_kang2000_improperGeo100_hydrated_baseline_10ps100K_postmin/phase1A_hydrated_baseline.gro")
    return p.parse_args()

def dist(a, b):
    return math.sqrt(
        (a["x_nm"] - b["x_nm"])**2 +
        (a["y_nm"] - b["y_nm"])**2 +
        (a["z_nm"] - b["z_nm"])**2
    )

args = parse_args()

gro = Path(args.gro)
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

lines = gro.read_text().splitlines()
natoms = int(lines[1].strip())

atoms = []
for line in lines[2:2 + natoms]:
    name = line[10:15].strip()
    resname = line[5:10].strip()
    elem_letters = ''.join([c for c in name if c.isalpha()])
    if not elem_letters:
        continue
    elem = elem_letters[0]
    atoms.append({
        "resid": int(line[0:5]),
        "resname": resname,
        "name": name,
        "elem": elem,
        "x_nm": float(line[20:28]),
        "y_nm": float(line[28:36]),
        "z_nm": float(line[36:44]),
    })

pair = set(args.pair)
pyr = [a for a in atoms if a["resname"] == "PYR" and a["resid"] in pair]

if not pyr:
    raise RuntimeError(f"No PYR atoms found for pair {args.pair}")

waters = {}
for a in atoms:
    if a["resname"] == "SOL":
        waters.setdefault(a["resid"], []).append(a)

selected = list(pyr)
selected_waters = []

for wid, wat in waters.items():
    ow = [a for a in wat if a["name"] == "OW"]
    if not ow:
        continue
    dmin = min(dist(ow[0], p) for p in pyr)
    if dmin < args.cutoff_nm:
        selected_waters.append(wid)
        selected.extend(wat)

coords = []
for a in selected:
    if a["elem"] in {"M", "X"}:
        continue
    coords.append(f"{a['elem']:2s} {a['x_nm']*10:12.6f} {a['y_nm']*10:12.6f} {a['z_nm']*10:12.6f}")

xyz = outdir / f"{args.tag}.xyz"
inp = outdir / f"{args.tag}_serial_tight.inp"

xyz.write_text(
    f"{len(coords)}\nPYR{args.pair[0]}-PYR{args.pair[1]} + waters within {args.cutoff_nm:.2f} nm; water_resids={selected_waters}\n"
    + "\n".join(coords) + "\n"
)

inp.write_text("""! PBE0 def2-SVP TightSCF SlowConv

%tddft
  nroots 20
end

%scf
  MaxIter 500
  ConvForced true
end

* xyz 0 1
""" + "\n".join(coords) + "\n*\n")

print(xyz, "qm_atoms", len(coords), "waters", len(selected_waters))
print(inp)
