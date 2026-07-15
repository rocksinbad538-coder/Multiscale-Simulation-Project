from pathlib import Path
import math

GRO = Path("runs/phase1A/accepted/hybrid_hbnBonded_kang2000_improperGeo100_hydrated_baseline_10ps100K_postmin/phase1A_hydrated_baseline.gro")
OUTDIR = Path("runs/phase1A/day013_orca_hydrated_dimer_pilots/PYR4_PYR5_water0p50_serial_tight")
OUTDIR.mkdir(parents=True, exist_ok=True)

PAIR = [4, 5]
CUTOFF_NM = 0.50

lines = GRO.read_text().splitlines()
natoms = int(lines[1].strip())

atoms = []
for line in lines[2:2+natoms]:
    name = line[10:15].strip()
    resname = line[5:10].strip()
    elem = ''.join([c for c in name if c.isalpha()])[0]
    atoms.append({
        "resid": int(line[0:5]),
        "resname": resname,
        "name": name,
        "elem": elem,
        "x_nm": float(line[20:28]),
        "y_nm": float(line[28:36]),
        "z_nm": float(line[36:44]),
    })

pyr = [a for a in atoms if a["resname"] == "PYR" and a["resid"] in PAIR]
waters = {}
for a in atoms:
    if a["resname"] == "SOL":
        waters.setdefault(a["resid"], []).append(a)

def dist(a,b):
    return math.sqrt((a["x_nm"]-b["x_nm"])**2 + (a["y_nm"]-b["y_nm"])**2 + (a["z_nm"]-b["z_nm"])**2)

selected = list(pyr)
selected_waters = []

for wid, wat in waters.items():
    ow = [a for a in wat if a["name"] == "OW"]
    if not ow:
        continue
    dmin = min(dist(ow[0], p) for p in pyr)
    if dmin < CUTOFF_NM:
        selected_waters.append(wid)
        selected.extend(wat)

coords = []
for a in selected:
    if a["elem"] in {"M", "X"}:
        continue
    coords.append(f"{a['elem']:2s} {a['x_nm']*10:12.6f} {a['y_nm']*10:12.6f} {a['z_nm']*10:12.6f}")

xyz = OUTDIR / "PYR4_PYR5_water0p50.xyz"
xyz.write_text(f"{len(coords)}\nPYR4-PYR5 + waters within 0.50 nm; water_resids={selected_waters}\n" + "\n".join(coords) + "\n")

inp = OUTDIR / "PYR4_PYR5_water0p50_serial_tight.inp"
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
