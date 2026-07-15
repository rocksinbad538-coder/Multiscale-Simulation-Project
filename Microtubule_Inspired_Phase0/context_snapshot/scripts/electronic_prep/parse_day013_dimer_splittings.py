from pathlib import Path
import re, csv

pairs = {
    "PYR2_PYR3": "runs/phase1A/day013_orca_dimers/PYR2_PYR3_pair_serial_tight/PYR2_PYR3_pair_serial_tight.out",
    "PYR2_PYR4": "runs/phase1A/day012_orca_pilots/PYR2_PYR4_pair_serial_tight/PYR2_PYR4_pair_serial_tight.out",
    "PYR2_PYR5": "runs/phase1A/day013_orca_dimers/PYR2_PYR5_pair_serial_tight/PYR2_PYR5_pair_serial_tight.out",
    "PYR3_PYR4": "runs/phase1A/day013_orca_dimers/PYR3_PYR4_pair_serial_tight/PYR3_PYR4_pair_serial_tight.out",
    "PYR3_PYR5": "runs/phase1A/day012_orca_pilots/PYR3_PYR5_pair_serial_tight/PYR3_PYR5_pair_serial_tight.out",
    "PYR4_PYR5": "runs/phase1A/day013_orca_dimers/PYR4_PYR5_pair_serial_tight/PYR4_PYR5_pair_serial_tight.out",
}

rows = []
pat = re.compile(r"STATE\s+(\d+):\s+E=\s+([0-9.]+)\s+au\s+([0-9.]+)\s+eV\s+([0-9.]+)\s+cm\*\*-1")

for pair, fname in pairs.items():
    text = Path(fname).read_text(errors="ignore")
    states = []
    for line in text.splitlines():
        m = pat.search(line)
        if m:
            states.append(float(m.group(3)))

    e1, e2 = states[0], states[1]
    splitting = e2 - e1
    J = splitting / 2.0

    rows.append({
        "pair": pair.replace("_", "-"),
        "E1_eV": f"{e1:.6f}",
        "E2_eV": f"{e2:.6f}",
        "splitting_eV": f"{splitting:.6f}",
        "J_eff_eV": f"{J:.6f}",
        "J_eff_meV": f"{J*1000:.3f}",
        "terminated_normally": "YES" if "ORCA TERMINATED NORMALLY" in text else "NO",
    })

out = Path("runs/phase1A/day013_orca_dimers/parsed_summary/dimer_splittings_and_Jeff.csv")
with out.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys())
    w.writeheader()
    w.writerows(rows)

print(out.read_text())
