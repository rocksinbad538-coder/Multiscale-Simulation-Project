from pathlib import Path
import argparse
import csv
import re

parser = argparse.ArgumentParser()
parser.add_argument("--out", required=True, help="ORCA output file")
parser.add_argument("--csv", required=True, help="Output CSV")
args = parser.parse_args()

out = Path(args.out)
csv_path = Path(args.csv)
csv_path.parent.mkdir(parents=True, exist_ok=True)

state_pat = re.compile(
    r"STATE\s+(\d+):\s+E=\s+([0-9.]+)\s+au\s+([0-9.]+)\s+eV\s+([0-9.]+)\s+cm"
)

rows = []
for line in out.read_text(errors="ignore").splitlines():
    m = state_pat.search(line)
    if m:
        rows.append({
            "state": f"S{m.group(1)}",
            "energy_au": m.group(2),
            "energy_eV": m.group(3),
            "energy_cm-1": m.group(4),
        })

with csv_path.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["state", "energy_au", "energy_eV", "energy_cm-1"])
    w.writeheader()
    w.writerows(rows)

print(csv_path)
print(f"Parsed {len(rows)} states")
