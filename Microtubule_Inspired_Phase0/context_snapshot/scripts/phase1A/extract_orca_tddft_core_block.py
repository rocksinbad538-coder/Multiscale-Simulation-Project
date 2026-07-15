from pathlib import Path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--out", required=True)
parser.add_argument("--txt", required=True)
args = parser.parse_args()

out = Path(args.out)
txt = Path(args.txt)
txt.parent.mkdir(parents=True, exist_ok=True)

lines = out.read_text(errors="ignore").splitlines()

start = None
end = None

for i, line in enumerate(lines):
    if "TD-DFT/TDA EXCITED STATES" in line:
        start = max(0, i - 20)
    if start is not None and "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS" in line:
        end = i + 35
        break

if start is None:
    raise SystemExit("Could not find TDDFT excited-state block")

if end is None:
    end = min(len(lines), start + 250)

txt.write_text("\n".join(lines[start:end]) + "\n")
print(txt)
