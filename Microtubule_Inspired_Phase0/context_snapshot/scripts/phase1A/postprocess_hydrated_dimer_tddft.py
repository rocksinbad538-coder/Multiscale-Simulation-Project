from pathlib import Path
import argparse
import re
import csv

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--tag", required=True)
    args = ap.parse_args()

    out = Path(args.out)
    workdir = Path(args.workdir)
    tag = args.tag
    text = out.read_text(errors="ignore").splitlines()

    parsed = workdir / "parsed"
    parsed.mkdir(parents=True, exist_ok=True)

    state_pat = re.compile(
        r"STATE\s+(\d+):\s+E=\s+([0-9.]+)\s+au\s+([0-9.]+)\s+eV\s+([0-9.]+)\s+cm"
    )

    states = []
    for line in text:
        m = state_pat.search(line)
        if m:
            states.append([m.group(1), m.group(2), m.group(3), m.group(4)])

    states_csv = parsed / f"{tag}_tddft_states.csv"
    with states_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["state", "energy_au", "energy_eV", "energy_cm-1"])
        w.writerows(states)

    def extract_block(marker, nlines, filename):
        for i, line in enumerate(text):
            if marker in line:
                block = "\n".join(text[i:i+nlines]) + "\n"
                p = parsed / filename
                p.write_text(block)
                return p
        return None

    core = extract_block("TD-DFT/TDA EXCITED STATES", 140, f"{tag}_tddft_core_block.txt")
    absb = extract_block("ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS", 45, f"{tag}_absorption_block.txt")

    print(states_csv)
    if core:
        print(core)
    if absb:
        print(absb)

if __name__ == "__main__":
    main()
