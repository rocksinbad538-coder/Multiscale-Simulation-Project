from pathlib import Path
import pandas as pd

framedir = Path("runs/phase1A/day016_md_bath_extraction/frames_0to100ps_every5ps")
out = Path("runs/phase1A/day016_md_bath_extraction/extracted_frames_audit.csv")

rows = []
for f in sorted(framedir.glob("frame*.gro")):
    lines = f.read_text(errors="ignore").splitlines()
    natoms = int(lines[1].strip())
    atoms = lines[2:-1]
    pyr_count = 0
    hbn_count = 0
    sol_count = 0
    for line in atoms:
        resname = line[5:10].strip()
        if resname == "PYR":
            pyr_count += 1
        elif resname == "HBN":
            hbn_count += 1
        elif resname == "SOL":
            sol_count += 1
    rows.append({
        "frame_file": f.name,
        "natoms": natoms,
        "HBN_atoms": hbn_count,
        "PYR_atoms": pyr_count,
        "SOL_atoms": sol_count,
    })

df = pd.DataFrame(rows)
df.to_csv(out, index=False)
print(df.to_string(index=False))
print("Wrote:", out)
