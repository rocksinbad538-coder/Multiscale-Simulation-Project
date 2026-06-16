#!/usr/bin/env python3
from pathlib import Path

inp = Path("parameters/phase1A/pyrene_gaff2/pyrene_gaff2.mol2")
out = Path("parameters/phase1A/pyrene_gaff2/pyrene_gaff2_neutral.mol2")

lines = inp.read_text().splitlines()
atom_start = None
bond_start = None

for i, line in enumerate(lines):
    if line.startswith("@<TRIPOS>ATOM"):
        atom_start = i + 1
    if line.startswith("@<TRIPOS>BOND"):
        bond_start = i
        break

if atom_start is None or bond_start is None:
    raise SystemExit("Could not find MOL2 ATOM/BOND sections.")

atom_lines = lines[atom_start:bond_start]
charges = []
parsed = []

for line in atom_lines:
    parts = line.split()
    if len(parts) < 9:
        parsed.append((line, None))
        continue
    q = float(parts[8])
    charges.append(q)
    parsed.append((line, parts))

total_q = sum(charges)
n = len(charges)
corr = -total_q / n

new_atom_lines = []
for line, parts in parsed:
    if parts is None:
        new_atom_lines.append(line)
        continue

    parts[8] = f"{float(parts[8]) + corr:.6f}"

    # keep a readable mol2 line
    new_line = (
        f"{int(parts[0]):7d} {parts[1]:<8s}"
        f"{float(parts[2]):10.4f} {float(parts[3]):10.4f} {float(parts[4]):10.4f} "
        f"{parts[5]:<8s} {parts[6]:>3s} {parts[7]:<8s} {parts[8]:>10s}"
    )
    new_atom_lines.append(new_line)

new_lines = lines[:atom_start] + new_atom_lines + lines[bond_start:]
out.write_text("\n".join(new_lines) + "\n")

print("Input total charge:", f"{total_q:.8f}")
print("Correction per atom:", f"{corr:.10f}")
print("Wrote:", out)

# verify
q2 = []
in_atom = False
for line in out.read_text().splitlines():
    if line.startswith("@<TRIPOS>ATOM"):
        in_atom = True
        continue
    if line.startswith("@<TRIPOS>BOND"):
        in_atom = False
        continue
    if in_atom and line.strip():
        parts = line.split()
        if len(parts) >= 9:
            q2.append(float(parts[8]))
print("Output total charge:", f"{sum(q2):.8f}")
