#!/usr/bin/env python3
from pathlib import Path
import parmed as pmd
import csv

IN_TOP = Path("parameters/phase1A/accepted/pyrene_gaff2/pyrene_neutral.prmtop")
IN_CRD = Path("parameters/phase1A/accepted/pyrene_gaff2/pyrene_neutral.inpcrd")
OUTDIR = Path("parameters/phase1A/pyrene_gromacs")
OUTDIR.mkdir(parents=True, exist_ok=True)

gro = OUTDIR / "pyrene.gro"
top = OUTDIR / "pyrene.top"
itp = OUTDIR / "pyrene.itp"
audit = OUTDIR / "audit_pyrene_gromacs.csv"

amber = pmd.load_file(str(IN_TOP), str(IN_CRD))

# Save combined GROMACS topology and coordinate file
amber.save(str(gro), overwrite=True)
amber.save(str(top), overwrite=True)

# Create a reusable ITP by extracting moleculetype section from TOP
top_text = top.read_text()
marker = "[ moleculetype ]"
system_marker = "[ system ]"

if marker not in top_text:
    raise SystemExit("ERROR: [ moleculetype ] not found in pyrene.top")

start = top_text.index(marker)
end = top_text.index(system_marker) if system_marker in top_text else len(top_text)
itp_text = top_text[start:end].strip() + "\n"
itp.write_text(itp_text)

total_charge = sum(atom.charge for atom in amber.atoms)
atom_types = sorted(set(atom.type for atom in amber.atoms))

with audit.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric", "value"])
    w.writerow(["input_prmtop_exists", IN_TOP.exists()])
    w.writerow(["input_inpcrd_exists", IN_CRD.exists()])
    w.writerow(["gro_exists", gro.exists()])
    w.writerow(["top_exists", top.exists()])
    w.writerow(["itp_exists", itp.exists()])
    w.writerow(["atom_count", len(amber.atoms)])
    w.writerow(["residue_count", len(amber.residues)])
    w.writerow(["bond_count", len(amber.bonds)])
    w.writerow(["angle_count", len(amber.angles)])
    w.writerow(["dihedral_count", len(amber.dihedrals)])
    w.writerow(["unique_atom_types", ";".join(atom_types)])
    w.writerow(["total_charge_e", f"{total_charge:.8f}"])
    w.writerow(["abs_total_charge_e", f"{abs(total_charge):.8e}"])

print(audit.read_text())
