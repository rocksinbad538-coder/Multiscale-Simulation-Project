from pathlib import Path
import re
import pandas as pd

ROOT = Path("runs/phase1A/day016_md_bath_extraction/orca_embedding_pilot_inputs")
OUTDIR = Path("runs/phase1A/day016_md_bath_extraction/orca_embedding_analysis")
OUTDIR.mkdir(parents=True, exist_ok=True)

state_re = re.compile(
    r"STATE\s+(\d+):\s+E=\s+([-0-9.]+)\s+au\s+([-0-9.]+)\s+eV\s+([-0-9.]+)\s+cm\*\*-1"
)

abs_re = re.compile(
    r"0-1A\s+->\s+(\d+)-1A\s+([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.Ee+]+)"
)

def parse_frame_cluster(name: str):
    m = re.match(r"frame(\d+)_([^_]+)_embedding", name)
    if not m:
        return None, None
    return int(m.group(1)), m.group(2)

def parse_pc(pc_path: Path):
    if not pc_path.exists():
        return None, None
    lines = [x.strip() for x in pc_path.read_text(errors="ignore").splitlines() if x.strip()]
    if not lines:
        return None, None
    try:
        declared = int(float(lines[0].split()[0]))
    except Exception:
        declared = None
    qsum = 0.0
    n = 0
    for line in lines[1:]:
        p = line.split()
        if len(p) == 4:
            try:
                qsum += float(p[0])
                n += 1
            except Exception:
                pass
    return declared if declared is not None else n, qsum

rows = []

for out in sorted(ROOT.glob("frame*_PYR*_embedding/frame*_PYR*_embedding.out")):
    name = out.parent.name
    frame, cluster = parse_frame_cluster(name)
    txt = out.read_text(errors="ignore")

    ok = "ORCA TERMINATED NORMALLY" in txt
    tddft_ok = "ORCA-CIS/TD-DFT FINISHED WITHOUT ERROR" in txt
    scf_ok = "SCF CONVERGED" in txt
    has_error = bool(re.search(r"\bERROR\b|aborting the run|error termination", txt, re.I))

    final_energy = None
    for line in txt.splitlines():
        if "FINAL SINGLE POINT ENERGY" in line:
            try:
                final_energy = float(line.split()[-1])
            except Exception:
                pass

    total_minutes = None
    m = re.search(r"TOTAL RUN TIME:\s+(\d+)\s+days\s+(\d+)\s+hours\s+(\d+)\s+minutes\s+([0-9.]+)\s+seconds", txt)
    if m:
        d, h, mi, s = m.groups()
        total_minutes = int(d)*1440 + int(h)*60 + int(mi) + float(s)/60.0

    pc_file = out.parent / (out.stem + ".pc")
    n_pc_file, pc_total_charge = parse_pc(pc_file)

    pc_reads = re.findall(r"ok \((\d+) point charges\)", txt)
    n_pc_orca = int(pc_reads[-1]) if pc_reads else None

    row = {
        "frame": frame,
        "cluster": cluster,
        "job": name,
        "terminated_normally": ok,
        "scf_converged": scf_ok,
        "tddft_finished": tddft_ok,
        "has_error_flag": has_error,
        "final_single_point_energy_Eh": final_energy,
        "n_point_charges_file": n_pc_file,
        "n_point_charges_orca": n_pc_orca,
        "point_charge_total": pc_total_charge,
        "total_runtime_min": total_minutes,
        "out_file": str(out),
    }

    states = state_re.findall(txt)
    for state_id, e_au, e_ev, e_cm in states[:10]:
        k = int(state_id)
        row[f"S{k}_au"] = float(e_au)
        row[f"S{k}_eV"] = float(e_ev)
        row[f"S{k}_cm-1"] = float(e_cm)

    in_abs = False
    for line in txt.splitlines():
        if "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in line:
            in_abs = True
            continue
        if in_abs and "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY" in line:
            break
        if in_abs:
            m = abs_re.search(line)
            if m:
                state, e_ev, e_cm, wl_nm, fosc = m.groups()
                k = int(state)
                row[f"f{k}"] = float(fosc)
                row[f"lambda{k}_nm"] = float(wl_nm)

    dip = re.search(r"Magnitude \(Debye\)\s+:\s+([-0-9.]+)", txt)
    row["dipole_D"] = float(dip.group(1)) if dip else None

    rows.append(row)

df = pd.DataFrame(rows).sort_values(["frame", "cluster"])
csv_path = OUTDIR / "embedding_pilot_summary.csv"
df.to_csv(csv_path, index=False)

successful = int((df["terminated_normally"] & df["tddft_finished"]).sum()) if len(df) else 0

md = OUTDIR / "EMBEDDING_PILOT_AUDIT_DAY016.md"
with md.open("w") as f:
    f.write("# Day016 ORCA embedding pilot audit\n\n")
    f.write(f"- Jobs parsed: {len(df)}\n")
    f.write(f"- Successful embedded TDDFT jobs: {successful}/{len(df)}\n")
    if len(df):
        f.write(f"- Point charges read by ORCA: {df['n_point_charges_orca'].min()}–{df['n_point_charges_orca'].max()}\n")
        f.write(f"- S1 range: {df['S1_eV'].min():.3f}–{df['S1_eV'].max():.3f} eV\n")
        f.write(f"- Mean S1: {df['S1_eV'].mean():.3f} eV\n")
        f.write(f"- Std S1: {df['S1_eV'].std(ddof=1):.3f} eV\n")
    f.write("\n## Parsed jobs\n\n")
    cols = ["frame","cluster","terminated_normally","tddft_finished","n_point_charges_orca","S1_eV","S2_eV","S3_eV","f1","dipole_D","total_runtime_min"]
    f.write(df[cols].to_string(index=False))
    f.write("\n")

print(df[["frame","cluster","terminated_normally","tddft_finished","n_point_charges_orca","S1_eV","S2_eV","S3_eV","f1","dipole_D","total_runtime_min"]].to_string(index=False))
print("Wrote:", csv_path)
print("Wrote:", md)
