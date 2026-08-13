#!/usr/bin/env python3

from pathlib import Path
import subprocess
import shutil

ROOT = Path(__file__).resolve().parents[2]

CAMPAIGN = ROOT / "runs" / "phase2" / "campaign"

WORK = ROOT / "runs" / "phase2" / "day044_md_protocol"

ANALYSIS = ROOT / "runs" / "phase2" / "day045_md_analysis"

scripts = [
    "analyze_md_trajectory.py",
    "analyze_md_log.py",
    "analyze_rmsf.py",
    "analyze_shape.py",
    "analyze_aligned_rmsd.py",
    "build_md_scientific_report.py",
]

for folder in sorted(CAMPAIGN.iterdir()):

    if not folder.is_dir():
        continue

    print("="*80)
    print(folder.name)
    print("="*80)

    # limpiar directorio temporal
    for f in ANALYSIS.glob("*"):
        if f.is_file():
            f.unlink()

    # copiar resultados de la campaña
    shutil.copy2(folder/"production.xyz", WORK/"production.xyz")
    shutil.copy2(folder/"production.log", WORK/"in.production.log")

    # ejecutar todos los análisis
    for s in scripts:

        subprocess.run(
            ["python3", str(ROOT/"scripts"/"phase2"/s)],
            check=True,
        )

    # copiar resultados al directorio de la campaña
    dest = folder/"analysis"

    dest.mkdir(exist_ok=True)

    for f in ANALYSIS.iterdir():
        shutil.copy2(f, dest/f.name)

print()
print("="*80)
print("MULTI-TEMPERATURE ANALYSIS COMPLETED")
print("="*80)
