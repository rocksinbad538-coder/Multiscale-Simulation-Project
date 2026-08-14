#!/usr/bin/env python3

from pathlib import Path
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

RUNNING = ROOT/"runs"/"phase2"/"campaign"/"running_statistics"

rows=[]

for csv in sorted(RUNNING.glob("*K_running_statistics.csv")):

    T = csv.stem.replace("_running_statistics","")

    df = pd.read_csv(csv)

    i80 = int(0.8*len(df))

    result = {"Temperature":T}

    for var in ["Rg","PotEng","RMSD"]:

        run = df[f"{var}_RunningMean"]

        final = run.iloc[-1]

        early = run.iloc[i80]

        drift = abs(final-early)

        rel = 100*drift/max(abs(final),1e-12)

        result[f"{var}_DriftPercent"]=rel

    rows.append(result)

summary=pd.DataFrame(rows)

outfile=RUNNING/"convergence_assessment.csv"

summary.to_csv(outfile,index=False)

(RUNNING/"convergence_assessment.json").write_text(
    json.dumps(
        summary.to_dict(orient="records"),
        indent=2
    )
)

print("="*90)
print("DAY047 / PHASE2-A45")
print("CONVERGENCE ASSESSMENT")
print("="*90)
print(outfile)
