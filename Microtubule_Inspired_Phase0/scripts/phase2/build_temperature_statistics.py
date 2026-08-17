#!/usr/bin/env python3

from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

META = ROOT/"runs"/"phase2"/"campaign"/"campaign_meta_summary.csv"

OUT = ROOT/"runs"/"phase2"/"campaign"/"temperature_statistics"

OUT.mkdir(exist_ok=True)

df = pd.read_csv(META)

variables=[

("final_Rg_A","Radius of gyration (Å)"),

("mean_RMSF_A","Mean RMSF (Å)"),

("maximum_RMSF_A","Maximum RMSF (Å)"),

("mean_aligned_rmsd_A","Mean aligned RMSD (Å)"),

("mean_relative_shape_anisotropy","Relative shape anisotropy"),

("maximum_atomic_displacement_A","Maximum displacement (Å)")

]

summary=[]

for column,label in variables:

    x=df["Temperature_K"].values.astype(float)

    y=df[column].values.astype(float)

    p=np.polyfit(x,y,1)

    fit=np.polyval(p,x)

    ss_res=((y-fit)**2).sum()

    ss_tot=((y-y.mean())**2).sum()

    r2=1.0-ss_res/ss_tot if ss_tot>0 else 1.0

    summary.append({

        "descriptor":column,

        "slope_per_K":float(p[0]),

        "intercept":float(p[1]),

        "R2":float(r2),

        "minimum":float(y.min()),

        "maximum":float(y.max()),

        "mean":float(y.mean()),

        "std":float(y.std(ddof=1))

    })

    plt.figure(figsize=(6,4))

    plt.scatter(x,y,s=60)

    plt.plot(x,fit,lw=2)

    plt.grid(alpha=.3)

    plt.xlabel("Temperature (K)")

    plt.ylabel(label)

    plt.title(label)

    plt.tight_layout()

    plt.savefig(

        OUT/f"{column}.png",

        dpi=300

    )

    plt.close()

summary_df=pd.DataFrame(summary)

summary_df.to_csv(

    OUT/"temperature_statistics.csv",

    index=False

)

with open(OUT/"temperature_statistics.json","w") as f:

    json.dump(

        summary,

        f,

        indent=2

    )

print("="*90)
print("DAY049 / PHASE2-B10")
print("TEMPERATURE STATISTICS")
print("="*90)

print(summary_df)

print()

print(OUT)
