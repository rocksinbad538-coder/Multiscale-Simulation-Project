#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import json

ROOT = Path(__file__).resolve().parents[2]

analysis = ROOT/"runs"/"phase2"/"day045_md_analysis"

df = pd.read_csv(analysis/"incremental_rmsd.csv")

top = (
    df.sort_values("IncrementalRMSD", ascending=False)
      .head(20)
      .reset_index(drop=True)
)

top.to_csv(
    analysis/"incremental_rmsd_largest_events.csv",
    index=False
)

report = {
    "maximum_incremental_rmsd_A": float(top.iloc[0]["IncrementalRMSD"]),
    "timestep_of_maximum": int(top.iloc[0]["timestep"]),
    "events_above_1A": int((df["IncrementalRMSD"] > 1.0).sum()),
    "events_above_2A": int((df["IncrementalRMSD"] > 2.0).sum()),
    "events_above_3A": int((df["IncrementalRMSD"] > 3.0).sum()),
}

(analysis/"INCREMENTAL_EVENT_AUDIT.json").write_text(
    json.dumps(report, indent=2)
)

print("="*90)
print("DAY047 / PHASE2-A38")
print("INCREMENTAL EVENT AUDIT")
print("="*90)

print(analysis/"incremental_rmsd_largest_events.csv")
print(analysis/"INCREMENTAL_EVENT_AUDIT.json")
