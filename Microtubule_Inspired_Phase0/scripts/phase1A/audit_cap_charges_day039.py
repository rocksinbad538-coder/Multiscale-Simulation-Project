#!/usr/bin/env python3

from pathlib import Path
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

latest = (
    ROOT
    / "runs"
    / "phase1A"
    / "day038_resp_stage1_executions"
    / "LATEST_RESP_STAGE1_UPPER_V7A_R1_EXECUTION.txt"
).read_text().strip()

run = Path(latest)

csv = run / "QM_F06_UPPER_V7A_R1_RESP_STAGE1_TRANSFERABILITY.csv"

df = pd.read_csv(csv)

caps = df[df["artificial_cap_bool"]].copy()

print("=" * 100)
print("DAY039 / D039-A1")
print("CAP CHARGE INVENTORY")
print("=" * 100)
print()

print("Artificial caps :", len(caps))
print()

total = caps["RESP_stage1_charge_e_float"].sum()

print(f"Total cap charge = {total:.8f} e")
print()

cols = [
    "atom_index_1based",
    "element",
    "atom_role",
    "RESP_stage1_charge_e_float",
]

print(caps[cols].to_string(index=False))

out = run / "DAY039_CAP_CHARGE_INVENTORY.json"

json.dump(
    {
        "cap_count": int(len(caps)),
        "total_cap_charge": float(total),
    },
    out.open("w"),
    indent=2,
)

print()
print("written:", out)
