#!/usr/bin/env python3

from pathlib import Path
import csv
from collections import Counter

ROOT = Path(__file__).resolve().parents[2]

INPUT = (
ROOT
/"runs"
/
"phase1A"
/
"day038_resp_stage1_executions"
/
"resp_stage1_upper_v7a_r1_20260803T202335Z"
/
"QM_F06_UPPER_V7A_R1_RESP_STAGE1_TRANSFERABILITY.csv"
)

rows=[]

elements=[]

ids=[]

with open(INPUT) as f:

    reader=csv.DictReader(f)

    for r in reader:

        rows.append(r)

        elements.append(r["element"])

        ids.append(r["atom_id"])

print("="*100)

print("TRANSFERABILITY AUDIT")

print("="*100)

print()

print("rows =",len(rows))

print()

print("elements")

for k,v in Counter(elements).items():

    print(k,v)

print()

print("unique atom ids =",len(set(ids)))

print("duplicate ids =",len(ids)-len(set(ids)))
