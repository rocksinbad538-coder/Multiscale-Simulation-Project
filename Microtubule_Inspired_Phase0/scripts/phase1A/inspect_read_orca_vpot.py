#!/usr/bin/env python3

from pathlib import Path
from pprint import pprint

from resp_common import read_orca_vpot

ROOT = Path(__file__).resolve().parents[2]

VPOT = (
    ROOT
    / "runs/phase1A/day036_qm_f06_upper_v7a_r1_esp_executions"
    / "esp_upper_v7a_r1_20260731T174832Z"
    / "esp_upper_v7a_r1.vpot"
)

obj = read_orca_vpot(VPOT)

print("=" * 80)
print(type(obj))
print("=" * 80)

if isinstance(obj, dict):

    print("\nKEYS\n")
    pprint(sorted(obj.keys()))

    print("\n")

    for k, v in obj.items():

        print("-" * 80)
        print(k)
        print(type(v))

        try:
            print("len =", len(v))
        except Exception:
            pass

        try:
            print("shape =", v.shape)
        except Exception:
            pass

        print()

else:
    print(obj)
