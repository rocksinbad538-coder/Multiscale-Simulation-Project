#!/usr/bin/env python3

from pathlib import Path
import re
import numpy as np
import pandas as pd

BASE = Path("runs/phase2/day053_dense_electronic_pilot")
ORCA_ROOT = BASE / "orca_inputs"

HIST_ROOT = Path(
    "runs/phase1A/day016_md_bath_extraction/"
    "orca_embedding_pilot_inputs"
)

OUTDIR = BASE / "analysis"
OUTDIR.mkdir(parents=True, exist_ok=True)

HOMO = 52
LUMO = 53
MIN_WEIGHT = 0.70
SITES = ["PYR3", "PYR4"]

STATE_RE = re.compile(
    r"STATE\s+(\d+):\s+E=\s+"
    r"([-+0-9.Ee]+)\s+au\s+"
    r"([-+0-9.Ee]+)\s+eV"
)

TRANSITION_RE = re.compile(
    r"(\d+)a\s*->\s*(\d+)a\s*:\s*"
    r"([-+0-9.Ee]+)\s*"
    r"\(c=\s*([-+0-9.Ee]+)\)"
)

ABSORPTION_RE = re.compile(
    r"0-1A\s+->\s+(\d+)-1A\s+"
    r"([-+0-9.Ee]+)\s+"
    r"([-+0-9.Ee]+)\s+"
    r"([-+0-9.Ee]+)\s+"
    r"([-+0-9.Ee]+)"
)

ERROR_PATTERNS = [
    re.compile(r"ORCA finished by error termination", re.I),
    re.compile(r"error termination", re.I),
    re.compile(r"aborting the run", re.I),
    re.compile(r"SCF NOT CONVERGED", re.I),
    re.compile(r"SCF failed", re.I),
    re.compile(r"segmentation fault", re.I),
]


def parse_states(path):
    text = path.read_text(errors="ignore")
    lines = text.splitlines()

    states = {}
    current = None

    for line in lines:
        m = STATE_RE.search(line)
        if m:
            root = int(m.group(1))
            current = root
            states[root] = {
                "energy_eV": float(m.group(3)),
                "transitions": [],
                "fosc": np.nan,
            }
            continue

        m = TRANSITION_RE.search(line)
        if m and current is not None:
            occ, virt, weight, coeff = m.groups()
            states[current]["transitions"].append({
                "occupied": int(occ),
                "virtual": int(virt),
                "weight": float(weight),
                "coefficient": float(coeff),
            })

    in_abs = False

    for line in lines:
        if "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in line:
            in_abs = True
            continue

        if in_abs and "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY" in line:
            break

        if in_abs:
            m = ABSORPTION_RE.search(line)
            if m:
                root = int(m.group(1))
                if root in states:
                    states[root]["fosc"] = float(m.group(5))

    return text, states


def hl_weight(state):
    return float(sum(
        t["weight"]
        for t in state["transitions"]
        if t["occupied"] == HOMO and t["virtual"] == LUMO
    ))


def analyze_job(frame_id, time_ps, site, out_path, source):
    text, states = parse_states(out_path)

    normal = "ORCA TERMINATED NORMALLY" in text
    scf = "SCF CONVERGED" in text
    tddft = "ORCA-CIS/TD-DFT FINISHED WITHOUT ERROR" in text
    error = any(p.search(text) for p in ERROR_PATTERNS)

    if not (normal and scf and tddft and not error):
        raise RuntimeError(f"QC FAIL: {out_path}")

    for root in (1, 2):
        if root not in states:
            raise RuntimeError(f"Missing S{root}: {out_path}")

    candidates = []

    for root in (1, 2):
        s = states[root]
        candidates.append({
            "root": root,
            "energy_eV": s["energy_eV"],
            "fosc": s["fosc"],
            "hl_weight": hl_weight(s),
        })

    tracked = max(candidates, key=lambda x: x["hl_weight"])
    alternate = min(candidates, key=lambda x: x["hl_weight"])

    return {
        "frame_id": frame_id,
        "time_ps": time_ps,
        "site": site,
        "source": source,
        "tracked_root": tracked["root"],
        "tracked_energy_eV": tracked["energy_eV"],
        "tracked_fosc": tracked["fosc"],
        "tracked_HOMO_LUMO_weight": tracked["hl_weight"],
        "alternate_root": alternate["root"],
        "alternate_energy_eV": alternate["energy_eV"],
        "alternate_HOMO_LUMO_weight": alternate["hl_weight"],
        "state_separation_meV": 1000.0 * abs(
            candidates[1]["energy_eV"] - candidates[0]["energy_eV"]
        ),
        "S1_energy_eV": candidates[0]["energy_eV"],
        "S1_HOMO_LUMO_weight": candidates[0]["hl_weight"],
        "S2_energy_eV": candidates[1]["energy_eV"],
        "S2_HOMO_LUMO_weight": candidates[1]["hl_weight"],
    }


rows = []

# Historical regression anchor at exactly 100.000 ps.
for site in SITES:
    job = f"frame020_{site}_embedding"
    path = HIST_ROOT / job / f"{job}.out"

    if not path.is_file():
        raise SystemExit(f"Missing historical anchor: {path}")

    rows.append(
        analyze_job(
            frame_id=0,
            time_ps=100.0,
            site=site,
            out_path=path,
            source="M4_frame020_anchor",
        )
    )

# New dense jobs: 100.05 ... 100.50 ps.
for frame_id in range(1, 11):
    time_ps = 100.0 + 0.05 * frame_id

    for site in SITES:
        job = f"dense{frame_id:03d}_{site}"
        path = ORCA_ROOT / job / f"{job}.out"

        if not path.is_file():
            raise SystemExit(f"Missing new output: {path}")

        rows.append(
            analyze_job(
                frame_id=frame_id,
                time_ps=time_ps,
                site=site,
                out_path=path,
                source="Day053_dense",
            )
        )

df = pd.DataFrame(rows).sort_values(
    ["frame_id", "site"]
).reset_index(drop=True)

if len(df) != 22:
    raise SystemExit(f"Expected 22 rows, found {len(df)}")

if df.duplicated(["frame_id", "site"]).any():
    raise SystemExit("Duplicate frame/site rows found")

if (df["tracked_HOMO_LUMO_weight"] < MIN_WEIGHT).any():
    print(
        df.loc[
            df["tracked_HOMO_LUMO_weight"] < MIN_WEIGHT,
            [
                "frame_id",
                "time_ps",
                "site",
                "tracked_root",
                "tracked_HOMO_LUMO_weight",
            ],
        ].to_string(index=False)
    )
    raise SystemExit("STATE_TRACKING_GATE=FAIL")

tracking_path = OUTDIR / "DAY053_STATE_TRACKING.csv"
df.to_csv(tracking_path, index=False)

energy = df.pivot(
    index="time_ps",
    columns="site",
    values="tracked_energy_eV",
).reset_index()

energy = energy[["time_ps", "PYR3", "PYR4"]]

energy["gap_PYR3_minus_PYR4_eV"] = (
    energy["PYR3"] - energy["PYR4"]
)

energy["gap_PYR3_minus_PYR4_meV"] = (
    1000.0 * energy["gap_PYR3_minus_PYR4_eV"]
)

for col in ["PYR3", "PYR4", "gap_PYR3_minus_PYR4_eV"]:
    energy[f"delta_{col}"] = energy[col] - energy[col].mean()

energy_path = OUTDIR / "DAY053_DENSE_GAP_TRAJECTORY.csv"
energy.to_csv(energy_path, index=False)

roots = (
    df.groupby("site")["tracked_root"]
    .apply(lambda x: sorted(set(int(v) for v in x)))
    .to_dict()
)

weights = (
    df.groupby("site")["tracked_HOMO_LUMO_weight"]
    .agg(["min", "mean", "max"])
)

print("=" * 78)
print("DAY053 DENSE ELECTRONIC PILOT ANALYSIS")
print("=" * 78)
print("N_ROWS =", len(df))
print("N_TIMES =", energy["time_ps"].nunique())
print("TIME_MIN_PS =", energy["time_ps"].min())
print("TIME_MAX_PS =", energy["time_ps"].max())
print("ROOTS =", roots)
print("MIN_CHARACTER_WEIGHT =", df["tracked_HOMO_LUMO_weight"].min())
print("STATE_TRACKING_GATE=PASS")
print()
print("CHARACTER WEIGHTS")
print(weights.to_string())
print()
print("GAP_MEAN_MEV =", energy["gap_PYR3_minus_PYR4_meV"].mean())
print("GAP_SD_MEV =", energy["gap_PYR3_minus_PYR4_meV"].std(ddof=1))
print("GAP_MIN_MEV =", energy["gap_PYR3_minus_PYR4_meV"].min())
print("GAP_MAX_MEV =", energy["gap_PYR3_minus_PYR4_meV"].max())
print()
print("TRACKING =", tracking_path)
print("GAP_TRAJECTORY =", energy_path)
