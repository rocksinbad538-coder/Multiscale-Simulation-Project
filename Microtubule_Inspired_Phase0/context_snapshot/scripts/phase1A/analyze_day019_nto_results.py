#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from pathlib import Path
from statistics import mean, pstdev

PROJECT_ROOT = Path(__file__).resolve().parents[2]
NTO_ROOT = PROJECT_ROOT / "runs/phase1A/day016_md_bath_extraction/day019_nto_inputs"
MANIFEST = NTO_ROOT / "NTO_INPUT_MANIFEST_DAY019.csv"
OUTPUT_ROOT = PROJECT_ROOT / "runs/phase1A/day016_md_bath_extraction/day019_nto_analysis"

STATE_METRICS_CSV = OUTPUT_ROOT / "nto_state_metrics.csv"
PAIR_OCCUPATIONS_CSV = OUTPUT_ROOT / "nto_pair_occupations.csv"
ENVIRONMENT_COMPARISON_CSV = OUTPUT_ROOT / "nto_environment_comparison.csv"
MODEL_ASSESSMENT_CSV = OUTPUT_ROOT / "model_space_assessment.csv"
REPORT_MD = OUTPUT_ROOT / "DAY019_NTO_ANALYSIS.md"

STATE_RE = re.compile(r"NATURAL TRANSITION ORBITALS FOR STATE\s+(\d+)")
ENERGY_RE = re.compile(
    r"^\s*E=\s*([-+0-9.Ee]+)\s+au\s+([-+0-9.Ee]+)\s+eV\s+([-+0-9.Ee]+)\s+cm\*\*-1"
)
OCC_RE = re.compile(
    r"^\s*(\d+)([ab])\s*->\s*(\d+)([ab])\s*:\s*n=\s*([-+0-9.Ee]+)"
)
ABS_RE = re.compile(
    r"^\s*\d+-\S+\s*->\s*(\d+)-\S+\s+"
    r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)"
)

NORMAL = "ORCA TERMINATED NORMALLY"
SCF = "SCF CONVERGED"
TDDFT = "ORCA-CIS/TD-DFT FINISHED WITHOUT ERROR"

SINGLE_PAIR_THRESHOLD = 0.75
TWO_PAIR_N1_MAX = 0.70
TWO_PAIR_N2_MIN = 0.25
LOW_GAP_LIMIT_MEV = 120.0


def read_manifest() -> list[dict[str, str]]:
    if not MANIFEST.is_file():
        raise SystemExit(f"Missing manifest: {MANIFEST}")
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 8:
        raise SystemExit(f"Expected 8 manifest rows, found {len(rows)}")
    return rows


def parse_fosc(lines: list[str]) -> dict[int, float]:
    values: dict[int, float] = {}
    inside = False
    for line in lines:
        if "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in line:
            inside = True
            continue
        if inside and "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY" in line:
            break
        if inside:
            match = ABS_RE.search(line)
            if match:
                values[int(match.group(1))] = float(match.group(5))
    return values


def classify(n1: float, n2: float) -> str:
    if n1 >= SINGLE_PAIR_THRESHOLD:
        return "single_pair_dominated"
    if n1 <= TWO_PAIR_N1_MAX and n2 >= TWO_PAIR_N2_MIN:
        return "two_pair_mixed"
    return "intermediate"


def parse_output(path: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    if not path.is_file():
        raise RuntimeError(f"Missing output: {path}")
    text = path.read_text(errors="ignore")
    for marker in (NORMAL, SCF, TDDFT):
        if marker not in text:
            raise RuntimeError(f"Missing marker '{marker}' in {path}")

    lines = text.splitlines()
    fosc = parse_fosc(lines)
    states: dict[int, dict[str, object]] = {}
    current: int | None = None

    for line in lines:
        match = STATE_RE.search(line)
        if match:
            current = int(match.group(1))
            states[current] = {"pairs": []}
            continue
        if current is None:
            continue
        match = ENERGY_RE.match(line)
        if match:
            states[current]["energy_au"] = float(match.group(1))
            states[current]["energy_eV"] = float(match.group(2))
            states[current]["energy_cm1"] = float(match.group(3))
            continue
        match = OCC_RE.match(line)
        if match:
            states[current]["pairs"].append(
                {
                    "root": current,
                    "hole_orbital": int(match.group(1)),
                    "hole_spin": match.group(2),
                    "particle_orbital": int(match.group(3)),
                    "particle_spin": match.group(4),
                    "occupation": float(match.group(5)),
                }
            )
            continue
        if "TD-DFT/TDA-EXCITATION SPECTRA" in line:
            current = None

    if sorted(states) != [1, 2]:
        raise RuntimeError(f"Expected states S1/S2 in {path}; found {sorted(states)}")

    state_rows: list[dict[str, object]] = []
    pair_rows: list[dict[str, object]] = []

    for root in (1, 2):
        state = states[root]
        pairs = list(state["pairs"])
        if not pairs or "energy_eV" not in state:
            raise RuntimeError(f"Incomplete NTO block for S{root} in {path}")
        pairs = sorted(pairs, key=lambda row: float(row["occupation"]), reverse=True)
        occupations = [float(row["occupation"]) for row in pairs]
        total = sum(occupations)
        normalized = [value / total for value in occupations if value > 0]
        pr = 1.0 / sum(value * value for value in normalized)
        entropy = -sum(value * math.log(value) for value in normalized)
        n1 = occupations[0]
        n2 = occupations[1] if len(occupations) > 1 else 0.0
        first, second = pairs[0], pairs[1]

        state_rows.append(
            {
                "root": root,
                "energy_au": state["energy_au"],
                "energy_eV": state["energy_eV"],
                "energy_cm1": state["energy_cm1"],
                "fosc": fosc.get(root, float("nan")),
                "n_printed_pairs": len(pairs),
                "sum_printed_occupations": total,
                "unprinted_weight_upper_bound": max(0.0, 1.0 - total),
                "n1": n1,
                "n2": n2,
                "n1_plus_n2": n1 + n2,
                "n2_over_n1": n2 / n1 if n1 else float("nan"),
                "participation_ratio_printed": pr,
                "shannon_entropy_printed": entropy,
                "character_class": classify(n1, n2),
                "dominant_pair": f"{first['hole_orbital']}{first['hole_spin']}->{first['particle_orbital']}{first['particle_spin']}",
                "second_pair": f"{second['hole_orbital']}{second['hole_spin']}->{second['particle_orbital']}{second['particle_spin']}",
            }
        )
        pair_rows.extend(pairs)

    return state_rows, pair_rows


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def f(value: object, digits: int = 6) -> str:
    if isinstance(value, float):
        return "nan" if math.isnan(value) else f"{value:.{digits}f}"
    return str(value)


def main() -> None:
    manifest = read_manifest()
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    state_rows: list[dict[str, object]] = []
    pair_rows: list[dict[str, object]] = []

    for item in manifest:
        job = item["target_job"]
        job_dir = NTO_ROOT / job
        output = job_dir / f"{job}.out"
        states, pairs = parse_output(output)
        tracked_root = int(item["tracked_root"])
        by_root = {int(row["root"]): row for row in states}
        gap_meV = (float(by_root[2]["energy_eV"]) - float(by_root[1]["energy_eV"])) * 1000.0
        brightest_root = max((1, 2), key=lambda root: float(by_root[root]["fosc"]))

        for state in states:
            root = int(state["root"])
            state_rows.append(
                {
                    "calculation_type": item["calculation_type"],
                    "frame": int(item["frame"]),
                    "cluster": item["cluster"],
                    "job": job,
                    "reason": item["reason"],
                    "tracked_root": tracked_root,
                    "brightest_root_S1_S2": brightest_root,
                    "tracked_equals_brightest": tracked_root == brightest_root,
                    "root": root,
                    "is_tracked_root": root == tracked_root,
                    "is_brightest_root": root == brightest_root,
                    "S1_S2_gap_meV": gap_meV,
                    **state,
                    "output": str(output.relative_to(PROJECT_ROOT)),
                    "nto_file": str((job_dir / f"{job}.s{root}.nto").relative_to(PROJECT_ROOT)),
                }
            )

        for pair in pairs:
            root = int(pair["root"])
            pair_rows.append(
                {
                    "calculation_type": item["calculation_type"],
                    "frame": int(item["frame"]),
                    "cluster": item["cluster"],
                    "job": job,
                    "tracked_root": tracked_root,
                    "is_tracked_root": root == tracked_root,
                    **pair,
                }
            )

    if len(state_rows) != 16:
        raise RuntimeError(f"Expected 16 states, found {len(state_rows)}")

    tracked = [row for row in state_rows if bool(row["is_tracked_root"])]
    alternate = [row for row in state_rows if not bool(row["is_tracked_root"])]
    if len(tracked) != 8 or len(alternate) != 8:
        raise RuntimeError("Expected 8 tracked and 8 alternate states")
    if not all(bool(row["tracked_equals_brightest"]) for row in tracked):
        bad = sorted({str(row["job"]) for row in tracked if not bool(row["tracked_equals_brightest"])})
        raise RuntimeError("Tracked root is not brightest for: " + ", ".join(bad))

    state_fields = [
        "calculation_type", "frame", "cluster", "job", "reason",
        "tracked_root", "brightest_root_S1_S2", "tracked_equals_brightest",
        "root", "is_tracked_root", "is_brightest_root", "S1_S2_gap_meV",
        "energy_au", "energy_eV", "energy_cm1", "fosc",
        "n_printed_pairs", "sum_printed_occupations", "unprinted_weight_upper_bound",
        "n1", "n2", "n1_plus_n2", "n2_over_n1",
        "participation_ratio_printed", "shannon_entropy_printed",
        "character_class", "dominant_pair", "second_pair", "output", "nto_file",
    ]
    pair_fields = [
        "calculation_type", "frame", "cluster", "job", "tracked_root",
        "is_tracked_root", "root", "hole_orbital", "hole_spin",
        "particle_orbital", "particle_spin", "occupation",
    ]
    write_csv(STATE_METRICS_CSV, state_rows, state_fields)
    write_csv(PAIR_OCCUPATIONS_CSV, pair_rows, pair_fields)

    vacuum = {
        (str(row["cluster"]), int(row["root"])): row
        for row in state_rows if row["calculation_type"] == "vacuum_reference"
    }
    environment_rows: list[dict[str, object]] = []
    for row in state_rows:
        if row["calculation_type"] != "embedded":
            continue
        ref = vacuum.get((str(row["cluster"]), int(row["root"])))
        if ref is None:
            continue
        environment_rows.append(
            {
                "frame": row["frame"], "cluster": row["cluster"], "root": row["root"],
                "is_tracked_root": row["is_tracked_root"],
                "embedded_energy_eV": row["energy_eV"], "vacuum_energy_eV": ref["energy_eV"],
                "delta_energy_meV": (float(row["energy_eV"]) - float(ref["energy_eV"])) * 1000.0,
                "embedded_n1": row["n1"], "vacuum_n1": ref["n1"],
                "delta_n1": float(row["n1"]) - float(ref["n1"]),
                "embedded_n2": row["n2"], "vacuum_n2": ref["n2"],
                "delta_n2": float(row["n2"]) - float(ref["n2"]),
                "embedded_participation_ratio": row["participation_ratio_printed"],
                "vacuum_participation_ratio": ref["participation_ratio_printed"],
                "delta_participation_ratio": float(row["participation_ratio_printed"]) - float(ref["participation_ratio_printed"]),
                "embedded_character_class": row["character_class"],
                "vacuum_character_class": ref["character_class"],
                "character_class_preserved": row["character_class"] == ref["character_class"],
            }
        )
    env_fields = list(environment_rows[0].keys()) if environment_rows else []
    if environment_rows:
        write_csv(ENVIRONMENT_COMPARISON_CSV, environment_rows, env_fields)

    gaps = sorted({round(float(row["S1_S2_gap_meV"]), 12) for row in state_rows})
    tracked_single_fraction = sum(row["character_class"] == "single_pair_dominated" for row in tracked) / 8
    alternate_two_pair_fraction = sum(row["character_class"] == "two_pair_mixed" for row in alternate) / 8
    min_gap, max_gap = min(gaps), max(gaps)

    primary_model = "8-state" if (
        max_gap <= LOW_GAP_LIMIT_MEV
        and tracked_single_fraction == 1.0
        and alternate_two_pair_fraction >= 0.75
    ) else "4-state"

    recommendation = (
        "Use the 8-state manifold as the primary low-energy excitonic model because both S1 and S2 remain within 120 meV and have systematically distinct NTO composition. Retain the 4-state tracked-bright manifold as a reduced control model."
        if primary_model == "8-state"
        else
        "Use the 4-state tracked-bright manifold as the primary model and retain the 8-state construction as a sensitivity test."
    )

    assessment_rows = [
        {"metric": "n_jobs", "value": 8, "interpretation": "Representative NTO calculations analyzed"},
        {"metric": "n_states", "value": 16, "interpretation": "S1 and S2 for all representative jobs"},
        {"metric": "tracked_equals_brightest_fraction", "value": 1.0, "interpretation": "Tracked root is the brighter S1/S2 root in every case"},
        {"metric": "tracked_single_pair_fraction", "value": tracked_single_fraction, "interpretation": "Tracked roots classified as single-pair dominated"},
        {"metric": "alternate_two_pair_fraction", "value": alternate_two_pair_fraction, "interpretation": "Alternate roots classified as two-pair mixed"},
        {"metric": "minimum_S1_S2_gap_meV", "value": min_gap, "interpretation": "Smallest representative low-state separation"},
        {"metric": "maximum_S1_S2_gap_meV", "value": max_gap, "interpretation": "Largest representative low-state separation"},
        {"metric": "mean_tracked_n1", "value": mean(float(row["n1"]) for row in tracked), "interpretation": "Mean dominant NTO occupation of tracked roots"},
        {"metric": "mean_alternate_n1", "value": mean(float(row["n1"]) for row in alternate), "interpretation": "Mean dominant NTO occupation of alternate roots"},
        {"metric": "primary_model_recommendation", "value": primary_model, "interpretation": recommendation},
    ]
    write_csv(MODEL_ASSESSMENT_CSV, assessment_rows, ["metric", "value", "interpretation"])

    with REPORT_MD.open("w", encoding="utf-8") as handle:
        handle.write("# Day019 NTO analysis and model-space assessment\n\n")
        handle.write("## Validation\n\n")
        handle.write("- Representative jobs analyzed: 8/8\n")
        handle.write("- NTO states analyzed: 16/16\n")
        handle.write("- Tracked root equals brightest S1/S2 root: 8/8 jobs\n")
        handle.write("- Every output passed normal termination, SCF, and TDDFT checks.\n\n")
        handle.write("## State-resolved metrics\n\n")
        handle.write("| Type | Frame | Site | Root | Tracked | Energy (eV) | fosc | n1 | n2 | n1+n2 | PR | Character |\n")
        handle.write("|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|\n")
        for row in sorted(state_rows, key=lambda x: (str(x["calculation_type"]), int(x["frame"]), str(x["cluster"]), int(x["root"]))):
            handle.write(
                f"| {row['calculation_type']} | {row['frame']} | {row['cluster']} | S{row['root']} | "
                f"{row['is_tracked_root']} | {f(row['energy_eV'],3)} | {f(row['fosc'],6)} | "
                f"{f(row['n1'],6)} | {f(row['n2'],6)} | {f(row['n1_plus_n2'],6)} | "
                f"{f(row['participation_ratio_printed'],4)} | {row['character_class']} |\n"
            )
        handle.write("\n## Aggregate findings\n\n")
        handle.write(f"- Representative S1-S2 gaps span {min_gap:.1f}-{max_gap:.1f} meV.\n")
        handle.write(f"- Tracked roots are single-pair dominated in {tracked_single_fraction*100:.1f}% of cases.\n")
        handle.write(f"- Alternate roots are two-pair mixed in {alternate_two_pair_fraction*100:.1f}% of cases.\n")
        handle.write(f"- Mean tracked-root dominant occupation: {mean(float(row['n1']) for row in tracked):.6f} (SD {pstdev(float(row['n1']) for row in tracked):.6f}).\n")
        handle.write(f"- Mean alternate-root dominant occupation: {mean(float(row['n1']) for row in alternate):.6f} (SD {pstdev(float(row['n1']) for row in alternate):.6f}).\n\n")
        handle.write("## Physical interpretation\n\n")
        handle.write(
            "For PYR2-PYR4 vacuum references, S2 is the bright tracked root and is dominated by the 52a->53a NTO pair, while S1 contains two comparably weighted NTO pairs. For PYR5, the ordering is reversed: S1 is the bright single-pair-dominated state and S2 is the two-pair-mixed state. Representative embedded cases preserve this distinction while changing the degree of mixing.\n\n"
        )
        handle.write("## Preliminary model-space decision\n\n")
        handle.write(f"**Primary recommendation: {primary_model}.**\n\n")
        handle.write(recommendation + "\n\n")
        handle.write(
            "This decision uses energies, oscillator strengths, and NTO occupation spectra. Orbital-shape visualization remains the next validation step because matching orbital indices do not by themselves prove spatial equivalence across sites and environments.\n"
        )

    print("Day019 NTO analysis completed.")
    print("Jobs analyzed: 8/8")
    print("States analyzed: 16/16")
    print("Tracked root equals brightest S1/S2 root: 8/8")
    print(f"Representative S1-S2 gap range: {min_gap:.1f}-{max_gap:.1f} meV")
    print(f"Tracked single-pair fraction: {tracked_single_fraction:.3f}")
    print(f"Alternate two-pair fraction: {alternate_two_pair_fraction:.3f}")
    print(f"Primary model recommendation: {primary_model}")
    print(f"Wrote: {OUTPUT_ROOT.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
