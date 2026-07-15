#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path
from statistics import mean, pstdev


PROJECT_ROOT = Path(__file__).resolve().parents[2]

DEFAULT_SOURCE_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction"
)

NTO_ANALYSIS_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/day019_nto_analysis"
)

VACUUM_STATE_METRICS = NTO_ANALYSIS_ROOT / "nto_state_metrics.csv"

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_full_two_state_manifold_analysis"
)

SOURCE_SELECTION_CSV = OUTPUT_ROOT / "SOURCE_OUTPUT_SELECTION_DAY019.csv"
OBSERVATIONS_CSV = OUTPUT_ROOT / "two_state_observations.csv"
FRAME_SUMMARY_CSV = OUTPUT_ROOT / "two_state_frame_summary.csv"
SITE_STATISTICS_CSV = OUTPUT_ROOT / "two_state_site_statistics.csv"
CORRELATIONS_CSV = OUTPUT_ROOT / "two_state_energy_correlations.csv"
REPORT_MD = OUTPUT_ROOT / "FULL_TWO_STATE_MANIFOLD_DAY019.md"

NORMAL_MARKER = "ORCA TERMINATED NORMALLY"
SCF_MARKER = "SCF CONVERGED"
TDDFT_MARKER = "ORCA-CIS/TD-DFT FINISHED WITHOUT ERROR"

FRAME_SITE_RE = re.compile(
    r"frame(?P<frame>\d{3})[_-](?P<site>PYR[2-5])",
    re.IGNORECASE,
)

STATE_ENERGY_RE = re.compile(
    r"STATE\s+(?P<root>\d+)\s*:\s*"
    r"E\s*=\s*(?P<energy_au>[-+0-9.Ee]+)\s+au\s+"
    r"(?P<energy_eV>[-+0-9.Ee]+)\s+eV",
    re.IGNORECASE,
)

ABSORPTION_RE = re.compile(
    r"^\s*\d+-\S+\s*->\s*(\d+)-\S+\s+"
    r"([-+0-9.Ee]+)\s+"
    r"([-+0-9.Ee]+)\s+"
    r"([-+0-9.Ee]+)\s+"
    r"([-+0-9.Ee]+)"
)

EXPECTED_SITES = ("PYR2", "PYR3", "PYR4", "PYR5")
EXPECTED_FRAMES = tuple(range(21))
EXPECTED_BRIGHT_ROOT = {
    "PYR2": 2,
    "PYR3": 2,
    "PYR4": 2,
    "PYR5": 1,
}


def log(message: str = "") -> None:
    print(message, flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Analyze S1/S2 over all 84 embedded TDDFT calculations, "
            "assign bright-like and alternate-like local-state families, "
            "and build the full two-state-per-site time series."
        )
    )
    parser.add_argument(
        "--source-root",
        type=Path,
        default=DEFAULT_SOURCE_ROOT,
        help="Root below which embedded ORCA .out files are discovered.",
    )
    return parser.parse_args()


def candidate_score(path: Path) -> int:
    text = str(path).lower()
    score = 0

    if "embedding" in text:
        score += 100
    if "tddft" in text:
        score += 20
    if "day018" in text:
        score += 10
    if "day017" in text:
        score += 5

    return score


def excluded_candidate(path: Path) -> bool:
    text = str(path).lower()

    excluded_tokens = (
        "day019_nto_inputs",
        "day019_nto_analysis",
        "day019_nto_cube",
        "day019_full_two_state",
        "vacuum_reference",
        ".incomplete_",
        "/archive",
        "/archives",
    )

    return any(token in text for token in excluded_tokens)


def identify_frame_site(path: Path) -> tuple[int, str] | None:
    match = FRAME_SITE_RE.search(str(path))

    if not match:
        return None

    frame = int(match.group("frame"))
    site = match.group("site").upper()

    if frame not in EXPECTED_FRAMES:
        return None

    return frame, site


def basic_output_validity(path: Path) -> tuple[bool, str]:
    try:
        text = path.read_text(errors="ignore")
    except OSError as exc:
        return False, f"read_error:{exc}"

    missing: list[str] = []

    if NORMAL_MARKER not in text:
        missing.append("normal_termination")
    if SCF_MARKER not in text:
        missing.append("scf_convergence")
    if TDDFT_MARKER not in text:
        missing.append("tddft_completion")

    if missing:
        return False, "missing:" + ",".join(missing)

    return True, "valid"


def discover_outputs(
    source_root: Path,
) -> tuple[
    dict[tuple[int, str], Path],
    list[dict[str, object]],
]:
    if not source_root.is_dir():
        raise SystemExit(f"Source root not found: {source_root}")

    grouped: dict[tuple[int, str], list[Path]] = defaultdict(list)

    for path in source_root.rglob("*.out"):
        if excluded_candidate(path):
            continue

        key = identify_frame_site(path)
        if key is None:
            continue

        valid, _ = basic_output_validity(path)
        if valid:
            grouped[key].append(path.resolve())

    expected_keys = {
        (frame, site)
        for frame in EXPECTED_FRAMES
        for site in EXPECTED_SITES
    }

    missing = sorted(expected_keys - set(grouped))

    if missing:
        preview = ", ".join(
            f"frame{frame:03d}_{site}"
            for frame, site in missing[:20]
        )
        raise RuntimeError(
            f"Missing valid embedded outputs for {len(missing)} keys: "
            f"{preview}"
        )

    selected: dict[tuple[int, str], Path] = {}
    audit_rows: list[dict[str, object]] = []

    for key in sorted(expected_keys):
        candidates = grouped[key]
        ranked = sorted(
            candidates,
            key=lambda path: (
                candidate_score(path),
                -len(str(path)),
                str(path),
            ),
            reverse=True,
        )

        top_score = candidate_score(ranked[0])
        top = [
            path
            for path in ranked
            if candidate_score(path) == top_score
        ]

        if len(top) > 1:
            raise RuntimeError(
                f"Ambiguous source outputs for frame{key[0]:03d}_{key[1]} "
                f"at score {top_score}:\n"
                + "\n".join(str(path) for path in top)
            )

        chosen = ranked[0]
        selected[key] = chosen

        for path in ranked:
            audit_rows.append(
                {
                    "frame": key[0],
                    "site": key[1],
                    "selected": path == chosen,
                    "score": candidate_score(path),
                    "size_bytes": path.stat().st_size,
                    "path": str(path.relative_to(PROJECT_ROOT)),
                }
            )

    if len(selected) != 84:
        raise RuntimeError(
            f"Expected 84 selected outputs, found {len(selected)}."
        )

    return selected, audit_rows


def parse_absorption_fosc(lines: list[str]) -> dict[int, float]:
    values: dict[int, float] = {}
    in_electric_absorption = False

    for line in lines:
        if (
            "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC "
            "DIPOLE MOMENTS" in line
        ):
            in_electric_absorption = True
            continue

        if (
            in_electric_absorption
            and "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY"
            in line
        ):
            break

        if in_electric_absorption:
            match = ABSORPTION_RE.search(line)
            if match:
                root = int(match.group(1))
                values[root] = float(match.group(5))

    return values


def parse_s1_s2(path: Path) -> dict[int, dict[str, float]]:
    text = path.read_text(errors="ignore")
    lines = text.splitlines()

    valid, reason = basic_output_validity(path)
    if not valid:
        raise RuntimeError(f"Invalid ORCA output {path}: {reason}")

    energies: dict[int, dict[str, float]] = {}

    for match in STATE_ENERGY_RE.finditer(text):
        root = int(match.group("root"))

        if root not in (1, 2):
            continue

        energies[root] = {
            "energy_au": float(match.group("energy_au")),
            "energy_eV": float(match.group("energy_eV")),
        }

    fosc = parse_absorption_fosc(lines)

    result: dict[int, dict[str, float]] = {}

    for root in (1, 2):
        if root not in energies:
            raise RuntimeError(
                f"Missing S{root} energy in {path}"
            )
        if root not in fosc:
            raise RuntimeError(
                f"Missing S{root} oscillator strength in {path}"
            )

        result[root] = {
            **energies[root],
            "fosc": fosc[root],
        }

    return result


def read_vacuum_references() -> dict[tuple[str, int], dict[str, float]]:
    if not VACUUM_STATE_METRICS.is_file():
        raise RuntimeError(
            f"Missing vacuum state metrics: {VACUUM_STATE_METRICS}"
        )

    with VACUUM_STATE_METRICS.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(csv.DictReader(handle))

    vacuum: dict[tuple[str, int], dict[str, float]] = {}

    for row in rows:
        if row["calculation_type"] != "vacuum_reference":
            continue

        site = row["cluster"]
        root = int(row["root"])

        vacuum[(site, root)] = {
            "energy_eV": float(row["energy_eV"]),
            "fosc": float(row["fosc"]),
        }

    expected = {
        (site, root)
        for site in EXPECTED_SITES
        for root in (1, 2)
    }

    if set(vacuum) != expected:
        raise RuntimeError(
            "Vacuum-reference table does not contain all eight "
            "site/root combinations."
        )

    return vacuum


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
    fieldnames: list[str] | None = None,
) -> None:
    if not rows:
        raise RuntimeError(f"No rows to write: {path}")

    if fieldnames is None:
        fieldnames = list(rows[0].keys())

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(rows)


def pearson(values_a: list[float], values_b: list[float]) -> float:
    if len(values_a) != len(values_b):
        raise RuntimeError("Correlation vectors have unequal lengths.")

    mean_a = mean(values_a)
    mean_b = mean(values_b)

    centered_a = [value - mean_a for value in values_a]
    centered_b = [value - mean_b for value in values_b]

    denom_a = math.sqrt(sum(value * value for value in centered_a))
    denom_b = math.sqrt(sum(value * value for value in centered_b))

    if denom_a == 0.0 or denom_b == 0.0:
        return float("nan")

    return (
        sum(a * b for a, b in zip(centered_a, centered_b))
        / (denom_a * denom_b)
    )


def summarize_values(values: list[float]) -> dict[str, float]:
    return {
        "mean": mean(values),
        "sd_population": pstdev(values),
        "minimum": min(values),
        "maximum": max(values),
        "range": max(values) - min(values),
    }


def main() -> None:
    args = parse_args()
    source_root = args.source_root.expanduser().resolve()

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    log("Day019 full two-state embedded-manifold analysis")
    log(f"Source root: {source_root}")

    selected, source_audit = discover_outputs(source_root)
    write_csv(SOURCE_SELECTION_CSV, source_audit)

    log(f"Selected embedded outputs: {len(selected)}/84")

    vacuum = read_vacuum_references()

    observations: list[dict[str, object]] = []
    frame_rows: list[dict[str, object]] = []

    for index, ((frame, site), path) in enumerate(
        sorted(selected.items()),
        start=1,
    ):
        states = parse_s1_s2(path)

        root_1 = states[1]
        root_2 = states[2]

        bright_root = max(
            (1, 2),
            key=lambda root: float(states[root]["fosc"]),
        )
        alternate_root = 1 if bright_root == 2 else 2

        expected_bright_root = EXPECTED_BRIGHT_ROOT[site]
        mapping_preserved = bright_root == expected_bright_root

        bright = states[bright_root]
        alternate = states[alternate_root]

        gap_meV = (
            float(root_2["energy_eV"])
            - float(root_1["energy_eV"])
        ) * 1000.0

        fosc_sum = float(root_1["fosc"]) + float(root_2["fosc"])
        brightness_fraction = (
            float(bright["fosc"]) / fosc_sum
            if fosc_sum > 0.0
            else float("nan")
        )

        frame_rows.append(
            {
                "frame": frame,
                "site": site,
                "S1_energy_eV": root_1["energy_eV"],
                "S1_fosc": root_1["fosc"],
                "S2_energy_eV": root_2["energy_eV"],
                "S2_fosc": root_2["fosc"],
                "S1_S2_gap_meV": gap_meV,
                "bright_root": bright_root,
                "alternate_root": alternate_root,
                "expected_bright_root": expected_bright_root,
                "mapping_preserved": mapping_preserved,
                "bright_energy_eV": bright["energy_eV"],
                "bright_fosc": bright["fosc"],
                "alternate_energy_eV": alternate["energy_eV"],
                "alternate_fosc": alternate["fosc"],
                "bright_fraction_of_S1_S2_fosc": brightness_fraction,
                "source_output": str(path.relative_to(PROJECT_ROOT)),
            }
        )

        for root in (1, 2):
            family = (
                "bright_like"
                if root == bright_root
                else "alternate_like"
            )

            reference = vacuum[(site, root)]

            observations.append(
                {
                    "frame": frame,
                    "site": site,
                    "root": root,
                    "family": family,
                    "is_bright_root": root == bright_root,
                    "expected_bright_root": expected_bright_root,
                    "mapping_preserved": mapping_preserved,
                    "energy_au": states[root]["energy_au"],
                    "energy_eV": states[root]["energy_eV"],
                    "fosc": states[root]["fosc"],
                    "vacuum_energy_eV_same_root": reference["energy_eV"],
                    "solvent_shift_meV_same_root": (
                        float(states[root]["energy_eV"])
                        - float(reference["energy_eV"])
                    )
                    * 1000.0,
                    "source_output": str(path.relative_to(PROJECT_ROOT)),
                }
            )

        log(
            f"[{index:02d}/84] frame{frame:03d} {site} "
            f"bright=S{bright_root} "
            f"gap={gap_meV:.1f} meV "
            f"f_bright={float(bright['fosc']):.6f} "
            f"mapping={'PASS' if mapping_preserved else 'SWITCH'}"
        )

    if len(frame_rows) != 84:
        raise RuntimeError(
            f"Expected 84 frame rows, found {len(frame_rows)}."
        )

    if len(observations) != 168:
        raise RuntimeError(
            f"Expected 168 state observations, found {len(observations)}."
        )

    write_csv(FRAME_SUMMARY_CSV, frame_rows)
    write_csv(OBSERVATIONS_CSV, observations)

    mapping_switches = [
        row
        for row in frame_rows
        if not bool(row["mapping_preserved"])
    ]

    site_statistics: list[dict[str, object]] = []

    for site in EXPECTED_SITES:
        site_frames = [
            row for row in frame_rows
            if row["site"] == site
        ]

        if len(site_frames) != 21:
            raise RuntimeError(
                f"Expected 21 frames for {site}, found {len(site_frames)}."
            )

        for family, energy_key, fosc_key in (
            ("bright_like", "bright_energy_eV", "bright_fosc"),
            (
                "alternate_like",
                "alternate_energy_eV",
                "alternate_fosc",
            ),
        ):
            energies = [
                float(row[energy_key])
                for row in site_frames
            ]
            foscs = [
                float(row[fosc_key])
                for row in site_frames
            ]

            energy_stats = summarize_values(energies)
            fosc_stats = summarize_values(foscs)

            site_statistics.append(
                {
                    "site": site,
                    "family": family,
                    "n_frames": len(site_frames),
                    "root": (
                        EXPECTED_BRIGHT_ROOT[site]
                        if family == "bright_like"
                        else (
                            1
                            if EXPECTED_BRIGHT_ROOT[site] == 2
                            else 2
                        )
                    ),
                    "energy_mean_eV": energy_stats["mean"],
                    "energy_sd_eV": energy_stats["sd_population"],
                    "energy_minimum_eV": energy_stats["minimum"],
                    "energy_maximum_eV": energy_stats["maximum"],
                    "energy_range_meV": energy_stats["range"] * 1000.0,
                    "fosc_mean": fosc_stats["mean"],
                    "fosc_sd": fosc_stats["sd_population"],
                    "fosc_minimum": fosc_stats["minimum"],
                    "fosc_maximum": fosc_stats["maximum"],
                }
            )

        gaps = [
            float(row["S1_S2_gap_meV"])
            for row in site_frames
        ]
        gap_stats = summarize_values(gaps)

        site_statistics.append(
            {
                "site": site,
                "family": "S1_S2_gap",
                "n_frames": len(site_frames),
                "root": "",
                "energy_mean_eV": "",
                "energy_sd_eV": "",
                "energy_minimum_eV": "",
                "energy_maximum_eV": "",
                "energy_range_meV": "",
                "fosc_mean": gap_stats["mean"],
                "fosc_sd": gap_stats["sd_population"],
                "fosc_minimum": gap_stats["minimum"],
                "fosc_maximum": gap_stats["maximum"],
            }
        )

    write_csv(SITE_STATISTICS_CSV, site_statistics)

    frame_lookup = {
        (int(row["frame"]), str(row["site"])): row
        for row in frame_rows
    }

    correlations: list[dict[str, object]] = []

    for family, key in (
        ("bright_like", "bright_energy_eV"),
        ("alternate_like", "alternate_energy_eV"),
    ):
        for index_a, site_a in enumerate(EXPECTED_SITES):
            for site_b in EXPECTED_SITES[index_a:]:
                values_a = [
                    float(frame_lookup[(frame, site_a)][key])
                    for frame in EXPECTED_FRAMES
                ]
                values_b = [
                    float(frame_lookup[(frame, site_b)][key])
                    for frame in EXPECTED_FRAMES
                ]

                correlations.append(
                    {
                        "family": family,
                        "site_a": site_a,
                        "site_b": site_b,
                        "pearson_r": pearson(values_a, values_b),
                        "n_frames": 21,
                    }
                )

    for site in EXPECTED_SITES:
        bright_values = [
            float(frame_lookup[(frame, site)]["bright_energy_eV"])
            for frame in EXPECTED_FRAMES
        ]
        alternate_values = [
            float(frame_lookup[(frame, site)]["alternate_energy_eV"])
            for frame in EXPECTED_FRAMES
        ]

        correlations.append(
            {
                "family": "within_site_bright_vs_alternate",
                "site_a": site,
                "site_b": site,
                "pearson_r": pearson(
                    bright_values,
                    alternate_values,
                ),
                "n_frames": 21,
            }
        )

    write_csv(CORRELATIONS_CSV, correlations)

    all_gaps = [
        float(row["S1_S2_gap_meV"])
        for row in frame_rows
    ]

    bright_fractions = [
        float(row["bright_fraction_of_S1_S2_fosc"])
        for row in frame_rows
        if math.isfinite(
            float(row["bright_fraction_of_S1_S2_fosc"])
        )
    ]

    with REPORT_MD.open("w", encoding="utf-8") as handle:
        handle.write(
            "# Day019 full two-state embedded-manifold analysis\n\n"
        )

        handle.write("## Validation\n\n")
        handle.write("- Embedded outputs selected: 84/84\n")
        handle.write("- S1/S2 state observations parsed: 168/168\n")
        handle.write(
            f"- Bright-root mapping preserved: "
            f"{84 - len(mapping_switches)}/84\n"
        )
        handle.write(
            f"- Bright-root switches detected: {len(mapping_switches)}\n\n"
        )

        handle.write("## Electronic-character convention\n\n")
        handle.write(
            "| Site | Bright-like local state | Alternate-like local state |\n"
        )
        handle.write("|---|---:|---:|\n")

        for site in EXPECTED_SITES:
            bright_root = EXPECTED_BRIGHT_ROOT[site]
            alternate_root = 1 if bright_root == 2 else 2
            handle.write(
                f"| {site} | S{bright_root} | S{alternate_root} |\n"
            )

        handle.write("\n")
        handle.write(
            "The family labels are defined by oscillator-strength ordering "
            "within S1/S2 and are consistent with the Day019 NTO occupation "
            "and cross-site subspace analyses.\n\n"
        )

        handle.write("## Aggregate results\n\n")
        handle.write(
            f"- Full S1-S2 gap range: "
            f"{min(all_gaps):.3f} to {max(all_gaps):.3f} meV\n"
        )
        handle.write(
            f"- Mean S1-S2 gap: {mean(all_gaps):.3f} meV\n"
        )
        handle.write(
            f"- Mean bright-state share of S1/S2 oscillator strength: "
            f"{mean(bright_fractions):.6f}\n"
        )
        handle.write(
            f"- Minimum bright-state share of S1/S2 oscillator strength: "
            f"{min(bright_fractions):.6f}\n\n"
        )

        handle.write("## Site-resolved statistics\n\n")
        handle.write(
            "| Site | Family | Root | Mean energy (eV) | SD (meV) | "
            "Energy range (meV) | Mean fosc |\n"
        )
        handle.write(
            "|---|---|---:|---:|---:|---:|---:|\n"
        )

        for row in site_statistics:
            if row["family"] == "S1_S2_gap":
                continue

            handle.write(
                f"| {row['site']} "
                f"| {row['family']} "
                f"| S{row['root']} "
                f"| {float(row['energy_mean_eV']):.6f} "
                f"| {float(row['energy_sd_eV']) * 1000.0:.3f} "
                f"| {float(row['energy_range_meV']):.3f} "
                f"| {float(row['fosc_mean']):.6f} |\n"
            )

        handle.write("\n## Model-space conclusion\n\n")

        if not mapping_switches:
            handle.write(
                "**The two local electronic families are preserved across "
                "all 84 embedded calculations.** The eight-state basis can "
                "therefore be indexed by electronic character rather than "
                "by a globally fixed root number. The four-state model is "
                "the bright-like subset of this basis.\n\n"
            )
        else:
            handle.write(
                "**Root switches were detected.** The eight-state basis "
                "must be tracked frame by frame by electronic character "
                "rather than root number.\n\n"
            )

        handle.write(
            "This analysis provides the complete diagonal-energy and "
            "oscillator-strength time series for the two-state-per-site "
            "manifold. Interstate and intersite couplings remain the next "
            "required ingredient before constructing the final excitonic "
            "Hamiltonian.\n"
        )

    log("")
    log("Day019 full two-state analysis completed.")
    log("Embedded outputs: 84/84")
    log("State observations: 168/168")
    log(
        f"Bright-root mapping preserved: "
        f"{84 - len(mapping_switches)}/84"
    )
    log(
        f"Full S1-S2 gap range: "
        f"{min(all_gaps):.3f}-{max(all_gaps):.3f} meV"
    )
    log(
        f"Minimum bright fosc share: "
        f"{min(bright_fractions):.6f}"
    )
    log(f"Wrote: {OUTPUT_ROOT.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
