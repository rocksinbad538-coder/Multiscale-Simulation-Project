#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATA_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction"
)

VACUUM_ROOT = (
    DATA_ROOT
    / "day019_vacuum_reference_inputs"
)

EMBEDDED_TRACKING_CSV = (
    DATA_ROOT
    / "state_identity_analysis"
    / "low_state_identity_tracking.csv"
)

OUTPUT_ROOT = (
    DATA_ROOT
    / "day019_vacuum_reference_analysis"
)

VACUUM_STATE_CSV = (
    OUTPUT_ROOT
    / "vacuum_reference_state_tracking.csv"
)

SOLVENT_SHIFT_CSV = (
    OUTPUT_ROOT
    / "solvent_induced_site_energy_shifts.csv"
)

SOLVENT_STATS_CSV = (
    OUTPUT_ROOT
    / "solvent_shift_statistics.csv"
)

DECOMPOSITION_CSV = (
    OUTPUT_ROOT
    / "pyr5_offset_decomposition.csv"
)

NTO_SELECTION_CSV = (
    OUTPUT_ROOT
    / "NTO_SELECTION_DAY019.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "VACUUM_AND_SOLVENT_SHIFT_ANALYSIS_DAY019.md"
)

CHROMOPHORES = ["PYR2", "PYR3", "PYR4", "PYR5"]
CANDIDATE_ROOTS = [1, 2]

DEFAULT_OCCUPIED_ORBITAL = 52
DEFAULT_VIRTUAL_ORBITAL = 53

NORMAL_MARKER = "ORCA TERMINATED NORMALLY"
SCF_MARKER = "SCF CONVERGED"
TDDFT_MARKER = "ORCA-CIS/TD-DFT FINISHED WITHOUT ERROR"

STATE_RE = re.compile(
    r"^\s*STATE\s+(\d+)\s*:\s*E\s*=\s*"
    r"([-+0-9.Ee]+)\s+au\s+"
    r"([-+0-9.Ee]+)\s+eV",
    flags=re.IGNORECASE,
)

TRANSITION_RE = re.compile(
    r"^\s*(\d+)\s*([abAB]?)\s*->\s*"
    r"(\d+)\s*([abAB]?)\s*:\s*"
    r"([-+0-9.Ee]+)"
)

ABSORPTION_RE = re.compile(
    r"^\s*\d+-\S+\s*->\s*(\d+)-\S+\s+"
    r"([-+0-9.Ee]+)\s+"
    r"([-+0-9.Ee]+)\s+"
    r"([-+0-9.Ee]+)\s+"
    r"([-+0-9.Ee]+)"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Parse Day019 vacuum-reference ORCA calculations, "
            "track the HOMO-to-LUMO-dominated state, quantify "
            "embedded-minus-vacuum solvent shifts, decompose the "
            "PYR5 offset, and prepare representative NTO cases."
        )
    )
    parser.add_argument(
        "--occupied-orbital",
        type=int,
        default=DEFAULT_OCCUPIED_ORBITAL,
        help="ORCA occupied-orbital index used for state tracking.",
    )
    parser.add_argument(
        "--virtual-orbital",
        type=int,
        default=DEFAULT_VIRTUAL_ORBITAL,
        help="ORCA virtual-orbital index used for state tracking.",
    )
    return parser.parse_args()


def parse_absorption_fosc(lines: list[str]) -> dict[int, float]:
    oscillator_strengths: dict[int, float] = {}
    in_electric_absorption = False

    for line in lines:
        if (
            "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC "
            "DIPOLE MOMENTS"
            in line
        ):
            in_electric_absorption = True
            continue

        if (
            in_electric_absorption
            and "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY"
            in line
        ):
            break

        if not in_electric_absorption:
            continue

        match = ABSORPTION_RE.search(line)
        if not match:
            continue

        root = int(match.group(1))
        oscillator_strength = float(match.group(5))
        oscillator_strengths[root] = oscillator_strength

    return oscillator_strengths


def parse_orca_output(
    path: Path,
    occupied_orbital: int,
    virtual_orbital: int,
) -> dict:
    if not path.is_file():
        raise RuntimeError(f"Missing ORCA output: {path}")

    text = path.read_text(errors="ignore")
    lines = text.splitlines()

    if NORMAL_MARKER not in text:
        raise RuntimeError(f"ORCA did not terminate normally: {path}")

    if SCF_MARKER not in text:
        raise RuntimeError(f"SCF convergence marker missing: {path}")

    if TDDFT_MARKER not in text:
        raise RuntimeError(f"TDDFT/TDA completion marker missing: {path}")

    states: dict[int, dict] = {}
    current_root: int | None = None

    for line in lines:
        state_match = STATE_RE.match(line)

        if state_match:
            root = int(state_match.group(1))
            current_root = root
            states[root] = {
                "energy_au": float(state_match.group(2)),
                "energy_eV": float(state_match.group(3)),
                "transitions": [],
            }
            continue

        transition_match = TRANSITION_RE.match(line)

        if transition_match and current_root is not None:
            source_orbital = int(transition_match.group(1))
            source_spin = transition_match.group(2).lower()
            target_orbital = int(transition_match.group(3))
            target_spin = transition_match.group(4).lower()
            weight = float(transition_match.group(5))

            states[current_root]["transitions"].append(
                {
                    "source_orbital": source_orbital,
                    "source_spin": source_spin,
                    "target_orbital": target_orbital,
                    "target_spin": target_spin,
                    "weight": weight,
                }
            )

    oscillator_strengths = parse_absorption_fosc(lines)

    missing_roots = [
        root
        for root in CANDIDATE_ROOTS
        if root not in states
    ]

    if missing_roots:
        raise RuntimeError(
            f"Missing candidate roots {missing_roots} in {path}"
        )

    missing_fosc = [
        root
        for root in CANDIDATE_ROOTS
        if root not in oscillator_strengths
    ]

    if missing_fosc:
        raise RuntimeError(
            f"Missing oscillator strengths for roots {missing_fosc} "
            f"in {path}. Parsed roots: "
            f"{sorted(oscillator_strengths)}"
        )

    for root in CANDIDATE_ROOTS:
        target_weight = sum(
            transition["weight"]
            for transition in states[root]["transitions"]
            if (
                transition["source_orbital"] == occupied_orbital
                and transition["target_orbital"] == virtual_orbital
            )
        )

        states[root]["fosc"] = oscillator_strengths[root]
        states[root]["target_weight"] = target_weight

    tracked_root = max(
        CANDIDATE_ROOTS,
        key=lambda root: states[root]["target_weight"],
    )

    alternate_root = next(
        root
        for root in CANDIDATE_ROOTS
        if root != tracked_root
    )

    tracked = states[tracked_root]
    alternate = states[alternate_root]

    return {
        "tracked_root": tracked_root,
        "tracked_energy_eV": tracked["energy_eV"],
        "tracked_fosc": tracked["fosc"],
        "tracked_HOMO_LUMO_weight": tracked["target_weight"],
        "alternate_root": alternate_root,
        "alternate_energy_eV": alternate["energy_eV"],
        "alternate_fosc": alternate["fosc"],
        "alternate_HOMO_LUMO_weight": alternate["target_weight"],
        "HOMO_LUMO_weight_contrast": (
            tracked["target_weight"]
            - alternate["target_weight"]
        ),
        "state_separation_meV": (
            abs(
                states[2]["energy_eV"]
                - states[1]["energy_eV"]
            )
            * 1000.0
        ),
        "S1_energy_eV": states[1]["energy_eV"],
        "S1_fosc": states[1]["fosc"],
        "S1_HOMO_LUMO_weight": states[1]["target_weight"],
        "S2_energy_eV": states[2]["energy_eV"],
        "S2_fosc": states[2]["fosc"],
        "S2_HOMO_LUMO_weight": states[2]["target_weight"],
    }


def vacuum_output_path(chromophore: str) -> Path:
    job = f"frame000_{chromophore}_vacuum_reference"
    return VACUUM_ROOT / job / f"{job}.out"


def validate_embedded_tracking(
    dataframe: pd.DataFrame,
) -> pd.DataFrame:
    required_columns = {
        "frame",
        "cluster",
        "tracked_root",
        "tracked_energy_eV",
        "tracked_fosc",
        "tracked_HOMO_LUMO_weight",
        "state_separation_meV",
    }

    missing_columns = sorted(
        required_columns - set(dataframe.columns)
    )

    if missing_columns:
        raise RuntimeError(
            "Embedded tracking CSV is missing columns: "
            f"{missing_columns}"
        )

    work = dataframe.copy()

    work["frame"] = pd.to_numeric(
        work["frame"],
        errors="raise",
    ).astype(int)

    if len(work) != 84:
        raise RuntimeError(
            f"Expected 84 embedded tracking rows, found {len(work)}."
        )

    if work.duplicated(["frame", "cluster"]).any():
        raise RuntimeError(
            "Duplicate embedded frame/cluster rows detected."
        )

    expected_frames = list(range(21))
    observed_frames = sorted(
        work["frame"].unique().tolist()
    )

    if observed_frames != expected_frames:
        raise RuntimeError(
            f"Unexpected embedded frame set: {observed_frames}"
        )

    for frame in expected_frames:
        observed_clusters = set(
            work.loc[
                work["frame"] == frame,
                "cluster",
            ].tolist()
        )

        if observed_clusters != set(CHROMOPHORES):
            raise RuntimeError(
                f"Frame {frame:03d} has incorrect chromophore set."
            )

    return (
        work.sort_values(["frame", "cluster"])
        .reset_index(drop=True)
    )


def build_nto_selection(
    embedded: pd.DataFrame,
    solvent: pd.DataFrame,
    vacuum: pd.DataFrame,
) -> pd.DataFrame:
    candidates: list[dict] = []

    for row in vacuum.itertuples(index=False):
        candidates.append(
            {
                "calculation_type": "vacuum_reference",
                "frame": 0,
                "cluster": row.cluster,
                "job": (
                    f"frame000_{row.cluster}_vacuum_reference"
                ),
                "tracked_root": int(row.tracked_root),
                "reason": (
                    "vacuum_reference_for_each_chromophore"
                ),
                "metric_value": float(
                    row.tracked_HOMO_LUMO_weight
                ),
            }
        )

    minimum_character = embedded.loc[
        embedded[
            "tracked_HOMO_LUMO_weight"
        ].idxmin()
    ]

    candidates.append(
        {
            "calculation_type": "embedded",
            "frame": int(minimum_character["frame"]),
            "cluster": minimum_character["cluster"],
            "job": (
                f"frame{int(minimum_character['frame']):03d}_"
                f"{minimum_character['cluster']}_embedding"
            ),
            "tracked_root": int(
                minimum_character["tracked_root"]
            ),
            "reason": "minimum_tracked_character_weight",
            "metric_value": float(
                minimum_character[
                    "tracked_HOMO_LUMO_weight"
                ]
            ),
        }
    )

    minimum_gap = embedded.loc[
        embedded["state_separation_meV"].idxmin()
    ]

    candidates.append(
        {
            "calculation_type": "embedded",
            "frame": int(minimum_gap["frame"]),
            "cluster": minimum_gap["cluster"],
            "job": (
                f"frame{int(minimum_gap['frame']):03d}_"
                f"{minimum_gap['cluster']}_embedding"
            ),
            "tracked_root": int(
                minimum_gap["tracked_root"]
            ),
            "reason": "minimum_S1_S2_separation_meV",
            "metric_value": float(
                minimum_gap["state_separation_meV"]
            ),
        }
    )

    pyr5 = solvent.loc[
        solvent["cluster"] == "PYR5"
    ].copy()

    minimum_shift = pyr5.loc[
        pyr5["solvent_shift_meV"].idxmin()
    ]

    candidates.append(
        {
            "calculation_type": "embedded",
            "frame": int(minimum_shift["frame"]),
            "cluster": "PYR5",
            "job": (
                f"frame{int(minimum_shift['frame']):03d}_"
                "PYR5_embedding"
            ),
            "tracked_root": int(
                minimum_shift["embedded_tracked_root"]
            ),
            "reason": "minimum_PYR5_solvent_shift",
            "metric_value": float(
                minimum_shift["solvent_shift_meV"]
            ),
        }
    )

    maximum_shift = pyr5.loc[
        pyr5["solvent_shift_meV"].idxmax()
    ]

    candidates.append(
        {
            "calculation_type": "embedded",
            "frame": int(maximum_shift["frame"]),
            "cluster": "PYR5",
            "job": (
                f"frame{int(maximum_shift['frame']):03d}_"
                "PYR5_embedding"
            ),
            "tracked_root": int(
                maximum_shift["embedded_tracked_root"]
            ),
            "reason": "maximum_PYR5_solvent_shift",
            "metric_value": float(
                maximum_shift["solvent_shift_meV"]
            ),
        }
    )

    result = pd.DataFrame(candidates)

    result = (
        result.groupby(
            [
                "calculation_type",
                "frame",
                "cluster",
                "job",
                "tracked_root",
            ],
            as_index=False,
        )
        .agg(
            reason=(
                "reason",
                lambda values: ";".join(
                    sorted(set(values))
                ),
            ),
            metric_value=("metric_value", "first"),
        )
        .sort_values(
            ["calculation_type", "frame", "cluster"]
        )
        .reset_index(drop=True)
    )

    return result


def main() -> None:
    args = parse_args()

    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    vacuum_rows: list[dict] = []

    for chromophore in CHROMOPHORES:
        output_path = vacuum_output_path(chromophore)

        parsed = parse_orca_output(
            output_path,
            occupied_orbital=args.occupied_orbital,
            virtual_orbital=args.virtual_orbital,
        )

        vacuum_rows.append(
            {
                "cluster": chromophore,
                "output": str(
                    output_path.relative_to(PROJECT_ROOT)
                ),
                **parsed,
            }
        )

    vacuum = pd.DataFrame(vacuum_rows)

    vacuum.to_csv(
        VACUUM_STATE_CSV,
        index=False,
    )

    if not EMBEDDED_TRACKING_CSV.is_file():
        raise SystemExit(
            f"Missing embedded tracking CSV: "
            f"{EMBEDDED_TRACKING_CSV}"
        )

    embedded = validate_embedded_tracking(
        pd.read_csv(EMBEDDED_TRACKING_CSV)
    )

    vacuum_lookup = vacuum.set_index("cluster")

    solvent_rows: list[dict] = []

    for row in embedded.itertuples(index=False):
        vacuum_row = vacuum_lookup.loc[row.cluster]

        solvent_shift_eV = (
            float(row.tracked_energy_eV)
            - float(
                vacuum_row["tracked_energy_eV"]
            )
        )

        solvent_rows.append(
            {
                "frame": int(row.frame),
                "time_ps": float(row.frame) * 5.0,
                "cluster": row.cluster,
                "embedded_tracked_root": int(
                    row.tracked_root
                ),
                "embedded_tracked_energy_eV": float(
                    row.tracked_energy_eV
                ),
                "embedded_tracked_fosc": float(
                    row.tracked_fosc
                ),
                "embedded_HOMO_LUMO_weight": float(
                    row.tracked_HOMO_LUMO_weight
                ),
                "vacuum_tracked_root": int(
                    vacuum_row["tracked_root"]
                ),
                "vacuum_tracked_energy_eV": float(
                    vacuum_row["tracked_energy_eV"]
                ),
                "vacuum_tracked_fosc": float(
                    vacuum_row["tracked_fosc"]
                ),
                "vacuum_HOMO_LUMO_weight": float(
                    vacuum_row[
                        "tracked_HOMO_LUMO_weight"
                    ]
                ),
                "solvent_shift_eV": solvent_shift_eV,
                "solvent_shift_meV": (
                    1000.0 * solvent_shift_eV
                ),
            }
        )

    solvent = pd.DataFrame(solvent_rows)

    solvent.to_csv(
        SOLVENT_SHIFT_CSV,
        index=False,
    )

    solvent_statistics = (
        solvent.groupby("cluster")
        .agg(
            n_frames=("frame", "count"),
            vacuum_energy_eV=(
                "vacuum_tracked_energy_eV",
                "first",
            ),
            mean_embedded_energy_eV=(
                "embedded_tracked_energy_eV",
                "mean",
            ),
            mean_solvent_shift_meV=(
                "solvent_shift_meV",
                "mean",
            ),
            std_solvent_shift_meV=(
                "solvent_shift_meV",
                "std",
            ),
            min_solvent_shift_meV=(
                "solvent_shift_meV",
                "min",
            ),
            max_solvent_shift_meV=(
                "solvent_shift_meV",
                "max",
            ),
        )
        .reset_index()
    )

    solvent_statistics.to_csv(
        SOLVENT_STATS_CSV,
        index=False,
    )

    vacuum_reference_mean = float(
        vacuum.loc[
            vacuum["cluster"].isin(
                ["PYR2", "PYR3", "PYR4"]
            ),
            "tracked_energy_eV",
        ].mean()
    )

    vacuum_pyr5 = float(
        vacuum.loc[
            vacuum["cluster"] == "PYR5",
            "tracked_energy_eV",
        ].iloc[0]
    )

    vacuum_offset_meV = (
        1000.0
        * (
            vacuum_reference_mean
            - vacuum_pyr5
        )
    )

    embedded_means = (
        solvent.groupby("cluster")[
            "embedded_tracked_energy_eV"
        ]
        .mean()
    )

    embedded_reference_mean = float(
        embedded_means.loc[
            ["PYR2", "PYR3", "PYR4"]
        ].mean()
    )

    embedded_pyr5 = float(
        embedded_means.loc["PYR5"]
    )

    embedded_offset_meV = (
        1000.0
        * (
            embedded_reference_mean
            - embedded_pyr5
        )
    )

    mean_solvent_shifts = (
        solvent.groupby("cluster")[
            "solvent_shift_meV"
        ]
        .mean()
    )

    solvent_offset_contribution_meV = float(
        mean_solvent_shifts.loc[
            ["PYR2", "PYR3", "PYR4"]
        ].mean()
        - mean_solvent_shifts.loc["PYR5"]
    )

    closure_error_meV = (
        embedded_offset_meV
        - vacuum_offset_meV
        - solvent_offset_contribution_meV
    )

    decomposition = pd.DataFrame(
        [
            {
                "quantity": (
                    "vacuum_geometry_offset_"
                    "PYR5_vs_PYR2_PYR4"
                ),
                "value_meV": vacuum_offset_meV,
            },
            {
                "quantity": (
                    "mean_solvent_contribution_"
                    "to_PYR5_offset"
                ),
                "value_meV": (
                    solvent_offset_contribution_meV
                ),
            },
            {
                "quantity": (
                    "mean_embedded_PYR5_offset"
                ),
                "value_meV": embedded_offset_meV,
            },
            {
                "quantity": (
                    "decomposition_closure_error"
                ),
                "value_meV": closure_error_meV,
            },
        ]
    )

    decomposition.to_csv(
        DECOMPOSITION_CSV,
        index=False,
    )

    nto_selection = build_nto_selection(
        embedded=embedded,
        solvent=solvent,
        vacuum=vacuum,
    )

    nto_selection.to_csv(
        NTO_SELECTION_CSV,
        index=False,
    )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day019 vacuum-reference and "
            "solvent-shift analysis\n\n"
        )

        handle.write(
            "## Vacuum state tracking\n\n"
        )

        handle.write(
            vacuum[
                [
                    "cluster",
                    "tracked_root",
                    "tracked_energy_eV",
                    "tracked_fosc",
                    "tracked_HOMO_LUMO_weight",
                    "alternate_root",
                    "alternate_energy_eV",
                    "alternate_HOMO_LUMO_weight",
                    "state_separation_meV",
                ]
            ].to_string(index=False)
        )

        handle.write("\n\n")

        handle.write(
            "## Solvent-shift statistics\n\n"
        )

        handle.write(
            solvent_statistics.to_string(index=False)
        )

        handle.write("\n\n")

        handle.write(
            "## PYR5 offset decomposition\n\n"
        )

        handle.write(
            decomposition.to_string(index=False)
        )

        handle.write("\n\n")

        handle.write(
            "The decomposition is defined as:\n\n"
        )

        handle.write(
            "`embedded PYR5 offset = vacuum geometry "
            "offset + mean differential solvent shift`.\n\n"
        )

        handle.write(
            "The vacuum calculation removes the "
            "TIP4P/2005 point-charge environment while "
            "preserving the frozen chromophore geometry. "
            "The vacuum offset therefore measures the "
            "fixed-geometry baseline difference, subject "
            "to residual numerical orientation or "
            "integration-grid effects. The embedded-minus-"
            "vacuum difference measures the electrostatic "
            "effect of the water point charges for each "
            "fixed chromophore geometry.\n\n"
        )

        handle.write(
            "## NTO cases selected\n\n"
        )

        handle.write(
            nto_selection.to_string(index=False)
        )

        handle.write("\n")

    print(
        "Day019 vacuum and solvent-shift "
        "analysis completed."
    )

    print(
        f"Vacuum calculations parsed: "
        f"{len(vacuum)}/4"
    )

    print(
        f"Embedded observations processed: "
        f"{len(solvent)}/84"
    )

    print("\nVacuum tracked states:")

    print(
        vacuum[
            [
                "cluster",
                "tracked_root",
                "tracked_energy_eV",
                "tracked_fosc",
                "tracked_HOMO_LUMO_weight",
                "state_separation_meV",
            ]
        ].to_string(index=False)
    )

    print("\nSolvent-shift statistics:")

    print(
        solvent_statistics.to_string(index=False)
    )

    print("\nPYR5 offset decomposition:")

    print(
        decomposition.to_string(index=False)
    )

    print(
        f"\nNTO cases selected: "
        f"{len(nto_selection)}"
    )

    print(
        f"Wrote: "
        f"{OUTPUT_ROOT.relative_to(PROJECT_ROOT)}"
    )

    if not math.isclose(
        closure_error_meV,
        0.0,
        abs_tol=1.0e-8,
    ):
        raise SystemExit(
            "Offset decomposition closure failed: "
            f"{closure_error_meV:.12e} meV"
        )

    print("Offset decomposition closure: PASS")


if __name__ == "__main__":
    main()
