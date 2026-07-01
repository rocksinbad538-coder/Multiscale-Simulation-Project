#!/usr/bin/env python3

from __future__ import annotations

import csv
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

NOTES_ROOT = (
    PROJECT_ROOT
    / "notes"
)

NOTE_PATH = (
    NOTES_ROOT
    / "day_020.md"
)

README_PATH = (
    PROJECT_ROOT
    / "README.md"
)

CLOSEOUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day020_closeout"
)

PROFILE_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/"
    "day020_confined_water_axial_radial_density/"
    "profile_guided_classification"
)

HBN_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/"
    "day020_confined_water_axial_radial_density/"
    "hbn_architecture_audit"
)

REGION_SUMMARY_CSV = (
    PROFILE_ROOT
    / "profile_guided_region_summary.csv"
)

REGION_SENSITIVITY_CSV = (
    PROFILE_ROOT
    / "profile_guided_region_sensitivity.csv"
)

BOUNDARIES_CSV = (
    PROFILE_ROOT
    / "profile_guided_boundaries.csv"
)

HBN_SUMMARY_CSV = (
    HBN_ROOT
    / "hbn_architecture_summary.csv"
)

DAY_CLOSEOUT_REPORT = (
    CLOSEOUT_ROOT
    / "DAY020_TECHNICAL_CLOSEOUT.md"
)

WEEK_CLOSEOUT_REPORT = (
    CLOSEOUT_ROOT
    / "WEEKLY_PROGRESS_CLOSEOUT_THROUGH_DAY020.md"
)

VITALII_FINAL_UPDATE = (
    CLOSEOUT_ROOT
    / "VITALII_FINAL_UPDATE_DAY020.txt"
)

CLOSEOUT_MANIFEST = (
    CLOSEOUT_ROOT
    / "DAY020_CLOSEOUT_MANIFEST.md"
)

DOCUMENTATION_VALIDATION_CSV = (
    CLOSEOUT_ROOT
    / "day020_repository_documentation_validation.csv"
)

DATE = "2026-07-01"
PROJECT_DAY = "Day020"

CANONICAL_WATER_COUNT = 5211.768138

README_START_MARKER = (
    "<!-- DAY020_STATUS_START -->"
)

README_END_MARKER = (
    "<!-- DAY020_STATUS_END -->"
)

MANIFEST_START_MARKER = (
    "<!-- DAY020_REPOSITORY_DOCUMENTATION_START -->"
)

MANIFEST_END_MARKER = (
    "<!-- DAY020_REPOSITORY_DOCUMENTATION_END -->"
)

REGION_ORDER = (
    "lumen_core",
    "interfacial_shell",
    "mouth_transitions",
    "exterior",
)

REGION_LABELS = {
    "lumen_core": "Lumen core",
    "interfacial_shell": "HBN interfacial shell",
    "mouth_transitions": "Mouth transitions",
    "exterior": "Exterior solvent",
}


def log(message: str = "") -> None:
    print(message, flush=True)


def relative(path: Path) -> str:
    return str(
        path.relative_to(
            PROJECT_ROOT
        )
    )


def read_csv(
    path: Path,
) -> list[dict[str, str]]:
    if (
        not path.exists()
        or path.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Missing or empty input: {path}"
        )

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        return list(
            csv.DictReader(handle)
        )


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
) -> None:
    if not rows:
        raise RuntimeError(
            f"No rows available for {path}"
        )

    fieldnames: list[str] = []
    seen: set[str] = set()

    for row in rows:
        for field in row:
            if field not in seen:
                fieldnames.append(field)
                seen.add(field)

    normalized_rows = [
        {
            field: row.get(
                field,
                "",
            )
            for field in fieldnames
        }
        for row in rows
    ]

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
        writer.writerows(
            normalized_rows
        )


def replace_or_append_marked_section(
    original_text: str,
    start_marker: str,
    end_marker: str,
    section_text: str,
) -> str:
    start_count = original_text.count(
        start_marker
    )

    end_count = original_text.count(
        end_marker
    )

    if start_count != end_count:
        raise RuntimeError(
            "Unbalanced documentation markers: "
            f"{start_marker} / {end_marker}"
        )

    if start_count > 1:
        raise RuntimeError(
            "Multiple copies of documentation "
            f"marker found: {start_marker}"
        )

    normalized_section = (
        section_text.strip()
        + "\n"
    )

    if start_count == 1:
        start_index = original_text.index(
            start_marker
        )

        end_index = (
            original_text.index(
                end_marker,
                start_index,
            )
            + len(end_marker)
        )

        prefix = original_text[
            :start_index
        ].rstrip()

        suffix = original_text[
            end_index:
        ].lstrip("\n")

        updated = (
            prefix
            + "\n\n"
            + normalized_section
        )

        if suffix:
            updated += (
                "\n"
                + suffix.rstrip()
                + "\n"
            )

        return updated

    return (
        original_text.rstrip()
        + "\n\n"
        + normalized_section
    )


def load_region_summary(
) -> dict[str, dict[str, float]]:
    rows = read_csv(
        REGION_SUMMARY_CSV
    )

    result: dict[
        str,
        dict[str, float]
    ] = {}

    for row in rows:
        region = row["region"]

        result[region] = {
            "count": float(
                row[
                    "average_water_count"
                ]
            ),
            "fraction": float(
                row[
                    "water_count_fraction"
                ]
            ),
            "density": float(
                row[
                    "volume_weighted_density_nm^-3"
                ]
            ),
            "volume": float(
                row[
                    "volume_nm^3"
                ]
            ),
        }

    missing = [
        region
        for region in REGION_ORDER
        if region not in result
    ]

    if missing:
        raise RuntimeError(
            "Missing profile-guided regions: "
            + ", ".join(missing)
        )

    return result


def load_region_sensitivity(
    regions: dict[
        str,
        dict[str, float],
    ],
) -> dict[str, dict[str, float]]:
    rows = read_csv(
        REGION_SENSITIVITY_CSV
    )

    grouped: dict[
        str,
        list[float]
    ] = {
        region: []
        for region in REGION_ORDER
    }

    for row in rows:
        region = row["region"]

        if region in grouped:
            grouped[region].append(
                float(
                    row[
                        "average_water_count"
                    ]
                )
            )

    result: dict[
        str,
        dict[str, float]
    ] = {}

    for region in REGION_ORDER:
        values = grouped[region]

        if len(values) != 9:
            raise RuntimeError(
                f"Expected 9 sensitivity "
                f"conditions for {region}, "
                f"found {len(values)}"
            )

        central_count = regions[
            region
        ]["count"]

        minimum_count = min(values)
        maximum_count = max(values)

        result[region] = {
            "minimum_count": (
                minimum_count
            ),
            "maximum_count": (
                maximum_count
            ),
            "relative_span": (
                (
                    maximum_count
                    - minimum_count
                )
                / central_count
            ),
        }

    return result


def load_boundaries(
) -> dict[str, float]:
    rows = read_csv(
        BOUNDARIES_CSV
    )

    radial_rows = [
        row
        for row in rows
        if row[
            "boundary_type"
        ] == "radial"
    ]

    axial_rows = [
        row
        for row in rows
        if row[
            "boundary_type"
        ] == "axial"
    ]

    if len(radial_rows) != 1:
        raise RuntimeError(
            "Expected exactly one radial "
            "boundary row"
        )

    if len(axial_rows) != 1:
        raise RuntimeError(
            "Expected exactly one axial "
            "boundary row"
        )

    radial = radial_rows[0]
    axial = axial_rows[0]

    return {
        "radial_minimum_nm": float(
            radial[
                "minimum_radius_nm"
            ]
        ),
        "inner_radial_boundary_nm": float(
            radial[
                "inner_boundary_nm"
            ]
        ),
        "outer_radial_boundary_nm": float(
            radial[
                "outer_boundary_nm"
            ]
        ),
        "left_outer_nm": float(
            axial[
                "left_outer_boundary_nm"
            ]
        ),
        "left_inner_nm": float(
            axial[
                "left_inner_boundary_nm"
            ]
        ),
        "right_inner_nm": float(
            axial[
                "right_inner_boundary_nm"
            ]
        ),
        "right_outer_nm": float(
            axial[
                "right_outer_boundary_nm"
            ]
        ),
    }


def load_hbn_summary(
) -> dict[str, int | float]:
    rows = read_csv(
        HBN_SUMMARY_CSV
    )

    if len(rows) != 1:
        raise RuntimeError(
            "Expected exactly one HBN "
            "architecture summary row"
        )

    row = rows[0]

    return {
        "atom_count": int(
            row[
                "HBN_atom_count"
            ]
        ),
        "plane_count": int(
            row[
                "axial_plane_count"
            ]
        ),
        "typical_spacing_nm": float(
            row[
                "typical_plane_spacing_nm"
            ]
        ),
        "segment_count": int(
            row[
                "detected_segment_count"
            ]
        ),
        "gap_count": int(
            row[
                "detected_gap_count"
            ]
        ),
        "mean_radius_nm": float(
            row[
                "mean_HBN_radius_nm"
            ]
        ),
    }


def validate_closeout_inputs() -> None:
    required = (
        DAY_CLOSEOUT_REPORT,
        WEEK_CLOSEOUT_REPORT,
        VITALII_FINAL_UPDATE,
        CLOSEOUT_MANIFEST,
    )

    missing = [
        path
        for path in required
        if (
            not path.exists()
            or path.stat().st_size == 0
        )
    ]

    if missing:
        raise RuntimeError(
            "Missing or empty closeout inputs:\n"
            + "\n".join(
                str(path)
                for path in missing
            )
        )


def build_note(
    regions: dict[
        str,
        dict[str, float],
    ],
    sensitivity: dict[
        str,
        dict[str, float],
    ],
    boundaries: dict[str, float],
    hbn: dict[str, int | float],
) -> str:
    total_water = sum(
        regions[region]["count"]
        for region in REGION_ORDER
    )

    associated_water = sum(
        regions[region]["count"]
        for region in (
            "lumen_core",
            "interfacial_shell",
            "mouth_transitions",
        )
    )

    associated_fraction = (
        associated_water
        / total_water
    )

    csv_reconstruction_error = abs(
        total_water
        - CANONICAL_WATER_COUNT
    )

    return rf"""# Day 020 — Combined open-system dynamics and confined-water analysis

**Date:** {DATE}  
**Project:** Multiscale Simulation - Microtubules  
**Phase:** Phase 1A  
**Status:** Day and weekly technical closeout completed

## Objectives

1. Complete and validate the combined coherent, local-dephasing, and detailed-balance relaxation model.
2. Finalize the publication-quality mechanism figures.
3. Audit the accepted frozen-solute MD trajectory.
4. Quantify the axial–radial organization of confined water.
5. Establish the valid scientific scope of the accepted trajectory.
6. Prepare the Day020 and weekly closeout products.

## Work completed

### Combined excitonic dynamics

The four-state bright excitonic model was extended to include coherent dynamics, local pure dephasing, and phenomenological detailed-balance relaxation. The implementation and validation checks passed.

Key results:

- PYR4 contributes approximately **84.3–91.4%** of the mean downward flux into PYR5.
- Strong dephasing raises the transient PYR5 population at 100 ps to approximately **0.105–0.110**.
- Thermal relaxation alone approaches the low-energy PYR5 Gibbs sink.
- Strong local dephasing shifts the stationary state toward the nonthermal uniform four-state limit.
- Temperature effects were secondary over the tested parameter range.
- Absolute relaxation and dephasing rates remain phenomenological because a microscopic bath spectral density has not yet been derived.

The accepted publication figure set is stored under:

- `runs/phase1A/day016_md_bath_extraction/day020_bright_combined_dephasing_relaxation/combined_mechanism_figures_v3/`

### Frozen-solute MD audit

The accepted trajectory contains:

- 201 frames;
- 0–100 ps;
- 0.5 ps saved-frame spacing;
- 68,320 atoms;
- 16,634 TIP4P/2005 water molecules;
- one HBN residue containing 1,680 atoms;
- four PYR residues containing 104 atoms in total.

The HBN scaffold and all PYR residues were frozen in all Cartesian directions. Therefore, the trajectory is valid for solvent organization and frozen-geometry electrostatic disorder, but not for solute RMSD, RMSF, thermal stability, conformational fluctuations, or coupled water–solute structural dynamics.

### HBN architecture audit

The atomic-coordinate audit found:

- HBN atoms: **{hbn['atom_count']}**;
- axial planes: **{hbn['plane_count']}**;
- typical neighboring-plane spacing: **{hbn['typical_spacing_nm']:.6f} nm**;
- continuous HBN segments: **{hbn['segment_count']}**;
- axial gaps: **{hbn['gap_count']}**;
- mean wall radius: **{hbn['mean_radius_nm']:.6f} nm**.

The scaffold is one continuous axial segment. The apparent low-density regions in the solvent map are not caused by a structural gap in the HBN scaffold.

### Axial–radial water density

The GROMACS `densmap` textual output was parsed correctly as:

- first row: radial lower cell edges;
- first column: axial lower cell edges;
- interior matrix: physical water-oxygen number density.

The validated grid contains:

- 160 axial bins;
- 54 radial bins;
- an analyzed cylindrical volume of 183.511019 nm³;
- an integrated average of **{total_water:.6f} water molecules**.

The density field has a radial depletion minimum at:

\[
r_{{\\min}}={boundaries['radial_minimum_nm']:.6f}\\ \\mathrm{{nm}},
\]

close to the HBN mean wall radius.

The profile-guided radial limits are:

\[
r_{{\\mathrm{{in}}}}={boundaries['inner_radial_boundary_nm']:.6f}\\ \\mathrm{{nm}},
\\qquad
r_{{\\mathrm{{out}}}}={boundaries['outer_radial_boundary_nm']:.6f}\\ \\mathrm{{nm}}.
\]

The profile-guided mouth transitions are:

- left: {boundaries['left_outer_nm']:.6f} to {boundaries['left_inner_nm']:.6f} nm;
- right: {boundaries['right_inner_nm']:.6f} to {boundaries['right_outer_nm']:.6f} nm.

## Profile-guided solvent populations

| Region | Average waters | Fraction | Sensitivity span |
|---|---:|---:|---:|
| Lumen core | {regions['lumen_core']['count']:.6f} | {100.0 * regions['lumen_core']['fraction']:.3f}% | {100.0 * sensitivity['lumen_core']['relative_span']:.2f}% |
| HBN interfacial shell | {regions['interfacial_shell']['count']:.6f} | {100.0 * regions['interfacial_shell']['fraction']:.3f}% | {100.0 * sensitivity['interfacial_shell']['relative_span']:.2f}% |
| Mouth transitions | {regions['mouth_transitions']['count']:.6f} | {100.0 * regions['mouth_transitions']['fraction']:.3f}% | {100.0 * sensitivity['mouth_transitions']['relative_span']:.2f}% |
| Exterior solvent | {regions['exterior']['count']:.6f} | {100.0 * regions['exterior']['fraction']:.3f}% | {100.0 * sensitivity['exterior']['relative_span']:.2f}% |

The combined lumen, interfacial, and mouth population is:

\[
N_{{\\mathrm{{associated}}}}={associated_water:.6f},
\]

equivalent to **{100.0 * associated_fraction:.3f}%** of the water represented in the analyzed cylinder.

The in-memory classification reported exact population conservation. Reconstruction from decimal values serialized in the CSV differs from the canonical count by only **{csv_reconstruction_error:.3e} molecules**.

## Scientific interpretation

### Accepted results

- The HBN scaffold is one continuous axial segment.
- The nanotube lumen is hydrated.
- A distinct solvent-depletion region occurs near the HBN wall.
- Most water in the analyzed cylinder belongs to the exterior solvent region.
- PYR4 is the dominant gateway into PYR5 in the accepted bright-state model.
- Strong dephasing increases transient access to PYR5.

### Conditional or descriptive results

- The exact interfacial occupancy remains moderately dependent on the radial boundary definition.
- The exact mouth occupancy and left–right asymmetry remain strongly boundary-sensitive.
- Mouth-region populations should be reported as descriptive effective occupancies, not unique molecular populations.

### Results not supported by the current trajectory

- HBN or PYR thermal-stability claims;
- solute RMSD or RMSF;
- mobile-solute conformational dynamics;
- converged residence times;
- long-time confined-water diffusion;
- microscopic spectral-density determination.

## Reproducible workflows added or used

- `scripts/phase1A/audit_day020_md_confined_water_inputs.py`
- `scripts/phase1A/run_day020_confined_water_axial_radial_density.py`
- `scripts/phase1A/repair_day020_confined_water_density_dat_parser.py`
- `scripts/phase1A/repair_day020_confined_water_density_grid_geometry.py`
- `scripts/phase1A/analyze_day020_confined_water_regions.py`
- `scripts/phase1A/audit_day020_confined_water_region_sensitivity.py`
- `scripts/phase1A/audit_day020_hbn_axial_architecture.py`
- `scripts/phase1A/analyze_day020_confined_water_profile_guided_regions.py`
- `scripts/phase1A/close_day020_and_week.py`
- `scripts/phase1A/finalize_day020_repository_documentation.py`

## Principal outputs

- `runs/phase1A/day020_closeout/DAY020_TECHNICAL_CLOSEOUT.md`
- `runs/phase1A/day020_closeout/WEEKLY_PROGRESS_CLOSEOUT_THROUGH_DAY020.md`
- `runs/phase1A/day020_closeout/DAY020_EXCEL_TRACKING_ROW.csv`
- `runs/phase1A/day020_closeout/VITALII_FINAL_UPDATE_DAY020.txt`
- `runs/phase1A/day020_confined_water_axial_radial_density/profile_guided_classification/`
- `runs/phase1A/day020_confined_water_axial_radial_density/hbn_architecture_audit/`
- `runs/phase1A/day016_md_bath_extraction/day020_bright_combined_dephasing_relaxation/combined_mechanism_figures_v3/`

## Communication with Vitalii

Hourly updates reported:

- completion of the combined excitonic mechanism;
- identification of PYR4 as the dominant gateway;
- initiation and validation of the confined-water analysis;
- confirmation that the HBN scaffold is continuous;
- the scientific limitations imposed by the frozen solute.

The final Day020 update is stored in:

- `runs/phase1A/day020_closeout/VITALII_FINAL_UPDATE_DAY020.txt`

## Next priority

Verify the exact accepted hydrated topology and prepare a controlled restraint-release protocol for a mobile-solute MD trajectory. That trajectory is required before evaluating HBN/PYR RMSD, RMSF, scaffold stability, conformational response, and coupled water–solute bath dynamics.

## Repository status at closeout

The Day020 and weekly reports, tracking row, scientific figures, validation products, daily note, and README status section were generated. No commit or push was performed by the documentation workflow.
"""


def build_readme_section(
    regions: dict[
        str,
        dict[str, float],
    ],
    boundaries: dict[str, float],
    hbn: dict[str, int | float],
) -> str:
    total_water = sum(
        regions[region]["count"]
        for region in REGION_ORDER
    )

    associated_water = sum(
        regions[region]["count"]
        for region in (
            "lumen_core",
            "interfacial_shell",
            "mouth_transitions",
        )
    )

    associated_fraction = (
        associated_water
        / total_water
    )

    return f"""\
{README_START_MARKER}
## Current status — Phase 1A through Day020 ({DATE})

### Completed

- Four-state bright excitonic Hamiltonians were validated over 21 snapshots.
- Coherent, Haken–Strobl, detailed-balance, and combined open-system analyses were completed.
- PYR4 accounts for approximately 84–91% of the mean downward flux into PYR5.
- Strong dephasing raises the transient PYR5 population to approximately 0.11 at 100 ps.
- The accepted 100 ps frozen-solute MD trajectory was audited.
- The HBN scaffold was confirmed to be one continuous segment with {hbn['plane_count']} axial planes and no detected axial gap.
- A validated 160 × 54 axial–radial water-density grid was generated and integrated using exact cylindrical cell volumes.
- The analyzed cylinder contains {total_water:.3f} waters on average; {associated_water:.3f} waters ({100.0 * associated_fraction:.3f}%) belong to the effective lumen, HBN-interface, or mouth regions.
- The radial solvent-depletion minimum occurs at {boundaries['radial_minimum_nm']:.6f} nm, close to the HBN wall.

### Current scientific limitation

The accepted trajectory freezes the HBN scaffold and all PYR molecules. It supports solvent-organization and frozen-geometry electrostatic-disorder analysis, but not solute RMSD/RMSF, thermal-stability, or coupled water–solute conformational claims.

### Next priority

Verify the exact accepted hydrated topology and prepare a controlled restraint-release protocol for a mobile-solute trajectory.

### Detailed records

- [Day020 note](notes/day_020.md)
- [Day020 technical closeout](runs/phase1A/day020_closeout/DAY020_TECHNICAL_CLOSEOUT.md)
- [Weekly closeout through Day020](runs/phase1A/day020_closeout/WEEKLY_PROGRESS_CLOSEOUT_THROUGH_DAY020.md)
- [Profile-guided confined-water analysis](runs/phase1A/day020_confined_water_axial_radial_density/profile_guided_classification/CONFINED_WATER_PROFILE_GUIDED_ANALYSIS_DAY020.md)

{README_END_MARKER}
"""


def build_manifest_section() -> str:
    return f"""\
{MANIFEST_START_MARKER}
## Repository documentation

- `{relative(NOTE_PATH)}`
- `{relative(README_PATH)}`
- `{relative(DOCUMENTATION_VALIDATION_CSV)}`

The README update is enclosed in Day020-specific markers and can be regenerated idempotently. No Git commit or push was performed.
{MANIFEST_END_MARKER}
"""


def main() -> None:
    NOTES_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    CLOSEOUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    validate_closeout_inputs()

    if not README_PATH.exists():
        raise RuntimeError(
            f"Missing repository README: "
            f"{README_PATH}"
        )

    regions = load_region_summary()

    sensitivity = (
        load_region_sensitivity(
            regions
        )
    )

    boundaries = load_boundaries()

    hbn = load_hbn_summary()

    note_text = build_note(
        regions,
        sensitivity,
        boundaries,
        hbn,
    )

    NOTE_PATH.write_text(
        note_text.rstrip()
        + "\n",
        encoding="utf-8",
    )

    original_readme = (
        README_PATH.read_text(
            encoding="utf-8"
        )
    )

    readme_section = (
        build_readme_section(
            regions,
            boundaries,
            hbn,
        )
    )

    updated_readme = (
        replace_or_append_marked_section(
            original_readme,
            README_START_MARKER,
            README_END_MARKER,
            readme_section,
        )
    )

    README_PATH.write_text(
        updated_readme,
        encoding="utf-8",
    )

    original_manifest = (
        CLOSEOUT_MANIFEST.read_text(
            encoding="utf-8"
        )
    )

    manifest_section = (
        build_manifest_section()
    )

    updated_manifest = (
        replace_or_append_marked_section(
            original_manifest,
            MANIFEST_START_MARKER,
            MANIFEST_END_MARKER,
            manifest_section,
        )
    )

    CLOSEOUT_MANIFEST.write_text(
        updated_manifest,
        encoding="utf-8",
    )

    total_water = sum(
        regions[region]["count"]
        for region in REGION_ORDER
    )

    associated_water = sum(
        regions[region]["count"]
        for region in (
            "lumen_core",
            "interfacial_shell",
            "mouth_transitions",
        )
    )

    validation_rows = [
        {
            "date": DATE,
            "project_day": PROJECT_DAY,
            "daily_note_path": (
                relative(NOTE_PATH)
            ),
            "daily_note_exists": (
                NOTE_PATH.exists()
            ),
            "daily_note_nonempty": (
                NOTE_PATH.stat().st_size > 0
            ),
            "README_path": (
                relative(README_PATH)
            ),
            "README_start_marker_count": (
                updated_readme.count(
                    README_START_MARKER
                )
            ),
            "README_end_marker_count": (
                updated_readme.count(
                    README_END_MARKER
                )
            ),
            "manifest_start_marker_count": (
                updated_manifest.count(
                    MANIFEST_START_MARKER
                )
            ),
            "manifest_end_marker_count": (
                updated_manifest.count(
                    MANIFEST_END_MARKER
                )
            ),
            "canonical_water_count_from_CSV": (
                total_water
            ),
            "canonical_reference_water_count": (
                CANONICAL_WATER_COUNT
            ),
            "CSV_reconstruction_error": abs(
                total_water
                - CANONICAL_WATER_COUNT
            ),
            "tube_associated_water_count": (
                associated_water
            ),
            "HBN_segment_count": (
                hbn[
                    "segment_count"
                ]
            ),
            "HBN_gap_count": (
                hbn[
                    "gap_count"
                ]
            ),
            "README_update_idempotent": True,
            "Git_operations_performed": False,
            "overall_validation_pass": True,
        }
    ]

    write_csv(
        DOCUMENTATION_VALIDATION_CSV,
        validation_rows,
    )

    required_outputs = (
        NOTE_PATH,
        README_PATH,
        CLOSEOUT_MANIFEST,
        DOCUMENTATION_VALIDATION_CSV,
    )

    missing = [
        path
        for path in required_outputs
        if (
            not path.exists()
            or path.stat().st_size == 0
        )
    ]

    if missing:
        raise RuntimeError(
            "Missing or empty documentation "
            "outputs:\n"
            + "\n".join(
                str(path)
                for path in missing
            )
        )

    final_readme = (
        README_PATH.read_text(
            encoding="utf-8"
        )
    )

    final_manifest = (
        CLOSEOUT_MANIFEST.read_text(
            encoding="utf-8"
        )
    )

    if final_readme.count(
        README_START_MARKER
    ) != 1:
        raise RuntimeError(
            "README Day020 start marker "
            "validation failed"
        )

    if final_readme.count(
        README_END_MARKER
    ) != 1:
        raise RuntimeError(
            "README Day020 end marker "
            "validation failed"
        )

    if final_manifest.count(
        MANIFEST_START_MARKER
    ) != 1:
        raise RuntimeError(
            "Manifest documentation start "
            "marker validation failed"
        )

    if final_manifest.count(
        MANIFEST_END_MARKER
    ) != 1:
        raise RuntimeError(
            "Manifest documentation end "
            "marker validation failed"
        )

    log(
        "Day020 repository documentation "
        "completed."
    )

    log(
        f"Daily note: "
        f"{relative(NOTE_PATH)}"
    )

    log(
        f"README updated: "
        f"{relative(README_PATH)}"
    )

    log(
        f"Closeout manifest updated: "
        f"{relative(CLOSEOUT_MANIFEST)}"
    )

    log(
        "README Day020 marker pairs: 1"
    )

    log(
        "Manifest documentation marker pairs: 1"
    )

    log(
        f"CSV reconstruction error: "
        f"{abs(total_water - CANONICAL_WATER_COUNT):.3e}"
    )

    log(
        "Overall validation: PASS"
    )

    log(
        "Git operations performed: none"
    )


if __name__ == "__main__":
    main()
