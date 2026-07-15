#!/usr/bin/env python3

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

AUDIT = (
    ROOT
    / "scripts/phase1A/"
    "audit_day024_r2_parent_rim_and_chemical_constraints.py"
)

CLASSIFICATION = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "01_r2_parent_rim_chemical_audit/"
    "topology_terminal_classification/"
    "hbn_terminal_coordination_classification_summary.csv"
)

EXPECTED_CLASSIFICATION_DECISION = (
    "HBN_EXPLICIT_TOPOLOGY_AND_TERMINAL_COORDINATION_CLASSIFIED"
)


def read_single_csv_row(
    path: Path,
) -> dict[str, str]:
    if (
        not path.exists()
        or not path.is_file()
        or path.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Missing or empty required file: {path}"
        )

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(
            csv.DictReader(handle)
        )

    if len(rows) != 1:
        raise RuntimeError(
            f"Expected one row in {path}; found {len(rows)}"
        )

    return rows[0]


def require_int(
    row: dict[str, str],
    key: str,
    expected: int,
) -> None:
    observed = int(
        float(
            row[key]
        )
    )

    if observed != expected:
        raise RuntimeError(
            f"Unexpected {key}: {observed}/{expected}"
        )


def replace_once(
    text: str,
    old: str,
    new: str,
    label: str,
) -> str:
    count = text.count(old)

    if count != 1:
        raise RuntimeError(
            f"Repair target {label!r} was found "
            f"{count} times instead of exactly once."
        )

    return text.replace(
        old,
        new,
        1,
    )


def main() -> None:
    classification = read_single_csv_row(
        CLASSIFICATION
    )

    if (
        classification.get(
            "decision",
            "",
        )
        != EXPECTED_CLASSIFICATION_DECISION
    ):
        raise RuntimeError(
            "The explicit-topology terminal classification "
            "is not in the approved state."
        )

    require_int(
        classification,
        "HBN_topology_atoms",
        1680,
    )

    require_int(
        classification,
        "explicit_bonds",
        2460,
    )

    require_int(
        classification,
        "geometry_bonds",
        2460,
    )

    require_int(
        classification,
        "explicit_only_bonds",
        0,
    )

    require_int(
        classification,
        "geometry_only_bonds",
        0,
    )

    require_int(
        classification,
        "explicit_degree1_atoms",
        60,
    )

    require_int(
        classification,
        "explicit_degree2_atoms",
        0,
    )

    require_int(
        classification,
        "explicit_degree3_atoms",
        1620,
    )

    require_int(
        classification,
        "terminal_B_atoms_total",
        30,
    )

    require_int(
        classification,
        "terminal_N_atoms_total",
        30,
    )

    require_int(
        classification,
        "lower_terminal_atoms",
        30,
    )

    require_int(
        classification,
        "lower_terminal_B_atoms",
        30,
    )

    require_int(
        classification,
        "lower_terminal_N_atoms",
        0,
    )

    require_int(
        classification,
        "upper_terminal_atoms",
        30,
    )

    require_int(
        classification,
        "upper_terminal_B_atoms",
        0,
    )

    require_int(
        classification,
        "upper_terminal_N_atoms",
        30,
    )

    text = AUDIT.read_text(
        encoding="utf-8",
    )

    text = replace_once(
        text,
        '''STATIC_SUMMARY = (
    DAY023_ROOT
    / "13_r2_topology_static_scan"
    / "r2_topology_static_scan_summary.csv"
)

SYSTEM_GRO = (
''',
        '''STATIC_SUMMARY = (
    DAY023_ROOT
    / "13_r2_topology_static_scan"
    / "r2_topology_static_scan_summary.csv"
)

TERMINAL_CLASSIFICATION_SUMMARY = (
    OUTPUT_ROOT
    / "topology_terminal_classification"
    / "hbn_terminal_coordination_classification_summary.csv"
)

SYSTEM_GRO = (
''',
        "terminal-classification source path",
    )

    text = replace_once(
        text,
        '''EXPECTED_SELECTION_DECISION = (
    "R2_SELECTED_AS_PRIMARY_PARTIAL_CAP_SCREENING_ARCHITECTURE"
)

EXPECTED_ATOMS = 68332
''',
        '''EXPECTED_SELECTION_DECISION = (
    "R2_SELECTED_AS_PRIMARY_PARTIAL_CAP_SCREENING_ARCHITECTURE"
)

EXPECTED_TERMINAL_CLASSIFICATION_DECISION = (
    "HBN_EXPLICIT_TOPOLOGY_AND_TERMINAL_COORDINATION_CLASSIFIED"
)

EXPECTED_ATOMS = 68332
''',
        "classification decision constant",
    )

    text = replace_once(
        text,
        '''EXPECTED_EDGE_ATOMS = 120
EXPECTED_EDGE_ATOMS_PER_END = 60
EXPECTED_EDGE_B_PER_END = 30
EXPECTED_EDGE_N_PER_END = 30

EXPECTED_INTERIOR_ATOMS = (
    EXPECTED_HBN_ATOMS
    - EXPECTED_EDGE_ATOMS
)

EXPECTED_GEOMETRY_BONDS = (
    (
        EXPECTED_EDGE_ATOMS
        * 2
        + EXPECTED_INTERIOR_ATOMS
        * 3
    )
    // 2
)
''',
        '''EXPECTED_EDGE_ATOMS = 60
EXPECTED_EDGE_ATOMS_PER_END = 30

EXPECTED_LOWER_EDGE_B = 30
EXPECTED_LOWER_EDGE_N = 0

EXPECTED_UPPER_EDGE_B = 0
EXPECTED_UPPER_EDGE_N = 30

EXPECTED_INTERIOR_ATOMS = (
    EXPECTED_HBN_ATOMS
    - EXPECTED_EDGE_ATOMS
)

EXPECTED_GEOMETRY_BONDS = (
    (
        EXPECTED_EDGE_ATOMS
        * 1
        + EXPECTED_INTERIOR_ATOMS
        * 3
    )
    // 2
)

EXPECTED_TERMINAL_COORDINATION_DEFICIT_TO_DEGREE3 = (
    EXPECTED_EDGE_ATOMS
    * 2
)
''',
        "terminal coordination constants",
    )

    text = replace_once(
        text,
        '''MIN_EDGE_ALTERNATION_FRACTION = 0.95
''',
        '''MIN_TERMINAL_ELEMENT_PURITY_FRACTION = 0.999999
''',
        "terminal purity threshold",
    )

    text = replace_once(
        text,
        '''    alternation_fraction = float(
        np.mean(
            alternation
        )
    )

    atom_rows = []
''',
        '''    alternation_fraction = float(
        np.mean(
            alternation
        )
    )

    b_count = int(
        np.count_nonzero(
            elements[
                indices
            ]
            == "B"
        )
    )

    n_count = int(
        np.count_nonzero(
            elements[
                indices
            ]
            == "N"
        )
    )

    dominant_element = (
        "B"
        if b_count > n_count
        else (
            "N"
            if n_count > b_count
            else "BALANCED"
        )
    )

    element_purity_fraction = (
        max(
            b_count,
            n_count,
        )
        / len(indices)
    )

    atom_rows = []
''',
        "terminal composition metrics",
    )

    text = replace_once(
        text,
        '''        "B_count": int(
            np.count_nonzero(
                elements[
                    indices
                ]
                == "B"
            )
        ),
        "N_count": int(
            np.count_nonzero(
                elements[
                    indices
                ]
                == "N"
            )
        ),
''',
        '''        "B_count": b_count,
        "N_count": n_count,
''',
        "terminal B/N metric reuse",
    )

    text = replace_once(
        text,
        '''        "element_alternation_fraction": (
            alternation_fraction
        ),
''',
        '''        "element_alternation_fraction": (
            alternation_fraction
        ),
        "dominant_element": (
            dominant_element
        ),
        "element_purity_fraction": (
            element_purity_fraction
        ),
''',
        "terminal purity output",
    )

    text = replace_once(
        text,
        '''        STATIC_SUMMARY,
        SYSTEM_GRO,
        SYSTEM_TPR,
''',
        '''        STATIC_SUMMARY,
        TERMINAL_CLASSIFICATION_SUMMARY,
        SYSTEM_GRO,
        SYSTEM_TPR,
''',
        "required classification input",
    )

    text = replace_once(
        text,
        '''    static_summary = (
        read_single_csv_row(
            STATIC_SUMMARY
        )
    )

    if (
        selection.get(
''',
        '''    static_summary = (
        read_single_csv_row(
            STATIC_SUMMARY
        )
    )

    terminal_classification = (
        read_single_csv_row(
            TERMINAL_CLASSIFICATION_SUMMARY
        )
    )

    if (
        selection.get(
''',
        "classification CSV read",
    )

    text = replace_once(
        text,
        '''    edge_indices = np.flatnonzero(
        degrees == 2
    )

    interior_indices = np.flatnonzero(
        degrees == 3
    )

    anomalous_indices = np.flatnonzero(
        (degrees < 2)
        | (degrees > 3)
    )
''',
        '''    edge_indices = np.flatnonzero(
        degrees == 1
    )

    degree2_indices = np.flatnonzero(
        degrees == 2
    )

    interior_indices = np.flatnonzero(
        degrees == 3
    )

    anomalous_indices = np.flatnonzero(
        (degrees != 1)
        & (degrees != 3)
    )
''',
        "degree-1 terminal selection",
    )

    text = replace_once(
        text,
        '''        {
            "constraint": (
                "terminal_anchor_atoms_per_end"
            ),
            "value": (
                EXPECTED_EDGE_ATOMS_PER_END
            ),
            "unit": "atoms/end",
            "basis": (
                "Geometry-derived degree-2 rim"
            ),
            "status": "FIXED",
        },
        {
            "constraint": (
                "terminal_B_anchor_atoms_per_end"
            ),
            "value": (
                EXPECTED_EDGE_B_PER_END
            ),
            "unit": "atoms/end",
            "basis": (
                "Balanced armchair-like BN rim"
            ),
            "status": "FIXED",
        },
        {
            "constraint": (
                "terminal_N_anchor_atoms_per_end"
            ),
            "value": (
                EXPECTED_EDGE_N_PER_END
            ),
            "unit": "atoms/end",
            "basis": (
                "Balanced armchair-like BN rim"
            ),
            "status": "FIXED",
        },
''',
        '''        {
            "constraint": (
                "terminal_anchor_atoms_per_end"
            ),
            "value": (
                EXPECTED_EDGE_ATOMS_PER_END
            ),
            "unit": "atoms/end",
            "basis": (
                "Explicit topology and geometry: "
                "degree-1 terminal sites"
            ),
            "status": "FIXED",
        },
        {
            "constraint": (
                "lower_terminal_B_anchor_atoms"
            ),
            "value": (
                EXPECTED_LOWER_EDGE_B
            ),
            "unit": "atoms/end",
            "basis": (
                "Lower end is purely B-terminated"
            ),
            "status": "FIXED",
        },
        {
            "constraint": (
                "lower_terminal_N_anchor_atoms"
            ),
            "value": (
                EXPECTED_LOWER_EDGE_N
            ),
            "unit": "atoms/end",
            "basis": (
                "Lower end is purely B-terminated"
            ),
            "status": "FIXED",
        },
        {
            "constraint": (
                "upper_terminal_B_anchor_atoms"
            ),
            "value": (
                EXPECTED_UPPER_EDGE_B
            ),
            "unit": "atoms/end",
            "basis": (
                "Upper end is purely N-terminated"
            ),
            "status": "FIXED",
        },
        {
            "constraint": (
                "upper_terminal_N_anchor_atoms"
            ),
            "value": (
                EXPECTED_UPPER_EDGE_N
            ),
            "unit": "atoms/end",
            "basis": (
                "Upper end is purely N-terminated"
            ),
            "status": "FIXED",
        },
        {
            "constraint": (
                "parent_terminal_coordination_deficit_to_degree3"
            ),
            "value": (
                EXPECTED_TERMINAL_COORDINATION_DEFICIT_TO_DEGREE3
            ),
            "unit": "missing neighbor incidences",
            "basis": (
                "60 degree-1 termini relative to "
                "three-coordinate interior sites"
            ),
            "status": "REQUIRES_CHEMICAL_RESOLUTION",
        },
''',
        "polar terminal design constraints",
    )

    text = replace_once(
        text,
        '''                "Passivate the existing BNNT rim "
                "using only terminal H, OH, or NHx groups."
''',
        '''                "Attempt to passivate the existing "
                "degree-1 polar BN termini using only "
                "small terminal groups."
''',
        "C0 structural concept",
    )

    text = replace_once(
        text,
        '''                "BN edge passivation is chemically "
                "documented."
''',
        '''                "Small-group BN edge passivation is "
                "chemically plausible in general, but the "
                "accepted parent sites are degree 1 rather "
                "than conventional degree-2 edge sites."
''',
        "C0 chemical basis",
    )

    text = replace_once(
        text,
        '''                "Insufficient as a direct R2-cap "
                "replacement because the required "
                f"inward radial closure is "
                f"{required_radial_occlusion_nm:.6f} nm."
''',
        '''                "Insufficient as a direct R2-cap "
                "replacement because it neither reproduces "
                "the required inward radial closure of "
                f"{required_radial_occlusion_nm:.6f} nm "
                "nor by itself resolves the full parent-side "
                "coordination deficit."
''',
        "C0 geometric assessment",
    )

    text = replace_once(
        text,
        '''                "Leaves a pore much wider than the "
                "validated 0.84 nm aperture."
''',
        '''                "Leaves a pore much wider than the "
                "validated 0.84 nm aperture and may leave "
                "the degree-1 polar termini chemically "
                "under-resolved."
''',
        "C0 principal risk",
    )

    text = replace_once(
        text,
        '''                "Coaxial atomically thin BN annular "
                "nanoflake placed at the validated cap "
                "plane and attached to the terminal rim "
                "through a chemically defined linker "
                "network."
''',
        '''                "Coaxial atomically thin BN annular "
                "nanoflake placed at the validated cap "
                "plane and attached through end-specific "
                "junction networks: one for the B-terminated "
                "lower rim and one for the N-terminated "
                "upper rim."
''',
        "C1 polar junction concept",
    )

    text = replace_once(
        text,
        '''                "Outer-rim linker chemistry and junction "
                "strain must be resolved; a direct "
                "90-degree seamless sp2 seam is not "
                "assumed."
''',
        '''                "B-end and N-end linker chemistry, "
                "charge compensation, and junction strain "
                "must be resolved independently; a direct "
                "90-degree seamless sp2 seam is not assumed."
''',
        "C1 principal risk",
    )

    text = replace_once(
        text,
        '''        "R2_architecture_selection_is_valid": (
            selection.get(
                "decision",
                "",
            )
            == EXPECTED_SELECTION_DECISION
        ),
''',
        '''        "R2_architecture_selection_is_valid": (
            selection.get(
                "decision",
                "",
            )
            == EXPECTED_SELECTION_DECISION
        ),
        "explicit_topology_terminal_classification_is_valid": (
            terminal_classification.get(
                "decision",
                "",
            )
            == EXPECTED_TERMINAL_CLASSIFICATION_DECISION
        ),
''',
        "classification validity gate",
    )

    text = replace_once(
        text,
        '''        "geometry_has_120_degree2_edge_atoms": (
            len(
                edge_indices
            )
            == EXPECTED_EDGE_ATOMS
        ),
        "geometry_has_1560_degree3_interior_atoms": (
            len(
                interior_indices
            )
            == EXPECTED_INTERIOR_ATOMS
        ),
''',
        '''        "geometry_has_60_degree1_terminal_atoms": (
            len(
                edge_indices
            )
            == EXPECTED_EDGE_ATOMS
        ),
        "geometry_has_zero_degree2_atoms": (
            len(
                degree2_indices
            )
            == 0
        ),
        "geometry_has_1620_degree3_interior_atoms": (
            len(
                interior_indices
            )
            == EXPECTED_INTERIOR_ATOMS
        ),
''',
        "coordination gates",
    )

    text = replace_once(
        text,
        '''        "lower_end_has_30B_and_30N": (
            int(
                lower_metrics[
                    "B_count"
                ]
            )
            == EXPECTED_EDGE_B_PER_END
            and int(
                lower_metrics[
                    "N_count"
                ]
            )
            == EXPECTED_EDGE_N_PER_END
        ),
        "upper_end_has_30B_and_30N": (
            int(
                upper_metrics[
                    "B_count"
                ]
            )
            == EXPECTED_EDGE_B_PER_END
            and int(
                upper_metrics[
                    "N_count"
                ]
            )
            == EXPECTED_EDGE_N_PER_END
        ),
''',
        '''        "lower_end_is_30B_and_0N": (
            int(
                lower_metrics[
                    "B_count"
                ]
            )
            == EXPECTED_LOWER_EDGE_B
            and int(
                lower_metrics[
                    "N_count"
                ]
            )
            == EXPECTED_LOWER_EDGE_N
        ),
        "upper_end_is_0B_and_30N": (
            int(
                upper_metrics[
                    "B_count"
                ]
            )
            == EXPECTED_UPPER_EDGE_B
            and int(
                upper_metrics[
                    "N_count"
                ]
            )
            == EXPECTED_UPPER_EDGE_N
        ),
''',
        "polar end composition gates",
    )

    text = replace_once(
        text,
        '''        "lower_end_element_sequence_is_alternating": (
            float(
                lower_metrics[
                    "element_alternation_fraction"
                ]
            )
            >= MIN_EDGE_ALTERNATION_FRACTION
        ),
        "upper_end_element_sequence_is_alternating": (
            float(
                upper_metrics[
                    "element_alternation_fraction"
                ]
            )
            >= MIN_EDGE_ALTERNATION_FRACTION
        ),
''',
        '''        "lower_end_is_elementally_pure_B": (
            lower_metrics[
                "dominant_element"
            ]
            == "B"
            and float(
                lower_metrics[
                    "element_purity_fraction"
                ]
            )
            >= MIN_TERMINAL_ELEMENT_PURITY_FRACTION
        ),
        "upper_end_is_elementally_pure_N": (
            upper_metrics[
                "dominant_element"
            ]
            == "N"
            and float(
                upper_metrics[
                    "element_purity_fraction"
                ]
            )
            >= MIN_TERMINAL_ELEMENT_PURITY_FRACTION
        ),
        "terminal_to_terminal_bond_count_is_zero": (
            edge_edge_bonds == 0
        ),
        "terminal_to_interior_bond_count_is_60": (
            edge_interior_bonds
            == EXPECTED_EDGE_ATOMS
        ),
''',
        "terminal polarity and connectivity gates",
    )

    text = replace_once(
        text,
        '''        "degree2_terminal_atoms": (
            len(
                edge_indices
            )
        ),
        "degree3_interior_atoms": (
''',
        '''        "degree1_terminal_atoms": (
            len(
                edge_indices
            )
        ),
        "degree2_atoms": (
            len(
                degree2_indices
            )
        ),
        "degree3_interior_atoms": (
''',
        "summary coordination fields",
    )

    text = replace_once(
        text,
        '''        "lower_end_alternation_fraction": (
            lower_metrics[
                "element_alternation_fraction"
            ]
        ),
''',
        '''        "lower_end_dominant_element": (
            lower_metrics[
                "dominant_element"
            ]
        ),
        "lower_end_element_purity_fraction": (
            lower_metrics[
                "element_purity_fraction"
            ]
        ),
        "lower_end_alternation_fraction": (
            lower_metrics[
                "element_alternation_fraction"
            ]
        ),
''',
        "lower summary composition",
    )

    text = replace_once(
        text,
        '''        "upper_end_alternation_fraction": (
            upper_metrics[
                "element_alternation_fraction"
            ]
        ),
''',
        '''        "upper_end_dominant_element": (
            upper_metrics[
                "dominant_element"
            ]
        ),
        "upper_end_element_purity_fraction": (
            upper_metrics[
                "element_purity_fraction"
            ]
        ),
        "upper_end_alternation_fraction": (
            upper_metrics[
                "element_alternation_fraction"
            ]
        ),
        "terminal_polarity": (
            "LOWER_B_TERMINATED_UPPER_N_TERMINATED"
        ),
        "terminal_coordination_deficit_to_degree3": (
            EXPECTED_TERMINAL_COORDINATION_DEFICIT_TO_DEGREE3
        ),
''',
        "upper summary composition",
    )

    text = replace_once(
        text,
        '''        "R2_static_summary"
            ),
            "file": relative(
                STATIC_SUMMARY
            ),
            "sha256": sha256(
                STATIC_SUMMARY
            ),
        },
        {
            "role": (
                "R2_parent_coordinates"
''',
        '''        "R2_static_summary"
            ),
            "file": relative(
                STATIC_SUMMARY
            ),
            "sha256": sha256(
                STATIC_SUMMARY
            ),
        },
        {
            "role": (
                "HBN_terminal_topology_classification"
            ),
            "file": relative(
                TERMINAL_CLASSIFICATION_SUMMARY
            ),
            "sha256": sha256(
                TERMINAL_CLASSIFICATION_SUMMARY
            ),
        },
        {
            "role": (
                "R2_parent_coordinates"
''',
        "classification source manifest",
    )

    text = replace_once(
        text,
        '''The degree-2 atoms represent unsaturated terminal-rim sites in the
current unpassivated parent model. Frozen-coordinate stability does
not itself establish chemical stability of these sites.
''',
        '''The accepted topology and the independent geometric graph identify
60 degree-1 terminal sites and 1620 degree-3 interior sites. The lower
end contains 30 B-only terminal atoms, whereas the upper end contains
30 N-only terminal atoms. These are strongly undercoordinated polar
termini in the current parent model. Frozen-coordinate stability does
not establish their chemical stability.
''',
        "report terminal interpretation",
    )

    text = replace_once(
        text,
        '''- Element-alternation fraction:
  **{lower_metrics['element_alternation_fraction']:.6f}**
''',
        '''- Dominant element/purity:
  **{lower_metrics['dominant_element']}/
  {lower_metrics['element_purity_fraction']:.6f}**
- Element-alternation fraction:
  **{lower_metrics['element_alternation_fraction']:.6f}**
''',
        "lower report polarity",
    )

    text = replace_once(
        text,
        '''- Element-alternation fraction:
  **{upper_metrics['element_alternation_fraction']:.6f}**
''',
        '''- Dominant element/purity:
  **{upper_metrics['dominant_element']}/
  {upper_metrics['element_purity_fraction']:.6f}**
- Element-alternation fraction:
  **{upper_metrics['element_alternation_fraction']:.6f}**
''',
        "upper report polarity",
    )

    text = replace_once(
        text,
        '''Simple H/OH/NHx passivation of the existing terminal atoms cannot by
itself reproduce the validated aperture because the required inward
radial closure is approximately
**{required_radial_occlusion_nm:.6f} nm**.
''',
        '''Simple H/OH/NHx termination cannot be promoted directly as the R2
chemical replacement. It cannot reproduce the validated aperture, and
the parent atoms are degree-1 polar termini rather than conventional
degree-2 edge sites. The chemical design must therefore resolve an
aggregate parent-side coordination deficit of
**{EXPECTED_TERMINAL_COORDINATION_DEFICIT_TO_DEGREE3} missing neighbor
incidences**, while also providing an inward radial closure of
**{required_radial_occlusion_nm:.6f} nm**.
''',
        "report passivation conclusion",
    )

    text = replace_once(
        text,
        '''The leading geometric analogue is a separately defined annular
nanoflake or rigid macrocycle positioned at the validated cap plane
and attached through an explicit chemically valid junction. A direct
unstrained 90-degree seamless sp2 junction is not assumed.
''',
        '''The leading geometric analogue remains a separately defined annular
nanoflake or rigid macrocycle positioned at the validated cap plane.
However, the lower B-terminated rim and upper N-terminated rim require
distinct junction chemistry or an explicitly justified compensation
strategy. A symmetric identical linker assignment and a direct
unstrained 90-degree seamless sp2 junction are not assumed.
''',
        "report candidate interpretation",
    )

    text = replace_once(
        text,
        '''        "DEFINE_AND_RANK_R2_EXPLICIT_END_RIM_CHEMISTRY_CANDIDATES"
''',
        '''        "DEFINE_AND_RANK_R2_POLAR_END_SPECIFIC_VALENCE_COMPLETION_CANDIDATES"
''',
        "next-step decision",
    )

    text = replace_once(
        text,
        '''        "Degree-2 edge / degree-3 interior / "
        "anomalous atoms: "
''',
        '''        "Degree-1 terminal / degree-2 / "
        "degree-3 interior / anomalous atoms: "
''',
        "coordination print label",
    )

    text = replace_once(
        text,
        '''        f"{len(edge_indices)}/"
        f"{len(interior_indices)}/"
        f"{len(anomalous_indices)}"
''',
        '''        f"{len(edge_indices)}/"
        f"{len(degree2_indices)}/"
        f"{len(interior_indices)}/"
        f"{len(anomalous_indices)}"
''',
        "coordination print values",
    )

    text = replace_once(
        text,
        '''        "Lower rim atoms / B / N / radius / "
        "axial std / alternation: "
''',
        '''        "Lower rim atoms / B / N / radius / "
        "axial std / dominant / purity: "
''',
        "lower print label",
    )

    text = replace_once(
        text,
        '''        f"{lower_metrics['radius_mean_nm']:.6f}/"
        f"{lower_metrics['axial_standard_deviation_nm']:.6f}/"
        f"{lower_metrics['element_alternation_fraction']:.6f}"
''',
        '''        f"{lower_metrics['radius_mean_nm']:.6f}/"
        f"{lower_metrics['axial_standard_deviation_nm']:.6f}/"
        f"{lower_metrics['dominant_element']}/"
        f"{lower_metrics['element_purity_fraction']:.6f}"
''',
        "lower print values",
    )

    text = replace_once(
        text,
        '''        "Upper rim atoms / B / N / radius / "
        "axial std / alternation: "
''',
        '''        "Upper rim atoms / B / N / radius / "
        "axial std / dominant / purity: "
''',
        "upper print label",
    )

    text = replace_once(
        text,
        '''        f"{upper_metrics['radius_mean_nm']:.6f}/"
        f"{upper_metrics['axial_standard_deviation_nm']:.6f}/"
        f"{upper_metrics['element_alternation_fraction']:.6f}"
''',
        '''        f"{upper_metrics['radius_mean_nm']:.6f}/"
        f"{upper_metrics['axial_standard_deviation_nm']:.6f}/"
        f"{upper_metrics['dominant_element']}/"
        f"{upper_metrics['element_purity_fraction']:.6f}"
''',
        "upper print values",
    )

    stale_tokens = (
        "EXPECTED_EDGE_B_PER_END",
        "EXPECTED_EDGE_N_PER_END",
        "MIN_EDGE_ALTERNATION_FRACTION",
        "geometry_has_120_degree2_edge_atoms",
        "geometry_has_1560_degree3_interior_atoms",
        "lower_end_has_30B_and_30N",
        "upper_end_has_30B_and_30N",
    )

    remaining = [
        token
        for token in stale_tokens
        if token in text
    ]

    if remaining:
        raise RuntimeError(
            "Stale pre-classification assumptions remain: "
            + " | ".join(
                remaining
            )
        )

    AUDIT.write_text(
        text,
        encoding="utf-8",
    )

    print(
        "Day024 parent-rim auditor repaired for "
        "degree-1 polar termini."
    )

    print(
        "Validated terminal model: "
        "60 degree-1 sites and 1620 degree-3 sites."
    )

    print(
        "Lower end: 30 B / 0 N."
    )

    print(
        "Upper end: 0 B / 30 N."
    )

    print(
        "Explicit and geometric connectivity prerequisite: PASS"
    )

    print(
        "No geometry, minimization, MD, or QM was executed."
    )


if __name__ == "__main__":
    main()
