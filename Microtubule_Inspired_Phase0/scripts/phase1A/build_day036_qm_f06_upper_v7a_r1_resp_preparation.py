#!/usr/bin/env python3
"""
Build the Day036 RESP-preparation gate for QM_F06 UPPER V7-A R1.

This stage:
- validates the adopted geometry and RESP-design authorization;
- preserves atom indexing and provenance;
- classifies real atoms and artificial QM caps;
- creates conservative equivalence candidates without enforcing them;
- defines the artificial-cap transfer policy;
- prepares a candidate ORCA ESP/CHELPG single-point input;
- writes a reproducible execution manifest.

This stage DOES NOT execute ORCA, ESP fitting, RESP, charge adoption,
force-field adoption, or molecular dynamics.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
import csv
import hashlib
import json


ROOT = Path(__file__).resolve().parents[2]

ADOPTION_DIR = (
    ROOT
    / "runs/phase1A/"
    / "day035_qm_f06_upper_v7a_r1_coordinate_adoption"
)

ADOPTION_REPORT = (
    ADOPTION_DIR
    / "QM_F06_UPPER_V7A_COORDINATE_ADOPTION.json"
)

RESP_DESIGN_REPORT = (
    ROOT
    / "runs/phase1A/"
    / "day035_qm_f06_upper_v7a_r1_resp_input_design/"
    / "QM_F06_UPPER_V7A_R1_RESP_INPUT_DESIGN.json"
)

SOURCE_DIAGNOSTIC_INPUT = (
    ROOT
    / "runs/phase1A/"
    / "day028_qm_f06_lower_boundary_v2b_electronic_diagnostic/"
    / "QM_F06_LOWER_BOUNDARY_V2B_ELECTRONIC_DIAGNOSTIC.inp"
)

OUTPUT_DIR = (
    ROOT
    / "runs/phase1A/"
    / "day036_qm_f06_upper_v7a_r1_resp_preparation"
)

README = OUTPUT_DIR / "README.md"

OUTPUT_PREPARATION = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_RESP_PREPARATION.json"
)

OUTPUT_ATOM_CLASSES = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_RESP_ATOM_CLASSES.csv"
)

OUTPUT_EQUIVALENCE = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_RESP_EQUIVALENCE_GROUPS.csv"
)

OUTPUT_CAP_POLICY = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_RESP_CAP_POLICY.csv"
)

OUTPUT_MANIFEST = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_RESP_EXECUTION_MANIFEST.json"
)

OUTPUT_PROTOCOL = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_RESP_PROTOCOL.md"
)

OUTPUT_ORCA_INPUT = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_ESP_ORCA.inp"
)

EXPECTED_ADOPTION_DECISION = (
    "QM_F06_UPPER_V7A_"
    "FINAL_COORDINATES_ADOPTED_"
    "RESP_INPUT_DESIGN_AUTHORIZED"
)

EXPECTED_RESP_DESIGN_DECISION = (
    "QM_F06_UPPER_V7A_R1_"
    "RESP_INPUT_DESIGN_PASS_"
    "RESP_INPUT_PREPARATION_AUTHORIZED"
)

EXPECTED_ATOM_COUNT = 52
EXPECTED_COMPOSITION = {
    "B": 17,
    "N": 14,
    "H": 21,
}
EXPECTED_CAP_COUNT = 15
EXPECTED_REAL_ATOM_COUNT = 37

NET_CHARGE = 0
MULTIPLICITY = 1

ORCA_METHOD_LINE = (
    "! PBE0 D4 def2-TZVP def2/J RIJCOSX "
    "TightSCF DefGrid3 CHELPG"
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def read_json(path: Path) -> dict:
    if not path.is_file():
        raise RuntimeError(
            f"Required JSON is missing: {path}"
        )

    return json.loads(
        path.read_text(encoding="utf-8")
    )


def read_csv(path: Path) -> list[dict]:
    if not path.is_file():
        raise RuntimeError(
            f"Required CSV is missing: {path}"
        )

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        return list(csv.DictReader(handle))


def write_csv(
    path: Path,
    fieldnames: list[str],
    rows: list[dict],
) -> None:
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


def read_xyz(path: Path) -> list[dict]:
    if not path.is_file():
        raise RuntimeError(
            f"Required XYZ is missing: {path}"
        )

    lines = path.read_text(
        encoding="utf-8"
    ).splitlines()

    if len(lines) < 2:
        raise RuntimeError(
            f"Incomplete XYZ file: {path}"
        )

    atom_count = int(lines[0].strip())
    atoms = []

    for index, line in enumerate(
        lines[2:2 + atom_count]
    ):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Incomplete XYZ row {index}: {path}"
            )

        atoms.append({
            "index_0based": index,
            "element": fields[0],
            "x_A": float(fields[1]),
            "y_A": float(fields[2]),
            "z_A": float(fields[3]),
        })

    if len(atoms) != atom_count:
        raise RuntimeError(
            f"XYZ atom-count mismatch: {path}"
        )

    return atoms


def parse_bool(value: str) -> bool:
    normalized = str(value).strip().lower()

    if normalized == "true":
        return True

    if normalized == "false":
        return False

    raise RuntimeError(
        f"Invalid Boolean value: {value!r}"
    )


def extract_chelpg_block(path: Path) -> list[str]:
    """
    Reuse the exact CHELPG block already exercised in the project.
    """
    if not path.is_file():
        raise RuntimeError(
            "Reference CHELPG input is missing: "
            f"{path}"
        )

    lines = path.read_text(
        encoding="utf-8"
    ).splitlines()

    start = None

    for index, line in enumerate(lines):
        if line.strip().lower() == "%chelpg":
            start = index
            break

    if start is None:
        raise RuntimeError(
            "No %chelpg block found in reference input."
        )

    block = []

    for line in lines[start:]:
        block.append(line)

        if (
            len(block) > 1
            and line.strip().lower() == "end"
        ):
            return block

    raise RuntimeError(
        "Reference %chelpg block has no terminating end."
    )


def build_atom_class(
    row: dict,
) -> tuple[str, str, str]:
    artificial_cap = parse_bool(
        row["artificial_cap"]
    )

    element = row["element"]
    atom_role = row["atom_role"]
    node_type = row["node_type"]

    if artificial_cap:
        transfer_status = (
            "QM_ONLY_EXCLUDED_FROM_DIRECT_TRANSFER"
        )
        class_id = (
            f"CAP__{element}__{node_type}__{atom_role}"
        )
        class_basis = (
            "ARTIFICIAL_CAP_ELEMENT_NODE_TYPE_ROLE"
        )
    else:
        transfer_status = "TRANSFERABLE_REAL_ATOM"
        class_id = (
            f"REAL__{element}__{node_type}__{atom_role}"
        )
        class_basis = (
            "REAL_ATOM_ELEMENT_NODE_TYPE_ROLE"
        )

    return (
        class_id,
        class_basis,
        transfer_status,
    )


def main() -> None:
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    adoption = read_json(
        ADOPTION_REPORT
    )

    resp_design = read_json(
        RESP_DESIGN_REPORT
    )

    gates = {}

    gates["adoption_decision_matches"] = (
        adoption.get("decision")
        == EXPECTED_ADOPTION_DECISION
    )

    gates["RESP_input_design_authorized"] = (
        adoption.get(
            "authorizations",
            {},
        ).get(
            "RESP_input_design_authorized"
        )
        is True
    )

    gates["RESP_design_decision_matches"] = (
        resp_design.get("decision")
        == EXPECTED_RESP_DESIGN_DECISION
    )

    gates[
        "RESP_input_preparation_authorized"
    ] = (
        resp_design.get(
            "authorizations",
            {},
        ).get(
            "RESP_input_preparation_authorized"
        )
        is True
    )

    adopted_xyz = Path(
        adoption["outputs"][
            "adopted_final_XYZ"
        ]
    )

    adopted_map = Path(
        adoption["outputs"][
            "adopted_atom_role_provenance_map"
        ]
    )

    adopted_edges = Path(
        adoption["outputs"][
            "adopted_nominal_edges"
        ]
    )

    atoms = read_xyz(adopted_xyz)
    map_rows = read_csv(adopted_map)
    edge_rows = read_csv(adopted_edges)

    map_rows.sort(
        key=lambda row: int(
            row["v7a_index_0based"]
        )
    )

    gates["atom_count_52"] = (
        len(atoms)
        == len(map_rows)
        == EXPECTED_ATOM_COUNT
    )

    composition = Counter(
        atom["element"]
        for atom in atoms
    )

    gates["composition_B17_N14_H21"] = (
        dict(composition)
        == EXPECTED_COMPOSITION
    )

    gates["map_element_order_matches_XYZ"] = (
        len(atoms) == len(map_rows)
        and all(
            atom["element"]
            == row["element"]
            for atom, row in zip(
                atoms,
                map_rows,
            )
        )
    )

    gates["nominal_edge_count_57"] = (
        len(edge_rows) == 57
    )

    gates["geometry_hash_matches_design"] = (
        sha256(adopted_xyz)
        == resp_design.get(
            "geometry_sha256"
        )
    )

    gates["geometry_hash_matches_adoption"] = (
        sha256(adopted_xyz)
        == adoption.get(
            "sha256",
            {},
        ).get(
            "adopted_final_XYZ"
        )
    )

    artificial_caps = [
        row
        for row in map_rows
        if parse_bool(
            row["artificial_cap"]
        )
    ]

    real_atoms = [
        row
        for row in map_rows
        if not parse_bool(
            row["artificial_cap"]
        )
    ]

    gates["artificial_cap_count_15"] = (
        len(artificial_caps)
        == EXPECTED_CAP_COUNT
    )

    gates["real_atom_count_37"] = (
        len(real_atoms)
        == EXPECTED_REAL_ATOM_COUNT
    )

    gates["charge_and_multiplicity_defined"] = (
        NET_CHARGE == 0
        and MULTIPLICITY == 1
    )

    chelpg_block = extract_chelpg_block(
        SOURCE_DIAGNOSTIC_INPUT
    )

    gates["reference_CHELPG_block_available"] = (
        len(chelpg_block) >= 2
    )

    if not all(gates.values()):
        failed = [
            name
            for name, passed in gates.items()
            if not passed
        ]

        raise RuntimeError(
            "RESP preparation upstream gate failed: "
            + " | ".join(failed)
        )

    atom_class_rows = []
    equivalence_rows = []
    cap_policy_rows = []

    candidate_groups: dict[
        str,
        list[int],
    ] = defaultdict(list)

    class_counts: Counter = Counter()

    for atom, row in zip(
        atoms,
        map_rows,
    ):
        index = int(
            row["v7a_index_0based"]
        )

        class_id, class_basis, transfer_status = (
            build_atom_class(row)
        )

        artificial_cap = parse_bool(
            row["artificial_cap"]
        )

        candidate_group_key = (
            f"{transfer_status}"
            f"|{row['element']}"
            f"|{row['node_type']}"
            f"|{row['atom_role']}"
        )

        candidate_groups[
            candidate_group_key
        ].append(index)

        class_counts[class_id] += 1

        atom_class_rows.append({
            "atom_index_0based": index,
            "atom_index_1based": index + 1,
            "atom_id": row["atom_id"],
            "element": row["element"],
            "atom_role": row["atom_role"],
            "node_type": row["node_type"],
            "artificial_cap": artificial_cap,
            "transfer_status": transfer_status,
            "RESP_atom_class": class_id,
            "classification_basis": class_basis,
            "candidate_equivalence_key": (
                candidate_group_key
            ),
            "equivalence_enforced": False,
            "x_A": f"{atom['x_A']:.12f}",
            "y_A": f"{atom['y_A']:.12f}",
            "z_A": f"{atom['z_A']:.12f}",
        })

        equivalence_rows.append({
            "atom_index_0based": index,
            "atom_index_1based": index + 1,
            "atom_id": row["atom_id"],
            "element": row["element"],
            "candidate_group_key": (
                candidate_group_key
            ),
            "candidate_group_size": 0,
            "enforced_group_id": (
                f"SINGLETON_{index:03d}"
            ),
            "equivalence_enforced": False,
            "review_status": (
                "TOPOLOGY_AND_LOCAL_GEOMETRY_"
                "REVIEW_REQUIRED"
            ),
            "scientific_basis": (
                "NO_EQUIVALENCE_ASSUMED_FROM_"
                "ELEMENT_OR_LABEL_ALONE"
            ),
        })

        if artificial_cap:
            cap_policy_rows.append({
                "atom_index_0based": index,
                "atom_index_1based": index + 1,
                "atom_id": row["atom_id"],
                "element": row["element"],
                "atom_role": row["atom_role"],
                "node_type": row["node_type"],
                "included_in_QM_ESP_fit": True,
                "directly_transferable_to_full_scaffold": (
                    False
                ),
                "transfer_policy": (
                    "EXCLUDE_FROM_TRANSFERABLE_"
                    "ATOM_TYPE_AVERAGES"
                ),
                "charge_redistribution_policy": (
                    "PENDING_JOINT_LOWER_UPPER_"
                    "PROTOCOL_VALIDATION"
                ),
                "charge_adoption_authorized": False,
            })

    group_sizes = {
        key: len(indices)
        for key, indices
        in candidate_groups.items()
    }

    for row in equivalence_rows:
        row["candidate_group_size"] = (
            group_sizes[
                row["candidate_group_key"]
            ]
        )

    write_csv(
        OUTPUT_ATOM_CLASSES,
        list(atom_class_rows[0].keys()),
        atom_class_rows,
    )

    write_csv(
        OUTPUT_EQUIVALENCE,
        list(equivalence_rows[0].keys()),
        equivalence_rows,
    )

    write_csv(
        OUTPUT_CAP_POLICY,
        list(cap_policy_rows[0].keys()),
        cap_policy_rows,
    )

    orca_input_lines = [
        ORCA_METHOD_LINE,
        "",
        "%pal",
        "  nprocs 4",
        "end",
        "",
        "%maxcore 2500",
        "",
        "%scf",
        "  MaxIter 500",
        "  AutoTRAH false",
        "end",
        "",
        "# Day036 candidate ESP/CHELPG single-point input",
        "# Uses the formally adopted V7-A R1 geometry",
        "# Preparation only; execution is not authorized",
        "",
        *chelpg_block,
        "",
        f"* xyz {NET_CHARGE} {MULTIPLICITY}",
    ]

    for atom in atoms:
        orca_input_lines.append(
            f"{atom['element']:<2s} "
            f"{atom['x_A']: .12f} "
            f"{atom['y_A']: .12f} "
            f"{atom['z_A']: .12f}"
        )

    orca_input_lines.append("*")

    OUTPUT_ORCA_INPUT.write_text(
        "\n".join(orca_input_lines) + "\n",
        encoding="utf-8",
    )

    preparation_decision = (
        "QM_F06_UPPER_V7A_R1_"
        "RESP_PREPARATION_PASS_"
        "ESP_INPUT_PREFLIGHT_AUTHORIZED"
    )

    preparation_report = {
        "schema_version": 1,
        "generated_at_utc": utc_now(),
        "model": "QM_F06_UPPER_V7A_R1",
        "stage": "DAY036_RESP_PREPARATION",
        "decision": preparation_decision,
        "upstream": {
            "coordinate_adoption_report": str(
                ADOPTION_REPORT
            ),
            "RESP_input_design_report": str(
                RESP_DESIGN_REPORT
            ),
            "execution_directory": adoption[
                "source_execution_directory"
            ],
        },
        "electronic_structure": {
            "method": "PBE0-D4",
            "basis": "def2-TZVP",
            "Coulomb_auxiliary_basis": "def2/J",
            "exchange_approximation": "RIJCOSX",
            "SCF": "TightSCF",
            "grid": "DefGrid3",
            "net_charge": NET_CHARGE,
            "multiplicity": MULTIPLICITY,
            "ESP_source": "ORCA_CHELPG",
        },
        "geometry": {
            "path": str(adopted_xyz),
            "sha256": sha256(adopted_xyz),
            "atom_count": len(atoms),
            "composition": dict(
                sorted(composition.items())
            ),
        },
        "classification": {
            "real_atom_count": len(real_atoms),
            "artificial_cap_count": len(
                artificial_caps
            ),
            "RESP_atom_class_count": len(
                class_counts
            ),
            "candidate_equivalence_group_count": (
                len(candidate_groups)
            ),
            "enforced_nonsingleton_equivalence_groups": 0,
        },
        "cap_policy": {
            "caps_included_in_QM_ESP_fit": True,
            "caps_directly_transferable": False,
            "policy": (
                "EXCLUDE_ARTIFICIAL_CAPS_FROM_"
                "TRANSFERABLE_ATOM_TYPE_AVERAGES"
            ),
            "redistribution_policy": (
                "PENDING_JOINT_LOWER_UPPER_"
                "PROTOCOL_VALIDATION"
            ),
        },
        "gates": gates,
        "outputs": {
            "atom_classes": str(
                OUTPUT_ATOM_CLASSES
            ),
            "equivalence_groups": str(
                OUTPUT_EQUIVALENCE
            ),
            "cap_policy": str(
                OUTPUT_CAP_POLICY
            ),
            "candidate_ORCA_ESP_input": str(
                OUTPUT_ORCA_INPUT
            ),
            "RESP_protocol": str(
                OUTPUT_PROTOCOL
            ),
            "execution_manifest": str(
                OUTPUT_MANIFEST
            ),
        },
        "authorizations": {
            "ESP_input_preflight_authorized": True,
            "ESP_execution_authorized": False,
            "RESP_input_generation_authorized": False,
            "RESP_execution_authorized": False,
            "RESP_validation_authorized": False,
            "charge_adoption_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    manifest = {
        "schema_version": 1,
        "generated_at_utc": utc_now(),
        "model": "QM_F06_UPPER_V7A_R1",
        "candidate_ORCA_executable": (
            "/Users/alejandro/projects/"
            "orca_6_1_1_macosx_intel_openmpi411/"
            "orca"
        ),
        "candidate_input": str(
            OUTPUT_ORCA_INPUT
        ),
        "working_directory": str(
            OUTPUT_DIR
        ),
        "net_charge": NET_CHARGE,
        "multiplicity": MULTIPLICITY,
        "nprocs": 4,
        "maxcore_MB_per_process": 2500,
        "input_sha256": sha256(
            OUTPUT_ORCA_INPUT
        ),
        "geometry_sha256": sha256(
            adopted_xyz
        ),
        "execution_authorized": False,
        "required_next_gate": (
            "INDEPENDENT_ESP_INPUT_PREFLIGHT"
        ),
        "prohibited_actions": [
            "DO_NOT_EXECUTE_ORCA_ESP",
            "DO_NOT_RUN_RESP",
            "DO_NOT_ADOPT_CHARGES",
            "DO_NOT_UPDATE_FORCE_FIELD",
            "DO_NOT_RUN_MD",
        ],
    }

    protocol_text = f"""# QM_F06 UPPER V7-A R1 RESP Protocol — Day036

## Current formal state

Decision:

`{preparation_decision}`

The UPPER V7-A R1 geometry has been optimized, independently audited,
coordinate-validated, and formally adopted.

This stage prepares, but does not execute, the ESP/RESP workflow.

## Adopted electronic-structure specification

- Method: PBE0-D4
- Basis: def2-TZVP
- Coulomb auxiliary basis: def2/J
- Exchange approximation: RIJCOSX
- SCF convergence: TightSCF
- Integration grid: DefGrid3
- Net charge: 0
- Multiplicity: 1
- Candidate ESP source: ORCA CHELPG

## Geometry and atom inventory

- Total atoms: {len(atoms)}
- Real transferable atoms: {len(real_atoms)}
- Artificial QM caps: {len(artificial_caps)}
- Composition: B17N14H21
- Adopted geometry SHA256: `{sha256(adopted_xyz)}`

## Artificial-cap policy

Artificial caps remain included in the QM electrostatic calculation
because they are part of the finite electronic-structure model.

Their fitted charges must not be transferred directly into the full
R2 scaffold or included in transferable atom-type averages.

Final cap-charge removal or redistribution remains pending a joint
LOWER/UPPER protocol validation.

## Equivalence policy

No non-singleton RESP equivalence restraint is enforced at this stage.

Candidate groups are reported from element, node type, atom role, and
transfer status. They require independent topology and local-geometry
validation before any equality constraint may be applied.

Element identity alone is not accepted as evidence of charge
equivalence.

## Conformational scope

The current UPPER V7-A R1 geometry is an accepted reference geometry.

A transferable production charge model must later evaluate consistency
against the accepted LOWER reference and determine whether additional
conformers are required.

## Authorization state

- RESP preparation gate: PASS
- ESP input preflight: AUTHORIZED
- ESP execution: NOT AUTHORIZED
- RESP input generation: NOT AUTHORIZED
- RESP execution: NOT AUTHORIZED
- Charge adoption: NOT AUTHORIZED
- Force-field adoption: NOT AUTHORIZED
- Molecular dynamics: NOT AUTHORIZED

## Required next step

Run an independent preflight of the candidate ORCA ESP input, including:

1. executable-path validation;
2. geometry-hash validation;
3. charge and multiplicity validation;
4. CHELPG-grid validation;
5. method and basis validation;
6. output-file and working-directory isolation;
7. explicit execution authorization.
"""

    OUTPUT_PREPARATION.write_text(
        json.dumps(
            preparation_report,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    OUTPUT_MANIFEST.write_text(
        json.dumps(
            manifest,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    OUTPUT_PROTOCOL.write_text(
        protocol_text,
        encoding="utf-8",
    )

    README.write_text(
        """# Day036 – QM_F06 UPPER V7-A R1 RESP Preparation

This directory contains the reproducible preparation gate for the
UPPER V7-A R1 electrostatic-parameterization workflow.

The adopted geometry, atom classification, candidate equivalence
groups, artificial-cap policy, ESP protocol and execution manifest
are generated by:

`scripts/phase1A/build_day036_qm_f06_upper_v7a_r1_resp_preparation.py`

This stage does not execute ORCA ESP, RESP, charge adoption,
force-field modification or molecular dynamics.
""",
        encoding="utf-8",
    )

    print("=" * 104)
    print("QM_F06 UPPER V7-A R1 RESP PREPARATION — DAY036")
    print("=" * 104)

    for name, value in gates.items():
        print(
            f"{name:58s}: "
            f"{'PASS' if value else 'FAIL'}"
        )

    print()
    print(f"Atoms: {len(atoms)}")
    print(f"Real atoms: {len(real_atoms)}")
    print(
        f"Artificial caps: "
        f"{len(artificial_caps)}"
    )
    print(
        f"RESP atom classes: "
        f"{len(class_counts)}"
    )
    print(
        "Enforced nonsingleton equivalence groups: 0"
    )

    print()
    print("Decision:", preparation_decision)
    print("Preparation report:", OUTPUT_PREPARATION)
    print("Atom classes:", OUTPUT_ATOM_CLASSES)
    print("Equivalence groups:", OUTPUT_EQUIVALENCE)
    print("Cap policy:", OUTPUT_CAP_POLICY)
    print("Protocol:", OUTPUT_PROTOCOL)
    print("Candidate ORCA ESP input:", OUTPUT_ORCA_INPUT)
    print("Execution manifest:", OUTPUT_MANIFEST)

    print()
    print("ESP input preflight authorized: True")
    print("ESP execution authorized: False")
    print("RESP execution authorized: False")
    print("Charge adoption authorized: False")
    print("Force-field adoption authorized: False")
    print("MD authorized: False")


if __name__ == "__main__":
    main()
