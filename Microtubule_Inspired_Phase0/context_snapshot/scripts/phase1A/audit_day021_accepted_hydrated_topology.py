#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import re
import shutil
import subprocess
from collections import Counter
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

ACCEPTED_RUN_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/accepted/"
    "hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute"
)

CANDIDATE_TOPOLOGY_ROOT = (
    PROJECT_ROOT
    / "parameters/phase1A/accepted/"
    "hybrid_hbnBonded_kang2000_improperGeo100_validated"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day021_accepted_hydrated_topology_audit"
)

TPR_PATH = (
    ACCEPTED_RUN_ROOT
    / "nvt_100ps_frozenSolute.tpr"
)

GRO_PATH = (
    ACCEPTED_RUN_ROOT
    / "nvt_100ps_frozenSolute.gro"
)

MDP_PATH = (
    ACCEPTED_RUN_ROOT
    / "nvt_100ps_frozenSolute.mdp"
)

TOP_PATH = (
    CANDIDATE_TOPOLOGY_ROOT
    / (
        "hbn_pyrene_4_hydratable_gap45_"
        "pyr5shift_clean032_hbnBonded_"
        "kang2000_improperGeo100.top"
    )
)

REBUILT_TPR = (
    OUTPUT_ROOT
    / "candidate_rebuilt_from_accepted_gro.tpr"
)

PROCESSED_TOP = (
    OUTPUT_ROOT
    / "candidate_processed.top"
)

MDOUT_MDP = (
    OUTPUT_ROOT
    / "candidate_mdout.mdp"
)

GROMPP_LOG = (
    OUTPUT_ROOT
    / "grompp_candidate_rebuild.log"
)

ACCEPTED_TPR_DUMP = (
    OUTPUT_ROOT
    / "accepted_tpr_dump.txt"
)

REBUILT_TPR_DUMP = (
    OUTPUT_ROOT
    / "rebuilt_tpr_dump.txt"
)

TPR_COMPARISON_LOG = (
    OUTPUT_ROOT
    / "accepted_vs_rebuilt_tpr_check.log"
)

FILE_INVENTORY_CSV = (
    OUTPUT_ROOT
    / "topology_file_inventory.csv"
)

INCLUDES_CSV = (
    OUTPUT_ROOT
    / "topology_include_hierarchy.csv"
)

MOLECULE_TYPES_CSV = (
    OUTPUT_ROOT
    / "topology_molecule_types.csv"
)

MOLECULES_CSV = (
    OUTPUT_ROOT
    / "topology_system_molecules.csv"
)

GRO_COMPOSITION_CSV = (
    OUTPUT_ROOT
    / "accepted_gro_composition.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "accepted_hydrated_topology_audit_summary.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "ACCEPTED_HYDRATED_TOPOLOGY_AUDIT_DAY021.md"
)

EXPECTED_RESIDUE_COUNTS = {
    "HBN": 1,
    "PYR": 4,
    "SOL": 16634,
}

EXPECTED_ATOM_COUNTS = {
    "HBN": 1680,
    "PYR": 104,
    "SOL": 66536,
}

EXPECTED_TOTAL_ATOMS = 68320

INCLUDE_PATTERN = re.compile(
    r'^\s*#include\s+[<"]([^>"]+)[>"]'
)

SECTION_PATTERN = re.compile(
    r"^\s*\[\s*([^\]]+)\s*\]\s*$"
)

TPR_NATOMS_PATTERN = re.compile(
    r"\bnatoms\s*=\s*(\d+)"
)


def log(message: str = "") -> None:
    print(message, flush=True)


def relative(path: Path) -> str:
    try:
        return str(
            path.relative_to(PROJECT_ROOT)
        )
    except ValueError:
        return str(path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            block = handle.read(
                1024 * 1024
            )

            if not block:
                break

            digest.update(block)

    return digest.hexdigest()


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

    normalized = [
        {
            field: row.get(field, "")
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
        writer.writerows(normalized)


def require_inputs() -> None:
    required = (
        TPR_PATH,
        GRO_PATH,
        MDP_PATH,
        TOP_PATH,
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
            "Missing or empty required inputs:\n"
            + "\n".join(
                str(path)
                for path in missing
            )
        )


def find_gmx() -> Path:
    preferred = Path(
        "/usr/local/gromacs/bin/gmx"
    )

    if preferred.exists():
        return preferred

    discovered = shutil.which("gmx")

    if discovered is None:
        raise RuntimeError(
            "Could not locate the GROMACS "
            "gmx executable"
        )

    return Path(discovered)


def run_command(
    command: list[str],
    cwd: Path,
    log_path: Path,
) -> subprocess.CompletedProcess[str]:
    environment = {
        **dict(__import__("os").environ),
        "GMX_MAXBACKUP": "-1",
    }

    completed = subprocess.run(
        command,
        cwd=cwd,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )

    log_path.write_text(
        completed.stdout,
        encoding="utf-8",
    )

    return completed


def strip_semicolon_comment(
    line: str,
) -> str:
    return line.split(
        ";",
        1,
    )[0].strip()


def resolve_include(
    token: str,
    parent_file: Path,
) -> Path | None:
    candidates = (
        parent_file.parent / token,
        CANDIDATE_TOPOLOGY_ROOT / token,
        PROJECT_ROOT / token,
    )

    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()

    return None


def collect_include_hierarchy(
    root_topology: Path,
) -> tuple[
    list[Path],
    list[dict[str, object]],
]:
    visited: set[Path] = set()
    ordered_files: list[Path] = []
    rows: list[dict[str, object]] = []

    def recurse(
        current: Path,
        depth: int,
    ) -> None:
        resolved_current = (
            current.resolve()
        )

        if resolved_current in visited:
            return

        visited.add(resolved_current)
        ordered_files.append(
            resolved_current
        )

        lines = resolved_current.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines()

        for line_number, line in enumerate(
            lines,
            start=1,
        ):
            match = INCLUDE_PATTERN.match(
                line
            )

            if not match:
                continue

            include_token = match.group(1)

            resolved_include = (
                resolve_include(
                    include_token,
                    resolved_current,
                )
            )

            rows.append(
                {
                    "parent_file": relative(
                        resolved_current
                    ),
                    "line_number": line_number,
                    "depth": depth,
                    "include_token": (
                        include_token
                    ),
                    "resolved_path": (
                        relative(
                            resolved_include
                        )
                        if resolved_include
                        else ""
                    ),
                    "resolved_locally": (
                        resolved_include
                        is not None
                    ),
                }
            )

            if resolved_include is not None:
                recurse(
                    resolved_include,
                    depth + 1,
                )

    recurse(
        root_topology,
        0,
    )

    return ordered_files, rows


def parse_molecule_types(
    topology_files: list[Path],
) -> tuple[
    dict[str, dict[str, object]],
    list[dict[str, object]],
]:
    definitions: dict[
        str,
        dict[str, object]
    ] = {}

    duplicate_rows: list[
        dict[str, object]
    ] = []

    for path in topology_files:
        current_section = ""
        current_molecule: str | None = None
        awaiting_molecule_name = False

        lines = path.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines()

        for line_number, raw_line in enumerate(
            lines,
            start=1,
        ):
            line = strip_semicolon_comment(
                raw_line
            )

            if not line:
                continue

            if line.startswith("#"):
                continue

            section_match = (
                SECTION_PATTERN.match(
                    line
                )
            )

            if section_match:
                current_section = (
                    section_match.group(1)
                    .strip()
                    .lower()
                )

                if (
                    current_section
                    == "moleculetype"
                ):
                    awaiting_molecule_name = (
                        True
                    )
                    current_molecule = None

                continue

            tokens = line.split()

            if (
                current_section
                == "moleculetype"
                and awaiting_molecule_name
            ):
                molecule_name = tokens[0]

                if molecule_name in definitions:
                    duplicate_rows.append(
                        {
                            "molecule_type": (
                                molecule_name
                            ),
                            "first_file": (
                                definitions[
                                    molecule_name
                                ]["source_file"]
                            ),
                            "duplicate_file": (
                                relative(path)
                            ),
                            "duplicate_line": (
                                line_number
                            ),
                        }
                    )
                else:
                    definitions[
                        molecule_name
                    ] = {
                        "molecule_type": (
                            molecule_name
                        ),
                        "source_file": (
                            relative(path)
                        ),
                        "definition_line": (
                            line_number
                        ),
                        "nrexcl": (
                            tokens[1]
                            if len(tokens) > 1
                            else ""
                        ),
                        "atom_count": 0,
                    }

                current_molecule = (
                    molecule_name
                )

                awaiting_molecule_name = (
                    False
                )

                continue

            if (
                current_section == "atoms"
                and current_molecule
                is not None
            ):
                try:
                    int(tokens[0])
                except (
                    ValueError,
                    IndexError,
                ):
                    continue

                definitions[
                    current_molecule
                ]["atom_count"] = (
                    int(
                        definitions[
                            current_molecule
                        ]["atom_count"]
                    )
                    + 1
                )

    return definitions, duplicate_rows


def parse_system_molecules(
    topology_path: Path,
) -> list[dict[str, object]]:
    current_section = ""
    rows: list[dict[str, object]] = []

    lines = topology_path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    for line_number, raw_line in enumerate(
        lines,
        start=1,
    ):
        line = strip_semicolon_comment(
            raw_line
        )

        if not line:
            continue

        if line.startswith("#"):
            continue

        section_match = (
            SECTION_PATTERN.match(line)
        )

        if section_match:
            current_section = (
                section_match.group(1)
                .strip()
                .lower()
            )
            continue

        if current_section != "molecules":
            continue

        tokens = line.split()

        if len(tokens) < 2:
            continue

        try:
            count = int(tokens[1])
        except ValueError:
            continue

        rows.append(
            {
                "molecule_type": tokens[0],
                "molecule_count": count,
                "source_file": (
                    relative(topology_path)
                ),
                "line_number": line_number,
            }
        )

    if not rows:
        raise RuntimeError(
            "No [ molecules ] entries were "
            "parsed from the candidate topology"
        )

    return rows


def parse_gro(
    path: Path,
) -> dict[str, object]:
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    if len(lines) < 3:
        raise RuntimeError(
            f"Invalid GRO file: {path}"
        )

    atom_count = int(
        lines[1].strip()
    )

    atom_lines = lines[
        2 : 2 + atom_count
    ]

    if len(atom_lines) != atom_count:
        raise RuntimeError(
            "GRO atom-count mismatch"
        )

    residue_counts: Counter[str] = Counter()
    atom_counts: Counter[str] = Counter()

    previous_residue_key: (
        tuple[int, str] | None
    ) = None

    velocity_line_count = 0

    for atom_line in atom_lines:
        residue_number = int(
            atom_line[0:5]
        )

        residue_name = (
            atom_line[5:10].strip()
        )

        residue_key = (
            residue_number,
            residue_name,
        )

        atom_counts[
            residue_name
        ] += 1

        if residue_key != (
            previous_residue_key
        ):
            residue_counts[
                residue_name
            ] += 1

            previous_residue_key = (
                residue_key
            )

        if len(atom_line) >= 68:
            velocity_line_count += 1

    return {
        "atom_count": atom_count,
        "residue_counts": residue_counts,
        "atom_counts": atom_counts,
        "velocity_line_count": (
            velocity_line_count
        ),
        "all_atoms_have_velocities": (
            velocity_line_count
            == atom_count
        ),
    }


def parse_tpr_natoms(
    dump_text: str,
) -> int:
    match = TPR_NATOMS_PATTERN.search(
        dump_text
    )

    if match is None:
        raise RuntimeError(
            "Could not parse natoms from "
            "gmx dump output"
        )

    return int(match.group(1))


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    require_inputs()

    gmx = find_gmx()

    topology_files, include_rows = (
        collect_include_hierarchy(
            TOP_PATH
        )
    )

    (
        molecule_types,
        duplicate_molecule_types,
    ) = parse_molecule_types(
        topology_files
    )

    system_molecules = (
        parse_system_molecules(
            TOP_PATH
        )
    )

    gro = parse_gro(
        GRO_PATH
    )

    inventory_paths = [
        TPR_PATH,
        GRO_PATH,
        MDP_PATH,
        TOP_PATH,
        *topology_files,
    ]

    unique_inventory_paths: list[Path] = []
    seen_inventory_paths: set[Path] = set()

    for path in inventory_paths:
        resolved = path.resolve()

        if resolved in seen_inventory_paths:
            continue

        seen_inventory_paths.add(
            resolved
        )

        unique_inventory_paths.append(
            resolved
        )

    inventory_rows = [
        {
            "path": relative(path),
            "suffix": path.suffix,
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in unique_inventory_paths
    ]

    write_csv(
        FILE_INVENTORY_CSV,
        inventory_rows,
    )

    write_csv(
        INCLUDES_CSV,
        (
            include_rows
            if include_rows
            else [
                {
                    "parent_file": (
                        relative(TOP_PATH)
                    ),
                    "line_number": "",
                    "depth": "",
                    "include_token": "",
                    "resolved_path": "",
                    "resolved_locally": "",
                }
            ]
        ),
    )

    molecule_type_rows = [
        molecule_types[name]
        for name in sorted(
            molecule_types
        )
    ]

    if duplicate_molecule_types:
        for row in duplicate_molecule_types:
            row[
                "molecule_type"
            ] = (
                "DUPLICATE:"
                + str(
                    row[
                        "molecule_type"
                    ]
                )
            )

        molecule_type_rows.extend(
            duplicate_molecule_types
        )

    write_csv(
        MOLECULE_TYPES_CSV,
        molecule_type_rows,
    )

    topology_total_atoms = 0
    molecule_rows: list[
        dict[str, object]
    ] = []

    missing_molecule_types: list[str] = []

    for molecule_entry in system_molecules:
        molecule_name = str(
            molecule_entry[
                "molecule_type"
            ]
        )

        molecule_count = int(
            molecule_entry[
                "molecule_count"
            ]
        )

        definition = molecule_types.get(
            molecule_name
        )

        if definition is None:
            atoms_per_molecule = None
            expected_atoms = None

            missing_molecule_types.append(
                molecule_name
            )
        else:
            atoms_per_molecule = int(
                definition[
                    "atom_count"
                ]
            )

            expected_atoms = (
                molecule_count
                * atoms_per_molecule
            )

            topology_total_atoms += (
                expected_atoms
            )

        gro_residue_count = int(
            gro[
                "residue_counts"
            ].get(
                molecule_name,
                0,
            )
        )

        gro_atom_count = int(
            gro[
                "atom_counts"
            ].get(
                molecule_name,
                0,
            )
        )

        molecule_rows.append(
            {
                **molecule_entry,
                "atoms_per_molecule": (
                    atoms_per_molecule
                ),
                "topology_atom_total": (
                    expected_atoms
                ),
                "GRO_residue_count": (
                    gro_residue_count
                ),
                "GRO_atom_count": (
                    gro_atom_count
                ),
                "molecule_count_matches_GRO": (
                    molecule_count
                    == gro_residue_count
                ),
                "atom_total_matches_GRO": (
                    expected_atoms
                    == gro_atom_count
                    if expected_atoms
                    is not None
                    else False
                ),
            }
        )

    write_csv(
        MOLECULES_CSV,
        molecule_rows,
    )

    gro_composition_rows = []

    all_residue_names = sorted(
        set(
            gro[
                "residue_counts"
            ].keys()
        )
        | set(
            gro[
                "atom_counts"
            ].keys()
        )
    )

    for residue_name in all_residue_names:
        gro_composition_rows.append(
            {
                "residue_name": residue_name,
                "residue_count": int(
                    gro[
                        "residue_counts"
                    ][residue_name]
                ),
                "atom_count": int(
                    gro[
                        "atom_counts"
                    ][residue_name]
                ),
                "expected_residue_count": (
                    EXPECTED_RESIDUE_COUNTS.get(
                        residue_name,
                        "",
                    )
                ),
                "expected_atom_count": (
                    EXPECTED_ATOM_COUNTS.get(
                        residue_name,
                        "",
                    )
                ),
            }
        )

    write_csv(
        GRO_COMPOSITION_CSV,
        gro_composition_rows,
    )

    accepted_dump_result = (
        run_command(
            [
                str(gmx),
                "dump",
                "-s",
                str(TPR_PATH),
            ],
            PROJECT_ROOT,
            ACCEPTED_TPR_DUMP,
        )
    )

    if accepted_dump_result.returncode != 0:
        raise RuntimeError(
            "gmx dump failed for the "
            "accepted TPR"
        )

    accepted_tpr_natoms = (
        parse_tpr_natoms(
            accepted_dump_result.stdout
        )
    )

    grompp_result = run_command(
        [
            str(gmx),
            "grompp",
            "-f",
            str(MDP_PATH),
            "-c",
            str(GRO_PATH),
            "-p",
            str(TOP_PATH),
            "-o",
            str(REBUILT_TPR),
            "-po",
            str(MDOUT_MDP),
            "-pp",
            str(PROCESSED_TOP),
            "-maxwarn",
            "0",
        ],
        CANDIDATE_TOPOLOGY_ROOT,
        GROMPP_LOG,
    )

    grompp_pass = (
        grompp_result.returncode == 0
        and REBUILT_TPR.exists()
        and REBUILT_TPR.stat().st_size > 0
    )

    rebuilt_tpr_natoms: int | None = None
    tpr_check_return_code: (
        int | None
    ) = None

    if grompp_pass:
        rebuilt_dump_result = (
            run_command(
                [
                    str(gmx),
                    "dump",
                    "-s",
                    str(REBUILT_TPR),
                ],
                PROJECT_ROOT,
                REBUILT_TPR_DUMP,
            )
        )

        if (
            rebuilt_dump_result.returncode
            == 0
        ):
            rebuilt_tpr_natoms = (
                parse_tpr_natoms(
                    rebuilt_dump_result.stdout
                )
            )

        tpr_check_result = run_command(
            [
                str(gmx),
                "check",
                "-s1",
                str(TPR_PATH),
                "-s2",
                str(REBUILT_TPR),
            ],
            PROJECT_ROOT,
            TPR_COMPARISON_LOG,
        )

        tpr_check_return_code = (
            tpr_check_result.returncode
        )
    else:
        REBUILT_TPR_DUMP.write_text(
            (
                "Not generated because "
                "grompp failed.\n"
            ),
            encoding="utf-8",
        )

        TPR_COMPARISON_LOG.write_text(
            (
                "Not generated because "
                "grompp failed.\n"
            ),
            encoding="utf-8",
        )

    gro_expected_composition_pass = all(
        int(
            gro[
                "residue_counts"
            ].get(
                residue_name,
                0,
            )
        )
        == expected_count
        for (
            residue_name,
            expected_count,
        )
        in EXPECTED_RESIDUE_COUNTS.items()
    ) and all(
        int(
            gro[
                "atom_counts"
            ].get(
                residue_name,
                0,
            )
        )
        == expected_count
        for (
            residue_name,
            expected_count,
        )
        in EXPECTED_ATOM_COUNTS.items()
    )

    topology_molecule_counts_pass = all(
        bool(
            row[
                "molecule_count_matches_GRO"
            ]
        )
        for row in molecule_rows
    )

    topology_atom_counts_pass = all(
        bool(
            row[
                "atom_total_matches_GRO"
            ]
        )
        for row in molecule_rows
    )

    unresolved_local_includes = [
        row
        for row in include_rows
        if not bool(
            row[
                "resolved_locally"
            ]
        )
    ]

    gro_tpr_atom_count_pass = (
        int(gro["atom_count"])
        == accepted_tpr_natoms
        == EXPECTED_TOTAL_ATOMS
    )

    topology_total_atom_count_pass = (
        topology_total_atoms
        == int(gro["atom_count"])
        == EXPECTED_TOTAL_ATOMS
    )

    rebuilt_atom_count_pass = (
        rebuilt_tpr_natoms
        == EXPECTED_TOTAL_ATOMS
        if rebuilt_tpr_natoms
        is not None
        else False
    )

    no_duplicate_molecule_types_pass = (
        len(
            duplicate_molecule_types
        )
        == 0
    )

    no_missing_molecule_types_pass = (
        len(
            missing_molecule_types
        )
        == 0
    )

    overall_pass = all(
        (
            gro_expected_composition_pass,
            gro_tpr_atom_count_pass,
            topology_total_atom_count_pass,
            topology_molecule_counts_pass,
            topology_atom_counts_pass,
            no_duplicate_molecule_types_pass,
            no_missing_molecule_types_pass,
            grompp_pass,
            rebuilt_atom_count_pass,
        )
    )

    summary_row = {
        "accepted_TPR": relative(
            TPR_PATH
        ),
        "accepted_GRO": relative(
            GRO_PATH
        ),
        "accepted_MDP": relative(
            MDP_PATH
        ),
        "candidate_TOP": relative(
            TOP_PATH
        ),
        "GROMACS_executable": str(gmx),
        "GRO_atom_count": int(
            gro["atom_count"]
        ),
        "accepted_TPR_atom_count": (
            accepted_tpr_natoms
        ),
        "candidate_topology_atom_count": (
            topology_total_atoms
        ),
        "rebuilt_TPR_atom_count": (
            rebuilt_tpr_natoms
        ),
        "GRO_all_atoms_have_velocities": (
            gro[
                "all_atoms_have_velocities"
            ]
        ),
        "local_topology_file_count": (
            len(topology_files)
        ),
        "include_directive_count": (
            len(include_rows)
        ),
        "unresolved_include_count": (
            len(
                unresolved_local_includes
            )
        ),
        "molecule_type_count": (
            len(molecule_types)
        ),
        "duplicate_molecule_type_count": (
            len(
                duplicate_molecule_types
            )
        ),
        "missing_molecule_types": (
            ",".join(
                missing_molecule_types
            )
        ),
        "grompp_return_code": (
            grompp_result.returncode
        ),
        "grompp_pass": grompp_pass,
        "gmx_check_return_code": (
            tpr_check_return_code
        ),
        "GRO_expected_composition_pass": (
            gro_expected_composition_pass
        ),
        "GRO_vs_TPR_atom_count_pass": (
            gro_tpr_atom_count_pass
        ),
        "topology_total_atom_count_pass": (
            topology_total_atom_count_pass
        ),
        "topology_molecule_counts_pass": (
            topology_molecule_counts_pass
        ),
        "topology_atom_counts_pass": (
            topology_atom_counts_pass
        ),
        "rebuilt_TPR_atom_count_pass": (
            rebuilt_atom_count_pass
        ),
        "no_duplicate_molecule_types_pass": (
            no_duplicate_molecule_types_pass
        ),
        "no_missing_molecule_types_pass": (
            no_missing_molecule_types_pass
        ),
        "structural_reconstruction_pass": (
            overall_pass
        ),
    }

    write_csv(
        SUMMARY_CSV,
        [summary_row],
    )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day021 Accepted Hydrated "
            "Topology Audit\n\n"
        )

        handle.write(
            "## Purpose\n\n"
        )

        handle.write(
            "Audit the strongest accepted "
            "hydrated-topology candidate against "
            "the accepted GRO and TPR before any "
            "mobile-solute simulation is prepared.\n\n"
        )

        handle.write(
            "## Accepted system\n\n"
        )

        handle.write(
            f"- GRO atoms: "
            f"{gro['atom_count']}.\n"
        )

        handle.write(
            f"- Accepted TPR atoms: "
            f"{accepted_tpr_natoms}.\n"
        )

        handle.write(
            f"- HBN residues/atoms: "
            f"{gro['residue_counts'].get('HBN', 0)}/"
            f"{gro['atom_counts'].get('HBN', 0)}.\n"
        )

        handle.write(
            f"- PYR residues/atoms: "
            f"{gro['residue_counts'].get('PYR', 0)}/"
            f"{gro['atom_counts'].get('PYR', 0)}.\n"
        )

        handle.write(
            f"- SOL residues/atoms: "
            f"{gro['residue_counts'].get('SOL', 0)}/"
            f"{gro['atom_counts'].get('SOL', 0)}.\n"
        )

        handle.write(
            f"- GRO velocity records: "
            f"{gro['velocity_line_count']}/"
            f"{gro['atom_count']}.\n\n"
        )

        handle.write(
            "## Candidate topology reconstruction\n\n"
        )

        handle.write(
            f"- Candidate topology: "
            f"`{relative(TOP_PATH)}`.\n"
        )

        handle.write(
            f"- Parsed topology atom total: "
            f"{topology_total_atoms}.\n"
        )

        handle.write(
            f"- Local topology files parsed: "
            f"{len(topology_files)}.\n"
        )

        handle.write(
            f"- Molecule types parsed: "
            f"{len(molecule_types)}.\n"
        )

        handle.write(
            f"- Duplicate molecule types: "
            f"{len(duplicate_molecule_types)}.\n"
        )

        handle.write(
            f"- Missing molecule definitions: "
            f"{len(missing_molecule_types)}.\n"
        )

        handle.write(
            f"- `gmx grompp` return code: "
            f"{grompp_result.returncode}.\n"
        )

        handle.write(
            f"- Rebuilt TPR atoms: "
            f"{rebuilt_tpr_natoms}.\n\n"
        )

        handle.write(
            "## Validation\n\n"
        )

        for key, value in (
            (
                "Expected GRO composition",
                gro_expected_composition_pass,
            ),
            (
                "GRO versus accepted TPR atom count",
                gro_tpr_atom_count_pass,
            ),
            (
                "Candidate topology total atom count",
                topology_total_atom_count_pass,
            ),
            (
                "Topology molecule counts versus GRO",
                topology_molecule_counts_pass,
            ),
            (
                "Topology atom totals versus GRO",
                topology_atom_counts_pass,
            ),
            (
                "No duplicate molecule types",
                no_duplicate_molecule_types_pass,
            ),
            (
                "No missing molecule types",
                no_missing_molecule_types_pass,
            ),
            (
                "Candidate grompp reconstruction",
                grompp_pass,
            ),
            (
                "Rebuilt TPR atom count",
                rebuilt_atom_count_pass,
            ),
        ):
            handle.write(
                f"- {key}: "
                f"{'PASS' if value else 'FAIL'}.\n"
            )

        handle.write(
            "\n## Interpretation\n\n"
        )

        if overall_pass:
            handle.write(
                "The candidate topology is "
                "compositionally and structurally "
                "consistent with the accepted GRO "
                "and TPR and can be advanced to the "
                "parameter-level TPR identity audit.\n"
            )
        else:
            handle.write(
                "The candidate topology must not be "
                "used for mobile-solute MD until the "
                "failed checks are resolved.\n"
            )

        handle.write(
            "\nThe `gmx check` comparison is stored "
            "as a diagnostic because differences in "
            "coordinates, velocities, run state, or "
            "generated metadata do not by themselves "
            "prove a topology mismatch.\n"
        )

    required_outputs = (
        FILE_INVENTORY_CSV,
        INCLUDES_CSV,
        MOLECULE_TYPES_CSV,
        MOLECULES_CSV,
        GRO_COMPOSITION_CSV,
        SUMMARY_CSV,
        REPORT_MD,
        GROMPP_LOG,
        ACCEPTED_TPR_DUMP,
        REBUILT_TPR_DUMP,
        TPR_COMPARISON_LOG,
    )

    missing_outputs = [
        path
        for path in required_outputs
        if (
            not path.exists()
            or path.stat().st_size == 0
        )
    ]

    if missing_outputs:
        raise RuntimeError(
            "Missing or empty audit outputs:\n"
            + "\n".join(
                str(path)
                for path in missing_outputs
            )
        )

    log(
        "Day021 accepted hydrated-topology "
        "audit completed."
    )

    log(
        f"GRO atoms: "
        f"{gro['atom_count']}"
    )

    log(
        f"Accepted TPR atoms: "
        f"{accepted_tpr_natoms}"
    )

    log(
        f"Candidate topology atoms: "
        f"{topology_total_atoms}"
    )

    log(
        f"Rebuilt TPR atoms: "
        f"{rebuilt_tpr_natoms}"
    )

    for molecule_row in molecule_rows:
        log(
            f"{molecule_row['molecule_type']}: "
            f"{molecule_row['molecule_count']} molecules, "
            f"{molecule_row['atoms_per_molecule']} "
            f"atoms/molecule, "
            f"{molecule_row['topology_atom_total']} atoms"
        )

    log(
        f"Local topology files: "
        f"{len(topology_files)}"
    )

    log(
        f"Unresolved includes: "
        f"{len(unresolved_local_includes)}"
    )

    log(
        f"Duplicate molecule types: "
        f"{len(duplicate_molecule_types)}"
    )

    log(
        f"Missing molecule types: "
        f"{len(missing_molecule_types)}"
    )

    log(
        f"grompp return code: "
        f"{grompp_result.returncode}"
    )

    log(
        "Structural reconstruction: "
        f"{'PASS' if overall_pass else 'FAIL'}"
    )

    log(
        f"Wrote: "
        f"{relative(OUTPUT_ROOT)}"
    )


if __name__ == "__main__":
    main()
