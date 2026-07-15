#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import re
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/"
    "day021_accepted_hydrated_topology_audit"
)

OUTPUT_CSV = (
    OUTPUT_ROOT
    / "hydrated_topology_candidate_inventory.csv"
)

OUTPUT_REPORT = (
    OUTPUT_ROOT
    / "HYDRATED_TOPOLOGY_CANDIDATE_INVENTORY_DAY021.md"
)

EXCLUDED_DIRECTORY_NAMES = {
    ".git",
    "__pycache__",
}

EXCLUDED_PATH_PREFIXES = (
    OUTPUT_ROOT.resolve(),
)

EXPECTED_MOLECULES = {
    "HBN": 1,
    "PYR": 4,
    "SOL": 16634,
}

SECTION_PATTERN = re.compile(
    r"^\s*\[\s*([^\]]+)\s*\]\s*$"
)

INCLUDE_PATTERN = re.compile(
    r'^\s*#include\s+[<"]([^>"]+)[>"]'
)

HBN_SECTIONS = (
    "atoms",
    "bonds",
    "angles",
    "dihedrals",
    "pairs",
    "pairs_nb",
    "constraints",
    "exclusions",
    "settles",
    "virtual_sites2",
    "virtual_sites3",
    "virtual_sites4",
)


def log(message: str = "") -> None:
    print(message, flush=True)


def relative(path: Path) -> str:
    try:
        return str(
            path.resolve().relative_to(
                PROJECT_ROOT
            )
        )
    except ValueError:
        return str(path.resolve())


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


def strip_comment(line: str) -> str:
    return line.split(
        ";",
        1,
    )[0].strip()


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
) -> None:
    if not rows:
        rows = [
            {
                "status": (
                    "NO_MATCHING_TOPOLOGIES_FOUND"
                )
            }
        ]

    fieldnames: list[str] = []
    seen: set[str] = set()

    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)

    normalized = [
        {
            key: row.get(key, "")
            for key in fieldnames
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


def is_excluded(path: Path) -> bool:
    if any(
        part in EXCLUDED_DIRECTORY_NAMES
        for part in path.parts
    ):
        return True

    resolved = path.resolve()

    return any(
        resolved == prefix
        or prefix in resolved.parents
        for prefix in EXCLUDED_PATH_PREFIXES
    )


def find_topology_files() -> list[Path]:
    result = []

    for path in PROJECT_ROOT.rglob("*.top"):
        if not path.is_file():
            continue

        if is_excluded(path):
            continue

        result.append(path.resolve())

    return sorted(
        result,
        key=lambda path: relative(path)
    )


def parse_molecules(
    topology_path: Path,
) -> dict[str, int]:
    section = ""
    molecules: dict[str, int] = {}

    lines = topology_path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    for raw_line in lines:
        line = strip_comment(raw_line)

        if not line:
            continue

        if line.startswith("#"):
            continue

        section_match = (
            SECTION_PATTERN.match(line)
        )

        if section_match:
            section = (
                section_match.group(1)
                .strip()
                .lower()
            )

            continue

        if section != "molecules":
            continue

        fields = line.split()

        if len(fields) < 2:
            continue

        try:
            count = int(fields[1])
        except ValueError:
            continue

        molecule_name = fields[0]

        molecules[molecule_name] = (
            molecules.get(
                molecule_name,
                0,
            )
            + count
        )

    return molecules


def resolve_include(
    parent: Path,
    include_token: str,
) -> Path | None:
    candidates = (
        parent.parent / include_token,
        PROJECT_ROOT / include_token,
    )

    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()

    return None


def collect_local_include_closure(
    root_topology: Path,
) -> tuple[list[Path], list[str]]:
    visited: set[Path] = set()
    ordered: list[Path] = []
    unresolved: list[str] = []

    def recurse(path: Path) -> None:
        resolved = path.resolve()

        if resolved in visited:
            return

        visited.add(resolved)
        ordered.append(resolved)

        lines = resolved.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines()

        for line in lines:
            match = INCLUDE_PATTERN.match(
                line
            )

            if match is None:
                continue

            token = match.group(1)

            included = resolve_include(
                resolved,
                token,
            )

            if included is None:
                unresolved.append(
                    f"{relative(resolved)}::{token}"
                )

                continue

            recurse(included)

    recurse(root_topology)

    return ordered, unresolved


def parse_hbn_definitions(
    topology_files: list[Path],
) -> list[dict[str, object]]:
    definitions: list[
        dict[str, object]
    ] = []

    for path in topology_files:
        lines = path.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines()

        section = ""
        current_molecule: str | None = None
        waiting_for_molecule_name = False
        current_definition: (
            dict[str, object] | None
        ) = None

        for line_number, raw_line in enumerate(
            lines,
            start=1,
        ):
            line = strip_comment(raw_line)

            if not line:
                continue

            if line.startswith("#"):
                continue

            section_match = (
                SECTION_PATTERN.match(line)
            )

            if section_match:
                section = (
                    section_match.group(1)
                    .strip()
                    .lower()
                )

                if section == "moleculetype":
                    waiting_for_molecule_name = (
                        True
                    )

                    current_molecule = None
                    current_definition = None

                continue

            fields = line.split()

            if (
                section == "moleculetype"
                and waiting_for_molecule_name
            ):
                current_molecule = fields[0]

                waiting_for_molecule_name = (
                    False
                )

                if current_molecule == "HBN":
                    current_definition = {
                        "source_file": (
                            relative(path)
                        ),
                        "definition_line": (
                            line_number
                        ),
                        **{
                            f"HBN_{name}_entries": 0
                            for name in HBN_SECTIONS
                        },
                    }

                    definitions.append(
                        current_definition
                    )

                continue

            if (
                current_molecule != "HBN"
                or current_definition is None
            ):
                continue

            if section not in HBN_SECTIONS:
                continue

            try:
                int(fields[0])
            except (
                IndexError,
                ValueError,
            ):
                continue

            key = (
                f"HBN_{section}_entries"
            )

            current_definition[key] = (
                int(
                    current_definition[key]
                )
                + 1
            )

    return definitions


def exact_molecular_composition(
    molecules: dict[str, int],
) -> bool:
    for name, count in (
        EXPECTED_MOLECULES.items()
    ):
        if molecules.get(name, 0) != count:
            return False

    extra_nonzero = {
        name: count
        for name, count
        in molecules.items()
        if (
            count != 0
            and name
            not in EXPECTED_MOLECULES
        )
    }

    return not extra_nonzero


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    topology_files = (
        find_topology_files()
    )

    exact_candidates: list[
        dict[str, object]
    ] = []

    composition_candidates = 0

    for top_path in topology_files:
        molecules = parse_molecules(
            top_path
        )

        if not molecules:
            continue

        if not exact_molecular_composition(
            molecules
        ):
            continue

        composition_candidates += 1

        include_files, unresolved = (
            collect_local_include_closure(
                top_path
            )
        )

        hbn_definitions = (
            parse_hbn_definitions(
                include_files
            )
        )

        if not hbn_definitions:
            exact_candidates.append(
                {
                    "topology_path": (
                        relative(top_path)
                    ),
                    "topology_sha256": (
                        sha256(top_path)
                    ),
                    "HBN_count": (
                        molecules.get(
                            "HBN",
                            0,
                        )
                    ),
                    "PYR_count": (
                        molecules.get(
                            "PYR",
                            0,
                        )
                    ),
                    "SOL_count": (
                        molecules.get(
                            "SOL",
                            0,
                        )
                    ),
                    "local_include_file_count": (
                        len(include_files)
                    ),
                    "unresolved_include_count": (
                        len(unresolved)
                    ),
                    "HBN_definition_found": (
                        False
                    ),
                    "HBN_definition_index": "",
                    "likely_original_frozen_HBN": (
                        False
                    ),
                }
            )

            continue

        for definition_index, definition in enumerate(
            hbn_definitions,
            start=1,
        ):
            bond_count = int(
                definition[
                    "HBN_bonds_entries"
                ]
            )

            angle_count = int(
                definition[
                    "HBN_angles_entries"
                ]
            )

            dihedral_count = int(
                definition[
                    "HBN_dihedrals_entries"
                ]
            )

            likely_frozen = (
                int(
                    definition[
                        "HBN_atoms_entries"
                    ]
                )
                == 1680
                and bond_count == 0
                and angle_count == 0
                and dihedral_count == 0
            )

            exact_candidates.append(
                {
                    "topology_path": (
                        relative(top_path)
                    ),
                    "topology_sha256": (
                        sha256(top_path)
                    ),
                    "HBN_count": (
                        molecules["HBN"]
                    ),
                    "PYR_count": (
                        molecules["PYR"]
                    ),
                    "SOL_count": (
                        molecules["SOL"]
                    ),
                    "local_include_file_count": (
                        len(include_files)
                    ),
                    "unresolved_include_count": (
                        len(unresolved)
                    ),
                    "HBN_definition_found": (
                        True
                    ),
                    "HBN_definition_index": (
                        definition_index
                    ),
                    **definition,
                    "likely_original_frozen_HBN": (
                        likely_frozen
                    ),
                }
            )

    write_csv(
        OUTPUT_CSV,
        exact_candidates,
    )

    likely_originals = [
        row
        for row in exact_candidates
        if bool(
            row.get(
                "likely_original_frozen_HBN",
                False,
            )
        )
    ]

    with OUTPUT_REPORT.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day021 Hydrated Topology "
            "Candidate Inventory\n\n"
        )

        handle.write(
            f"- `.top` files inspected: "
            f"{len(topology_files)}.\n"
        )

        handle.write(
            f"- Exact HBN/PYR/SOL composition "
            f"candidates: "
            f"{composition_candidates}.\n"
        )

        handle.write(
            f"- Candidate HBN definitions recorded: "
            f"{len(exact_candidates)}.\n"
        )

        handle.write(
            f"- Unbonded 1680-atom HBN candidates: "
            f"{len(likely_originals)}.\n\n"
        )

        handle.write(
            "## Candidate table\n\n"
        )

        if not exact_candidates:
            handle.write(
                "No exact current-working-tree "
                "candidate was found.\n"
            )

        for index, row in enumerate(
            exact_candidates,
            start=1,
        ):
            handle.write(
                f"### Candidate {index}\n\n"
            )

            handle.write(
                f"- Topology: "
                f"`{row['topology_path']}`\n"
            )

            handle.write(
                f"- HBN definition: "
                f"`{row.get('source_file', 'not found')}`\n"
            )

            handle.write(
                f"- HBN atoms: "
                f"{row.get('HBN_atoms_entries', '')}\n"
            )

            handle.write(
                f"- HBN bonds: "
                f"{row.get('HBN_bonds_entries', '')}\n"
            )

            handle.write(
                f"- HBN angles: "
                f"{row.get('HBN_angles_entries', '')}\n"
            )

            handle.write(
                f"- HBN dihedrals: "
                f"{row.get('HBN_dihedrals_entries', '')}\n"
            )

            handle.write(
                f"- Likely original frozen HBN: "
                f"{row.get('likely_original_frozen_HBN', False)}\n\n"
            )

    log(
        "Day021 hydrated-topology "
        "candidate inventory completed."
    )

    log(
        f"TOP files inspected: "
        f"{len(topology_files)}"
    )

    log(
        f"Exact composition candidates: "
        f"{composition_candidates}"
    )

    log(
        f"HBN definitions recorded: "
        f"{len(exact_candidates)}"
    )

    log(
        f"Likely original frozen-HBN candidates: "
        f"{len(likely_originals)}"
    )

    for index, row in enumerate(
        exact_candidates,
        start=1,
    ):
        log(
            f"Candidate {index}: "
            f"{row['topology_path']}"
        )

        log(
            "  HBN source: "
            f"{row.get('source_file', 'NOT FOUND')}"
        )

        log(
            "  HBN atoms/bonds/angles/dihedrals: "
            f"{row.get('HBN_atoms_entries', '')}/"
            f"{row.get('HBN_bonds_entries', '')}/"
            f"{row.get('HBN_angles_entries', '')}/"
            f"{row.get('HBN_dihedrals_entries', '')}"
        )

        log(
            "  likely original frozen HBN: "
            f"{row.get('likely_original_frozen_HBN', False)}"
        )

    log(
        f"Wrote: "
        f"{relative(OUTPUT_CSV)}"
    )

    log(
        f"Wrote: "
        f"{relative(OUTPUT_REPORT)}"
    )


if __name__ == "__main__":
    main()
