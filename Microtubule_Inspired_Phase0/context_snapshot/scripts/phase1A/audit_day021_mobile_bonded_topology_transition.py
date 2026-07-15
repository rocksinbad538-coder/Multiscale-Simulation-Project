#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import os
import re
import shutil
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

FROZEN_ROOT = (
    PROJECT_ROOT
    / "parameters/phase1A/accepted/"
    "hybrid_hydrated_gap45_pyr5shift_clean032"
)

FROZEN_TOP = (
    FROZEN_ROOT
    / "hbn_pyrene_4_hydratable_gap45_"
    "pyr5shift_clean032.top"
)

FROZEN_HBN = (
    FROZEN_ROOT
    / "hbn_fixed_dummy.itp"
)

FROZEN_PYR = (
    FROZEN_ROOT
    / "pyrene.itp"
)

FROZEN_WATER = (
    FROZEN_ROOT
    / "tip4p2005.itp"
)

MOBILE_ROOT = (
    PROJECT_ROOT
    / "parameters/phase1A/accepted/"
    "hybrid_hbnBonded_kang2000_"
    "improperGeo100_validated"
)

MOBILE_TOP = (
    MOBILE_ROOT
    / "hbn_pyrene_4_hydratable_gap45_"
    "pyr5shift_clean032_hbnBonded_"
    "kang2000_improperGeo100.top"
)

MOBILE_HBN = (
    MOBILE_ROOT
    / "hbn_bonded_candidate_"
    "kang2000_improperGeo100.itp"
)

MOBILE_PYR = (
    MOBILE_ROOT
    / "pyrene.itp"
)

MOBILE_WATER = (
    MOBILE_ROOT
    / "tip4p2005.itp"
)

ACCEPTED_RUN_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/accepted/"
    "hybrid_hydrated_gap45_pyr5shift_clean032_"
    "nvt_100ps_frozenSolute"
)

FINAL_GRO = (
    ACCEPTED_RUN_ROOT
    / "nvt_100ps_frozenSolute.gro"
)

AUDIT_MDP = (
    PROJECT_ROOT
    / "runs/phase1A/"
    "day021_accepted_hydrated_topology_audit/"
    "original_frozen_topology_identity/"
    "nvt_100ps_frozenSolute_accepted_seed.mdp"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/"
    "day021_accepted_hydrated_topology_audit/"
    "mobile_bonded_transition_audit"
)

REBUILT_TPR = (
    OUTPUT_ROOT
    / "mobile_bonded_candidate_static_check.tpr"
)

PROCESSED_TOP = (
    OUTPUT_ROOT
    / "mobile_bonded_candidate_processed.top"
)

MDOUT_MDP = (
    OUTPUT_ROOT
    / "mobile_bonded_candidate_mdout.mdp"
)

GROMPP_LOG = (
    OUTPUT_ROOT
    / "grompp_mobile_bonded_candidate.log"
)

TPR_DUMP = (
    OUTPUT_ROOT
    / "mobile_bonded_candidate_dump.txt"
)

TPR_DUMP_STDERR = (
    OUTPUT_ROOT
    / "mobile_bonded_candidate_dump_stderr.log"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "mobile_bonded_transition_audit_summary.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "MOBILE_BONDED_TRANSITION_AUDIT_DAY021.md"
)

EXPECTED_MOLECULES = {
    "HBN": 1,
    "PYR": 4,
    "SOL": 16634,
}

EXPECTED_HBN_COUNTS = {
    "atoms": 1680,
    "bonds": 2460,
    "angles": 4860,
    "dihedrals": 1620,
}

EXPECTED_TOTAL_ATOMS = 68320

SECTION_PATTERN = re.compile(
    r"^\s*\[\s*([^\]]+)\s*\]\s*$"
)

NATOMS_PATTERN = re.compile(
    r"\bnatoms\s*=\s*(\d+)"
)

NTYPES_PATTERN = re.compile(
    r"\bntypes\s*=\s*(\d+)"
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


def normalized_semantic_lines(
    path: Path,
) -> list[str]:
    result: list[str] = []

    for raw_line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        line = strip_comment(raw_line)

        if not line:
            continue

        result.append(
            " ".join(line.split())
        )

    return result


def parse_sections(
    path: Path,
) -> dict[str, list[list[str]]]:
    sections: dict[
        str,
        list[list[str]]
    ] = {}

    current_section = ""

    for raw_line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        line = strip_comment(raw_line)

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

            sections.setdefault(
                current_section,
                [],
            )

            continue

        if not current_section:
            continue

        sections.setdefault(
            current_section,
            [],
        ).append(
            line.split()
        )

    return sections


def first_data_row(
    sections: dict[
        str,
        list[list[str]],
    ],
    section_name: str,
) -> list[str]:
    rows = sections.get(
        section_name,
        [],
    )

    if not rows:
        raise RuntimeError(
            f"Missing [{section_name}] data"
        )

    return rows[0]


def parse_atomtypes(
    top_path: Path,
) -> dict[str, list[str]]:
    sections = parse_sections(
        top_path
    )

    result: dict[
        str,
        list[str]
    ] = {}

    for fields in sections.get(
        "atomtypes",
        [],
    ):
        if not fields:
            continue

        result[fields[0]] = fields

    return result


def parse_molecules(
    top_path: Path,
) -> dict[str, int]:
    sections = parse_sections(
        top_path
    )

    result: dict[str, int] = {}

    for fields in sections.get(
        "molecules",
        [],
    ):
        if len(fields) < 2:
            continue

        try:
            count = int(fields[1])
        except ValueError:
            continue

        result[fields[0]] = (
            result.get(
                fields[0],
                0,
            )
            + count
        )

    return result


def parse_hbn_itp(
    path: Path,
) -> dict[str, object]:
    sections = parse_sections(
        path
    )

    moleculetype = first_data_row(
        sections,
        "moleculetype",
    )

    atom_rows = [
        tuple(fields)
        for fields in sections.get(
            "atoms",
            [],
        )
        if (
            fields
            and fields[0].isdigit()
        )
    ]

    counts: dict[str, int] = {}

    for section_name in (
        "atoms",
        "bonds",
        "angles",
        "dihedrals",
        "pairs",
        "constraints",
        "exclusions",
    ):
        counts[section_name] = sum(
            1
            for fields in sections.get(
                section_name,
                [],
            )
            if (
                fields
                and fields[0].isdigit()
            )
        )

    return {
        "molecule_name": (
            moleculetype[0]
        ),
        "nrexcl": (
            int(moleculetype[1])
            if len(moleculetype) > 1
            else None
        ),
        "atom_rows": atom_rows,
        "counts": counts,
    }


def find_gmx() -> Path:
    preferred = Path(
        "/usr/local/gromacs/bin/gmx"
    )

    if preferred.exists():
        return preferred

    discovered = shutil.which("gmx")

    if discovered is None:
        raise RuntimeError(
            "Could not locate GROMACS gmx"
        )

    return Path(discovered)


def run_combined(
    command: list[str],
    cwd: Path,
    log_path: Path,
) -> subprocess.CompletedProcess[str]:
    environment = {
        **os.environ,
        "GMX_MAXBACKUP": "-1",
    }

    result = subprocess.run(
        command,
        cwd=cwd,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )

    log_path.write_text(
        result.stdout,
        encoding="utf-8",
    )

    return result


def run_dump(
    gmx: Path,
    tpr_path: Path,
) -> None:
    environment = {
        **os.environ,
        "GMX_MAXBACKUP": "-1",
    }

    result = subprocess.run(
        [
            str(gmx),
            "dump",
            "-s",
            str(tpr_path),
        ],
        cwd=PROJECT_ROOT,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )

    TPR_DUMP.write_text(
        result.stdout,
        encoding="utf-8",
    )

    TPR_DUMP_STDERR.write_text(
        result.stderr,
        encoding="utf-8",
    )

    if result.returncode != 0:
        raise RuntimeError(
            "gmx dump failed for mobile "
            "bonded candidate"
        )


def parse_first_integer(
    text: str,
    pattern: re.Pattern[str],
    label: str,
) -> int:
    match = pattern.search(text)

    if match is None:
        raise RuntimeError(
            f"Could not parse {label}"
        )

    return int(match.group(1))


def write_csv(
    path: Path,
    row: dict[str, object],
) -> None:
    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(row.keys()),
        )

        writer.writeheader()
        writer.writerow(row)


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    required = (
        FROZEN_TOP,
        FROZEN_HBN,
        FROZEN_PYR,
        FROZEN_WATER,
        MOBILE_TOP,
        MOBILE_HBN,
        MOBILE_PYR,
        MOBILE_WATER,
        FINAL_GRO,
        AUDIT_MDP,
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

    frozen_top_sections = parse_sections(
        FROZEN_TOP
    )

    mobile_top_sections = parse_sections(
        MOBILE_TOP
    )

    defaults_equal = (
        frozen_top_sections.get(
            "defaults",
            [],
        )
        == mobile_top_sections.get(
            "defaults",
            [],
        )
    )

    frozen_atomtypes = parse_atomtypes(
        FROZEN_TOP
    )

    mobile_atomtypes = parse_atomtypes(
        MOBILE_TOP
    )

    required_atomtypes = (
        "B0",
        "N0",
    )

    atomtype_equalities = {
        name: (
            frozen_atomtypes.get(name)
            == mobile_atomtypes.get(name)
        )
        for name in required_atomtypes
    }

    frozen_b0 = frozen_atomtypes.get(
        "B0"
    )

    frozen_n0 = frozen_atomtypes.get(
        "N0"
    )

    mobile_b0 = mobile_atomtypes.get(
        "B0"
    )

    mobile_n0 = mobile_atomtypes.get(
        "N0"
    )

    def nonzero_lj(
        fields: list[str] | None,
    ) -> bool:
        if (
            fields is None
            or len(fields) < 2
        ):
            return False

        try:
            sigma = float(fields[-2])
            epsilon = float(fields[-1])
        except ValueError:
            return False

        return (
            sigma > 0.0
            and epsilon > 0.0
        )

    mobile_hbn_lj_nonzero = all(
        (
            nonzero_lj(mobile_b0),
            nonzero_lj(mobile_n0),
        )
    )

    frozen_molecules = parse_molecules(
        FROZEN_TOP
    )

    mobile_molecules = parse_molecules(
        MOBILE_TOP
    )

    molecule_composition_equal = (
        frozen_molecules
        == mobile_molecules
        == EXPECTED_MOLECULES
    )

    frozen_hbn = parse_hbn_itp(
        FROZEN_HBN
    )

    mobile_hbn = parse_hbn_itp(
        MOBILE_HBN
    )

    hbn_atom_table_equal = (
        frozen_hbn["atom_rows"]
        == mobile_hbn["atom_rows"]
    )

    pyrene_semantic_equal = (
        normalized_semantic_lines(
            FROZEN_PYR
        )
        == normalized_semantic_lines(
            MOBILE_PYR
        )
    )

    water_semantic_equal = (
        normalized_semantic_lines(
            FROZEN_WATER
        )
        == normalized_semantic_lines(
            MOBILE_WATER
        )
    )

    mobile_counts = mobile_hbn[
        "counts"
    ]

    expected_bonded_counts_pass = all(
        int(
            mobile_counts.get(
                section_name,
                -1,
            )
        )
        == expected_count
        for (
            section_name,
            expected_count,
        )
        in EXPECTED_HBN_COUNTS.items()
    )

    frozen_hbn_unbonded_pass = all(
        int(
            frozen_hbn["counts"].get(
                section_name,
                -1,
            )
        )
        == expected_count
        for (
            section_name,
            expected_count,
        )
        in {
            "atoms": 1680,
            "bonds": 0,
            "angles": 0,
            "dihedrals": 0,
        }.items()
    )

    gmx = find_gmx()

    grompp_result = run_combined(
        [
            str(gmx),
            "grompp",
            "-f",
            str(AUDIT_MDP),
            "-c",
            str(FINAL_GRO),
            "-p",
            str(MOBILE_TOP),
            "-o",
            str(REBUILT_TPR),
            "-po",
            str(MDOUT_MDP),
            "-pp",
            str(PROCESSED_TOP),
            "-maxwarn",
            "0",
        ],
        MOBILE_ROOT,
        GROMPP_LOG,
    )

    grompp_pass = (
        grompp_result.returncode == 0
        and REBUILT_TPR.exists()
        and REBUILT_TPR.stat().st_size > 0
    )

    if not grompp_pass:
        raise RuntimeError(
            "Mobile bonded candidate failed "
            f"grompp. See {GROMPP_LOG}"
        )

    run_dump(
        gmx,
        REBUILT_TPR,
    )

    dump_text = TPR_DUMP.read_text(
        encoding="utf-8",
        errors="replace",
    )

    rebuilt_natoms = parse_first_integer(
        dump_text,
        NATOMS_PATTERN,
        "natoms",
    )

    rebuilt_ntypes = parse_first_integer(
        dump_text,
        NTYPES_PATTERN,
        "ntypes",
    )

    atom_count_pass = (
        rebuilt_natoms
        == EXPECTED_TOTAL_ATOMS
    )

    intended_transition_pass = all(
        (
            defaults_equal,
            all(
                atomtype_equalities.values()
            ),
            mobile_hbn_lj_nonzero,
            molecule_composition_equal,
            hbn_atom_table_equal,
            pyrene_semantic_equal,
            water_semantic_equal,
            frozen_hbn_unbonded_pass,
            expected_bonded_counts_pass,
            grompp_pass,
            atom_count_pass,
        )
    )

    summary = {
        "frozen_TOP": (
            relative(FROZEN_TOP)
        ),
        "mobile_TOP": (
            relative(MOBILE_TOP)
        ),
        "frozen_HBN_ITP": (
            relative(FROZEN_HBN)
        ),
        "mobile_HBN_ITP": (
            relative(MOBILE_HBN)
        ),
        "frozen_TOP_sha256": (
            sha256(FROZEN_TOP)
        ),
        "mobile_TOP_sha256": (
            sha256(MOBILE_TOP)
        ),
        "defaults_equal": (
            defaults_equal
        ),
        "B0_atomtype_equal": (
            atomtype_equalities["B0"]
        ),
        "N0_atomtype_equal": (
            atomtype_equalities["N0"]
        ),
        "mobile_B0_record": (
            " ".join(mobile_b0 or [])
        ),
        "mobile_N0_record": (
            " ".join(mobile_n0 or [])
        ),
        "mobile_HBN_LJ_nonzero": (
            mobile_hbn_lj_nonzero
        ),
        "molecule_composition_equal": (
            molecule_composition_equal
        ),
        "HBN_atom_table_equal": (
            hbn_atom_table_equal
        ),
        "frozen_HBN_nrexcl": (
            frozen_hbn["nrexcl"]
        ),
        "mobile_HBN_nrexcl": (
            mobile_hbn["nrexcl"]
        ),
        "frozen_HBN_atoms": (
            frozen_hbn["counts"]["atoms"]
        ),
        "frozen_HBN_bonds": (
            frozen_hbn["counts"]["bonds"]
        ),
        "frozen_HBN_angles": (
            frozen_hbn["counts"]["angles"]
        ),
        "frozen_HBN_dihedrals": (
            frozen_hbn["counts"]["dihedrals"]
        ),
        "mobile_HBN_atoms": (
            mobile_counts["atoms"]
        ),
        "mobile_HBN_bonds": (
            mobile_counts["bonds"]
        ),
        "mobile_HBN_angles": (
            mobile_counts["angles"]
        ),
        "mobile_HBN_dihedrals": (
            mobile_counts["dihedrals"]
        ),
        "mobile_HBN_pairs": (
            mobile_counts["pairs"]
        ),
        "mobile_HBN_constraints": (
            mobile_counts["constraints"]
        ),
        "mobile_HBN_exclusions": (
            mobile_counts["exclusions"]
        ),
        "pyrene_semantic_equal": (
            pyrene_semantic_equal
        ),
        "water_semantic_equal": (
            water_semantic_equal
        ),
        "grompp_return_code": (
            grompp_result.returncode
        ),
        "grompp_pass": (
            grompp_pass
        ),
        "rebuilt_TPR_atoms": (
            rebuilt_natoms
        ),
        "rebuilt_TPR_ntypes": (
            rebuilt_ntypes
        ),
        "intended_model_transition_validation": (
            "PASS"
            if intended_transition_pass
            else "FAIL"
        ),
    }

    write_csv(
        SUMMARY_CSV,
        summary,
    )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day021 Mobile Bonded-Topology "
            "Transition Audit\n\n"
        )

        handle.write(
            "## Preserved model components\n\n"
        )

        handle.write(
            f"- Defaults equal: "
            f"{defaults_equal}.\n"
        )

        handle.write(
            f"- B0 atom type equal: "
            f"{atomtype_equalities['B0']}.\n"
        )

        handle.write(
            f"- N0 atom type equal: "
            f"{atomtype_equalities['N0']}.\n"
        )

        handle.write(
            f"- Mobile HBN LJ nonzero: "
            f"{mobile_hbn_lj_nonzero}.\n"
        )

        handle.write(
            f"- HBN atom table equal: "
            f"{hbn_atom_table_equal}.\n"
        )

        handle.write(
            f"- Pyrene model equal: "
            f"{pyrene_semantic_equal}.\n"
        )

        handle.write(
            f"- TIP4P/2005 model equal: "
            f"{water_semantic_equal}.\n"
        )

        handle.write(
            f"- Molecular composition equal: "
            f"{molecule_composition_equal}.\n\n"
        )

        handle.write(
            "## Intended HBN changes\n\n"
        )

        handle.write(
            f"- Frozen HBN bonded counts: "
            f"{frozen_hbn['counts']}.\n"
        )

        handle.write(
            f"- Mobile HBN bonded counts: "
            f"{mobile_counts}.\n"
        )

        handle.write(
            f"- Expected bonded counts pass: "
            f"{expected_bonded_counts_pass}.\n\n"
        )

        handle.write(
            "## GROMACS reconstruction\n\n"
        )

        handle.write(
            f"- `grompp`: "
            f"{'PASS' if grompp_pass else 'FAIL'}.\n"
        )

        handle.write(
            f"- Rebuilt atoms: "
            f"{rebuilt_natoms}.\n"
        )

        handle.write(
            f"- Rebuilt interaction parameter "
            f"types: {rebuilt_ntypes}.\n\n"
        )

        handle.write(
            "## Decision\n\n"
        )

        handle.write(
            "- Intended frozen-to-bonded model "
            "transition: "
            f"{'PASS' if intended_transition_pass else 'FAIL'}.\n"
        )

        handle.write(
            "- This test validates model continuity "
            "and topology construction only. It does "
            "not yet authorize unrestrained production "
            "MD.\n"
        )

    log(
        "Day021 mobile bonded-topology "
        "transition audit completed."
    )

    log(
        f"Defaults preserved: "
        f"{'PASS' if defaults_equal else 'FAIL'}"
    )

    log(
        f"B0 atomtype preserved: "
        f"{'PASS' if atomtype_equalities['B0'] else 'FAIL'}"
    )

    log(
        f"N0 atomtype preserved: "
        f"{'PASS' if atomtype_equalities['N0'] else 'FAIL'}"
    )

    log(
        f"Mobile HBN LJ nonzero: "
        f"{'PASS' if mobile_hbn_lj_nonzero else 'FAIL'}"
    )

    log(
        f"HBN atom table/order preserved: "
        f"{'PASS' if hbn_atom_table_equal else 'FAIL'}"
    )

    log(
        f"PYR model preserved: "
        f"{'PASS' if pyrene_semantic_equal else 'FAIL'}"
    )

    log(
        f"TIP4P/2005 model preserved: "
        f"{'PASS' if water_semantic_equal else 'FAIL'}"
    )

    log(
        f"Molecular composition preserved: "
        f"{'PASS' if molecule_composition_equal else 'FAIL'}"
    )

    log(
        "Frozen HBN atoms/bonds/angles/dihedrals: "
        f"{frozen_hbn['counts']['atoms']}/"
        f"{frozen_hbn['counts']['bonds']}/"
        f"{frozen_hbn['counts']['angles']}/"
        f"{frozen_hbn['counts']['dihedrals']}"
    )

    log(
        "Mobile HBN atoms/bonds/angles/dihedrals: "
        f"{mobile_counts['atoms']}/"
        f"{mobile_counts['bonds']}/"
        f"{mobile_counts['angles']}/"
        f"{mobile_counts['dihedrals']}"
    )

    log(
        f"Expected mobile bonded counts: "
        f"{'PASS' if expected_bonded_counts_pass else 'FAIL'}"
    )

    log(
        f"grompp reconstruction: "
        f"{'PASS' if grompp_pass else 'FAIL'}"
    )

    log(
        f"Rebuilt TPR atoms/types: "
        f"{rebuilt_natoms}/{rebuilt_ntypes}"
    )

    log(
        "Intended model transition validation: "
        f"{'PASS' if intended_transition_pass else 'FAIL'}"
    )

    log(
        f"Wrote: {relative(OUTPUT_ROOT)}"
    )


if __name__ == "__main__":
    main()
