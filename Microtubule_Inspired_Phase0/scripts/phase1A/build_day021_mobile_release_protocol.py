#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import os
import shutil
import subprocess
from collections import OrderedDict
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

SOURCE_ROOT = (
    PROJECT_ROOT
    / "parameters/phase1A/accepted/"
    "hybrid_hbnBonded_kang2000_improperGeo100_validated"
)

SOURCE_TOP = (
    SOURCE_ROOT
    / "hbn_pyrene_4_hydratable_gap45_pyr5shift_clean032_"
    "hbnBonded_kang2000_improperGeo100.top"
)

SOURCE_HBN_ITP = (
    SOURCE_ROOT
    / "hbn_bonded_candidate_kang2000_improperGeo100.itp"
)

SOURCE_PYR_ITP = (
    SOURCE_ROOT
    / "pyrene.itp"
)

SOURCE_WATER_ITP = (
    SOURCE_ROOT
    / "tip4p2005.itp"
)

ACCEPTED_RUN_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/accepted/"
    "hybrid_hydrated_gap45_pyr5shift_clean032_"
    "nvt_100ps_frozenSolute"
)

START_GRO = (
    ACCEPTED_RUN_ROOT
    / "nvt_100ps_frozenSolute.gro"
)

BASE_MDP = (
    ACCEPTED_RUN_ROOT
    / "nvt_100ps_frozenSolute.mdp"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol"
)

PROTOCOL_ROOT = (
    OUTPUT_ROOT
    / "protocol_inputs"
)

TOPOLOGY_ROOT = (
    PROTOCOL_ROOT
    / "topology"
)

MDP_ROOT = (
    PROTOCOL_ROOT
    / "mdp"
)

STATIC_ROOT = (
    OUTPUT_ROOT
    / "static_validation"
)

DERIVED_TOP = (
    TOPOLOGY_ROOT
    / "hbn_pyrene_mobile_release.top"
)

DERIVED_HBN_ITP = (
    TOPOLOGY_ROOT
    / "hbn_bonded_mobile_release.itp"
)

DERIVED_PYR_ITP = (
    TOPOLOGY_ROOT
    / "pyrene_mobile_release.itp"
)

DERIVED_WATER_ITP = (
    TOPOLOGY_ROOT
    / "tip4p2005.itp"
)

INDEX_NDX = (
    PROTOCOL_ROOT
    / "mobile_release_index.ndx"
)

PLAN_CSV = (
    PROTOCOL_ROOT
    / "mobile_release_protocol_plan.csv"
)

MANIFEST_CSV = (
    PROTOCOL_ROOT
    / "mobile_release_protocol_manifest.csv"
)

STATIC_SUMMARY_CSV = (
    STATIC_ROOT
    / "mobile_release_static_validation.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "MOBILE_RELEASE_PROTOCOL_DAY021.md"
)

TOTAL_ATOMS = 68320

GROUPS = OrderedDict(
    [
        ("System", (1, 68320)),
        ("HBN", (1, 1680)),
        ("PYR", (1681, 1784)),
        ("HBN_PYR", (1, 1784)),
        ("SOL", (1785, 68320)),
    ]
)

VELOCITY_SEED = 20260706

RESTRAINT_STRENGTHS = (
    100000,
    10000,
    1000,
    100,
)

STAGES = (
    {
        "name": "00_em_k100000",
        "kind": "em",
        "restraint": 100000,
        "nsteps": 50000,
        "dt": "",
        "time_ps": "",
        "continuation": "no",
        "gen_vel": "no",
    },
    {
        "name": "01_em_k10000",
        "kind": "em",
        "restraint": 10000,
        "nsteps": 50000,
        "dt": "",
        "time_ps": "",
        "continuation": "no",
        "gen_vel": "no",
    },
    {
        "name": "02_nvt_k10000_1ps",
        "kind": "nvt",
        "restraint": 10000,
        "nsteps": 4000,
        "dt": "0.00025",
        "time_ps": "1.0",
        "continuation": "no",
        "gen_vel": "yes",
    },
    {
        "name": "03_nvt_k1000_2ps",
        "kind": "nvt",
        "restraint": 1000,
        "nsteps": 8000,
        "dt": "0.00025",
        "time_ps": "2.0",
        "continuation": "yes",
        "gen_vel": "no",
    },
    {
        "name": "04_nvt_k100_2ps",
        "kind": "nvt",
        "restraint": 100,
        "nsteps": 4000,
        "dt": "0.0005",
        "time_ps": "2.0",
        "continuation": "yes",
        "gen_vel": "no",
    },
    {
        "name": "05_nvt_unrestrained_2ps",
        "kind": "nvt",
        "restraint": 0,
        "nsteps": 4000,
        "dt": "0.0005",
        "time_ps": "2.0",
        "continuation": "yes",
        "gen_vel": "no",
    },
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


def require_inputs() -> None:
    required = (
        SOURCE_TOP,
        SOURCE_HBN_ITP,
        SOURCE_PYR_ITP,
        SOURCE_WATER_ITP,
        START_GRO,
        BASE_MDP,
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


def parse_mdp(
    path: Path,
) -> OrderedDict[str, str]:
    settings: OrderedDict[
        str,
        str
    ] = OrderedDict()

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    for raw_line in lines:
        line = raw_line.split(
            ";",
            1,
        )[0].strip()

        if not line or "=" not in line:
            continue

        key, value = line.split(
            "=",
            1,
        )

        normalized_key = (
            key.strip()
            .lower()
            .replace("_", "-")
        )

        settings[
            normalized_key
        ] = value.strip()

    return settings


def restraint_define(
    force_constant: int,
) -> str:
    if force_constant <= 0:
        return ""

    return (
        f"-DPOSRES_HBN_K{force_constant} "
        f"-DPOSRES_PYR_K{force_constant}"
    )


def append_restraint_blocks(
    source_path: Path,
    destination_path: Path,
    atom_count: int,
    macro_prefix: str,
) -> None:
    source_text = source_path.read_text(
        encoding="utf-8",
        errors="replace",
    ).rstrip()

    output_lines = [
        source_text,
        "",
        "; ============================================================",
        "; Day021 staged position restraints for mobile release",
        "; Accepted source file preserved without modification",
        "; ============================================================",
    ]

    for force_constant in RESTRAINT_STRENGTHS:
        macro = (
            f"{macro_prefix}_K{force_constant}"
        )

        output_lines.extend(
            [
                "",
                f"#ifdef {macro}",
                "[ position_restraints ]",
                "; ai  funct       fcx       fcy       fcz",
            ]
        )

        for atom_index in range(
            1,
            atom_count + 1,
        ):
            output_lines.append(
                f"{atom_index:6d}"
                f"{1:7d}"
                f"{force_constant:12.1f}"
                f"{force_constant:12.1f}"
                f"{force_constant:12.1f}"
            )

        output_lines.append("#endif")

    destination_path.write_text(
        "\n".join(output_lines)
        + "\n",
        encoding="utf-8",
    )


def create_derived_topology() -> None:
    source_text = SOURCE_TOP.read_text(
        encoding="utf-8",
        errors="replace",
    )

    old_hbn_include = (
        '#include '
        '"hbn_bonded_candidate_'
        'kang2000_improperGeo100.itp"'
    )

    new_hbn_include = (
        '#include '
        '"hbn_bonded_mobile_release.itp"'
    )

    old_pyr_include = (
        '#include "pyrene.itp"'
    )

    new_pyr_include = (
        '#include "pyrene_mobile_release.itp"'
    )

    if old_hbn_include not in source_text:
        raise RuntimeError(
            "Expected HBN include was not found "
            "in the accepted source topology"
        )

    if old_pyr_include not in source_text:
        raise RuntimeError(
            "Expected PYR include was not found "
            "in the accepted source topology"
        )

    derived_text = source_text.replace(
        old_hbn_include,
        new_hbn_include,
        1,
    )

    derived_text = derived_text.replace(
        old_pyr_include,
        new_pyr_include,
        1,
    )

    DERIVED_TOP.write_text(
        derived_text,
        encoding="utf-8",
    )

    append_restraint_blocks(
        SOURCE_HBN_ITP,
        DERIVED_HBN_ITP,
        atom_count=1680,
        macro_prefix="POSRES_HBN",
    )

    append_restraint_blocks(
        SOURCE_PYR_ITP,
        DERIVED_PYR_ITP,
        atom_count=26,
        macro_prefix="POSRES_PYR",
    )

    shutil.copy2(
        SOURCE_WATER_ITP,
        DERIVED_WATER_ITP,
    )


def write_index_group(
    handle,
    name: str,
    first_atom: int,
    last_atom: int,
) -> None:
    handle.write(
        f"[ {name} ]\n"
    )

    current_line: list[str] = []

    for atom_index in range(
        first_atom,
        last_atom + 1,
    ):
        current_line.append(
            str(atom_index)
        )

        if len(current_line) == 15:
            handle.write(
                " ".join(current_line)
                + "\n"
            )

            current_line = []

    if current_line:
        handle.write(
            " ".join(current_line)
            + "\n"
        )

    handle.write("\n")


def create_index() -> None:
    with INDEX_NDX.open(
        "w",
        encoding="utf-8",
    ) as handle:
        for (
            group_name,
            atom_range,
        ) in GROUPS.items():
            write_index_group(
                handle,
                group_name,
                atom_range[0],
                atom_range[1],
            )


def count_index_group(
    path: Path,
    group_name: str,
) -> int:
    active_group = ""
    count = 0

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    for raw_line in lines:
        line = raw_line.strip()

        if (
            line.startswith("[")
            and line.endswith("]")
        ):
            active_group = (
                line[1:-1].strip()
            )

            continue

        if (
            active_group != group_name
            or not line
        ):
            continue

        count += len(
            line.split()
        )

    return count


def render_stage_mdp(
    stage: dict[str, object],
) -> Path:
    settings = parse_mdp(
        BASE_MDP
    )

    remove_keys = {
        "freezegrps",
        "freezedim",
        "define",
        "ld-seed",
        "gen-seed",
        "gen-temp",
    }

    for key in remove_keys:
        settings.pop(
            key,
            None,
        )

    force_constant = int(
        stage["restraint"]
    )

    define_value = restraint_define(
        force_constant
    )

    if define_value:
        settings["define"] = (
            define_value
        )

    settings["nsteps"] = str(
        stage["nsteps"]
    )

    settings["continuation"] = str(
        stage["continuation"]
    )

    settings["gen-vel"] = str(
        stage["gen_vel"]
    )

    settings["constraints"] = "none"
    settings["pcoupl"] = "no"

    if stage["kind"] == "em":
        settings["integrator"] = "steep"
        settings.pop("dt", None)

        settings["emtol"] = "1000.0"
        settings["emstep"] = "0.0005"

        settings["tcoupl"] = "no"
        settings["nstxout"] = "0"
        settings["nstvout"] = "0"
        settings["nstfout"] = "0"
        settings[
            "nstxout-compressed"
        ] = "0"

        settings["nstenergy"] = "100"
        settings["nstlog"] = "100"

    else:
        settings["integrator"] = "md"
        settings["dt"] = str(
            stage["dt"]
        )

        if int(stage["restraint"]) > 0:
            settings["comm-mode"] = "None"
            settings.pop(
                "comm-grps",
                None,
            )
            settings.pop(
                "nstcomm",
                None,
            )
        else:
            settings["comm-mode"] = "Linear"
            settings["comm-grps"] = "System"
            settings["nstcomm"] = "100"

        settings["tcoupl"] = "V-rescale"
        settings["tc-grps"] = "System"
        settings["tau-t"] = "1.0"
        settings["ref-t"] = "300"

        settings["nstxout"] = "0"
        settings["nstvout"] = "0"
        settings["nstfout"] = "0"

        if str(stage["dt"]) == "0.00025":
            output_stride = "400"
        else:
            output_stride = "200"

        settings[
            "nstxout-compressed"
        ] = output_stride

        settings["nstenergy"] = (
            output_stride
        )

        settings["nstlog"] = (
            output_stride
        )

        settings[
            "nstcalcenergy"
        ] = "100"

        if stage["gen_vel"] == "yes":
            settings["gen-temp"] = "300"
            settings["gen-seed"] = str(
                VELOCITY_SEED
            )

    destination = (
        MDP_ROOT
        / f"{stage['name']}.mdp"
    )

    with destination.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "; Day021 mobile-release protocol\n"
        )

        handle.write(
            f"; Derived from {relative(BASE_MDP)}\n"
        )

        handle.write(
            "; Accepted source files were not modified\n\n"
        )

        for key, value in settings.items():
            handle.write(
                f"{key:<28s} = {value}\n"
            )

    return destination


def count_position_restraints(
    processed_topology: Path,
) -> int:
    current_section = ""
    count = 0

    lines = processed_topology.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    for raw_line in lines:
        line = raw_line.split(
            ";",
            1,
        )[0].strip()

        if not line:
            continue

        if (
            line.startswith("[")
            and line.endswith("]")
        ):
            current_section = (
                line[1:-1]
                .strip()
                .lower()
            )

            continue

        if (
            current_section
            != "position_restraints"
        ):
            continue

        fields = line.split()

        try:
            int(fields[0])
        except (
            IndexError,
            ValueError,
        ):
            continue

        count += 1

    return count


def run_static_grompp(
    gmx: Path,
    stage: dict[str, object],
) -> dict[str, object]:
    stage_name = str(
        stage["name"]
    )

    stage_root = (
        STATIC_ROOT
        / stage_name
    )

    stage_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    mdp_path = (
        MDP_ROOT
        / f"{stage_name}.mdp"
    )

    tpr_path = (
        stage_root
        / f"{stage_name}.tpr"
    )

    processed_topology = (
        stage_root
        / f"{stage_name}_processed.top"
    )

    mdout_path = (
        stage_root
        / f"{stage_name}_mdout.mdp"
    )

    log_path = (
        stage_root
        / f"{stage_name}_grompp.log"
    )

    for path in (
        tpr_path,
        processed_topology,
        mdout_path,
    ):
        if path.exists():
            path.unlink()

    command = [
        str(gmx),
        "grompp",
        "-f",
        str(mdp_path),
        "-c",
        str(START_GRO),
        "-r",
        str(START_GRO),
        "-p",
        str(DERIVED_TOP),
        "-n",
        str(INDEX_NDX),
        "-o",
        str(tpr_path),
        "-po",
        str(mdout_path),
        "-pp",
        str(processed_topology),
        "-maxwarn",
        "0",
    ]

    environment = {
        **os.environ,
        "GMX_MAXBACKUP": "-1",
    }

    result = subprocess.run(
        command,
        cwd=TOPOLOGY_ROOT,
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

    outputs_exist = all(
        path.exists()
        and path.stat().st_size > 0
        for path in (
            tpr_path,
            processed_topology,
            mdout_path,
        )
    )

    restraint_count = (
        count_position_restraints(
            processed_topology
        )
        if processed_topology.exists()
        else -1
    )

    expected_restraints = (
        1706
        if int(stage["restraint"]) > 0
        else 0
    )

    grompp_pass = (
        result.returncode == 0
        and outputs_exist
    )

    restraint_pass = (
        restraint_count
        == expected_restraints
    )

    return {
        "stage": stage_name,
        "kind": stage["kind"],
        "restraint_k": stage["restraint"],
        "grompp_return_code": (
            result.returncode
        ),
        "grompp_pass": grompp_pass,
        "position_restraint_entries": (
            restraint_count
        ),
        "expected_position_restraint_entries": (
            expected_restraints
        ),
        "restraint_entry_validation": (
            restraint_pass
        ),
        "tpr_path": relative(
            tpr_path
        ),
        "processed_topology_path": (
            relative(
                processed_topology
            )
        ),
        "grompp_log": relative(
            log_path
        ),
    }


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
        for key in row:
            if key in seen:
                continue

            seen.add(key)
            fieldnames.append(key)

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

        for row in rows:
            writer.writerow(
                {
                    key: row.get(
                        key,
                        "",
                    )
                    for key in fieldnames
                }
            )


def create_plan() -> None:
    rows: list[
        dict[str, object]
    ] = []

    for stage in STAGES:
        force_constant = int(
            stage["restraint"]
        )

        rows.append(
            {
                **stage,
                "define": (
                    restraint_define(
                        force_constant
                    )
                ),
                "HBN_restrained_atoms": (
                    1680
                    if force_constant > 0
                    else 0
                ),
                "PYR_restrained_atoms_per_molecule": (
                    26
                    if force_constant > 0
                    else 0
                ),
                "PYR_molecule_count": 4,
                "velocity_policy": (
                    "regenerate_all_velocities"
                    if stage["gen_vel"] == "yes"
                    else (
                        "not_applicable"
                        if stage["kind"] == "em"
                        else "continue_from_previous_stage"
                    )
                ),
            }
        )

    write_csv(
        PLAN_CSV,
        rows,
    )


def create_manifest() -> None:
    paths = [
        SOURCE_TOP,
        SOURCE_HBN_ITP,
        SOURCE_PYR_ITP,
        SOURCE_WATER_ITP,
        START_GRO,
        BASE_MDP,
        DERIVED_TOP,
        DERIVED_HBN_ITP,
        DERIVED_PYR_ITP,
        DERIVED_WATER_ITP,
        INDEX_NDX,
        PLAN_CSV,
    ]

    paths.extend(
        sorted(
            MDP_ROOT.glob(
                "*.mdp"
            )
        )
    )

    rows = [
        {
            "path": relative(path),
            "size_bytes": (
                path.stat().st_size
            ),
            "sha256": sha256(path),
        }
        for path in paths
    ]

    write_csv(
        MANIFEST_CSV,
        rows,
    )


def main() -> None:
    require_inputs()

    TOPOLOGY_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    MDP_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    STATIC_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    create_derived_topology()
    create_index()

    for stage in STAGES:
        render_stage_mdp(stage)

    create_plan()

    index_counts = {
        group_name: count_index_group(
            INDEX_NDX,
            group_name,
        )
        for group_name in GROUPS
    }

    expected_index_counts = {
        group_name: (
            atom_range[1]
            - atom_range[0]
            + 1
        )
        for (
            group_name,
            atom_range,
        ) in GROUPS.items()
    }

    index_validation = (
        index_counts
        == expected_index_counts
    )

    gmx = find_gmx()

    static_rows = [
        run_static_grompp(
            gmx,
            stage,
        )
        for stage in STAGES
    ]

    write_csv(
        STATIC_SUMMARY_CSV,
        static_rows,
    )

    create_manifest()

    all_grompp_pass = all(
        bool(
            row["grompp_pass"]
        )
        for row in static_rows
    )

    all_restraint_counts_pass = all(
        bool(
            row[
                "restraint_entry_validation"
            ]
        )
        for row in static_rows
    )

    protocol_static_validation = all(
        (
            index_validation,
            all_grompp_pass,
            all_restraint_counts_pass,
        )
    )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day021 Mobile-Release Protocol\n\n"
        )

        handle.write(
            "## Provenance\n\n"
        )

        handle.write(
            f"- Source topology: "
            f"`{relative(SOURCE_TOP)}`\n"
        )

        handle.write(
            f"- Starting GRO: "
            f"`{relative(START_GRO)}`\n"
        )

        handle.write(
            "- Accepted source files modified: no.\n\n"
        )

        handle.write(
            "## Index groups\n\n"
        )

        for (
            group_name,
            count,
        ) in index_counts.items():
            handle.write(
                f"- {group_name}: "
                f"{count} atoms.\n"
            )

        handle.write(
            f"- Index validation: "
            f"{'PASS' if index_validation else 'FAIL'}.\n\n"
        )

        handle.write(
            "## Static GROMACS validation\n\n"
        )

        for row in static_rows:
            handle.write(
                f"- {row['stage']}: "
                f"grompp="
                f"{'PASS' if row['grompp_pass'] else 'FAIL'}, "
                f"position restraints="
                f"{row['position_restraint_entries']}/"
                f"{row['expected_position_restraint_entries']}.\n"
            )

        handle.write(
            "\n## Decision\n\n"
        )

        handle.write(
            "- Protocol static validation: "
            f"{'PASS' if protocol_static_validation else 'FAIL'}.\n"
        )

        handle.write(
            "- Scientific calculation started: no.\n"
        )

    log(
        "Day021 mobile-release protocol "
        "construction completed."
    )

    log(
        "Accepted source files modified: NO"
    )

    log(
        f"Derived topology: "
        f"{relative(DERIVED_TOP)}"
    )

    log("Index counts:")

    for (
        group_name,
        count,
    ) in index_counts.items():
        log(
            f"  {group_name}: {count}"
        )

    log(
        f"Index validation: "
        f"{'PASS' if index_validation else 'FAIL'}"
    )

    log(
        "Static stage validation:"
    )

    for row in static_rows:
        log(
            f"  {row['stage']}: "
            f"grompp="
            f"{'PASS' if row['grompp_pass'] else 'FAIL'}, "
            f"posres="
            f"{row['position_restraint_entries']}/"
            f"{row['expected_position_restraint_entries']}"
        )

    log(
        f"All grompp checks: "
        f"{'PASS' if all_grompp_pass else 'FAIL'}"
    )

    log(
        f"All restraint-count checks: "
        f"{'PASS' if all_restraint_counts_pass else 'FAIL'}"
    )

    log(
        "Protocol static validation: "
        f"{'PASS' if protocol_static_validation else 'FAIL'}"
    )

    log(
        "Scientific calculation started: NO"
    )

    log(
        f"Wrote: {relative(OUTPUT_ROOT)}"
    )


if __name__ == "__main__":
    main()
