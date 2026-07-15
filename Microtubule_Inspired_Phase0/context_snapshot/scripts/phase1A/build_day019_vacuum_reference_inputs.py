#!/usr/bin/env python3
from __future__ import annotations

import csv
import difflib
import hashlib
import re
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

SOURCE_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/orca_embedding_pilot_inputs"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/day019_vacuum_reference_inputs"
)

MANIFEST_CSV = OUTPUT_ROOT / "VACUUM_REFERENCE_INPUT_MANIFEST_DAY019.csv"
REPORT_MD = OUTPUT_ROOT / "VACUUM_REFERENCE_INPUT_AUDIT_DAY019.md"

CHROMOPHORES = ["PYR2", "PYR3", "PYR4", "PYR5"]
EXPECTED_METHOD_LINE = "! wB97X-D3 def2-SVP def2/J RIJCOSX TightSCF"
EXPECTED_CHARGE = 0
EXPECTED_MULTIPLICITY = 1
EXPECTED_ATOMS = 26
EXPECTED_C = 16
EXPECTED_H = 10


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def parse_xyz_block(lines: list[str]) -> tuple[int, int, int, int, list[str]]:
    start = None
    stop = None
    charge = None
    multiplicity = None

    for index, line in enumerate(lines):
        match = re.match(
            r"^\s*\*\s+xyz\s+(-?\d+)\s+(\d+)\s*$",
            line,
            flags=re.IGNORECASE,
        )
        if match:
            start = index
            charge = int(match.group(1))
            multiplicity = int(match.group(2))
            continue

        if start is not None and re.match(r"^\s*\*\s*$", line):
            stop = index
            break

    if start is None or stop is None:
        raise RuntimeError("Could not locate complete '* xyz ... *' block.")

    atom_lines = [line for line in lines[start + 1 : stop] if line.strip()]
    return start, stop, charge, multiplicity, atom_lines


def validate_source_input(text: str, source_path: Path) -> dict:
    lines = text.splitlines()
    if not lines:
        raise RuntimeError(f"Empty input: {source_path}")

    method_line = lines[0].strip()
    if method_line != EXPECTED_METHOD_LINE:
        raise RuntimeError(
            f"Unexpected method line in {source_path}:\n"
            f"Observed: {method_line}\n"
            f"Expected: {EXPECTED_METHOD_LINE}"
        )

    if not re.search(r"(?mi)^\s*%maxcore\s+4096\s*$", text):
        raise RuntimeError(f"Missing or altered %maxcore in {source_path}")

    if not re.search(r"(?mi)^\s*%tddft\s*$", text):
        raise RuntimeError(f"Missing %tddft block in {source_path}")

    if not re.search(r"(?mi)^\s*nroots\s+10\s*$", text):
        raise RuntimeError(f"Missing 'nroots 10' in {source_path}")

    if not re.search(r"(?mi)^\s*tda\s+true\s*$", text):
        raise RuntimeError(f"Missing 'tda true' in {source_path}")

    pointcharge_lines = [
        line
        for line in lines
        if re.match(r"^\s*%pointcharges\b", line, flags=re.IGNORECASE)
    ]
    if len(pointcharge_lines) != 1:
        raise RuntimeError(
            f"Expected exactly one %pointcharges line in {source_path}; "
            f"found {len(pointcharge_lines)}."
        )

    start, stop, charge, multiplicity, atom_lines = parse_xyz_block(lines)
    elements = [line.split()[0] for line in atom_lines]

    if len(atom_lines) != EXPECTED_ATOMS:
        raise RuntimeError(
            f"{source_path}: expected {EXPECTED_ATOMS} atoms, "
            f"found {len(atom_lines)}."
        )

    if elements.count("C") != EXPECTED_C or elements.count("H") != EXPECTED_H:
        raise RuntimeError(
            f"{source_path}: expected C{EXPECTED_C}H{EXPECTED_H}, "
            f"found C{elements.count('C')}H{elements.count('H')}."
        )

    if charge != EXPECTED_CHARGE or multiplicity != EXPECTED_MULTIPLICITY:
        raise RuntimeError(
            f"{source_path}: expected charge/multiplicity "
            f"{EXPECTED_CHARGE}/{EXPECTED_MULTIPLICITY}, "
            f"found {charge}/{multiplicity}."
        )

    return {
        "lines": lines,
        "pointcharge_line": pointcharge_lines[0],
        "xyz_start": start,
        "xyz_stop": stop,
        "charge": charge,
        "multiplicity": multiplicity,
        "atom_lines": atom_lines,
        "elements": elements,
    }


def build_vacuum_input(source_text: str) -> str:
    output_lines = []

    for line in source_text.splitlines():
        if re.match(r"^\s*%pointcharges\b", line, flags=re.IGNORECASE):
            continue
        output_lines.append(line)

    output_text = "\n".join(output_lines).rstrip() + "\n"

    if re.search(r"(?i)%pointcharges", output_text):
        raise RuntimeError("Vacuum input still contains %pointcharges.")

    return output_text


def write_xyz(path: Path, atom_lines: list[str], comment: str) -> None:
    content = (
        f"{len(atom_lines)}\n"
        f"{comment}\n"
        + "\n".join(atom_lines)
        + "\n"
    )
    path.write_text(content, encoding="utf-8")


def validate_vacuum_input(text: str, target_path: Path) -> dict:
    if re.search(r"(?i)%pointcharges", text):
        raise RuntimeError(f"%pointcharges remains in {target_path}")

    lines = text.splitlines()

    if lines[0].strip() != EXPECTED_METHOD_LINE:
        raise RuntimeError(f"Method line changed in {target_path}")

    if not re.search(r"(?mi)^\s*nroots\s+10\s*$", text):
        raise RuntimeError(f"nroots changed in {target_path}")

    if not re.search(r"(?mi)^\s*tda\s+true\s*$", text):
        raise RuntimeError(f"TDA setting changed in {target_path}")

    start, stop, charge, multiplicity, atom_lines = parse_xyz_block(lines)
    elements = [line.split()[0] for line in atom_lines]

    if len(atom_lines) != EXPECTED_ATOMS:
        raise RuntimeError(
            f"{target_path}: expected {EXPECTED_ATOMS} atoms, "
            f"found {len(atom_lines)}."
        )

    if elements.count("C") != EXPECTED_C or elements.count("H") != EXPECTED_H:
        raise RuntimeError(
            f"{target_path}: unexpected composition "
            f"C{elements.count('C')}H{elements.count('H')}."
        )

    if charge != EXPECTED_CHARGE or multiplicity != EXPECTED_MULTIPLICITY:
        raise RuntimeError(
            f"{target_path}: unexpected charge/multiplicity "
            f"{charge}/{multiplicity}."
        )

    return {
        "atom_lines": atom_lines,
        "elements": elements,
        "charge": charge,
        "multiplicity": multiplicity,
        "xyz_start": start,
        "xyz_stop": stop,
    }


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    rows = []

    for chromophore in CHROMOPHORES:
        source_job = f"frame000_{chromophore}_embedding"
        target_job = f"frame000_{chromophore}_vacuum_reference"

        source_dir = SOURCE_ROOT / source_job
        source_input = source_dir / f"{source_job}.inp"

        if not source_input.is_file():
            raise SystemExit(f"Missing source input: {source_input}")

        source_text = source_input.read_text(encoding="utf-8")
        source_info = validate_source_input(source_text, source_input)

        target_dir = OUTPUT_ROOT / target_job
        target_dir.mkdir(parents=True, exist_ok=True)

        target_input = target_dir / f"{target_job}.inp"
        target_xyz = target_dir / f"{target_job}_qm.xyz"
        diff_path = target_dir / f"{target_job}_vs_embedding.diff"

        vacuum_text = build_vacuum_input(source_text)
        target_input.write_text(vacuum_text, encoding="utf-8")

        target_info = validate_vacuum_input(vacuum_text, target_input)

        if source_info["atom_lines"] != target_info["atom_lines"]:
            raise RuntimeError(
                f"Geometry changed while generating {target_input}"
            )

        write_xyz(
            target_xyz,
            target_info["atom_lines"],
            (
                f"{target_job}; frozen geometry copied from "
                f"{source_job}; no electrostatic point charges"
            ),
        )

        diff_text = "".join(
            difflib.unified_diff(
                source_text.splitlines(keepends=True),
                vacuum_text.splitlines(keepends=True),
                fromfile=str(source_input.relative_to(PROJECT_ROOT)),
                tofile=str(target_input.relative_to(PROJECT_ROOT)),
            )
        )
        diff_path.write_text(diff_text, encoding="utf-8")

        changed_lines = [
            line
            for line in diff_text.splitlines()
            if (
                (line.startswith("+") or line.startswith("-"))
                and not line.startswith("+++")
                and not line.startswith("---")
            )
        ]

        expected_removed = f'-{source_info["pointcharge_line"]}'
        if changed_lines != [expected_removed]:
            raise RuntimeError(
                f"Unexpected input differences for {chromophore}:\n"
                + "\n".join(changed_lines)
            )

        rows.append(
            {
                "chromophore": chromophore,
                "source_job": source_job,
                "target_job": target_job,
                "source_input": str(source_input.relative_to(PROJECT_ROOT)),
                "target_input": str(target_input.relative_to(PROJECT_ROOT)),
                "target_xyz": str(target_xyz.relative_to(PROJECT_ROOT)),
                "diff_file": str(diff_path.relative_to(PROJECT_ROOT)),
                "method_line": EXPECTED_METHOD_LINE,
                "nroots": 10,
                "tda": True,
                "charge": target_info["charge"],
                "multiplicity": target_info["multiplicity"],
                "n_atoms": len(target_info["atom_lines"]),
                "n_C": target_info["elements"].count("C"),
                "n_H": target_info["elements"].count("H"),
                "pointcharges_present": False,
                "source_sha256": sha256_text(source_text),
                "target_sha256": sha256_text(vacuum_text),
                "only_removed_line": source_info["pointcharge_line"],
            }
        )

        print(f"Generated and validated: {target_input.relative_to(PROJECT_ROOT)}")

    with MANIFEST_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    with REPORT_MD.open("w", encoding="utf-8") as handle:
        handle.write("# Day019 vacuum-reference input audit\n\n")
        handle.write("## Scope\n\n")
        handle.write(
            "Four vacuum-reference ORCA inputs were generated from the "
            "frame000 embedded inputs for PYR2âPYR5.\n\n"
        )
        handle.write("## Preserved settings\n\n")
        handle.write(f"- Method line: `{EXPECTED_METHOD_LINE}`\n")
        handle.write("- `%maxcore 4096`\n")
        handle.write("- `nroots 10`\n")
        handle.write("- `tda true`\n")
        handle.write("- Charge/multiplicity: `0 1`\n")
        handle.write("- Geometry: exactly 26 atoms, C16H10\n\n")
        handle.write("## Deliberate modification\n\n")
        handle.write(
            "The only removed input line is the corresponding "
            "`%pointcharges` directive. No coordinates or electronic-"
            "structure settings were changed.\n\n"
        )
        handle.write("## Generated jobs\n\n")
        handle.write("| Chromophore | Vacuum-reference job |\n")
        handle.write("|---|---|\n")
        for row in rows:
            handle.write(
                f"| {row['chromophore']} | `{row['target_job']}` |\n"
            )
        handle.write("\n## Validation status\n\n")
        handle.write("- Inputs generated: 4/4\n")
        handle.write("- Geometry identity checks: 4/4 passed\n")
        handle.write("- Composition checks: 4/4 passed\n")
        handle.write("- Charge/multiplicity checks: 4/4 passed\n")
        handle.write("- Residual `%pointcharges` directives: 0\n")
        handle.write("- Unexpected input differences: 0\n")

    print(f"Wrote manifest: {MANIFEST_CSV.relative_to(PROJECT_ROOT)}")
    print(f"Wrote report: {REPORT_MD.relative_to(PROJECT_ROOT)}")
    print("Day019 vacuum-reference input generation: PASS")


if __name__ == "__main__":
    main()
