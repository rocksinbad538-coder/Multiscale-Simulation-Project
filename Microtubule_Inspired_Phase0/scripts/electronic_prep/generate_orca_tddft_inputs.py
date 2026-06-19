from pathlib import Path

SETS = [
    (
        Path("runs/phase1A/day012_electronic_prep/qm_extracts_pyr_only"),
        Path("runs/phase1A/day012_electronic_prep/orca_inputs/pyr_only_tddft")
    ),
    (
        Path("runs/phase1A/day012_electronic_prep/qm_extracts"),
        Path("runs/phase1A/day012_electronic_prep/orca_inputs/pyr_water0p50_tddft")
    ),
]

HEADER = """! PBE0 def2-SVP D3BJ RIJCOSX TightSCF TDDFT

%pal
  nprocs 4
end

%tddft
  nroots 10
  maxdim 50
end

%scf
  MaxIter 300
end

* xyz 0 1
"""

VALID_QM_ELEMENTS = {"H", "C", "N", "O", "B"}

def clean_coord_line(line):
    parts = line.split()
    if len(parts) != 4:
        return None

    elem = parts[0]

    # Exclude TIP4P/2005 virtual site.
    if elem in {"M", "MW", "X"}:
        return None

    if elem not in VALID_QM_ELEMENTS:
        raise ValueError(f"Unexpected QM element '{elem}' in line: {line}")

    return line

def main():
    for srcdir, outdir in SETS:
        outdir.mkdir(parents=True, exist_ok=True)

        for xyz in sorted(srcdir.glob("*.xyz")):
            lines = xyz.read_text().splitlines()
            coords_raw = lines[2:]

            coords = []
            skipped = 0

            for line in coords_raw:
                cleaned = clean_coord_line(line)
                if cleaned is None:
                    skipped += 1
                    continue
                coords.append(cleaned)

            inp = outdir / (xyz.stem + ".inp")

            with inp.open("w") as f:
                f.write(HEADER)
                for line in coords:
                    f.write(line + "\n")
                f.write("*\n")

            print(f"{inp}: qm_atoms={len(coords)} skipped_virtual_sites={skipped}")

if __name__ == "__main__":
    main()
