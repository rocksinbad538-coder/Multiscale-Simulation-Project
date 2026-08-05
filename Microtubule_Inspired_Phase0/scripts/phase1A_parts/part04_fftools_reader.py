sections = {

    "ATOMS",

    "BONDS",

    "ANGLES",

    "IMPROPER",

    "DIHEDRALS",

}

for line in FF_FILE.read_text().splitlines():

    line = line.strip()

    if not line:

        continue

    if line.startswith("#"):

        continue

    upper = line.upper()

    if upper in sections:

        current_section = upper

        continue

    if current_section is None:

        continue

    entries.append(

        FFEntry(

            section=current_section,

            tokens=line.split(),

            raw=line,

        )

    )

print("[2] TOTAL ENTRIES")

print(len(entries))

print()
