summary = {}

for e in entries:

    summary.setdefault(

        e.section,

        0,

    )

    summary[e.section] += 1

print("[3] SECTION COUNTS")

for k in (

    "ATOMS",

    "BONDS",

    "ANGLES",

    "IMPROPER",

    "DIHEDRALS",

):

    print(

        k,

        summary.get(k,0),

    )

print()
