coverage=sum(required.values())

total=len(required)

coverage_fraction=coverage/total

print("[3] COVERAGE")

print(f"{coverage}/{total}")

print()

decision=(
    "PASS"
    if coverage==total
    else
    "REVIEW"
)

csv_file=RUN_DIR/"SEMANTIC_PARAMETER_AUDIT.csv"

with open(csv_file,"w",newline="") as f:

    w=csv.writer(f)

    w.writerow(

        [

            "requirement",

            "status",

            "evidence",

        ]

    )

    for k in required:

        w.writerow(

            [

                k,

                required[k],

                evidence.get(k,""),

            ]

        )

certificate={

"timestamp":utc(),

"coverage":coverage,

"total":total,

"coverage_fraction":coverage_fraction,

"decision":decision,

"HB_atom_type_detected":required["HB atom type"],

"HN_atom_type_detected":required["HN atom type"],

"scientific_conclusion":
"Hydrogen type assignment is explicitly separated into HB and HN according to parent atom identity in the authoritative force field.",

"phase1B_parameterization_ready":
coverage==total,

}

json_file=RUN_DIR/"PHASE1B_PARAMETER_COMPLETENESS_CERTIFICATE.json"

json_file.write_text(

    json.dumps(

        certificate,

        indent=2,

    )

)

print("[4] OUTPUTS")

print(csv_file)

print(json_file)

print()

print("[5] DECISION")

print(decision)
