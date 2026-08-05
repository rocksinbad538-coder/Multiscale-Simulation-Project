csv_file = RUN_DIR / "FFTOOLS_PARAMETER_INVENTORY.csv"

with open(

    csv_file,

    "w",

    newline="",

) as f:

    w = csv.writer(f)

    w.writerow(

        [

            "section",

            "raw",

        ]

    )

    for e in entries:

        w.writerow(

            [

                e.section,

                e.raw,

            ]

        )

report = {

    "timestamp": utc(),

    "sections": summary,

    "entry_count": len(entries),

    "source_sha256": sha256(FF_FILE),

}

json_file = RUN_DIR / "FFTOOLS_PARAMETER_AUDIT.json"

json_file.write_text(

    json.dumps(

        report,

        indent=2,

    )

)

print("[4] OUTPUTS")

print(csv_file)

print(json_file)

print()

print("DECISION")

print(

    "AUTHORITATIVE_FORCE_FIELD_PARSED"

)
