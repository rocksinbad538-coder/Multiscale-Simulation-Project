csv_file=RUN/"HYDROGEN_TAXONOMY.csv"

with open(csv_file,"w",newline="") as f:

    writer=csv.DictWriter(

        f,

        fieldnames=rows[0].keys()

    )

    writer.writeheader()

    writer.writerows(rows)

report={

"timestamp":utc(),

"functional_edge":counts["FUNCTIONAL_EDGE"],

"original_fragment":counts["ORIGINAL_FRAGMENT"],

"resp_temporary":counts["RESP_TEMPORARY"],

"decision":

"D040_A10_HYDROGEN_TAXONOMY_COMPLETE",

"phase1A_scientific_taxonomy_complete":True

}

json_file=RUN/"HYDROGEN_TAXONOMY.json"

json_file.write_text(

json.dumps(

report,

indent=2,

)

)

print("[2] OUTPUTS")

print(csv_file)

print(json_file)

print()

print("[3] DECISION")

print(report["decision"])
