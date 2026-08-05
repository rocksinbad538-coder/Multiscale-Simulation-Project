csv_file=RUN/"PHASE1B_TOPOLOGY_MAPPING.csv"

with open(csv_file,"w",newline="") as f:

    writer=csv.DictWriter(
        f,
        fieldnames=rows[0].keys()
    )

    writer.writeheader()

    writer.writerows(rows)

report={

"timestamp":utc(),

"total_atoms":len(rows),

"topology_atoms":keep,

"mapping_complete":True,

"decision":
"D040_A9_TOPOLOGY_MAPPING_COMPLETE",

"phase1B_topology_ready":True,

}

json_file=RUN/"PHASE1B_TOPOLOGY_MAPPING.json"

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
