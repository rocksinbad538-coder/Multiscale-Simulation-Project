csv_file=RUN/"PHASE1B_PARAMETER_MAPPING.csv"

with open(csv_file,"w",newline="") as f:

    writer=csv.DictWriter(
        f,
        fieldnames=list(mapping[0].keys())
    )

    writer.writeheader()

    writer.writerows(mapping)

report={

"timestamp":utc(),

"mapped_atoms":len(mapping),

"mapping_generated":True,

"HB_detected":any(
x["proposed_forcefield_type"]=="HB"
for x in mapping
),

"HN_detected":any(
x["proposed_forcefield_type"]=="HN"
for x in mapping
),

"decision":
"D040_A8_PARAMETER_MAPPING_COMPLETE",

"phase1B_mapping_ready":True

}

json_file=RUN/"PHASE1B_PARAMETER_MAPPING.json"

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
