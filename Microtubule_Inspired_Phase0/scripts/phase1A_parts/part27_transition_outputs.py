phys_csv=RUN/"PHASE1B_PHYSICAL_HYDROGENS.csv"

with open(phys_csv,"w",newline="") as f:

    writer=csv.DictWriter(
        f,
        fieldnames=physical[0].keys()
    )

    writer.writeheader()
    writer.writerows(physical)

caps_csv=RUN/"PHASE1B_QM_CAPS.csv"

with open(caps_csv,"w",newline="") as f:

    writer=csv.DictWriter(
        f,
        fieldnames=qm_caps[0].keys()
    )

    writer.writeheader()
    writer.writerows(qm_caps)

report={

"timestamp":utc(),

"physical_hydrogens":len(physical),

"qm_caps":len(qm_caps),

"phase1A_complete":True,

"phase1B_ready":True,

"decision":"D040_A11_PHASE1B_TRANSITION_READY"

}

json_file=RUN/"PHASE1B_TRANSITION_SPECIFICATION.json"

json_file.write_text(

json.dumps(report,indent=2)

)

print("[2] OUTPUTS")

print(phys_csv)
print(caps_csv)
print(json_file)

print()

print("[3] DECISION")

print(report["decision"])
