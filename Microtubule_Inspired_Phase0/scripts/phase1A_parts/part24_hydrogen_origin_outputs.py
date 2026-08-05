csv_file=RUN/"RESP_HYDROGEN_ORIGIN_TABLE.csv"

with open(csv_file,"w",newline="") as f:

    writer=csv.DictWriter(
        f,
        fieldnames=rows[0].keys()
    )

    writer.writeheader()

    writer.writerows(rows)

report={

"timestamp":utc(),

"hydrogen_count":len(rows),

"role_distribution":dict(role_counter),

"node_distribution":dict(node_counter),

"artificial_cap_distribution":dict(artificial_counter),

"transfer_status_distribution":dict(transfer_counter),

"decision":"D040_A10C_RESP_HYDROGEN_ORIGIN_AUDIT_COMPLETE"

}

json_file=RUN/"RESP_HYDROGEN_ORIGIN_AUDIT.json"

json_file.write_text(

json.dumps(

report,

indent=2,

)

)

print("[6] OUTPUTS")

print(csv_file)

print(json_file)

print()

print("[7] DECISION")

print(report["decision"])
