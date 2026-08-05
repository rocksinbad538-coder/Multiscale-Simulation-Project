report={

"timestamp":utc(),

"total_hydrogen_rows":len(hydrogen_ids),

"unique_hydrogen_atom_ids":len(counter),

"duplicate_count":len(duplicates),

"duplicates":duplicates,

"decision":

"D040_A10A_HYDROGEN_MAPPING_INTEGRITY_COMPLETE"

}

json_file=RUN/"HYDROGEN_MAPPING_INTEGRITY.json"

json_file.write_text(

json.dumps(

report,

indent=2,

)

)

print("[3] OUTPUT")

print(json_file)

print()

print("[4] DECISION")

print(report["decision"])
