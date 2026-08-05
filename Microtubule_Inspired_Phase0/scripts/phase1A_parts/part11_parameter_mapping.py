mapping = []

parent_lookup = {}

with open(PARENT_MAPPING) as f:

    reader = csv.DictReader(f)

    for row in reader:

        parent_lookup[
            row["hydrogen_atom_id"]
        ] = row

with open(TRANSFERABILITY) as f:

    reader = csv.DictReader(f)

    for row in reader:

        atom_id = row["atom_id"]

        element = row["element"]

        role = row["atom_role"]

        parent_element = ""

        proposed_type = element

        if atom_id in parent_lookup:

            parent_element = parent_lookup[
                atom_id
            ]["parent_element"]

            if parent_element == "B":

                proposed_type = "HB"

            elif parent_element == "N":

                proposed_type = "HN"

        mapping.append(

            {

                "atom_id": atom_id,

                "element": element,

                "atom_role": role,

                "parent_element": parent_element,

                "proposed_forcefield_type": proposed_type,

                "RESP_stage1_charge_e": row[
                    "RESP_stage1_charge_e_float"
                ],

            }

        )

print("[1] TRACEABILITY TABLE")

print("rows =", len(mapping))

print()

counts = {}

for m in mapping:

    t = m["proposed_forcefield_type"]

    counts[t] = counts.get(t, 0) + 1

print("[2] PROPOSED FORCE-FIELD TYPES")

for k in sorted(counts):

    print(f"{k:>4s} : {counts[k]}")

print()
