rows=[]

role_counter=Counter()

node_counter=Counter()

artificial_counter=Counter()

transfer_counter=Counter()

with open(TRANSFER) as f:

    reader=csv.DictReader(f)

    for row in reader:

        if row["element"]!="H":

            continue

        rows.append(row)

        role_counter[row["atom_role"]]+=1

        node_counter[row["node_type"]]+=1

        artificial_counter[row["artificial_cap"]]+=1

        transfer_counter[row["transfer_status"]]+=1

print("[1] TOTAL HYDROGENS")

print(len(rows))

print()

print("[2] ROLE DISTRIBUTION")

for k,v in sorted(role_counter.items()):

    print(f"{v:3d}   {k}")

print()

print("[3] NODE TYPES")

for k,v in sorted(node_counter.items()):

    print(f"{v:3d}   {k}")

print()

print("[4] ARTIFICIAL CAP")

for k,v in sorted(artificial_counter.items()):

    print(f"{k:8s} {v}")

print()

print("[5] TRANSFER STATUS")

for k,v in sorted(transfer_counter.items()):

    print(f"{v:3d}   {k}")

print()
