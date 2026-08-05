rows=[]

hydrogen_ids=[]

with open(INPUT) as f:

    reader=csv.DictReader(f)

    for row in reader:

        rows.append(row)

        if row["element"]=="H":

            hydrogen_ids.append(row["atom_id"])

counter=Counter(hydrogen_ids)

duplicates={

k:v

for k,v in counter.items()

if v>1

}

print("[1] HYDROGEN COUNTS")

print()

print("total_H_rows =",len(hydrogen_ids))

print("unique_H =",len(counter))

print("duplicate_atom_ids =",len(duplicates))

print()

if duplicates:

    print("[2] DUPLICATED HYDROGENS")

    print()

    for atom_id,n in sorted(duplicates.items()):

        print(atom_id,"count=",n)

print()
