rows=[]

counts={

"FUNCTIONAL_EDGE":0,

"ORIGINAL_FRAGMENT":0,

"RESP_TEMPORARY":0,

}

with open(INPUT) as f:

    reader=csv.DictReader(f)

    for row in reader:

        if row["element"]!="H":

            continue

        ff=row["proposed_forcefield_type"]

        role=row["atom_role"].upper()

        if ff in ("HB","HN"):

            taxonomy="FUNCTIONAL_EDGE"

        elif "ORIGINAL_FRAGMENT" in role:

            taxonomy="ORIGINAL_FRAGMENT"

        else:

            taxonomy="RESP_TEMPORARY"

        row["hydrogen_taxonomy"]=taxonomy

        counts[taxonomy]+=1

        rows.append(row)

print("[1] HYDROGEN TAXONOMY")

print()

for k,v in counts.items():

    print(f"{k:20s} {v}")

print()

print("total =",len(rows))

print()
