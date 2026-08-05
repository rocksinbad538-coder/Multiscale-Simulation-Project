rows=[]

with open(INPUT) as f:

    reader=csv.DictReader(f)

    for r in reader:

        keep=False

        reason=""

        ff=r["proposed_forcefield_type"]

        role=r["atom_role"]

        if ff in ("HB","HN"):

            keep=True
            reason="FUNCTIONAL_EDGE"

        elif r["element"]!="H":

            keep=True
            reason="FRAMEWORK"

        else:

            keep=False
            reason="NON_FUNCTIONAL_H"

        r["topology_keep"]=keep
        r["topology_reason"]=reason

        rows.append(r)

print("[1] TOPOLOGY CLASSIFICATION")

print("atoms =",len(rows))

print()

keep=sum(
r["topology_keep"]
for r in rows
)

print("kept =",keep)

print("discarded =",len(rows)-keep)

print()
