physical=[]

qm_caps=[]

with open(INPUT) as f:

    reader=csv.DictReader(f)

    for row in reader:

        if row["transfer_status"]=="TRANSFERABLE_REAL_ATOM":

            physical.append(row)

        else:

            qm_caps.append(row)

print("[1] TRANSITION SUMMARY")

print()

print("physical hydrogens =",len(physical))

print("QM caps =",len(qm_caps))

print()
