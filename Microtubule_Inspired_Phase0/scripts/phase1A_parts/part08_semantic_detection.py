required={

"HB atom type":False,
"HN atom type":False,

"B-HB bond":False,
"N-HN bond":False,
"B-N bond":False,

"B-N-HN angle":False,
"N-B-HB angle":False,

"B improper":False,
"N improper":False,

"HB dihedral":False,
"HN dihedral":False,

}

evidence={}

for row in rows:

    text=row["raw"]

    tokens=text.split()

    sec=row["section"]

    if sec=="ATOMS":

        if text.startswith("HB "):

            required["HB atom type"]=True
            evidence["HB atom type"]=text

        if text.startswith("HN "):

            required["HN atom type"]=True
            evidence["HN atom type"]=text

    elif sec=="BONDS":

        if text.startswith("B  HB"):

            required["B-HB bond"]=True
            evidence["B-HB bond"]=text

        elif text.startswith("N  HN"):

            required["N-HN bond"]=True
            evidence["N-HN bond"]=text

        elif text.startswith("B  N"):

            required["B-N bond"]=True
            evidence["B-N bond"]=text

    elif sec=="ANGLES":

        if len(tokens)>=3:

            if tokens[:3]==["B","N","HN"]:

                required["B-N-HN angle"]=True
                evidence["B-N-HN angle"]=text

            elif tokens[:3]==["N","B","HB"]:

                required["N-B-HB angle"]=True
                evidence["N-B-HB angle"]=text

    elif sec=="IMPROPER":

        if text.startswith("B"):

            required["B improper"]=True
            evidence["B improper"]=text

        elif text.startswith("N"):

            required["N improper"]=True
            evidence["N improper"]=text

    elif sec=="DIHEDRALS":

        if "HN" in text:

            required["HN dihedral"]=True
            evidence["HN dihedral"]=text

        if "HB" in text:

            required["HB dihedral"]=True
            evidence["HB dihedral"]=text

print("[2] SEMANTIC AUDIT")

for k,v in required.items():

    print(f"{k:25s}",v)

print()
