#!/usr/bin/env python3
from pathlib import Path

src = Path("parameters/phase1A/accepted/hybrid_dry_gromacs")
water = Path("parameters/phase1A/accepted/water_tip4p2005")
out = Path("parameters/phase1A/hybrid_hydrated_gromacs")
out.mkdir(parents=True, exist_ok=True)

# Copy coordinate and component ITP files
(out/"hbn_pyrene_4_dry_for_solvation.gro").write_text((src/"hbn_pyrene_4_dry.gro").read_text())
(out/"hbn_fixed_dummy.itp").write_text((src/"hbn_fixed_dummy.itp").read_text())
(out/"pyrene.itp").write_text((src/"pyrene.itp").read_text())
(out/"tip4p2005.itp").write_text((water/"tip4p2005.itp").read_text())
(out/"tip4p2005_atomtypes.itp").write_text((water/"tip4p2005_atomtypes.itp").read_text())

top = out/"hbn_pyrene_4_hydratable.top"

top.write_text("""; Hydratable h-BN + 4 pyrene topology, Phase 1A.9b
; h-BN currently dummy/fixed/non-interacting: not physical h-BN FF

[ defaults ]
; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ
1 2 yes 0.5 0.83333333

[ atomtypes ]
; name at.num mass charge ptype sigma epsilon
B0 5 10.810000 0.000000 A 0.000000 0.000000
N0 7 14.007000 0.000000 A 0.000000 0.000000

; GAFF2 atomtypes used by pyrene
ca 6 12.010000 0.000000 A 0.33152123 0.4133792
ha 1 1.008000 0.000000 A 0.26254785 0.0673624

; TIP4P/2005 atomtypes
OW 8 15.9994 0.0000 A 0.315890 0.774900
HW 1 1.0080 0.0000 A 0.000000 0.000000
MW 0 0.0000 0.0000 D 0.000000 0.000000

#include "hbn_fixed_dummy.itp"
#include "pyrene.itp"
#include "tip4p2005.itp"

[ system ]
h-BN fixed dummy scaffold + 4 pyrenes + TIP4P/2005 water

[ molecules ]
; molname count
HBN 1
PYR 4
""")

print("Wrote", top)
print("Prepared hydration directory:", out)
