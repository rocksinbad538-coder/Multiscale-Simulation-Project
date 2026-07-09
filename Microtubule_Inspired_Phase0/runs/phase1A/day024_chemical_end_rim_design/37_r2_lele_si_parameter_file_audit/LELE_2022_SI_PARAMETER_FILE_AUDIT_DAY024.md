# Lele 2022 SI and ReaxFF Parameter-File Audit

## Verified artifacts

- Primary article:
  `runs/phase1A/day024_chemical_end_rim_design/37_r2_lele_si_parameter_file_audit/raw/jp1c09648.pdf`
- Supporting Information:
  `runs/phase1A/day024_chemical_end_rim_design/37_r2_lele_si_parameter_file_audit/raw/jp1c09648_si_001.pdf`
- ReaxFF parameter file:
  `runs/phase1A/day024_chemical_end_rim_design/37_r2_lele_si_parameter_file_audit/raw/jp1c09648_si_002.txt`

All three artifacts match the expected SHA-256 hashes.

## Supporting Information

- Pages: **2**
- Extracted characters: **2365**
- Quality-factor evidence rows:
  **8**
- Borazine-snapshot evidence rows:
  **6**
- Parameter-table hits:
  **0**
- Water hits:
  **0**

The SI contains the quality-factor derivation and the 4 ns
borazine simulation snapshot. It does not supply a separate
R2-relevant validation set.

## ReaxFF parameter-file structure

- General parameters:
  **39**
- Elements:
  **6**
- Bonds:
  **21**
- Off-diagonal terms:
  **15**
- Angles:
  **105**
- Torsions:
  **58**
- Hydrogen bonds:
  **4**
- Element order:
  **C | H | O | N | B | Al**
- B/N/H-relevant parameter records:
  **166**

## Scientific conclusion

The file is a genuine composite ReaxFF parameter set containing
C, H, O, N, B and Al. It includes reactive records involving
B-N, B-H and N-H chemistry. Its demonstrated target domain is
gas-phase B/N/H chemistry and high-temperature BN nanostructure
formation.

It does not establish validated transferability to:

- equilibrium mechanics of the selected R2 BNNT;
- reconstructed annulus environments;
- four-atom B-N-B-N bridges;
- confined water;
- anisotropic scaffold-water polarization.

## Decision

- Decision:
  **R2_LELE_2022_SI_AND_REAXFF_PARAMETER_FILE_AUDITED_NOT_AUTHORIZED_FOR_R2**
- Failed gates:
  **NONE**
- R2 force-field coverage established:
  **NO**
- Parameter adoption authorized:
  **NO**
- Topology generation authorized:
  **NO**
- Charge assignment authorized:
  **NO**
- Force-field parameterization authorized:
  **NO**
- Energy minimization authorized:
  **NO**
- MD authorized:
  **NO**
- QM calculation authorized:
  **NO**
- Required next step:
  `COMPARE_LELE_REAXFF_DOMAIN_WITH_FUNCTIONALIZED_HBN_FIXED_TOPOLOGY_AND_QM_REFERENCE_ROUTES`
