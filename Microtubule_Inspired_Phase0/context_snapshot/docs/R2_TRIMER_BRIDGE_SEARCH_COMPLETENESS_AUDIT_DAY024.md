# R2 Trimer-Bridge Search-Completeness Audit

## Scope

This gate tests whether the Gate 3K angle failures result from the
trimer topology itself or from the restricted equal-angle conformer
library previously used.

The three internal trimer angles are varied independently. Candidate
coordinates are diagnostic only and are not applied to the accepted
structure.

## Search space

- Independent angle values:
  **7**
- Torsion combinations:
  **144**
- Total library conformers:
  **49392**
- Base conformers retained per failed path:
  **600**
- Azimuths per mirror:
  **72**
- Mirrors:
  **2**

## Original failed paths

- Total:
  **12**
- Lower/upper:
  **6/6**

## Expanded-search results

- `LOWER:BRIDGE:01`: passing candidates=0; best angle range=67.169–142.270°; best max bond deviation=0.000206 nm; best clashes=0; solved=False
- `LOWER:BRIDGE:02`: passing candidates=0; best angle range=67.231–142.391°; best max bond deviation=0.000233 nm; best clashes=0; solved=False
- `LOWER:BRIDGE:06`: passing candidates=0; best angle range=66.295–142.967°; best max bond deviation=0.000120 nm; best clashes=0; solved=False
- `LOWER:BRIDGE:07`: passing candidates=0; best angle range=66.791–142.648°; best max bond deviation=0.000188 nm; best clashes=0; solved=False
- `LOWER:BRIDGE:11`: passing candidates=0; best angle range=66.614–142.663°; best max bond deviation=0.000100 nm; best clashes=0; solved=False
- `LOWER:BRIDGE:12`: passing candidates=0; best angle range=66.158–143.133°; best max bond deviation=0.000138 nm; best clashes=0; solved=False
- `UPPER:BRIDGE:02`: passing candidates=0; best angle range=66.439–142.886°; best max bond deviation=0.000188 nm; best clashes=0; solved=False
- `UPPER:BRIDGE:03`: passing candidates=0; best angle range=66.765–142.650°; best max bond deviation=0.000120 nm; best clashes=0; solved=False
- `UPPER:BRIDGE:07`: passing candidates=0; best angle range=67.703–141.929°; best max bond deviation=0.000233 nm; best clashes=0; solved=False
- `UPPER:BRIDGE:08`: passing candidates=0; best angle range=67.264–142.178°; best max bond deviation=0.000206 nm; best clashes=0; solved=False
- `UPPER:BRIDGE:12`: passing candidates=0; best angle range=66.413–142.963°; best max bond deviation=0.000138 nm; best clashes=0; solved=False
- `UPPER:BRIDGE:13`: passing candidates=0; best angle range=66.270–142.896°; best max bond deviation=0.000100 nm; best clashes=0; solved=False

## Aggregate result

- Locally solved paths:
  **0**
- Locally unsolved paths:
  **12**
- Lower solved:
  **0/6**
- Upper solved:
  **0/6**

## Audit gates

- `Gate3I_graph_is_accepted`: **PASS**
- `Gate3K_has_expected_topology_review_decision`: **PASS**
- `all_Gate3K_angle_failed_paths_were_identified`: **PASS**
- `independent_angle_library_is_nonempty`: **PASS**
- `every_failed_path_received_expanded_search`: **PASS**
- `all_reported_search_metrics_are_finite`: **PASS**
- `candidate_coordinates_are_search_only`: **PASS**

## Decision

- Decision:
  **R2_TRIMER_BRIDGE_TOPOLOGY_REDESIGN_CONFIRMED_BY_EXPANDED_CONFORMER_SEARCH**
- Failed audit-integrity gates:
  **NONE**
- Current trimer graph retained:
  **NO**
- Global coordinate refinement authorized:
  **NO**
- Candidate coordinates applied:
  **NO**
- Molecular topology generation authorized:
  **NO**
- Energy minimization authorized:
  **NO**
- MD authorized:
  **NO**
- QM authorized:
  **NO**
- Required next step:
  `SCREEN_R2_SPARSE_AND_LONGER_BRIDGE_TOPOLOGIES`
