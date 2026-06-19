#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/pyr_only_tddft"

INPUT="PYR2_only.inp"
OUTPUT="PYR2_only.out"

orca "$INPUT" > "$OUTPUT"

grep -Ei "ORCA TERMINATED NORMALLY|TOTAL RUN TIME|EXCITED STATES|STATE" "$OUTPUT" || true
