#!/usr/bin/env bash
set -u
set -o pipefail

cd "$(dirname "$0")"

ORCA="/Users/alejandro/projects/orca_6_1_1_macosx_intel_openmpi411/orca"

echo "$$" > v4.supervisor_pid

START_UTC="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf '%s\n' "$START_UTC" > v4.start_utc

"$ORCA" v4.inp > v4.out 2> v4.stderr &
ORCA_PID="$!"

echo "$ORCA_PID" > v4.orca_pid

wait "$ORCA_PID"
STATUS="$?"

END_UTC="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

printf '%s\n' "$STATUS" > v4.exit_status
printf '%s\n' "$END_UTC" > v4.end_utc

exit "$STATUS"
