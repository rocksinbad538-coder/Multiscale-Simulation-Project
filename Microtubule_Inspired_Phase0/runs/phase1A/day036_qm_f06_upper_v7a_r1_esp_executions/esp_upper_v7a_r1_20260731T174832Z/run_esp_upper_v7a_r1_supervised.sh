#!/usr/bin/env bash
set -u

ORCA_BIN='/Users/alejandro/projects/orca_6_1_1_macosx_intel_openmpi411/orca'
INPUT_FILE='esp_upper_v7a_r1.inp'
OUTPUT_FILE='esp_upper_v7a_r1.out'
STDERR_FILE='esp_upper_v7a_r1.stderr'
STATUS_FILE='esp_upper_v7a_r1.exit_status'
PID_FILE='esp_upper_v7a_r1.orca_pid'
SUPERVISOR_LOG='esp_upper_v7a_r1.supervisor.log'

{
    echo "Start UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "Working directory: $PWD"
    echo "ORCA binary: $ORCA_BIN"
    echo "Input: $INPUT_FILE"
} > "$SUPERVISOR_LOG"

"$ORCA_BIN" "$INPUT_FILE"     > "$OUTPUT_FILE"     2> "$STDERR_FILE" &

ORCA_PID=$!
printf '%s\n' "$ORCA_PID" > "$PID_FILE"

echo "ORCA PID: $ORCA_PID" >> "$SUPERVISOR_LOG"

wait "$ORCA_PID"
STATUS=$?

printf '%s\n' "$STATUS" > "$STATUS_FILE"

{
    echo "Exit status: $STATUS"
    echo "End UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
} >> "$SUPERVISOR_LOG"

exit "$STATUS"
