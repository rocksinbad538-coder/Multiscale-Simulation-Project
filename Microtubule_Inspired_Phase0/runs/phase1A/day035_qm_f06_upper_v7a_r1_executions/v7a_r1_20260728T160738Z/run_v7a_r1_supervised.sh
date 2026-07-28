#!/bin/bash

set -uo pipefail

EXEC_DIR='/Users/alejandro/projects/Multiscale-Simulation-Project/Microtubule_Inspired_Phase0/runs/phase1A/day035_qm_f06_upper_v7a_r1_executions/v7a_r1_20260728T160738Z'
ORCA_BIN='/Users/alejandro/projects/orca_6_1_1_macosx_intel_openmpi411/orca'
OUTPUT_FILE="v7a_r1.out"
STDERR_FILE="v7a_r1.stderr"

cd "$EXEC_DIR" || {
    printf '%s\n' 90 > v7a_r1.exit_status
    exit 90
}

printf '%s\n' "$$" > v7a_r1.supervisor_pid

{
    echo "========================================================="
    echo "QM_F06 UPPER V7-A R1 SUPERVISED EXECUTION"
    echo "========================================================="
    echo "Started: $(date)"
    echo "Execution directory: $EXEC_DIR"
    echo "ORCA binary: $ORCA_BIN"
} > v7a_r1.supervisor.log

"$ORCA_BIN" v7a_r1.inp > "$OUTPUT_FILE" 2> "$STDERR_FILE" &
ORCA_PID="$!"

printf '%s\n' "$ORCA_PID" > v7a_r1.orca_pid
echo "ORCA PID: $ORCA_PID" >> v7a_r1.supervisor.log

wait "$ORCA_PID"
SHELL_STATUS="$?"

printf '%s\n' "$SHELL_STATUS" > v7a_r1.orca_shell_status

NORMAL_COUNT="$(grep -c 'ORCA TERMINATED NORMALLY' "$OUTPUT_FILE" 2>/dev/null || true)"
ERROR_COUNT="$(grep -c 'ORCA finished by error termination' "$OUTPUT_FILE" 2>/dev/null || true)"
OPT_COUNT="$(grep -c 'THE OPTIMIZATION HAS CONVERGED' "$OUTPUT_FILE" 2>/dev/null || true)"

if [ "$ERROR_COUNT" -gt 0 ]; then
    FINAL_STATUS=91
    CLASSIFICATION="ORCA_ERROR_TERMINATION"
elif [ "$SHELL_STATUS" -ne 0 ]; then
    FINAL_STATUS="$SHELL_STATUS"
    CLASSIFICATION="NONZERO_SHELL_STATUS"
elif [ "$NORMAL_COUNT" -eq 0 ]; then
    FINAL_STATUS=92
    CLASSIFICATION="NORMAL_TERMINATION_MARKER_ABSENT"
else
    FINAL_STATUS=0
    CLASSIFICATION="ORCA_NORMAL_TERMINATION"
fi

printf '%s\n' "$FINAL_STATUS" > v7a_r1.exit_status
printf '%s\n' "$CLASSIFICATION" > v7a_r1.termination_classification

{
    echo "Finished: $(date)"
    echo "ORCA shell status: $SHELL_STATUS"
    echo "Normal termination markers: $NORMAL_COUNT"
    echo "Error termination markers: $ERROR_COUNT"
    echo "Optimization convergence markers: $OPT_COUNT"
    echo "Final classified status: $FINAL_STATUS"
    echo "Classification: $CLASSIFICATION"
} >> v7a_r1.supervisor.log

exit "$FINAL_STATUS"
