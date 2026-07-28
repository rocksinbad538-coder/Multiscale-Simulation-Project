#!/bin/bash

set -uo pipefail

EXEC_DIR='/Users/alejandro/projects/Multiscale-Simulation-Project/Microtubule_Inspired_Phase0/runs/phase1A/day035_qm_f06_upper_v7a_executions/v7a_20260727T193301Z'
ORCA_BIN='/Users/alejandro/projects/orca_6_1_1_macosx_intel_openmpi411/orca'

cd "$EXEC_DIR" || {
    echo "ERROR: cannot enter execution directory" > v7a.supervisor.log
    printf '%s\n' 90 > v7a.exit_status
    exit 90
}

printf '%s\n' "$$" > v7a.supervisor_pid

{
    echo "========================================================="
    echo "QM_F06 UPPER V7-A SUPERVISED EXECUTION"
    echo "========================================================="
    echo "Started: $(date)"
    echo "Execution directory: $EXEC_DIR"
    echo "ORCA binary: $ORCA_BIN"
} > v7a.supervisor.log

"$ORCA_BIN" v7a.inp > v7a.out 2> v7a.stderr &

ORCA_PID="$!"
printf '%s\n' "$ORCA_PID" > v7a.orca_pid

{
    echo "ORCA PID: $ORCA_PID"
    echo "Waiting for ORCA..."
} >> v7a.supervisor.log

wait "$ORCA_PID"
STATUS="$?"

printf '%s\n' "$STATUS" > v7a.exit_status

{
    echo "Finished: $(date)"
    echo "Exit status: $STATUS"
} >> v7a.supervisor.log

exit "$STATUS"
