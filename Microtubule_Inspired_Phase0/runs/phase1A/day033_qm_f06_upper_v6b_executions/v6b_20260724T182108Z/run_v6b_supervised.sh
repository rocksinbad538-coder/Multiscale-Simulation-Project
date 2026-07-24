#!/bin/bash

set -uo pipefail

EXEC_DIR='/Users/alejandro/projects/Multiscale-Simulation-Project/Microtubule_Inspired_Phase0/runs/phase1A/day033_qm_f06_upper_v6b_executions/v6b_20260724T182108Z'
ORCA_BIN='/Users/alejandro/projects/orca_6_1_1_macosx_intel_openmpi411/orca'

cd "$EXEC_DIR" || {
    echo "ERROR: cannot enter execution directory" > v6b.supervisor.log
    printf '%s\n' 90 > v6b.exit_status
    exit 90
}

printf '%s\n' "$$" > v6b.supervisor_pid

{
    echo "========================================================="
    echo "QM_F06 UPPER V6-B SUPERVISED EXECUTION"
    echo "========================================================="
    echo "Started: $(date)"
    echo "Execution directory: $EXEC_DIR"
    echo "ORCA binary: $ORCA_BIN"
} > v6b.supervisor.log

"$ORCA_BIN" v6b.inp > v6b.out 2> v6b.stderr &

ORCA_PID="$!"

printf '%s\n' "$ORCA_PID" > v6b.orca_pid

{
    echo "ORCA PID: $ORCA_PID"
    echo "Waiting for ORCA..."
} >> v6b.supervisor.log

wait "$ORCA_PID"
STATUS="$?"

printf '%s\n' "$STATUS" > v6b.exit_status

{
    echo "Finished: $(date)"
    echo "Exit status: $STATUS"
} >> v6b.supervisor.log

exit "$STATUS"
