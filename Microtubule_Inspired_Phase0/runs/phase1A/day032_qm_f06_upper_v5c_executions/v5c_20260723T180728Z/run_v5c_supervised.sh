#!/bin/bash

set -uo pipefail

EXEC_DIR='/Users/alejandro/projects/Multiscale-Simulation-Project/Microtubule_Inspired_Phase0/runs/phase1A/day032_qm_f06_upper_v5c_executions/v5c_20260723T180728Z'
ORCA_BIN='/Users/alejandro/projects/orca_6_1_1_macosx_intel_openmpi411/orca'

cd "$EXEC_DIR" || {
    echo "ERROR: cannot enter execution directory" > v5c.supervisor.log
    printf '%s\n' 90 > v5c.exit_status
    exit 90
}

printf '%s\n' "$$" > v5c.supervisor_pid

{
    echo "========================================================="
    echo "QM_F06 UPPER V5-C SUPERVISED EXECUTION"
    echo "========================================================="
    echo "Started: $(date)"
    echo "Execution directory: $EXEC_DIR"
    echo "ORCA binary: $ORCA_BIN"
} > v5c.supervisor.log

"$ORCA_BIN" v5c.inp > v5c.out 2> v5c.stderr &

ORCA_PID="$!"

printf '%s\n' "$ORCA_PID" > v5c.orca_pid

{
    echo "ORCA PID: $ORCA_PID"
    echo "Waiting for ORCA..."
} >> v5c.supervisor.log

wait "$ORCA_PID"
STATUS="$?"

printf '%s\n' "$STATUS" > v5c.exit_status

{
    echo "Finished: $(date)"
    echo "Exit status: $STATUS"
} >> v5c.supervisor.log

exit "$STATUS"
