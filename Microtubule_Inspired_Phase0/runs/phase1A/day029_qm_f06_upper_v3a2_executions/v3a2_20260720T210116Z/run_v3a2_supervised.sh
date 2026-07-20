#!/usr/bin/env bash

set -uo pipefail

cd "$(dirname "$0")" || exit 1

ORCA_BIN="$(command -v orca)"

if [ -z "${ORCA_BIN:-}" ]; then
    echo "ORCA executable not found." >&2
    echo "127" > v3a2.exit_status
    exit 127
fi

if pgrep -af "[o]rca.*v3a2.inp" >/dev/null 2>&1; then
    echo "Another ORCA process using v3a2.inp appears active." >&2
    echo "98" > v3a2.exit_status
    exit 98
fi

{
    echo "ORCA executable: $ORCA_BIN"
    echo "Execution directory: $PWD"
    echo "Supervisor PID: $$"
    echo "Start UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "Start local: $(date '+%Y-%m-%d %H:%M:%S %Z')"
} > v3a2.execution_metadata

"$ORCA_BIN" v3a2.inp \
    > v3a2.out \
    2> v3a2.stderr &

orca_pid=$!

echo "$orca_pid" > v3a2.orca_pid

wait "$orca_pid"
status=$?

echo "$status" > v3a2.exit_status

{
    echo "End UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "End local: $(date '+%Y-%m-%d %H:%M:%S %Z')"
    echo "Exit status: $status"
} >> v3a2.execution_metadata

exit "$status"
