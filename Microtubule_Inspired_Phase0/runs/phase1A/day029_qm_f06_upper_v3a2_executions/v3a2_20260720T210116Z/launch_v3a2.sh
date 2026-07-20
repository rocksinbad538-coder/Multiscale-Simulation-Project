#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

if pgrep -af "[o]rca.*v3a2.inp" >/dev/null 2>&1; then
    echo "Another v3a2 ORCA process appears to be active."
    exit 1
fi

ORCA_BIN="$(command -v orca)"

if [ -z "${ORCA_BIN:-}" ]; then
    echo "ORCA executable not found."
    exit 1
fi

echo "ORCA executable: $ORCA_BIN"
echo "Execution directory: $PWD"
echo "Start UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"

"$ORCA_BIN" v3a2.inp > v3a2.out 2> v3a2.stderr

status=$?

echo "$status" > v3a2.exit_status

echo "End UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Exit status: $status"

exit "$status"
