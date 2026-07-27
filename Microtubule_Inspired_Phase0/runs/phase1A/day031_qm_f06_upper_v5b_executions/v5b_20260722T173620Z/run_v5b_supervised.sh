#!/bin/bash
set -u

cd "$(dirname "$0")" || exit 1

echo "$$" > v5b.supervisor_pid

"/Users/alejandro/projects/orca_6_1_1_macosx_intel_openmpi411/orca" v5b.inp > v5b.out 2> v5b.stderr &
orca_pid="$!"

echo "$orca_pid" > v5b.orca_pid

wait "$orca_pid"
status="$?"

echo "$status" > v5b.exit_status
exit "$status"
