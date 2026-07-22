#!/bin/bash
set -u

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

LATEST_POINTER="$ROOT/runs/phase1A/day031_qm_f06_upper_v5b_executions/LATEST_V5B_EXECUTION.txt"
LATEST_REL="$(cat "$LATEST_POINTER")"
EXEC_DIR="$ROOT/$LATEST_REL"

OUT="$EXEC_DIR/v5b.out"
EXIT_FILE="$EXEC_DIR/v5b.exit_status"
WATCH_LOG="$EXEC_DIR/v5b.overnight_watch.log"
SUMMARY="$EXEC_DIR/v5b.overnight_completion_summary.txt"

{
    echo "========================================================="
    echo "QM_F06 UPPER V5-B OVERNIGHT WATCH"
    echo "Started: $(date)"
    echo "Execution directory: $EXEC_DIR"
    echo "========================================================="
} >> "$WATCH_LOG"

while true
do
    if [ -f "$EXIT_FILE" ]; then
        status="$(cat "$EXIT_FILE")"

        {
            echo
            echo "ORCA supervisor finished: $(date)"
            echo "Exit status: $status"
        } >> "$WATCH_LOG"

        {
            echo "========================================================="
            echo "QM_F06 UPPER V5-B COMPLETION SUMMARY"
            echo "========================================================="
            echo "Generated: $(date)"
            echo "Execution directory: $EXEC_DIR"
            echo "Exit status: $status"
            echo
            echo "Geometry cycles:"
            grep -c "GEOMETRY OPTIMIZATION CYCLE" "$OUT" || true
            echo
            echo "SCF convergence records:"
            grep -c "SCF CONVERGED AFTER" "$OUT" || true
            echo
            grep -E \
            "THE OPTIMIZATION HAS CONVERGED|FINAL SINGLE POINT ENERGY|ORCA TERMINATED NORMALLY|ORCA finished by error termination|SCF NOT CONVERGED" \
            "$OUT" \
            | tail -40 || true
        } > "$SUMMARY"

        break
    fi

    if ! pgrep -f "orca v5b.inp" >/dev/null 2>&1; then
        {
            echo
            echo "WARNING: ORCA process absent but exit-status file not present."
            echo "Detected: $(date)"
        } >> "$WATCH_LOG"
    fi

    {
        printf "%s | size=" "$(date '+%Y-%m-%d %H:%M:%S')"
        stat -f "%z" "$OUT" 2>/dev/null || echo "missing"
    } >> "$WATCH_LOG"

    sleep 300
done
