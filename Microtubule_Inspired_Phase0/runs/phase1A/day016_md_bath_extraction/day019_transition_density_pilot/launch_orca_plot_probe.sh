#!/usr/bin/env bash
set -euo pipefail

SITE="${1:-PYR2}"
ROOT="runs/phase1A/day016_md_bath_extraction/day019_transition_density_pilot"
SITE_DIR="$ROOT/inputs/$SITE"

case "$SITE" in
  PYR2|PYR3|PYR4|PYR5) ;;
  *)
    echo "ERROR: site must be PYR2, PYR3, PYR4, or PYR5" >&2
    exit 2
    ;;
esac

if [[ ! -e "$SITE_DIR/pilot.gbw" || ! -e "$SITE_DIR/pilot.cis" ]]; then
  echo "ERROR: missing pilot.gbw or pilot.cis under $SITE_DIR" >&2
  exit 3
fi

cd "$SITE_DIR"

echo "Launching ORCA transition-density probe for $SITE"
echo "Transcript: orca_plot_probe_${SITE}.typescript"
echo "Select option 7 (CIS/TD-DFT transition densities)."
echo "This is a menu/protocol discovery run; do not generate production cubes yet."

script -q "orca_plot_probe_${SITE}.typescript" \
  orca_plot pilot.gbw -i
