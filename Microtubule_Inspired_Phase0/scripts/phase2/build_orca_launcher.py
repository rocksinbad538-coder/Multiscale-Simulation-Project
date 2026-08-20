#!/usr/bin/env python3

from pathlib import Path
import json

ROOT = Path("runs/phase2/campaign_phase5_corrected/orca_campaign")

manifest = json.loads(
    (ROOT/"ORCA_CAMPAIGN_MANIFEST.json").read_text()
)

script = ROOT/"run_all_orca.sh"

with script.open("w") as f:

    f.write("#!/bin/bash\n\n")
    f.write("set -e\n\n")

    for job in manifest:

        folder = Path(job["input"]).parent

        f.write(f'echo "=================================================="\n')
        f.write(f'echo "{folder.name}"\n')
        f.write(f'echo "=================================================="\n')
        f.write(f"cd {folder}\n")
        f.write("orca job.inp > job.out\n")
        f.write("cd - >/dev/null\n\n")

script.chmod(0o755)

print("="*90)
print("PHASE4-A07")
print("ORCA LAUNCHER CREATED")
print("="*90)
print(script)
