from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import subprocess
import shutil

BASE = Path("runs/phase2/day053_dense_electronic_pilot/orca_inputs")
ORCA = shutil.which("orca")

if ORCA is None:
    raise SystemExit("ERROR: ORCA not found in PATH")

jobs = []

for d in sorted(BASE.glob("dense*_PYR*")):
    stem = d.name

    # dense000 == validated historical M4 frame020 at 100 ps.
    if stem.startswith("dense000_"):
        continue

    inp = d / f"{stem}.inp"
    out = d / f"{stem}.out"

    if not inp.exists():
        raise SystemExit(f"ERROR: missing input {inp}")

    if out.exists() and "ORCA TERMINATED NORMALLY" in out.read_text(errors="ignore"):
        print(f"SKIP_COMPLETED {stem}")
        continue

    jobs.append((stem, d, inp.name, out.name))


def run_job(job):
    stem, workdir, inp, out = job

    with open(workdir / out, "w") as fh:
        p = subprocess.run(
            [ORCA, inp],
            cwd=workdir,
            stdout=fh,
            stderr=subprocess.STDOUT,
        )

    text = (workdir / out).read_text(errors="ignore")
    normal = "ORCA TERMINATED NORMALLY" in text

    return stem, p.returncode, normal


print("=" * 72)
print("DAY053 DENSE ELECTRONIC PILOT")
print("ORCA =", ORCA)
print("MAX_CONCURRENT = 2")
print("NEW_JOBS =", len(jobs))
print("ANCHOR_SKIPPED = dense000_PYR3,dense000_PYR4")
print("=" * 72)

fails = []

with ThreadPoolExecutor(max_workers=2) as ex:
    futures = {ex.submit(run_job, j): j[0] for j in jobs}

    for fut in as_completed(futures):
        stem, rc, normal = fut.result()
        status = "PASS" if rc == 0 and normal else "FAIL"
        print(f"{stem} RC={rc} NORMAL={normal} STATUS={status}", flush=True)

        if status != "PASS":
            fails.append(stem)

print("=" * 72)
print("FAILED_JOBS =", len(fails))
print("PRODUCTION_GATE =", "PASS" if not fails else "FAIL")
if fails:
    print("FAIL_LIST =", ",".join(fails))
