from pathlib import Path
import argparse
import re
import subprocess
import time

ROOT = Path("runs/phase1A/day016_md_bath_extraction/orca_embedding_pilot_inputs")
PAT = re.compile(r"frame(\d{3})_(PYR[2-5])_embedding\.inp$")

def output_is_done(out: Path) -> bool:
    if not out.exists():
        return False
    txt = out.read_text(errors="ignore")
    return (
        "ORCA TERMINATED NORMALLY" in txt
        and "ORCA-CIS/TD-DFT FINISHED WITHOUT ERROR" in txt
    )

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frame-min", type=int, default=0)
    ap.add_argument("--frame-max", type=int, default=20)
    ap.add_argument("--max-jobs", type=int, default=None)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    jobs = []
    for inp in sorted(ROOT.glob("frame*_PYR*_embedding/frame*_PYR*_embedding.inp")):
        m = PAT.match(inp.name)
        if not m:
            continue
        frame = int(m.group(1))
        chrom = m.group(2)
        if frame < args.frame_min or frame > args.frame_max:
            continue
        out = inp.with_suffix(".out")
        jobs.append((frame, chrom, inp, out))

    pending = []
    done = []
    for frame, chrom, inp, out in jobs:
        if output_is_done(out):
            done.append((frame, chrom, inp, out))
        else:
            pending.append((frame, chrom, inp, out))

    if args.max_jobs is not None:
        pending = pending[:args.max_jobs]

    print("ORCA embedding batch runner")
    print("frame_min:", args.frame_min)
    print("frame_max:", args.frame_max)
    print("jobs_in_range:", len(jobs))
    print("already_done:", len(done))
    print("pending_selected:", len(pending))
    print()

    for frame, chrom, inp, out in pending:
        print(f"RUN frame={frame:03d} chrom={chrom} workdir={inp.parent}")
        if args.dry_run:
            continue

        t0 = time.time()
        with out.open("w") as fout:
            ret = subprocess.run(
                ["orca", inp.name],
                cwd=inp.parent,
                stdout=fout,
                stderr=subprocess.STDOUT,
                check=False,
            )
        dt_min = (time.time() - t0) / 60.0

        ok = output_is_done(out)
        print(
            f"DONE frame={frame:03d} chrom={chrom} "
            f"returncode={ret.returncode} ok={ok} runtime_min={dt_min:.2f}"
        )

        if not ok:
            print(f"WARNING: calculation did not pass completion check: {out}")

    print()
    print("Batch finished.")

if __name__ == "__main__":
    main()
