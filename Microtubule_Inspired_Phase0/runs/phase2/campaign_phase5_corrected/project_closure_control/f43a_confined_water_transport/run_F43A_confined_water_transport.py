#!/usr/bin/env python3

from pathlib import Path
import csv
import json
import math
import sys
import numpy as np

print("=" * 148)
print("PHASE5-F43A — CONFINED-WATER RESIDENCE AND TRANSLATIONAL DIFFUSION")
print("DELIVERABLE_TARGET=M2.4_RESIDENCE + M2.5_DIFFUSION")
print("SOURCE=EXISTING_ACCEPTED_MOBILE_MD_TRAJECTORY")
print("=" * 148)

ROOT = Path.cwd()
P0 = ROOT / "Microtubule_Inspired_Phase0"

STAGE = (
    P0 /
    "runs/phase1A/day021_mobile_restraint_protocol/"
    "execution/08_nvt_mobile_100ps"
)

DAY020 = (
    P0 /
    "runs/phase1A/day020_confined_water_axial_radial_density"
)

NDX = DAY020 / "confined_water_analysis_groups.ndx"

BOUND = (
    DAY020 /
    "profile_guided_classification/"
    "profile_guided_boundaries.csv"
)

OUT = (
    P0 /
    "runs/phase2/campaign_phase5_corrected/"
    "project_closure_control/"
    "f43a_confined_water_transport"
)

OUT.mkdir(parents=True, exist_ok=True)

# ================================================================
# 0. Dependency check
# ================================================================

try:
    import MDAnalysis as mda
except Exception as e:
    print(f"F43A_ABORT=MDANALYSIS_IMPORT_FAILED {e}")
    sys.exit(2)

print("\n[1] INPUT DISCOVERY")
print("-" * 120)

tprs = sorted(STAGE.glob("*08_nvt_mobile_100ps*.tpr"))
xtcs = sorted(STAGE.glob("*08_nvt_mobile_100ps*.xtc"))

if len(tprs) != 1:
    print("TPR_CANDIDATES=")
    for x in tprs:
        print(x)
    raise SystemExit(
        f"F43A_ABORT=EXPECTED_ONE_TPR_FOUND_{len(tprs)}"
    )

if len(xtcs) != 1:
    print("XTC_CANDIDATES=")
    for x in xtcs:
        print(x)
    raise SystemExit(
        f"F43A_ABORT=EXPECTED_ONE_XTC_FOUND_{len(xtcs)}"
    )

TPR = tprs[0]
XTC = xtcs[0]

for label, p in [
    ("TPR", TPR),
    ("XTC", XTC),
    ("NDX", NDX),
    ("PROFILE_BOUNDARIES", BOUND),
]:
    print(
        f"{label}="
        + ("FOUND" if p.exists() else "NOT_FOUND")
        + f" path={p.relative_to(ROOT)}"
    )

if not all(p.exists() for p in [TPR, XTC, NDX, BOUND]):
    raise SystemExit("F43A_ABORT=MISSING_CANONICAL_INPUT")

# ================================================================
# 1. NDX parsing
# ================================================================

def read_ndx(path):
    groups = {}
    current = None

    for raw in path.read_text(errors="ignore").splitlines():
        line = raw.strip()

        if not line:
            continue

        if line.startswith("[") and line.endswith("]"):
            current = line[1:-1].strip()
            groups[current] = []
            continue

        if current is not None:
            for token in line.split():
                try:
                    groups[current].append(int(token) - 1)
                except ValueError:
                    pass

    return groups


groups = read_ndx(NDX)

print("\n[2] INDEX-GROUP AUDIT")
print("-" * 120)

for name, ids in groups.items():
    print(f"GROUP={name} N={len(ids)}")

# Resolve scaffold group conservatively.
hbn_candidates = [
    name for name in groups
    if "hbn" in name.lower()
    or "scaffold" in name.lower()
]

if not hbn_candidates:
    raise SystemExit("F43A_ABORT=NO_HBN_SCAFFOLD_GROUP")

# Prefer exact-looking HBN group, then largest candidate.
hbn_candidates.sort(
    key=lambda n: (
        0 if n.lower() in {"hbn", "scaffold", "hbn_scaffold"} else 1,
        -len(groups[n])
    )
)

HBN_GROUP = hbn_candidates[0]
HBN_IDX = np.asarray(groups[HBN_GROUP], dtype=int)

# Resolve water oxygen group if available.
ow_candidates = [
    name for name in groups
    if (
        "water" in name.lower()
        and (
            "oxygen" in name.lower()
            or "_o" in name.lower()
        )
    )
]

print(f"SELECTED_HBN_GROUP={HBN_GROUP}")
print(f"HBN_ATOMS={len(HBN_IDX)}")

# ================================================================
# 2. Load trajectory
# ================================================================

print("\n[3] TRAJECTORY AUDIT")
print("-" * 120)

u = mda.Universe(str(TPR), str(XTC))

print(f"N_ATOMS={u.atoms.n_atoms}")
print(f"N_FRAMES={len(u.trajectory)}")

if HBN_IDX.max() >= u.atoms.n_atoms:
    raise SystemExit(
        "F43A_ABORT=NDX_TOPOLOGY_ATOM_ORDER_MISMATCH"
    )

# Water oxygen resolution:
if ow_candidates:
    OW_GROUP = ow_candidates[0]
    OW_IDX = np.asarray(groups[OW_GROUP], dtype=int)

    if OW_IDX.max() >= u.atoms.n_atoms:
        raise SystemExit(
            "F43A_ABORT=WATER_OXYGEN_NDX_INDEX_MISMATCH"
        )

    print(f"SELECTED_WATER_O_GROUP={OW_GROUP}")
else:
    # Conservative topology-derived fallback.
    possible = u.select_atoms(
        "resname SOL HOH WAT TIP3 TIP3P SPC SPCE "
        "and name O OW OWT OH2"
    )

    if len(possible) == 0:
        # second fallback: oxygen-like names + water-like residues
        water_resnames = {
            "SOL", "HOH", "WAT", "TIP3",
            "TIP3P", "SPC", "SPCE"
        }

        idx = []
        for atom in u.atoms:
            if (
                atom.resname.upper() in water_resnames
                and atom.name.upper() in {
                    "O", "OW", "OWT", "OH2"
                }
            ):
                idx.append(atom.index)

        OW_IDX = np.asarray(idx, dtype=int)
    else:
        OW_IDX = possible.indices.copy()

    print("SELECTED_WATER_O_GROUP=TOPOLOGY_FALLBACK")

if len(OW_IDX) == 0:
    raise SystemExit("F43A_ABORT=NO_WATER_OXYGENS")

print(f"WATER_OXYGENS={len(OW_IDX)}")

# ================================================================
# 3. Resolve radial confinement boundary from Day020 product
# ================================================================

print("\n[4] CONFINED-REGION CONTRACT")
print("-" * 120)

with BOUND.open(newline="", errors="ignore") as fh:
    rows = list(csv.DictReader(fh))

print(f"BOUNDARY_ROWS={len(rows)}")

for row in rows:
    print("BOUNDARY_ROW=" + repr(row))

def numeric_values(row):
    out = []
    for k, v in row.items():
        if v is None:
            continue
        try:
            x = float(v)
            if math.isfinite(x):
                out.append((k, x))
        except Exception:
            pass
    return out


radial_candidates = []

for row in rows:
    text = " ".join(
        str(v).lower()
        for v in row.values()
        if v is not None
    )

    # Require radial + physically interpretable profile feature.
    if (
        "radial" in text
        and (
            "depletion" in text
            or "minimum" in text
            or "wall" in text
            or "lumen" in text
        )
    ):
        for key, val in numeric_values(row):
            kl = key.lower()

            if (
                (
                    "nm" in kl
                    or "coordinate" in kl
                    or "position" in kl
                    or "boundary" in kl
                    or "value" in kl
                )
                and 0.5 < abs(val) < 5.0
            ):
                radial_candidates.append(
                    (text, key, abs(val))
                )

# unique approximately
unique = []
for rec in radial_candidates:
    val = rec[2]
    if not any(abs(val - x[2]) < 1e-6 for x in unique):
        unique.append(rec)

print("RADIAL_BOUNDARY_CANDIDATES=")

for x in unique:
    print(
        f"  VALUE_NM={x[2]:.9f} "
        f"FIELD={x[1]} TEXT={x[0][:160]}"
    )

# Preferred validated Day020 profile minimum ≈1.174620 nm.
preferred = [
    x for x in unique
    if 1.0 <= x[2] <= 1.35
]

if len(preferred) == 1:
    R_CUT_NM = preferred[0][2]

elif len(preferred) > 1:
    # choose candidate closest to previously validated depletion minimum
    preferred.sort(
        key=lambda x: abs(x[2] - 1.174620)
    )

    if (
        len(preferred) > 1
        and abs(
            abs(preferred[0][2] - 1.174620)
            - abs(preferred[1][2] - 1.174620)
        ) < 1e-5
    ):
        raise SystemExit(
            "F43A_ABORT=AMBIGUOUS_RADIAL_BOUNDARY"
        )

    R_CUT_NM = preferred[0][2]

else:
    # Explicit documented fallback from validated Day020 closeout.
    R_CUT_NM = 1.174620
    print(
        "RADIAL_BOUNDARY_SOURCE="
        "DOCUMENTED_DAY020_DEPLETION_MINIMUM_FALLBACK"
    )

print(f"R_CUT_NM={R_CUT_NM:.9f}")

# ================================================================
# 4. Reference scaffold geometry
# ================================================================

u.trajectory[0]

ref_hbn_A = u.atoms[HBN_IDX].positions.astype(float)
ref_center_A = ref_hbn_A.mean(axis=0)

X = ref_hbn_A - ref_center_A
cov = X.T @ X / len(X)

evals, evecs = np.linalg.eigh(cov)
axis = evecs[:, np.argmax(evals)]
axis = axis / np.linalg.norm(axis)

# Stable sign convention.
major = np.argmax(np.abs(axis))
if axis[major] < 0:
    axis *= -1

z_ref_nm = (
    (ref_hbn_A - ref_center_A) @ axis
) / 10.0

Z_MIN_NM = float(z_ref_nm.min())
Z_MAX_NM = float(z_ref_nm.max())

print(
    "AXIS_REFERENCE="
    f"({axis[0]:+.8f},{axis[1]:+.8f},{axis[2]:+.8f})"
)

print(
    "REFERENCE_CENTER_A="
    f"({ref_center_A[0]:.6f},"
    f"{ref_center_A[1]:.6f},"
    f"{ref_center_A[2]:.6f})"
)

print(f"Z_MIN_NM={Z_MIN_NM:.9f}")
print(f"Z_MAX_NM={Z_MAX_NM:.9f}")
print(f"Z_SPAN_NM={Z_MAX_NM-Z_MIN_NM:.9f}")

print(
    "CONFINED_REGION="
    "CORE_LUMEN: "
    "r<=profile_guided_radial_wall_boundary AND "
    "HBN_reference_zmin<=z<=HBN_reference_zmax"
)

# ================================================================
# 5. Kabsch rigid alignment
# ================================================================

def kabsch(current, reference):
    c0 = current.mean(axis=0)
    r0 = reference.mean(axis=0)

    C = current - c0
    R = reference - r0

    H = C.T @ R

    U, s, Vt = np.linalg.svd(H)

    Q = U @ Vt

    if np.linalg.det(Q) < 0:
        Vt[-1, :] *= -1
        Q = U @ Vt

    # row-vector convention
    return c0, r0, Q


def align_points(points, c0, r0, Q):
    return (points - c0) @ Q + r0


# ================================================================
# 6. First pass — occupancy matrix
# ================================================================

print("\n[5] FIRST PASS — CONFINED-WATER OCCUPANCY")
print("-" * 120)

nF = len(u.trajectory)
nW = len(OW_IDX)

occupancy = np.zeros(
    (nF, nW),
    dtype=np.bool_
)

times_ps = np.zeros(nF, dtype=float)

for iframe, ts in enumerate(u.trajectory):

    cur_hbn = u.atoms[HBN_IDX].positions.astype(float)

    c0, r0, Q = kabsch(cur_hbn, ref_hbn_A)

    water = u.atoms[OW_IDX].positions.astype(float)
    water_aligned = align_points(
        water,
        c0,
        r0,
        Q
    )

    d_nm = (
        water_aligned - ref_center_A
    ) / 10.0

    zz = d_nm @ axis
    perp = d_nm - np.outer(zz, axis)
    rr = np.linalg.norm(perp, axis=1)

    inside = (
        (rr <= R_CUT_NM)
        & (zz >= Z_MIN_NM)
        & (zz <= Z_MAX_NM)
    )

    occupancy[iframe] = inside

    try:
        times_ps[iframe] = float(ts.time)
    except Exception:
        times_ps[iframe] = float(iframe)

    if (
        iframe == 0
        or iframe == nF - 1
        or iframe % max(1, nF // 10) == 0
    ):
        print(
            f"FRAME={iframe}/{nF-1} "
            f"TIME_PS={times_ps[iframe]:.6f} "
            f"N_CONFINED={inside.sum()}"
        )

if nF < 3:
    raise SystemExit("F43A_ABORT=TRAJECTORY_TOO_SHORT")

dt = float(np.median(np.diff(times_ps)))

print(f"FRAME_DT_PS={dt:.9f}")
print(
    f"MEAN_CONFINED_WATERS="
    f"{occupancy.sum(axis=1).mean():.6f}"
)
print(
    f"MIN_CONFINED_WATERS="
    f"{occupancy.sum(axis=1).min()}"
)
print(
    f"MAX_CONFINED_WATERS="
    f"{occupancy.sum(axis=1).max()}"
)

ever = occupancy.any(axis=0)
ever_idx_local = np.where(ever)[0]

print(
    f"WATERS_EVER_IN_CORE_LUMEN="
    f"{len(ever_idx_local)}/{nW}"
)

# ================================================================
# 7. Residence correlation
# ================================================================

print("\n[6] RESIDENCE CORRELATION")
print("-" * 120)

max_lag = nF // 2

cres = np.full(max_lag + 1, np.nan)
denoms = np.zeros(max_lag + 1, dtype=np.int64)

cres[0] = 1.0
denoms[0] = int(occupancy.sum())

for lag in range(1, max_lag + 1):

    a = occupancy[:-lag]
    b = occupancy[lag:]

    denom = int(a.sum())
    denoms[lag] = denom

    if denom > 0:
        cres[lag] = (
            np.logical_and(a, b).sum()
            / denom
        )

lags_ps = np.arange(max_lag + 1) * dt

# Finite-window integral.
valid = np.isfinite(cres)

tau_res_window_ps = float(
    np.trapz(
        cres[valid],
        lags_ps[valid]
    )
)

# First 1/e crossing
cross = np.where(cres <= 1.0 / math.e)[0]

if len(cross):
    tau_1e_ps = float(lags_ps[cross[0]])
else:
    tau_1e_ps = None

print(
    f"RESIDENCE_INTEGRAL_WINDOW_PS="
    f"{tau_res_window_ps:.9f}"
)

print(
    "RESIDENCE_1E_CROSSING_PS="
    + (
        f"{tau_1e_ps:.9f}"
        if tau_1e_ps is not None
        else "NOT_REACHED"
    )
)

print(
    f"C_RES_FINAL="
    f"{cres[-1]:.9f}"
)

# Continuous run-length residence statistics
run_lengths = []

for j in ever_idx_local:

    arr = occupancy[:, j]

    starts = np.where(
        arr & np.r_[True, ~arr[:-1]]
    )[0]

    ends = np.where(
        arr & np.r_[~arr[1:], True]
    )[0]

    for s, e in zip(starts, ends):
        run_lengths.append(e - s + 1)

run_lengths = np.asarray(
    run_lengths,
    dtype=int
)

run_ps = run_lengths * dt

print(f"CONTINUOUS_RESIDENCE_EVENTS={len(run_ps)}")

if len(run_ps):

    print(
        f"CONTINUOUS_RESIDENCE_MEAN_PS="
        f"{run_ps.mean():.9f}"
    )

    print(
        f"CONTINUOUS_RESIDENCE_MEDIAN_PS="
        f"{np.median(run_ps):.9f}"
    )

    print(
        f"CONTINUOUS_RESIDENCE_P95_PS="
        f"{np.percentile(run_ps,95):.9f}"
    )

    print(
        f"CONTINUOUS_RESIDENCE_MAX_PS="
        f"{run_ps.max():.9f}"
    )

# ================================================================
# 8. Second pass — aligned coordinates only for waters that enter
# ================================================================

print("\n[7] SECOND PASS — CONFINED-WATER COORDINATES")
print("-" * 120)

OW_EVER_IDX = OW_IDX[ever_idx_local]

coords = np.empty(
    (nF, len(OW_EVER_IDX), 3),
    dtype=np.float32
)

occ_e = occupancy[:, ever_idx_local]

for iframe, ts in enumerate(u.trajectory):

    cur_hbn = u.atoms[HBN_IDX].positions.astype(float)
    c0, r0, Q = kabsch(cur_hbn, ref_hbn_A)

    water = (
        u.atoms[OW_EVER_IDX]
        .positions.astype(float)
    )

    water_aligned = align_points(
        water,
        c0,
        r0,
        Q
    )

    coords[iframe] = (
        (water_aligned - ref_center_A)
        / 10.0
    ).astype(np.float32)

    if (
        iframe == 0
        or iframe == nF - 1
        or iframe % max(1, nF // 10) == 0
    ):
        print(
            f"FRAME={iframe}/{nF-1}"
        )

# ================================================================
# 9. Continuous-survival conditioned MSD
# ================================================================

print("\n[8] CONTINUOUS-CONFINEMENT MSD")
print("-" * 120)

# cumulative count of outside states per water
outside = (~occ_e).astype(np.int32)

cs = np.vstack([
    np.zeros(
        (1, outside.shape[1]),
        dtype=np.int32
    ),
    np.cumsum(
        outside,
        axis=0,
        dtype=np.int32
    )
])

# Analyze a dense lag grid, but not more than half trajectory.
lag_indices = np.arange(
    1,
    max_lag + 1,
    dtype=int
)

msd_rows = []

for lag in lag_indices:

    # Number of outside frames in inclusive interval [t,t+lag]
    exits = (
        cs[lag+1:]
        - cs[:-(lag+1)]
    )

    continuous = (exits == 0)

    if not np.any(continuous):
        msd_rows.append(
            (
                lag,
                lag * dt,
                0,
                np.nan,
                np.nan,
                np.nan,
                np.nan
            )
        )
        continue

    disp = (
        coords[lag:]
        - coords[:-lag]
    )

    # Match shapes:
    # continuous is (nF-lag, nW)
    dr2 = np.sum(
        disp * disp,
        axis=2
    )

    z0 = coords[:-lag] @ axis
    z1 = coords[lag:] @ axis
    dz2 = (z1 - z0) ** 2

    radial_vec0 = (
        coords[:-lag]
        - z0[..., None] * axis
    )

    radial_vec1 = (
        coords[lag:]
        - z1[..., None] * axis
    )

    dperp = radial_vec1 - radial_vec0
    dperp2 = np.sum(
        dperp * dperp,
        axis=2
    )

    vals = dr2[continuous]
    vz = dz2[continuous]
    vp = dperp2[continuous]

    msd_rows.append(
        (
            lag,
            lag * dt,
            int(len(vals)),
            float(vals.mean()),
            float(vz.mean()),
            float(vp.mean()),
            float(vals.std(ddof=1))
            if len(vals) > 1
            else 0.0
        )
    )

# ================================================================
# 10. Diffusion fit diagnostics
# ================================================================

print("\n[9] DIFFUSION FIT DIAGNOSTICS")
print("-" * 120)

arr = np.asarray(
    [
        row
        for row in msd_rows
        if (
            row[2] >= 100
            and np.isfinite(row[3])
        )
    ],
    dtype=float
)

fit_result = {
    "accepted": False,
    "reason": None,
}

if len(arr) < 10:

    fit_result["reason"] = (
        "INSUFFICIENT_CONTINUOUS_SURVIVOR_MSD_POINTS"
    )

else:
    t = arr[:,1]
    msd3 = arr[:,3]
    msdz = arr[:,4]
    msdp = arr[:,5]

    # Candidate windows: latter 25–50%, 33–66%, and 50–80%
    fractions = [
        (0.25, 0.50),
        (0.33, 0.66),
        (0.50, 0.80),
    ]

    fits = []

    for f0, f1 in fractions:

        lo = t.min() + f0 * (t.max()-t.min())
        hi = t.min() + f1 * (t.max()-t.min())

        mask = (
            (t >= lo)
            & (t <= hi)
        )

        if mask.sum() < 5:
            continue

        def fit(y):
            x = t[mask]
            yy = y[mask]

            p = np.polyfit(x, yy, 1)
            pred = np.polyval(p, x)

            ss_res = np.sum(
                (yy - pred)**2
            )

            ss_tot = np.sum(
                (yy - yy.mean())**2
            )

            r2 = (
                1.0 - ss_res/ss_tot
                if ss_tot > 0
                else np.nan
            )

            return float(p[0]), float(p[1]), float(r2)

        s3, b3, r3 = fit(msd3)
        sz, bz, rz = fit(msdz)
        sp, bp, rp = fit(msdp)

        # nm^2/ps -> 1e-6 m^2/s
        #
        # 1 nm^2/ps = 1e-6 m^2/s
        D3 = s3 / 6.0
        Dz = sz / 2.0
        Dperp = sp / 4.0

        fits.append({
            "fraction_start": f0,
            "fraction_end": f1,
            "t_start_ps": float(t[mask][0]),
            "t_end_ps": float(t[mask][-1]),
            "n_points": int(mask.sum()),

            "slope_msd3_nm2_ps": s3,
            "R2_msd3": r3,
            "D3_nm2_ps": D3,
            "D3_m2_s": D3 * 1e-6,

            "slope_msd_z_nm2_ps": sz,
            "R2_msd_z": rz,
            "Dz_nm2_ps": Dz,
            "Dz_m2_s": Dz * 1e-6,

            "slope_msd_perp_nm2_ps": sp,
            "R2_msd_perp": rp,
            "Dperp_nm2_ps": Dperp,
            "Dperp_m2_s": Dperp * 1e-6,
        })

    print(f"CANDIDATE_FIT_WINDOWS={len(fits)}")

    for x in fits:
        print(
            "FIT "
            f"{x['t_start_ps']:.4f}-"
            f"{x['t_end_ps']:.4f} ps "
            f"D3={x['D3_m2_s']:.6e} m2/s "
            f"R2={x['R2_msd3']:.5f} "
            f"Dz={x['Dz_m2_s']:.6e} "
            f"Dperp={x['Dperp_m2_s']:.6e}"
        )

    acceptable = [
        x for x in fits
        if (
            x["D3_nm2_ps"] > 0
            and x["R2_msd3"] >= 0.95
            and x["R2_msd_z"] >= 0.90
            and x["R2_msd_perp"] >= 0.90
        )
    ]

    if acceptable:
        # choose best 3D linearity
        acceptable.sort(
            key=lambda x: x["R2_msd3"],
            reverse=True
        )

        fit_result = {
            "accepted": True,
            "reason":
                "LINEAR_MSD_WINDOW_IDENTIFIED",
            "selected":
                acceptable[0],
            "all_windows":
                fits,
        }

    else:
        fit_result = {
            "accepted": False,
            "reason":
                "NO_ROBUST_LINEAR_DIFFUSIVE_WINDOW",
            "all_windows":
                fits,
        }

print(
    "DIFFUSION_FIT_ACCEPTED="
    + str(fit_result["accepted"])
)

print(
    f"DIFFUSION_FIT_REASON="
    f"{fit_result['reason']}"
)

if fit_result["accepted"]:
    s = fit_result["selected"]

    print(
        f"D_3D_M2_S="
        f"{s['D3_m2_s']:.9e}"
    )

    print(
        f"D_Z_M2_S="
        f"{s['Dz_m2_s']:.9e}"
    )

    print(
        f"D_PERP_M2_S="
        f"{s['Dperp_m2_s']:.9e}"
    )

# ================================================================
# 11. Write outputs
# ================================================================

print("\n[10] WRITE DELIVERABLE PRODUCTS")
print("-" * 120)

cres_csv = OUT / "M2_4_CONFINED_WATER_RESIDENCE_CORRELATION.csv"

with cres_csv.open("w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow([
        "lag_frames",
        "lag_ps",
        "C_res_intermittent",
        "occupancy_denominator"
    ])

    for i in range(len(cres)):
        w.writerow([
            i,
            lags_ps[i],
            cres[i],
            denoms[i],
        ])


runs_csv = OUT / "M2_4_CONTINUOUS_RESIDENCE_EVENTS.csv"

with runs_csv.open("w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow([
        "event_id",
        "duration_frames",
        "duration_ps"
    ])

    for i, n in enumerate(run_lengths):
        w.writerow([
            i + 1,
            int(n),
            float(n * dt)
        ])


msd_csv = OUT / "M2_5_CONFINED_WATER_MSD.csv"

with msd_csv.open("w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow([
        "lag_frames",
        "lag_ps",
        "continuous_survivor_samples",
        "MSD_3D_nm2",
        "MSD_axial_nm2",
        "MSD_transverse_nm2",
        "MSD_3D_sd_nm2",
    ])
    w.writerows(msd_rows)


summary = {
    "phase": "PHASE5-F43A",

    "source_tpr": str(
        TPR.relative_to(ROOT)
    ),

    "source_xtc": str(
        XTC.relative_to(ROOT)
    ),

    "trajectory_frames": nF,
    "frame_dt_ps": dt,

    "region_definition": {
        "type":
            "profile-guided core lumen",
        "radial_cutoff_nm":
            R_CUT_NM,
        "axial_min_nm":
            Z_MIN_NM,
        "axial_max_nm":
            Z_MAX_NM,
        "reference_frame":
            "HBN rigid-body aligned",
    },

    "water_population": {
        "total_water_oxygens":
            int(nW),
        "mean_confined":
            float(
                occupancy.sum(axis=1).mean()
            ),
        "minimum_confined":
            int(
                occupancy.sum(axis=1).min()
            ),
        "maximum_confined":
            int(
                occupancy.sum(axis=1).max()
            ),
        "waters_ever_confined":
            int(len(ever_idx_local)),
    },

    "residence": {
        "definition":
            "<h(0)h(t)>/<h(0)^2>",
        "finite_window_integral_ps":
            tau_res_window_ps,
        "one_over_e_crossing_ps":
            tau_1e_ps,
        "C_res_at_max_lag":
            float(cres[-1]),
        "continuous_event_count":
            int(len(run_ps)),
        "continuous_mean_ps":
            float(run_ps.mean())
            if len(run_ps) else None,
        "continuous_median_ps":
            float(np.median(run_ps))
            if len(run_ps) else None,
        "continuous_p95_ps":
            float(np.percentile(run_ps,95))
            if len(run_ps) else None,
        "continuous_max_ps":
            float(run_ps.max())
            if len(run_ps) else None,
    },

    "diffusion": fit_result,

    "scientific_guardrails": [
        "Residence is defined by occupancy of the geometric core lumen, not chromophore contact.",
        "MSD is conditioned on continuous confinement over each lag interval.",
        "The coordinate frame is rigid-body aligned to the HBN scaffold.",
        "A diffusion coefficient is accepted only if a robust linear MSD regime is detected.",
        "If no linear regime is detected, the correct deliverable is an upper/lower mobility characterization rather than an invented D.",
    ],

    "new_md_simulations": 0,
    "source_trajectory_modified": False,
}


summary_json = OUT / "PHASE5_F43A_CONFINED_WATER_TRANSPORT.json"

summary_json.write_text(
    json.dumps(
        summary,
        indent=2
    )
)

for p in [
    cres_csv,
    runs_csv,
    msd_csv,
    summary_json,
]:
    print(f"WROTE={p.relative_to(ROOT)}")

# ================================================================
# 12. Final scientific gate
# ================================================================

print("\n[11] FINAL SCIENTIFIC GATE")
print("-" * 120)

print("M2_4_RESIDENCE_ANALYSIS=COMPUTED")

if (
    tau_1e_ps is not None
    and cres[-1] < 0.5
):
    print(
        "M2_4_RESIDENCE_STATUS="
        "QUANTITATIVELY_RESOLVED_WITHIN_WINDOW"
    )
else:
    print(
        "M2_4_RESIDENCE_STATUS="
        "FINITE_WINDOW_LOWER_BOUND_OR_INCOMPLETE_DECAY"
    )

if fit_result["accepted"]:
    print(
        "M2_5_DIFFUSION_STATUS="
        "QUANTITATIVELY_RESOLVED"
    )
else:
    print(
        "M2_5_DIFFUSION_STATUS="
        "NO_CONVERGED_DIFFUSIVE_REGIME_IN_CURRENT_WINDOW"
    )

print(
    "M2_6_HYDROGEN_BONDS="
    "DEFERRED_TO_F43B_CHEMICAL_DEFINITION_GATE"
)

print("F43A_STATUS=0")
print("NEW_MD_SIMULATIONS=0")
print("SOURCE_FILES_MODIFIED=NO")
print(
    "NEXT_ACTION="
    "INTERPRET_RESIDENCE_AND_DIFFUSION_THEN_RUN_F43B_HBOND_GATE"
)

print("\n" + "=" * 148)
print("PHASE5-F43A COMPLETE")
print("TERMINAL_REMAINS_OPEN=YES")
print("=" * 148)

