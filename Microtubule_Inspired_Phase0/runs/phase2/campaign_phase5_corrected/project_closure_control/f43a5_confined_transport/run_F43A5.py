#!/usr/bin/env python3

from pathlib import Path
import csv
import json
import math
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import minimize_vectors

print("=" * 152)
print("PHASE5-F43A5 — VALIDATED CONFINED-WATER RESIDENCE / MSD / DIFFUSION")
print("DELIVERABLE_TARGET=M2.4_RESIDENCE + M2.5_DIFFUSION")
print("COORDINATE_CONTRACT=VALIDATED_F43A4_DENSMAP_AXIS_PBC")
print("=" * 152)

ROOT = Path.cwd()
P0 = ROOT / "Microtubule_Inspired_Phase0"

MOBILE = (
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
    "profile_guided_classification/profile_guided_boundaries.csv"
)

OUT = (
    P0 /
    "runs/phase2/campaign_phase5_corrected/"
    "project_closure_control/f43a5_confined_transport"
)

OUT.mkdir(parents=True, exist_ok=True)

TPR = MOBILE / "08_nvt_mobile_100ps.tpr"
XTC = MOBILE / "08_nvt_mobile_100ps.xtc"


# =============================================================================
# Helpers
# =============================================================================

def read_ndx(path):

    groups = {}
    current = None

    for raw in path.read_text(errors="ignore").splitlines():

        s = raw.strip()

        if not s:
            continue

        if s.startswith("[") and s.endswith("]"):

            current = s[1:-1].strip()
            groups[current] = []

        elif current is not None:

            groups[current].extend(
                int(x) - 1
                for x in s.split()
            )

    return groups


def pbc_group_com(atomgroup, box):

    pos = atomgroup.positions.astype(np.float64)

    anchor = pos[0].copy()

    whole = (
        anchor
        + minimize_vectors(
            pos - anchor,
            box
        )
    )

    masses = atomgroup.masses.astype(np.float64)

    if (
        len(masses) != len(whole)
        or not np.all(np.isfinite(masses))
        or masses.sum() <= 0
    ):
        return whole.mean(axis=0)

    return np.average(
        whole,
        axis=0,
        weights=masses
    )


def densmap_geometry(u, idx_minus, idx_plus):

    box = u.dimensions.copy()

    cminus = pbc_group_com(
        u.atoms[idx_minus],
        box
    )

    cplus = pbc_group_com(
        u.atoms[idx_plus],
        box
    )

    avec = minimize_vectors(
        np.asarray(
            [cplus - cminus],
            dtype=np.float64
        ),
        box
    )[0]

    L = np.linalg.norm(avec)

    axis = avec / L

    center = cminus + 0.5 * avec

    return center, axis, L / 10.0


def relative_pbc_positions(atomgroup, center, box):

    return minimize_vectors(
        atomgroup.positions.astype(np.float64) - center,
        box
    )


def kabsch(current, reference):

    """
    Row-vector convention:
        current @ Q ~= reference
    Both arrays must already be centered.
    """

    H = current.T @ reference

    U, s, Vt = np.linalg.svd(H)

    Q = U @ Vt

    if np.linalg.det(Q) < 0:

        Vt[-1, :] *= -1

        Q = U @ Vt

    return Q


def linear_fit(x, y):

    p = np.polyfit(x, y, 1)

    pred = np.polyval(p, x)

    ss_res = float(
        np.sum((y - pred) ** 2)
    )

    ss_tot = float(
        np.sum((y - np.mean(y)) ** 2)
    )

    r2 = (
        1.0 - ss_res / ss_tot
        if ss_tot > 0
        else float("nan")
    )

    return (
        float(p[0]),
        float(p[1]),
        r2
    )


# =============================================================================
# [1] Input audit
# =============================================================================

print("\n[1] INPUT AUDIT")
print("-" * 120)

for label, p in [
    ("TPR", TPR),
    ("XTC", XTC),
    ("NDX", NDX),
    ("BOUNDARIES", BOUND),
]:

    print(
        f"{label}="
        + ("FOUND" if p.exists() else "NOT_FOUND")
        + f" path={p.relative_to(ROOT)}"
    )

    if not p.exists():
        raise SystemExit(
            f"F43A5_ABORT=MISSING_{label}"
        )


groups = read_ndx(NDX)

for g in [
    "AxisMinus",
    "AxisPlus",
    "HBN",
    "Water_O",
]:

    if g not in groups:
        raise SystemExit(
            f"F43A5_ABORT=MISSING_GROUP_{g}"
        )

    print(
        f"GROUP={g} N={len(groups[g])}"
    )


IDX_MINUS = np.asarray(
    groups["AxisMinus"],
    dtype=int
)

IDX_PLUS = np.asarray(
    groups["AxisPlus"],
    dtype=int
)

IDX_HBN = np.asarray(
    groups["HBN"],
    dtype=int
)

IDX_WATER = np.asarray(
    groups["Water_O"],
    dtype=int
)


# =============================================================================
# [2] Boundary contract
# =============================================================================

print("\n[2] PROFILE-GUIDED REGION CONTRACT")
print("-" * 120)

with BOUND.open(
    newline="",
    errors="ignore"
) as fh:

    rows = list(csv.DictReader(fh))


radial = next(
    x for x in rows
    if x["boundary_type"] == "radial"
)

axial = next(
    x for x in rows
    if x["boundary_type"] == "axial"
)


R_INNER = float(
    radial["inner_boundary_nm"]
)

R_MINIMUM = float(
    radial["minimum_radius_nm"]
)

R_OUTER = float(
    radial["outer_boundary_nm"]
)

ZLI = float(
    axial["left_inner_boundary_nm"]
)

ZRI = float(
    axial["right_inner_boundary_nm"]
)

ZHO = float(
    axial["HBN_lower_boundary_nm"]
)

ZHP = float(
    axial["HBN_upper_boundary_nm"]
)


print(f"PRIMARY_R_MAX_NM={R_INNER:.9f}")
print(f"PRIMARY_Z_MIN_NM={ZLI:.9f}")
print(f"PRIMARY_Z_MAX_NM={ZRI:.9f}")

print()
print(
    "PRIMARY_REGION="
    "PROFILE_GUIDED_LUMEN_CORE"
)

print(
    "SECONDARY_SENSITIVITY_REGION="
    "RADIAL_DEPLETION_MINIMUM_OVER_HBN_SPAN"
)


# =============================================================================
# [3] Trajectory
# =============================================================================

print("\n[3] TRAJECTORY")
print("-" * 120)

u = mda.Universe(
    str(TPR),
    str(XTC)
)

NF = len(u.trajectory)
NW = len(IDX_WATER)

print(f"N_FRAMES={NF}")
print(f"N_WATER_O={NW}")

times = np.zeros(
    NF,
    dtype=float
)

# Primary and sensitivity occupancy
occ_primary = np.zeros(
    (NF, NW),
    dtype=bool
)

occ_sensitivity = np.zeros(
    (NF, NW),
    dtype=bool
)

# Full body-frame coordinates for waters that ever occupy
# the primary region will be collected after first pass.


# =============================================================================
# [4] First frame reference scaffold
# =============================================================================

u.trajectory[0]

center0, axis0, axisL0 = densmap_geometry(
    u,
    IDX_MINUS,
    IDX_PLUS
)

box0 = u.dimensions.copy()

hbn_ref = (
    relative_pbc_positions(
        u.atoms[IDX_HBN],
        center0,
        box0
    )
)

hbn_ref -= hbn_ref.mean(axis=0)

print("\n[4] REFERENCE FRAME")
print("-" * 120)

print(
    "CENTER0_A="
    f"({center0[0]:.6f},"
    f"{center0[1]:.6f},"
    f"{center0[2]:.6f})"
)

print(
    "AXIS0="
    f"({axis0[0]:+.9f},"
    f"{axis0[1]:+.9f},"
    f"{axis0[2]:+.9f})"
)

print(
    f"AXIS_LENGTH0_NM="
    f"{axisL0:.9f}"
)


# =============================================================================
# [5] Occupancy pass
# =============================================================================

print("\n[5] OCCUPANCY PASS")
print("-" * 120)

axis_series = []
center_series = []

for iframe, ts in enumerate(u.trajectory):

    center, axis, axisL = densmap_geometry(
        u,
        IDX_MINUS,
        IDX_PLUS
    )

    box = u.dimensions.copy()

    dr = relative_pbc_positions(
        u.atoms[IDX_WATER],
        center,
        box
    ) / 10.0

    z = dr @ axis

    perp = (
        dr
        - np.outer(z, axis)
    )

    r = np.linalg.norm(
        perp,
        axis=1
    )

    primary = (
        (r <= R_INNER)
        &
        (z >= ZLI)
        &
        (z <= ZRI)
    )

    sensitivity = (
        (r <= R_MINIMUM)
        &
        (z >= ZHO)
        &
        (z <= ZHP)
    )

    occ_primary[iframe] = primary
    occ_sensitivity[iframe] = sensitivity

    times[iframe] = float(ts.time)

    axis_series.append(axis)
    center_series.append(center)

    if (
        iframe == 0
        or iframe == NF - 1
        or iframe % 20 == 0
    ):

        print(
            f"FRAME={iframe:3d} "
            f"T={times[iframe]:8.3f}ps "
            f"Ncore={primary.sum():4d} "
            f"Nsensitivity={sensitivity.sum():4d}"
        )


DT = float(
    np.median(
        np.diff(times)
    )
)

Ncore = occ_primary.sum(axis=1)
Nsens = occ_sensitivity.sum(axis=1)

print()
print(f"FRAME_DT_PS={DT:.9f}")

print(
    f"PRIMARY_MEAN_OCCUPANCY="
    f"{Ncore.mean():.6f}"
)

print(
    f"PRIMARY_SD_OCCUPANCY="
    f"{Ncore.std(ddof=1):.6f}"
)

print(
    f"PRIMARY_MIN_OCCUPANCY="
    f"{Ncore.min()}"
)

print(
    f"PRIMARY_MAX_OCCUPANCY="
    f"{Ncore.max()}"
)

print(
    f"SENSITIVITY_MEAN_OCCUPANCY="
    f"{Nsens.mean():.6f}"
)

# Regression against F43A4
if abs(Ncore.mean() - 329.482587) <= 0.05:

    print(
        "F43A4_OCCUPANCY_REGRESSION=PASS"
    )

else:

    print(
        "F43A4_OCCUPANCY_REGRESSION=FAIL"
    )

    raise SystemExit(
        "F43A5_ABORT=OCCUPANCY_REGRESSION_FAILURE"
    )


# =============================================================================
# [6] Residence correlation
# =============================================================================

print("\n[6] RESIDENCE CORRELATION")
print("-" * 120)

MAXLAG = NF // 2

lags = np.arange(
    MAXLAG + 1
)

lag_ps = lags * DT

Cint = np.full(
    MAXLAG + 1,
    np.nan
)

Scont = np.full(
    MAXLAG + 1,
    np.nan
)

denom = np.zeros(
    MAXLAG + 1,
    dtype=np.int64
)

Cint[0] = 1.0
Scont[0] = 1.0

denom[0] = int(
    occ_primary.sum()
)


# Continuous survival computed through cumulative outside counts.
outside = (
    ~occ_primary
).astype(np.int32)

outside_cs = np.vstack([
    np.zeros(
        (1, NW),
        dtype=np.int32
    ),
    np.cumsum(
        outside,
        axis=0,
        dtype=np.int32
    )
])


for lag in range(1, MAXLAG + 1):

    h0 = occ_primary[:-lag]
    ht = occ_primary[lag:]

    d = int(h0.sum())
    denom[lag] = d

    if d == 0:
        continue

    # Intermittent correlation
    Cint[lag] = (
        np.logical_and(
            h0,
            ht
        ).sum()
        / d
    )

    # Continuous occupancy of every sampled frame from t0..t0+lag.
    outside_count = (
        outside_cs[lag+1:]
        - outside_cs[:-(lag+1)]
    )

    survived = (
        (outside_count == 0)
        &
        h0
    )

    Scont[lag] = (
        survived.sum()
        / d
    )


def first_crossing(y, threshold):

    idx = np.where(
        np.isfinite(y)
        & (y <= threshold)
    )[0]

    return (
        float(lag_ps[idx[0]])
        if len(idx)
        else None
    )


tau_int = float(
    np.trapezoid(
        Cint[np.isfinite(Cint)],
        lag_ps[np.isfinite(Cint)]
    )
)

tau_cont = float(
    np.trapezoid(
        Scont[np.isfinite(Scont)],
        lag_ps[np.isfinite(Scont)]
    )
)

C_1e = first_crossing(
    Cint,
    1.0 / math.e
)

S_1e = first_crossing(
    Scont,
    1.0 / math.e
)


print(
    f"C_RES_INTEGRAL_0_{lag_ps[-1]:.1f}PS="
    f"{tau_int:.9f} ps"
)

print(
    f"S_CONT_INTEGRAL_0_{lag_ps[-1]:.1f}PS="
    f"{tau_cont:.9f} ps"
)

print(
    "C_RES_1E_PS="
    + (
        f"{C_1e:.9f}"
        if C_1e is not None
        else "NOT_REACHED"
    )
)

print(
    "S_CONT_1E_PS="
    + (
        f"{S_1e:.9f}"
        if S_1e is not None
        else "NOT_REACHED"
    )
)

print(
    f"C_RES_50PS="
    f"{Cint[-1]:.9f}"
)

print(
    f"S_CONT_50PS="
    f"{Scont[-1]:.9f}"
)


# =============================================================================
# [7] Continuous residence events
# =============================================================================

print("\n[7] CONTINUOUS RESIDENCE EVENTS")
print("-" * 120)

event_durations = []

for j in range(NW):

    h = occ_primary[:,j]

    if not np.any(h):
        continue

    starts = np.where(
        h
        &
        np.r_[
            True,
            ~h[:-1]
        ]
    )[0]

    ends = np.where(
        h
        &
        np.r_[
            ~h[1:],
            True
        ]
    )[0]

    for s, e in zip(starts, ends):

        event_durations.append(
            (e - s + 1) * DT
        )


event_durations = np.asarray(
    event_durations,
    dtype=float
)

print(
    f"CONTINUOUS_EVENTS="
    f"{len(event_durations)}"
)

if len(event_durations):

    print(
        f"EVENT_MEAN_PS="
        f"{event_durations.mean():.9f}"
    )

    print(
        f"EVENT_MEDIAN_PS="
        f"{np.median(event_durations):.9f}"
    )

    print(
        f"EVENT_P95_PS="
        f"{np.percentile(event_durations,95):.9f}"
    )

    print(
        f"EVENT_MAX_PS="
        f"{event_durations.max():.9f}"
    )


# =============================================================================
# [8] Prepare body-frame coordinates for waters ever in core
# =============================================================================

print("\n[8] BODY-FRAME COORDINATE PASS")
print("-" * 120)

ever_local = np.where(
    occ_primary.any(axis=0)
)[0]

ever_global = IDX_WATER[
    ever_local
]

NEVER = len(ever_local)

print(
    f"WATERS_EVER_IN_PRIMARY_CORE="
    f"{NEVER}/{NW}"
)

coords = np.empty(
    (NF, NEVER, 3),
    dtype=np.float32
)

occ_e = occ_primary[
    :,
    ever_local
]


for iframe, ts in enumerate(u.trajectory):

    center, axis, axisL = densmap_geometry(
        u,
        IDX_MINUS,
        IDX_PLUS
    )

    box = u.dimensions.copy()

    # Reconstruct scaffold relative to instantaneous densmap center
    hbn = relative_pbc_positions(
        u.atoms[IDX_HBN],
        center,
        box
    )

    hbn -= hbn.mean(axis=0)

    Q = kabsch(
        hbn,
        hbn_ref
    )

    water = relative_pbc_positions(
        u.atoms[ever_global],
        center,
        box
    )

    # Remove scaffold rigid-body rotation.
    water_body = water @ Q

    coords[iframe] = (
        water_body / 10.0
    ).astype(np.float32)

    if (
        iframe == 0
        or iframe == NF - 1
        or iframe % 40 == 0
    ):

        rms = np.sqrt(
            np.mean(
                np.sum(
                    (
                        hbn @ Q
                        - hbn_ref
                    ) ** 2,
                    axis=1
                )
            )
        )

        print(
            f"FRAME={iframe:3d} "
            f"HBN_ALIGNMENT_RMS_A={rms:.6f}"
        )


# Reference axis after rigid alignment.
axis_ref = axis0.copy()


# =============================================================================
# [9] Continuous-confinement MSD
# =============================================================================

print("\n[9] CONTINUOUS-CONFINEMENT MSD")
print("-" * 120)

outside_e = (
    ~occ_e
).astype(np.int32)

cs = np.vstack([
    np.zeros(
        (1, NEVER),
        dtype=np.int32
    ),
    np.cumsum(
        outside_e,
        axis=0,
        dtype=np.int32
    )
])


msd_rows = []

for lag in range(
    1,
    MAXLAG + 1
):

    continuous = (
        cs[lag+1:]
        - cs[:-(lag+1)]
    ) == 0

    n_samples = int(
        continuous.sum()
    )

    if n_samples == 0:

        msd_rows.append(
            {
                "lag_frames":lag,
                "lag_ps":lag*DT,
                "samples":0,
                "MSD_3D_nm2":None,
                "MSD_axial_nm2":None,
                "MSD_transverse_nm2":None,
            }
        )

        continue

    d = (
        coords[lag:]
        - coords[:-lag]
    )

    d2 = np.sum(
        d*d,
        axis=2
    )

    dz = (
        d @ axis_ref
    )

    dz2 = dz * dz

    dperp2 = d2 - dz2

    msd_rows.append(
        {
            "lag_frames":lag,
            "lag_ps":lag*DT,
            "samples":n_samples,

            "MSD_3D_nm2":
                float(
                    d2[continuous].mean()
                ),

            "MSD_axial_nm2":
                float(
                    dz2[continuous].mean()
                ),

            "MSD_transverse_nm2":
                float(
                    dperp2[continuous].mean()
                ),
        }
    )


valid = [
    x for x in msd_rows
    if (
        x["samples"] >= 100
        and x["MSD_3D_nm2"] is not None
    )
]

print(
    f"MSD_VALID_LAGS="
    f"{len(valid)}/{MAXLAG}"
)

if valid:

    print(
        f"MSD_MAX_VALID_LAG_PS="
        f"{valid[-1]['lag_ps']:.6f}"
    )

    print(
        f"MSD_MAX_VALID_LAG_SAMPLES="
        f"{valid[-1]['samples']}"
    )


# =============================================================================
# [10] Diffusive-regime scan
# =============================================================================

print("\n[10] DIFFUSIVE-REGIME SCAN")
print("-" * 120)

fit_windows = []

if len(valid) >= 10:

    t = np.asarray(
        [x["lag_ps"] for x in valid]
    )

    y3 = np.asarray(
        [x["MSD_3D_nm2"] for x in valid]
    )

    yz = np.asarray(
        [x["MSD_axial_nm2"] for x in valid]
    )

    yp = np.asarray(
        [x["MSD_transverse_nm2"] for x in valid]
    )

    candidate_windows = [
        (0.10,0.30),
        (0.20,0.40),
        (0.25,0.50),
        (0.33,0.66),
        (0.50,0.80),
    ]

    for f0, f1 in candidate_windows:

        lo = (
            t.min()
            + f0 * (t.max()-t.min())
        )

        hi = (
            t.min()
            + f1 * (t.max()-t.min())
        )

        m = (
            (t >= lo)
            &
            (t <= hi)
        )

        if m.sum() < 5:
            continue

        s3,b3,r3 = linear_fit(
            t[m],
            y3[m]
        )

        sz,bz,rz = linear_fit(
            t[m],
            yz[m]
        )

        sp,bp,rp = linear_fit(
            t[m],
            yp[m]
        )

        rec = {
            "fraction_start":f0,
            "fraction_end":f1,

            "t_start_ps":
                float(t[m][0]),

            "t_end_ps":
                float(t[m][-1]),

            "n_points":
                int(m.sum()),

            "slope_3D_nm2_ps":s3,
            "R2_3D":r3,

            "D_3D_nm2_ps":
                s3 / 6.0,

            "D_3D_m2_s":
                s3 / 6.0 * 1e-6,

            "slope_axial_nm2_ps":sz,
            "R2_axial":rz,

            "D_axial_nm2_ps":
                sz / 2.0,

            "D_axial_m2_s":
                sz / 2.0 * 1e-6,

            "slope_transverse_nm2_ps":sp,
            "R2_transverse":rp,

            "D_transverse_nm2_ps":
                sp / 4.0,

            "D_transverse_m2_s":
                sp / 4.0 * 1e-6,
        }

        fit_windows.append(rec)

        print(
            f"WINDOW="
            f"{rec['t_start_ps']:.3f}-"
            f"{rec['t_end_ps']:.3f}ps "
            f"N={rec['n_points']} "
            f"D3={rec['D_3D_m2_s']:.6e} "
            f"R2_3D={rec['R2_3D']:.5f} "
            f"Dz={rec['D_axial_m2_s']:.6e} "
            f"R2z={rec['R2_axial']:.5f} "
            f"Dperp={rec['D_transverse_m2_s']:.6e} "
            f"R2perp={rec['R2_transverse']:.5f}"
        )


acceptable = [
    x for x in fit_windows
    if (
        x["D_3D_m2_s"] > 0
        and x["D_axial_m2_s"] > 0
        and x["D_transverse_m2_s"] > 0
        and x["R2_3D"] >= 0.95
        and x["R2_axial"] >= 0.90
        and x["R2_transverse"] >= 0.90
    )
]


if acceptable:

    acceptable.sort(
        key=lambda x: (
            x["R2_3D"]
            + x["R2_axial"]
            + x["R2_transverse"]
        ),
        reverse=True
    )

    selected = acceptable[0]

    diffusion_status = (
        "LINEAR_REGIME_DETECTED"
    )

    print()
    print(
        "DIFFUSION_WINDOW_SELECTED="
        f"{selected['t_start_ps']:.3f}-"
        f"{selected['t_end_ps']:.3f} ps"
    )

    print(
        f"D_3D_M2_S="
        f"{selected['D_3D_m2_s']:.9e}"
    )

    print(
        f"D_PARALLEL_M2_S="
        f"{selected['D_axial_m2_s']:.9e}"
    )

    print(
        f"D_PERP_M2_S="
        f"{selected['D_transverse_m2_s']:.9e}"
    )

else:

    selected = None

    diffusion_status = (
        "NO_ROBUST_LINEAR_REGIME"
    )


print(
    f"DIFFUSION_STATUS="
    f"{diffusion_status}"
)


# =============================================================================
# [11] Sensitivity of occupancy
# =============================================================================

print("\n[11] REGION-SENSITIVITY SUMMARY")
print("-" * 120)

print(
    f"PROFILE_CORE_MEAN="
    f"{Ncore.mean():.6f}"
)

print(
    f"DEPLETION_MINIMUM_HBN_SPAN_MEAN="
    f"{Nsens.mean():.6f}"
)

print(
    f"POPULATION_RATIO_SENSITIVITY_TO_PRIMARY="
    f"{Nsens.mean()/Ncore.mean():.6f}"
)


# =============================================================================
# [12] Write outputs
# =============================================================================

print("\n[12] OUTPUT")
print("-" * 120)

RESCSV = (
    OUT /
    "M2_4_CONFINED_WATER_RESIDENCE.csv"
)

with RESCSV.open(
    "w",
    newline=""
) as fh:

    w = csv.writer(fh)

    w.writerow([
        "lag_frames",
        "lag_ps",
        "C_res_intermittent",
        "S_res_continuous",
        "origin_occupancy_denominator",
    ])

    for i in range(
        MAXLAG + 1
    ):

        w.writerow([
            int(lags[i]),
            float(lag_ps[i]),
            float(Cint[i])
                if np.isfinite(Cint[i])
                else "",
            float(Scont[i])
                if np.isfinite(Scont[i])
                else "",
            int(denom[i]),
        ])


EVENTCSV = (
    OUT /
    "M2_4_CONTINUOUS_RESIDENCE_EVENTS.csv"
)

with EVENTCSV.open(
    "w",
    newline=""
) as fh:

    w = csv.writer(fh)

    w.writerow([
        "event_id",
        "duration_ps"
    ])

    for i, x in enumerate(
        event_durations,
        1
    ):

        w.writerow([
            i,
            float(x)
        ])


MSDCSV = (
    OUT /
    "M2_5_CONFINED_WATER_MSD.csv"
)

with MSDCSV.open(
    "w",
    newline=""
) as fh:

    fields = [
        "lag_frames",
        "lag_ps",
        "samples",
        "MSD_3D_nm2",
        "MSD_axial_nm2",
        "MSD_transverse_nm2",
    ]

    w = csv.DictWriter(
        fh,
        fieldnames=fields
    )

    w.writeheader()

    for x in msd_rows:
        w.writerow(x)


FITCSV = (
    OUT /
    "M2_5_DIFFUSION_WINDOW_SCAN.csv"
)

with FITCSV.open(
    "w",
    newline=""
) as fh:

    if fit_windows:

        fields = list(
            fit_windows[0].keys()
        )

        w = csv.DictWriter(
            fh,
            fieldnames=fields
        )

        w.writeheader()
        w.writerows(fit_windows)


summary = {

    "phase":
        "PHASE5-F43A5",

    "coordinate_contract":
        "F43A4 validated AxisMinus/AxisPlus densmap PBC",

    "primary_confined_region": {
        "name":
            "profile-guided lumen core",

        "r_max_nm":
            R_INNER,

        "z_min_nm":
            ZLI,

        "z_max_nm":
            ZRI,
    },

    "trajectory": {
        "frames":NF,
        "dt_ps":DT,
        "duration_ps":
            float(
                times[-1]
                - times[0]
            ),
    },

    "occupancy": {
        "mean":
            float(Ncore.mean()),

        "sd":
            float(
                Ncore.std(ddof=1)
            ),

        "min":
            int(Ncore.min()),

        "max":
            int(Ncore.max()),

        "waters_ever_in_core":
            int(NEVER),

        "sensitivity_mean":
            float(Nsens.mean()),
    },

    "residence": {
        "C_res_definition":
            "<h(t0)h(t0+t)>/<h(t0)>",

        "S_cont_definition":
            "continuous occupancy from t0 through t0+t",

        "integration_window_ps":
            float(lag_ps[-1]),

        "C_res_integral_ps":
            tau_int,

        "S_cont_integral_ps":
            tau_cont,

        "C_res_1e_ps":
            C_1e,

        "S_cont_1e_ps":
            S_1e,

        "continuous_event_count":
            int(
                len(event_durations)
            ),

        "event_mean_ps":
            float(
                event_durations.mean()
            ) if len(event_durations)
            else None,

        "event_median_ps":
            float(
                np.median(
                    event_durations
                )
            ) if len(event_durations)
            else None,

        "event_p95_ps":
            float(
                np.percentile(
                    event_durations,
                    95
                )
            ) if len(event_durations)
            else None,

        "event_max_ps":
            float(
                event_durations.max()
            ) if len(event_durations)
            else None,
    },

    "diffusion": {
        "conditioning":
            "continuous confinement over each lag interval",

        "status":
            diffusion_status,

        "selected":
            selected,

        "windows":
            fit_windows,
    },

    "limitations": [
        "Trajectory duration is 100 ps.",
        "Trajectory output interval is 0.5 ps; residence events shorter than the output interval are unresolved.",
        "Diffusion is conditioned on continuous occupancy of the confined region and therefore characterizes retained confined water rather than unrestricted bulk water.",
        "Residence and diffusion estimates are not declared converged if the required decay/linear regime is absent.",
    ],

    "new_md_simulations":
        0,

    "source_files_modified":
        False,
}


JOUT = (
    OUT /
    "PHASE5_F43A5_CONFINED_WATER_TRANSPORT.json"
)

JOUT.write_text(
    json.dumps(
        summary,
        indent=2
    )
)


for p in [
    RESCSV,
    EVENTCSV,
    MSDCSV,
    FITCSV,
    JOUT,
]:

    print(
        f"WROTE="
        f"{p.relative_to(ROOT)}"
    )


# =============================================================================
# [13] Final gate
# =============================================================================

print("\n[13] FINAL SCIENTIFIC GATE")
print("-" * 120)

print(
    "M2_4_RESIDENCE_CORRELATION="
    "COMPUTED"
)

if (
    C_1e is not None
    and S_1e is not None
):

    print(
        "M2_4_RESIDENCE_STATUS="
        "DECAY_RESOLVED_WITHIN_TRAJECTORY"
    )

else:

    print(
        "M2_4_RESIDENCE_STATUS="
        "FINITE_WINDOW_ONLY"
    )


if selected is not None:

    print(
        "M2_5_DIFFUSION_STATUS="
        "LINEAR_REGIME_CANDIDATE_IDENTIFIED"
    )

else:

    print(
        "M2_5_DIFFUSION_STATUS="
        "NO_DEFENSIBLE_LONG_TIME_D_FOUND"
    )


print(
    "M2_6_HYDROGEN_BONDS="
    "PENDING_F43B"
)

print(
    "NEW_MD_SIMULATIONS=0"
)

print(
    "SOURCE_FILES_MODIFIED=NO"
)

print(
    "F43A5_STATUS=0"
)

print(
    "NEXT_ACTION="
    "SCIENTIFIC_INTERPRETATION_AND_F43B_HBOND_GATE"
)

print("\n" + "=" * 152)
print("PHASE5-F43A5 COMPLETE")
print("TERMINAL_REMAINS_OPEN=YES")
print("=" * 152)

