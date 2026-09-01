#!/usr/bin/env python3

from pathlib import Path
import csv
import json
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import minimize_vectors

print("=" * 148)
print("PHASE5-F43A4 — EXACT DENSMAP AXIS / PBC COORDINATE-CONTRACT VALIDATION")
print("PURPOSE=REPRODUCE_DAY020_GEOMETRY_BEFORE_RESIDENCE_OR_DIFFUSION")
print("=" * 148)

ROOT = Path.cwd()
P0 = ROOT / "Microtubule_Inspired_Phase0"

DAY020 = (
    P0 /
    "runs/phase1A/day020_confined_water_axial_radial_density"
)

NDX = DAY020 / "confined_water_analysis_groups.ndx"

BOUND = (
    DAY020 /
    "profile_guided_classification/profile_guided_boundaries.csv"
)

PROFILE_SUMMARY = (
    DAY020 /
    "profile_guided_classification/profile_guided_region_summary.csv"
)

REGIONAL_SUMMARY = (
    DAY020 /
    "regional_classification/confined_water_region_summary.csv"
)

FROZEN = (
    P0 /
    "runs/phase1A/accepted/"
    "hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute"
)

MOBILE = (
    P0 /
    "runs/phase1A/day021_mobile_restraint_protocol/"
    "execution/08_nvt_mobile_100ps"
)

OUT = (
    P0 /
    "runs/phase2/campaign_phase5_corrected/"
    "project_closure_control/"
    "f43a4_densmap_coordinate_contract"
)

OUT.mkdir(parents=True, exist_ok=True)

# =============================================================================
# Helpers
# =============================================================================

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

        elif current is not None:

            groups[current].extend(
                int(x) - 1
                for x in line.split()
            )

    return groups


def print_csv(path, title):

    print()
    print(title)
    print("-" * 120)

    if not path.exists():
        print("FILE=NOT_FOUND")
        return

    print(f"FILE={path.relative_to(ROOT)}")

    with path.open(newline="", errors="ignore") as fh:

        rows = list(csv.DictReader(fh))

    print(f"ROWS={len(rows)}")

    for row in rows:
        print(row)


def discover_one(directory, pattern):

    hits = sorted(directory.glob(pattern))

    if len(hits) != 1:

        print(f"DISCOVERY_DIR={directory}")
        print(f"PATTERN={pattern}")
        print(f"HITS={len(hits)}")

        for p in hits:
            print(p)

        raise SystemExit(
            "F43A4_ABORT=AMBIGUOUS_INPUT_DISCOVERY"
        )

    return hits[0]


def pbc_group_com(atomgroup, box):

    """
    PBC-aware COM for a spatially compact reference group.

    Positions are first reconstructed around one anchor using
    the minimum-image convention, then mass weighted.
    """

    pos = atomgroup.positions.astype(np.float64)

    anchor = pos[0].copy()

    vec = pos - anchor

    vec = minimize_vectors(
        vec,
        box
    )

    whole = anchor + vec

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


def densmap_coordinates(u, idx_minus, idx_plus, idx_water):

    """
    Reproduce gmx densmap axial-radial coordinate definition:

      COM(AxisMinus) -> COM(AxisPlus)
      midpoint = origin
      positive axis = Minus -> Plus

    Water displacements relative to the midpoint are reduced
    with the minimum-image convention.
    """

    box = u.dimensions.copy()

    minus = u.atoms[idx_minus]
    plus  = u.atoms[idx_plus]
    water = u.atoms[idx_water]

    cminus = pbc_group_com(
        minus,
        box
    )

    cplus = pbc_group_com(
        plus,
        box
    )

    axis_vec = minimize_vectors(
        np.asarray(
            [cplus - cminus],
            dtype=np.float64
        ),
        box
    )[0]

    axis_length_A = np.linalg.norm(axis_vec)

    if axis_length_A <= 0:
        raise RuntimeError(
            "Zero AxisMinus-AxisPlus separation"
        )

    axis = axis_vec / axis_length_A

    center = (
        cminus
        + 0.5 * axis_vec
    )

    dr = (
        water.positions.astype(np.float64)
        - center
    )

    dr = minimize_vectors(
        dr,
        box
    )

    z_A = dr @ axis

    perp = (
        dr
        - np.outer(z_A, axis)
    )

    r_A = np.linalg.norm(
        perp,
        axis=1
    )

    return {
        "r_nm": r_A / 10.0,
        "z_nm": z_A / 10.0,
        "axis": axis,
        "center_A": center,
        "axis_length_nm":
            axis_length_A / 10.0,
        "com_minus_A": cminus,
        "com_plus_A": cplus,
    }


# =============================================================================
# [1] Input audit
# =============================================================================

print("\n[1] INPUT AUDIT")
print("-" * 120)

groups = read_ndx(NDX)

required = [
    "AxisMinus",
    "AxisPlus",
    "Water_O",
    "HBN",
]

for g in required:

    if g not in groups:
        raise SystemExit(
            f"F43A4_ABORT=MISSING_GROUP_{g}"
        )

    print(
        f"GROUP={g} "
        f"N={len(groups[g])}"
    )


FROZEN_TPR = discover_one(
    FROZEN,
    "*frozenSolute.tpr"
)

FROZEN_XTC = discover_one(
    FROZEN,
    "*frozenSolute.xtc"
)

MOBILE_TPR = discover_one(
    MOBILE,
    "*08_nvt_mobile_100ps*.tpr"
)

MOBILE_XTC = discover_one(
    MOBILE,
    "*08_nvt_mobile_100ps*.xtc"
)

for label, p in [
    ("FROZEN_TPR", FROZEN_TPR),
    ("FROZEN_XTC", FROZEN_XTC),
    ("MOBILE_TPR", MOBILE_TPR),
    ("MOBILE_XTC", MOBILE_XTC),
]:

    print(
        f"{label}="
        f"{p.relative_to(ROOT)}"
    )


# =============================================================================
# [2] Canonical Day020 products
# =============================================================================

print_csv(
    PROFILE_SUMMARY,
    "[2A] DAY020 PROFILE-GUIDED REGION SUMMARY"
)

print_csv(
    REGIONAL_SUMMARY,
    "[2B] DAY020 OPERATIONAL REGION SUMMARY"
)


# =============================================================================
# [3] Parse profile-guided boundaries
# =============================================================================

print("\n[3] PROFILE-GUIDED BOUNDARY CONTRACT")
print("-" * 120)

with BOUND.open(
    newline="",
    errors="ignore"
) as fh:

    rows = list(csv.DictReader(fh))

radial = next(
    r for r in rows
    if r["boundary_type"] == "radial"
)

axial = next(
    r for r in rows
    if r["boundary_type"] == "axial"
)

R_MIN = float(
    radial["minimum_radius_nm"]
)

R_INNER = float(
    radial["inner_boundary_nm"]
)

R_OUTER = float(
    radial["outer_boundary_nm"]
)

Z_LEFT_OUTER = float(
    axial["left_outer_boundary_nm"]
)

Z_LEFT_INNER = float(
    axial["left_inner_boundary_nm"]
)

Z_RIGHT_INNER = float(
    axial["right_inner_boundary_nm"]
)

Z_RIGHT_OUTER = float(
    axial["right_outer_boundary_nm"]
)

Z_HBN_LOW = float(
    axial["HBN_lower_boundary_nm"]
)

Z_HBN_HIGH = float(
    axial["HBN_upper_boundary_nm"]
)

for k, v in [
    ("R_MIN_NM", R_MIN),
    ("R_INNER_NM", R_INNER),
    ("R_OUTER_NM", R_OUTER),
    ("Z_LEFT_OUTER_NM", Z_LEFT_OUTER),
    ("Z_LEFT_INNER_NM", Z_LEFT_INNER),
    ("Z_RIGHT_INNER_NM", Z_RIGHT_INNER),
    ("Z_RIGHT_OUTER_NM", Z_RIGHT_OUTER),
    ("Z_HBN_LOW_NM", Z_HBN_LOW),
    ("Z_HBN_HIGH_NM", Z_HBN_HIGH),
]:

    print(f"{k}={v:.9f}")


# =============================================================================
# Selectors
# =============================================================================

def selectors(r, z):

    return {

        # Strict profile-guided interior.
        "PROFILE_INTERIOR_CORE":
            (
                (r <= R_INNER)
                &
                (z >= Z_LEFT_INNER)
                &
                (z <= Z_RIGHT_INNER)
            ),

        # Interior radius over physical HBN axial segment.
        "INNER_RADIUS_HBN_SPAN":
            (
                (r <= R_INNER)
                &
                (z >= Z_HBN_LOW)
                &
                (z <= Z_HBN_HIGH)
            ),

        # Density-minimum radius over central profile region.
        "MINIMUM_RADIUS_CENTRAL":
            (
                (r <= R_MIN)
                &
                (z >= Z_LEFT_INNER)
                &
                (z <= Z_RIGHT_INNER)
            ),

        # Density-minimum radius over full HBN span.
        "MINIMUM_RADIUS_HBN_SPAN":
            (
                (r <= R_MIN)
                &
                (z >= Z_HBN_LOW)
                &
                (z <= Z_HBN_HIGH)
            ),

        # Interfacial radial shell over central profile region.
        "PROFILE_INTERFACIAL_SHELL":
            (
                (r > R_INNER)
                &
                (r <= R_OUTER)
                &
                (z >= Z_LEFT_INNER)
                &
                (z <= Z_RIGHT_INNER)
            ),

        # Entire analyzed profile cylinder.
        "PROFILE_ANALYSIS_CYLINDER":
            (
                (r <= R_OUTER)
                &
                (z >= Z_LEFT_OUTER)
                &
                (z <= Z_RIGHT_OUTER)
            ),
    }


# =============================================================================
# [4] Frozen trajectory — validation dataset
# =============================================================================

def analyze_trajectory(label, tpr, xtc):

    print()
    print(
        f"[TRAJECTORY] {label}"
    )
    print("-" * 120)

    u = mda.Universe(
        str(tpr),
        str(xtc)
    )

    print(
        f"N_FRAMES={len(u.trajectory)}"
    )

    idx_minus = np.asarray(
        groups["AxisMinus"],
        dtype=int
    )

    idx_plus = np.asarray(
        groups["AxisPlus"],
        dtype=int
    )

    idx_water = np.asarray(
        groups["Water_O"],
        dtype=int
    )

    all_counts = {}
    axes = []
    centers = []
    axis_lengths = []

    first = None

    for iframe, ts in enumerate(u.trajectory):

        geo = densmap_coordinates(
            u,
            idx_minus,
            idx_plus,
            idx_water
        )

        masks = selectors(
            geo["r_nm"],
            geo["z_nm"]
        )

        if first is None:
            first = {
                "geo": geo,
                "counts": {
                    k: int(v.sum())
                    for k, v in masks.items()
                }
            }

        for name, mask in masks.items():

            all_counts.setdefault(
                name,
                []
            ).append(
                int(mask.sum())
            )

        axes.append(
            geo["axis"]
        )

        centers.append(
            geo["center_A"]
        )

        axis_lengths.append(
            geo["axis_length_nm"]
        )

    axes = np.asarray(axes)
    centers = np.asarray(centers)
    axis_lengths = np.asarray(
        axis_lengths
    )

    print()
    print("FIRST_FRAME_DENSMAP_GEOMETRY")

    g = first["geo"]

    print(
        "CENTER_A="
        f"({g['center_A'][0]:.6f},"
        f"{g['center_A'][1]:.6f},"
        f"{g['center_A'][2]:.6f})"
    )

    print(
        "AXIS="
        f"({g['axis'][0]:+.9f},"
        f"{g['axis'][1]:+.9f},"
        f"{g['axis'][2]:+.9f})"
    )

    print(
        f"AXIS_LENGTH_NM="
        f"{g['axis_length_nm']:.9f}"
    )

    print()
    print("FIRST_FRAME_COUNTS")

    for name, n in first["counts"].items():
        print(f"{name}={n}")

    print()
    print("TIME_AVERAGED_COUNTS")

    summary = {}

    for name, values in all_counts.items():

        a = np.asarray(
            values,
            dtype=float
        )

        summary[name] = {
            "mean": float(a.mean()),
            "sd": float(
                a.std(ddof=1)
                if len(a) > 1
                else 0.0
            ),
            "min": int(a.min()),
            "max": int(a.max()),
        }

        print(
            f"{name} "
            f"MEAN={a.mean():.6f} "
            f"SD="
            f"{a.std(ddof=1):.6f} "
            f"MIN={int(a.min())} "
            f"MAX={int(a.max())}"
        )

    print()
    print("AXIS_STABILITY")

    print(
        "MEAN_AXIS="
        f"({axes[:,0].mean():+.9f},"
        f"{axes[:,1].mean():+.9f},"
        f"{axes[:,2].mean():+.9f})"
    )

    print(
        f"AXIS_LENGTH_MEAN_NM="
        f"{axis_lengths.mean():.9f}"
    )

    print(
        f"AXIS_LENGTH_SD_NM="
        f"{axis_lengths.std(ddof=1):.9f}"
    )

    print(
        "CENTER_MEAN_A="
        f"({centers[:,0].mean():.6f},"
        f"{centers[:,1].mean():.6f},"
        f"{centers[:,2].mean():.6f})"
    )

    return {
        "label": label,
        "n_frames":
            len(u.trajectory),
        "first_frame_counts":
            first["counts"],
        "counts": summary,
        "mean_axis":
            axes.mean(axis=0).tolist(),
        "mean_center_A":
            centers.mean(axis=0).tolist(),
        "mean_axis_length_nm":
            float(axis_lengths.mean()),
    }


print("\n[4] FROZEN — DAY020 RECONSTRUCTION")
print("=" * 120)

frozen_result = analyze_trajectory(
    "FROZEN_DAY020_REFERENCE",
    FROZEN_TPR,
    FROZEN_XTC
)


# =============================================================================
# [5] Mobile trajectory using identical coordinate contract
# =============================================================================

print("\n[5] MOBILE — SAME DENSMAP CONTRACT")
print("=" * 120)

mobile_result = analyze_trajectory(
    "MOBILE_DAY021",
    MOBILE_TPR,
    MOBILE_XTC
)


# =============================================================================
# [6] Scientific comparison
# =============================================================================

print("\n[6] CROSS-TRAJECTORY COMPARISON")
print("-" * 120)

for name in frozen_result["counts"]:

    f = frozen_result[
        "counts"
    ][name]["mean"]

    m = mobile_result[
        "counts"
    ][name]["mean"]

    print(
        f"{name} "
        f"FROZEN_MEAN={f:.6f} "
        f"MOBILE_MEAN={m:.6f} "
        f"DELTA={m-f:+.6f}"
    )


# =============================================================================
# [7] Sanity gate
# =============================================================================

print("\n[7] COORDINATE-CONTRACT SANITY GATE")
print("-" * 120)

f_core = frozen_result[
    "counts"
]["PROFILE_INTERIOR_CORE"]["mean"]

m_core = mobile_result[
    "counts"
]["PROFILE_INTERIOR_CORE"]["mean"]

print(
    f"FROZEN_PROFILE_CORE_MEAN="
    f"{f_core:.6f}"
)

print(
    f"MOBILE_PROFILE_CORE_MEAN="
    f"{m_core:.6f}"
)

# Previous broken selector was ~2.55 waters/frame.
# A correctly reconstructed lumen population should be
# orders of magnitude larger.

if f_core >= 100 and m_core >= 100:

    contract_pass = True

    print(
        "DENSMAP_COORDINATE_CONTRACT="
        "PLAUSIBLE_PASS"
    )

else:

    contract_pass = False

    print(
        "DENSMAP_COORDINATE_CONTRACT="
        "FAIL"
    )


# =============================================================================
# [8] Write audit
# =============================================================================

payload = {
    "phase":
        "PHASE5-F43A4",

    "coordinate_contract":
        "GROMACS densmap AxisMinus/AxisPlus COM midpoint axis",

    "pbc":
        "minimum-image displacement relative to dynamic densmap midpoint",

    "frozen":
        frozen_result,

    "mobile":
        mobile_result,

    "contract_pass":
        contract_pass,

    "new_md_simulations":
        0,

    "source_files_modified":
        False,
}

jout = (
    OUT /
    "PHASE5_F43A4_DENSMAP_COORDINATE_CONTRACT.json"
)

jout.write_text(
    json.dumps(
        payload,
        indent=2
    )
)

print()
print(f"WROTE={jout.relative_to(ROOT)}")


# =============================================================================
# [9] Final decision
# =============================================================================

print("\n[9] FINAL DECISION")
print("-" * 120)

if contract_pass:

    print(
        "F43A4_COORDINATE_CONTRACT=PASS"
    )

    print(
        "OLD_PCA_SELECTOR=REJECTED"
    )

    print(
        "F43A_RESIDENCE_DIFFUSION="
        "AUTHORIZED_FOR_REIMPLEMENTATION_WITH_DENSMAP_CONTRACT"
    )

    print(
        "NEXT_ACTION="
        "F43A5_RESIDENCE_AND_MSD_WITH_VALIDATED_COORDINATES"
    )

else:

    print(
        "F43A4_COORDINATE_CONTRACT=FAIL"
    )

    print(
        "RESIDENCE_DIFFUSION=REMAIN_BLOCKED"
    )

    print(
        "NEXT_ACTION="
        "COMPARE_PYTHON_COORDINATES_DIRECTLY_WITH_GROMACS_DENSMAP"
    )

print("NEW_MD_SIMULATIONS=0")
print("SOURCE_FILES_MODIFIED=NO")
print("F43A4_STATUS=0")

print("\n" + "=" * 148)
print("PHASE5-F43A4 COMPLETE")
print("TERMINAL_REMAINS_OPEN=YES")
print("=" * 148)

