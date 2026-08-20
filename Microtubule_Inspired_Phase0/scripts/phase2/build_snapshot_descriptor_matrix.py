#!/usr/bin/env python3

from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

CAMPAIGN = (
    ROOT
    / "runs"
    / "phase2"
    / "campaign_phase5_corrected"
)

ENSEMBLE = CAMPAIGN / "representative_ensemble"
OUT = CAMPAIGN / "snapshot_descriptors"

OUT.mkdir(parents=True, exist_ok=True)


def read_dump_snapshot(path):

    lines = path.read_text().splitlines()

    atom_header = None

    for i,line in enumerate(lines):

        if line.startswith("ITEM: ATOMS"):
            atom_header = i
            break

    if atom_header is None:
        raise RuntimeError(
            f"No ITEM: ATOMS header in {path}"
        )

    columns = lines[atom_header].split()[2:]

    required = ["id","xu","yu","zu"]

    missing = [
        c for c in required
        if c not in columns
    ]

    if missing:
        raise RuntimeError(
            f"{path}: missing required unwrapped "
            f"columns {missing}; found {columns}"
        )

    index = {
        name: columns.index(name)
        for name in required
    }

    rows = []

    for line in lines[atom_header+1:]:

        if not line.strip():
            continue

        if line.startswith("ITEM:"):
            break

        fields = line.split()

        rows.append(
            (
                int(fields[index["id"]]),
                float(fields[index["xu"]]),
                float(fields[index["yu"]]),
                float(fields[index["zu"]]),
            )
        )

    rows.sort(key=lambda row: row[0])

    ids = np.asarray(
        [r[0] for r in rows],
        dtype=int
    )

    xyz = np.asarray(
        [[r[1],r[2],r[3]] for r in rows],
        dtype=float
    )

    return ids, xyz


def center_coordinates(xyz):

    centroid = xyz.mean(axis=0)

    return xyz-centroid, centroid


def kabsch_align(mobile, reference):

    # Both coordinate sets must already be centered.

    covariance = mobile.T @ reference

    U, S, Vt = np.linalg.svd(covariance)

    rotation = U @ Vt

    # Prevent reflection.
    if np.linalg.det(rotation) < 0.0:

        U[:,-1] *= -1.0

        rotation = U @ Vt

    aligned = mobile @ rotation

    return aligned, rotation


snapshots = []

for folder in sorted(ENSEMBLE.glob("*K")):

    T = int(folder.name[:-1])

    for snap in sorted(folder.glob("snapshot_*.dump")):

        ids, xyz = read_dump_snapshot(snap)

        snapshots.append(
            {
                "temperature_K":T,
                "snapshot":snap.stem,
                "path":snap,
                "ids":ids,
                "xyz":xyz,
            }
        )


if not snapshots:
    raise RuntimeError(
        f"No snapshots found under {ENSEMBLE}"
    )


# ------------------------------------------------------------
# Common atom ordering
# ------------------------------------------------------------

reference_ids = snapshots[0]["ids"]

for item in snapshots[1:]:

    if not np.array_equal(
        item["ids"],
        reference_ids
    ):
        raise RuntimeError(
            f"Atom-ID mismatch in {item['path']}"
        )


# ------------------------------------------------------------
# Common centered reference
# ------------------------------------------------------------

reference_xyz, reference_centroid = (
    center_coordinates(
        snapshots[0]["xyz"]
    )
)


rows = []
alignment_records = []

for item in snapshots:

    centered, centroid = center_coordinates(
        item["xyz"]
    )

    aligned, rotation = kabsch_align(
        centered,
        reference_xyz
    )

    diff = aligned-reference_xyz

    rmsd = float(
        np.sqrt(
            np.mean(
                np.sum(diff*diff,axis=1)
            )
        )
    )

    vector = aligned.reshape(-1)

    row = {
        "temperature_K":
            item["temperature_K"],
        "snapshot":
            item["snapshot"],
    }

    for i,value in enumerate(vector):

        row[f"f{i:03d}"] = float(value)

    rows.append(row)

    alignment_records.append(
        {
            "temperature_K":
                item["temperature_K"],
            "snapshot":
                item["snapshot"],
            "source_file":
                str(item["path"].resolve()),
            "centroid_before_alignment_A":
                centroid.tolist(),
            "rmsd_to_reference_A":
                rmsd,
            "rotation_determinant":
                float(np.linalg.det(rotation)),
        }
    )


df = pd.DataFrame(rows)

descriptor_csv = (
    OUT
    / "snapshot_descriptor_matrix.csv"
)

df.to_csv(
    descriptor_csv,
    index=False
)


audit = {
    "status":
        "COORDINATE_PREPROCESSING_VALIDATED",
    "coordinate_source":
        "LAMMPS xu yu zu",
    "periodic_boundary_handling":
        "UNWRAPPED_AT_DUMP_TIME",
    "translation_handling":
        "CENTROID_REMOVAL",
    "rotation_handling":
        "KABSCH_ALIGNMENT",
    "alignment_reference":
        {
            "temperature_K":
                snapshots[0]["temperature_K"],
            "snapshot":
                snapshots[0]["snapshot"],
            "source_file":
                str(
                    snapshots[0]["path"].resolve()
                ),
        },
    "atom_count":
        int(len(reference_ids)),
    "snapshot_count":
        int(len(snapshots)),
    "reflection_allowed":
        False,
    "alignment_records":
        alignment_records,
}

(
    OUT
    / "coordinate_preprocessing_audit.json"
).write_text(
    json.dumps(
        audit,
        indent=2
    )
    + "\n"
)


print("="*90)
print("PHASE5-D22")
print("UNWRAPPED + CENTERED + KABSCH-ALIGNED DESCRIPTORS")
print("="*90)

print("Snapshots :",len(snapshots))
print("Atoms     :",len(reference_ids))
print("Features  :",3*len(reference_ids))
print()
print(descriptor_csv)
print()
print("COORDINATE_PREPROCESSING_VALIDATED")
