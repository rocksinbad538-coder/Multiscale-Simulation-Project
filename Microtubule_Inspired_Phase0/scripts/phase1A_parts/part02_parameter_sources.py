#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import tarfile
import zipfile

ROOT = pathlib.Path(__file__).resolve().parents[2]

DOWNLOADS = pathlib.Path.home() / "Downloads"

candidate_files = []

for extension in (
    "*.zip",
    "*.tar",
    "*.tar.gz",
    "*.tgz",
    "*.pdf",
):

    candidate_files.extend(
        DOWNLOADS.glob(extension)
    )

print()
print("[4] DISCOVER PRIMARY LITERATURE FILES")
print()

for f in sorted(candidate_files):

    print(f.name)

print()

supplement_archive = None

for f in candidate_files:

    name = f.name.lower()

    if "supp" in name:
        supplement_archive = f

    elif "support" in name:
        supplement_archive = f

if supplement_archive is None:

    print("supplement_archive = NOT FOUND")

else:

    print(
        "supplement_archive =",
        supplement_archive,
    )

print()

hbn_archive = None

for f in candidate_files:

    if "hbn" in f.name.lower():

        hbn_archive = f

if hbn_archive is None:

    print("hbn_archive = NOT FOUND")

else:

    print(
        "hbn_archive =",
        hbn_archive,
    )

print()

archive_members = []

if supplement_archive is not None:

    if zipfile.is_zipfile(supplement_archive):

        with zipfile.ZipFile(supplement_archive) as z:

            archive_members = z.namelist()

elif hbn_archive is not None:

    with tarfile.open(hbn_archive) as tar:

        archive_members = tar.getnames()

print("[5] ARCHIVE CONTENTS")

for member in archive_members:

    print(member)
