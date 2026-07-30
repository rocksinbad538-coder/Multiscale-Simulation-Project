#!/usr/bin/env python3
"""
resp_common.py

Shared utilities for the RESP pipeline.

Version: 0.1 (foundation)

Purpose
-------
Centralize all reusable validation logic for the RESP pipeline so that
all gates share identical provenance, hashing, JSON handling and report
generation.

Current utilities
-----------------
- sha256()
- load_json()
- save_json()
- require_file()
- require_authorization()
- require_decision()
- validate_execution_binding()
- utc_now()
"""

from __future__ import annotations

from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def require_file(path: Path) -> None:
    if (not path.exists()) or path.stat().st_size == 0:
        raise RuntimeError(f"Required file missing: {path}")


def load_json(path: Path) -> dict:
    require_file(path)
    return json.loads(path.read_text(encoding="utf-8"))


def save_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(data, indent=2) + "\n",
        encoding="utf-8",
    )


def sha256(path: Path) -> str:
    require_file(path)

    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def require_decision(report: dict, expected: str) -> None:
    decision = report.get("decision")

    if decision != expected:
        raise RuntimeError(
            f"Decision mismatch.\n"
            f"Expected: {expected}\n"
            f"Observed: {decision}"
        )


def require_authorization(
    report: dict,
    authorization_name: str,
) -> None:

    authorizations = report.get("authorizations", {})

    if authorizations.get(authorization_name) is not True:
        raise RuntimeError(
            f"Authorization '{authorization_name}' not granted."
        )


def validate_execution_binding(
    previous_report: dict,
    execution_directory: str,
) -> None:

    observed = previous_report.get(
        "source_execution_directory",
        previous_report.get("execution_directory"),
    )

    if observed != execution_directory:
        raise RuntimeError(
            "Execution binding failed.\n"
            f"Expected: {execution_directory}\n"
            f"Observed: {observed}"
        )


def build_report(
    *,
    decision: str,
    execution_directory: str,
    authorizations: dict,
    sha256_items: dict | None = None,
    summary: dict | None = None,
    metadata: dict | None = None,
) -> dict:

    report = {
        "generated_utc": utc_now(),
        "decision": decision,
        "execution_directory": execution_directory,
        "authorizations": authorizations,
    }

    if sha256_items:
        report["sha256"] = sha256_items

    if summary:
        report["summary"] = summary

    if metadata:
        report["metadata"] = metadata

    return report


def print_gate_banner(title: str) -> None:
    print("=" * 100)
    print(title)
    print("=" * 100)
