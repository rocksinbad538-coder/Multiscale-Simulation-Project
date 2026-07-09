#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/phase1A/day024_chemical_end_rim_design"
SCRIPTS = ROOT / "scripts/phase1A"
G3I = BASE / "11_r2_alternating_bn_trimer_bridge_graph"
G3K = BASE / "13_r2_trimer_bridge_conformer_and_h_refinement"
G3K1 = BASE / "14_r2_trimer_bridge_search_completeness_audit"
OUT = BASE / "15_r2_full_density_longer_bn_bridge_screen"

HELPER_GEOMETRY = SCRIPTS / "audit_day024_r2_trimer_bridge_search_completeness.py"
HELPER_GRAPH = SCRIPTS / "build_and_validate_day024_r2_alternating_bn_trimer_bridge_graph.py"
GRAPH_NODES = G3I / "r2_alternating_bn_trimer_bridge_graph_nodes.csv"
GRAPH_EDGES = G3I / "r2_alternating_bn_trimer_bridge_graph_edges.csv"
GRAPH_SUMMARY = G3I / "r2_alternating_bn_trimer_bridge_graph_summary.csv"
FIXED_COORDINATES = G3K / "r2_trimer_bridge_refined_coordinates.csv"
REDESIGN_SUMMARY = G3K1 / "r2_trimer_search_completeness_summary.csv"

LIBRARY_CSV = OUT / "r2_full_density_longer_bridge_library_summary.csv"
PAIR_CSV = OUT / "r2_full_density_longer_bridge_pair_screen.csv"
MAPPING_CSV = OUT / "r2_full_density_longer_bridge_mapping_screen.csv"
CLASS_CSV = OUT / "r2_full_density_longer_bridge_class_summary.csv"
SELECTED_CSV = OUT / "r2_full_density_longer_bridge_selected_candidate.csv"
SUMMARY_CSV = OUT / "r2_full_density_longer_bridge_screen_summary.csv"
GATES_CSV = OUT / "r2_full_density_longer_bridge_screen_gates.csv"
JSON_OUT = OUT / "r2_full_density_longer_bridge_screen.json"
MANIFEST_CSV = OUT / "r2_full_density_longer_bridge_source_manifest.csv"
REPORT_MD = OUT / "R2_FULL_DENSITY_LONGER_BN_BRIDGE_SCREEN_DAY024.md"

EXPECTED_GRAPH_DECISION = "R2_ALTERNATING_BN_TRIMER_BRIDGE_GRAPH_VALIDATED"
EXPECTED_REDESIGN_DECISION = (
    "R2_TRIMER_BRIDGE_TOPOLOGY_REDESIGN_CONFIRMED_"
    "BY_EXPANDED_CONFORMER_SEARCH"
)
PASS_DECISION = "R2_FULL_DENSITY_LONGER_BN_BRIDGE_CLASS_IDENTIFIED"
SPARSE_DECISION = "R2_FULL_DENSITY_LONGER_BN_BRIDGES_REQUIRE_SPARSE_ATTACHMENT_SCREEN"

BN = 0.144973
BRIDGE_CLASSES = (4, 5, 6)
ATTACHMENTS_PER_END = 15
LIBRARY_SAMPLES = 40000
NEAREST_CONFORMERS = 100
AZIMUTHS_DEG = tuple(float(value) for value in range(0, 360, 15))
PAIR_POOL_SIZE = 12
PAIR_PASS_TARGET = 36
GLOBAL_SWEEPS = 8
LOCAL_RADIUS_NM = 0.95
ANGLE_MIN_LIBRARY_DEG = 105.0
ANGLE_MAX_LIBRARY_DEG = 135.0
ANGLE_MIN_GATE_DEG = 70.0
ANGLE_MAX_GATE_DEG = 175.0
MAX_BOND_DEVIATION_NM = 0.003
MIN_HEAVY_CLEARANCE_NM = 0.120
SOFT_HEAVY_CLEARANCE_NM = 0.140


def load_module(path: Path, name: str):
    if not path.is_file():
        raise RuntimeError(f"Missing helper module: {path}")
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load helper module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty required file: {path}")


def read_rows(path: Path) -> list[dict[str, str]]:
    require_file(path)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError(f"No rows found in {path}")
    return rows


def read_one(path: Path) -> dict[str, str]:
    rows = read_rows(path)
    if len(rows) != 1:
        raise RuntimeError(f"Expected one row in {path}; found {len(rows)}")
    return rows[0]


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"No rows available for {path}")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def parse_float(row: dict[str, str], key: str) -> float:
    value = float(row[key])
    if not math.isfinite(value):
        raise RuntimeError(f"Non-finite numeric field {key!r}")
    return value


def parse_int(row: dict[str, str], key: str) -> int:
    return int(float(row[key]))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def normalized_batch(vectors: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(vectors, axis=1)
    if np.any(norms <= 1.0e-12):
        raise RuntimeError("Zero vector in batch normalization.")
    return vectors / norms[:, None]


def place_next_batch(
    first: np.ndarray,
    second: np.ndarray,
    third: np.ndarray,
    angle_rad: np.ndarray,
    torsion_rad: np.ndarray,
) -> np.ndarray:
    previous = normalized_batch(third - second)
    normal = np.cross(second - first, previous)
    bad = np.linalg.norm(normal, axis=1) <= 1.0e-12
    if np.any(bad):
        reference = np.tile(np.asarray([0.0, 0.0, 1.0]), (int(np.sum(bad)), 1))
        previous_bad = previous[bad]
        parallel = np.abs(np.sum(previous_bad * reference, axis=1)) > 0.90
        reference[parallel] = np.asarray([0.0, 1.0, 0.0])
        normal[bad] = np.cross(previous_bad, reference)
    normal = normalized_batch(normal)
    in_plane = normalized_batch(np.cross(normal, previous))
    direction = (
        -np.cos(angle_rad)[:, None] * previous
        + np.sin(angle_rad)[:, None]
        * (
            np.cos(torsion_rad)[:, None] * in_plane
            + np.sin(torsion_rad)[:, None] * normal
        )
    )
    return third + BN * direction


def build_library(bridge_atoms: int) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(35000 + bridge_atoms)
    angles = np.deg2rad(
        rng.uniform(
            ANGLE_MIN_LIBRARY_DEG,
            ANGLE_MAX_LIBRARY_DEG,
            size=(LIBRARY_SAMPLES, bridge_atoms),
        )
    )
    torsions = rng.uniform(
        0.0,
        2.0 * math.pi,
        size=(LIBRARY_SAMPLES, bridge_atoms - 1),
    )
    p0 = np.zeros((LIBRARY_SAMPLES, 3), dtype=float)
    p1 = np.tile(np.asarray([BN, 0.0, 0.0]), (LIBRARY_SAMPLES, 1))
    p2 = p1 + BN * np.column_stack(
        (-np.cos(angles[:, 0]), np.sin(angles[:, 0]), np.zeros(LIBRARY_SAMPLES))
    )
    points = [p0, p1, p2]
    for index in range(1, bridge_atoms):
        points.append(
            place_next_batch(
                points[-3],
                points[-2],
                points[-1],
                angles[:, index],
                torsions[:, index - 1],
            )
        )
    chains = np.stack(points, axis=1)
    distances = np.linalg.norm(chains[:, -1] - chains[:, 0], axis=1)
    return {"chains": chains, "distances": distances}




def map_chain_generic(
    chain: np.ndarray,
    start: np.ndarray,
    finish: np.ndarray,
    azimuth_deg: float,
    mirror: bool,
    geometry,
) -> np.ndarray:
    working = np.array(chain, dtype=float, copy=True)
    if mirror:
        working[:, 2] *= -1.0
    target = finish - start
    alignment = geometry.rotation_from_vectors(working[-1] - working[0], target)
    mapped = (working - working[0]) @ alignment.T
    mapped = mapped @ geometry.rotation_about_axis(
        geometry.normalized(target),
        math.radians(azimuth_deg),
    ).T
    correction = target - mapped[-1]
    denominator = mapped.shape[0] - 1
    for index in range(1, mapped.shape[0]):
        mapped[index] += (index / denominator) * correction
    mapped += start
    return mapped

def opposite(element: str) -> str:
    if element == "B":
        return "N"
    if element == "N":
        return "B"
    raise RuntimeError(f"Unexpected BN element: {element}")


def required_annulus_element(seed_element: str, bridge_atoms: int) -> str:
    return seed_element if (bridge_atoms + 1) % 2 == 0 else opposite(seed_element)


def bridge_sequence(seed_element: str, bridge_atoms: int) -> list[str]:
    values = []
    current = seed_element
    for _ in range(bridge_atoms):
        current = opposite(current)
        values.append(current)
    return values


def soft_penalty(distances: np.ndarray) -> float:
    return float(
        np.sum(np.maximum(SOFT_HEAVY_CLEARANCE_NM - distances, 0.0) ** 2)
    )


def evaluate_chain(
    chain: np.ndarray,
    start_neighbors: list[np.ndarray],
    finish_neighbors: list[np.ndarray],
    local_fixed: np.ndarray,
    geometry,
) -> dict[str, Any]:
    bonds = np.linalg.norm(np.diff(chain, axis=0), axis=1)
    max_bond_deviation = float(np.max(np.abs(bonds - BN)))
    angles = [
        geometry.angle_degrees(chain[index - 1], chain[index], chain[index + 1])
        for index in range(1, chain.shape[0] - 1)
    ]
    angles.extend(
        geometry.angle_degrees(neighbor, chain[0], chain[1])
        for neighbor in start_neighbors
    )
    angles.extend(
        geometry.angle_degrees(neighbor, chain[-1], chain[-2])
        for neighbor in finish_neighbors
    )
    angle_values = np.asarray(angles, dtype=float)
    angle_violations = int(
        np.sum(
            (angle_values < ANGLE_MIN_GATE_DEG)
            | (angle_values > ANGLE_MAX_GATE_DEG)
        )
    )
    self_distances = []
    for first in range(chain.shape[0]):
        for second in range(first + 2, chain.shape[0]):
            self_distances.append(float(np.linalg.norm(chain[first] - chain[second])))
    internal = chain[1:-1]
    fixed_distances = (
        np.linalg.norm(
            internal[:, None, :] - local_fixed[None, :, :],
            axis=2,
        ).reshape(-1)
        if local_fixed.size
        else np.asarray([], dtype=float)
    )
    nonbonded = np.concatenate((np.asarray(self_distances), fixed_distances))
    clashes = int(np.sum(nonbonded < MIN_HEAVY_CLEARANCE_NM))
    minimum_clearance = float(np.min(nonbonded))
    passes = (
        angle_violations == 0
        and clashes == 0
        and max_bond_deviation <= MAX_BOND_DEVIATION_NM
    )
    return {
        "passes": passes,
        "angle_violations": angle_violations,
        "minimum_angle_deg": float(np.min(angle_values)),
        "maximum_angle_deg": float(np.max(angle_values)),
        "maximum_bond_deviation_nm": max_bond_deviation,
        "clash_count": clashes,
        "minimum_clearance_nm": minimum_clearance,
        "clearance_penalty": soft_penalty(nonbonded),
    }


def local_score(candidate: dict[str, Any]) -> tuple[Any, ...]:
    return (
        int(candidate["clash_count"]),
        int(candidate["angle_violations"]),
        float(candidate["clearance_penalty"]),
        -float(candidate["minimum_angle_deg"]),
        -float(candidate["minimum_clearance_nm"]),
        float(candidate["maximum_bond_deviation_nm"]),
        float(candidate["library_distance_error_nm"]),
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for required in (
        HELPER_GEOMETRY,
        HELPER_GRAPH,
        GRAPH_NODES,
        GRAPH_EDGES,
        GRAPH_SUMMARY,
        FIXED_COORDINATES,
        REDESIGN_SUMMARY,
    ):
        require_file(required)

    geometry = load_module(HELPER_GEOMETRY, "day024_geometry_helper")
    graph_tools = load_module(HELPER_GRAPH, "day024_graph_helper")

    graph_summary = read_one(GRAPH_SUMMARY)
    redesign_summary = read_one(REDESIGN_SUMMARY)
    node_rows = read_rows(GRAPH_NODES)
    edge_rows = read_rows(GRAPH_EDGES)
    coordinate_rows = read_rows(FIXED_COORDINATES)

    if graph_summary.get("decision") != EXPECTED_GRAPH_DECISION:
        raise RuntimeError("Gate 3I graph is not accepted.")
    if redesign_summary.get("decision") != EXPECTED_REDESIGN_DECISION:
        raise RuntimeError("Gate 3K.1 does not contain the expected redesign decision.")

    nodes = {row["node_id"]: row for row in node_rows}
    coordinates = {
        row["node_id"]: np.asarray(
            [parse_float(row, "x_nm"), parse_float(row, "y_nm"), parse_float(row, "z_nm")],
            dtype=float,
        )
        for row in coordinate_rows
    }
    if set(nodes) != set(coordinates):
        raise RuntimeError("Gate 3I nodes and Gate 3K coordinates disagree.")

    current_adjacency = {node_id: set() for node_id in nodes}
    for row in edge_rows:
        current_adjacency[row["source_node"]].add(row["target_node"])
        current_adjacency[row["target_node"]].add(row["source_node"])

    old_bridge_ids = {
        node_id
        for node_id, row in nodes.items()
        if row["node_type"] == "ALTERNATING_BN_TRIMER_BRIDGE"
    }
    fixed_heavy_ids = [
        node_id
        for node_id, row in nodes.items()
        if row["element"] != "H" and node_id not in old_bridge_ids
    ]
    fixed_heavy_set = set(fixed_heavy_ids)
    base_adjacency = {node_id: set() for node_id in fixed_heavy_ids}
    for row in edge_rows:
        first = row["source_node"]
        second = row["target_node"]
        if first in fixed_heavy_set and second in fixed_heavy_set:
            base_adjacency[first].add(second)
            base_adjacency[second].add(first)

    libraries: dict[int, dict[str, np.ndarray]] = {}
    library_rows = []
    print("Generating deterministic full-density longer-bridge libraries...")
    for bridge_atoms in BRIDGE_CLASSES:
        library = build_library(bridge_atoms)
        libraries[bridge_atoms] = library
        library_rows.append(
            {
                "bridge_atoms_per_attachment": bridge_atoms,
                "bonds_per_path": bridge_atoms + 1,
                "sample_count": int(library["chains"].shape[0]),
                "minimum_endpoint_distance_nm": float(np.min(library["distances"])),
                "median_endpoint_distance_nm": float(np.median(library["distances"])),
                "maximum_endpoint_distance_nm": float(np.max(library["distances"])),
                "contour_length_nm": (bridge_atoms + 1) * BN,
                "random_seed": 35000 + bridge_atoms,
            }
        )
        print(
            f"  m={bridge_atoms}: samples={LIBRARY_SAMPLES}; "
            f"span={np.min(library['distances']):.6f}-{np.max(library['distances']):.6f} nm"
        )
    write_rows(LIBRARY_CSV, library_rows)

    pair_cache: dict[tuple[str, int, str, str], dict[str, Any]] = {}
    pair_rows: list[dict[str, Any]] = []

    def screen_pair(end: str, bridge_atoms: int, seed_id: str, annulus_id: str) -> dict[str, Any]:
        key = (end, bridge_atoms, seed_id, annulus_id)
        if key in pair_cache:
            return pair_cache[key]
        start = coordinates[seed_id]
        finish = coordinates[annulus_id]
        start_neighbors = [coordinates[node_id] for node_id in base_adjacency[seed_id]]
        finish_neighbors = [coordinates[node_id] for node_id in base_adjacency[annulus_id]]
        if len(start_neighbors) != 2 or len(finish_neighbors) != 2:
            raise RuntimeError(f"{key}: endpoint base-heavy neighbors are not 2/2.")
        midpoint = 0.5 * (start + finish)
        local_ids = [
            node_id
            for node_id in fixed_heavy_ids
            if node_id not in {seed_id, annulus_id}
            and float(np.linalg.norm(coordinates[node_id] - midpoint)) <= LOCAL_RADIUS_NM
        ]
        local_fixed = (
            np.asarray([coordinates[node_id] for node_id in local_ids], dtype=float)
            if local_ids
            else np.empty((0, 3), dtype=float)
        )
        library = libraries[bridge_atoms]
        target_distance = float(np.linalg.norm(finish - start))
        errors = np.abs(library["distances"] - target_distance)
        count = min(NEAREST_CONFORMERS, errors.size)
        indices = np.argpartition(errors, count - 1)[:count]
        indices = indices[np.argsort(errors[indices])]
        feasible: list[dict[str, Any]] = []
        best_failed: dict[str, Any] | None = None
        tested = 0
        for library_index in indices:
            base_chain = library["chains"][int(library_index)]
            for mirror in (False, True):
                for azimuth in AZIMUTHS_DEG:
                    tested += 1
                    mapped = map_chain_generic(
                        base_chain,
                        start,
                        finish,
                        azimuth,
                        mirror,
                        geometry,
                    )
                    metrics = evaluate_chain(
                        mapped,
                        start_neighbors,
                        finish_neighbors,
                        local_fixed,
                        geometry,
                    )
                    candidate = {
                        "chain": mapped,
                        "library_index": int(library_index),
                        "library_distance_error_nm": float(errors[int(library_index)]),
                        "mirror": mirror,
                        "azimuth_deg": azimuth,
                        **metrics,
                    }
                    if metrics["passes"]:
                        feasible.append(candidate)
                        if len(feasible) >= PAIR_PASS_TARGET:
                            break
                    elif best_failed is None or local_score(candidate) < local_score(best_failed):
                        best_failed = candidate
                if len(feasible) >= PAIR_PASS_TARGET:
                    break
            if len(feasible) >= PAIR_PASS_TARGET:
                break
        feasible.sort(key=local_score)
        feasible = feasible[:PAIR_POOL_SIZE]
        best = feasible[0] if feasible else best_failed
        if best is None:
            raise RuntimeError(f"No candidates tested for {key}")
        result = {
            "end": end,
            "bridge_atoms": bridge_atoms,
            "seed_node": seed_id,
            "annulus_node": annulus_id,
            "target_distance_nm": target_distance,
            "local_fixed_heavy_atoms": len(local_ids),
            "candidates_tested": tested,
            "pair_feasible": bool(feasible),
            "retained_feasible_candidates": len(feasible),
            "best_library_distance_error_nm": float(best["library_distance_error_nm"]),
            "best_minimum_angle_deg": float(best["minimum_angle_deg"]),
            "best_maximum_angle_deg": float(best["maximum_angle_deg"]),
            "best_maximum_bond_deviation_nm": float(best["maximum_bond_deviation_nm"]),
            "best_minimum_clearance_nm": float(best["minimum_clearance_nm"]),
            "best_angle_violations": int(best["angle_violations"]),
            "best_clash_count": int(best["clash_count"]),
            "pool": feasible,
        }
        pair_cache[key] = result
        pair_rows.append({key: value for key, value in result.items() if key != "pool"})
        return result

    def refine_mapping(pair_results: list[dict[str, Any]]) -> dict[str, Any]:
        if any(not result["pair_feasible"] for result in pair_results):
            return {
                "feasible": False,
                "interbridge_clashes": -1,
                "minimum_interbridge_clearance_nm": 0.0,
                "minimum_angle_deg": min(float(result["best_minimum_angle_deg"]) for result in pair_results),
                "maximum_bond_deviation_nm": max(float(result["best_maximum_bond_deviation_nm"]) for result in pair_results),
                "minimum_local_clearance_nm": min(float(result["best_minimum_clearance_nm"]) for result in pair_results),
                "indices": [],
            }
        indices = [0] * len(pair_results)
        for sweep in range(GLOBAL_SWEEPS):
            changed = 0
            order = range(len(pair_results)) if sweep % 2 == 0 else reversed(range(len(pair_results)))
            for path_index in order:
                other_points = []
                for other_index, result in enumerate(pair_results):
                    if other_index != path_index:
                        other_points.extend(result["pool"][indices[other_index]]["chain"][1:-1])
                other_array = np.asarray(other_points, dtype=float)
                scores = []
                for candidate in pair_results[path_index]["pool"]:
                    internal = candidate["chain"][1:-1]
                    distances = np.linalg.norm(
                        internal[:, None, :] - other_array[None, :, :],
                        axis=2,
                    ).reshape(-1)
                    scores.append(
                        (
                            int(np.sum(distances < MIN_HEAVY_CLEARANCE_NM)),
                            soft_penalty(distances),
                            -min(float(candidate["minimum_clearance_nm"]), float(np.min(distances))),
                            -float(candidate["minimum_angle_deg"]),
                            float(candidate["maximum_bond_deviation_nm"]),
                        )
                    )
                best_index = min(range(len(scores)), key=lambda index: scores[index])
                if best_index != indices[path_index]:
                    indices[path_index] = best_index
                    changed += 1
            if changed == 0:
                break
        inter_distances = []
        for first in range(len(pair_results)):
            first_points = pair_results[first]["pool"][indices[first]]["chain"][1:-1]
            for second in range(first + 1, len(pair_results)):
                second_points = pair_results[second]["pool"][indices[second]]["chain"][1:-1]
                inter_distances.extend(
                    np.linalg.norm(
                        first_points[:, None, :] - second_points[None, :, :],
                        axis=2,
                    ).reshape(-1)
                )
        inter_values = np.asarray(inter_distances, dtype=float)
        clashes = int(np.sum(inter_values < MIN_HEAVY_CLEARANCE_NM))
        selected = [result["pool"][indices[index]] for index, result in enumerate(pair_results)]
        return {
            "feasible": clashes == 0,
            "interbridge_clashes": clashes,
            "minimum_interbridge_clearance_nm": float(np.min(inter_values)),
            "minimum_angle_deg": min(float(candidate["minimum_angle_deg"]) for candidate in selected),
            "maximum_bond_deviation_nm": max(float(candidate["maximum_bond_deviation_nm"]) for candidate in selected),
            "minimum_local_clearance_nm": min(float(candidate["minimum_clearance_nm"]) for candidate in selected),
            "indices": indices,
        }

    mapping_rows: list[dict[str, Any]] = []
    mapping_objects: dict[str, dict[str, Any]] = {}
    class_rows: list[dict[str, Any]] = []
    best_by_class_end: dict[tuple[int, str], dict[str, Any]] = {}

    print("Screening 15-attachment longer-bridge mappings...")
    for bridge_atoms in BRIDGE_CLASSES:
        feasible_counts: dict[str, int] = {}
        for end in ("LOWER", "UPPER"):
            seeds = sorted(
                [
                    row
                    for row in node_rows
                    if row["end"] == end
                    and row["node_type"] == "HEXAGONAL_EDGE_COMPLETION_SEED"
                ],
                key=lambda row: parse_int(row, "circumferential_index"),
            )
            seed_element = seeds[0]["element"]
            if len(seeds) != 30 or any(row["element"] != seed_element for row in seeds):
                raise RuntimeError(f"{end}: invalid seed population.")
            annulus_element = required_annulus_element(seed_element, bridge_atoms)
            outer = sorted(
                [
                    row
                    for row in node_rows
                    if row["end"] == end
                    and row["node_type"] == "ANNULUS_OUTER_BOUNDARY"
                    and row["element"] == annulus_element
                ],
                key=lambda row: parse_float(row, "angle_turns"),
            )
            if len(outer) != 15:
                raise RuntimeError(f"{end}: invalid outer-annulus endpoint population.")
            feasible_mappings = []
            for parity in (0, 1):
                selected_seeds = [
                    row
                    for row in seeds
                    if parse_int(row, "circumferential_index") % 2 == parity
                ]
                for orientation in (1, -1):
                    for rotation in range(15):
                        mapped_outer = [
                            outer[(orientation * index + rotation) % 15]
                            for index in range(15)
                        ]
                        pair_results = []
                        for seed_row, outer_row in zip(selected_seeds, mapped_outer):
                            pair_results.append(
                                screen_pair(
                                    end,
                                    bridge_atoms,
                                    seed_row["node_id"],
                                    outer_row["node_id"],
                                )
                            )
                            if not pair_results[-1]["pair_feasible"]:
                                break
                        result = (
                            refine_mapping(pair_results)
                            if len(pair_results) == 15
                            else {
                                "feasible": False,
                                "interbridge_clashes": -1,
                                "minimum_interbridge_clearance_nm": 0.0,
                                "minimum_angle_deg": min(float(item["best_minimum_angle_deg"]) for item in pair_results),
                                "maximum_bond_deviation_nm": max(float(item["best_maximum_bond_deviation_nm"]) for item in pair_results),
                                "minimum_local_clearance_nm": min(float(item["best_minimum_clearance_nm"]) for item in pair_results),
                                "indices": [],
                            }
                        )
                        mapping_id = f"M{bridge_atoms}:{end}:P{parity}:O{orientation}:R{rotation}"
                        row = {
                            "mapping_id": mapping_id,
                            "bridge_atoms_per_attachment": bridge_atoms,
                            "attachments_per_end": 15,
                            "end": end,
                            "seed_element": seed_element,
                            "bridge_element_sequence": "-".join(bridge_sequence(seed_element, bridge_atoms)),
                            "annulus_endpoint_element": annulus_element,
                            "seed_parity": parity,
                            "orientation": orientation,
                            "rotation": rotation,
                            "endpoint_pairs_evaluated": len(pair_results),
                            "all_endpoint_pairs_locally_feasible": len(pair_results) == 15 and all(item["pair_feasible"] for item in pair_results),
                            "mapping_geometrically_feasible": bool(result["feasible"]),
                            "interbridge_clashes": int(result["interbridge_clashes"]),
                            "minimum_interbridge_clearance_nm": float(result["minimum_interbridge_clearance_nm"]),
                            "minimum_angle_deg": float(result["minimum_angle_deg"]),
                            "maximum_bond_deviation_nm": float(result["maximum_bond_deviation_nm"]),
                            "minimum_local_clearance_nm": float(result["minimum_local_clearance_nm"]),
                            "seed_nodes": " | ".join(item["node_id"] for item in selected_seeds),
                            "annulus_nodes": " | ".join(item["node_id"] for item in mapped_outer),
                        }
                        mapping_rows.append(row)
                        mapping_object = {
                            "row": row,
                            "seeds": selected_seeds,
                            "outer": mapped_outer,
                            "pair_results": pair_results,
                            "indices": result["indices"],
                        }
                        mapping_objects[mapping_id] = mapping_object
                        if result["feasible"]:
                            feasible_mappings.append(mapping_object)
            feasible_counts[end] = len(feasible_mappings)
            if feasible_mappings:
                feasible_mappings.sort(
                    key=lambda item: (
                        -float(item["row"]["minimum_angle_deg"]),
                        -min(
                            float(item["row"]["minimum_local_clearance_nm"]),
                            float(item["row"]["minimum_interbridge_clearance_nm"]),
                        ),
                        float(item["row"]["maximum_bond_deviation_nm"]),
                    )
                )
                best_by_class_end[(bridge_atoms, end)] = feasible_mappings[0]
            print(f"  m={bridge_atoms}, {end}: feasible mappings={len(feasible_mappings)}/60")
        class_feasible = feasible_counts["LOWER"] > 0 and feasible_counts["UPPER"] > 0
        if class_feasible:
            lower = best_by_class_end[(bridge_atoms, "LOWER")]["row"]
            upper = best_by_class_end[(bridge_atoms, "UPPER")]["row"]
            minimum_angle = min(float(lower["minimum_angle_deg"]), float(upper["minimum_angle_deg"]))
            minimum_clearance = min(
                float(lower["minimum_local_clearance_nm"]),
                float(lower["minimum_interbridge_clearance_nm"]),
                float(upper["minimum_local_clearance_nm"]),
                float(upper["minimum_interbridge_clearance_nm"]),
            )
            maximum_bond_deviation = max(
                float(lower["maximum_bond_deviation_nm"]),
                float(upper["maximum_bond_deviation_nm"]),
            )
        else:
            minimum_angle = 0.0
            minimum_clearance = 0.0
            maximum_bond_deviation = math.inf
        class_rows.append(
            {
                "bridge_atoms_per_attachment": bridge_atoms,
                "bonds_per_path": bridge_atoms + 1,
                "attachments_per_end": 15,
                "lower_feasible_mapping_count": feasible_counts["LOWER"],
                "upper_feasible_mapping_count": feasible_counts["UPPER"],
                "uniform_class_feasible": class_feasible,
                "class_minimum_angle_deg": minimum_angle,
                "class_minimum_clearance_nm": minimum_clearance,
                "class_maximum_bond_deviation_nm": maximum_bond_deviation,
                "added_bridge_heavy_atoms_per_end": 15 * bridge_atoms,
                "total_added_heavy_atoms_per_end": 156 + 15 * bridge_atoms,
                "total_H_atoms_per_end": 72 + 15 * (bridge_atoms - 2),
                "candidate_is_final_chemistry": False,
            }
        )

    write_rows(PAIR_CSV, pair_rows)
    write_rows(MAPPING_CSV, mapping_rows)
    write_rows(CLASS_CSV, class_rows)

    feasible_classes = [row for row in class_rows if row["uniform_class_feasible"]]
    selected_class = None
    selected_rows: list[dict[str, Any]] = []
    topology_metrics: dict[str, Any] = {}
    if feasible_classes:
        feasible_classes.sort(
            key=lambda row: (
                int(row["bridge_atoms_per_attachment"]),
                -float(row["class_minimum_angle_deg"]),
                -float(row["class_minimum_clearance_nm"]),
            )
        )
        selected_class = feasible_classes[0]
        selected_m = int(selected_class["bridge_atoms_per_attachment"])
        selected_objects = [best_by_class_end[(selected_m, end)] for end in ("LOWER", "UPPER")]

        candidate_adjacency = {node_id: set(neighbors) for node_id, neighbors in base_adjacency.items()}
        candidate_elements = {node_id: nodes[node_id]["element"] for node_id in fixed_heavy_ids}
        candidate_paths = []
        for mapping in selected_objects:
            end = mapping["row"]["end"]
            selected_rows.append({"classification": "SELECTED_END_MAPPING", **mapping["row"]})
            for path_index, (seed_row, outer_row) in enumerate(zip(mapping["seeds"], mapping["outer"])):
                sequence = bridge_sequence(seed_row["element"], selected_m)
                bridge_ids = []
                for position, element in enumerate(sequence, start=1):
                    node_id = f"CAND:{end}:{path_index:02d}:{position}"
                    candidate_adjacency[node_id] = set()
                    candidate_elements[node_id] = element
                    bridge_ids.append(node_id)
                path_nodes = [seed_row["node_id"], *bridge_ids, outer_row["node_id"]]
                path_edges = set()
                for index in range(len(path_nodes) - 1):
                    first = path_nodes[index]
                    second = path_nodes[index + 1]
                    candidate_adjacency[first].add(second)
                    candidate_adjacency[second].add(first)
                    path_edges.add(tuple(sorted((first, second))))
                candidate_paths.append(
                    {
                        "seed": seed_row["node_id"],
                        "outer": outer_row["node_id"],
                        "edges": path_edges,
                        "edge_count": selected_m + 1,
                    }
                )

        nonheteropolar = 0
        for first, neighbors in candidate_adjacency.items():
            for second in neighbors:
                if first < second and {candidate_elements[first], candidate_elements[second]} != {"B", "N"}:
                    nonheteropolar += 1
        components = graph_tools.connected_components(candidate_adjacency)
        is_bipartite, _ = graph_tools.bipartite_coloring(candidate_adjacency)
        four_cycles = graph_tools.count_four_cycles(candidate_adjacency)
        paths_without_cycle = 0
        cycle_lengths = []
        for path in candidate_paths:
            alternative = graph_tools.shortest_path_length(
                candidate_adjacency,
                path["seed"],
                path["outer"],
                path["edges"],
            )
            if alternative is None:
                paths_without_cycle += 1
            else:
                cycle_lengths.append(alternative + path["edge_count"])
        topology_metrics = {
            "selected_heavy_nodes": len(candidate_adjacency),
            "selected_components": len(components),
            "selected_bipartite": is_bipartite,
            "selected_four_cycles": four_cycles,
            "selected_nonheteropolar_edges": nonheteropolar,
            "selected_degree_over3": sum(len(neighbors) > 3 for neighbors in candidate_adjacency.values()),
            "selected_degree_below2": sum(len(neighbors) < 2 for neighbors in candidate_adjacency.values()),
            "selected_H_required": sum(max(0, 3 - len(neighbors)) for neighbors in candidate_adjacency.values()),
            "selected_paths_without_cycle": paths_without_cycle,
            "selected_cycle_minimum": min(cycle_lengths) if cycle_lengths else 0,
            "selected_cycle_maximum": max(cycle_lengths) if cycle_lengths else 0,
        }
        selected_rows.append({"classification": "SELECTED_CLASS", **selected_class, **topology_metrics})
        write_rows(SELECTED_CSV, selected_rows)

    audit_gates = {
        "Gate3I_graph_is_accepted": graph_summary.get("decision") == EXPECTED_GRAPH_DECISION,
        "Gate3K1_trimer_redesign_is_confirmed": redesign_summary.get("decision") == EXPECTED_REDESIGN_DECISION,
        "bridge_classes_4_5_6_were_screened": len(class_rows) == 3,
        "60_mappings_per_end_and_class_were_screened": len(mapping_rows) == 3 * 2 * 60,
        "all_pair_metrics_are_finite": all(
            all(
                math.isfinite(float(row[field]))
                for field in (
                    "target_distance_nm",
                    "best_library_distance_error_nm",
                    "best_minimum_angle_deg",
                    "best_maximum_angle_deg",
                    "best_maximum_bond_deviation_nm",
                    "best_minimum_clearance_nm",
                )
            )
            for row in pair_rows
        ),
    }
    failed_audit_gates = [name for name, passed in audit_gates.items() if not passed]
    integrity_pass = not failed_audit_gates
    topology_pass = False
    if selected_class is not None:
        topology_pass = (
            topology_metrics["selected_components"] == 1
            and topology_metrics["selected_bipartite"]
            and topology_metrics["selected_four_cycles"] == 0
            and topology_metrics["selected_nonheteropolar_edges"] == 0
            and topology_metrics["selected_degree_over3"] == 0
            and topology_metrics["selected_degree_below2"] == 0
            and topology_metrics["selected_paths_without_cycle"] == 0
            and topology_metrics["selected_cycle_minimum"] >= 6
        )
    accepted = integrity_pass and selected_class is not None and topology_pass
    decision = PASS_DECISION if accepted else SPARSE_DECISION
    next_step = (
        "BUILD_AND_VALIDATE_R2_SELECTED_FULL_DENSITY_LONGER_BN_BRIDGE_GRAPH"
        if accepted
        else "SCREEN_R2_SPARSE_LONGER_BN_BRIDGE_TOPOLOGIES"
    )

    summary = {
        "decision": decision,
        "bridge_classes_screened": len(class_rows),
        "endpoint_pairs_screened": len(pair_rows),
        "mappings_screened": len(mapping_rows),
        "feasible_class_count": len(feasible_classes),
        "selected_bridge_atoms_per_attachment": "" if selected_class is None else selected_class["bridge_atoms_per_attachment"],
        "selected_attachments_per_end": "" if selected_class is None else 15,
        "selected_minimum_angle_deg": "" if selected_class is None else selected_class["class_minimum_angle_deg"],
        "selected_minimum_clearance_nm": "" if selected_class is None else selected_class["class_minimum_clearance_nm"],
        "selected_maximum_bond_deviation_nm": "" if selected_class is None else selected_class["class_maximum_bond_deviation_nm"],
        **topology_metrics,
        "audit_integrity_pass": integrity_pass,
        "selected_topology_graph_checks_pass": topology_pass,
        "candidate_is_final_chemistry": False,
        "selected_graph_generation_authorized": accepted,
        "coordinate_generation_authorized": False,
        "molecular_topology_generation_authorized": False,
        "formal_charge_assignment_authorized": False,
        "force_field_parameterization_authorized": False,
        "energy_minimization_authorized": False,
        "MD_authorized": False,
        "QM_authorized": False,
        "failed_audit_gates": " | ".join(failed_audit_gates),
        "required_next_step": next_step,
    }
    write_rows(SUMMARY_CSV, [summary])
    write_rows(GATES_CSV, [{"gate": name, "pass": passed} for name, passed in audit_gates.items()])
    JSON_OUT.write_text(
        json.dumps(
            {
                "summary": summary,
                "class_summaries": class_rows,
                "selected_rows": selected_rows,
                "audit_gates": audit_gates,
                "limitations": [
                    "This is a deterministic geometric screen, not an energy calculation.",
                    "Diagnostic bridge conformers are not applied to the accepted structure.",
                    "A selected class must still be rebuilt and validated as a complete graph.",
                    "No molecular topology, formal charges, force-field parameters, minimization, MD, or QM calculation was generated.",
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    write_rows(
        MANIFEST_CSV,
        [
            {"role": "geometry_helper", "file": relative(HELPER_GEOMETRY), "sha256": sha256(HELPER_GEOMETRY)},
            {"role": "graph_helper", "file": relative(HELPER_GRAPH), "sha256": sha256(HELPER_GRAPH)},
            {"role": "Gate3I_graph_nodes", "file": relative(GRAPH_NODES), "sha256": sha256(GRAPH_NODES)},
            {"role": "Gate3I_graph_edges", "file": relative(GRAPH_EDGES), "sha256": sha256(GRAPH_EDGES)},
            {"role": "Gate3I_graph_summary", "file": relative(GRAPH_SUMMARY), "sha256": sha256(GRAPH_SUMMARY)},
            {"role": "Gate3K_fixed_coordinates", "file": relative(FIXED_COORDINATES), "sha256": sha256(FIXED_COORDINATES)},
            {"role": "Gate3K1_redesign_summary", "file": relative(REDESIGN_SUMMARY), "sha256": sha256(REDESIGN_SUMMARY)},
        ],
    )

    class_lines = "\n".join(
        f"- m={row['bridge_atoms_per_attachment']}: lower/upper feasible mappings="
        f"{row['lower_feasible_mapping_count']}/{row['upper_feasible_mapping_count']}; "
        f"feasible={row['uniform_class_feasible']}"
        for row in class_rows
    )
    gate_lines = "\n".join(
        f"- `{name}`: **{'PASS' if passed else 'FAIL'}**"
        for name, passed in audit_gates.items()
    )
    REPORT_MD.write_text(
        f"""# R2 Full-Density Longer BN Bridge Screen

## Scope

This gate screens alternating BN bridges containing four, five or six
heavy atoms while retaining 15 attachments per end.

No accepted coordinates were replaced. No molecular topology, formal
charges, force-field parameters, minimization, MD or QM calculation was
generated.

## Results

{class_lines}

## Selection

- Selected bridge atoms per attachment: **{summary['selected_bridge_atoms_per_attachment'] or 'NONE'}**
- Selected attachments per end: **{summary['selected_attachments_per_end'] or 'NONE'}**
- Selected minimum angle: **{summary['selected_minimum_angle_deg'] or 'N/A'}**
- Selected minimum clearance: **{summary['selected_minimum_clearance_nm'] or 'N/A'}**
- Selected graph checks pass: **{topology_pass}**

## Audit gates

{gate_lines}

## Decision

- Decision: **{decision}**
- Failed audit-integrity gates: **{'NONE' if not failed_audit_gates else ' | '.join(failed_audit_gates)}**
- Selected graph generation authorized: **{'YES' if accepted else 'NO'}**
- Coordinate generation authorized: **NO**
- Molecular topology generation authorized: **NO**
- Energy minimization authorized: **NO**
- MD authorized: **NO**
- QM authorized: **NO**
- Required next step: `{next_step}`
""",
        encoding="utf-8",
    )

    print("Day024 R2 full-density longer BN bridge screen completed.")
    print(f"Bridge classes / endpoint pairs / mappings screened: {len(class_rows)}/{len(pair_rows)}/{len(mapping_rows)}")
    for row in class_rows:
        print(
            "Bridge atoms / lower feasible / upper feasible / class feasible / minimum angle / minimum clearance: "
            f"{row['bridge_atoms_per_attachment']}/"
            f"{row['lower_feasible_mapping_count']}/"
            f"{row['upper_feasible_mapping_count']}/"
            f"{row['uniform_class_feasible']}/"
            f"{row['class_minimum_angle_deg']}/"
            f"{row['class_minimum_clearance_nm']}"
        )
    print(
        "Selected bridge atoms / attachments per end / minimum angle / minimum clearance / maximum bond deviation: "
        f"{summary['selected_bridge_atoms_per_attachment'] or 'NONE'}/"
        f"{summary['selected_attachments_per_end'] or 'NONE'}/"
        f"{summary['selected_minimum_angle_deg'] or 'N/A'}/"
        f"{summary['selected_minimum_clearance_nm'] or 'N/A'}/"
        f"{summary['selected_maximum_bond_deviation_nm'] or 'N/A'}"
    )
    if topology_metrics:
        print(
            "Selected graph components/bipartite/four-cycles/nonheteropolar/degree>3/degree<2: "
            f"{topology_metrics['selected_components']}/"
            f"{topology_metrics['selected_bipartite']}/"
            f"{topology_metrics['selected_four_cycles']}/"
            f"{topology_metrics['selected_nonheteropolar_edges']}/"
            f"{topology_metrics['selected_degree_over3']}/"
            f"{topology_metrics['selected_degree_below2']}"
        )
        print(
            "Selected paths without cycle / cycle min-max / H required: "
            f"{topology_metrics['selected_paths_without_cycle']}/"
            f"{topology_metrics['selected_cycle_minimum']}-"
            f"{topology_metrics['selected_cycle_maximum']}/"
            f"{topology_metrics['selected_H_required']}"
        )
    print(f"Decision: {decision}")
    print("Failed audit-integrity gates: " + ("NONE" if not failed_audit_gates else " | ".join(failed_audit_gates)))
    print(f"Selected graph generation authorized: {'YES' if accepted else 'NO'}")
    print("Coordinate generation authorized: NO")
    print("Molecular topology generation authorized: NO")
    print("Formal charge assignment authorized: NO")
    print("Force-field parameterization authorized: NO")
    print("Energy minimization authorized: NO")
    print("MD authorized: NO")
    print("QM authorized: NO")
    print(f"Required next step: {next_step}")
    for path in (
        LIBRARY_CSV,
        PAIR_CSV,
        MAPPING_CSV,
        CLASS_CSV,
        SUMMARY_CSV,
        GATES_CSV,
        JSON_OUT,
        MANIFEST_CSV,
        REPORT_MD,
    ):
        print(f"Wrote: {relative(path)}")
    if selected_class is not None:
        print(f"Wrote: {relative(SELECTED_CSV)}")


if __name__ == "__main__":
    main()
