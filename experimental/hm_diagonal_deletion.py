"""
Hertel-Mehlhorn-style diagonal deletion from an existing triangulation.

This is a project-local convex pre-merge backend.  It consumes triangulation
pieces, builds the face-adjacency graph induced by shared diagonals, and then
removes a shared diagonal whenever the two incident current faces have a convex
simple union.

The first implementation is intentionally conservative: it reuses the existing
boundary-chain union and convexity predicates from star_merge.py.  The important
algorithmic difference from the previous greedy pre-merge is that this is
edge/deletion-driven rather than score-driven all-candidate greedy merging.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import Deque, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from geometry_utils import Polygon
from star_merge import (
    Region,
    EPS,
    _edge_key,
    _edges,
    are_adjacent,
    clean_polygon,
    is_convex_polygon,
    shared_boundary_length,
    union_if_simple,
)

Pair = Tuple[int, int]


@dataclass
class HMDiagonalDeletionStats:
    initial_faces: int = 0
    initial_adjacency_edges: int = 0
    diagonal_tests: int = 0
    convexity_tests: int = 0
    accepted_deletions: int = 0
    union_failures: int = 0
    rejected_nonconvex: int = 0
    local_adjacency_updates: int = 0
    final_faces: int = 0


def _pair(i: int, j: int) -> Pair:
    return (i, j) if i < j else (j, i)


def _region_from_convex_piece(poly: Polygon, source_count: int = 1) -> Region:
    poly = clean_polygon(poly)
    if not is_convex_polygon(poly):
        raise ValueError("Expected a convex piece for HM diagonal-deletion initialization/output.")
    # For a convex polygon, K(P)=P.
    return Region(polygon=poly, kernel=poly[:], source_count=source_count)


def _build_initial_adjacency(active: Dict[int, Region]) -> Dict[int, Set[int]]:
    """Build adjacency from exact shared edge keys, with robust fallback if needed."""
    neighbors: Dict[int, Set[int]] = {i: set() for i in active}
    edge_to_faces: Dict[Tuple[Tuple[float, float], Tuple[float, float]], List[int]] = {}

    for face_id, region in active.items():
        for a, b in _edges(region.polygon):
            edge_to_faces.setdefault(_edge_key(a, b), []).append(face_id)

    for face_ids in edge_to_faces.values():
        if len(face_ids) == 2:
            i, j = face_ids
            neighbors[i].add(j)
            neighbors[j].add(i)

    # Fallback for rare split-edge cases.  Triangulations should not need this,
    # but it keeps the backend robust if the source triangulator changes.
    ids = sorted(active)
    for pos, i in enumerate(ids):
        for j in ids[pos + 1 :]:
            if j in neighbors[i]:
                continue
            if shared_boundary_length(active[i].polygon, active[j].polygon) > EPS:
                neighbors[i].add(j)
                neighbors[j].add(i)

    return neighbors


def _remove_region(neighbors: Dict[int, Set[int]], region_id: int) -> None:
    for nb in list(neighbors.get(region_id, set())):
        if nb in neighbors:
            neighbors[nb].discard(region_id)
    neighbors.pop(region_id, None)


def _enqueue_pair(queue: Deque[Pair], queued: Set[Pair], i: int, j: int) -> None:
    p = _pair(i, j)
    if p not in queued:
        queue.append(p)
        queued.add(p)


def hm_diagonal_deletion_from_triangles(
    triangles: Sequence[Polygon],
    log_prefix: str = "HM diagonal deletion",
) -> Tuple[List[Region], List[str], HMDiagonalDeletionStats]:
    """
    Convert triangulation pieces into convex regions by diagonal deletion.

    Returns:
        final convex Region objects, log messages, and diagnostic stats.
    """
    active: Dict[int, Region] = {
        i: _region_from_convex_piece(tri, source_count=1)
        for i, tri in enumerate(triangles)
    }
    neighbors = _build_initial_adjacency(active)

    stats = HMDiagonalDeletionStats(
        initial_faces=len(active),
        initial_adjacency_edges=sum(len(v) for v in neighbors.values()) // 2,
    )
    log: List[str] = []

    queue: Deque[Pair] = deque()
    queued: Set[Pair] = set()
    for i in sorted(neighbors):
        for j in sorted(neighbors[i]):
            if i < j:
                _enqueue_pair(queue, queued, i, j)

    next_id = len(active)
    step = 0

    while queue:
        i, j = queue.popleft()
        queued.discard((i, j))

        if i not in active or j not in active:
            continue
        if j not in neighbors.get(i, set()):
            continue

        stats.diagonal_tests += 1
        union_poly = union_if_simple(active[i].polygon, active[j].polygon)
        if union_poly is None:
            stats.union_failures += 1
            continue

        stats.convexity_tests += 1
        if not is_convex_polygon(union_poly):
            stats.rejected_nonconvex += 1
            continue

        # Delete the shared diagonal and merge the two incident faces.
        merged = Region(
            polygon=clean_polygon(union_poly),
            kernel=clean_polygon(union_poly),
            source_count=active[i].source_count + active[j].source_count,
        )
        if not is_convex_polygon(merged.polygon):
            # Safety guard against numerical or boundary-reconstruction issues.
            stats.rejected_nonconvex += 1
            continue

        step += 1
        stats.accepted_deletions += 1
        old_count = len(active)
        log.append(
            f"{log_prefix} {step}: deleted shared diagonal between faces {i} and {j}; "
            f"pieces {active[i].source_count}+{active[j].source_count} -> {merged.source_count}; "
            f"face count {old_count}->{old_count - 1}."
        )

        candidate_neighbors = (neighbors.get(i, set()) | neighbors.get(j, set())) - {i, j}

        _remove_region(neighbors, i)
        _remove_region(neighbors, j)
        active.pop(i)
        active.pop(j)

        new_id = next_id
        next_id += 1
        active[new_id] = merged
        neighbors[new_id] = set()

        for nb in sorted(candidate_neighbors):
            if nb not in active:
                continue
            stats.local_adjacency_updates += 1
            if are_adjacent(merged.polygon, active[nb].polygon):
                neighbors[new_id].add(nb)
                neighbors.setdefault(nb, set()).add(new_id)
                _enqueue_pair(queue, queued, new_id, nb)

    final_regions = [active[i] for i in sorted(active)]
    stats.final_faces = len(final_regions)
    return final_regions, log, stats


def hm_diagonal_premerge(
    initial_pieces: Sequence[Polygon],
    log_prefix: str = "HM diagonal deletion",
) -> Tuple[List[Region], List[str]]:
    """Compatibility wrapper returning only regions and log."""
    regions, log, _stats = hm_diagonal_deletion_from_triangles(
        initial_pieces,
        log_prefix=log_prefix,
    )
    return regions, log

