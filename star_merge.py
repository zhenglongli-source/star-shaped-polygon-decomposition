"""
Greedy star-shaped merging for the star-shaped polygon decomposition project.

Pipeline idea:
1. Start from a triangulation of the input simple polygon.
2. Optional HM-inspired convex pre-merge: greedily merge adjacent pieces when
   their union remains convex.  This plays the same practical role as removing
   inessential diagonals before the star-shaped stage.
3. Greedily merge adjacent pieces when their union remains star-shaped.
   Star-shapedness is checked by the kernel computation kernel computation.
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple
import copy

from shapely.geometry import Polygon as ShapelyPolygon
from shapely.geometry.base import BaseGeometry

from geometry_utils import Point, Polygon, ensure_ccw, polygon_area, cross
from kernel_adapter import compute_kernel, kernel_area
# from kernel_baseline import compute_kernel


EPS = 1e-8


# ---------------------------------------------------------------------------
# Lightweight instrumentation
# ---------------------------------------------------------------------------
# The main algorithm can run without looking at these counters.  They are used
# to diagnose where the merge stage spends work after replacing global pair
# enumeration with adjacency-driven local updates.
_MERGE_STATS = {"total": {}, "by_stage": {}}
_CURRENT_STAGE = "global"


def reset_merge_stats() -> None:
    """Reset diagnostic counters for merge experiments."""
    global _MERGE_STATS
    _MERGE_STATS = {"total": {}, "by_stage": {}}


def get_merge_stats():
    """Return a deep copy of diagnostic counters."""
    return copy.deepcopy(_MERGE_STATS)


def _set_stage(stage: str) -> None:
    global _CURRENT_STAGE
    _CURRENT_STAGE = stage
    _MERGE_STATS["by_stage"].setdefault(stage, {})


def _stat(name: str, amount: int = 1) -> None:
    """Increment a diagnostic counter for both total and current stage."""
    total = _MERGE_STATS["total"]
    total[name] = total.get(name, 0) + amount

    stage_dict = _MERGE_STATS["by_stage"].setdefault(_CURRENT_STAGE, {})
    stage_dict[name] = stage_dict.get(name, 0) + amount


@dataclass
class Region:
    """A polygonal region, together with its kernel when star-shaped."""
    polygon: Polygon
    kernel: Polygon
    source_count: int = 1


def clean_polygon(poly: Polygon, eps: float = EPS) -> Polygon:
    """Remove repeated closing point and nearly repeated consecutive vertices."""
    if not poly:
        return []

    cleaned: Polygon = []
    for p in poly:
        if not cleaned:
            cleaned.append(p)
            continue
        if abs(p[0] - cleaned[-1][0]) > eps or abs(p[1] - cleaned[-1][1]) > eps:
            cleaned.append(p)

    if len(cleaned) > 1:
        first = cleaned[0]
        last = cleaned[-1]
        if abs(first[0] - last[0]) <= eps and abs(first[1] - last[1]) <= eps:
            cleaned.pop()

    return ensure_ccw(cleaned)


def to_shapely(poly: Polygon) -> ShapelyPolygon:
    """Convert a point-list polygon to a Shapely polygon."""
    return ShapelyPolygon(poly)


def from_shapely_polygon(geom: BaseGeometry) -> Optional[Polygon]:
    """
    Convert a Shapely polygon union back to a simple exterior boundary.

    Returns None if the geometry is not a single simple polygon without holes.
    This keeps the project focused on simple polygon pieces.
    """
    if geom.is_empty or geom.geom_type != "Polygon":
        return None
    if len(getattr(geom, "interiors", [])) > 0:
        return None

    coords = list(geom.exterior.coords)
    poly = [(float(x), float(y)) for x, y in coords[:-1]]
    poly = clean_polygon(poly)
    if len(poly) < 3 or abs(polygon_area(poly)) < EPS:
        return None
    return poly


def _point_key(p: Point, ndigits: int = 10) -> Tuple[float, float]:
    """Rounded coordinate key used for combinatorial edge matching."""
    return (round(float(p[0]), ndigits), round(float(p[1]), ndigits))


def _edge_key(a: Point, b: Point) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """Undirected rounded edge key."""
    ka = _point_key(a)
    kb = _point_key(b)
    return (ka, kb) if ka <= kb else (kb, ka)


def _edges(poly: Polygon) -> List[Tuple[Point, Point]]:
    """Directed boundary edges of a polygon."""
    poly = clean_polygon(poly)
    return [(poly[i], poly[(i + 1) % len(poly)]) for i in range(len(poly))]


def _dist(a: Point, b: Point) -> float:
    return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) ** 0.5


def _same_point(a: Point, b: Point, eps: float = EPS) -> bool:
    return abs(a[0] - b[0]) <= eps and abs(a[1] - b[1]) <= eps


def _collinear(a: Point, b: Point, c: Point, eps: float = EPS) -> bool:
    return abs(cross(a, b, c)) <= eps


def _segment_overlap_length(a: Point, b: Point, c: Point, d: Point, eps: float = EPS) -> float:
    """
    Positive length of overlap between two collinear segments.

    This supports the adjacency predicate even when a shared boundary is split
    into smaller collinear pieces on one side.
    """
    if not (_collinear(a, b, c, eps) and _collinear(a, b, d, eps)):
        return 0.0

    dx = b[0] - a[0]
    dy = b[1] - a[1]
    use_x = abs(dx) >= abs(dy)

    if use_x:
        a1, b1 = sorted((a[0], b[0]))
        c1, d1 = sorted((c[0], d[0]))
    else:
        a1, b1 = sorted((a[1], b[1]))
        c1, d1 = sorted((c[1], d[1]))

    lo = max(a1, c1)
    hi = min(b1, d1)
    if hi <= lo + eps:
        return 0.0

    # Convert projected overlap back to Euclidean length.
    seg_len = _dist(a, b)
    proj_len = abs((b[0] - a[0]) if use_x else (b[1] - a[1]))
    if proj_len <= eps:
        return 0.0
    return (hi - lo) * seg_len / proj_len


def shared_boundary_length(poly_a: Polygon, poly_b: Polygon, eps: float = EPS) -> float:
    """Total positive-length boundary overlap between two polygonal regions."""
    length = 0.0
    edges_a = _edges(poly_a)
    edges_b = _edges(poly_b)

    # Fast path: exact matching undirected edges.
    b_by_key = {_edge_key(c, d): (c, d) for c, d in edges_b}
    for a, b in edges_a:
        if _edge_key(a, b) in b_by_key:
            length += _dist(a, b)

    if length > eps:
        return length

    # Robust path: collinear partial overlaps.
    for a, b in edges_a:
        for c, d in edges_b:
            length += _segment_overlap_length(a, b, c, d, eps)

    return length


def are_adjacent_shapely(poly_a: Polygon, poly_b: Polygon, eps: float = EPS) -> bool:
    """Original Shapely-based adjacency test, kept for debugging/validation."""
    a = to_shapely(poly_a)
    b = to_shapely(poly_b)
    inter = a.boundary.intersection(b.boundary)
    return inter.length > eps


def are_adjacent(poly_a: Polygon, poly_b: Polygon, eps: float = EPS) -> bool:
    """
    Two pieces are adjacent if their boundaries share a positive-length segment.

    This replaces the previous Shapely boundary-intersection query with a
    DCEL-style shared-edge predicate.  It still allows collinear split edges.
    """
    _stat("adjacency_tests")
    result = shared_boundary_length(poly_a, poly_b, eps) > eps
    if result:
        _stat("adjacency_true")
    return result


def _union_by_boundary_edges(poly_a: Polygon, poly_b: Polygon, eps: float = EPS) -> Optional[Polygon]:
    """
    Merge two adjacent simple pieces by deleting their shared internal edges.

    This is a lightweight DCEL-style boundary update. It works when the shared
    boundary is represented by matching directed edges with the same endpoints
    in the two polygons. If the boundary is split differently or the result is
    not a single simple exterior cycle, return None and let the Shapely fallback
    handle it.
    """
    poly_a = clean_polygon(poly_a)
    poly_b = clean_polygon(poly_b)
    if len(poly_a) < 3 or len(poly_b) < 3:
        return None

    all_edges = _edges(poly_a) + _edges(poly_b)
    counts: Dict[Tuple[Tuple[float, float], Tuple[float, float]], int] = {}
    for u, v in all_edges:
        key = _edge_key(u, v)
        counts[key] = counts.get(key, 0) + 1

    shared_keys = {key for key, count in counts.items() if count == 2}
    if not shared_keys:
        return None

    boundary_edges = [(u, v) for u, v in all_edges if _edge_key(u, v) not in shared_keys]
    if len(boundary_edges) < 3:
        return None

    outgoing: Dict[Tuple[float, float], List[Tuple[Point, Point]]] = {}
    for u, v in boundary_edges:
        outgoing.setdefault(_point_key(u), []).append((u, v))

    # A simple exterior cycle has exactly one outgoing edge per boundary vertex.
    if any(len(items) != 1 for items in outgoing.values()):
        return None

    start = boundary_edges[0][0]
    start_key = _point_key(start)
    current_key = start_key
    result: Polygon = []
    used = set()

    for _ in range(len(boundary_edges) + 1):
        if current_key not in outgoing:
            return None
        u, v = outgoing[current_key][0]
        edge_id = (_point_key(u), _point_key(v))
        if edge_id in used:
            return None
        used.add(edge_id)
        result.append(u)
        current_key = _point_key(v)
        if current_key == start_key:
            break
    else:
        return None

    if len(used) != len(boundary_edges):
        return None

    result = clean_polygon(result, eps)
    if len(result) < 3 or abs(polygon_area(result)) < eps:
        return None

    # Cheap validity checks.  Use Shapely only as a validator, not as the union
    # constructor.  The area should equal the sum of the two non-overlapping
    # pieces because adjacent regions share only boundary.
    geom = ShapelyPolygon(result)
    if geom.is_empty or not geom.is_valid or len(getattr(geom, "interiors", [])) > 0:
        return None

    area_result = abs(polygon_area(result))
    area_expected = abs(polygon_area(poly_a)) + abs(polygon_area(poly_b))
    if abs(area_result - area_expected) > max(1.0, area_expected) * 1e-7:
        return None

    return result


def union_if_simple_shapely(poly_a: Polygon, poly_b: Polygon) -> Optional[Polygon]:
    """Original Shapely union, kept as a robust fallback/reference."""
    geom = to_shapely(poly_a).union(to_shapely(poly_b))
    return from_shapely_polygon(geom)


def union_if_simple(poly_a: Polygon, poly_b: Polygon) -> Optional[Polygon]:
    """
    Return the simple polygon union of two pieces, or None if invalid.

    Prefer a lightweight boundary-chain merge.  Fall back to Shapely only when
    the shared boundary is not represented in a clean DCEL-like way.
    """
    _stat("union_attempts")
    _stat("boundary_union_attempts")
    merged = _union_by_boundary_edges(poly_a, poly_b)
    if merged is not None:
        _stat("boundary_union_success")
        return merged

    _stat("shapely_union_fallback")
    merged = union_if_simple_shapely(poly_a, poly_b)
    if merged is not None:
        _stat("shapely_union_success")
    return merged


def is_convex_polygon(poly: Polygon, eps: float = EPS) -> bool:
    """Return True if a polygon is convex or weakly convex."""
    poly = clean_polygon(poly)
    n = len(poly)
    if n < 3:
        return False

    seen_positive = False
    seen_negative = False
    for i in range(n):
        a = poly[i]
        b = poly[(i + 1) % n]
        c = poly[(i + 2) % n]
        z = cross(a, b, c)
        if z > eps:
            seen_positive = True
        elif z < -eps:
            seen_negative = True
        if seen_positive and seen_negative:
            return False
    return True


def make_region(poly: Polygon) -> Optional[Region]:
    """Create a Region if poly is star-shaped; otherwise return None."""
    poly = clean_polygon(poly)
    kernel = compute_kernel(poly)
    if not kernel:
        return None
    return Region(polygon=poly, kernel=kernel, source_count=1)


def make_regions(initial_pieces: List[Polygon]) -> List[Region]:
    """Convert initial pieces into star-shaped Region objects."""
    regions: List[Region] = []
    for piece in initial_pieces:
        region = make_region(piece)
        if region is None:
            raise ValueError("Initial piece is not star-shaped. This should not happen for triangles/convex pieces.")
        regions.append(region)
    return regions


def try_convex_merge(a: Region, b: Region, assume_adjacent: bool = False) -> Optional[Region]:
    """Merge two adjacent regions if their union is convex."""
    if not assume_adjacent and not are_adjacent(a.polygon, b.polygon):
        return None
    union_poly = union_if_simple(a.polygon, b.polygon)
    if union_poly is None:
        return None
    _stat("convexity_tests")
    if not is_convex_polygon(union_poly):
        return None

    # For a convex polygon, the kernel is the polygon itself.  We still compute
    # it through the same kernel computation routine for consistency with the rest of the code.
    _stat("kernel_tests")
    kernel = compute_kernel(union_poly)
    if not kernel:
        return None
    return Region(union_poly, kernel, a.source_count + b.source_count)


def try_star_merge(a: Region, b: Region, assume_adjacent: bool = False) -> Optional[Region]:
    """Merge two adjacent regions if their union remains star-shaped."""
    if not assume_adjacent and not are_adjacent(a.polygon, b.polygon):
        return None

    union_poly = union_if_simple(a.polygon, b.polygon)
    if union_poly is None:
        return None

    _stat("kernel_tests")
    kernel = compute_kernel(union_poly)
    if not kernel:
        return None

    return Region(
        polygon=union_poly,
        kernel=kernel,
        source_count=a.source_count + b.source_count,
    )


def _score_region(region: Region, score: str) -> float:
    """Scoring rule used to choose among valid merge candidates."""
    region_area = abs(polygon_area(region.polygon))
    kernel_area = abs(polygon_area(region.kernel)) if region.kernel else 0.0

    if score == "largest_area":
        return region_area
    if score == "largest_kernel":
        return kernel_area
    if score == "kernel_ratio":
        return kernel_area / region_area if region_area > EPS else 0.0
    if score == "fewest_boundary_vertices":
        # Larger is better in the generic loop, so negate vertex count.
        return -float(len(region.polygon))
    return region_area


def _select_try_merge(merge_type: str):
    """Return the merge predicate for the requested stage."""
    if merge_type == "convex":
        return try_convex_merge
    if merge_type == "star":
        return try_star_merge
    raise ValueError("merge_type must be 'convex' or 'star'.")


def _initial_adjacency_graph(regions: List[Region]) -> Dict[int, Set[int]]:
    """
    Build the initial region adjacency graph.

    This is still an O(k^2) initialization, but it is paid once. The older
    merge loop effectively rediscovered adjacency from scratch in every merge
    round, producing an O(k^3)-style candidate enumeration pattern.
    """
    neighbors: Dict[int, Set[int]] = {i: set() for i in range(len(regions))}
    for i in range(len(regions)):
        for j in range(i + 1, len(regions)):
            if are_adjacent(regions[i].polygon, regions[j].polygon):
                neighbors[i].add(j)
                neighbors[j].add(i)
    return neighbors


def _adjacency_edges(neighbors: Dict[int, Set[int]], active: Dict[int, Region]) -> List[Tuple[int, int]]:
    """Return all active undirected adjacency edges in deterministic order."""
    edges: List[Tuple[int, int]] = []
    active_ids = set(active.keys())
    for i in sorted(active_ids):
        for j in sorted(neighbors.get(i, set())):
            if i < j and j in active_ids:
                edges.append((i, j))
    return edges


def _remove_active_region(neighbors: Dict[int, Set[int]], region_id: int) -> None:
    """Remove one region id from the adjacency graph."""
    for nb in list(neighbors.get(region_id, set())):
        if nb in neighbors:
            neighbors[nb].discard(region_id)
    neighbors.pop(region_id, None)


def greedy_merge_regions_global(
    regions: List[Region],
    merge_type: str,
    strategy: str = "largest_area",
    log_prefix: str = "merge",
) -> Tuple[List[Region], List[str]]:
    """
    Original global all-pairs greedy merge loop.

    Kept as a reference baseline. In each merge round it scans all O(r^2)
    region pairs, even though only adjacent regions can possibly merge.
    """
    _set_stage(log_prefix.replace(" ", "_"))
    try_merge = _select_try_merge(merge_type)
    regions = list(regions)
    log: List[str] = []
    step = 0

    while True:
        best = None  # (score, i, j, merged_region)
        for i in range(len(regions)):
            for j in range(i + 1, len(regions)):
                _stat("candidate_pairs_examined")
                merged = try_merge(regions[i], regions[j])
                if merged is None:
                    continue
                _stat("valid_merge_candidates")
                value = _score_region(merged, strategy)
                if strategy == "first":
                    best = (value, i, j, merged)
                    break
                if best is None or value > best[0]:
                    best = (value, i, j, merged)
            if strategy == "first" and best is not None:
                break

        if best is None:
            break

        _, i, j, merged = best
        step += 1
        _stat("accepted_merges")
        old_count = len(regions)
        log.append(
            f"{log_prefix} {step}: merged regions {i} and {j}; "
            f"pieces {regions[i].source_count}+{regions[j].source_count} -> {merged.source_count}; "
            f"region count {old_count}->{old_count - 1}."
        )
        regions = [region for k, region in enumerate(regions) if k not in (i, j)] + [merged]

    return regions, log


def greedy_merge_regions(
    regions: List[Region],
    merge_type: str,
    strategy: str = "largest_area",
    log_prefix: str = "merge",
    use_adjacency_graph: bool = True,
) -> Tuple[List[Region], List[str]]:
    """
    Generic greedy merge loop.

    With use_adjacency_graph=True, this is an adjacency-driven local update
    version: only adjacent region pairs are tested, and after a merge only the
    neighbors of the new region are recomputed. This is the first step toward a
    DCEL implementation while preserving the old greedy strategy.
    """
    stage_name = log_prefix.replace(" ", "_")
    _set_stage(stage_name)

    if not use_adjacency_graph:
        return greedy_merge_regions_global(regions, merge_type, strategy, log_prefix)

    try_merge = _select_try_merge(merge_type)
    active: Dict[int, Region] = {i: region for i, region in enumerate(regions)}
    neighbors = _initial_adjacency_graph(list(regions))
    _stat("initial_adjacency_edges", sum(len(v) for v in neighbors.values()) // 2)
    next_id = len(regions)
    log: List[str] = []
    step = 0

    while True:
        best = None  # (score, i, j, merged_region)
        for i, j in _adjacency_edges(neighbors, active):
            _stat("candidate_pairs_examined")
            merged = try_merge(active[i], active[j], assume_adjacent=True)
            if merged is None:
                continue
            _stat("valid_merge_candidates")
            value = _score_region(merged, strategy)
            if strategy == "first":
                best = (value, i, j, merged)
                break
            if best is None or value > best[0]:
                best = (value, i, j, merged)

        if best is None:
            break

        _, i, j, merged = best
        step += 1
        _stat("accepted_merges")
        old_count = len(active)
        log.append(
            f"{log_prefix} {step}: merged regions {i} and {j}; "
            f"pieces {active[i].source_count}+{active[j].source_count} -> {merged.source_count}; "
            f"region count {old_count}->{old_count - 1}."
        )

        # Only former neighbors of i or j can become neighbors of the merged region.
        candidate_neighbors = (neighbors.get(i, set()) | neighbors.get(j, set())) - {i, j}

        _remove_active_region(neighbors, i)
        _remove_active_region(neighbors, j)
        active.pop(i)
        active.pop(j)

        new_id = next_id
        next_id += 1
        active[new_id] = merged
        neighbors[new_id] = set()

        for nb in list(candidate_neighbors):
            if nb not in active:
                continue
            _stat("local_adjacency_updates")
            if are_adjacent(merged.polygon, active[nb].polygon):
                neighbors[new_id].add(nb)
                neighbors.setdefault(nb, set()).add(new_id)

    return [active[i] for i in sorted(active.keys())], log

def hm_inspired_convex_premerge(
    initial_pieces: List[Polygon],
    strategy: str = "largest_area",
) -> Tuple[List[Region], List[str]]:
    """
    A practical HM-inspired preprocessing step.

    It starts with triangles and greedily removes a diagonal whenever the union
    of the two adjacent current pieces is convex.  This is not a full DCEL-based
    O(n) implementation of Hertel-Mehlhorn, but it captures the same algorithmic
    idea for this project: simplify the initial triangulation before using the
    star-shaped merge test.
    """
    regions = make_regions(initial_pieces)
    return greedy_merge_regions(regions, merge_type="convex", strategy=strategy, log_prefix="convex pre-merge")


def greedy_star_merge(
    initial_pieces: List[Polygon],
    strategy: str = "largest_area",
) -> Tuple[List[Region], List[str]]:
    """
    Original baseline: start from pieces and greedily merge adjacent pieces when
    their union remains star-shaped.
    """
    regions = make_regions(initial_pieces)
    return greedy_merge_regions(regions, merge_type="star", strategy=strategy, log_prefix="star merge")


def optimized_star_decomposition(
    initial_pieces: List[Polygon],
    convex_strategy: str = "largest_area",
    star_strategy: str = "largest_area",
    use_convex_premerge: bool = True,
) -> Tuple[List[Region], List[str], List[Region]]:
    """
    Optimized pipeline:
        triangulation pieces
        -> optional HM-inspired convex pre-merge
        -> star-shaped merge.

    Returns final regions, combined log, and the intermediate regions after the
    convex pre-merge stage.
    """
    if use_convex_premerge:
        convex_regions, convex_log = hm_inspired_convex_premerge(initial_pieces, strategy=convex_strategy)
    else:
        convex_regions = make_regions(initial_pieces)
        convex_log = []

    final_regions, star_log = greedy_merge_regions(
        convex_regions,
        merge_type="star",
        strategy=star_strategy,
        log_prefix="star merge",
    )
    return final_regions, convex_log + star_log, convex_regions


def compare_star_strategies(
    initial_pieces: List[Polygon],
    strategies: List[str],
    use_convex_premerge: bool = True,
) -> List[Tuple[str, int, int, int]]:
    """
    Small experimental helper.

    Returns rows of:
        (strategy, convex_piece_count, final_region_count, accepted_merge_count)
    """
    rows = []
    for strategy in strategies:
        final_regions, log, convex_regions = optimized_star_decomposition(
            initial_pieces,
            convex_strategy="largest_area",
            star_strategy=strategy,
            use_convex_premerge=use_convex_premerge,
        )
        rows.append((strategy, len(convex_regions), len(final_regions), len(log)))
    return rows


