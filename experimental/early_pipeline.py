"""
Early-pipeline backends for the star-shaped decomposition star-shaped decomposition project.

Backends:
    ear_greedy:
        existing Python ear clipping + old greedy convex pre-merge.

    polypp_mono_hm:
        PolyPartition Triangulate_MONO + project-local HM-style diagonal deletion.

    polypp_hm_reference:
        PolyPartition ConvexPartition_HM direct output.
        This is kept only as an external reference backend, not recommended as main.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

from geometry_utils import Polygon, polygon_area
from triangulation import triangulate_polygon
from star_merge import (
    Region,
    clean_polygon,
    greedy_merge_regions,
    hm_inspired_convex_premerge,
    is_convex_polygon,
    make_regions,
)
from hm_seeded_kernel_merge_optimized import (
    HybridHMSeededStats,
    hm_seeded_star_merge_hybrid,
)
from hm_seeded_kernel_merge import KernelStateRegion
from polypartition_adapter import (
    convex_partition_polypartition_hm,
    partition_area_error,
    triangulate_polypartition_mono,
)
from hm_diagonal_deletion import (
    HMDiagonalDeletionStats,
    hm_diagonal_deletion_from_triangles,
)

EPS = 1e-8


@dataclass
class EarlyPipelineResult:
    backend: str
    triangles: List[Polygon]
    convex_regions: List[Region]
    log: List[str]
    area_error: float = 0.0
    hm_stats: Optional[HMDiagonalDeletionStats] = None


@dataclass
class FullPipelineResult:
    early: EarlyPipelineResult
    final_regions: List[KernelStateRegion]
    full_log: List[str]
    hybrid_stats: HybridHMSeededStats


def _convex_piece_to_region(poly: Polygon) -> Region:
    poly = clean_polygon(poly)
    if len(poly) < 3 or abs(polygon_area(poly)) < EPS:
        raise ValueError("Degenerate convex partition piece.")
    if not is_convex_polygon(poly):
        raise ValueError(f"Expected convex piece, got non-convex polygon: {poly}")
    return Region(polygon=poly, kernel=poly[:], source_count=1)


def build_early_pipeline(
    poly: Polygon,
    backend: str = "ear_greedy",
    convex_strategy: str = "largest_area",
) -> EarlyPipelineResult:
    """
    Build the early decomposition before star-shaped merging.

    Returns convex regions that can seed the HM-seeded portal-kernel backend.
    """
    if backend == "ear_greedy":
        triangles = triangulate_polygon(poly)
        convex_regions, convex_log = hm_inspired_convex_premerge(
            triangles,
            strategy=convex_strategy,
        )
        return EarlyPipelineResult(
            backend=backend,
            triangles=triangles,
            convex_regions=convex_regions,
            log=convex_log,
            area_error=0.0,
            hm_stats=None,
        )

    if backend == "polypp_mono_hm":
        triangles = triangulate_polypartition_mono(poly)
        area_error = partition_area_error(poly, triangles)

        convex_regions, hm_log, hm_stats = hm_diagonal_deletion_from_triangles(
            triangles,
            log_prefix="polypp-mono HM deletion",
        )
        return EarlyPipelineResult(
            backend=backend,
            triangles=triangles,
            convex_regions=convex_regions,
            log=hm_log,
            area_error=area_error,
            hm_stats=hm_stats,
        )

    if backend == "polypp_hm_reference":
        convex_pieces = convex_partition_polypartition_hm(poly)
        area_error = partition_area_error(poly, convex_pieces)
        convex_regions = [_convex_piece_to_region(piece) for piece in convex_pieces]
        return EarlyPipelineResult(
            backend=backend,
            triangles=[],
            convex_regions=convex_regions,
            log=[f"PolyPartition ConvexPartition_HM returned {len(convex_regions)} convex pieces."],
            area_error=area_error,
            hm_stats=None,
        )

    raise ValueError(
        "Unknown early backend. Use one of: "
        "'ear_greedy', 'polypp_mono_hm', 'polypp_hm_reference'."
    )


def run_hybrid_pipeline(
    poly: Polygon,
    early_backend: str = "ear_greedy",
    convex_strategy: str = "largest_area",
    star_strategy: str = "largest_area",
) -> FullPipelineResult:
    """
    Full pipeline:
        early backend -> convex regions -> HM-seeded hybrid star merge.
    """
    early = build_early_pipeline(
        poly,
        backend=early_backend,
        convex_strategy=convex_strategy,
    )

    final_regions, star_log, hybrid_stats = hm_seeded_star_merge_hybrid(
        early.convex_regions,
        original_polygon=poly,
        strategy=star_strategy,
        log_prefix=f"{early_backend} hybrid star merge",
    )

    return FullPipelineResult(
        early=early,
        final_regions=final_regions,
        full_log=early.log + star_log,
        hybrid_stats=hybrid_stats,
    )


def run_standard_star_from_early(
    poly: Polygon,
    early_backend: str = "ear_greedy",
    convex_strategy: str = "largest_area",
    star_strategy: str = "largest_area",
):
    """
    Optional compatibility path:
        early backend -> old standard star merge.

    This is useful for visual comparison with the older main.py pipeline.
    """
    early = build_early_pipeline(
        poly,
        backend=early_backend,
        convex_strategy=convex_strategy,
    )

    final_regions, star_log = greedy_merge_regions(
        early.convex_regions,
        merge_type="star",
        strategy=star_strategy,
        log_prefix=f"{early_backend} standard star merge",
    )

    return final_regions, early.log + star_log, early
