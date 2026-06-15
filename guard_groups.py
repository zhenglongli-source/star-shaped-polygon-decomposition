"""
Optional post-processing for the star-shaped decomposition star-shaped decomposition project.

After the main algorithm finishes, we have final star-shaped regions

    R1, R2, ..., Rk

and each region Ri has a kernel Ki. Placing one guard anywhere in Ki guards Ri.

This module tries to reduce the number of actual guard points by checking
whether several final regions have kernels with a nonempty common intersection.
If

    K1 ∩ K2 ∩ ... ∩ Km != empty,

then one guard placed in that common intersection guards all of those regions.

Important:
    This module does NOT merge polygonal regions.
    It only groups final regions for guard placement.

Why this is optional:
    Sharing one guard can reduce the number of guards, but it also restricts
    the feasible guard-placement area from Ki and Kj to Ki ∩ Kj, which may be
    much smaller, or even only a point/line segment.
"""

from dataclasses import dataclass
from typing import List, Sequence, Tuple, Optional

from shapely.geometry import Polygon as ShapelyPolygon
from shapely.geometry import Point as ShapelyPoint
from shapely.geometry import LineString, MultiPoint, MultiLineString, GeometryCollection
from shapely.geometry.base import BaseGeometry


Point = Tuple[float, float]
Polygon = List[Point]


@dataclass
class GuardGroup:
    """
    A group of final star-shaped regions that can share one guard.

    region_indices:
        Indices of the final regions included in this group.

    common_kernel:
        A representative coordinate list for the common kernel intersection.
        This can be:
            - a polygon boundary if the common intersection has positive area;
            - a line segment if the common intersection is one-dimensional;
            - a single point if the common intersection is zero-dimensional.

    geom_type:
        Shapely geometry type of the common kernel, such as:
            "Polygon", "LineString", or "Point".

    area:
        Area of the common kernel. This is 0 for point/line intersections.
    """

    region_indices: List[int]
    common_kernel: Polygon
    geom_type: str
    area: float


def _to_shapely_kernel(kernel: Sequence[Point]) -> Optional[BaseGeometry]:
    """
    Convert a kernel coordinate list into a Shapely geometry.

    Kernels from compute_kernel are usually polygons. However, after intersection,
    a common kernel can become a point or line segment. For the initial region
    kernels, we expect at least 3 points.
    """
    if kernel is None or len(kernel) == 0:
        return None

    try:
        if len(kernel) == 1:
            x, y = kernel[0]
            return ShapelyPoint(x, y)

        if len(kernel) == 2:
            return LineString(kernel)

        geom = ShapelyPolygon(kernel)

        if geom.is_empty:
            return None

        # buffer(0) often fixes minor numerical/self-touching issues.
        geom = geom.buffer(0)

        if geom.is_empty:
            return None

        return geom

    except Exception:
        return None


def _coords_from_geometry(geom: BaseGeometry) -> Polygon:
    """
    Convert a Shapely geometry into a representative coordinate list.

    This does not try to preserve all components of a Multi geometry. For this
    project, we only need a concrete representative guardable set to report.
    """
    if geom is None or geom.is_empty:
        return []

    if isinstance(geom, ShapelyPoint):
        return [(geom.x, geom.y)]

    if isinstance(geom, LineString):
        return list(geom.coords)

    if isinstance(geom, ShapelyPolygon):
        return list(geom.exterior.coords)[:-1]

    if isinstance(geom, MultiPoint):
        pts = list(geom.geoms)
        if not pts:
            return []
        p = pts[0]
        return [(p.x, p.y)]

    if isinstance(geom, MultiLineString):
        lines = list(geom.geoms)
        if not lines:
            return []
        return list(lines[0].coords)

    if isinstance(geom, GeometryCollection):
        # Prefer polygon, then line, then point.
        for target_type in ("Polygon", "LineString", "Point"):
            for g in geom.geoms:
                if g.geom_type == target_type:
                    return _coords_from_geometry(g)

    # Fallback for unusual geometry types.
    if hasattr(geom, "geoms"):
        for g in geom.geoms:
            result = _coords_from_geometry(g)
            if result:
                return result

    return []


def _valid_intersection(geom: BaseGeometry, allow_touch: bool, eps: float) -> bool:
    """
    Decide whether an intersection is considered usable for sharing one guard.

    If allow_touch=True:
        A positive-area polygon, line segment, or point intersection is accepted.

    If allow_touch=False:
        Only positive-area polygon intersection is accepted.
    """
    if geom is None or geom.is_empty:
        return False

    if not allow_touch:
        return geom.geom_type in ("Polygon", "MultiPolygon") and geom.area > eps

    # With allow_touch=True, a point/segment intersection is enough in theory.
    return True


def intersect_kernel_geometries(
    kernel_a: Sequence[Point],
    kernel_b: Sequence[Point],
    allow_touch: bool = True,
    eps: float = 1e-9,
) -> Optional[BaseGeometry]:
    """
    Compute the intersection of two kernels.

    Returns:
        A Shapely geometry if the intersection is usable.
        None otherwise.
    """
    A = _to_shapely_kernel(kernel_a)
    B = _to_shapely_kernel(kernel_b)

    if A is None or B is None:
        return None

    inter = A.intersection(B)

    if not _valid_intersection(inter, allow_touch=allow_touch, eps=eps):
        return None

    return inter


def greedy_kernel_guard_groups(
    regions,
    allow_touch: bool = True,
    eps: float = 1e-9,
) -> List[GuardGroup]:
    """
    Greedily group final regions whose kernels have a common intersection.

    Input:
        regions:
            List of StarRegion objects from star_merge.py.
            Each object should have a .kernel field.

        allow_touch:
            If True, point or line-segment kernel intersections are accepted.
            If False, only positive-area kernel intersections are accepted.

    Output:
        A list of GuardGroup objects.

    Correctness idea:
        For each group, we maintain the common intersection of all kernels
        currently inside that group. A new region can join the group only if
        its kernel intersects this maintained common kernel. Therefore every
        returned group has a nonempty common kernel intersection.
    """
    groups: List[GuardGroup] = []

    for idx, region in enumerate(regions):
        kernel = region.kernel

        kernel_geom = _to_shapely_kernel(kernel)
        if kernel_geom is None or kernel_geom.is_empty:
            groups.append(
                GuardGroup(
                    region_indices=[idx],
                    common_kernel=[],
                    geom_type="Empty",
                    area=0.0,
                )
            )
            continue

        placed = False

        for group in groups:
            if not group.common_kernel:
                continue

            inter = intersect_kernel_geometries(
                group.common_kernel,
                kernel,
                allow_touch=allow_touch,
                eps=eps,
            )

            if inter is not None:
                group.region_indices.append(idx)
                group.common_kernel = _coords_from_geometry(inter)
                group.geom_type = inter.geom_type
                group.area = inter.area
                placed = True
                break

        if not placed:
            groups.append(
                GuardGroup(
                    region_indices=[idx],
                    common_kernel=_coords_from_geometry(kernel_geom),
                    geom_type=kernel_geom.geom_type,
                    area=kernel_geom.area,
                )
            )

    return groups


def print_guard_group_summary(groups: List[GuardGroup]) -> None:
    """
    Print a compact summary of guard groups.
    """
    print(f"Optional shared-guard groups: {len(groups)}")

    for i, group in enumerate(groups):
        print(
            f"  Guard group {i}: regions={group.region_indices}, "
            f"common kernel type={group.geom_type}, "
            f"common kernel vertices={len(group.common_kernel)}, "
            f"area={group.area:.6g}"
        )

