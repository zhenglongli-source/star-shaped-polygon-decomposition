"""
Ear-clipping triangulation for a simple polygon.

This file is used as the baseline decomposition stage for the star-shaped polygon decomposition
project.  The triangulation is intentionally simple and readable: it assumes the
input polygon is simple and returns n-2 triangles for a polygon with n vertices.
"""

from typing import List
from geometry_utils import Point, Polygon, cross, ensure_ccw, polygon_area

EPS = 1e-9


def point_in_triangle(p: Point, a: Point, b: Point, c: Point, eps: float = EPS) -> bool:
    """Return True if p lies inside or on the boundary of triangle abc."""
    c1 = cross(a, b, p)
    c2 = cross(b, c, p)
    c3 = cross(c, a, p)
    return c1 >= -eps and c2 >= -eps and c3 >= -eps


def is_convex_vertex(prev: Point, curr: Point, nxt: Point, eps: float = EPS) -> bool:
    """For a CCW polygon, a vertex is convex if the left turn is positive."""
    return cross(prev, curr, nxt) > eps


def is_ear(poly: Polygon, indices: List[int], pos: int) -> bool:
    """Check whether indices[pos] is an ear tip in the current polygon."""
    n = len(indices)
    i_prev = indices[(pos - 1) % n]
    i_curr = indices[pos]
    i_next = indices[(pos + 1) % n]

    a = poly[i_prev]
    b = poly[i_curr]
    c = poly[i_next]

    if not is_convex_vertex(a, b, c):
        return False

    # No other remaining vertex may lie inside the candidate ear triangle.
    for idx in indices:
        if idx in (i_prev, i_curr, i_next):
            continue
        if point_in_triangle(poly[idx], a, b, c):
            return False

    return True


def triangulate_polygon(poly: Polygon) -> List[Polygon]:
    """
    Triangulate a simple polygon using ear clipping.

    Parameters
    ----------
    poly:
        List of polygon vertices in clockwise or counterclockwise order.

    Returns
    -------
    List[Polygon]
        A list of triangles, each represented by three points in CCW order.
    """
    poly = ensure_ccw(poly)
    n = len(poly)
    if n < 3:
        raise ValueError("A polygon needs at least three vertices.")
    if abs(polygon_area(poly)) < EPS:
        raise ValueError("Degenerate polygon with near-zero area.")
    if n == 3:
        return [poly[:]]

    indices = list(range(n))
    triangles: List[Polygon] = []

    # A safety bound prevents infinite loops on degenerate inputs.
    attempts_without_clip = 0
    max_attempts = n * n

    while len(indices) > 3:
        clipped = False
        m = len(indices)
        for pos in range(m):
            if is_ear(poly, indices, pos):
                i_prev = indices[(pos - 1) % m]
                i_curr = indices[pos]
                i_next = indices[(pos + 1) % m]
                triangles.append([poly[i_prev], poly[i_curr], poly[i_next]])
                del indices[pos]
                clipped = True
                attempts_without_clip = 0
                break

        if not clipped:
            attempts_without_clip += 1
            if attempts_without_clip > max_attempts:
                raise RuntimeError(
                    "Ear clipping failed. The input may be degenerate or not a simple polygon."
                )

    triangles.append([poly[indices[0]], poly[indices[1]], poly[indices[2]]])
    return triangles


