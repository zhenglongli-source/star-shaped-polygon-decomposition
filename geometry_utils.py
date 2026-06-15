from typing import List, Tuple, Optional

Point = Tuple[float, float]
Polygon = List[Point]


def cross(a: Point, b: Point, c: Point) -> float:
    """
    Cross product of AB x AC.
    Positive if C is to the left of directed line AB.
    """
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def polygon_area(poly: Polygon) -> float:
    """
    Signed area. Positive if polygon is CCW.
    """
    n = len(poly)
    area = 0.0
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        area += x1 * y2 - x2 * y1
    return area / 2.0


def ensure_ccw(poly: Polygon) -> Polygon:
    """
    Return polygon in CCW order.
    """
    if polygon_area(poly) < 0:
        return list(reversed(poly))
    return poly[:]


def line_intersection(p1: Point, p2: Point, q1: Point, q2: Point) -> Optional[Point]:
    """
    Intersection point of lines p1-p2 and q1-q2.
    Returns None if lines are parallel.
    """
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = q1
    x4, y4 = q2

    den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if abs(den) < 1e-12:
        return None

    px = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / den
    py = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / den
    return (px, py)


def inside_halfplane(a: Point, b: Point, p: Point, eps: float = 1e-9) -> bool:
    """
    Check whether point p lies in the left closed half-plane of directed edge a->b.
    For CCW polygons, interior is on the left side of each edge.
    """
    return cross(a, b, p) >= -eps

