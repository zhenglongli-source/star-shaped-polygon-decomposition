from typing import List, Tuple
from geometry_utils import Point, Polygon, ensure_ccw, inside_halfplane, line_intersection

BIG = 10**6


def clip_polygon_with_halfplane(poly: Polygon, a: Point, b: Point) -> Polygon:
    """
    Clip polygon poly by the left closed half-plane of directed line a->b
    using Sutherland-Hodgman style clipping.
    """
    if not poly:
        return []

    clipped = []
    n = len(poly)

    for i in range(n):
        curr = poly[i]
        prev = poly[i - 1]

        curr_in = inside_halfplane(a, b, curr)
        prev_in = inside_halfplane(a, b, prev)

        if prev_in and curr_in:
            clipped.append(curr)
        elif prev_in and not curr_in:
            inter = line_intersection(prev, curr, a, b)
            if inter is not None:
                clipped.append(inter)
        elif not prev_in and curr_in:
            inter = line_intersection(prev, curr, a, b)
            if inter is not None:
                clipped.append(inter)
            clipped.append(curr)

    return clipped


def initial_bounding_box(poly: Polygon, margin: float = 10.0) -> Polygon:
    xs = [p[0] for p in poly]
    ys = [p[1] for p in poly]
    minx, maxx = min(xs), max(xs)
    miny, maxy = min(ys), max(ys)

    dx = maxx - minx
    dy = maxy - miny
    pad = max(dx, dy, 1.0) * margin

    return [
        (minx - pad, miny - pad),
        (maxx + pad, miny - pad),
        (maxx + pad, maxy + pad),
        (minx - pad, maxy + pad),
    ]


def compute_kernel(poly: Polygon) -> Polygon:
    """
    Baseline method:
    kernel = intersection of all interior half-planes defined by polygon edges.
    """
    poly = ensure_ccw(poly)
    kernel = initial_bounding_box(poly)

    n = len(poly)
    for i in range(n):
        a = poly[i]
        b = poly[(i + 1) % n]
        kernel = clip_polygon_with_halfplane(kernel, a, b)
        if not kernel:
            return []

    return kernel


def is_star_shaped(poly: Polygon) -> bool:
    return len(compute_kernel(poly)) > 0

