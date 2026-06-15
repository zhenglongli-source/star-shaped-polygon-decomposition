"""
decomposition-aware Lee-Preparata-inspired support-scan kernel computation.

This module is independent from the baseline implementation.

Core idea:
    Maintain the current kernel K_i and two boundary handles h_i, l_i.
    When a new interior half-plane is added, try to find the affected
    boundary chain by scanning from the previous handles instead of always
    restarting from vertex 0.

Important:
    This is a bounded support-scan adaptation, not a full exact reproduction
    of every unbounded case in the Lee-Preparata paper.

    Full Lee-Preparata:
        - handles unbounded intermediate kernels
        - uses exact F_i, L_i support points
        - has detailed reflex/convex cases

    This version:
        - keeps a bounded polygonal kernel
        - maintains engineering handles h_idx, l_idx
        - stops once the affected chain is found
        - falls back to full clipping only for zero-crossing / degenerate cases
"""

from dataclasses import dataclass
from typing import List, Tuple, Optional

from geometry_utils import (
    Point,
    Polygon,
    ensure_ccw,
    polygon_area,
    inside_halfplane,
    line_intersection,
)

EPS = 1e-9


@dataclass
class SupportState:
    h_idx: int = 0
    l_idx: int = 0


@dataclass
class SupportScanStep:
    edge_index: int
    edge: Tuple[Point, Point]

    h_before: int
    l_before: int
    h_after: int
    l_after: int

    before_vertices: int
    after_vertices: int

    scanned_edges: int
    used_support_scan: bool
    fallback_used: bool

    changed: bool
    empty: bool
    note: str


@dataclass
class Crossing:
    edge_curr_index: int
    point: Point
    entering: bool
    scanned_to_find: int


def remove_nearly_duplicate_points(poly: Polygon, eps: float = EPS) -> Polygon:
    if not poly:
        return []

    cleaned: Polygon = []
    for p in poly:
        if not cleaned:
            cleaned.append(p)
            continue

        q = cleaned[-1]
        if abs(p[0] - q[0]) > eps or abs(p[1] - q[1]) > eps:
            cleaned.append(p)

    if len(cleaned) > 1:
        first = cleaned[0]
        last = cleaned[-1]
        if abs(first[0] - last[0]) <= eps and abs(first[1] - last[1]) <= eps:
            cleaned.pop()

    return cleaned


def is_degenerate_polygon(poly: Polygon, eps: float = EPS) -> bool:
    return len(poly) < 3 or abs(polygon_area(poly)) < eps


def full_clip_kernel(kernel: Polygon, a: Point, b: Point, eps: float = EPS) -> Tuple[Polygon, bool]:
    """
    Robust full clipping fallback.

    This is essentially the baseline Sutherland-Hodgman style update.
    It is used only when the support-scan path cannot safely classify
    the update.
    """
    if not kernel:
        return [], False

    clipped: Polygon = []
    changed = False
    n = len(kernel)

    for i in range(n):
        curr = kernel[i]
        prev = kernel[i - 1]

        curr_in = inside_halfplane(a, b, curr, eps)
        prev_in = inside_halfplane(a, b, prev, eps)

        if prev_in and curr_in:
            clipped.append(curr)

        elif prev_in and not curr_in:
            inter = line_intersection(prev, curr, a, b)
            if inter is not None:
                clipped.append(inter)
            changed = True

        elif not prev_in and curr_in:
            inter = line_intersection(prev, curr, a, b)
            if inter is not None:
                clipped.append(inter)
            clipped.append(curr)
            changed = True

        else:
            changed = True

    clipped = remove_nearly_duplicate_points(clipped, eps)

    if is_degenerate_polygon(clipped, eps):
        return [], True

    return clipped, changed


def find_next_crossing(
    kernel: Polygon,
    a: Point,
    b: Point,
    start_curr_index: int,
    eps: float = EPS,
) -> Optional[Crossing]:
    """
    Scan forward on the current kernel boundary starting from start_curr_index.

    The edge tested at index i is:
        kernel[i-1] -> kernel[i]

    A crossing occurs when the endpoints are on different sides of the
    half-plane boundary.
    """
    n = len(kernel)
    if n == 0:
        return None

    start_curr_index %= n

    for step in range(n):
        curr_idx = (start_curr_index + step) % n
        prev_idx = (curr_idx - 1) % n

        prev = kernel[prev_idx]
        curr = kernel[curr_idx]

        prev_in = inside_halfplane(a, b, prev, eps)
        curr_in = inside_halfplane(a, b, curr, eps)

        if prev_in != curr_in:
            inter = line_intersection(prev, curr, a, b)
            if inter is None:
                return None

            # entering means outside -> inside as we follow kernel boundary
            entering = (not prev_in) and curr_in

            return Crossing(
                edge_curr_index=curr_idx,
                point=inter,
                entering=entering,
                scanned_to_find=step + 1,
            )

    return None


def cyclic_vertices_between(
    kernel: Polygon,
    start_curr_index: int,
    end_curr_index: int,
) -> Polygon:
    """
    Return vertices from start_curr_index up to end_curr_index - 1 cyclically.

    Used after an entering crossing and before an exiting crossing.
    """
    n = len(kernel)
    result: Polygon = []

    i = start_curr_index % n
    stop = end_curr_index % n

    while i != stop:
        result.append(kernel[i])
        i = (i + 1) % n

    return result


def support_scan_clip(
    kernel: Polygon,
    a: Point,
    b: Point,
    state: SupportState,
    eps: float = EPS,
) -> Tuple[Polygon, SupportState, int, bool, bool, str]:
    """
    Clip current kernel by the left half-plane of directed edge a -> b.

    Returns:
        new_kernel
        new_state
        scanned_edges
        changed
        fallback_used
        note

    Main support-scan path:
        1. Start from h_idx and find first crossing.
        2. Continue from after first crossing and find second crossing.
        3. If two crossings are found, keep the inside chain and stop.

    Fallback:
        Used when there is no crossing, only one crossing, degeneracy,
        or inconsistent entering/exiting order.
    """
    if not kernel:
        return [], state, 0, False, False, "Kernel already empty."

    n = len(kernel)
    h_start = state.h_idx % n
    l_start = state.l_idx % n

    # First try from h_idx.
    first = find_next_crossing(kernel, a, b, h_start, eps)
    scanned_edges = 0

    if first is not None:
        scanned_edges += first.scanned_to_find

        second_start = (first.edge_curr_index + 1) % n
        second = find_next_crossing(kernel, a, b, second_start, eps)

        if second is not None:
            scanned_edges += second.scanned_to_find

            c1 = first
            c2 = second

            # We need one entering crossing and one exiting crossing.
            if c1.entering == c2.entering:
                new_kernel, changed = full_clip_kernel(kernel, a, b, eps)
                return (
                    new_kernel,
                    SupportState(0, max(0, len(new_kernel) - 1)),
                    scanned_edges + n,
                    changed,
                    True,
                    "Degenerate/inconsistent crossings; used full clipping fallback.",
                )

            if c1.entering:
                enter = c1
                exit_ = c2
            else:
                enter = c2
                exit_ = c1

            kept_vertices = cyclic_vertices_between(
                kernel,
                enter.edge_curr_index,
                exit_.edge_curr_index,
            )

            new_kernel: Polygon = [enter.point] + kept_vertices + [exit_.point]
            new_kernel = remove_nearly_duplicate_points(new_kernel, eps)

            if is_degenerate_polygon(new_kernel, eps):
                return (
                    [],
                    SupportState(0, 0),
                    scanned_edges,
                    True,
                    False,
                    "Support scan produced empty or degenerate kernel.",
                )

            # New handles are placed near the two new intersection points.
            # Since the new polygon starts with enter.point and ends with exit_.point:
            # h_idx = 0, l_idx = last vertex.
            new_state = SupportState(
                h_idx=0,
                l_idx=len(new_kernel) - 1,
            )

            return (
                new_kernel,
                new_state,
                scanned_edges,
                True,
                False,
                "Support scan found two crossings and clipped only the affected chain.",
            )

    # If h_idx did not find crossing, try l_idx.
    # This still uses the maintained handle pair, rather than vertex 0.
    first = find_next_crossing(kernel, a, b, l_start, eps)

    if first is not None:
        scanned_edges += first.scanned_to_find

        second_start = (first.edge_curr_index + 1) % n
        second = find_next_crossing(kernel, a, b, second_start, eps)

        if second is not None:
            scanned_edges += second.scanned_to_find

            c1 = first
            c2 = second

            if c1.entering == c2.entering:
                new_kernel, changed = full_clip_kernel(kernel, a, b, eps)
                return (
                    new_kernel,
                    SupportState(0, max(0, len(new_kernel) - 1)),
                    scanned_edges + n,
                    changed,
                    True,
                    "Degenerate/inconsistent crossings; used full clipping fallback.",
                )

            if c1.entering:
                enter = c1
                exit_ = c2
            else:
                enter = c2
                exit_ = c1

            kept_vertices = cyclic_vertices_between(
                kernel,
                enter.edge_curr_index,
                exit_.edge_curr_index,
            )

            new_kernel = [enter.point] + kept_vertices + [exit_.point]
            new_kernel = remove_nearly_duplicate_points(new_kernel, eps)

            if is_degenerate_polygon(new_kernel, eps):
                return (
                    [],
                    SupportState(0, 0),
                    scanned_edges,
                    True,
                    False,
                    "Support scan produced empty or degenerate kernel.",
                )

            new_state = SupportState(
                h_idx=0,
                l_idx=len(new_kernel) - 1,
            )

            return (
                new_kernel,
                new_state,
                scanned_edges,
                True,
                False,
                "Support scan from l_idx found two crossings.",
            )

    # Zero-crossing case:
    # The new half-plane either contains the whole current kernel,
    # excludes the whole current kernel, or we missed a degenerate boundary case.
    # For correctness, certify it by full clipping.
    new_kernel, changed = full_clip_kernel(kernel, a, b, eps)

    if not new_kernel:
        note = "Zero-crossing case certified by full clipping: kernel becomes empty."
    elif changed:
        note = "Zero-crossing/degenerate case required full clipping fallback."
    else:
        note = "Zero-crossing case certified: half-plane did not change kernel."

    return (
        new_kernel,
        SupportState(0, max(0, len(new_kernel) - 1)),
        scanned_edges + n,
        changed,
        True,
        note,
    )


def ordered_edges(poly: Polygon) -> List[Tuple[Point, Point]]:
    poly = ensure_ccw(poly)
    n = len(poly)
    return [(poly[i], poly[(i + 1) % n]) for i in range(n)]


def compute_kernel_support_scan(
    poly: Polygon,
    return_trace: bool = False,
) -> Tuple[Polygon, List[SupportScanStep]]:
    """
    Compute kernel using the bounded support-scan adaptation.

    Initialization:
        We start from the polygon itself as a bounded candidate.
        This is safe because K(P) is a subset of P.

    Then each boundary edge contributes an interior half-plane.
    """
    poly = ensure_ccw(poly)
    kernel: Polygon = poly[:]

    state = SupportState(h_idx=0, l_idx=max(0, len(kernel) - 1))
    trace: List[SupportScanStep] = []

    for edge_index, (a, b) in enumerate(ordered_edges(poly)):
        before_vertices = len(kernel)
        h_before = state.h_idx
        l_before = state.l_idx

        new_kernel, new_state, scanned_edges, changed, fallback_used, note = support_scan_clip(
            kernel=kernel,
            a=a,
            b=b,
            state=state,
        )

        after_vertices = len(new_kernel)
        empty = after_vertices == 0

        trace.append(
            SupportScanStep(
                edge_index=edge_index,
                edge=(a, b),
                h_before=h_before,
                l_before=l_before,
                h_after=new_state.h_idx,
                l_after=new_state.l_idx,
                before_vertices=before_vertices,
                after_vertices=after_vertices,
                scanned_edges=scanned_edges,
                used_support_scan=not fallback_used,
                fallback_used=fallback_used,
                changed=changed,
                empty=empty,
                note=note,
            )
        )

        kernel = new_kernel
        state = new_state

        if empty:
            break

    if return_trace:
        return kernel, trace

    return kernel, []


def kernel_area(kernel: Polygon) -> float:
    if not kernel:
        return 0.0
    return abs(polygon_area(kernel))


def is_star_shaped_support_scan(poly: Polygon) -> bool:
    kernel, _ = compute_kernel_support_scan(poly, return_trace=False)
    return len(kernel) > 0

