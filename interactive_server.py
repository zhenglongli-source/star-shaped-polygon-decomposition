"""
Local backend for the Computational Geometry star-shaped decomposition interactive viewer.

Place this file either:
  A) in the project root that contains submission/star-shaped decomposition_star_decomposition_code_v2/, or
  B) directly next to the star-shaped decomposition_star_decomposition_code_v2/ folder, or
  C) inside star-shaped decomposition_star_decomposition_code_v2/ itself.

Run:
  python interactive_server_updated.py
Open:
  http://localhost:8010/interactive.html
"""

from __future__ import annotations

import json
import os
import sys
import traceback
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

ROOT = os.path.dirname(os.path.abspath(__file__))


def find_code_dir(root: str) -> str:
    candidates = [
        os.path.join(root, "submission", "star-shaped decomposition_star_decomposition_code_v2"),
        os.path.join(root, "star-shaped decomposition_star_decomposition_code_v2"),
        root if os.path.basename(root) == "star-shaped decomposition_star_decomposition_code_v2" else "",
    ]
    for path in candidates:
        if path and os.path.isfile(os.path.join(path, "star_merge.py")):
            return path
    raise RuntimeError(
        "Could not find star-shaped decomposition_star_decomposition_code_v2. Put this server in the project root, "
        "next to star-shaped decomposition_star_decomposition_code_v2, or inside star-shaped decomposition_star_decomposition_code_v2."
    )


CODE_DIR = ROOT
if CODE_DIR not in sys.path:
    sys.path.insert(0, CODE_DIR)

from examples import ALL_EXAMPLES  # type: ignore
from geometry_utils import polygon_area  # type: ignore
from guard_groups import greedy_kernel_guard_groups  # type: ignore
from star_merge import (  # type: ignore
    compare_star_strategies,
    get_merge_stats,
    optimized_star_decomposition,
    reset_merge_stats,
)
from triangulation import triangulate_polygon  # type: ignore
from visualize import polygon_centroid  # type: ignore

Point = Tuple[float, float]
Polygon = List[Point]
EPS = 1e-9


def signed_area(poly: Sequence[Point]) -> float:
    area = 0.0
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        area += x1 * y2 - x2 * y1
    return area / 2.0


def point_to_list(p: Sequence[float]) -> List[float]:
    return [float(p[0]), float(p[1])]


def poly_to_list(poly: Optional[Sequence[Sequence[float]]]) -> List[List[float]]:
    if not poly:
        return []
    return [point_to_list(p) for p in poly]


def same_point(a: Point, b: Point, eps: float = EPS) -> bool:
    return abs(a[0] - b[0]) <= eps and abs(a[1] - b[1]) <= eps


def cross(a: Point, b: Point, c: Point) -> float:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def on_segment(a: Point, b: Point, p: Point, eps: float = EPS) -> bool:
    if abs(cross(a, b, p)) > eps:
        return False
    return (
        min(a[0], b[0]) - eps <= p[0] <= max(a[0], b[0]) + eps
        and min(a[1], b[1]) - eps <= p[1] <= max(a[1], b[1]) + eps
    )


def segment_intersection_kind(a: Point, b: Point, c: Point, d: Point, eps: float = EPS) -> str:
    """Return none, touch, or cross/overlap for closed segments ab and cd."""
    o1 = cross(a, b, c)
    o2 = cross(a, b, d)
    o3 = cross(c, d, a)
    o4 = cross(c, d, b)

    # Proper crossing.
    if ((o1 > eps and o2 < -eps) or (o1 < -eps and o2 > eps)) and (
        (o3 > eps and o4 < -eps) or (o3 < -eps and o4 > eps)
    ):
        return "cross"

    touches = 0
    for p, x, y in ((c, a, b), (d, a, b), (a, c, d), (b, c, d)):
        if on_segment(x, y, p, eps):
            touches += 1

    if touches == 0:
        return "none"
    if touches == 1:
        return "touch"
    return "overlap"


def validate_simple_polygon(poly: Sequence[Point]) -> None:
    n = len(poly)
    if n < 3:
        raise ValueError("A polygon needs at least 3 vertices.")
    if abs(signed_area(poly)) <= EPS:
        raise ValueError("Degenerate polygon: area is too close to zero.")

    for i in range(n):
        if same_point(poly[i], poly[(i + 1) % n]):
            raise ValueError(f"Degenerate polygon: edge {i} has repeated endpoints.")

    # Non-adjacent edges may not intersect. Adjacent edges may share only their endpoint.
    for i in range(n):
        a, b = poly[i], poly[(i + 1) % n]
        for j in range(i + 1, n):
            if i == j:
                continue
            # Skip adjacent edge pairs and the first/last edge pair.
            if (i + 1) % n == j or (j + 1) % n == i:
                continue
            c, d = poly[j], poly[(j + 1) % n]
            kind = segment_intersection_kind(a, b, c, d)
            if kind != "none":
                raise ValueError(
                    f"Input is not a simple polygon: edge {i} intersects edge {j} ({kind})."
                )


def normalize_polygon(raw_poly: Iterable[Sequence[Any]]) -> Polygon:
    clean: Polygon = []
    for raw in raw_poly:
        if len(raw) < 2:
            raise ValueError("Each vertex must contain x and y.")
        p = (float(raw[0]), float(raw[1]))
        if not clean or not same_point(clean[-1], p):
            clean.append(p)

    if len(clean) >= 2 and same_point(clean[0], clean[-1]):
        clean.pop()

    validate_simple_polygon(clean)

    # The algorithm routines use CCW polygons.
    if signed_area(clean) < 0:
        clean.reverse()
    return clean


def region_to_dict(region: Any, index: int) -> Dict[str, Any]:
    guard = polygon_centroid(region.kernel) if region.kernel else None
    return {
        "index": index,
        "polygon": poly_to_list(region.polygon),
        "kernel": poly_to_list(region.kernel),
        "guard": point_to_list(guard) if guard is not None else None,
        "source_count": int(region.source_count),
        "boundary_vertices": len(region.polygon),
        "kernel_vertices": len(region.kernel) if region.kernel else 0,
    }


def guard_group_to_dict(group: Any, index: int) -> Dict[str, Any]:
    guard = polygon_centroid(group.common_kernel) if group.common_kernel else None
    return {
        "index": index,
        "regions": list(group.region_indices),
        "common_kernel": poly_to_list(group.common_kernel),
        "guard": point_to_list(guard) if guard is not None else None,
        "geom_type": group.geom_type,
        "area": float(group.area),
    }


def build_case(name: str, raw_poly: Iterable[Sequence[Any]], include_strategy: bool = True) -> Dict[str, Any]:
    poly = normalize_polygon(raw_poly)

    reset_merge_stats()
    triangles = triangulate_polygon(poly)
    final_regions, merge_log, convex_regions = optimized_star_decomposition(
        triangles,
        convex_strategy="largest_area",
        star_strategy="largest_area",
        use_convex_premerge=True,
    )
    merge_stats = get_merge_stats()
    guard_groups = greedy_kernel_guard_groups(final_regions, allow_touch=True)

    strategy_comparison = []
    if include_strategy:
        rows = compare_star_strategies(
            triangles,
            strategies=[
                "first",
                "largest_area",
                "largest_kernel",
                "kernel_ratio",
                "fewest_boundary_vertices",
            ],
            use_convex_premerge=True,
        )
        strategy_comparison = [
            {
                "strategy": row[0],
                "convex_pieces": row[1],
                "final_star_regions": row[2],
                "accepted_merges": row[3],
            }
            for row in rows
        ]

    return {
        "name": name,
        "polygon": poly_to_list(poly),
        "input_area": abs(float(polygon_area(poly))),
        "triangles": [poly_to_list(t) for t in triangles],
        "convex_regions": [region_to_dict(r, i) for i, r in enumerate(convex_regions)],
        "final_regions": [region_to_dict(r, i) for i, r in enumerate(final_regions)],
        "guard_groups": [guard_group_to_dict(g, i) for i, g in enumerate(guard_groups)],
        "merge_log": list(merge_log),
        "merge_stats": merge_stats,
        "summary": {
            "vertices": len(poly),
            "initial_triangles": len(triangles),
            "convex_regions": len(convex_regions),
            "final_regions": len(final_regions),
            "guard_groups": len(guard_groups),
            "accepted_merges": len(merge_log),
        },
        "strategy_comparison": strategy_comparison,
    }


def json_response(handler: SimpleHTTPRequestHandler, payload: Dict[str, Any], status: int = 200) -> None:
    data = json.dumps(payload, indent=2).encode("utf-8")
    handler.send_response(status)
    handler.send_header("Content-Type", "application/json; charset=utf-8")
    handler.send_header("Content-Length", str(len(data)))
    handler.end_headers()
    handler.wfile.write(data)


class Handler(SimpleHTTPRequestHandler):
    def do_GET(self) -> None:
        if self.path == "/health":
            json_response(self, {"ok": True, "root": ROOT, "code_dir": CODE_DIR})
            return

        if self.path == "/examples":
            examples = [
                {"name": name, "polygon": poly_to_list(poly)}
                for name, poly in ALL_EXAMPLES.items()
            ]
            json_response(self, {"ok": True, "examples": examples})
            return

        return super().do_GET()

    def do_POST(self) -> None:
        if self.path != "/run_decomposition":
            self.send_error(404)
            return

        try:
            length = int(self.headers.get("Content-Length", "0"))
            body = self.rfile.read(length).decode("utf-8")
            payload = json.loads(body) if body else {}
            poly = payload.get("polygon", [])
            name = payload.get("name", "interactive_custom_polygon")
            include_strategy = bool(payload.get("include_strategy", True))

            result = build_case(name, poly, include_strategy=include_strategy)
            json_response(self, {"ok": True, "case": result})

        except Exception as exc:
            json_response(
                self,
                {
                    "ok": False,
                    "error": str(exc),
                    "traceback": traceback.format_exc(limit=8),
                },
                status=200,
            )


def main() -> None:
    os.chdir(ROOT)
    port = int(os.environ.get("star-shaped decomposition_INTERACTIVE_PORT", "8010"))
    print(f"Code directory: {CODE_DIR}")
    print(f"Serving interactive visualizer at http://localhost:{port}/interactive.html")
    server = ThreadingHTTPServer(("localhost", port), Handler)
    server.serve_forever()


if __name__ == "__main__":
    main()


