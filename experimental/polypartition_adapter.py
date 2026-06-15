"""
Adapter for the optional PolyPartition command-line backend.

This module assumes that polypp_cli.exe has already been built in the project
root.  The CLI is a thin wrapper around Ivan Fratric's PolyPartition library.

Supported modes:
    triangulate_mono  -> PolyPartition Triangulate_MONO
    convex_hm         -> PolyPartition ConvexPartition_HM

The rest of the star-shaped decomposition pipeline only sees ordinary Python polygons:
    List[List[Tuple[float, float]]]
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from geometry_utils import Point, Polygon, ensure_ccw, polygon_area

EPS = 1e-8


def _clean_polygon(poly: Sequence[Point], eps: float = EPS) -> Polygon:
    """Small local cleaner to avoid importing star_merge from this adapter."""
    cleaned: Polygon = []
    for x, y in poly:
        p = (float(x), float(y))
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

    return ensure_ccw(cleaned)


def _default_exe_path() -> Path:
    return Path(__file__).resolve().with_name("polypp_cli.exe")


def _write_input_polygon(poly: Polygon, path: Path) -> Polygon:
    poly = _clean_polygon(poly)
    if len(poly) < 3 or abs(polygon_area(poly)) < EPS:
        raise ValueError("Input polygon is degenerate or has fewer than 3 vertices.")

    with path.open("w", encoding="utf-8") as f:
        f.write(f"{len(poly)}\n")
        for x, y in poly:
            f.write(f"{x:.17g} {y:.17g}\n")
    return poly


def _read_output_json(path: Path) -> Dict:
    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    if not data.get("ok", False):
        raise RuntimeError(f"PolyPartition returned ok=false in {path}")
    if "pieces" not in data:
        raise RuntimeError(f"PolyPartition output has no pieces field: {path}")
    return data


def call_polypartition(
    poly: Polygon,
    mode: str,
    exe_path: Optional[str | Path] = None,
    keep_temp: bool = False,
) -> List[Polygon]:
    """
    Call polypp_cli.exe and return output pieces as CCW Python polygons.

    Parameters
    ----------
    poly:
        Input simple polygon in boundary order.
    mode:
        "triangulate_mono" or "convex_hm".
    exe_path:
        Optional path to polypp_cli.exe.  Defaults to project-root executable.
    keep_temp:
        If True, keep temporary input/output files under _polypp_tmp for debugging.
    """
    if mode not in {"triangulate_mono", "convex_hm"}:
        raise ValueError("mode must be 'triangulate_mono' or 'convex_hm'.")

    exe = Path(exe_path) if exe_path is not None else _default_exe_path()
    if not exe.exists():
        raise FileNotFoundError(
            f"Cannot find {exe}. Build polypp_cli.exe first in the project root."
        )

    if keep_temp:
        tmp_root = Path(__file__).resolve().with_name("_polypp_tmp")
        tmp_root.mkdir(exist_ok=True)
        input_path = tmp_root / f"{mode}_input.txt"
        output_path = tmp_root / f"{mode}_output.json"
        _write_input_polygon(poly, input_path)
        result = subprocess.run(
            [str(exe), mode, str(input_path), str(output_path)],
            cwd=str(Path(__file__).resolve().parent),
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(
                "PolyPartition CLI failed.\n"
                f"returncode={result.returncode}\n"
                f"stdout={result.stdout}\n"
                f"stderr={result.stderr}"
            )
        data = _read_output_json(output_path)
    else:
        with tempfile.TemporaryDirectory(prefix="polypp_") as d:
            tmp = Path(d)
            input_path = tmp / "input.txt"
            output_path = tmp / "output.json"
            _write_input_polygon(poly, input_path)
            result = subprocess.run(
                [str(exe), mode, str(input_path), str(output_path)],
                cwd=str(Path(__file__).resolve().parent),
                capture_output=True,
                text=True,
            )
            if result.returncode != 0:
                raise RuntimeError(
                    "PolyPartition CLI failed.\n"
                    f"returncode={result.returncode}\n"
                    f"stdout={result.stdout}\n"
                    f"stderr={result.stderr}"
                )
            data = _read_output_json(output_path)

    pieces: List[Polygon] = []
    for raw_piece in data["pieces"]:
        piece = _clean_polygon([(float(x), float(y)) for x, y in raw_piece])
        if len(piece) >= 3 and abs(polygon_area(piece)) > EPS:
            pieces.append(piece)
    return pieces


def triangulate_polypartition_mono(
    poly: Polygon,
    exe_path: Optional[str | Path] = None,
) -> List[Polygon]:
    """Return PolyPartition monotone-triangulation pieces."""
    pieces = call_polypartition(poly, "triangulate_mono", exe_path=exe_path)
    bad = [p for p in pieces if len(p) != 3]
    if bad:
        raise RuntimeError(f"Triangulate_MONO returned non-triangle pieces: {bad[:1]}")
    return pieces


def convex_partition_polypartition_hm(
    poly: Polygon,
    exe_path: Optional[str | Path] = None,
) -> List[Polygon]:
    """Return PolyPartition Hertel-Mehlhorn convex-partition pieces."""
    return call_polypartition(poly, "convex_hm", exe_path=exe_path)


def partition_area_error(original: Polygon, pieces: Sequence[Polygon]) -> float:
    """Absolute area mismatch between original polygon and sum of piece areas."""
    original_area = abs(polygon_area(ensure_ccw(list(original))))
    piece_area = sum(abs(polygon_area(p)) for p in pieces)
    return abs(piece_area - original_area)


def summarize_partition(original: Polygon, pieces: Sequence[Polygon]) -> Dict[str, float]:
    original_area = abs(polygon_area(ensure_ccw(list(original))))
    piece_area = sum(abs(polygon_area(p)) for p in pieces)
    return {
        "pieces": len(pieces),
        "original_area": original_area,
        "piece_area_sum": piece_area,
        "area_error": abs(piece_area - original_area),
        "relative_area_error": abs(piece_area - original_area) / original_area if original_area > EPS else 0.0,
    }

