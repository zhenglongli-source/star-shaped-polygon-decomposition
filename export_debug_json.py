import json
import os

from examples import ALL_EXAMPLES
from triangulation import triangulate_polygon
from star_merge import optimized_star_decomposition, compare_star_strategies
from guard_groups import greedy_kernel_guard_groups
from visualize import polygon_centroid


def point_to_list(p):
    return [float(p[0]), float(p[1])]


def poly_to_list(poly):
    if poly is None:
        return []
    return [point_to_list(p) for p in poly]


def region_to_dict(region, index):
    guard = polygon_centroid(region.kernel) if region.kernel else None
    return {
        "index": index,
        "polygon": poly_to_list(region.polygon),
        "kernel": poly_to_list(region.kernel),
        "guard": point_to_list(guard) if guard is not None else None,
        "source_count": region.source_count,
        "boundary_vertices": len(region.polygon),
        "kernel_vertices": len(region.kernel) if region.kernel else 0,
    }


def guard_group_to_dict(group, index):
    guard = polygon_centroid(group.common_kernel) if group.common_kernel else None
    return {
        "index": index,
        "regions": group.region_indices,
        "common_kernel": poly_to_list(group.common_kernel),
        "guard": point_to_list(guard) if guard is not None else None,
        "geom_type": group.geom_type,
        "area": float(group.area),
    }


def build_case(name, poly):
    triangles = triangulate_polygon(poly)

    final_regions, merge_log, convex_regions = optimized_star_decomposition(
        triangles,
        convex_strategy="largest_area",
        star_strategy="largest_area",
        use_convex_premerge=True,
    )

    guard_groups = greedy_kernel_guard_groups(final_regions, allow_touch=True)

    strategy_rows = compare_star_strategies(
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

    return {
        "name": name,
        "polygon": poly_to_list(poly),
        "triangles": [poly_to_list(t) for t in triangles],
        "convex_regions": [
            region_to_dict(r, i) for i, r in enumerate(convex_regions)
        ],
        "final_regions": [
            region_to_dict(r, i) for i, r in enumerate(final_regions)
        ],
        "guard_groups": [
            guard_group_to_dict(g, i) for i, g in enumerate(guard_groups)
        ],
        "merge_log": merge_log,
        "summary": {
            "vertices": len(poly),
            "initial_triangles": len(triangles),
            "convex_regions": len(convex_regions),
            "final_regions": len(final_regions),
            "guard_groups": len(guard_groups),
            "accepted_merges": len(merge_log),
        },
        "strategy_comparison": [
            {
                "strategy": row[0],
                "convex_pieces": row[1],
                "final_star_regions": row[2],
                "accepted_merges": row[3],
            }
            for row in strategy_rows
        ],
    }


def main():
    cases = []

    for name, poly in ALL_EXAMPLES.items():
        print(f"Exporting {name}...")
        cases.append(build_case(name, poly))

    data = {
        "project": "Computational Geometry Star-Shaped Polygon Decomposition",
        "cases": cases,
    }

    out_path = "debug_output.json"

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)

    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()

