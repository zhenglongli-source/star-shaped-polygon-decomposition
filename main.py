"""
Star-shaped polygon decomposition demo:
Guarding a simple polygon through star-shaped decomposition.

This version has two algorithmic stages:
1. Triangulation baseline.
2. HM-inspired convex pre-merge, followed by star-shaped merging.

Each final star-shaped region contributes a kernel.  One guard placed anywhere
inside each kernel guards the corresponding region.
"""
from guard_groups import greedy_kernel_guard_groups, print_guard_group_summary
from triangulation import triangulate_polygon
from star_merge import optimized_star_decomposition, compare_star_strategies
from visualize_decomposition import plot_decomposition
from examples import notched_polygon, c_like_polygon, zigzag_polygon, two_room_corridor, comb_polygon, staircase_polygon, deep_notch_polygon, spiral_like_polygon 


def run_demo(poly, name: str, star_strategy: str = "largest_area", use_convex_premerge: bool = True, save_path=None):
    print("=" * 70)
    print(name)
    print(f"Number of polygon vertices: {len(poly)}")

    triangles = triangulate_polygon(poly)
    print(f"Initial triangulation pieces: {len(triangles)}")
    safe_name = name.lower().replace(" ", "_").replace("-", "_")

    plot_decomposition(
        original=poly,
        regions=triangles,
        kernels=None,
        title=f"{name}: initial ear-clipping triangulation",
        save_path=f"{safe_name}_stage1_triangulation.png",
    )
    regions, log, convex_regions = optimized_star_decomposition(
        triangles,
        convex_strategy="largest_area",
        star_strategy=star_strategy,
        use_convex_premerge=use_convex_premerge,
    )
    plot_decomposition(
        original=poly,
        regions=[r.polygon for r in convex_regions],
        kernels=[r.kernel for r in convex_regions],
        title=f"{name}: after HM-inspired convex pre-merge",
        save_path=f"{safe_name}_stage2_hm_premerge.png",
    )
    if use_convex_premerge:
        print(f"After HM-inspired convex pre-merge: {len(convex_regions)} pieces")
    else:
        print("HM-inspired convex pre-merge: skipped")

    print(f"Accepted merges: {len(log)}")
    for line in log:
        print("  " + line)

    print(f"Final star-shaped regions: {len(regions)}")
    print(f"Kernel regions before optional sharing: {len(regions)}")

    guard_groups = greedy_kernel_guard_groups(
        regions,
        allow_touch=True,
    )
    print_guard_group_summary(guard_groups)
    print(f"Optional guard count after kernel sharing: {len(guard_groups)}")

    for i, region in enumerate(regions):
        print(
            f"  Region {i}: boundary vertices={len(region.polygon)}, "
            f"kernel vertices={len(region.kernel)}, source triangles={region.source_count}"
        )

    print("\nStrategy comparison with convex pre-merge:")
    rows = compare_star_strategies(
        triangles,
        strategies=["first", "largest_area", "largest_kernel", "kernel_ratio", "fewest_boundary_vertices"],
        use_convex_premerge=True,
    )
    for strategy, convex_count, final_count, merge_count in rows:
        print(
            f"  {strategy:24s} convex pieces={convex_count:2d}, "
            f"final star regions={final_count:2d}, accepted merges={merge_count:2d}"
        )

    plot_decomposition(
        original=poly,
        regions=[r.polygon for r in regions],
        kernels=[r.kernel for r in regions],
        title=(
            f"{name}: {len(regions)} star-shaped regions, "
            f"{len(guard_groups)} optional guard groups"
        ),
        save_path=save_path,
    )


def main():
    # Change this polygon to test a different example.
    # run_demo(c_like_polygon, "C-like polygon", save_path="c_like_result.png")
    # Other examples:
    # run_demo(notched_polygon, "Notched polygon", save_path="star-shaped decomposition_notched_result.png")
    # run_demo(zigzag_polygon, "Zig-zag polygon", save_path="star-shaped decomposition_zigzag_result.png")
    run_demo(two_room_corridor, "Two room corridor", save_path="two_room_corridor_result.png")
    # run_demo(comb_polygon, "Comb polygon", save_path="star-shaped decomposition_comb_polygon_result.png")
    # run_demo(staircase_polygon, "Staircase polygon", save_path="star-shaped decomposition_staircase_polygon_result.png")
    # run_demo(deep_notch_polygon, "Deep notch polygon", save_path="star-shaped decomposition_deep_notch_polygon_result.png")
    # run_demo(spiral_like_polygon, "Spiral like polygon", save_path="star-shaped decomposition_spiral_like_polygon_result.png")

if __name__ == "__main__":
    main()


