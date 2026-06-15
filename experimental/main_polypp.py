"""
Demo driver for the PolyPartition-MONO + project-local HM early pipeline.

This does not replace main.py yet.
It tests:

    PolyPartition Triangulate_MONO
    -> project-local HM-style diagonal deletion
    -> HM-seeded hybrid star merge
"""

from guard_groups import greedy_kernel_guard_groups, print_guard_group_summary
from examples import (
    c_like_polygon,
    notched_polygon,
    zigzag_polygon,
    two_room_corridor,
    comb_polygon,
    staircase_polygon,
    deep_notch_polygon,
    spiral_like_polygon,
)
from visualize_decomposition import plot_decomposition
from early_pipeline import run_hybrid_pipeline


def run_demo(
    poly,
    name: str,
    early_backend: str = "polypp_mono_hm",
    save_path=None,
):
    print("=" * 70)
    print(name)
    print(f"Early backend: {early_backend}")
    print(f"Number of polygon vertices: {len(poly)}")

    result = run_hybrid_pipeline(
        poly,
        early_backend=early_backend,
        convex_strategy="largest_area",
        star_strategy="largest_area",
    )

    early = result.early
    final_regions = result.final_regions
    hybrid_stats = result.hybrid_stats

    print(f"Initial triangles: {len(early.triangles)}")
    print(f"Early partition area error: {early.area_error:.6g}")
    print(f"Convex regions after early pipeline: {len(early.convex_regions)}")

    if early.hm_stats is not None:
        print(
            "HM diagonal deletion stats: "
            f"initial_faces={early.hm_stats.initial_faces}, "
            f"initial_adjacency_edges={early.hm_stats.initial_adjacency_edges}, "
            f"diagonal_tests={early.hm_stats.diagonal_tests}, "
            f"accepted_deletions={early.hm_stats.accepted_deletions}, "
            f"rejected_nonconvex={early.hm_stats.rejected_nonconvex}, "
            f"final_faces={early.hm_stats.final_faces}"
        )

    print(f"Final star-shaped regions: {len(final_regions)}")
    print(
        "Hybrid merge stats: "
        f"candidates={hybrid_stats.candidates}, "
        f"valid_candidates={hybrid_stats.valid_candidates}, "
        f"accepted={hybrid_stats.accepted_merges}, "
        f"active_fast_accepts={hybrid_stats.active_fast_accepts}, "
        f"exact_record_candidates={hybrid_stats.exact_record_candidates}, "
        f"full_kernel_fallbacks={hybrid_stats.full_kernel_fallbacks}, "
        f"side_constraints={hybrid_stats.side_constraints_computed}"
    )

    print("\nAccepted merge log:")
    for line in result.full_log:
        print("  " + line)

    guard_groups = greedy_kernel_guard_groups(final_regions, allow_touch=True)
    print_guard_group_summary(guard_groups)
    print(f"Optional guard count after kernel sharing: {len(guard_groups)}")

    for i, region in enumerate(final_regions):
        print(
            f"  Region {i}: boundary vertices={len(region.polygon)}, "
            f"kernel vertices={len(region.kernel)}, "
            f"source pieces={region.source_count}"
        )

    safe_name = name.lower().replace(" ", "_").replace("-", "_")

    plot_decomposition(
        original=poly,
        regions=[r.polygon for r in early.convex_regions],
        kernels=[r.kernel for r in early.convex_regions],
        title=f"{name}: {early_backend} convex pre-processing",
        save_path=f"{safe_name}_{early_backend}_convex.png",
    )

    plot_decomposition(
        original=poly,
        regions=[r.polygon for r in final_regions],
        kernels=[r.kernel for r in final_regions],
        title=f"{name}: {len(final_regions)} star-shaped regions via {early_backend}",
        save_path=save_path,
    )


def main():
    # Recommended current demo.
    run_demo(
        two_room_corridor,
        "Two room corridor",
        early_backend="polypp_mono_hm",
        save_path="two_room_corridor_polypp_mono_hm_result.png",
    )

    # Other examples:
    # run_demo(c_like_polygon, "C-like polygon", early_backend="polypp_mono_hm")
    # run_demo(notched_polygon, "Notched polygon", early_backend="polypp_mono_hm")
    # run_demo(zigzag_polygon, "Zig-zag polygon", early_backend="polypp_mono_hm")
    # run_demo(comb_polygon, "Comb polygon", early_backend="polypp_mono_hm")
    # run_demo(staircase_polygon, "Staircase polygon", early_backend="polypp_mono_hm")
    # run_demo(deep_notch_polygon, "Deep notch polygon", early_backend="polypp_mono_hm")
    # run_demo(spiral_like_polygon, "Spiral like polygon", early_backend="polypp_mono_hm")


if __name__ == "__main__":
    main()
