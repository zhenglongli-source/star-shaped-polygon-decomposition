"""Visualization utilities for star-shaped decomposition and kernels."""

from typing import List, Optional
import matplotlib.pyplot as plt
from geometry_utils import Polygon
from visualize import polygon_centroid


def closed_xy(poly: Polygon):
    xs = [p[0] for p in poly] + [poly[0][0]]
    ys = [p[1] for p in poly] + [poly[0][1]]
    return xs, ys


def plot_decomposition(
    original: Polygon,
    regions: List[Polygon],
    kernels: Optional[List[Polygon]] = None,
    title: str = "Star-shaped decomposition",
    show_indices: bool = True,
    save_path: Optional[str] = None,
):
    fig, ax = plt.subplots(figsize=(8, 6))

    # Draw original polygon boundary more heavily.
    ox, oy = closed_xy(original)
    ax.plot(ox, oy, linewidth=2.0, marker="o", label="Original polygon")

    if show_indices:
        for i, (x, y) in enumerate(original):
            ax.text(x, y, str(i), fontsize=9, ha="right", va="bottom")

    # Draw final regions.
    for idx, region in enumerate(regions):
        rx, ry = closed_xy(region)
        ax.fill(rx, ry, alpha=0.18)
        ax.plot(rx, ry, linewidth=1.2)
        c = polygon_centroid(region)
        if c is not None:
            ax.text(c[0], c[1], f"R{idx}", fontsize=10, ha="center", va="center")

    # Draw kernels and one representative guard point in each kernel.
    if kernels:
        for idx, kernel in enumerate(kernels):
            if not kernel:
                continue
            kx, ky = closed_xy(kernel)
            ax.fill(kx, ky, alpha=0.35, hatch="//", label="Kernel" if idx == 0 else None)
            kc = polygon_centroid(kernel)
            if kc is not None:
                ax.plot([kc[0]], [kc[1]], marker="x", markersize=10)
                ax.text(kc[0], kc[1], f"g{idx}", fontsize=9, ha="left", va="bottom")

    ax.set_aspect("equal")
    ax.set_title(title)
    ax.legend(loc="best")
    ax.grid(True, alpha=0.25)

    if save_path:
        plt.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.show()


