import matplotlib.pyplot as plt
from geometry_utils import Polygon


def polygon_centroid(poly: Polygon):
    if not poly:
        return None

    area2 = 0.0
    cx = 0.0
    cy = 0.0
    n = len(poly)

    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        cross = x1 * y2 - x2 * y1
        area2 += cross
        cx += (x1 + x2) * cross
        cy += (y1 + y2) * cross

    if abs(area2) < 1e-12:
        # fallback: average of vertices
        xs = [p[0] for p in poly]
        ys = [p[1] for p in poly]
        return (sum(xs) / len(xs), sum(ys) / len(ys))

    cx /= (3.0 * area2)
    cy /= (3.0 * area2)
    return (cx, cy)

def draw_polygon(ax, poly: Polygon, color="blue", label="Polygon"):
    xs = [p[0] for p in poly] + [poly[0][0]]
    ys = [p[1] for p in poly] + [poly[0][1]]
    ax.plot(xs, ys, marker="o", label=label)

    for i, (x, y) in enumerate(poly):
        ax.text(x, y, str(i), fontsize=9, ha="right", va="bottom")

def draw_filled_polygon(ax, poly: Polygon, alpha=0.3, label="Kernel"):
    if not poly:
        return
    xs = [p[0] for p in poly]
    ys = [p[1] for p in poly]
    ax.fill(xs, ys, alpha=alpha, label=label)


def draw_kernel_point(ax, kernel: Polygon):
    if not kernel:
        return
    c = polygon_centroid(kernel)
    if c is not None:
        x, y = c
        ax.plot([x], [y], marker="x", markersize=10, label="Kernel point")


def show_result(poly: Polygon, kernel: Polygon, title: str):
    fig, ax = plt.subplots()
    draw_polygon(ax, poly, label="Polygon")
    if kernel:
        draw_filled_polygon(ax, kernel, alpha=0.35, label="Kernel")
        draw_kernel_point(ax, kernel)
    ax.set_title(title)
    ax.set_aspect("equal")
    ax.legend()
    plt.show()

