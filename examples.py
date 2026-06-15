"""
Example polygons for the star-shaped decomposition star-shaped decomposition project.

Each polygon is given as a list of vertices in boundary order.
The examples are designed to test different behaviors of the algorithm:

1. C-like polygon:
   A simple non-star-shaped polygon that usually decomposes into a small
   number of star-shaped regions.

2. Notched polygon:
   A rectangle with several inward notches. This tests whether convex
   pre-merging reduces the triangulation before star-shaped merging.

3. Zigzag polygon:
   A monotone zigzag-like polygon with several reflex vertices.

4. Two-room corridor:
   Two large rooms connected by a narrow corridor. This is a good example
   where one global guard region is usually not enough.

5. Comb polygon:
   A long rectangle with many inward teeth. This is a harder example and
   often leaves several final star-shaped regions.

6. Staircase polygon:
   A mostly monotone polygon with many reflex vertices but relatively simple
   global structure.

7. Deep-notch polygon:
   A polygon with several deep alternating notches.

8. Spiral-like polygon:
   A corridor-like polygon that bends inward. This is useful as a difficult
   test case for star-shaped merging.
"""


# ------------------------------------------------------------
# Original-style examples
# ------------------------------------------------------------

# A C-like polygon.
# This is a good medium example: simple, non-convex, and usually not
# globally star-shaped.
c_like_polygon = [
    (0, 0),
    (6, 0),
    (6, 1),
    (2, 1),
    (2, 4),
    (6, 4),
    (6, 5),
    (0, 5),
]


# A rectangle with two top notches.
# This is useful for testing whether the convex pre-merge step improves
# over raw triangulation.
notched_polygon = [
    (0, 0),
    (9, 0),
    (9, 5),
    (7, 5),
    (7, 3),
    (6, 3),
    (6, 5),
    (3, 5),
    (3, 3),
    (2, 3),
    (2, 5),
    (0, 5),
]


# A simple zigzag / staircase-like polygon.
# It is not too hard, but it has enough reflex vertices to make the
# decomposition nontrivial.
zigzag_polygon = [
    (0, 0),
    (10, 0),
    (10, 1),
    (9, 1),
    (9, 2),
    (8, 2),
    (8, 3),
    (7, 3),
    (7, 4),
    (6, 4),
    (6, 5),
    (5, 5),
    (5, 6),
    (0, 6),
]


# ------------------------------------------------------------
# More typical / harder project examples
# ------------------------------------------------------------

# Two large rooms connected by a narrow corridor.
# This is one of the best presentation examples because the final result
# is naturally expected to have more than one guardable region.
two_room_corridor = [
    (0, 0),
    (4, 0),
    (4, 1.5),
    (8, 1.5),
    (8, 0),
    (12, 0),
    (12, 5),
    (8, 5),
    (8, 3.5),
    (4, 3.5),
    (4, 5),
    (0, 5),
]


# Comb-shaped polygon.
# This has many repeated inward notches. It is good for stress-testing
# the merge strategy because many local merges are possible, but global
# star-shapedness is limited.
comb_polygon = [
    (0, 0),
    (14, 0),
    (14, 5),
    (13, 5),
    (13, 1),
    (12, 1),
    (12, 5),
    (11, 5),
    (11, 1),
    (10, 1),
    (10, 5),
    (9, 5),
    (9, 1),
    (8, 1),
    (8, 5),
    (7, 5),
    (7, 1),
    (6, 1),
    (6, 5),
    (5, 5),
    (5, 1),
    (4, 1),
    (4, 5),
    (3, 5),
    (3, 1),
    (2, 1),
    (2, 5),
    (0, 5),
]


# Staircase polygon.
# This is useful for a case with many reflex vertices but a more ordered
# shape than the comb example.
staircase_polygon = [
    (0, 0),
    (12, 0),
    (12, 1),
    (11, 1),
    (11, 2),
    (10, 2),
    (10, 3),
    (9, 3),
    (9, 4),
    (8, 4),
    (8, 5),
    (7, 5),
    (7, 6),
    (6, 6),
    (6, 7),
    (0, 7),
]


# Deep alternating notches.
# This example is useful because the polygon is still easy to understand,
# but the inward notches can block large star-shaped merges.
deep_notch_polygon = [
    (0, 0),
    (12, 0),
    (12, 6),
    (10, 6),
    (10, 2),
    (9, 2),
    (9, 6),
    (7, 6),
    (7, 1),
    (6, 1),
    (6, 6),
    (4, 6),
    (4, 2),
    (3, 2),
    (3, 6),
    (0, 6),
]


# A spiral-like orthogonal polygon.
# This is a harder example. It is useful for showing that the method is
# heuristic and that some shapes naturally require several guardable regions.
spiral_like_polygon = [
    (0, 0),
    (12, 0),
    (12, 12),
    (0, 12),
    (0, 2),
    (10, 2),
    (10, 10),
    (2, 10),
    (2, 4),
    (8, 4),
    (8, 8),
    (4, 8),
    (4, 6),
    (6, 6),
    (6, 5),
    (3, 5),
    (3, 9),
    (9, 9),
    (9, 3),
    (1, 3),
    (1, 11),
    (11, 11),
    (11, 1),
    (0, 1),
]

# ------------------------------------------------------------
# Optional dictionary for convenient testing
# ------------------------------------------------------------

ALL_EXAMPLES = {
    "c_like": c_like_polygon,
    "notched": notched_polygon,
    "zigzag": zigzag_polygon,
    "two_room_corridor": two_room_corridor,
    "comb": comb_polygon,
    "staircase": staircase_polygon,
    "deep_notch": deep_notch_polygon,
    "spiral_like": spiral_like_polygon,
}


# A recommended default for algorithm testing.
# Change this string in main.py if you use ALL_EXAMPLES there.
DEFAULT_EXAMPLE_NAME = "two_room_corridor"
DEFAULT_EXAMPLE = ALL_EXAMPLES[DEFAULT_EXAMPLE_NAME]

