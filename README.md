# Star-Shaped Polygon Decomposition

This repository contains an experimental computational geometry implementation for decomposing a simple polygon into star-shaped pieces. The project was developed as an extension of an Computational Geometry polygon-guarding project.

The implementation includes triangulation, an HM-inspired convex pre-merge stage, greedy star-shaped merging, kernel computation, optional shared-guard grouping, and browser-based visualization tools.

## Main Pipeline

The main pipeline is:

```text
simple polygon
-> triangulation
-> HM-inspired convex pre-merge
-> adjacency-driven star-shaped merge
-> kernel computation
-> optional shared-guard grouping
```

The current implementation focuses on practical algorithm engineering rather than being a production geometry library.

## Repository Contents

```text
main.py                         Main demo script
export_debug_json.py               Exports debug_output.json for the static visualizer
examples.py                     Built-in polygon examples
triangulation.py                   Ear-clipping triangulation
star_merge.py                      Convex pre-merge and star-shaped merge logic
guard_groups.py                    Optional shared-guard grouping
kernel_*.py                        Kernel computation support
visualizer.html                    Static browser visualizer for exported examples
interactive.html                   Interactive browser demo for drawing custom polygons
interactive_server.py              Local backend for the interactive demo
debug_output.json                  Pre-exported example data
assets/                            Example output images
```

## Requirements

Python 3.10+ is recommended.

The core scripts use standard Python plus the geometry dependencies already used in the project environment. If Shapely is available, some fallback geometry operations can use it.

## Run the Main Demo

From the repository root:

```powershell
python main.py
```

This prints the decomposition summary for the main example, including triangulation count, convex pre-merge count, final star-shaped regions, accepted merges, and optional guard sharing.

## Export Debug Data

```powershell
python export_debug_json.py
```

This generates:

```text
debug_output.json
```

which is used by the static visualizer.

## Static Visualizer

Start a local server from the repository root:

```powershell
python -m http.server 8000
```

Open:

```text
http://localhost:8000/visualizer.html
```

This viewer loads the precomputed `debug_output.json` file and displays built-in examples, triangulation, convex pre-merge regions, final star-shaped regions, kernels, and guards.

## Interactive Demo

Start the interactive backend:

```powershell
python interactive_server.py
```

Open:

```text
http://localhost:8010/interactive.html
```

You can draw a simple polygon by clicking vertices in boundary order, load built-in examples, and run the decomposition locally.


## Algorithm Notes

The implementation uses an engineering-oriented pipeline:

1. Triangulate the input simple polygon.
2. Apply an HM-inspired convex pre-merge stage to reduce the number of initial pieces.
3. Run adjacency-driven star-shaped merging.
4. Test merged-region visibility through kernel computation.
5. Optionally group regions that can share a guard through a common kernel.

The merge stage avoids scanning all region pairs globally when possible and focuses on adjacent candidate pairs. The visualization tools expose the intermediate stages so the decomposition can be inspected case by case.

## Optional Experimental Backend

The default pipeline uses a pure-Python ear-clipping triangulation backend for reproducibility. An optional experimental backend is included under `experimental/`, using PolyPartition's MIT-licensed `Triangulate_MONO` routine together with an HM-style diagonal deletion stage.

PolyPartition is included under its own MIT License in `third_party/polypartition/`. The small `polypp_cli.cpp` wrapper and Python adapter code are project-local code.

This optional backend is not required for the main demo, static visualizer, or interactive visualizer.

## Current Limitations

- Input polygons are expected to be simple.
- Degenerate inputs such as repeated points, overlapping edges, or very small nearly collinear edges may fail.
- The default public pipeline uses a compact ear-clipping triangulation backend for reproducibility.
- An optional PolyPartition-based backend is included under experimental/, but it is kept separate from the default pipeline and may require local C++ build setup.
- The implementation is designed for educational and experimental computational geometry use, not for robust industrial CAD/GIS geometry processing.
- This public version focuses on the main reproducible pipeline; exploratory diagnostics are kept out of the repository to keep the project concise.

## Future Improvements

Possible next improvements include:

- Integrate the optional PolyPartition / monotone-triangulation backend more tightly into the main Python pipeline.
- Expand the HM-style convex decomposition stage and compare it systematically with the default pre-merge stage.
- Strengthen DCEL-based data structures for maintaining adjacency and boundary updates.
- Add more systematic benchmark examples.
- Add a polished hosted web demo.

## Project Status

The current version is a cleaned GitHub-ready implementation of the star-shaped decomposition pipeline with local visualization and interactive testing support.
