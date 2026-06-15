# Experimental PolyPartition Backend

This folder contains an optional experimental backend for the star-shaped polygon decomposition pipeline.

The default public pipeline uses the pure-Python ear-clipping triangulation backend in `triangulation.py`. This experimental backend explores a stronger early-stage pipeline:

```text
PolyPartition Triangulate_MONO
-> HM-style diagonal deletion / convex pre-merge
-> star-shaped merge pipeline
```

## Files

```text
early_pipeline.py              Experimental early-stage pipeline
hm_diagonal_deletion.py        HM-style diagonal deletion logic
polypartition_adapter.py       Python adapter for the PolyPartition CLI wrapper
main_polypp.py                 Demo entry point for the experimental backend
polypp_cli.cpp                 Small local C++ wrapper around PolyPartition
```

## Third-Party Dependency

This backend uses PolyPartition, an MIT-licensed polygon partitioning and triangulation library by Ivan Fratric and contributors.

The PolyPartition source and license are kept under:

```text
third_party/polypartition/
```

PolyPartition is not part of the original project code. It is included under its own MIT License. The local `polypp_cli.cpp` wrapper and Python adapter code are project-local code.

## Build Notes

The repository does not include a prebuilt executable for the PolyPartition wrapper. To use this backend, build `polypp_cli.cpp` locally together with the PolyPartition source files from `third_party/polypartition/`.

A typical build command may look like:

```powershell
g++ experimental/polypp_cli.cpp third_party/polypartition/src/polypartition.cpp -I third_party/polypartition/src -O2 -o polypp_cli
```

The exact command may vary depending on compiler setup and the local PolyPartition folder layout.

## Notes

This backend is kept separate from the default pipeline to keep the main repository easy to run and reproduce. It is intended as an experimental branch showing how the triangulation and convex pre-merge stage can be strengthened.
