# Overview

There are some requirements for a 3D mesh to be print-ready. For example, the objects must be watertight (have no holes), have no triangles intersecting with one another, have manifold geometry (every edge belongs to exactly two faces), etc. This project currently contains code for **removing self-intersections** and **filling holes** in a 3D mesh. These two can also be used in tandem, but there are some [known issues](https://github.com/nvdomidi/repair/issues/2#issuecomment-1798056218).

# Using the code

The project structure is like this: 1. main.go