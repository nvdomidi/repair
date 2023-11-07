<style type="text/css">
    ol { list-style-type: upper-alpha; }
</style>

# Overview

There are some requirements for a 3D mesh to be print-ready. For example, the objects must be watertight (have no holes), have no triangles intersecting with one another, have manifold geometry (every edge belongs to exactly two faces), etc. This project currently contains code for **removing self-intersections** and **filling holes** in a 3D mesh. These two can also be used in tandem, but there are some [known issues](https://github.com/nvdomidi/repair/issues/2#issuecomment-1798056218).

# Project Structure

The project structure is like this:

- File `main.go` contains code for running repair and calling the functions.
- The repair package contains:
  - File `intersect.go` used for finding self-intersections.
  - File `holedetector.go` (written by [@pouriamsz](https://github.com/pouriamsz)) used to find the holes.
  - File `liepa.go` used for filling the holes (has its own [repo](https://github.com/nvdomidi/HoleFillingGo)).
- The `octree.go` package is used to subdivide the 3D space into an octree structure, used to detect intersections faster. This also has a [known issue](https://github.com/nvdomidi/repair/issues/2#issuecomment-1798056218).
- The `io.go` package is used to work with OBJ and STL files.

# Using the Code

In order to repair a mesh using this project, read main.go and do the following:

- Use either `io.ReadBinarySTL` or `io.ReadObj` and specify the path to your 3D mesh.
- Use `repair.RemoveSelfIntersections` to remove self-intersecting triangles.
  - The first flag `useOctree` specifies whether to use the octree structure or not. Only use if mesh is very large.
  - The second flag `saveRemovedTriangles` can be set to save intersections found in folder "intrsct_obj". Use this to verify whether detected intersections are correct or not.
- Use `repair.FindBorders` to find the hole loops.
- Use `repair.FillHoleLiepa` to fill the loops found. You can choose the triangle weight calculation method to be either "area" or "angle", based on papers by: [Barequet & Sharir](https://dl.acm.org/doi/10.1016/0167-8396%2894%2900011-G) and [Peter Liepa](https://diglib.eg.org/handle/10.2312/SGP.SGP03.200-206), respectively.

# Removing Self Intersections

The objective is to develop an algorithm that given the vertex coordinates of a pair of triangles, can determine whether they intersect or not. Triangle-triangle intersections in 3D space can have different configurations. If two triangles $T_0$ and $T_1$ exist on planes $P_0$ and $P_1$, the following cases (depending on design choices) can be identified as intersections:

<ol> 
<li> Planes are equal, and one vertex from a triangle lies inside the other. </li>
<li> Planes are equal, the triangles share a common edge and their other edges intersect. </li>
<li> Planes are equal, and the triangles have three intersecting edges. </li>
<li> Planes are equal, and the triangles are equal. </li>
<li> Planes are equal, and a vertex from one triangle connects with the edge of the other. </li>
<li> Planes are equal, and the triangles share a common vertex with no other intersection. </li>
<li> Planes are equal, and the triangles share a common edge with no other intersection. </li>
<li> Planes are equal, and the triangles share part of an edge with no other intersection. </li>
<br>
<li> Planes arent equal, and the triangles intersect in a line spanning from one edge of a triangle to another edge of that triangle. </li>
<li> Planes arent equal, and the triangles intersect in a line spanning from one edge of each triangle to another point inside that triangle. </li>
<li> Planes arent equal, and the triangles intersect in a line spanning from one edge of each triangle to another edge of that triangle. </li>
<li> Planes arent equal, and a vertex from one triangle lies inside the area of the other. </li>
<li> Planes arent equal, and a vertex from one triangle connects to the edge of the other. </li>
<li> Planes arent equal, and an edge from each triangle intersects with an edge of the other. </li>
<li> Planes arent equal, and the triangles share a common vertex. </li>
<li> Planes arent equal, and the triangles share a common edge. </li>
<li> Planes arent equal, and the triangles share part of a common edge. </li>
</ol>

Figure 1 shows all the different configurations stated above for 2D intersections (A to H) and Figure 2 shows all the configurations for 3D intersections (I to Q).

![](/pics/2dintersections.png)
*Figure 1: Different configurations for 2D triangle-triangle intersection.*

![](/pics/3dintersections.png)
*Figure 2: Different configurations for 3D triangle-triangle intersection*
