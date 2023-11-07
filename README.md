# 1. Overview

There are some requirements for a 3D mesh to be print-ready. For example, the objects must be watertight (have no holes), have no triangles intersecting with one another, have manifold geometry (every edge belongs to exactly two faces), etc. This project currently contains code for **removing self-intersections** and **filling holes** in a 3D mesh. These two can also be used in tandem, but there are some [known issues](https://github.com/nvdomidi/repair/issues/2#issuecomment-1798056218).

# 2. Project Structure

The project structure is like this:

- File `main.go` contains code for running repair and calling the functions.
- The repair package contains:
  - File `intersect.go` used for finding self-intersections.
  - File `holedetector.go` (written by [@pouriamsz](https://github.com/pouriamsz)) used to find the holes.
  - File `liepa.go` used for filling the holes (has its own [repo](https://github.com/nvdomidi/HoleFillingGo)).
- The `octree.go` package is used to subdivide the 3D space into an octree structure, used to detect intersections faster. This also has a [known issue](https://github.com/nvdomidi/repair/issues/2#issuecomment-1798056218).
- The `io.go` package is used to work with OBJ and STL files.

# 3. Using the Code

In order to repair a mesh using this project, read main.go and do the following:

- Use either `io.ReadBinarySTL` or `io.ReadObj` and specify the path to your 3D mesh.
- Use `repair.RemoveSelfIntersections` to remove self-intersecting triangles.
  - The first flag `useOctree` specifies whether to use the octree structure or not. Only use if mesh is very large.
  - The second flag `saveRemovedTriangles` can be set to save intersections found in folder "intrsct_obj". Use this to verify whether detected intersections are correct or not.
- Use `repair.FindBorders` to find the hole loops.
- Use `repair.FillHoleLiepa` to fill the loops found. You can choose the triangle weight calculation method to be either "area" or "angle", based on papers by: [Barequet & Sharir](https://dl.acm.org/doi/10.1016/0167-8396%2894%2900011-G) and [Peter Liepa](https://diglib.eg.org/handle/10.2312/SGP.SGP03.200-206), respectively.

# 4. Removing Self Intersections

The objective is to develop an algorithm that given the vertex coordinates of a pair of triangles, can determine whether they intersect or not.

# 4.1. Intersection Configurations

Triangle-triangle intersections in 3D space can have different configurations. If two triangles $T_0$ and $T_1$ exist on planes $P_0$ and $P_1$, the following cases (depending on design choices) can be identified as intersections:

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

[Figure 1](#fig1) shows all the different configurations stated above for 2D intersections (1 to 8) and [Figure 2](#fig2) shows all the configurations for 3D intersections (9 to 17).

![](/pics/2dintersections.png)
<br>
<span name = "fig1">*Figure 1: Different configurations for 2D triangle-triangle intersection.*</span>
<br>
<br>

![](/pics/3dintersections.png)
<br>
<span name = "fig2">*Figure 2: Different configurations for 3D triangle-triangle intersection.*</span>
<br>
<br>

# 4.2. 2D Triangle-triangle Intersection

Separating Axis Theorem (SAT) is a method used to determine whether two convex objects are intersecting. The idea is that if there exists an axis for which the intervals of projection of the two convex objects onto it do not intersect, then the objects do not intersect. The axis - if it exists - is known as a separating axis. In order to search for a separating axis, and consequently determine that the objects do not intersect, we only need to test edge normals for both objects. The reason for this is that it can be proven that if there exists at least one separating axis, then there also exists an edge normal that is a separating axis [[1]](#geo). [Figure 3](#fig3) shows how non-intersecting shapes can be separated along the found axis.

![](/pics/sat.png)
<br>
<span name = "fig3">*Figure 3: Left: Separation of convex polygons along the axis - Right: No possible separation [[1]](#geo)*. </span>
<br>

A straightforward approach to implementing an algorithm for this method would be to consider each edge normal as a candidate axis, project every vertex onto it (simply by calculating the dot product of normal vector and vertex), find the maximum projection value from the left polygon and the minimum from the right one, compare them and return true only if there is a separation. Function `Test2DIntersection` in `intersect.go` implements this method.

# 4.3. 3D Triangle-triangle Intersection

An easy approach to solving this problem would be to reduce it to edge-face intersection checks between each edge of one triangle and the face of the other. However, Moller and Haines [[2]](#render) propose a method which is also mentioned in [[1]](#geo). In this method they first reject pairs of triangles whose vertices exist on one side of each other’s planes. This is done by computing the signed distance of each point to the other triangle’s plane.

Given three points in 3D space $v_1,v_2,v_3$, we must first find the plane equation for this triangle. To represent a plane, only the normal vector and a point on the plane is needed. The normal vector $\hat{n}$ of the plane can be calculated as the normalized cross product of vectors formed by subtracting $v_2$ from $v_1$ and $v_3$ from $v_1$:
$$\hat{n} = \frac{v_1 - v_2}{|| v_1 - v_2 ||} \times \frac{v_1 - v_3}{|| v_1 - v_3 ||}$$

And the point can be simply chosen to be $v_1$. The signed distance from a point $x_0$ to the plane is given by $D = \hat{n} \cdot (x_0 - v_1)$ [[3]](#mathworld).

If the signed distances for each vertex of a triangle to the plane of the other triangle all have the same sign, then the vertices of that triangle exist on one side of the other triangle’s plane. Function `TestOnOneSideOfPlane` implements this. [Figure 4](#fig4) shows what is meant by existing on one side of one another’s planes.

![](/pics/oneside.png)
<br>
<span name = "fig4">*Figure 4: Triangle vertices exist on one side of each other’s planes.* </span>
<br>

By rejecting these pair of triangles, we can confirm that the line $L$ at which the planes intersect will also intersect both triangles. The line will be clipped by both triangles into two intervals. If these two intervals overlap, then the triangles intersect; otherwise, they don’t.

In order to find the intersection intervals, we must first find the equation for $L$:
- By calculating the dot product between the normal vector of each plane and the point on it, we can find the Hessian Normal Form for each plane:
$$P_1:\vec{n_1} \cdot P = s_1$$
$$P_2:\vec{n_2} \cdot P = s_2$$
- Then the line equation can be found as:

$$ a = \frac{s_2 \vec{n_1} \cdot \vec{n_2} - s_1 ||\vec{n_2}||^2}{(\vec{n_1} \cdot \vec{n_2})^2 - ||\vec{n_1}||^2||\vec{n_2}||^2}$$

$$ b = \frac{s_1 \vec{n_1} \cdot \vec{n_2} - s_2 ||\vec{n_1}||^2}{(\vec{n_1} \cdot \vec{n_2})^2 - ||\vec{n_1}||^2||\vec{n_2}||^2}$$

$$ L = P + t(\vec{n_1} \times \vec{n_2}) = (a\vec{n_1} + b\vec{n_2}) + t(\vec{n_1} \times \vec{n_2})$$

Function `LineFromTwoPlanes` uses these equations to find the line at which two planes intersect.

Now in order to calculate the intervals, the signed distances previously computed can be put to use. We know that one vertex of each triangle lies on the opposite side of $L$ from the other two. Assume $V_{0,0}, V_{0,1},V_{0,2}$ belong to triangle $0$ and $V_{1,0}, V_{1,1},V_{1,2}$ belong to triangle $1$. First the triangle vertices will be projected onto $L$:

$$ V_{0,i}^{\prime} = \vec{d} \cdot (V_{0,i} - P), i \in \{0,1,2\} $$

Where $\vec{d}$ is the direction vector and equal to $\vec{n_1} \times \vec{n_2}$. Then, the left and right bounds of the interval of intersection for triangle $0$ can be computed as:

$$ t_{0, i}=V_{0, i}^{\prime}+ (V_{0,2}^{\prime}-V_{0, i}^{\prime}) $$

By comparing $t_{0,0}$ and $t_{0,1}$ to $t_{1,0}$ and $t_{1,0}$, we can find if there is an intersection. [Figure 5](#fig5) shows examples of overlapping and non-overlapping intervals.

![](/pics/intervals.png)
<br>
<span name = "fig5">*Figure 5: Left: Overlapping intervals. - Right: Non-overlapping intervals. [[1]](#geo)* </span>
<br>

The 3D intersection check is implemented inside the `TestIntersection` function.


# References

1. <span name = "geo">D. E. Philip J. Schneider, "Geometric Tools for Computer Graphics," Elsevier Science Inc., 2002, pp. 265 - 271, 539 - 542, 347 - 376, 529 - 531. </span>
2. <span name = "render">E. H. Tomas Akenine-Moller, Real-Time Rendering, Fourth Edition. </span>
3. <span name = "mathworld"> W. MathWorld, "Point-Plane Distance," [Online]. Available: https://mathworld.wolfram.com/Point-PlaneDistance.html. </span>
