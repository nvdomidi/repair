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
- Use `repair.FillHoleLiepa` to fill the loops found. You can choose the triangle weight calculation method to be either "area" or "angle", based on papers by: Barequet & Sharir [[1]](#barq) and Peter Liepa [[2]](#liepa), respectively.

# 4. Removing Self Intersections

The objective is to develop an algorithm that given the vertex coordinates of a pair of triangles, can determine whether they intersect or not.

## 4.1. Intersection Configurations

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

<p align="center">
<img src = "/pics/2dintersections.png">
<br>
<span name = "fig1">*Figure 1: Different configurations for 2D triangle-triangle intersection.*</span>
<br>
<br>
</p>

![](/pics/3dintersections.png)
<br>
<span name = "fig2">*Figure 2: Different configurations for 3D triangle-triangle intersection.*</span>
<br>
<br>

## 4.2. 2D Triangle-triangle Intersection

Separating Axis Theorem (SAT) is a method used to determine whether two convex objects are intersecting. The idea is that if there exists an axis for which the intervals of projection of the two convex objects onto it do not intersect, then the objects do not intersect. The axis - if it exists - is known as a separating axis. In order to search for a separating axis, and consequently determine that the objects do not intersect, we only need to test edge normals for both objects. The reason for this is that it can be proven that if there exists at least one separating axis, then there also exists an edge normal that is a separating axis [[3]](#geo). [Figure 3](#fig3) shows how non-intersecting shapes can be separated along the found axis.

![](/pics/sat.png)
<br>
<span name = "fig3">*Figure 3: Left: Separation of convex polygons along the axis - Right: No possible separation [[3]](#geo)*. </span>
<br>

A straightforward approach to implementing an algorithm for this method would be to consider each edge normal as a candidate axis, project every vertex onto it (simply by calculating the dot product of normal vector and vertex), find the maximum projection value from the left polygon and the minimum from the right one, compare them and return true only if there is a separation. Function `Test2DIntersection` in `intersect.go` implements this method.

## 4.3. 3D Triangle-triangle Intersection

An easy approach to solving this problem would be to reduce it to edge-face intersection checks between each edge of one triangle and the face of the other. However, Moller and Haines [[4]](#render) propose a method which is also mentioned in [[3]](#geo). In this method they first reject pairs of triangles whose vertices exist on one side of each other’s planes. This is done by computing the signed distance of each point to the other triangle’s plane.

Given three points in 3D space $v_1,v_2,v_3$, we must first find the plane equation for this triangle. To represent a plane, only the normal vector and a point on the plane is needed. The normal vector $\hat{n}$ of the plane can be calculated as the normalized cross product of vectors formed by subtracting $v_2$ from $v_1$ and $v_3$ from $v_1$:
$$\hat{n} = \frac{v_1 - v_2}{|| v_1 - v_2 ||} \times \frac{v_1 - v_3}{|| v_1 - v_3 ||}$$

And the point can be simply chosen to be $v_1$. The signed distance from a point $x_0$ to the plane is given by $D = \hat{n} \cdot (x_0 - v_1)$ [[5]](#mathworld).

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

$$ t_{0, i}=V_{0, i}^{\prime}+ (V_{0,2}^{\prime}-V_{0, i}^{\prime}) \frac{dist_{V_{0,i}}}{dist_{V_{0,i}} - dist_{V_{0,2}}}, i \in\{0,1\} $$

By comparing $t_{0,0}$ and $t_{0,1}$ to $t_{1,0}$ and $t_{1,0}$, we can find if there is an intersection. [Figure 5](#fig5) shows examples of overlapping and non-overlapping intervals.

![](/pics/intervals.png)
<br>
<span name = "fig5">*Figure 5: Left: Overlapping intervals. - Right: Non-overlapping intervals. [[3]](#geo)* </span>
<br>

The 3D intersection check is implemented inside the `TestIntersection` function.

## 4.4. Intersection Edge Cases

We consider cases (5-8) to be 2D edge cases and cases (12-17) to be 3D edge cases. That means triangle vertices are touching the other triangle (getting very close to it) without penetrating it. In a valid mesh, adjacent triangles share vertices and edges. That means cases 6, 7, 15 and 16 can exist and should not be identified as intersections. The other edge cases will be considered to be intersections. So our goal is to identify whether two triangles are in an edge case configuration and if so, are they considered to be intersections or not.

For 2D configurations, after computing the projection intervals onto the separation axis, these intervals can be compared to each other. If these intervals are touching each other – that means minimum of one interval is very close to maximum of the other interval – there is a high probability (but it’s not guaranteed) that the configuration is an edge case. This can be further confirmed by performing a point-in-polygon test.

For 3D configurations, when testing whether one triangle is on one side of the other triangle’s plane, the signed distances of the vertices are computed. Vertices that have zero distance to the other triangle’s plane may create an edge case. Those vertices, along with the vertices of the other triangle can used in a 2D point-in-polygon test to find out whether they are touching the other triangle or not.

## 4.5. Triangle-triangle Intersection Algorithm

The entire triangle-triangle intersection pipeline is explained below. You can read the ``TestIntersection`` function for details.


1. Check if triangles are degenerate. If so, exist the program.
2. Compute the plane $P_0$ for the first triangle $T_0$.
3. Order the vertices of both triangles counter-clockwise rotation.  
4. Check if vertices from $T_1$ are on one side of $T_0$.
    * Calculate signed distances of each vertex from $T_1$ to plane $P_0$.
    * Store vertices with positive, negative and zero distances separately.
    * If there is at least one positive distance vertex and at least one negative, that means the triangle exists on both sides of the plane, so it is not on one side.
    * If there are three zero distance vertices, that means the triangle exists on the other triangle’s plane, so not on one side.
    * Return whether the vertices are on one side + positive distance vertices + negative + zeros.
5. If there are less than three zero distance vertices, that means the triangles are in 3D space, so perform a 3D edge case test:
    * Set edgeCase to False and valid to True.
    * If any of the vertices from $T_1$ are inside $T_0$, they must be equal to a vertex from $T_0$.
        * For each vertex from $T_1$, rotate it along with all vertices from $T_0$ to be parallel to XY plane.
        * Drop Z values and perform point-in-polygon test. Point-in-polygon is performed by first checking whether the point is almost equal to a vertex from the polygon. Then it checks whether it lies on a horizontal edge (constant Y value). Then it checks whether it lies on a non-horizontal edge. At last, it performs the check which is similar to shooting a ray to the right of the point, counting the number of times it intersects with an edge of the polygon. If the number is even, it is outside, if odd, it is inside. Lastly, return true if any of these conditions are satisfied.
        * If a point is inside the polygon, then edgeCase is set to true.
        * If a point is inside the polygon, it must be equal to one of the vertices of that polygon. If this condition is false, it means it is lying somewhere on the edge or inside the area of that triangle. Therefore, set valid to false.
    * Check if the edges between the positive distance and negative distance vertices are intersecting with the edges of $T_0$. If any of them are intersecting, then it is an edge case and it is not valid.
6. If there is an edge case and it is not valid, we assume there is an intersection. If there is an edge case and it is valid, we assume the triangles do not intersect.
7. If vertices from $T_1$ are not on one side of $P_0$, they must be ordered such that the two vertices on one side are separate from the vertex on the other side (Important for the input to 3D intersection test).
    * If two positive distance vertices and one negative, first return positives then negative.
    * If two negative distance vertices and one positive, first return negatives then positive.
    * If one positive, one negative and one zero, first return negative, then zero, then positive (not exactly like the previous two, but it still works).
8. If all vertices from $T_1$ are on one side of $T_0$, and the triangles have not been identified as an edge case, return “no intersection”.
9. If there are three vertices on the plane of the other triangle, perform 2D separating axis intersection check. If the intervals are touching each other, consider it an edge case, and check if it is a valid edge case.
    * For each vertex of $T_0$, perform point-in-polygon test against $T_1$. If it is inside the polygon, it must be equal to a vertex from $T_1$. If it is inside but not equal, then it is not valid.
    * Repeat the same for $T_1$.
10. Handle 2D intersection as such:
    * If the projected intervals are overlapping, there is an intersection in 2D.
    * If the intervals are not overlapping, but touching each other (edge case), then test if it is an acceptable edge case. If acceptable = no intersection, else = intersection.
    * If there exists a separating axis, then there is no intersection.
11. Find plane $P_1$ for triangle $T_1$. Repeat steps 5, 6, 7 and 8 for vertices of $T_0$ against $P_1$ and its triangle $T_1$.
12. If the algorithm reaches this part, the triangles are in 3D, they are not an edge case and they do not exist on one side of each other’s planes. Compute the intersection line of the two planes, and perform the 3D intersection test.
    * Calculate the intervals that the triangles intersect with the line ($t_{00}$ to $t_{01}$ and $t_{10}$ to $t_{11}$). This is done according to the equations in the previous report.
    * If one of the intervals has close to zero length and is almost equal to one end of the other interval, consider not intersecting (this is a special case that arises when testing many triangles).
    * If the intervals are nearly touching return no intersection (we previously handled edge cases).
    * If the intervals are overlapping, there is an intersection.
    * Otherwise, there is no intersection.

# 5. Octree Structure

To avoid the inefficiency of checking each triangle in a mesh against every other triangle (e.g., 89,075,844 checks for a 9,438-triangle mesh), an octree structure is employed. This approach reduces the number of intersection checks significantly (e.g., to 626,066 for the same mesh).

The octree divides 3D space into a box-shaped boundary. Each node is either a "leaf" node that contains objects, or it is an "inner" node that contains 8 smaller octrees. The objects within each node will be mesh triangles. You can read the implementation in ``octree.go``.

**Implementation steps:**

1. Creating the octree:
    * Initialize the root node as a "leaf" with a specified 3D boundary, maximum recursion level, and maximum objects per node.
    * Determine the boundary for each triangle and insert it into the octree.
2. Inserting objects:
    * For each object to be added, find "leaf" nodes corresponding to the boundary. Insert the triangles into these nodes. Upon reaching maximum capacity, subdivide the node into 8 smaller cubes and recursively insert triangles into the correct child nodes.
3. Retreiving objects:
    * Recursively search for "leaf" nodes within a specified boundary.
    * Return all objects within these nodes, facilitating efficient intersection checks with nearby triangles.

# 6. Hole-filling

The automatic hole-filling algorithm of Peter Liepa [[2]](#liepa) was employed for this part. The file ``liepa.go`` contains this code and separate repo can be found [here](https://github.com/nvdomidi/HoleFillingGo). The hole-filling algorithm depends on ``holedetector.go`` to find the boundary loops. The original paper contains the complete explanation of the algorithm, but an overview is also given here:

This algorithm is based on using dynamic programming and it computes "minimum-weights" for series of ordered vertices iteratively. Note: **_weights_** are assigned to potential triangles and **_minumum-weights_** are calculated for series of vertices. The weights look something like $\Omega(v_i,v_k,v_m)$ - determining how good a potential triangle would be - and minimum-weights look something like $W_{i,i+m}$ - starting from $m=1$ and eventually calculating $W_{0,n-1}$. When all the minimum-weights are calculated and stored, they can be used to retreive the best triangulation of the hole.

The weight function $\Omega$ can be chosen to be either: 
1. The area of the triangle (which is what Barequet and Sharir [[1]](#barq) use in their method).
2. An ordered pair of the largest dihedral angle of the triangle with adjacent triangles + the area of the triangle (what Peter Liepa [[2]](#liepa) uses).

## 6.1. The Hole-filling Algorithm

The hole-filling algorithm is performed as such:
1. Minimum-weights between adjacent vertices are set to be zero. Minimum-weights between vertices with 2 distance are set to be the weight of the triangle between them.
2. Find the minimum-weights for series of vertices with 3 distance, then 4, 5, ... until the minimum-weight for the entire hole is found. To find $W_{i,k}$, change $m$ between $i$ (the starting vertex) and $k$ (the last vertex), calculate $W_{i,m}$, $W_{m,k}$ and the weight of the triangle between these three $\Omega(v_i, v_m, v_k)$, and add these together. The index $m$ at which the minimum value for $W_{i,k}$ is found will be saved as $\lambda_{i,k}$.
3. Repeat step 2 until $W_{0,n-1}$ is found.
4. Recover the triangulation using the $\lambda_{i,k}$ values saved. The triangulation is performed by calling the _"Trace"_ method described by Barequet and Sharir [[1]](#barq).

## 6.2. Recovering the Triangles (Trace Algorithm)

The algorithm below is from Barequet and Sharir [[1]](#barq). The notation is modified to match the previous section.

Let $S := 0$. Invoke the recursive function _"Trace"_ with the parameters $(0,n-1)$.

Function Trace $(i,k)$:

if i + 2 = k then $S := S \cup \Delta v_iv_ov_k$

else do:

1. let $o := \lambda_{i,k}$
2. if $o \not ={i+1}$ then Trace $(i,o)$
3. $S := S \cup \Delta v_iv_ov_k$
4. if $o \not ={k-1}$ then Trace $(o,k)$

od

This recursive method uses the minimum-weights and the indices calculated to perform triangulation.





    


# References

1. <span name = "barq"> M. S. Gill Barequet, "Filling Gaps in the Boundary of a Polyhedron," Computer-Aided Geometric Design, pp. 207-229, 1995. </span>
2. <span name = "liepa"> P. Liepa, "Filling Holes in Meshes," Eurographics Symposium on Geometry Processing, 2003. </span>
3. <span name = "geo">D. E. Philip J. Schneider, "Geometric Tools for Computer Graphics," Elsevier Science Inc., 2002, pp. 265 - 271, 539 - 542, 347 - 376, 529 - 531. </span>
4. <span name = "render">E. H. Tomas Akenine-Moller, Real-Time Rendering, Fourth Edition. </span>
5. <span name = "mathworld"> W. MathWorld, "Point-Plane Distance," [Online]. Available: https://mathworld.wolfram.com/Point-PlaneDistance.html. </span>
