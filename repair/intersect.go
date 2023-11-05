package repair

import (
	"fmt"
	"math"
)

// Checks if the three points are colinear
func CheckTriangleDegenerate(v []Vertex) bool {
	if len(v) != 3 {
		return true
	}
	// if cross product of two edges is zero then they are on the same line
	vector1 := Subtract(v[1], v[0])
	vector2 := Subtract(v[2], v[0])
	crossProduct := CrossProduct(vector1, vector2)
	return crossProduct.X == 0 && crossProduct.Y == 0 && crossProduct.Z == 0
}

/* 2D Triangles Intersection */
func signedArea(a, b, c Vertex) float32 {
	return 0.5 * (a.X*(b.Y-c.Y) + b.X*(c.Y-a.Y) + c.X*(a.Y-b.Y))
}

// This function is used to order all vertices in consistent CCW order
func orderVertices(a, b, c Vertex, desiredOrder string) (a1, a2, a3 Vertex) {
	// default is CCW
	if desiredOrder == "" {
		desiredOrder = "CCW"
	}

	if (signedArea(a, b, c) > 0 && desiredOrder == "CW") || (signedArea(a, b, c) < 0 && desiredOrder == "CCW") {
		b, c = c, b
	}

	return a, b, c
}

// This function finds edge normal and returns it as a vertex
// A, B = input, vector normal to BA is returned
func perp(v0, v1 Vertex) Vertex {
	dx := v1.X - v0.X
	dy := v1.Y - v0.Y
	return Vertex{dy, -dx, 0}

}

// This function is used to project vertices of a triangle onto an axis and find minimum and maximum of projections
// Input: Vertex slice representing triangle in CCW order, d which is of type vertex and represents the axis
// Output: min value and max value
func computeInterval(c []Vertex, d Vertex) (float32, float32) {
	min := DotProduct(d, c[0])
	max := min

	for i := 1; i < len(c); i++ {
		value := DotProduct(d, c[i])
		if value < min {
			min = value
		} else if value > max {
			max = value
		}
	}

	return min, max
}

/* SAT Direct Implementation */
// first return: intersecting or not
// second return: is it edge case?
func Test2DIntersection(a1, b1, c1, a2, b2, c2 Vertex) (bool, bool) {
	triangle0 := []Vertex{a1, b1, c1}
	triangle1 := []Vertex{a2, b2, c2}
	epsilon := float32(1e-5)

	// test edge normals of triangle0 for separation
	i1 := len(triangle0) - 1
	for i0 := 0; i0 < len(triangle0); i0++ {
		D := perp(triangle0[i0], triangle0[i1])
		min0, max0 := computeInterval(triangle0, D)
		min1, max1 := computeInterval(triangle1, D)

		//fmt.Println("min0:", min0, "max0:", max0, "min1:", min1, "max1:", max1)

		if math.Abs(float64(max0-min1)) < float64(epsilon) || math.Abs(float64(max1-min0)) < float64(epsilon) {
			//fmt.Println("Edge case: very close")
			return false, true
		}

		if max0 < min1 || max1 < min0 {
			return false, false
		}

		i1 = i0
	}

	// test edge normals of triangle1 for separation
	i1 = len(triangle1) - 1
	for i0 := 0; i0 < len(triangle1); i0++ {
		D := perp(triangle1[i0], triangle1[i1])
		min0, max0 := computeInterval(triangle0, D)
		min1, max1 := computeInterval(triangle1, D)

		//fmt.Println("min0:", min0, "max0:", max0, "min1:", min1, "max1:", max1)

		if math.Abs(float64(max0-min1)) < float64(epsilon) || math.Abs(float64(max1-min0)) < float64(epsilon) {
			//fmt.Println("Edge case: very close")
			return false, true
		}

		if max0 < min1 || max1 < min0 {
			return false, false
		}

		i1 = i0
	}

	return true, false

}

// Takes three points as input and returns normal + one point on plane
// used for Hessian Normal Form equation of a plane
func FindPlaneFromPoints(v []Vertex) (p Plane) {
	n := CrossProduct(Subtract(v[0], v[1]), Subtract(v[0], v[2]))
	n = Multiply(n, 1/Length(n)) // normalized
	return Plane{normal: n, p: v[0]}
}

// D= n∙(x_0 - v_1)
func PointToPlaneSignedDistance(v Vertex, p Plane) float32 {
	return DotProduct(p.normal, Add(v, Multiply(p.p, -1)))
}

// Input: three vertices and plane of the other triangle
// Output: bool->Are all vertices on one side of the plane?
// three maps, for positive, negative, zero distances in order
func TestOnOneSideOfPlane(v1, v2, v3 Vertex, p Plane) (bool, map[Vertex]float32, map[Vertex]float32, map[Vertex]float32) {

	epsilon := float32(1e-5)

	distances := []float32{
		PointToPlaneSignedDistance(v1, p),
		PointToPlaneSignedDistance(v2, p),
		PointToPlaneSignedDistance(v3, p),
	}
	vertices := []Vertex{v1, v2, v3}

	positive := make(map[Vertex]float32)
	negative := make(map[Vertex]float32)
	zero := make(map[Vertex]float32)

	for i, d := range distances {
		if d > epsilon {
			positive[vertices[i]] = distances[i]
		} else if d < -epsilon {
			negative[vertices[i]] = distances[i]
		} else {
			zero[vertices[i]] = 0.0
		}
	}

	onOneSide := true
	if (len(positive) > 0 && len(negative) > 0) || len(zero) == 3 {
		onOneSide = false
	}

	return onOneSide, positive, negative, zero

}

// Finds line equation from two planes. Based on page 531 from "Geometry Tools for Computer Graphics"
func LineFromTwoPlanes(p1, p2 Plane) (Line, bool) {
	d := CrossProduct(p1.normal, p2.normal)
	if Length(d) == 0 {
		fmt.Println("Planes are parallel, no line equation found")
		return Line{}, false
	}

	var line Line
	line.dir = d
	var s1, s2, a, b float32
	s1 = DotProduct(p1.normal, p1.p)
	s2 = DotProduct(p2.normal, p2.p)
	n1n2dot := DotProduct(p1.normal, p2.normal)
	n1normsqr := DotProduct(p1.normal, p1.normal)
	n2normsqr := DotProduct(p2.normal, p2.normal)

	a = (s2*n1n2dot - s1*n2normsqr) / (n1n2dot*n1n2dot - n1normsqr*n2normsqr)
	b = (s1*n1n2dot - s2*n2normsqr) / (n1n2dot*n1n2dot - n1normsqr*n2normsqr)

	line.p = Add(Multiply(p1.normal, a), Multiply(p2.normal, b))

	return line, true
}

// This function takes a vertex loop on a plane and rotates them to have similar Z values
// Input: Vertex loop and normal vector of the plane
// Output: Vertex loop with similar Z values
func RotateVertices(vertices []Vertex, normal Vertex) []Vertex {

	var rotatedVertices []Vertex
	// Change this to the normal vector of the target plane you want to match
	zaxis := Vertex{0, 0, 1}

	k := CrossProduct(normal, zaxis)
	if Length(k) < 1e-7 {
		k = Vertex{1, 0, 0}
	}
	k = Normalize(k)

	cosTheta := DotProduct(normal, zaxis) / Length(normal)
	theta := math.Acos(float64(cosTheta))

	for _, vert := range vertices {

		term1 := Multiply(vert, cosTheta)
		term2 := Multiply(CrossProduct(k, vert), float32(math.Sin(theta)))
		term3 := Multiply(k, (DotProduct(k, vert) * (1 - cosTheta)))

		vrot := Add(term1, Add(term2, term3))

		rotatedVertices = append(rotatedVertices, vrot)
	}

	return rotatedVertices

}

// Vertex connectivity: each vertex connected to the next one
// Tests if the point is inside the polygon
// also returns true if its on edge or on a vertex or very close to it
func pointInPolygon(v []Vertex, p Vertex) bool {

	// first we rotate to match XY plane and drop Z value
	plane := FindPlaneFromPoints(v)
	v = RotateVertices(v, plane.normal)
	p = RotateVertices([]Vertex{p}, plane.normal)[0]

	isInside := false

	for i := 0; i < len(v); i++ {
		j := (i + 1) % len(v)

		// Check if point is very close or on a vertex
		if almostEqual(p.X, v[i].X) && almostEqual(p.Y, v[i].Y) {
			return true
		}

		// Check if point is on the horizontal boundary
		if almostEqual(v[i].Y, p.Y) && almostEqual(v[j].Y, p.Y) &&
			(v[i].X < p.X) != (v[j].X < p.X) {
			return true
		}

		// Check if point is on a non-horizontal boundary
		if almostEqual((v[j].X-v[i].X)*(p.Y-v[i].Y)/(v[j].Y-v[i].Y)+v[i].X, p.X) &&
			math.Min(float64(v[i].Y), float64(v[j].Y)) < float64(p.Y) && float64(p.Y) < math.Max(float64(v[i].Y), float64(v[j].Y)) {
			return true
		}

		if ((v[i].Y > p.Y) != (v[j].Y > p.Y)) &&
			(p.X < (v[j].X-v[i].X)*(p.Y-v[i].Y)/(v[j].Y-v[i].Y)+v[i].X) {
			isInside = !isInside
		}
	}

	return isInside
}

// check intersection of segments ab and cd in 3d space
func checkSegmentsIntersection(a, b, c, d Vertex) bool {

	r := Subtract(b, a)
	s := Subtract(d, c)
	q := Subtract(a, c)

	dotqr := DotProduct(q, r)
	dotqs := DotProduct(q, s)
	dotrs := DotProduct(r, s)
	dotrr := DotProduct(r, r)
	dotss := DotProduct(s, s)

	denom := dotrr*dotss - dotrs*dotrs
	numer := dotqs*dotrs - dotqr*dotss

	t := numer / denom
	u := (dotqs + t*dotrs) / dotss

	p0 := Add(a, Multiply(r, t))
	p1 := Add(c, Multiply(s, u))

	onSegment := false
	intersects := false
	if 0 <= t && t <= 1 && 0 <= u && u <= 1 {
		onSegment = true
	}

	if almostEqual(Length(Subtract(p0, p1)), 0.0) {
		intersects = true
	}

	return onSegment && intersects
}

// Input: Two triangles that have one or two points on each other's planes
// First we do point in polygon test -- IF it IS in polygon then:
// If zeroCount == 2 --> It can pass this part if they share an edge
// If zeroCount == 1 --> It's only acceptable if they share a vertex
// if zeroCount == 0 --> Then first two are on one side of other

func checkEdgeCase3D(positive, negative, zero map[Vertex]float32, a, b, c Vertex) (bool, bool) {
	edgeCase := false
	valid := true
	for v, _ := range zero {
		inside := pointInPolygon([]Vertex{a, b, c}, v)
		if inside {
			edgeCase = true
			if !(pointsEqual(v, a) || pointsEqual(v, b) || pointsEqual(v, c)) {
				valid = false
			}
		}
	}

	// now check edges
	if len(positive) == 1 && len(negative) == 1 {
		var pos1, neg1 Vertex
		for v, _ := range positive {
			pos1 = v
			break
		}
		for v, _ := range negative {
			neg1 = v
			break
		}

		if checkSegmentsIntersection(pos1, neg1, a, b) ||
			checkSegmentsIntersection(pos1, neg1, a, c) ||
			checkSegmentsIntersection(pos1, neg1, b, c) {
			edgeCase = true
			valid = false
			fmt.Println("3D segments intersecting")
		}

	} else if len(positive) == 2 && len(negative) == 1 {
		var pos1, pos2, neg1 Vertex
		i := 0
		for v, _ := range positive {
			i++
			if i == 1 {
				pos1 = v
			} else if i == 2 {
				pos2 = v
			} else {
				break
			}
		}
		for v, _ := range negative {
			neg1 = v
			break
		}

		if checkSegmentsIntersection(pos1, neg1, a, b) ||
			checkSegmentsIntersection(pos1, neg1, a, c) ||
			checkSegmentsIntersection(pos1, neg1, b, c) ||
			checkSegmentsIntersection(pos2, neg1, a, b) ||
			checkSegmentsIntersection(pos2, neg1, a, c) ||
			checkSegmentsIntersection(pos2, neg1, b, c) {
			edgeCase = true
			valid = false
			fmt.Println("3D segments intersecting")
		}

	} else if len(positive) == 1 && len(negative) == 2 {
		var neg1, neg2, pos1 Vertex
		i := 0
		for v, _ := range negative {
			i++
			if i == 1 {
				neg1 = v
			} else if i == 2 {
				neg2 = v
			} else {
				break
			}
		}
		for v, _ := range positive {
			pos1 = v
			break
		}

		if checkSegmentsIntersection(pos1, neg1, a, b) ||
			checkSegmentsIntersection(pos1, neg1, a, c) ||
			checkSegmentsIntersection(pos1, neg1, b, c) ||
			checkSegmentsIntersection(pos1, neg2, a, b) ||
			checkSegmentsIntersection(pos1, neg2, a, c) ||
			checkSegmentsIntersection(pos1, neg2, b, c) {
			fmt.Println("3D segments intersecting")
			edgeCase = true
			valid = false
		}
	}

	return edgeCase, valid
}

// This function is used to test whether the 2D edge  case
// If a point from one triangle is inside the other triangle, it can only be
// equal to one of the vertices
func isEdgeCase2DAcceptable(a1, b1, c1, a2, b2, c2 Vertex) bool {
	t1 := []Vertex{a1, b1, c1}
	t2 := []Vertex{a2, b2, c2}

	for _, vert := range t1 {
		if pointInPolygon([]Vertex{a2, b2, c2}, vert) {
			if pointsEqual(vert, a2) || pointsEqual(vert, b2) || pointsEqual(vert, c2) {
				continue
			} else {
				return false
			}
		}
	}

	for _, vert := range t2 {
		if pointInPolygon([]Vertex{a1, b1, c1}, vert) {
			if pointsEqual(vert, a1) || pointsEqual(vert, b1) || pointsEqual(vert, c1) {
				continue
			} else {
				return false
			}
		}
	}

	return true
}

// This function is called when there are vertices from one triangle on both sides of the other triangle.
// Vertices with the same sign of signed distance must be first and second, and the other one must be third
// This order is important later for the 3D intersection test.
func findVertexOrdering(positive, negative, zero map[Vertex]float32) (Vertex, Vertex, Vertex, float32, float32, float32) {
	var v0, v1, v2 Vertex
	var d0, d1, d2 float32
	if len(positive) == 2 && len(negative) == 1 {
		i := 0
		for v, d := range positive {
			if i == 0 {
				v0, d0 = v, d
			} else if i == 1 {
				v1, d1 = v, d
			}
			i++
		}

		for v, d := range negative {
			v2, d2 = v, d
			break
		}
	} else if len(negative) == 2 && len(positive) == 1 {
		i := 0
		for v, d := range negative {
			if i == 0 {
				v0, d0 = v, d
			} else if i == 1 {
				v1, d1 = v, d
			}
			i++
		}

		for v, d := range positive {
			v2, d2 = v, d
			break
		}
	} else if len(negative) == 1 && len(positive) == 1 {
		for v, d := range negative {
			v0, d0 = v, d
			break
		}
		for v, d := range zero {
			v1, d1 = v, d
			break
		}
		for v, d := range positive {
			v2, d2 = v, d
			break
		}
	}

	return v0, v1, v2, d0, d1, d2
}

func TestIntersection(m Mesh, t1, t2 uint32) (bool, error) {

	// STEP 1: Check if triangles are degenerate
	triangle0Degenerate := CheckTriangleDegenerate([]Vertex{m.V[m.T[t1*3]], m.V[m.T[t1*3+1]], m.V[m.T[t1*3+2]]})
	triangle1Degenerate := CheckTriangleDegenerate([]Vertex{m.V[m.T[t2*3]], m.V[m.T[t2*3+1]], m.V[m.T[t2*3+2]]})

	if triangle0Degenerate || triangle1Degenerate {
		return false, fmt.Errorf("error: triangle is degenerate - colinear points")
	}

	// STEP 2: Compute plane equation of T0
	p0 := FindPlaneFromPoints([]Vertex{m.V[m.T[t1*3]], m.V[m.T[t1*3+1]], m.V[m.T[t1*3+2]]})

	// STEP 3: Order the vertices of both triangles in counter-clockwise rotation.
	a1, b1, c1 := orderVertices(m.V[m.T[t1*3]], m.V[m.T[t1*3+1]], m.V[m.T[t1*3+2]], "CCW") //triangle 0
	a2, b2, c2 := orderVertices(m.V[m.T[t2*3]], m.V[m.T[t2*3+1]], m.V[m.T[t2*3+2]], "CCW") //triangle 1

	// STEP 4: Calculate signed distances from T1 to P0 and check if all vertices are on one side.
	t1OnOneSideT0, positive, negative, zero := TestOnOneSideOfPlane(a2, b2, c2, p0)
	var v10, v11, v12 Vertex
	var d10, d11, d12 float32

	// STEPS 5, 6, 7, 8: 3D edge case test. If edge case? lets check if its acceptable (acceptable means no intersection).
	// If not edge case and T1 is on one side of P0, that means no intersection.
	// If not edge case and T1 is not on one side of P0, that means further checks are needed, therefore, lets just order the vertices
	// such that the ones with the same signed distance are v10, v11 and the other is v12.
	if len(zero) < 3 {
		isEdgeCase3D, isAcceptable := checkEdgeCase3D(positive, negative, zero, a1, b1, c1)
		if isEdgeCase3D && !isAcceptable {
			return true, nil
		} else if isEdgeCase3D && isAcceptable {
			return false, nil
		} else if t1OnOneSideT0 {
			return false, nil // no intersection, all vertices of T1 exist on one side of P0
		} else {
			v10, v11, v12, d10, d11, d12 = findVertexOrdering(positive, negative, zero)
		}
	}

	// STEP 9: Perform 2D intersectoin test.
	if len(zero) == 3 {
		rotatedVerticesT0 := RotateVertices([]Vertex{a1, b1, c1}, p0.normal) // rotate vertices from triangle1 to align with XY plane
		rotatedVerticesT1 := RotateVertices([]Vertex{a2, b2, c2}, p0.normal) // rotate vertices from triangle1 to align with XY plane

		a1, b1, c1 = rotatedVerticesT0[0], rotatedVerticesT0[1], rotatedVerticesT0[2]
		a2, b2, c2 = rotatedVerticesT1[0], rotatedVerticesT1[1], rotatedVerticesT1[2]

		intersectingIn2D, edgeCase2D := Test2DIntersection(a1, b1, c1, a2, b2, c2)

		// STEP 10: If all projected intervals overlap = intersecting in 2D = return true
		// If there is one projection where intervals are touching eachother, test edge case and see if its acceptable
		// If edge case is acceptable = no intersection, O.W = intersection
		// If there is a separation in the projected intervals = no intersection
		if intersectingIn2D {
			return true, nil
		} else if edgeCase2D {
			if isEdgeCase2DAcceptable(a1, b1, c1, a2, b2, c2) {
				return false, nil
			} else {
				return true, nil
			}
		} else {
			return false, nil
		}
	}

	// STEP 11: Compute plane equation of T1 and repeat steps 5,6,7,8 for vertices of T0 against P1 and triangle T1.
	p1 := FindPlaneFromPoints([]Vertex{m.V[m.T[t2*3]], m.V[m.T[t2*3+1]], m.V[m.T[t2*3+2]]})

	t0OnOneSideT1, positive, negative, zero := TestOnOneSideOfPlane(a1, b1, c1, p1)
	var v00, v01, v02 Vertex
	var d00, d01, d02 float32

	if len(zero) < 3 {
		isEdgeCase3D, isAcceptable := checkEdgeCase3D(positive, negative, zero, a2, b2, c2)
		if isEdgeCase3D && !isAcceptable {
			return true, nil
		} else if isEdgeCase3D && isAcceptable {
			return false, nil
		} else if t0OnOneSideT1 {
			return false, nil // no intersection, all vertices of T1 exist on one side of P0
		} else {
			v00, v01, v02, d00, d01, d02 = findVertexOrdering(positive, negative, zero)
		}
	}

	// STEP 12: Compute intersection line and perform 3D intersection test
	l, b := LineFromTwoPlanes(p0, p1)
	if b { // this is just a check to see if theres no error finding the line
		v00proj := DotProduct(l.dir, Subtract(v00, l.p))
		v01proj := DotProduct(l.dir, Subtract(v01, l.p))
		v02proj := DotProduct(l.dir, Subtract(v02, l.p))
		v10proj := DotProduct(l.dir, Subtract(v10, l.p))
		v11proj := DotProduct(l.dir, Subtract(v11, l.p))
		v12proj := DotProduct(l.dir, Subtract(v12, l.p))

		t00 := v00proj + (v02proj-v00proj)*(d00/(d00-d02))
		t01 := v01proj + (v02proj-v01proj)*(d01/(d01-d02))
		t10 := v10proj + (v12proj-v10proj)*(d10/(d10-d12))
		t11 := v11proj + (v12proj-v11proj)*(d11/(d11-d12))

		if t00 > t01 {
			t00, t01 = t01, t00
		}
		if t10 > t11 {
			t10, t11 = t11, t10
		}

		// lets check if the interval is zero

		epsilon := float64(1e-7)

		if math.Abs(float64(t00-t01)) < epsilon {
			if math.Abs(float64(t00-t10)) < epsilon || math.Abs(float64(t00-t11)) < epsilon {
				// not intersecting
				return false, nil
			}
		}

		if math.Abs(float64(t10-t11)) < epsilon {
			if math.Abs(float64(t10-t00)) < epsilon || math.Abs(float64(t10-t01)) < epsilon {
				// not intersecting
				return false, nil
			}
		}

		if t00 > t10 {
			temp := t00
			t00 = t10
			t10 = temp

			temp = t01
			t01 = t11
			t11 = temp
		}

		if math.Abs(float64(t10-t01)) < epsilon {
			return false, nil
		} else if t10 < t01 {
			return true, nil
		} else {
			return false, nil
		}
	}

	return false, fmt.Errorf("something went wrong")

}

// This function takes a mesh as input, removes the self-intersections and returns new mesh
func RemoveSelfIntersections(m Mesh, saveRemovedTriangles bool) Mesh {

	octree := OctreeFromMesh(m)
	intersectingIndices := make(map[uint32]bool)

	numIntersections := 0
	for i := 0; i < len(m.T)/3; i++ {
		tri := []Vertex{m.V[m.T[i*3]], m.V[m.T[i*3+1]], m.V[m.T[i*3+2]]}
		_, nearbyTriangleIndices := octree.retrieve(triangleToBounds(tri))
		for j := 0; j < len(nearbyTriangleIndices); j++ {

			idx := nearbyTriangleIndices[j]
			if idx == uint32(i) {
				continue
			}

			isIntersecting, _ := TestIntersection(m, uint32(i), idx)
			if isIntersecting {
				numIntersections++
				// i is intersecting with idx. just remove idx
				intersectingIndices[idx] = true

				if saveRemovedTriangles {
					fmt.Println("triangle: ", i, " intersecting with: ", idx)
					m2 := Mesh{}
					m2verts := []Vertex{}
					m2verts = append(m2verts, tri...)

					m2verts = append(m2verts, m.V[m.T[idx*3]])
					m2verts = append(m2verts, m.V[m.T[idx*3+1]])
					m2verts = append(m2verts, m.V[m.T[idx*3+2]])
					m2.V = m2verts
					m2.T = []uint32{0, 1, 2, 3, 4, 5}
					filename := fmt.Sprintf("intersect/intersect%d.obj", numIntersections)
					m2.WriteToOBJ(filename)
					fmt.Println("saved: ", filename)

					fmt.Println("v0 := Vertex{", m2verts[0].X, ",", m2verts[0].Y, ",", m2verts[0].Z, "}")
					fmt.Println("v1 := Vertex{", m2verts[1].X, ",", m2verts[1].Y, ",", m2verts[1].Z, "}")
					fmt.Println("v2 := Vertex{", m2verts[2].X, ",", m2verts[2].Y, ",", m2verts[2].Z, "}")
					fmt.Println("v3 := Vertex{", m2verts[3].X, ",", m2verts[3].Y, ",", m2verts[3].Z, "}")
					fmt.Println("v4 := Vertex{", m2verts[4].X, ",", m2verts[4].Y, ",", m2verts[4].Z, "}")
					fmt.Println("v5 := Vertex{", m2verts[5].X, ",", m2verts[5].Y, ",", m2verts[5].Z, "}")
				}

			}
		}
	}

	newTriangles := []uint32{}

	for i := 0; i < len(m.T)/3; i++ {
		if intersectingIndices[uint32(i)] {
			fmt.Println("triangle: ", i)
			continue
		} else {
			newTriangles = append(newTriangles, []uint32{m.T[i*3], m.T[i*3+1], m.T[i*3+2]}...)
		}
	}

	fmt.Println(intersectingIndices)
	m.T = newTriangles
	return m
}
