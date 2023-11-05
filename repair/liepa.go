package repair

import (
	"fmt"
	"math"
	"sort"
)

func reverseInts(ints []int) {
	for i := len(ints)/2 - 1; i >= 0; i-- {
		opp := len(ints) - 1 - i
		ints[i], ints[opp] = ints[opp], ints[i]
	}
}

// Function used to find the boundary loops.
// Input: Faces of the mesh: slice of Faces (which is slice of ints)
// Output: Slice of slice of integers = Slice of Loops
func FindBoundaryLoops(faces []Face) [][]int {

	edges := [][2]int{}
	for i := range faces {
		edges = append(edges, [2]int{faces[i][0], faces[i][1]})
	}
	for i := range faces {
		edges = append(edges, [2]int{faces[i][1], faces[i][2]})
	}
	for i := range faces {
		edges = append(edges, [2]int{faces[i][2], faces[i][0]})
	}

	edgesSorted := make([][2]int, len(edges))

	for i := range edges {
		if edges[i][0] > edges[i][1] {
			edgesSorted[i][0], edgesSorted[i][1] = edges[i][1], edges[i][0]
		} else {
			edgesSorted[i][0], edgesSorted[i][1] = edges[i][0], edges[i][1]
		}
	}

	uniqueEdgesMap := make(map[[2]int]int)
	var uniqueEdges [][2]int

	for _, row := range edgesSorted {

		uniqueEdgesMap[row] = uniqueEdgesMap[row] + 1
	}

	for _, row := range edges {
		var sortedEdge [2]int
		if row[0] > row[1] {
			sortedEdge = [2]int{row[1], row[0]}
		} else {
			sortedEdge = [2]int{row[0], row[1]}
		}
		if uniqueEdgesMap[sortedEdge] == 1 {
			uniqueEdges = append(uniqueEdges, row)
		}
	}

	boundaryMap := make(map[int]int)

	for _, row := range uniqueEdges {
		boundaryMap[row[0]] = row[1]
	}

	boundaryLoops := [][]int{}
	boundaryLoop := []int{}
	vertex := -1

	for {
		if vertex == -1 {
			if len(boundaryMap) == 0 {
				break
			}
			for key, _ := range boundaryMap {
				vertex = key
				break
			}
			boundaryLoop = []int{vertex}
		} else {
			nextVertex := boundaryMap[vertex]
			delete(boundaryMap, vertex)
			if nextVertex == boundaryLoop[0] {
				reverseInts(boundaryLoop)
				boundaryLoops = append(boundaryLoops, boundaryLoop)
				vertex = -1
			} else {
				boundaryLoop = append(boundaryLoop, nextVertex)
				vertex = nextVertex

			}
		}
	}

	return boundaryLoops

}

func computeTriangleArea(v []Vertex, i, j, k int) float32 {
	AB := Add(v[j], Multiply(v[i], -1))
	AC := Add(v[k], Multiply(v[i], -1))

	area := Length(CrossProduct(AB, AC)) * 0.5
	return area
}

func computeTriangleNormal(v []Vertex, i, j, k int) Vertex {
	AB := Add(v[j], Multiply(v[i], -1))
	AC := Add(v[k], Multiply(v[i], -1))

	return Normalize(CrossProduct(AB, AC))
}

func cycle3Origins(b_face Face, n int) [2]int {
	sort.Ints(b_face)
	i := b_face[0]
	j := b_face[1]
	k := b_face[2]

	if i == -1 {
		if j == 0 && k == n-1 {
			return [2]int{n - 1, -1}
		}
		if j+1 == k {
			return [2]int{j, -1}
		}
		return [2]int{-1, -1}
	}
	if i == 0 && k == n-1 {
		if j == 1 {
			return [2]int{n - 1, 0}
		}
		if j == n-2 {
			return [2]int{n - 2, n - 1}
		}
	}
	return [2]int{i, j}

}

// Function used to fill holes of a triangular mesh.
// Input params: vertices and faces of original mesh. boundaryLoop = a slice of vertex numbers. method = "area" or "angle"
// Output: returns slice of new faces that fill the hole
func FillHoleLiepa(vertices []Vertex, faces []Face, boundaryLoop []int, method string) []Face {
	var holeTriangles []Face
	n := len(boundaryLoop)

	// index 0 = for adjacent vertices. 1 = for vertices with 2 distance. 2 = for vertices with 3 distance. ... n-2 = for between first and last vertex.
	areas := make([][]float32, n-1)
	lambdas := make([][]int, n-1)

	// initialize everything to zero
	for i := n - 1; i > 0; i-- {
		// zero vector with size of i
		zeroVector := make([]float32, i)
		for it := range zeroVector {
			zeroVector[it] = 0
		}

		areas[n-1-i] = zeroVector

		if i < n-2 {
			// zero vector with size of i
			zeroVector := make([]int, i)
			for it := range zeroVector {
				zeroVector[it] = 0
			}
			lambdas[n-1-i] = zeroVector
		} else {
			lambdas[n-1-i] = []int{}
		}
	}

	// for vertices with 2 distance, just calculate the triangle area
	for i := 0; i < n-2; i++ {
		areas[1][i] = computeTriangleArea(vertices, boundaryLoop[i], boundaryLoop[i+1], boundaryLoop[i+2])
	}

	if method == "area" {
		fmt.Println("Area-based approach à la Barequet and Sharir. Areas are used as weights.")
		// Area-based approach à la Barequet and Sharir. Areas are used as weights.
		for j := 3; j < n; j++ {
			for i := 0; i < n-j; i++ {
				min_area := float32(math.MaxFloat32)
				optimal_m := -1
				for m := 0; m < j-1; m++ {
					m1 := j - m - 2
					i1 := i + 1 + m
					area := areas[m][i] + areas[m1][i1]
					area += computeTriangleArea(vertices, boundaryLoop[i], boundaryLoop[i1], boundaryLoop[i+j])
					if area < min_area {
						min_area = area
						optimal_m = m
					}
				}
				areas[j-1][i] = min_area
				lambdas[j-1][i] = i + 1 + optimal_m
			}
		}

	} else if method == "angle" {
		fmt.Println("Dihedral-angle-based approach by Liepa. Angle-area pairs are used as weights.")
		// Dihedral-angle-based approach by Liepa. Angle-area pairs are used as weights.

		b := make([]int, len(vertices))
		for i := range b {
			b[i] = -1
		}

		for i := 0; i < len(boundaryLoop); i++ {
			b[boundaryLoop[i]] = i
		}

		b_faces := make([]Face, len(faces))
		for i := range faces {
			f := Face{}
			for j := range faces[i] {
				f = append(f, b[faces[i][j]])
			}
			b_faces[i] = f // faces that are inside the boundary loop look something like this: [22 23 -1] --> that means vertices 22, 23 from boundary loop + something not in boundary loop
		}

		edgeFaceNormals := make([][]Vertex, n-1)
		for i := n - 1; i > 0; i-- {

			var x int
			if i < n-1 {
				x = i
			} else {
				x = n
			}
			edgeFaceNormal := make([]Vertex, x)
			for j := range edgeFaceNormal {
				edgeFaceNormal[j] = Vertex{0, 0, 0}
			}

			edgeFaceNormals[n-1-i] = edgeFaceNormal
		}

		// edgeFaceNormals first element: matrix containing normals for triangles adjacent to each edge of the boundary loop
		for f := 0; f < len(faces); f++ {
			b_face := b_faces[f]

			// check if its loop edge
			num := 0
			for _, i := range b_face {
				if i == -1 {
					num++
				}
			}
			if num < 2 {
				face := faces[f]
				normal := computeTriangleNormal(vertices, face[0], face[1], face[2])
				ij := cycle3Origins(b_face, n)
				if ij[0] != -1 {
					edgeFaceNormals[0][ij[0]] = normal
				}
				if ij[1] != -1 {
					edgeFaceNormals[0][ij[1]] = normal
				}

			}
		}

		dotProducts := make([][]float64, n-1)
		for i := n - 1; i > 0; i-- {
			onesVector := make([]float64, i)
			for j := range onesVector {
				onesVector[j] = 1
			}
			dotProducts[n-1-i] = onesVector
		}

		// edgeFaceNormals second elemenet: triangles inside the loop
		for i := 0; i < n-2; i++ {
			edgeFaceNormals[1][i] = computeTriangleNormal(vertices, boundaryLoop[i], boundaryLoop[i+1], boundaryLoop[i+2])
		}

		// dot products are calculated to compare triangles adjacent to the loop, to the new triangles inside the loop
		for i := 0; i < n-2; i++ {
			dot0 := DotProduct(edgeFaceNormals[1][i], edgeFaceNormals[0][i])
			dot1 := DotProduct(edgeFaceNormals[1][i], edgeFaceNormals[0][i+1])
			dotProducts[1][i] = math.Min(float64(dot0), float64(dot1))
		}

		for j := 3; j < n; j++ {
			for i := 0; i < n-j; i++ {
				max_d := -math.MaxFloat32
				min_area := math.MaxFloat32
				optimal_m := -1
				var optimal_normal Vertex

				for m := 0; m < j-1; m++ {

					m1 := j - m - 2
					i1 := i + 1 + m
					triangle := Face{boundaryLoop[i], boundaryLoop[i1], boundaryLoop[i+j]}
					normal := computeTriangleNormal(vertices, triangle[0], triangle[1], triangle[2])
					// compare for created triangles inside the hole
					d := math.Min(float64(DotProduct(normal, edgeFaceNormals[m][i])), float64(DotProduct(normal, edgeFaceNormals[m1][i1])))
					if i == 0 && j == n-1 {
						// compare to triangle adjacent and outside the hole
						d = math.Min(d, float64(DotProduct(normal, edgeFaceNormals[0][n-1])))
					}
					d = math.Min(d, float64(dotProducts[m][i]))
					d = math.Min(d, float64(dotProducts[m1][i1]))
					area := areas[m][i] + areas[m1][i1] + computeTriangleArea(vertices, triangle[0], triangle[1], triangle[2])
					if max_d < d || (max_d == d && area < float32(min_area)) {
						max_d = d
						min_area = float64(area)
						optimal_m = m
						optimal_normal = normal
					}
				}
				dotProducts[j-1][i] = max_d
				areas[j-1][i] = float32(min_area)
				lambdas[j-1][i] = i + 1 + optimal_m
				edgeFaceNormals[j-1][i] = optimal_normal
			}
		}
	}

	// Reconstruct triangulation
	sections := [][2]int{{0, n - 1}}
	triangles := [][3]int{}

	num := 0

	for len(sections) > 0 {
		num++
		section := sections[len(sections)-1]
		d := section[0]
		b := section[1]
		sections = sections[:len(sections)-1]
		var m int
		if b-d == 2 {
			m = d + 1
		} else {
			m = lambdas[b-d-1][d]
		}
		triangles = append(triangles, [3]int{d, m, b})
		if 1 < m-d {
			sections = append(sections, [2]int{d, m})
		}
		if 1 < b-m {
			sections = append(sections, [2]int{m, b})
		}
	}

	/// add hole triangles ...
	for i := range triangles {
		var holeTriangle Face
		triangle := triangles[i]
		for j := 0; j < 3; j++ {
			holeTriangle = append(holeTriangle, boundaryLoop[triangle[j]])
		}
		holeTriangles = append(holeTriangles, holeTriangle)
	}

	return holeTriangles

}
