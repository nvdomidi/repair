package repair

import "main/geom"

type pair struct {
	frst, scnd uint32
}

func FindBorders(m geom.Mesh) [][]int {
	fxn := findAdjacencyMatrix(m)
	borders := facetBorders(m, fxn)

	loops := [][]int{}
	for _, border := range borders {
		loop := []int{}

		for i := len(border) - 2; i >= 0; i-- {
			loop = append(loop, int(border[i]))
		}

		loops = append(loops, loop)
	}

	return loops
}

// This function takes a mesh as input and returns an adjacenct matrix as output.
// Adjacency matrix is [][3]uint8 datatype
// Each row corresponds to a triangle
// Each entry corresponds to a triangle's edge.
// The value of each entry shows how many triangles are adjacent to that edge (if == 2 --> not boundary edge)
func findAdjacencyMatrix(m geom.Mesh) [][3]uint32 {

	// used to store how many times an edge appears in the mesh's triangles
	edgeCount := make(map[pair]uint32)

	for i := 0; i < len(m.T)/3; i++ {
		tri := []uint32{m.T[i*3], m.T[i*3+1], m.T[i*3+2]}
		// iterate over  triangle edges and count them in edgeCount
		for j := 0; j < len(tri); j++ {
			a := tri[j]
			b := tri[(j+1)%3]
			// order them such that left one is lower, this is done so (a,b) and (b,a) dont get confused
			if a > b {
				a, b = b, a
			}
			newEdge := pair{frst: a, scnd: b}
			edgeCount[newEdge] = edgeCount[newEdge] + 1
		}
	}

	fxn := make([][3]uint32, len(m.T)/3)
	for i := 0; i < len(m.T)/3; i++ {
		tri := []uint32{m.T[i*3], m.T[i*3+1], m.T[i*3+2]}
		for j := 0; j < len(tri); j++ {
			a := tri[j]
			b := tri[(j+1)%3]
			// order them such that left one is lower, this is done so (a,b) and (b,a) dont get confused
			if a > b {
				a, b = b, a
			}
			newEdge := pair{frst: a, scnd: b}
			fxn[i][j] = edgeCount[newEdge]
		}
	}

	return fxn
}

// void MeshAlgorithm::GetFacetBorders
func facetBorders(m geom.Mesh, fxn [][3]uint32) [][]uint32 {
	var openPointDegree = make(map[uint32]int8)

	var borders [][]uint32
	var nbr uint32
	// collect all boundary edges (unsorted)
	var edges []pair
	c := uint32(len(m.T) / 3)
	for i := uint32(0); i < c; i++ {
		for j := uint32(0); j < 3; j++ {
			nbr = fxn[i][j]
			if nbr != 2 {
				edges = append(edges, pair{frst: m.T[int(i*3+j)], scnd: m.T[int(i*3+((j+1)%3))]})
				openPointDegree[m.T[int(i*3+j)]] = openPointDegree[m.T[int(i*3+j)]] + 1
				openPointDegree[m.T[int(i*3+((j+1)%3))]] = openPointDegree[m.T[int(i*3+((j+1)%3))]] + 1
			}
		}
	}

	if len(edges) == 0 {
		return nil // no borders found (=> solid)
	}

	// search for edges in the unsorted list
	var first, last uint32
	var border []uint32
	first = edges[0].frst
	last = edges[0].scnd
	edges = append(edges[:0], edges[1:]...)
	border = append(border, first)
	border = append(border, last)
	// var saveI uint32
	var isEndOfEdgesAndIsAnotherHolePossible bool

	for len(edges) > 0 {
		isEndOfEdgesAndIsAnotherHolePossible = true
		// get adjacent edge
		for i := uint32(0); i < uint32(len(edges)); i++ {
			// saveI = i
			if edges[i].frst == last {
				last = edges[i].scnd
				border = append(border, last)
				edges = append(edges[:i], edges[i+1:]...)
				i = 0
				// saveI = i
				isEndOfEdgesAndIsAnotherHolePossible = false
				break
			} else if edges[i].scnd == first {
				first = edges[i].frst
				// push_front
				// https://stackoverflow.com/questions/53737435/how-to-prepend-int-to-slice
				border = append([]uint32{first}, border...)
				// border = append(border, first)
				edges = append(edges[:i], edges[i+1:]...)
				i = 0
				// saveI = i
				isEndOfEdgesAndIsAnotherHolePossible = false
				break
			} else if edges[i].scnd == last {
				last = edges[i].frst
				border = append(border, last)
				edges = append(edges[:i], edges[i+1:]...)
				i = 0
				// saveI = i
				isEndOfEdgesAndIsAnotherHolePossible = false
				break
			} else if edges[i].frst == first {
				first = edges[i].scnd
				border = append([]uint32{first}, border...)
				// border = append(border, first)
				edges = append(edges[:i], edges[i+1:]...)
				i = 0
				// saveI = i
				isEndOfEdgesAndIsAnotherHolePossible = false
				break
			}
		}

		// Note: Calling erase on list iterators doesn't force a re-allocation and
		// thus doesn't invalidate the iterator itself, only the referenced object
		// edges[saveI] == edges[len(edges)-1]
		if isEndOfEdgesAndIsAnotherHolePossible || len(edges) == 0 || last == first {
			// no further edge found or closed polyline, respectively
			borders = append(borders, border)
			border = nil
			if len(edges) > 0 {
				// start new boundary
				first = edges[0].frst
				last = edges[0].scnd
				edges = append(edges[:0], edges[1:]...)
				border = append(border, first)
				border = append(border, last)
			}
		}
	}
	// borders = splitBoundaryLoops(borders, openPointDegree)

	// split boundary loops if needed
	return splitBoundaryLoops(borders, openPointDegree)
}

func splitBoundaryLoops(borders [][]uint32, openPointDegree map[uint32]int8) [][]uint32 {
	// var nbr uint32
	// Count the number of open edges for each point
	// var openPointDegree = make(map[uint32]int8)
	// for i := uint32(0); i < fx.TriCount; i++ {
	// 	for j := 0; j < 3; j++ {
	// 		nbr = fx.n[i][j]
	// 		if nbr == math.MaxUint32 {
	// 			openPointDegree[fx.GetVertexID(i, j)] = openPointDegree[fx.GetVertexID(i, j)] + 1
	// 			openPointDegree[fx.GetVertexID(i, (j+1)%3)] = openPointDegree[fx.GetVertexID(i, (j+1)%3)] + 1
	// 		}
	// 	}
	// }

	// go through all boundaries and split them if needed
	var splitBorders [][]uint32
	for i := uint32(0); i < uint32(len(borders)); i++ {
		var split = false
		// for j := uint32(0); j < uint32(len(borders[i])); j++ {
		for _, j := range borders[i] {
			// two (or more) boundaries meet in one non-manifold point
			if openPointDegree[uint32(j)] > 2 {
				split = true
				break
			}
		}
		if !split {
			splitBorders = append(splitBorders, borders[i])
		} else {
			splited := splitBoundaryLoopsNnMnfldPnt(borders[i])
			splitBorders = append(splitBorders, splited...)

		}

	}
	return splitBorders
}

func splitBoundaryLoopsNnMnfldPnt(border []uint32) [][]uint32 {
	var splitBorders [][]uint32
	var ptDegree = make(map[uint32]int8)
	var bound []uint32
	for i := uint32(0); i < uint32(len(border)); i++ {
		ptDegree[border[i]]++
		deg := ptDegree[border[i]] + 1
		if deg > 0 {
			for j := uint32(0); j < uint32(len(bound)); j++ {
				if bound[j] == border[i] {
					var boundLoop []uint32
					boundLoop = append(boundLoop, bound[j:]...)
					boundLoop = append(boundLoop, border[i])
					bound = bound[:j]
					// bound = append(bound[:j], bound[j+1:]...)
					splitBorders = append(splitBorders, boundLoop)
					ptDegree[border[i]]--
					break
				}
			}
		}
		bound = append(bound, border[i])
	}

	return splitBorders
}
