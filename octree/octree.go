package octree

import (
	"fmt"
	"main/geom"
	"math"
)

/* subdivide into octree */

// This struct is used to specify the boundaries of space partitioning
type Bounds struct {
	origin geom.Vertex
	x      float32 // width
	y      float32 // height
	z      float32 // length
}

// This function is used to specify the boundaries of each triangle
func TriangleToBounds(triangle []geom.Vertex) Bounds {
	var minX, minY, minZ, maxX, maxY, maxZ float32
	minX, minY, minZ = math.MaxFloat32, math.MaxFloat32, math.MaxFloat32
	maxX, maxY, maxZ = -math.MaxFloat32, -math.MaxFloat32, -math.MaxFloat32

	for _, v := range triangle {
		if v.X < minX {
			minX = v.X
		} else if v.X > maxX {
			maxX = v.X
		}
		if v.Y < minY {
			minY = v.Y
		} else if v.Y > maxY {
			maxY = v.Y
		}
		if v.Z < minZ {
			minZ = v.Z
		} else if v.Z > maxZ {
			maxZ = v.Z
		}
	}

	return Bounds{origin: geom.Vertex{X: minX, Y: minY, Z: minZ}, x: maxX - minX, y: maxY - minY, z: maxZ - minZ}
}

type Node struct {
	nodetype      string          // "leaf" or "inner"
	objectBounds  []Bounds        // triangle bounds
	objects       [][]geom.Vertex // triangles
	objectIndices []uint32        // corresponding triangle indices from the mesh
	tfl           *Octree         // topfrontleft
	tfr           *Octree         // topfrontright
	tbl           *Octree         // topbackleft
	tbr           *Octree         // topbackright
	bfl           *Octree         // bottomfrontleft
	bfr           *Octree         // bottomfrontright
	bbl           *Octree         // bottombackleft
	bbr           *Octree         // bottombackright
}

// top, bottom = +z , -z
// front, back = +y, -y
// left, right = +x, -x
type Octree struct {
	node       Node
	bounds     Bounds
	maxObjects int
	maxLevels  int
	level      int
}

// This function returns child nodes of current node relevant to the bounds input
func (octree *Octree) GetRelevantNodes(limits Bounds) []*Octree {
	if octree.node.nodetype == "leaf" {
		return nil // error: cant get relevant nodes in leaf
	}

	relevantOctrees := []*Octree{}

	midX := octree.bounds.origin.X + octree.bounds.x/2
	midY := octree.bounds.origin.Y + octree.bounds.y/2
	midZ := octree.bounds.origin.Z + octree.bounds.z/2

	isLeft := limits.origin.X <= midX
	isRight := limits.origin.X+limits.x > midX

	isBack := limits.origin.Y <= midY
	isFront := limits.origin.Y+limits.y > midY

	isBot := limits.origin.Z <= midZ
	isTop := limits.origin.Z+limits.z > midZ

	if isLeft {
		if isFront {
			if isTop {
				relevantOctrees = append(relevantOctrees, octree.node.tfl)
			}
			if isBot {
				relevantOctrees = append(relevantOctrees, octree.node.bfl)
			}
		}

		if isBack {
			if isTop {
				relevantOctrees = append(relevantOctrees, octree.node.tbl)
			}
			if isBot {
				relevantOctrees = append(relevantOctrees, octree.node.bbl)
			}
		}
	}

	if isRight {
		if isFront {
			if isTop {
				relevantOctrees = append(relevantOctrees, octree.node.tfr)
			}
			if isBot {
				relevantOctrees = append(relevantOctrees, octree.node.bfr)
			}
		}

		if isBack {
			if isTop {
				relevantOctrees = append(relevantOctrees, octree.node.tbr)
			}
			if isBot {
				relevantOctrees = append(relevantOctrees, octree.node.bbr)
			}
		}
	}

	return relevantOctrees
}

// This function splits the octree into 8 smaller ones and inserts its objects into the children
func (octree *Octree) Split() {
	lvl := octree.level + 1
	x := octree.bounds.x / 2 // width
	y := octree.bounds.y / 2 // height
	z := octree.bounds.z / 2 // length

	if octree.node.nodetype != "leaf" {
		fmt.Println("Cannot split non-leaf")
		return
	}

	// right means x + width/2
	// front means y + height/2
	// top means z + length/2

	octree.node.nodetype = "inner"
	// top front left
	octree.node.tfl = new(Octree)
	octree.node.tfl.node = Node{}
	octree.node.tfl.bounds = Bounds{origin: geom.Vertex{X: octree.bounds.origin.X, Y: octree.bounds.origin.Y + y, Z: octree.bounds.origin.Z + z}, x: x, y: y, z: z}
	octree.node.tfl.maxObjects, octree.node.tfl.maxLevels, octree.node.tfl.level, octree.node.tfl.node.nodetype = octree.maxObjects, octree.maxLevels, lvl, "leaf"
	octree.node.tfl.node.objects, octree.node.tfl.node.objectBounds = [][]geom.Vertex{}, []Bounds{}
	octree.node.tfl.node.objectIndices = []uint32{}
	// top front right
	octree.node.tfr = new(Octree)
	octree.node.tfr.node = Node{}
	octree.node.tfr.bounds = Bounds{origin: geom.Vertex{X: octree.bounds.origin.X + x, Y: octree.bounds.origin.Y + y, Z: octree.bounds.origin.Z + z}, x: x, y: y, z: z}
	octree.node.tfr.maxObjects, octree.node.tfr.maxLevels, octree.node.tfr.level, octree.node.tfr.node.nodetype = octree.maxObjects, octree.maxLevels, lvl, "leaf"
	octree.node.tfr.node.objects, octree.node.tfr.node.objectBounds = [][]geom.Vertex{}, []Bounds{}
	octree.node.tfr.node.objectIndices = []uint32{}
	// top back left
	octree.node.tbl = new(Octree)
	octree.node.tbl.node = Node{}
	octree.node.tbl.bounds = Bounds{origin: geom.Vertex{X: octree.bounds.origin.X, Y: octree.bounds.origin.Y, Z: octree.bounds.origin.Z + z}, x: x, y: y, z: z}
	octree.node.tbl.maxObjects, octree.node.tbl.maxLevels, octree.node.tbl.level, octree.node.tbl.node.nodetype = octree.maxObjects, octree.maxLevels, lvl, "leaf"
	octree.node.tbl.node.objects, octree.node.tbl.node.objectBounds = [][]geom.Vertex{}, []Bounds{}
	octree.node.tbl.node.objectIndices = []uint32{}
	// top back right
	octree.node.tbr = new(Octree)
	octree.node.tbr.node = Node{}
	octree.node.tbr.bounds = Bounds{origin: geom.Vertex{X: octree.bounds.origin.X + x, Y: octree.bounds.origin.Y, Z: octree.bounds.origin.Z + z}, x: x, y: y, z: z}
	octree.node.tbr.maxObjects, octree.node.tbr.maxLevels, octree.node.tbr.level, octree.node.tbr.node.nodetype = octree.maxObjects, octree.maxLevels, lvl, "leaf"
	octree.node.tbr.node.objects, octree.node.tbr.node.objectBounds = [][]geom.Vertex{}, []Bounds{}
	octree.node.tbr.node.objectIndices = []uint32{}
	// bot front left
	octree.node.bfl = new(Octree)
	octree.node.bfl.node = Node{}
	octree.node.bfl.bounds = Bounds{origin: geom.Vertex{X: octree.bounds.origin.X, Y: octree.bounds.origin.Y + y, Z: octree.bounds.origin.Z}, x: x, y: y, z: z}
	octree.node.bfl.maxObjects, octree.node.bfl.maxLevels, octree.node.bfl.level, octree.node.bfl.node.nodetype = octree.maxObjects, octree.maxLevels, lvl, "leaf"
	octree.node.bfl.node.objects, octree.node.bfl.node.objectBounds = [][]geom.Vertex{}, []Bounds{}
	octree.node.bfl.node.objectIndices = []uint32{}
	// bot front right
	octree.node.bfr = new(Octree)
	octree.node.bfr.node = Node{}
	octree.node.bfr.bounds = Bounds{origin: geom.Vertex{X: octree.bounds.origin.X + x, Y: octree.bounds.origin.Y + y, Z: octree.bounds.origin.Z}, x: x, y: y, z: z}
	octree.node.bfr.maxObjects, octree.node.bfr.maxLevels, octree.node.bfr.level, octree.node.bfr.node.nodetype = octree.maxObjects, octree.maxLevels, lvl, "leaf"
	octree.node.bfr.node.objects, octree.node.bfr.node.objectBounds = [][]geom.Vertex{}, []Bounds{}
	octree.node.bfr.node.objectIndices = []uint32{}
	// bot back left
	octree.node.bbl = new(Octree)
	octree.node.bbl.node = Node{}
	octree.node.bbl.bounds = Bounds{origin: geom.Vertex{X: octree.bounds.origin.X, Y: octree.bounds.origin.Y, Z: octree.bounds.origin.Z}, x: x, y: y, z: z}
	octree.node.bbl.maxObjects, octree.node.bbl.maxLevels, octree.node.bbl.level, octree.node.bbl.node.nodetype = octree.maxObjects, octree.maxLevels, lvl, "leaf"
	octree.node.bbl.node.objects, octree.node.bbl.node.objectBounds = [][]geom.Vertex{}, []Bounds{}
	octree.node.bbl.node.objectIndices = []uint32{}
	// bot back right
	octree.node.bbr = new(Octree)
	octree.node.bbr.node = Node{}
	octree.node.bbr.bounds = Bounds{origin: geom.Vertex{X: octree.bounds.origin.X + x, Y: octree.bounds.origin.Y, Z: octree.bounds.origin.Z}, x: x, y: y, z: z}
	octree.node.bbr.maxObjects, octree.node.bbr.maxLevels, octree.node.bbr.level, octree.node.bbr.node.nodetype = octree.maxObjects, octree.maxLevels, lvl, "leaf"
	octree.node.bbr.node.objects, octree.node.bbr.node.objectBounds = [][]geom.Vertex{}, []Bounds{}
	octree.node.bbr.node.objectIndices = []uint32{}

	for i := 0; i < len(octree.node.objects); i++ {
		limits := octree.node.objectBounds[i]
		obj := octree.node.objects[i]
		idx := octree.node.objectIndices[i]

		relevantOctrees := octree.GetRelevantNodes(limits)
		for _, relevantOctree := range relevantOctrees {
			relevantOctree.insert(limits, obj, idx)
		}
	}

	octree.node.objectBounds = []Bounds{}
	octree.node.objects = [][]geom.Vertex{}
	octree.node.objectIndices = []uint32{}

}

// This function recursively finds child leaf nodes relavant to the bounds inserted,
// Then it inserts the object (triangle) along with the object index (triangle index from original mesh)
func (octree *Octree) insert(limits Bounds, obj []geom.Vertex, idx uint32) {
	if octree.node.nodetype != "leaf" {
		relevantOctrees := octree.GetRelevantNodes(limits)
		for _, relevantOctree := range relevantOctrees {
			relevantOctree.insert(limits, obj, idx)
		}
		return
	}

	octree.node.objectBounds = append(octree.node.objectBounds, limits)
	octree.node.objects = append(octree.node.objects, obj)
	octree.node.objectIndices = append(octree.node.objectIndices, idx)

	if len(octree.node.objects) > octree.maxObjects && octree.level < octree.maxLevels {
		octree.Split()
	}

}

// This function takes a bouinds as input and finds the neighboring objects from the relavant leaf node
func (octree *Octree) Retrieve(limits Bounds) ([][]geom.Vertex, []uint32) {
	if octree.node.nodetype == "leaf" {
		return octree.node.objects, octree.node.objectIndices
	}

	relevantTriangles := [][]geom.Vertex{}
	relevantIndices := []uint32{}

	relevantOctrees := octree.GetRelevantNodes(limits)
	for _, relevantOctree := range relevantOctrees {
		subOctreeTriangles, subOctreeIndices := relevantOctree.Retrieve(limits)
		relevantTriangles = append(relevantTriangles, subOctreeTriangles...)
		relevantIndices = append(relevantIndices, subOctreeIndices...)
	}

	return relevantTriangles, relevantIndices
}

// This function inserts all the triangles of a mesh into an octree structure
func OctreeFromMesh(m geom.Mesh) *Octree {
	octree := new(Octree)
	octree.node = Node{}
	octree.node.nodetype = "leaf"
	octree.node.objectBounds = []Bounds{}
	octree.node.objects = [][]geom.Vertex{}
	// finding the bounds
	//octree.bounds = Bounds{origin: geom.Vertex{-20, -20, -20}, x: 40, y: 40, z: 40}
	var minx, maxx, miny, maxy, minz, maxz float32
	for i, v := range m.V {
		if i == 0 || v.X < minx {
			minx = v.X
		}
		if i == 0 || v.X > maxx {
			maxx = v.X
		}
		if i == 0 || v.Y < miny {
			miny = v.Y
		}
		if i == 0 || v.X > maxy {
			maxy = v.Y
		}
		if i == 0 || v.Z < minz {
			minz = v.Z
		}
		if i == 0 || v.Z > maxz {
			maxz = v.Z
		}
	}

	// adding 10 percent to either side
	octree.bounds = Bounds{origin: geom.Vertex{X: minx - (maxx-minx)*0.1, Y: miny - (maxy-miny)*0.1, Z: minz - (maxz-minz)*0.1},
		x: (maxx - minx) * 1.2,
		y: (maxy - miny) * 1.2,
		z: (maxz - minz) * 1.2}

	octree.maxLevels = 4
	octree.maxObjects = 10
	octree.level = 0

	for i := 0; i < len(m.T)/3; i++ {
		triangleVerts := []geom.Vertex{m.V[m.T[i*3]], m.V[m.T[i*3+1]], m.V[m.T[i*3+2]]}
		triangleBounds := TriangleToBounds(triangleVerts)
		octree.insert(triangleBounds, triangleVerts, uint32(i))
	}

	return octree
}
