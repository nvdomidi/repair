package main

import (
	"fmt"

	"main/repair"
)

func main() {
	fmt.Println("hello world")
	m := repair.ReadBinarySTL("stls/teapot_hole.stl")

	m.WriteToOBJ("teapot_hole.obj")

	m = repair.RemoveSelfIntersections(m, true)

	m.WriteToOBJ("removed_self_intersections.obj")

	vertices := m.V
	faces := []repair.Face{}
	face := repair.Face{}
	for i, idx := range m.T {
		face = append(face, int(idx))
		if (i+1)%3 == 0 {
			faces = append(faces, face)
			face = repair.Face{}
		}
	}

	var mesh repair.Mesh
	mesh.V = vertices
	t := []uint32{}

	loops := repair.FindBorders(m)

	for _, loop := range loops {

		holeTriangles := repair.FillHoleLiepa(vertices, faces, loop, "area")

		faces = append(faces, holeTriangles...)

		for i := range faces {
			for j := range faces[i] {
				t = append(t, uint32(faces[i][j]))
			}
		}
	}

	mesh.T = t

	mesh.WriteToOBJ("hole_filled.obj")

	fmt.Println("succesfully filled hole")

}
