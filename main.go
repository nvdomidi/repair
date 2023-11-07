package main

import (
	"fmt"
	"main/geom"
	"main/io"
	"main/repair"
)

func main() {

	//m := io.ReadBinarySTL("stls/teapot_hole.stl")

	//io.WriteToOBJ(m, "objs/teapot_hole.obj")

	m, _ := io.ReadObj("objs/teapot_hole.obj")

	//repair.WriteToOBJ(m, "teapot_deleted.obj")

	m = repair.RemoveSelfIntersections(m, false, true)

	io.WriteToOBJ(m, "results/removed_self_intersections.obj")

	vertices := m.V
	faces := []geom.Face{}
	face := geom.Face{}
	for i, idx := range m.T {
		face = append(face, int(idx))
		if (i+1)%3 == 0 {
			faces = append(faces, face)
			face = geom.Face{}
		}
	}

	loops := repair.FindBorders(m)
	fmt.Println("number of loops: ", len(loops))

	newTriangles := []uint32{}

	for _, loop := range loops {

		holeTriangles := repair.FillHoleLiepa(vertices, faces, loop, "angle")

		faces = append(faces, holeTriangles...)

		for i := range faces {
			for j := range faces[i] {
				newTriangles = append(newTriangles, uint32(faces[i][j]))
			}
		}
	}

	m.T = newTriangles

	io.WriteToOBJ(m, "results/hole_filled.obj")

	fmt.Println("succesfully filled hole")

}
