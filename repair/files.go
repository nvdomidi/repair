package repair

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"fmt"
	"main/geom"
	"os"
	"strconv"
	"strings"
)

/* STL */
// VertexMap is used to identify unique vertices and assign them indices
type VertexMap struct {
	vrtToInt map[geom.Vertex]int
	intToVrt map[int]geom.Vertex
}

// Read binary STL file and return as mesh
func ReadBinarySTL(filepath string) geom.Mesh {
	data, err := os.ReadFile(filepath)
	if err != nil {
		panic(err)
	}

	buffer := bytes.NewBuffer(data[80:])
	var numTriangles uint32
	binary.Read(buffer, binary.LittleEndian, &numTriangles)

	fmt.Println("NumTriangles: ", numTriangles)

	triangles := make([][3]geom.Vertex, numTriangles)

	for i := 0; i < int(numTriangles); i++ {
		buffer.Next(12)
		for j := 0; j < 3; j++ {
			binary.Read(buffer, binary.LittleEndian, &triangles[i][j])
		}
		buffer.Next(2)
	}

	vMap := make(map[geom.Vertex]uint32)
	var m geom.Mesh
	for _, triangle := range triangles {
		for _, vertex := range triangle {
			if idx, exists := vMap[vertex]; !exists {
				vMap[vertex] = uint32(len(m.V))
				m.V = append(m.V, vertex)
				m.T = append(m.T, vMap[vertex])
			} else {
				m.T = append(m.T, idx)
			}
		}
	}

	return m
}

/* OBJ */

// This function is used when we have some vertices and no faces
// func WriteToOBJ(v []geom.Vertex, filepath string) error {
// 	file, err := os.Create(filepath)
// 	if err != nil {
// 		return err
// 	}
// 	defer file.Close()

// 	// Write the vertices
// 	for _, vert := range v {
// 		_, err := fmt.Fprintf(file, "v %f %f %f\n", vert.X, vert.Y, vert.Z)
// 		if err != nil {
// 			return err
// 		}
// 	}

// 	// Write the faces
// 	for i, _ := range v {
// 		_, err = fmt.Fprintf(file, "f %d %d\n", i+1, (i+1)%len(v)+1)
// 		if err != nil {
// 			return err
// 		}
// 	}

// 	return nil
// }

// This function is used to save a mesh to filepath
func WriteToOBJ(m geom.Mesh, filepath string) error {

	file, err := os.Create(filepath)
	if err != nil {
		return err
	}
	defer file.Close()

	// Write the vertices
	for _, v := range m.V {
		_, err := fmt.Fprintf(file, "v %f %f %f\n", v.X, v.Y, v.Z)
		if err != nil {
			return err
		}
	}

	// Write the faces
	for i := 0; i < len(m.T); i += 3 {
		// One has to be added to every triangle index
		_, err := fmt.Fprintf(file, "f %d %d %d\n", m.T[i]+1, m.T[i+1]+1, m.T[i+2]+1)
		if err != nil {
			return err
		}
	}

	return nil
}

// Utility function that reads OBJ into slice of vertices and faces
func ReadObj(filePath string) ([]geom.Vertex, []geom.Face, error) {
	var vertices []geom.Vertex
	var faces []geom.Face

	file, err := os.Open(filePath)
	if err != nil {
		return nil, nil, err
	}

	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		words := strings.Fields(line)
		// if line is empty continue
		if len(words) == 0 {
			continue
		}

		switch words[0] {
		case "v":
			var vertex geom.Vertex
			// if length of vertex line is not 4 there's a problem
			if len(words) != 4 {
				return nil, nil, fmt.Errorf("Vertices must be 3D")
			}

			value, _ := strconv.ParseFloat(words[1], 32)
			vertex.X = float32(value)
			value, _ = strconv.ParseFloat(words[2], 32)
			vertex.Y = float32(value)
			value, _ = strconv.ParseFloat(words[3], 32)
			vertex.Z = float32(value)

			vertices = append(vertices, vertex)

		case "f":
			var face geom.Face

			for _, word := range words[1:] {
				value, err := strconv.Atoi(word)
				if err != nil {
					return nil, nil, err
				}
				// Adjusting the index to be 0-based
				face = append(face, value-1)
			}
			faces = append(faces, face)
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, nil, err
	}

	return vertices, faces, nil
}
