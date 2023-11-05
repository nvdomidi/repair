package geom

import "math"

/* Data structures */

type Vertex struct {
	X, Y, Z float32
}

type Mesh struct {
	V []Vertex
	T []uint32
}

type Line struct {
	Dir Vertex // line direction
	P   Vertex // some point on the line
}

type Plane struct {
	Normal Vertex
	P      Vertex // some point on the plane
}

type Face []int

/* Vertex Operations */

// A + B
func Add(a, b Vertex) Vertex {
	return Vertex{a.X + b.X, a.Y + b.Y, a.Z + b.Z}
}

// A - B
func Subtract(a, b Vertex) Vertex {
	return Vertex{a.X - b.X, a.Y - b.Y, a.Z - b.Z}
}

// A x B
func CrossProduct(a, b Vertex) Vertex {
	return Vertex{
		X: a.Y*b.Z - a.Z*b.Y,
		Y: a.Z*b.X - a.X*b.Z,
		Z: a.X*b.Y - a.Y*b.X,
	}
}

// A . B
func DotProduct(a, b Vertex) float32 {
	return a.X*b.X + a.Y*b.Y + a.Z*b.Z
}

// pA
func Multiply(p Vertex, scalar float32) Vertex {
	return Vertex{p.X * scalar, p.Y * scalar, p.Z * scalar}
}

// |A|
func Length(a Vertex) float32 {
	return float32(math.Sqrt(float64(a.X)*float64(a.X) + float64(a.Y)*float64(a.Y) + float64(a.Z)*float64(a.Z)))
}

// A / |A|
func Normalize(a Vertex) Vertex {
	length := Length(a)
	return Vertex{a.X / length, a.Y / length, a.Z / length}
}

func AlmostEqual(a, b float32) bool {
	return math.Abs(float64(a-b)) < 1e-7
}

func PointsEqual(a, b Vertex) bool {
	return AlmostEqual(a.X, b.X) && AlmostEqual(a.Y, b.Y) && AlmostEqual(a.Z, b.Z)
}
