package main

import (
	"math"
	"math/rand"
)

type Double3 struct {
	x float64
	y float64
	z float64
}

func Double3Norm(d Double3) float64 {
	return d.x*d.x + d.y*d.y + d.z*d.z
}

type Pair struct {
	i int
	j int
}

type Particle struct {
	pos Double3
	vel Double3
}

type Variables struct {
	particles []Particle
	time      float64
	pairs     []Pair
	margin    float64
}

func makeParticle(x float64, y float64, z float64, v0 float64) Particle {
	pos := Double3{x, y, z}
	vz := 2.0*rand.Float64() - 1.0
	phi := 2.0 * rand.Float64() * math.Pi
	vel := Double3{v0 * math.Sqrt(1.0-vz*vz) * math.Cos(phi), v0 * math.Sqrt(1.0-vz*vz) * math.Sin(phi), v0 * vz}
	return Particle{pos, vel}
}

func MakeParticles() []Particle {
	rand.Seed(0)
	s := 1.0 / math.Pow(Density*0.25, 1.0/3.0)
	hs := 0.5 * s
	is := int(BoxLen / s)
	Particles := make([]Particle, 0)
	for iz := 0; iz < is; iz++ {
		for iy := 0; iy < is; iy++ {
			for ix := 0; ix < is; ix++ {
				ixD := float64(ix)
				iyD := float64(iy)
				izD := float64(iz)
				Particles = append(Particles, makeParticle(ixD*s, iyD*s, izD*s, 1.0))
				Particles = append(Particles, makeParticle(ixD*s+hs, iyD*s, izD*s, 1.0))
				Particles = append(Particles, makeParticle(ixD*s, iyD*s+hs, izD*s, 1.0))
				Particles = append(Particles, makeParticle(ixD*s, iyD*s, izD*s+hs, 1.0))
			}
		}
	}
	FixCMDrift(Particles)

	return Particles
}
