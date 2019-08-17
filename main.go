package main

import (
	"fmt"
	"math"
	"math/rand"
)

const (
	BoxLen    = 10.0
	BoxHalf   = BoxLen * 0.5
	dt        = 0.01
	cutoff    = 2.0
	cutoff2   = cutoff * cutoff
	density   = 0.5
	nSteps    = 10000
	nStepsObs = 100
)

type double3 struct {
	x float64
	y float64
	z float64
}

type particle struct {
	pos double3
	vel double3
}

func makeParticle(x float64, y float64, z float64, v0 float64) particle {
	pos := double3{x, y, z}
	vz := 2.0*rand.Float64() - 1.0
	phi := 2.0 * rand.Float64() * math.Pi
	vel := double3{v0 * math.Sqrt(1.0-vz*vz) * math.Cos(phi), v0 * math.Sqrt(1.0-vz*vz) * math.Sin(phi), v0 * vz}
	return particle{pos, vel}
}

func fixCMDrift(particles []particle) {
	cmVel := double3{0.0, 0.0, 0.0}
	for i := range particles {
		cmVel.x += particles[i].vel.x
		cmVel.y += particles[i].vel.y
		cmVel.z += particles[i].vel.z
	}

	cmVel.x /= float64(len(particles))
	cmVel.y /= float64(len(particles))
	cmVel.z /= float64(len(particles))

	for i := range particles {
		particles[i].vel.x -= cmVel.x
		particles[i].vel.y -= cmVel.y
		particles[i].vel.z -= cmVel.z
	}
}

func makeParticles() []particle {
	rand.Seed(0)
	s := 1.0 / math.Pow(density*0.25, 1.0/3.0)
	hs := 0.5 * s
	is := int(BoxLen / s)
	particles := make([]particle, 0)
	for iz := 0; iz < is; iz++ {
		for iy := 0; iy < is; iy++ {
			for ix := 0; ix < is; ix++ {
				ixD := float64(ix)
				iyD := float64(iy)
				izD := float64(iz)
				particles = append(particles, makeParticle(ixD*s, iyD*s, izD*s, 1.0))
				particles = append(particles, makeParticle(ixD*s+hs, iyD*s, izD*s, 1.0))
				particles = append(particles, makeParticle(ixD*s, iyD*s+hs, izD*s, 1.0))
				particles = append(particles, makeParticle(ixD*s, iyD*s, izD*s+hs, 1.0))
			}
		}
	}
	fixCMDrift(particles)

	return particles
}

func getMinimumImage(v double3) double3 {
	if v.x < -BoxHalf {
		v.x += BoxLen
	}
	if v.x > BoxHalf {
		v.x -= BoxLen
	}
	if v.y < -BoxHalf {
		v.y += BoxLen
	}
	if v.y > BoxHalf {
		v.y -= BoxLen
	}
	if v.z < -BoxHalf {
		v.z += BoxLen
	}
	if v.z > BoxHalf {
		v.z -= BoxLen
	}
	return v
}

func adjustPeriodic(particles []particle) {
	for i := range particles {
		if particles[i].pos.x < 0 {
			particles[i].pos.x += BoxLen
		}
		if particles[i].pos.x > BoxLen {
			particles[i].pos.x -= BoxLen
		}
		if particles[i].pos.y < 0 {
			particles[i].pos.y += BoxLen
		}
		if particles[i].pos.y > BoxLen {
			particles[i].pos.y -= BoxLen
		}
		if particles[i].pos.z < 0 {
			particles[i].pos.z += BoxLen
		}
		if particles[i].pos.z > BoxLen {
			particles[i].pos.z -= BoxLen
		}
	}
}

func calcForce(particles []particle) {
	nParticles := len(particles)
	for i := 0; i < nParticles-1; i++ {
		for j := i + 1; j < nParticles; j++ {
			qi := particles[i].pos
			qj := particles[j].pos
			dqij := double3{qj.x - qi.x, qj.y - qi.y, qj.z - qi.z}
			dqij = getMinimumImage(dqij)
			dq2 := dqij.x*dqij.x + dqij.y*dqij.y + dqij.z*dqij.z
			if dq2 > cutoff2 {
				continue
			}
			dq6 := dq2 * dq2 * dq2
			df := (24.0*dq6 - 48.0) / (dq6 * dq6 * dq2) * dt
			particles[i].vel.x += df * dqij.x
			particles[i].vel.y += df * dqij.y
			particles[i].vel.z += df * dqij.z
			particles[j].vel.x -= df * dqij.x
			particles[j].vel.y -= df * dqij.y
			particles[j].vel.z -= df * dqij.z
		}
	}
}

func kick(particles []particle) {
	for i := range particles {
		particles[i].pos.x += 0.5 * particles[i].vel.x * dt
		particles[i].pos.y += 0.5 * particles[i].vel.y * dt
		particles[i].pos.z += 0.5 * particles[i].vel.z * dt
	}
}

func calcKineticEnergy(particles []particle) float64 {
	k := 0.0
	for i := range particles {
		k += particles[i].vel.x * particles[i].vel.x
		k += particles[i].vel.y * particles[i].vel.y
		k += particles[i].vel.z * particles[i].vel.z
	}
	k /= float64(len(particles))
	return 0.5 * k
}

func main() {
	particles := makeParticles()

	for time := 0; time < nSteps; time++ {
		kick(particles)
		calcForce(particles)
		kick(particles)
		adjustPeriodic(particles)

		if time%nStepsObs == 0 {
			fmt.Println(time, calcKineticEnergy(particles))
		}
	}
}
