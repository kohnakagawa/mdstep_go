package main

import "fmt"

type MD struct {
	observer Observer
	vars     Variables
}

func MakeMdSystem() MD {
	return MD{Observer{}, Variables{MakeParticles(), 0.0}}
}

func (md *MD) calculate() {
	kick(md.vars.particles)
	CalcForceBruteForce(md.vars.particles)
	kick(md.vars.particles)
	adjustPeriodic(md.vars.particles)
	md.vars.time += Dt
}

func (md *MD) Run() {
	for i := 0; i < NSteps; i++ {
		if i%NStepsObs == 0 {
			k := md.observer.CalcKineticEnergy(md.vars.particles)
			v := md.observer.CalcPotentialEnergy(md.vars.particles)
			fmt.Println(md.vars.time, k, v, k+v)
		}
		md.calculate()
	}
}

func FixCMDrift(Particles []Particle) {
	cmVel := Double3{0.0, 0.0, 0.0}
	for i := range Particles {
		cmVel.x += Particles[i].vel.x
		cmVel.y += Particles[i].vel.y
		cmVel.z += Particles[i].vel.z
	}

	cmVel.x /= float64(len(Particles))
	cmVel.y /= float64(len(Particles))
	cmVel.z /= float64(len(Particles))

	for i := range Particles {
		Particles[i].vel.x -= cmVel.x
		Particles[i].vel.y -= cmVel.y
		Particles[i].vel.z -= cmVel.z
	}
}

func adjustPeriodic(Particles []Particle) {
	for i := range Particles {
		if Particles[i].pos.x < 0 {
			Particles[i].pos.x += BoxLen
		}
		if Particles[i].pos.x > BoxLen {
			Particles[i].pos.x -= BoxLen
		}
		if Particles[i].pos.y < 0 {
			Particles[i].pos.y += BoxLen
		}
		if Particles[i].pos.y > BoxLen {
			Particles[i].pos.y -= BoxLen
		}
		if Particles[i].pos.z < 0 {
			Particles[i].pos.z += BoxLen
		}
		if Particles[i].pos.z > BoxLen {
			Particles[i].pos.z -= BoxLen
		}
	}
}

func kick(Particles []Particle) {
	for i := range Particles {
		Particles[i].pos.x += 0.5 * Particles[i].vel.x * Dt
		Particles[i].pos.y += 0.5 * Particles[i].vel.y * Dt
		Particles[i].pos.z += 0.5 * Particles[i].vel.z * Dt
	}
}
