package main

import (
	"fmt"
	"math"
)

type MD struct {
	observer Observer
	vars     Variables
	meshlist MeshList
}

func NewMD() MD {
	vars := NewVariable()
	return MD{NewObserver(), vars, NewMeshList(vars)}
}

func (md *MD) makePairsBruteforce() {
	md.vars.pairs = md.vars.pairs[:0]
	nParticles := len(md.vars.particles)
	for i := 0; i < nParticles-1; i++ {
		for j := i + 1; j < nParticles; j++ {
			qi := md.vars.particles[i].pos
			qj := md.vars.particles[j].pos
			dqij := Double3{qj.x - qi.x, qj.y - qi.y, qj.z - qi.z}
			dqij = GetMinimumImage(dqij)
			dq2 := dqij.x*dqij.x + dqij.y*dqij.y + dqij.z*dqij.z
			if dq2 > SL2 {
				continue
			}
			md.vars.pairs = append(md.vars.pairs, Pair{i, j})
		}
	}
}

func (md *MD) checkPairs() {
	vmax2 := 0.0
	for _, particle := range md.vars.particles {
		v2 := Double3Norm(particle.vel)
		if vmax2 < Double3Norm(particle.vel) {
			vmax2 = v2
		}
	}
	vmax := math.Sqrt(vmax2)

	md.vars.margin -= vmax * 2.0 * Dt
	if md.vars.margin < 0.0 {
		md.vars.margin = Margin
		md.meshlist.makePairs(&md.vars)
	}
}

func (md *MD) calculate() {
	kick(md.vars.particles)
	md.checkPairs()
	CalcForcePair(md.vars.particles, md.vars.pairs)
	kick(md.vars.particles)
	adjustPeriodic(md.vars.particles)
	md.vars.time += Dt
}

func (md *MD) Run() {
	md.meshlist.makePairs(&md.vars)
	for i := 0; i < NSteps; i++ {
		if i%NStepsObs == 0 {
			k := md.observer.CalcKineticEnergy(md.vars.particles)
			v := md.observer.CalcPotentialEnergy(md.vars.particles, md.vars.pairs)
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

func showPairs(pairs []Pair) {
	for _, pair := range pairs {
		fmt.Println(pair.i, pair.j)
	}
}
