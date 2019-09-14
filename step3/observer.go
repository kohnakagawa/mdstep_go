package main

type Observer struct {
}

func NewObserver() Observer {
	return Observer{}
}

func (observer *Observer) CalcKineticEnergy(particles []Particle) float64 {
	k := 0.0
	for i := range particles {
		k += particles[i].vel.x * particles[i].vel.x
		k += particles[i].vel.y * particles[i].vel.y
		k += particles[i].vel.z * particles[i].vel.z
	}
	k /= float64(len(particles))
	return 0.5 * k
}

func (observer *Observer) CalcPotentialEnergy(particles []Particle, pairs []Pair) float64 {
	v := 0.0
	nParticles := len(particles)
	for _, pair := range pairs {
		qi := particles[pair.i].pos
		qj := particles[pair.j].pos
		dqij := Double3{qj.x - qi.x, qj.y - qi.y, qj.z - qi.z}
		dqij = GetMinimumImage(dqij)
		dq2 := dqij.x*dqij.x + dqij.y*dqij.y + dqij.z*dqij.z
		if dq2 > Cutoff2 {
			continue
		}
		dq6 := dq2 * dq2 * dq2
		dq12 := dq6 * dq6
		v += 4.0*(1.0/dq12-1.0/dq6) + C0
	}
	v /= float64(nParticles)
	return v
}
