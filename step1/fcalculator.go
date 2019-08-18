package main

func CalcForceBruteForce(particles []Particle) {
	nParticles := len(particles)
	for i := 0; i < nParticles-1; i++ {
		for j := i + 1; j < nParticles; j++ {
			qi := particles[i].pos
			qj := particles[j].pos
			dqij := Double3{qj.x - qi.x, qj.y - qi.y, qj.z - qi.z}
			dqij = GetMinimumImage(dqij)
			dq2 := dqij.x*dqij.x + dqij.y*dqij.y + dqij.z*dqij.z
			if dq2 > Cutoff2 {
				continue
			}
			dq6 := dq2 * dq2 * dq2
			df := (24.0*dq6 - 48.0) / (dq6 * dq6 * dq2) * Dt
			particles[i].vel.x += df * dqij.x
			particles[i].vel.y += df * dqij.y
			particles[i].vel.z += df * dqij.z
			particles[j].vel.x -= df * dqij.x
			particles[j].vel.y -= df * dqij.y
			particles[j].vel.z -= df * dqij.z
		}
	}
}
