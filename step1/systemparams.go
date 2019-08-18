package main

const (
	BoxLen    = 10.0
	BoxHalf   = BoxLen * 0.5
	Dt        = 0.01
	Cutoff    = 2.0
	Cutoff2   = Cutoff * Cutoff
	Density   = 0.5
	NSteps    = 10000
	NStepsObs = 100
	RC2       = 1.0 / Cutoff2
	C0        = -4.0 * RC2 * RC2 * RC2 * (RC2*RC2*RC2 - 1.0)
)

func GetMinimumImage(v Double3) Double3 {
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
