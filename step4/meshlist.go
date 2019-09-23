package main

import (
	"fmt"
	"math"
)

type MeshList struct {
	meshSize     float64
	meshDim      int
	numberOfMesh int
	count        []int
	indexes      []int
	sortedBuffer []int

	particlePosition []int
	meshPointer      []int

	neighborList    []int
	iPosition       []int
	jCount          []int
	particlePointer []int
}

func NewMeshList(vars Variables) MeshList {
	m := int(BoxLen/SL) - 1 // Assume: m > 1
	meshSize := BoxLen / float64(m)
	numberOfMesh := m * m * m
	numberOfParticles := len(vars.particles)
	initialBufferSize := int(float64(numberOfParticles) * Density * (4.0 / 3.0) * math.Pi * SL2 * SL)
	return MeshList{
		meshSize:         meshSize,
		meshDim:          m,
		numberOfMesh:     numberOfMesh,
		count:            make([]int, numberOfMesh),
		indexes:          make([]int, numberOfMesh),
		sortedBuffer:     make([]int, numberOfParticles),
		particlePosition: make([]int, numberOfParticles),
		meshPointer:      make([]int, numberOfMesh),
		neighborList:     make([]int, 0, initialBufferSize),
		iPosition:        make([]int, numberOfParticles+1),
		jCount:           make([]int, numberOfParticles),
		particlePointer:  make([]int, numberOfParticles),
	}
}

func getHash(pos Double3, im float64, m int) int {
	ix := int(pos.x * im)
	iy := int(pos.y * im)
	iz := int(pos.z * im)
	if ix < 0 {
		ix += m
	}
	if ix >= m {
		ix -= m
	}
	if iy < 0 {
		iy += m
	}
	if iy >= m {
		iy -= m
	}
	if iz < 0 {
		iz += m
	}
	if iz >= m {
		iz -= m
	}

	return ix + m*(iy+m*iz)
}

func (meshlist *MeshList) searchSameMesh(id int, ix int, iy int, iz int, particles []Particle, pairs []Pair) []Pair {
	ms2 := meshlist.meshSize * meshlist.meshSize

	for k := meshlist.indexes[id]; k < meshlist.indexes[id]+meshlist.count[id]-1; k++ {
		for l := k + 1; l < meshlist.indexes[id]+meshlist.count[id]; l++ {
			i := meshlist.sortedBuffer[k]
			j := meshlist.sortedBuffer[l]
			dx := particles[i].pos.x - particles[j].pos.x
			dy := particles[i].pos.y - particles[j].pos.y
			dz := particles[i].pos.z - particles[j].pos.z
			dq := Double3{dx, dy, dz}
			dq = GetMinimumImage(dq)
			dqNorm2 := Double3Norm(dq)
			if dqNorm2 > ms2 {
				continue
			}
			pairs = append(pairs, Pair{i, j})
		}
	}
	return pairs
}

func (meshlist *MeshList) searchDifferentMesh(id int, jx int, jy int, jz int, particles []Particle, pairs []Pair) []Pair {
	m := meshlist.meshDim
	ms2 := meshlist.meshSize * meshlist.meshSize

	if jx < 0 {
		jx += m
	}
	if jx >= m {
		jx -= m
	}
	if jy < 0 {
		jy += m
	}
	if jy >= m {
		jy -= m
	}
	if jz < 0 {
		jz += m
	}
	if jz >= m {
		jz -= m
	}

	id2 := jx + m*(jy+jz*m)

	for k := meshlist.indexes[id]; k < meshlist.indexes[id]+meshlist.count[id]; k++ {
		for l := meshlist.indexes[id2]; l < meshlist.indexes[id2]+meshlist.count[id2]; l++ {
			i := meshlist.sortedBuffer[k]
			j := meshlist.sortedBuffer[l]
			dx := particles[i].pos.x - particles[j].pos.x
			dy := particles[i].pos.y - particles[j].pos.y
			dz := particles[i].pos.z - particles[j].pos.z
			dq := Double3{dx, dy, dz}
			dq = GetMinimumImage(dq)
			dqNorm2 := Double3Norm(dq)
			if dqNorm2 > ms2 {
				continue
			}
			pairs = append(pairs, Pair{i, j})
		}
	}

	return pairs
}

func (meshlist *MeshList) search(ix int, iy int, iz int, vars *Variables) {
	m := meshlist.meshDim
	id := ix + m*(iy+m*iz)

	vars.pairs = meshlist.searchDifferentMesh(id, ix-1, iy-1, iz-1, vars.particles, vars.pairs)
	vars.pairs = meshlist.searchDifferentMesh(id, ix, iy-1, iz-1, vars.particles, vars.pairs)
	vars.pairs = meshlist.searchDifferentMesh(id, ix+1, iy-1, iz-1, vars.particles, vars.pairs)
	vars.pairs = meshlist.searchDifferentMesh(id, ix-1, iy, iz-1, vars.particles, vars.pairs)
	vars.pairs = meshlist.searchDifferentMesh(id, ix, iy, iz-1, vars.particles, vars.pairs)
	vars.pairs = meshlist.searchDifferentMesh(id, ix+1, iy, iz-1, vars.particles, vars.pairs)
	vars.pairs = meshlist.searchDifferentMesh(id, ix-1, iy+1, iz-1, vars.particles, vars.pairs)
	vars.pairs = meshlist.searchDifferentMesh(id, ix, iy+1, iz-1, vars.particles, vars.pairs)
	vars.pairs = meshlist.searchDifferentMesh(id, ix+1, iy+1, iz-1, vars.particles, vars.pairs)
	vars.pairs = meshlist.searchDifferentMesh(id, ix-1, iy-1, iz, vars.particles, vars.pairs)
	vars.pairs = meshlist.searchDifferentMesh(id, ix, iy-1, iz, vars.particles, vars.pairs)
	vars.pairs = meshlist.searchDifferentMesh(id, ix+1, iy-1, iz, vars.particles, vars.pairs)
	vars.pairs = meshlist.searchDifferentMesh(id, ix-1, iy, iz, vars.particles, vars.pairs)

	vars.pairs = meshlist.searchSameMesh(id, ix, iy, iz, vars.particles, vars.pairs)
}

func (meshlist *MeshList) MakePairs(vars *Variables) {
	vars.pairs = vars.pairs[:0]

	for i := range meshlist.particlePosition {
		meshlist.particlePosition[i] = 0
	}
	for i := range meshlist.count {
		meshlist.count[i] = 0
	}
	for i := range meshlist.meshPointer {
		meshlist.meshPointer[i] = 0
	}

	nParticles := len(vars.particles)
	im := 1.0 / meshlist.meshSize
	for i, p := range vars.particles {
		index := getHash(p.pos, im, meshlist.meshDim)
		if index < 0 || index > meshlist.meshDim*meshlist.meshDim*meshlist.meshDim {
			fmt.Println(index)
			fmt.Println(p.pos)
		}
		meshlist.count[index]++
		meshlist.particlePosition[i] = index
	}

	meshlist.indexes[0] = 0
	partialSum := 0
	for i := 0; i < meshlist.numberOfMesh-1; i++ {
		partialSum += meshlist.count[i]
		meshlist.indexes[i+1] = partialSum
	}

	for i := 0; i < nParticles; i++ {
		pos := meshlist.particlePosition[i]
		j := meshlist.indexes[pos] + meshlist.meshPointer[pos]
		meshlist.sortedBuffer[j] = i
		meshlist.meshPointer[pos]++
	}

	for ix := 0; ix < meshlist.meshDim; ix++ {
		for iy := 0; iy < meshlist.meshDim; iy++ {
			for iz := 0; iz < meshlist.meshDim; iz++ {
				meshlist.search(ix, iy, iz, vars)
			}
		}
	}
}

func (meshlist *MeshList) MakeNeighborList(vars *Variables) {
	numberOfPairs := len(vars.pairs)
	if cap(meshlist.neighborList) >= numberOfPairs {
		meshlist.neighborList = meshlist.neighborList[:numberOfPairs]
	} else {
		meshlist.neighborList = make([]int, numberOfPairs)
	}

	for i := range meshlist.jCount {
		meshlist.jCount[i] = 0
	}
	for _, p := range vars.pairs {
		meshlist.jCount[p.i]++
	}

	numberOfParticles := len(vars.particles)
	meshlist.iPosition[0] = 0
	for i := 0; i < numberOfParticles; i++ {
		meshlist.iPosition[i+1] = meshlist.iPosition[i] + meshlist.jCount[i]
	}

	for i := range meshlist.particlePointer {
		meshlist.particlePointer[i] = 0
	}
	for _, p := range vars.pairs {
		pos := meshlist.iPosition[p.i] + meshlist.particlePointer[p.i]
		meshlist.neighborList[pos] = p.j
		meshlist.particlePointer[p.i]++
	}
}
