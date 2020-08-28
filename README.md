# rOMP.jl
A recursive implementation of Orthogonal Matching Pursuit implemented in Julia.
This method uses the block inversion formulas to incrementally update the least squares solution as atoms are selected. Assumes D and the augmented dictionaries schur complement are both invertible which should always hold since the consecutive atoms should be linearly independent (if the new atom were in the span of the old ones the first step of OMP would see 0 inner product between the residual and the new atom since the residual is in the perp space of the the subspace that the new atom lies in.)
