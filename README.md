# rOMP.jl
A recursive implementation of Orthogonal Matching Pursuit implemented in Julia.
This method uses the block inversion formulas to incrementally update the least squares solution as atoms are selected. 


I assume that D_I'D_I = G the gram matrix of the sub dictionary corresponding to indices I takes the form G=[A,B;C,D] and we assume the 1x1 inner product D, and G's schur complement are both invertible which should always hold since the consecutive atoms should be linearly independent (if a potential "new" atom were in the span of the old ones, the argmax step of OMP would see 0 inner product between the residual and this potential "new" atom since the residual is in the perp space of the subspace that this atom lies in... so it wouldn't get picked to be the new atom)
