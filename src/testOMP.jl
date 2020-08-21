include("rOMP.jl")

using Main.rOMP
using LinearAlgebra

N=200

D=zeros(N,N)
# Generating DCT-I basis... There are definitely faster ways to do this...
for kk in 0:N-1
    for nn in 0:N-1
        if kk == 0
            ll=1/sqrt(N)
        else
            ll=sqrt(2/N)*cos((2*nn+1)*(kk)*pi/(2*N)) # don't forget the root N terms so it's orthonormal not just orthogonal
        end
        D[kk+1,nn+1]=ll
    end
end

y= sum(D[:,1:2:5],dims=2) # so y is the sum of the 1st, 3rd, and 5th frequency components

@time x= rOMP.omp(y,D,3,1e-6) # Note we give ourselves the best shot by setting Tau = 3 the number of components we already know are in y

err=y-D*x
display(err'*err)
