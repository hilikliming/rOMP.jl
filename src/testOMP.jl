include("rOMP.jl")

using Main.rOMP
using LinearAlgebra

N=200

D=zeros(N,N)
# Generating DCT-I basis.
for kk in 0:N-1
    for nn in 0:N-1
        if kk == 0
            ll=1/sqrt(N)
        else
            ll=sqrt(2/N)*cos((2*nn+1)*(kk)*pi/(2*N))
        end
        D[kk+1,nn+1]=ll
    end
end

y= sum(D[:,1:2:5],dims=2) # so y is the sum of the 3rd and fourth frequency components

@time x= rOMP.omp(y,D,3,1e-6)

err=y-D*x
display(err'*err)
