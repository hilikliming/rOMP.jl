module rOMP
export omp
using Base.Threads, Random, SparseArrays, LinearAlgebra

Random.seed!(1234)  # for stability of tests

function omp(data::Array{Float64,2}, dictionary::Array{Float64,2}, tau::Int, tolerance::Float64)
    K= size(dictionary,2)
    X_sparse=zeros(K,size(data,2))
    for nn=1:size(data,2)
        x = deepcopy(data[:,nn]) # a single data sample to be represented in terms of atoms in D
        r = deepcopy(x) # residual vector
        D = dictionary # Dictionary
        # Note that the o (optimal) variables are used to avoid an uncommon
        # scenario (that does occur) where a lower sparsity solution may have
        # had lower error than the final solution (with tau non zeros) but
        # wasn't low enough to break out of the coefficient solver via the error
        # tolerance. A litte more memory for significantly better solutions,
        # thanks to CR for the tip (JJH)
        γ       = 0 # this will be the growing coefficient vector
        γₒ      = 0 # this will store whatever the minimum error solution was during computation of the coefficients
        av_err  = 1e6 # norm of the error vector.
        best_err= 1e6 # will store lowest error vector norm
        ii      = 1   # while loop index
        DI      = []  # This holds the atoms selected via OMP as its columns (it grows along 2nd dimension)
        DIGI    = []  # Inverse of DI's gram matrix
        DIdag   = []  # PseudoInverse of DI
        I       = []  # set of indices corresponding to atoms selected in reconstruction
        Iₒ      = []  # I think you get the deal with these guys now (best set of indices lul)
        while (length(I)<tau) && (av_err > tolerance)
            k = argmax(broadcast(abs,D'*r))
            dk= D[:,k]
            if ii==1
                I = k
                #display("we made it")
                DI=dk
                DIGI=(DI'*DI)^(-1)
            else
                I = cat(dims=1,I,k)
                rho=DI'*dk
                DI=cat(dims=2,DI,dk)
                ipk=dk'*dk
                DIGI=blockMatrixInv(DIGI,rho,rho,ipk)
            end
            DIdag   = DIGI*DI'
            γ   = DIdag*x
            r       = x-DI*γ
            av_err  = norm(r)
            if av_err<= best_err
                best_err=av_err
                γₒ= γ
                Iₒ=I
            end
            X_sparse[I,nn]=γ
            ii+=1
        end
        if av_err > best_err
            X_sparse[I,nn]= 0*X_sparse[I,nn]
            X_sparse[Iₒ,nn]=γₒ
        end

    end
    return X_sparse
end

function blockMatrixInv(Ai::Array{Float64,2}, B::Array{Float64,1}, C::Array{Float64,1}, D::Float64)
    C=C'
    DCABi= (D-C*Ai*B)^(-1)
    return [Ai+Ai*B*DCABi*C*Ai -Ai*B*DCABi; -DCABi*C*Ai DCABi]
end

function blockMatrixInv(Ai::Float64, B::Float64, C::Float64, D::Float64)
    DCABi= (D-C*Ai*B)^(-1)
    return [Ai+Ai*B*DCABi*C*Ai -Ai*B*DCABi; -DCABi*C*Ai DCABi]
end

end # module
