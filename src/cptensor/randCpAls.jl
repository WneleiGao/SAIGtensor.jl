"""
    sidx2sub!(I, s)

convert linear indices to sub-index

# Arguments
* `I`: Array{Int64,1}, dimensions
* `s`: Array{Int64,1}, linear index

# Returns
* `sidx`: Array{Tuple,1} with length(I)

"""
function sidx2sub(I::Array{Int64,1}, s::Array{Int64,1})
    N = length(I)
    sidx = Array{Tuple}(N)
    In = zeros(Int64, N-1)
    for n = 1 : N
        tmp = setdiff(collect(1:N), [n])
        for i = 1 : N-1
            In[i] = I[tmp[i]]
        end
        sidx[n] = ind2sub((In...), s)
    end
    return sidx
end

"""
    randSpl!(Zs, Xs, A, T, ns, n)

random sampling certain rows from Z[n]

# Arguments
* `Zs`: Array{Float64,2}, sampled Z[n]
* `Xst`: Array{Array{Float64,2},1}, transpose of sampled X[n]
* `A`: Array{Array{Float64,2},1}, factor matrix
* `ns`: Int64, number of samples
* `n`: Int64, specify which dimension
"""
function randSpl!(Zs::Array{Float64,2}, Xst::Array{Array{Float64,2},1}, A::Array{Array{Float64,2},1}, T::tensor, ns::Int64, n::Int64)
    N = T.N
    length(A) == N || throw(DimensionMismatch())
    (l1, l2) = size(Zs)
    (l1, l2) == (ns, size(A[1],2)) || throw(DimensionMismatch())
    dim = setdiff(collect(1:N), [n])
    # generate random multi-linear index
    idx = Array{Array{Int64,1}}(N-1)
    for i = 1 : N-1
        In = size(A[dim[i]],1)
        In == T.I[dim[i]] || throw(DimensionMismatch())
        idx[i] = rand(collect(1:In), ns)
    end
    formZs!(Zs, A, idx, dim);
    formXst!(Xst[n], T.D, idx, dim, n)
    return idx
end

"""
   random sampling X[n]^T
"""
function formXst!(Xst::Array{Float64,2}, D::Array{Float64,4}, idx::Array{Array{Int64,1},1}, dim::Array{Int64,1}, n::Int64)
    (l1,l2) = size(Xst)
    N = ndims(D)
    ind = zeros(Int64,N)
    for i1 = 1 : l1
        for i = 1 : N-1
            ind[dim[i]] = idx[i][i1]
        end
        for i2 = 1 : l2
            ind[n] = i2
            Xst[i1, i2] = D[ind[1],ind[2],ind[3],ind[4]]
        end
    end
    return nothing
end

"""
  random sampling Z[n]
"""
function formZs!(Zs::Array{Float64,2}, A::Array{Array{Float64,2},1}, idx::Array{Array{Int64,1},1}, dim::Array{Int64,1})
    (l1,l2) = size(Zs)
    N = length(A)
    for i2 = 1 : l2
        for i1 = 1 : l1
            Zs[i1,i2] = 1.
        end
    end
    for i = 1 : N-1
        for i2 = 1: l2
            for i1 = 1 : l1
                Zs[i1,i2] = Zs[i1,i2] * A[dim[i]][idx[i][i1],i2]
            end
        end
    end
    return nothing
end

function updateAn!(lambda::Array{Float64,1}, A::Array{Array{Float64,2},1}, Zs::Array{Float64,2}, Xst::Array{Float64,2}, n::Int64)
    R = size(A[1],2)
    In = size(A[n],1)
    tmp = (Zs \ Xst)'
    for i2 = 1 : R
        lambda[i2] = maximum(abs(tmp))
        for i1 = 1 : In
            tmp[i1,i2] = tmp[i1,i2] / lambda[i2]
        end
    end
    A[n] = tmp
    return nothing
end

function randCpAls(X::tensor, R::Int64; maxiter::Int64=50, fitchangetol=1e-5, pflag=false)
    N = X.N; I = X.I;
    # minimum(X.I) > R || error("too large rank")
    normX = tnorm(X)
    fit = 0.
    # initialize factor matrix
    lambda = zeros(R)
    A = Array{Array{Float64,2}}(N)
    for m = 1 : N
        A[m] = randn(I[m], R)
    end
    ns = ceil(Int64, 10*R*log(R))
    Zs = zeros(ns, R)
    Xst = Array{Array{Float64,2}}(N)
    for i = 1 : N
        Xst[i] = zeros(ns, I[i])
    end
    println("random CP_ALS")  # alternating Least square
    for iter = 1 : maxiter
        fitold = fit
        # updating ns factor
        for n = 1 : N
            randSpl!(Zs, Xst, A, X, ns, n)
            updateAn!(lambda, A, Zs, Xst[n], n)
        end
        # check the fitness
        normresidual = sqrt(normX^2 + cpnorm(lambda,A)^2 - 2*innerprod(lambda, A, X))
        fit = 1 - normresidual/normX
        fitchange = abs(fitold-fit)
        if (iter > 1) && fitchange < fitchangetol
           flag = 0
        else
           flag = 1
        end
        if pflag
           println("iter: $iter, fit: $fit, fitchange: $fitchange")
        end
        if flag == 0
           break
        end
    end
    cpnormalize!(lambda, A)
    return lambda, A
end
