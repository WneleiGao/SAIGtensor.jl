"""
composite type for tucker tensor
"""
type tucker{Tv<:AbstractFloat}
     N  :: Int64
     rk :: Vector{Int64}
     I  :: Vector{Int64}
     core :: Array{Tv}
     U  :: Vector{Matrix{Tv}}
end

"""
  initialize a tucker tensor
"""
function initTucker(T::tucker)
    T1 = Tucker(T.N, copy(T.rk), copy(T.I), copy(T.core), deepcopy(U))
    return T1
end

"""
etp = Float32
rk = [13, 14, 37, 31]
I  = [33, 34, 35, 38]
core = rand(etp, rk...)
U  = Vector{Matrix{Float32}}(ndims(core))
for i = 1 : length(rk)
    U[i] = rand(etp, I[i], rk[i])
end
T = InitTucker(core, U)
"""
function initTucker{Tv<:AbstractFloat}(core::Array{Tv}, U::Vector{Matrix{Tv}})
    # check dimensions
    N = ndims(core)
    rk = collect(size(core))
    N == length(U) || error("core dimension does not match factors")
    I = zeros(Int64,N)
    for i = 1 : N
        rk[i] == size(U[i], 2) || error("core dimension does not match factors")
        I[i] = size(U[i],1)
    end
    T = Tucker(N, rk, I, copy(core), deepcopy(U))
    return T
end

"""
core = reshape(collect(1.:81.),3,3,3,3)
U = Vector{Matrix{Float64}}(4)
for i = 1 : 4
    U[i] = reshape(collect(1.:6.),2,3)
end
T = initTucker(core, U)
C = tucker2Tensor(T)
"""
function tucker2Tensor(T::tucker)
    etp = eltype(T.core)
    N = T.N
    C = initTensor(T.core)
    for i = 1 : N
        C = ttm(C, T.U[i], i)
    end
    return C
end

"""
nth mode factor matrix
"""
function nvecs(X::tensor, n::Ti, rk::Ti) where {Ti<:Int64}
    Xn = matricization(X, n)
    Y  = Xn * Xn'
    (u, s, v) = svd(Y)
    return u[:,1:rk]
end

"""
tucker decomposition by alternative minimization
"""
function tuckerAls(X::tensor, rk::Union{Vector{Ti},Ti};
                   tol=1e-4, maxIter=50, init="rand") where {Ti<:Int64}
    etp = eltype(X.D)
    N = X.N
    I = X.I
    fitchangetol = tol
    normX = vecnorm(X.D)

    if typeof(rk) <: Integer
       rk = rk*ones(Int64, N)
    end
    N == length(rk) || error("dimension of D does not match length of rank")

    # initialization factors
    U = Vector{Matrix{etp}}(N)
    if init == "rand"
       for i = 1 : N
           U[i] = randn(etp, I[i], rk[i])
       end
    elseif init == "eigvec"
       for i = 1 : N
           U[i] = nvecs(X, i, rk[i])
       end
    else
       error("unsupported initialization method")
    end

    fit = 0.0; fitold = 0.0;
    core = zeros(etp, rk...)
    for iter = 1 : maxIter
        fitold = fit
        Utilde = zeros()
        for n = 1 : N
            Utilde = ttm(X, U, -n, tflag=true)
            U[n]   = nvecs(Utilde, n, rk[n])
        end
        core = ttm(Utilde, U[N], N, tflag=true)
        normresidual = sqrt(abs(normX^2 - vecnorm(core.D)^2))
        fit = 1 - normresidual/normX
        fitchange = fit - fitold

        if fitchange < tol
           break
        end
        # println("iter $iter, fit $fit, fitchange $fitchange")
    end
    return tucker(N, I, rk, core.D, U)
end

"""
   Wrap tucker decomposition, the input is multi-dimensional array
"""
function WrapTuckerAls(dp::Array{Tv}, rk::Int64) where Tv <: AbstractFloat
    N    = ndims(dp)

    X    = tensor(N, collect(size(dp)), convert(Array{Float64,N}, dp))
    T    = tuckerAls(X, rk, tol=tol)
    dp   = tucker2Tensor(T); dp = convert(Array{Float32,1}, vec(dp.D))
    return dp
    return nothing
end

"""
   Wrap tucker decomposition, the input is a tuple
"""
function WrapTuckerAls(par::Tuple{String, Vector{Int64}, Float64})
    path = par[1]
    println("$path")
    rk   = par[2]
    tol  = par[3]
    dp   = getOnePatch(path)
    N    = ndims(dp)
    X    = tensor(N, collect(size(dp)), convert(Array{Float64,N}, dp))
    T    = tuckerAls(X, rk, tol=tol)
    dp   = tucker2Tensor(T); dp = convert(Array{Float32,1}, vec(dp.D))
    fid  = open(path, "r+"); pos = sizeof(Int32)*9;
    seek(fid, pos); write(fid, dp); close(fid);
    return nothing
end
