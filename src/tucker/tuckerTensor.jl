type tucker{Tv<:AbstractFloat}
     N  :: Int64
     rk :: Vector{Int64}
     I  :: Vector{Int64}
     core :: Array{Tv}
     U  :: Vector{Matrix{Tv}}
end

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

function nvecs{Tv<:AbstractFloat, Ti<:Int64}(X::tensor{Tv}, n::Ti, rk::Ti)
    Xn = matricization(X, n)
    Y  = Xn * Xn'
    (u, s, v) = svd(Y)
    return u[:,1:rk]
end

function tuckerAls{Tv<:AbstractFloat, Ti<:Int64}(X::tensor{Tv}, rk::Union{Vector{Ti},Ti};
                   tol=1e-4, maxIter=50, init="rand")
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




# construct a low tucker-rank 4D tensor
# N = 4
# rk = [11, 13, 15, 8];
# I  = [42, 22, 39, 88];
# core1 = randn(rk...)
# U1  = Vector{Matrix{Float64}}(N)
# for i = 1 : N
#     U1[i] = randn(I[i], rk[i])
# end
# T = tucker(N, I, rk, core1, U1)
# X = tucker2Tensor(T)
# dobs = SeisAddNoise(X.D, 2.0)
# Xn = initTensor(dobs)
#
# T1 = tuckerAls(Xn, rk);
# X1 = tucker2Tensor(T1);
#
# d = X.D; dr = X1.D
#
# SeisPlot(d[1,1,:,:], pclip=100)
# SeisPlot(dr[1,1,:,:], pclip=100)
# SeisPlot(dobs[1,1,:,:], pclip=100)
# vecnorm(d[1,1,:,:]-dr[1,1,:,:])
