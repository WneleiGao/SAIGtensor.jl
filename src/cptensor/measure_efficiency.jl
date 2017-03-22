function randCpAls_fit(X::tensor, R::Int64; maxiter::Int64=50, fitchangetol=1e-5)
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
    his = zeros(maxiter)
    for iter = 1 : maxiter
        fitold = fit
        # updating ns factor
        for n = 1 : N
            randSpl!(Zs, Xst, A, X, ns, n)
            updateAn!(lambda, A, Zs, Xst[n], n)
        end
        # check the fitness
        normresidual = sqrt(normX^2 + cpnorm(lambda,A)^2 - 2*innerprod(lambda, A, X))
        fit = 1. - normresidual/normX
        his[iter] = 1. - fit
    end
    cpnormalize!(lambda, A)
    return his
end

function randCpAls_time(X::tensor, R::Int64; maxiter::Int64=50, fitchangetol=1e-5)
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
    his = zeros(maxiter)
    for iter = 1 : maxiter
        tic()
        for n = 1 : N
            randSpl!(Zs, Xst, A, X, ns, n)
            updateAn!(lambda, A, Zs, Xst[n], n)
        end
        his[iter] = toq()
    end
    return his
end

function cpals_fit(X::tensor, R::Int64; maxiter::Int64=50)
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
    AtA = Array{Array{Float64,2}}(N)
    for m = 1 : N
        AtA[m] = At_mul_B(A[m], A[m])
    end
    # determine the updating order
    ns = cptfirst(I)
    his = zeros(maxiter)
    for iter = 1 : maxiter
        fitold = fit
        # updating ns factor
        Rn = cptRight(X, A, ns)
        Gn = cp_gradient(Rn, A, I, ns, dir="L")
        updateAn!(lambda, A, AtA, Gn, ns)
        # updating ns-1:-1:1 factors
        for m = ns-1 : -1 : 1
            Rn = cptRight(Rn, A[m+1], I, m)
            Gn = cp_gradient(Rn, A, I, m, dir="L")
            updateAn!(lambda, A, AtA, Gn, m)
        end
        # updating ns+1 factor
        Rn = cptLeft(X, A, ns+1)
        Gn = cp_gradient(Rn, A, I, ns+1, dir="R")
        updateAn!(lambda, A, AtA, Gn, ns+1)
        # updating ns+2:N factors
        for m = ns+2: N
            Rn = cptLeft(Rn, A[m-1], I, m)
            Gn = cp_gradient(Rn, A, I, m, dir="R")
            updateAn!(lambda, A, AtA, Gn, m)
        end
        # check the fitness
        normresidual = sqrt(normX^2 + cpnorm(lambda,A)^2 - 2*innerprod(lambda, A, X))
        fit = 1 - normresidual/normX
        his[iter] = 1.- fit
        println("$iter, $fit")
    end
    return his
end

function cpals_time(X::tensor, R::Int64; maxiter::Int64=50)
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
    AtA = Array{Array{Float64,2}}(N)
    for m = 1 : N
        AtA[m] = At_mul_B(A[m], A[m])
    end
    # determine the updating order
    ns = cptfirst(I)
    his = zeros(maxiter)
    for iter = 1 : maxiter
        tic()
        # updating ns factor
        Rn = cptRight(X, A, ns)
        Gn = cp_gradient(Rn, A, I, ns, dir="L")
        updateAn!(lambda, A, AtA, Gn, ns)
        # updating ns-1:-1:1 factors
        for m = ns-1 : -1 : 1
            Rn = cptRight(Rn, A[m+1], I, m)
            Gn = cp_gradient(Rn, A, I, m, dir="L")
            updateAn!(lambda, A, AtA, Gn, m)
        end
        # updating ns+1 factor
        Rn = cptLeft(X, A, ns+1)
        Gn = cp_gradient(Rn, A, I, ns+1, dir="R")
        updateAn!(lambda, A, AtA, Gn, ns+1)
        # updating ns+2:N factors
        for m = ns+2: N
            Rn = cptLeft(Rn, A[m-1], I, m)
            Gn = cp_gradient(Rn, A, I, m, dir="R")
            updateAn!(lambda, A, AtA, Gn, m)
        end
        his[iter] = toq()
    end
    return his
end
