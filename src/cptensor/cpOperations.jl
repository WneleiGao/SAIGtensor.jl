type cptensor
     N :: Int64
     I :: Array{Int64,1}
     R :: Int64
     lambda :: Array{Float64,1}
     A :: Array{Array{Float64,2}}
end

"""
    initCptensor(N, I, R, lambda, A)

create a cptensor object, right to test, it is a great work. fantastic packages.

# Arguments
* `N`: Int64, dimensions
* `I`: Array{Int64}, size of each dimension
* `R`: Int64, rank of tensor
* `lambda`: eigenvalue of tensor
* `A`: Array{Array{Float64,2}}, basis

# Returns
* `X`: CP tensor

# Example
```julia
julia> N = 5; I = [23, 42, 13, 14, 17];
julia> R = 10; lambda = rand(N);
julia> A = Array{Array{Float64,2}}(N)
julia> for m = 1 : N
julia>     A[m] = rand(I[m], R)
julia> end
julia> X = initCptensor(N, I, R, lambda, A);
```
"""
function initCptensor(N::Int64, I::Array{Int64}, R::Int64, lambda::Array{Float64,1}, A::Array{Array{Float64,2}})
    N == length(I) || throw(DimensionMismatch("length(I) != N"))
    for m = 1 : N
        I[m] == size(A[m],1) || throw(DimensionMismatch("I[m] != size(A[m],1)"))
        R == size(A[m],2) || throw(DimensionMismatch("R != size(A[m],2)"))
    end
    return cptensor(N, I, R, lambda, A)
end

"""
    l = cpnorm(X)

efficient way of computing the Frobenius norm of cp tensor

# Arguments
* `X`: cptensor

# Returns
* `l`: the frobenius norm of X

# Example
```julia
julia> I = [21, 32, 43, 24]; R = 3;
julia> N = length(I);
julia> A = Array{Array{Float64}}(N)
julia> for m = 1 : N
julia>     A[m] = rand(I[m], R)
julia> end
julia> X = initCptensor(N, I, R, lambda, A); l = cpnorm(X);
```
"""
function cpnorm(X::cptensor)
    coef = X.lambda * X.lambda'
    N = X.N
    for m = 1 : N
        coef = coef .* (X.A[m]'*X.A[m])
    end
    return sqrt(sum(coef))
end
function cpnorm(lambda::Array{Float64,1}, A::Array{Array{Float64,2},1})
    coef = lambda * lambda'
    N = length(A)
    for m = 1 : N
        coef = coef .* (A[m]'*A[m])
    end
    return sqrt(sum(coef))
end

function innerprod(X::cptensor, Y::tensor)
    X.N == Y.N || throw(DimensionMismatch("X.N != Y.N"))
    N = X.N
    for m = 1 : N
        if X.I[m] != Y.I[m]
           throw(DimensionMismatch("X.I[m] != Y.I[m]"))
        end
    end
    I = X.I
    tmp = 0.
    for r = 1 : X.R
        T = reshape(Y.D, prod(I[1:N-1]), I[N]) * X.A[N][:,r]
        for m = N-1 : -1 : 1
            T = reshape(T, prod(I[1:m-1]), I[m]) * X.A[m][:,r]
        end
        tmp = tmp + T[1] * X.lambda[r]
    end
    return tmp
end
function innerprod(lambda::Array{Float64,1}, A::Array{Array{Float64,2},1}, Y::tensor)
    R = length(lambda)
    N = length(A)
    I = zeros(Int64,N)
    for i = 1 : N
        I[i] = size(A[i],1)
    end
    tmp = 0.
    for r = 1 : R
        T = reshape(Y.D, prod(I[1:N-1]), I[N]) * A[N][:,r]
        for m = N-1 : -1 : 1
            T = reshape(T, prod(I[1:m-1]), I[m]) * A[m][:,r]
        end
        tmp = tmp + T[1] * lambda[r]
    end
    return tmp
end

"""
    f = fitness(X, Y)

compute the relative fitness of lowrank cptensor to a full tensor

# Arguments
* `X`: cptensor
* `Y`: tensor

# Returns
* `f`: relative fitness
"""
function fitness(X::cptensor, Y::tensor)
    lx = cpnorm(X);
    ly = tnorm(Y)
    ipt = innerprod(X, Y)
    fit = sqrt(lx^2 + ly^2 -2*ipt) / ly
    return fit
end

"""
    updateAn!(A, AtA, Gn, n)

update the nth entries in A and AtA

# Arguments
* `A`  : Array{Array{Float64,2}}, Array of factor matrix
* `AtA`: Array{Array{Float64,2}}, Array of A[n]'×A[n]
* `Gn` : Array{Float64,2}, Xn×(⊗A[k]) for all k!=n
* `n`  : Int64, nth factor
"""
function updateAn!(lambda::Array{Float64,1}, A::Array{Array{Float64,2}}, AtA::Array{Array{Float64,2}}, Gn::Array{Float64,2}, n::Int64)
    N = length(A); R = size(A[1],2);
    L = ones(R, R); In = size(A[n],1);
    I = setdiff(collect(1:N), [n])
    for i in I
        L = L .* AtA[i]
    end
    tmp = Gn / L
    for i2 = 1 : R
        # lambda[i2] = sqrt(sum(tmp[:,i2].^2))
        lambda[i2] = maximum(abs(tmp))
        for i1 = 1 : In
            tmp[i1,i2] = tmp[i1,i2] / lambda[i2]
        end
    end
    A[n] = tmp
    AtA[n] = A[n]'*A[n]
    return nothing
end

"""
    cptfirst(I)

decide which factor should be update first, the updating order is
ns, ns-1:-1:1, ns+1, ns+2:N

# Arguments
* `I`  : Array{Int64,1}, dimensions of tensor
"""
function cptfirst(I::Array{Int64})
    N = length(I)
    ns = 1; Jn = prod(I[1:ns]); Kn = prod(I[ns+1:N]);
    while Jn <= Kn
        ns = ns + 1
        Jn = Jn * I[ns]
        Kn = round(Int64, Kn/I[ns])
    end
    return ns - 1
end

"""
    cptRight(X, A, n)

compute X × A[n+1] × A[n+2] × ... × A[N]

# Arguments
* `X`: tensor
* `A`: Array{Array{Float64,2}}, Array of factor matrix
* `n`: Int64

# Returns
* `Rn` : Array{Float64,2}
"""
function cptRight(X::tensor, A::Array{Array{Float64,2}}, n::Int64)
    N  = X.N
    n < N ? nothing : error("n must less than N")
    B = A[n+1:N]
    I = X.I
    tmp = recursiveKRP(B)
    szr = prod(I[1:n])
    szc = prod(I[n+1:N])
    Xn  = reshape(X.D, szr, szc)
    Rn  = Xn * tmp
    return Rn
end

"When applied matrix, using the efficient way to compute them based on the
property X × A[n+1] × A[n+2] × ... × A[N] = (X × A[n+2] × ... × A[N]) × A[n+1]
* `Rnp1`: Array{Float64,2} with size prod(I[1:n+1]) * R
* `Anp1`: Array{Array{Float64,2}} with size I[n+1] * R
* `I`   : Array{Int64}
* `n`   : Int64                                                               "
function cptRight(Rnp1::Array{Float64,2}, Anp1::Array{Float64,2}, I::Array{Int64,1}, n::Int64)
    R = size(Rnp1, 2)
    szr = prod(I[1:n])
    szc = I[n+1]
    Rnew = zeros(szr, R)
    tmp  = zeros(szr)
    for i = 1 : R
        T = reshape(Rnp1[:,i], szr, szc)
        A_mul_B!(tmp, T, Anp1[:,i])
        for k = 1 : szr
            Rnew[k,i] = tmp[k]
        end
    end
    return Rnew
end

"""
    cptLeft(X, A, n)

compute X × A[1] × A[2] × ... × A[n-1]

# Arguments
* `X`  : tensor
* `A`: Array{Array{Float64,2}}, Array of factor matrix
* `n`  : Int64

# Returns
* `Ln` : Array{Float64,2}
"""
function cptLeft(X::tensor, A::Array{Array{Float64,2}}, n::Int64)
    N = X.N
    n > 1 ? nothing : error("n must be larger than 1")
    I = X.I
    R = size(A[1],2)
    B = A[1:n-1]
    tmp = recursiveKRP(B)
    szr = prod(I[1:n-1])
    szc = prod(I[n:N])
    Xn  = reshape(X.D, szr, szc)
    Ln  = Array{Float64}(szc, R)
    At_mul_B!(Ln, Xn, tmp)
    return Ln
end

"efficient way to compute Ln based on the property
X × A[1] ... ×A[n-1] × A[n] = (X × A[1] ... ×A[n-1]) × A[n]
* `Lnm1`: Array{Float64,2} with size prod(I[1:n+1]) * R
* `Anm1`: Array{Array{Float64,2}} with size I[n+1] * R
* `I`   : Array{Int64}
* `n`   : Int64                                             "
function cptLeft(Lnm1::Array{Float64,2}, Anm1::Array{Float64,2}, I::Array{Int64,1}, n::Int64)
    R = size(Lnm1, 2)
    szr = I[n-1]
    szc = prod(I[n:end])
    Lnew = zeros(szc, R)
    tmp  = zeros(szc)
    for i = 1 : R
        T = reshape(Lnm1[:,i], szr, szc)
        At_mul_B!(tmp, T, Anm1[:,i])
        for k = 1 : szc
            Lnew[k,i] = tmp[k]
        end
    end
    return Lnew
end

"""
    cp_gradient(Z, A, I, n)

when dir="L", it compute X × A[-n] from X × A[n+1]...A[N]
when dir="R", it compute X × A[-n] from X × A[1]...A[n-1]

# Arguments
* `Z`: Array{Float64,2}
* `A`: Array{Array{Float64,2}}, Array of factor matrix
* `I`: Array{Int64,1}, dimension of tensor
* `n`: Int64.

# Returns
* `G` : Array{Float64,2} with size I[n] × R
"""
function cp_gradient(Z::Array{Float64,2}, A::Array{Array{Float64,2}}, I::Array{Int64,1}, n::Int64; dir="L")
    R = size(Z, 2)
    N = length(A)
    G = zeros(I[n], R)
    tv = zeros(I[n])
    if dir == "L"
       if n == 1
          G = copy(Z)
       else
          szr = prod(I[1:n-1])
          szc = I[n]
          B = A[1:n-1]
          tmp = recursiveKRP(B)
          for i = 1 : R
              S = reshape(Z[:,i], szr, szc)
              At_mul_B!(tv, S, tmp[:,i])
              for j = 1 : I[n]
                  G[j,i] = tv[j]
              end
          end
       end
    elseif dir == "R"
       if n == N
          G = copy(Z)
       else
          szr = I[n]
          szc = prod(I[n+1:end])
          B = A[n+1:end]
          tmp = recursiveKRP(B)
          for i = 1 : R
              S = reshape(Z[:,i], szr, szc)
              A_mul_B!(tv, S, tmp[:,i])
              for j = 1 : I[n]
                  G[j,i] = tv[j]
              end
          end
       end
    end
    return G
end

function cpals(X::tensor, R::Int64; maxiter::Int64=50, fitchangetol=1e-6, pflag=false)
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
    println("CP_ALS")  # alternating Least square
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

function cpals_simplify(X::tensor, R::Int64; maxiter::Int64=50)
    N = X.N; I = X.I;
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
    for iter = 1 : maxiter
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
    end
    cpnormalize!(lambda, A)
    return lambda, A
end

"""
    ttf(X, A, n)

efficiet way computing tensor times all factor matrix except nth factor

# Arguments
* `X`: tensor
* `A`: Array{Array{Float64,2}}, Array of factor matrix
* `n`: Int64, nth factor

# Returns
* `G` : Array{Float64,2} with size I[n] × R
"""
function ttf(X::tensor, A::Array{Array{Float64,2}}, n::Int64)
    I = X.I
    ns = cptfirst(I)
    if n <= ns
       Rn = cptRight(X, A, n);
       G  = cp_gradient(Rn, A, I, n, dir="L")
    else
       Ln = cptLeft(X, A, n);
       G  = cp_gradient(Ln, A, I, n, dir="R")
    end
    return G
end

"""
    ttf_slow(X, A, n)

tensor times all factor matrix except nth factor
the usage of this method is deprecated.

# Arguments
* `X`: tensor
* `A`: Array{Array{Float64,2}}, Array of factor matrix
* `n`: Int64, nth factor

# Returns
* `G` : Array{Float64,2} with size I[n] × R
"""
function ttf_slow(X::tensor, A::Array{Array{Float64,2},1}, n)
    R = size(A[1],2)
    N = X.N;
    I = X.I;
    G = zeros(X.I[n], R)
    if n == 1
       G = reshape(X.D, X.I[1], prod(X.I[2:end]))
       tmp = recursiveKRP(A[2:end])
       G = G * tmp
    elseif n == N
       G = reshape(X.D, prod(X.I[1:end-1]), X.I[end])
       tmp = recursiveKRP(A[1:end-1])
       G = G' * tmp
    else
       X1 = matricization(X, 1)
       X1 = X1' * A[1]
       Irem = setdiff(collect(1:N), [1, n])
       for r = 1 : R
           I = copy(X.I); I[1] = 1;
           T = unmatricization(X1[:,r], 1, I)
           for m in Irem
               T = matricization(T, m)
               T = T' * A[m][:,r]
               I[m] = 1
               T = unmatricization(T, m, I)
           end
           G[:,r] = vec(T.D)
       end
    end
    return G
end

"""
    cpnormalize!(lambda, A)

normalize factor matrix for CP decomposition and add weight to lambda

# Arguments
* `lambda`: Array{Float64,1}, non-sorted weight
* `A`: Array{Array{Float64,2}}, Array of factor matrix
"""
function cpnormalize!(lambda::Array{Float64,1}, A::Array{Array{Float64,2},1})
    R = length(lambda)
    N = length(A)
    for r = 1 : R
        for n = 1 : N
            tmp = norm(A[n][:,r])
            if tmp > 0
               for i = 1 : size(A[n],1)
                   A[n][i,r] = A[n][i,r] / tmp
               end
            end
            lambda[r] = lambda[r] * tmp
        end
        if lambda[r] < 0
           lambda[r] = -lambda[r]
           for i = 1 : size(A[1],1)
               A[1][i,r] = -A[1][i,r]
           end
        end
    end
    return nothing
end

"""
    cp2tensor(X)

convert 4D cptensor to tensor
"""
function cp2tensor(X::cptensor)
    N = X.N
    I = X.I
    R = X.R
    lambda = X.lambda
    T = zeros(I...)
    if N == 4
       for r = 1 : R
           for i4 = 1 : I[4]
               for i3 = 1 : I[3]
                   for i2 = 1 : I[2]
                       for i1 = 1 : I[1]
                           T[i1,i2,i3,i4] = T[i1,i2,i3,i4] + lambda[r]*X.A[4][i4,r]*X.A[3][i3,r]*X.A[2][i2,r]*X.A[1][i1,r]
                       end
                   end
               end
           end
       end
    elseif N == 3
       for r = 1 : R
           for i3 = 1 : I[3]
               for i2 = 1 : I[2]
                   for i1 = 1 : I[1]
                       T[i1,i2,i3] = T[i1,i2,i3] + lambda[r]*X.A[3][i3,r]*X.A[2][i2,r]*X.A[1][i1,r]
                   end
               end
           end
       end
    end
    return tensor(N, I, T)
end

function synthesis(lambda::Vector{Float64}, A::Array{Array{Float64,2},1})

end
