type tensor
     N :: Int64
     I :: Array{Int64,1}
     D :: Array{Float64}
end

"""
    InitTensor(N, I, D)

create a tensor object

# Arguments
* `N`: Int64, dimensions
* `I`: Array{Int64}, size of each dimension
* `D`: Array{Float64}, multi-dimension array

# Returns
* `X`: tensor object

# Example
```julia
julia> N = 5; I = [23, 42, 13, 14, 17]; D = rand(I...)
julia> X = InitTensor(N, I, D)
```
"""
function initTensor(N, I, D)
    N==length(I) || throw(DimensionMismatch("length(I) != N"))
    for m = 1 : N
        I[m]==size(D,m) || throw(DimensionMismatch("I[m] != size(D,m)"))
    end
    return tensor(N, I, D)
end

"""
    Xn = matricization(X::tensor, n::Int64)

unfolding tensor along n dimension

# Arguments
* `X`: tensor, tensor object
* `n`: Int64, unfolding along n dimension, 1<= n <= N

# Returns
* `Xn`: Array{Float64,2}, matrix

# Example
```julia
julia> N = 3; I = [3, 4, 2]; D = reshape(collect(1:prod(I)), I...)
julia> X = initTensor(N, I, D)
julia> n = 3
julia> Xn = matricization(X, n)
```
"""
function matricization(X::tensor, n::Int64)
    N = X.N
    I = X.I
    if n == 1 # first dimension
       Xn = copy(reshape(X.D, I[1], prod(I[2:N])))
    elseif n == N
       Xn = reshape(X.D, prod(I[1:N-1]), I[N])
       Xn = Xn' # last dimension
    else
       p = collect(1:N); p[n] = 1; p[1] = n;  # permute first and nth dimension
       Xn = permutedims(X.D, p)
       Xn = reshape(Xn, I[n], round(Int64, prod(I)/I[n]))
    end
    return Xn
end

"""
  X = unmatricization(Xn, n, I)

folding an Array back to tensor

# Arguments
* `Xn`: Array{Float64}, the array to be folded
* `n`: Int64, folding along n dimension, 1<= n <= N
* `I`: Array{Int64}, size of tensor

# Returns
* `X`: tensor

# Example
```julia
julia> N = 3; I = [3, 4, 2]; D = reshape(collect(1:prod(I)), I[1], prod(I[2:N]))
julia> n = 1; X = unmatricization(D, n, I)
```
"""
function unmatricization(D::Array{Float64}, n::Int64, I::Array{Int64,1})
    if n == 1
       D = reshape(D, I...)
    elseif n == length(I)
       D = D'
       D = reshape(D, I...)
    else
       I1 = copy(I); I1[1] = I[n]; I1[n] = I[1];
       D = reshape(D, I1...)
       p  = collect(1:length(I)); p[n] = 1; p[1] = n;
       D  = permutedims(D, p)
    end
    X = tensor(length(I), I, D)
    return X
end

"""
    r = KRP(a1, a2)

Khatri-Rao product a2 ⊙ a1

# Arguments
* `a1`: Array{Float64}, it can be a 1D or 2D array
* `a2`: Array{Float64}, it can be a 1D or 2D array

# Returns
* `r`: Array{Float64}, it can be a 1D or 2D array

# Example
```julia
julia> R = 20; m1 = 77; m2 = 33;
julia> a1 = rand(m1,R); a2 = rand(m2,R);
julia> r = KRP(a1, a2)
```
"""
function KRP(a1::Array{Float64,1}, a2::Array{Float64,1})
    m1 = length(a1)
    m2 = length(a2)
    m  = m1*m2
    r  = Array{Float64}(m)
    for i2 = 1 : m2
        for i1 = 1 : m1
            idx = i1 + (i2-1)*m1
            r[idx] = a2[i2] * a1[i1]
        end
    end
    return r
end
function KRP(a1::Array{Float64,2}, a2::Array{Float64,2})
    size(a1,2) == size(a2,2) || throw(DimensionMismatch("size(a1,2) !=  size(a2,2)"))
    m1 = size(a1,1)
    m2 = size(a2,1)
    m  = m1*m2
    n  = size(a1,2)
    r  = Array{Float64}(m, n)
    r = zeros(typeof(a1[1]), m, n)
    for j = 1 : n
        for i2 = 1 : m2
            for i1 = 1 : m1
                idx = i1 + (i2-1)*m1
                r[idx,j] = a2[i2,j] * a1[i1,j]
            end
        end
    end
    return r
end

"""
    r = recursiveKRP(A)

recursive Khatri-Rao product A[N]⊙...⊙A[1]

# Arguments
* `A`: Array{Array{Float64}}, each element can be 1D or 2D array

# Returns
* `r`: Array{Float64}, it may be a 1D or 2D array

# Example
```julia
julia> I = [21, 32, 43, 24]; R = 3;
julia> N = length(I);
julia> A = Array{Array{Float64}}(N)
julia> for m = 1 : N
julia>     A[m] = rand(I[m], R)
julia> end
julia> r1 = recursiveKRP(A);
```
"""
function recursiveKRP(A::Array{Array{Float64,1},1})
     N = length(A)
     if N == 1
        return A[1]
     elseif N == 2
        return KRP(A[1], A[2])
     else
        r = KRP(A[1], A[2])
        for m = 3 : N
            r = KRP(r, A[m])
        end
        return r
     end
end
function recursiveKRP(A::Array{Array{Float64,2},1})
    N = length(A)
    if N == 1
       return A[1]
    elseif N == 2
       return KRP(A[1], A[2])
    else
       r = KRP(A[1], A[2])
       for m = 3 : N
           r = KRP(r, A[m])
       end
       return r
    end
end

"""
    l = tnorm(X)

compute the frobenius norm of tensor

# Arguments
* `X`: tensor

# Returns
* `l`: the frobenius norm of X

# Example
```julia
julia> I = [21, 32, 43]; N = 3; D = rand(I...)
julia> l = tnorm(tensor(N, I, D))
```
"""
function tnorm(X::tensor)
    sqrt(sum(X.D.^2))
end

"""
    y = ttv(X, v, dims)

compute tensor times a set of vectors

# Arguments
* `X`: tensor
* `v`: Array{Array{Float64,1}}
* `dims`: specification the multiplication on which dimension

# Returns
* `y`: Array{Float64,1}, length is prod(remdims)
"""
function ttv(X::tensor, v::Array{Array{Float64,1}}, dims::Array{Int64})
    for i = 1 : length(dims)
        1 <= dims[i] <= X.N || throw(DimensionMismatch())
        X.I[dims[i]] == length(v[i]) || throw(DimensionMismatch())
    end
    n = X.N
    remdims = setdiff(collect(1:N), dims)
    pdims = vcat(remdims, dims)
    newI  = X.I[pdims]
    des = zeros(newI...)
    if length(dims) > 1
       permutedims!(des, X.D, pdims)
    end
    for i = length(dims) : -1 : 1
        des = reshape(des, prod(newI[1:n-1]), newI[n])
        des = des * v[i]
        n = n - 1
    end
    return des
end
