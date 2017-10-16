"""
tensor multiplication performed in frequency domain, support order-N real tensor
"""
function tprod{Tv<:AbstractFloat}(A::Array{Tv}, B::Array{Tv})
    etp= eltype(A)
    sa = collect(size(A))
    sb = collect(size(B))
    sa[2] == sb[1] || error("size dismatch")
    sa[3:end] == sb[3:end] || error("size dismatch")

    Ac = convert(Array{Complex{etp}}, A);
    Bc = convert(Array{Complex{etp}}, B);
    for i = 3 : length(sa)
        fft!(Ac, i)
        fft!(Bc, i)
    end

    nfaces = prod(sa[3:end])
    Ac = reshape(Ac, sa[1], sa[2], nfaces)
    Bc = reshape(Bc, sb[1], sb[2], nfaces)
    C  = zeros(Complex{etp}, sa[1], sb[2], nfaces)
    for i = 1 : nfaces
        A_mul_B!(view(C,:,:,i), view(Ac,:,:,i), view(Bc,:,:,i))
    end
    sc = copy(sa); sc[2] = sb[2];
    C = reshape(C, sc...)

    for i = 3 : length(sa)
        ifft!(C, i)
    end
    return real(C)
end

"""
tensor svd for N-order real tensor
"""
function tsvd{Tv<:AbstractFloat}(A::Array{Tv})
    etp = eltype(A)
    sa = size(A)
    N  = length(sa)
    n1 = sa[1]; n2 = sa[2]; n3 = prod(sa[3:end]);

    # for 3rd order tensor, use conjugate symmetric trick to save computation
    k = min(n1, n2)
    if N == 3
       Ac = convert(Array{Complex{etp}}, A)
       U  = zeros(Complex{etp}, n1, k, n3)
       S  = zeros(etp, k, n3)
       V  = zeros(Complex{etp}, n2, k, n3)
       fft!(Ac, 3)
       nh = floor(Int64, n3/2)
       for i = 1 : nh+1
           (U[:,:,i], S[:,i], V[:,:,i]) = svd(view(Ac, :, :, i))
       end
       if mod(n3, 2) == 0
          for i = 2 : nh
              idx = n3 - i + 2
              U[:,:,idx] = conj(U[:,:,i])
              S[:  ,idx] =      S[:  ,i]
              V[:,:,idx] = conj(V[:,:,i])
          end
       else
          for i = 2 : nh+1
              idx = n3 - i + 2
              U[:,:,idx] = conj(U[:,:,i])
              S[:  ,idx] =      S[:  ,i]
              V[:,:,idx] = conj(V[:,:,i])
          end
       end
    end
    ifft!(U, 3); U = real(U);
    ifft!(V, 3); V = real(V);
    return U, S, V
end

function ttranspose{Tv<:AbstractFloat}(A::Array{Tv})
    sa = size(A)
    N  = length(sa)
    if N == 3
       At = zeros(sa[2], sa[1], sa[3])
       nh = floor(Int64, sa[3]/2)
       if mod(sa[3], 2) == 0
          for i = 2 : nh
              idx = n3 - i + 2
              At[:,:,i]   = (A[:,:,idx])'
              At[:,:,idx] = (A[:,:,i]  )'
          end
          At[:,:,nh+1] = (A[:,:,nh+1])'
       else
          for i = 2 : nh+1
              idx = n3 - i + 2
              U[:,:,idx] = conj(U[:,:,i])
              S[:  ,idx] =      S[:  ,i]
              V[:,:,idx] = conj(V[:,:,i])
          end
       end
    end
end
# ////////////////////////////
# n1 = 15; n2 = 15; n3 =17; n4 =18;
# n1b = n2; n2b = n2;
# A = randn(n1 , n2 , n3, n4);
# B = randn(n1b, n2b, n3, n4);
# @time C = tprod(A, B);
# path = "/Users/wenlei/Desktop/test.bin"
# fid = open(path, "w"); write(fid, n1, n2, n3, n4);
# write(fid, vec(A)); write(fid, vec(B)); write(fid, vec(C)); close(fid);

# test in matlab;
# fid = fopen('test.bin', 'r')
# n = fread(fid, 4, 'int64')'
# A = fread(fid, prod(n), 'float64'); A = reshape(A, n);
# B = fread(fid, prod(n), 'float64'); B = reshape(B, n);
# C = fread(fid, prod(n), 'float64'); C = reshape(C, n);
# C1 = tprod(A, B);
# norm(C1(:)-C(:));

# ////////////////////////////
n1 = 13; n2 = 29; n3 =128;
A = randn(n1 , n2 , n3);
(U, S, V) = tsvd(A);
k = size(S,1)
Sm= zeros(k, k, n3);
for i = 1 : n3
    Sm[:,:,i] = diagm(S[:,i])
end
































#
