function GetFourierMtx(N::Int64)
    w = e^(-2.0*pi*im/N)
    F = zeros(Complex128, N, N)
    for j = 1 : N
        for i = 1 : N
            F[i,j] = w^((i-1)*(j-1))
        end
    end
    return F / sqrt(N)
end

function GetCircMtx{Tv<:AbstractFloat}(v::Union{Vector{Tv}, Vector{Complex{Tv}}})
    N = length(v)
    C = zeros(eltype(v), N, N)
    for j = 1 : N
        for i = N-j+2 : N
            C[i+j-N-1, j] = v[i]
        end
        for i = 1 : N-j+1
            C[i+j-1, j] = v[i]
        end
    end
    return C
end

function GetCircBlkMtx{Tv<:AbstractFloat}(a::Array{Tv,3})
    (n1,n2,n3) = size(a)
    cmtx = zeros(eltype(a), n1*n3, n2*n3)
    for j = 1 : n3
        for i = n3-j+2 : n3
            i1l = (i+j-n3-2)*n1+1; i1u = i1l+n1-1;
            i2l = (j-1     )*n2+1; i2u = i2l+n2-1;
            cmtx[i1l:i1u,i2l:i2u] = a[:,:,i]
        end
        for i = 1 : n3-j+1
            i1l = (i+j-2)*n1+1; i1u = i1l+n1-1;
            i2l = (j-1  )*n2+1; i2u = i2l+n2-1;
            cmtx[i1l:i1u,i2l:i2u] = a[:,:,i]
        end
    end
    return cmtx
end

"""
tensor multiplication performed in frequency domain
"""
function FTproduct{Tv<:AbstractFloat}(a::Array{Tv,3}, b::Array{Tv,3})
    (n1, n2, n3) = size(a)
    (m1, m2, m3) = size(b)
    n2 == m1 || error("size dismatch")
    n3 == m3 || error("size dismatch")
    fa = fft(a,3); fb = fft(b,3)
    fc = zeros(Complex{eltype(a)}, n1, m2, n3)
    for i = 1 : n3
        fc[:,:,i] = fa[:,:,i] * fb[:,:,i]
    end
    return real(ifft(fc,3))
end

"""
tensor multiplication in time domain
"""
function Tproduct{Tv<:AbstractFloat}(a::Array{Tv,3}, b::Array{Tv,3})
    (n1, n2, n3) = size(a)
    (m1, m2, m3) = size(b)
    n2 == m1 || error("size dismatch")
    n3 == m3 || error("size dismatch")
    Ma = GetCircBlkMtx(a)
    Vb = unfold(b)
    Vc = Ma * Vb
    return fold(Vc, n1)
end

"""
  unfold a tensor into block vector
"""
function unfold{Tv<:AbstractFloat}(b::Array{Tv,3})
    (n1,n2,n3) = size(b)
    r = zeros(eltype(b), n1*n3, n2)
    for i = 1 : n3
        indl = (i-1)*n1+1
        indu = indl + n1 - 1
        r[indl:indu,1:n2] = b[:,:,i]
    end
    return r
end

"""
  fold block vector back into tensor
"""
function fold{Tv<:AbstractFloat}(a::Matrix{Tv}, n1::Int64)
    (N, n2) = size(a)
    n3 = round(Int64, N/n1)
    b  = zeros(eltype(a), n1, n2, n3)
    for i = 1 : n3
        indl = (i-1)*n1+1
        indu = indl + n1 - 1
        b[:,:,i] = a[indl:indu,1:n2]
    end
    return b
end

function randomSample(T::Array{Float64,3}, rate::Float64)
    (n1, n2, n3) = size(T)
    PT = zeros(T)
    omega = zeros(Float64, n1, n2, n3)
    pb = 1.0 - rate
    for i3 = 1 : n3
        for i2 = 1 : n2
            for i1 = 1 : n1
                if rand() > pb
                   omega[i1,i2,i3] = 1.0
                   PT[i1,i2,i3] = T[i1,i2,i3]
                end
            end
        end
    end
    return PT, omega
end
"""
   Alternating minimization of third-order tensor, T is the observed data. P is sampling tensor
r is predefined rank, maxIter is the number of iterations.
"""
function TubalAltMin(T::Array{Float64,3},  P::Array{Float64,3}, r::Int64, maxIter::Int64)

    (n1, n2, n3) = size(T)
    X = zeros(Complex128, n1, r , n3)
    Y = zeros(Complex128, r , n2, n3)
    FT = fft(T,3)
    # initialize X from the left sigular eigenslices
    for i3 = 1 : n3
        (u, s, v) = svd(view(FT, :, :, i3))
        for i2 = 1 : r
            for i1 = 1 : n1
                X[i1,i2,i3] = u[i1,i2]
            end
        end
    end
    FP = fft(P, 3)


end

# generate synthetic data
m = 200; n = 200; k = 50; r = 10;
X = rand(m, r, k); Y = rand(r, n, k);
T = FTproduct(X, Y)
(PT, P) = randomSample(T, 0.5)

FPT = fft(PT, 3)/sqrt(k)
j = 17;
lsT = FPT[:,j,:]

FX = fft(X, 3) / sqrt(k)
FY = fft(Y, 3) / sqrt(k)
lsY = FY[:,j,:]

tmp = zeros(Complex128, m, k)
for i = 1 : k
    tmp[:,i] = FX[:,:,i] * lsY[:,i]
end

FP = fft(P,3) / sqrt(k)
lsP = FP[:,j,:]

result = zeros(Complex128, size(lsP,1), size(lsP,2))
for i = 1 : size(tmp,1)
    tube = vec(lsP[i,:])
    C = GetCircMtx(tube)
    result[i,:] = C * tmp[i,:]
end
rtrue = FPT[:,j,:]

TFX = zeros(Complex128, r, m, k)
TFY = zeros(Complex128, n, r, k)
TFPT = zeros(Complex128, n, m, k)
TFP  = zeros(Complex128, n, m, k)
for i = 1 : k
    TFX[:,:,i] = FX[:,:,i]'
    TFY[:,:,i] = FY[:,:,i]'
    TFPT[:,:,i] = FPT[:,:,i]'
    TFP[:,:,i] = FP[:,:,i]'
end
rtrue = TFPT[:,j,:]

tmp = zeros(Complex128, n, k)
for i = 1 : k
    tmp[:,i] = TFY[:,:,i] * TFX[:,j,i]
end

lsP = TFP[:,j,:]
result = zeros(Complex128, size(lsP,1), size(lsP,2))
for i = 1 : size(tmp,1)
    tube = vec(lsP[i,:])
    C = GetCircMtx(tube)
    result[i,:] = C * tmp[i,:]
end













# test diagonalize circulant matrix;
N = 1024
c = rand(N)
v = rand(N)
C = GetCircMtx(c);
F = GetFourierMtx(N);
d = diag(F * C * F')
fc = fft(c)

# test diagonalize block circulant matrix
a = reshape(collect(1.:64.),4,4,4)
Mc = GetCircBlkMtx(a)
F  = kron(GetFourierMtx(4), eye(4))

r = F * Mc * F'
rc = fft(a,3)

# test tensor product
n1 = 37; n2 = 77; n3 = 128;
m1 = 77; m2 = 59; m3 = 128;

A = rand(n1, n2, n3);
B = rand(m1, m2, m3);

C = Tproduct(A, B);
C1 = FTproduct(A, B)
vecnorm(C-C1)

# test associative property of tensor product
n1 = 37; n2 = 77; n3 = 128;
m1 = 77; m2 = 59; m3 = 128;
k1 = 59; k2 = 11; k3 = 128;
A = rand(n1, n2, n3);
B = rand(m1, m2, m3);
C = rand(k1, k2, k3);
tmp = FTproduct(A,B); R = FTproduct(tmp, C);
tmp = FTproduct(B,C); R1 = FTproduct(A, tmp);
vecnorm(R-R1)

# test tsvd
n1 = 37; n2 = 77; n3 = 128;
A = rand(n1,n2,n3);
fA = fft(A,3);
U  = zeros(Complex128, n1, n1, n3)
S  = zeros(Complex128, n1, n2, n3)
V  = zeros(Complex128, n2, n2, n3)
for i = 1 : 3
    (s, v, d) = svd(fa[:,:,i])

end








# =======
