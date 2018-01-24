# ==============================================================================
#                          test for 3D MSSA
# ==============================================================================
using PyPlot, Seismic
(d0, ext)=SeisLinearEvents(nt= 250, dt=0.008, nx1=32, nx2=32, dx1=30, dx2=30,
                         tau=[0.4, 0.8, 1.2, 1.6],
                         p1=[ 0.0003,  0.0004, -0.0004, 0.],
                         p2=[-0.0002, -0.0004,  0.0001, 0.],
                         p3=[0., 0., -0., -0.],
                         p4=[0., 0., -0., -0.],
                         amp=[1., 1.2, -0.5, -0.9])

snr = 1.0; d = SeisAddNoise(d0, snr, L=5);
(nt, n2, n3) = size(d)
df = fft(d, 1);

ncol2 = floor(Int64, n2/2) + 1
ncol3 = floor(Int64, n3/2) + 1
nrow2 = n2 - ncol2 + 1
nrow3 = n3 - ncol3 + 1

# create full Hankel matrix
f0 = 20.0; dt = 0.008;
iw = floor(Int64, f0*dt*nt) + 1; sr = 0.3;
(rowIdx, colIdx, dobs) = FormHankel2D(df, iw, sr, nrow2, ncol2, nrow3, ncol3);

# sampling part of the matrix
rk=30; maxIter = 100; nrow = nrow2 * nrow3; ncol = ncol2 * ncol3;
mu = 1.0;  (count, w1 ) = L1MC(dobs, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 10.0; (count, w10) = L1MC(dobs, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 20.0; (count, w20) = L1MC(dobs, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 30.0; (count, w30) = L1MC(dobs, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 40.0; (count, w40) = L1MC(dobs, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 50.0; (count, w50) = L1MC(dobs, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 60.0; (count, w60) = L1MC(dobs, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 70.0; (count, w70) = L1MC(dobs, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 80.0; (count, w80) = L1MC(dobs, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 90.0; (count, w90) = L1MC(dobs, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu =100.0; (count, w100)= L1MC(dobs, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);

# plot(collect(1:10), abs.(diff((abs.(w[1:11])))));
plot(collect(1:30), sort(abs.(w1[1:30]), rev=true), label="mu=1")
plot(collect(1:30), sort(abs.(w10[1:30]), rev=true), label="mu=10")
plot(collect(1:30), sort(abs.(w20[1:30]), rev=true), label="mu=20")
plot(collect(1:30), sort(abs.(w30[1:30]), rev=true), label="mu=30")
plot(collect(1:30), sort(abs.(w40[1:30]), rev=true), label="mu=40")
plot(collect(1:30), sort(abs.(w50[1:30]), rev=true), label="mu=50")
plot(collect(1:30), sort(abs.(w60[1:30]), rev=true), label="mu=60")
plot(collect(1:30), sort(abs.(w70[1:30]), rev=true), label="mu=70")
plot(collect(1:30), sort(abs.(w80[1:30]), rev=true), label="mu=80")
plot(collect(1:30), sort(abs.(w90[1:30]), rev=true), label="mu=90")
plot(collect(1:30), sort(abs.(w100[1:30]), rev=true), label="mu=100")

plot(collect(2:30), abs.(diff(sort(abs.(w1[1:30]), rev=true))), label="mu=1")
plot(collect(2:30), abs.(diff(sort(abs.(w10[1:30]), rev=true))), label="mu=10")
plot(collect(2:30), abs.(diff(sort(abs.(w20[1:30]), rev=true))), label="mu=20")
plot(collect(2:30), abs.(diff(sort(abs.(w30[1:30]), rev=true))), label="mu=30")
plot(collect(2:30), abs.(diff(sort(abs.(w40[1:30]), rev=true))), label="mu=40")
plot(collect(2:30), abs.(diff(sort(abs.(w50[1:30]), rev=true))), label="mu=50")
plot(collect(2:30), abs.(diff(sort(abs.(w60[1:30]), rev=true))), label="mu=60")
plot(collect(2:30), abs.(diff(sort(abs.(w70[1:30]), rev=true))), label="mu=70")
plot(collect(2:30), abs.(diff(sort(abs.(w80[1:30]), rev=true))), label="mu=80")
plot(collect(2:30), abs.(diff(sort(abs.(w90[1:30]), rev=true))), label="mu=90")
plot(collect(2:30), abs.(diff(sort(abs.(w100[1:30]), rev=true))), label="mu=100")

dobs = zeros(d)
for j = 1 : n3
    for i = 1 : n2
        if rand() < sr
           dobs[:,i,j] = d[:,i,j]
        end
    end
end

tmp="/Users/wenlei/Desktop/d0.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(d0))); close(tmpfid);
tmp="/Users/wenlei/Desktop/dn.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(d ))); close(tmpfid);
tmp="/Users/wenlei/Desktop/dobs.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dobs))); close(tmpfid);

pscube < /Users/wenlei/Desktop/d0.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=32 d2num=8 f2num=8 label2="X" n3=32 d3num=8 f3num=8 label3="Y" xbox=0.5 ybox=0.5 perc=0.96 > /Users/wenlei/Desktop/d0.eps
pscube < /Users/wenlei/Desktop/dn.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=32 d2num=8 f2num=8 label2="X" n3=32 d3num=8 f3num=8 label3="Y" xbox=0.5 ybox=0.5 perc=0.96 > /Users/wenlei/Desktop/dn.eps
pscube < /Users/wenlei/Desktop/dobs.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=32 d2num=8 f2num=8 label2="X" n3=32 d3num=8 f3num=8 label3="Y" xbox=0.5 ybox=0.5 perc=0.96 > /Users/wenlei/Desktop/dobs.eps
psimage < /Users/wenlei/Desktop/tscl.bin height=3.0 width=2.0 labelsize=12 n1=180 d1=1 d1num=40 f1num=40 label1="Crossline" n2=120 d2num=20 f2num=20 label2="Inline" xbox=0.5 ybox=0.5 bclip=-14387.4 wclip=13873.6 > /Users/wenlei/Desktop/tscl.eps

# ==============================================================================
#                         brute force to determine turning point (works)
# ==============================================================================
w = sort(abs.(w10[1:30]), rev=true);
function DetectTurningPoint(w::Vector{Float64})
    n = length(w)
    cost = zeros(n-2)
    for k = 2 : n-2
          b  = w[1:k]
          A  = hcat(collect(1:k), ones(k))
          lc = (A'*A) \ (A'*b)

          # fit the right line
          b  = w[k+1:n]
          A  = hcat(collect(k+1:n), ones(n-k))
          rc = (A'*A) \ (A'*b)
          for i = 1 : k
              cost[k] = cost[k] + (w[i]-lc[1]*i - lc[2])^2
          end
          for i = k+1 : n
              cost[k] = cost[k] + (w[i]-rc[1]*i - rc[2])^2
          end
     end
end

# ==============================================================================
#                          test for real data
# ==============================================================================
M = 377; N = 412; R= 27;
A = randn(M, R); B = randn(N, R);
C = A * B'; C = SeisAddNoise(C, 1.5, L=5);
SR = 0.35;
rowIdx = Int64[];
colIdx = Int64[];
d      = Float64[];
for j = 1 : N
    for i = 1 : M
        if rand() < SR
           push!(rowIdx, i)
           push!(colIdx, j)
           push!(d, C[i,j])
        end
    end
end

rk = 100; maxIter = 100; mu = 50.0;
(count, w) = L1MC(d, rowIdx, colIdx, M, N, rk, mu, maxIter);
plot(abs.(diff(abs.(w))));
plot(abs.(w))

# ==============================================================================
#                          test for complex data
# ==============================================================================
M = 377; N = 412; R= 27;
A = randn(M,R)+im*randn(M,R);
B = randn(N,R)+im*randn(N,R);
C = A * B';
SR = 0.35;
rowIdx = Int64[];
colIdx = Int64[];
d      = Complex128[];
for j = 1 : N
    for i = 1 : M
        if rand() < SR
           push!(rowIdx, i)
           push!(colIdx, j)
           push!(d, C[i,j])
        end
    end
end
rk=100; mu = 50.0; maxIter = 100;

# ==============================================================================
#                          test for complex data
# ==============================================================================
using PyPlot, Seismic
# generate clean data with 3 linear events
f0 = 20.0; dt = 0.004; tmax= 1.2; t=collect(0.0:dt:tmax); nt = length(t);
tau = [0.1  ,  0.4   , 0.55    ];
p   = [0.001, -0.0004, 0.000001]; nx =60; dx = 10.0; ph=zeros(length(p));
amp = [1.0  , -1.0   , 1.0     ];
(dobs, extent) = SeisLinearEvents(dt=dt, nt=nt, dx1=10, nx1=nx, tau=tau, p1=p, p2=ph, p3=ph, p4=ph, f0=f0, amp=amp)
snr = 1.5; dobs = SeisAddNoise(dobs, snr, L=5);

# add band-liminited gaussian noise to data, snr is defined as the ratio of maximum amplitude
(nt, nx) = size(dobs)
nf = nextpow2(nt);
nf = nf >= 512 ? nf : 512;
if nf > nt
   df = vcat(convert(Matrix{Complex{Float64}},dobs), zeros(Complex{Float64},nf-nt,nx))
else
   df = convert(Matrix{Complex{Float64}},dobs)
end
fft!(df, 1)  #in-place fft
ncol = floor(Int64, nx/2) + 1
nrow = nx - ncol + 1

iw = 40; SR=0.6
(rowIdx, colIdx, dobs) = Hankel(nrow, ncol, df, iw, SR)
rk=30; mu = 30.0; maxIter = 100;
(count, w) = L1MC(dobs, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
plot(abs.(diff((abs.(w)))));
plot(abs.(w))


# ==============================================================================
#                          Float64 cases
# ==============================================================================
function FormHankel2D(df::Array{Complex{Tv}, 3}, iw::Ti, SR::Tv,
                      nrow2::Ti, ncol2::Ti, nrow3::Ti, ncol3::Ti) where {Tv<:Float64, Ti<:Int64}
     nx1 = nrow2 + ncol2 - 1
     nx2 = nrow3 + ncol3 - 1
     idx1 = Int64[]
     idx2 = Int64[]
     rowIdx = Int64[]
     colIdx = Int64[]
     dobs   = Complex128[]
     for j = 1 : nx2
         for i = 1 : nx1
             if rand() < SR
                push!(idx1, i)
                push!(idx2, j)
             end
         end
     end
     for j = 1 : ncol3
         for i = 1 : nrow3
             tmp2 = i+j-1
             rnum = (i-1)*nrow2
             cnum = (j-1)*ncol2
             for k = 1 : ncol2
                 for l = 1 : nrow2
                     tmp1 = k + l - 1
                     for m = 1 : length(idx1)
                         if tmp1 == idx1[m] && tmp2 == idx2[m]
                            irh  = rnum + l
                            irc  = cnum + k
                            push!(rowIdx, irh)
                            push!(colIdx, irc)
                            push!(dobs, df[iw, tmp1, tmp2])
                         end
                     end
                 end
             end
         end
     end
     return rowIdx, colIdx, dobs
end


function Hankel(nrow::Ti, ncol::Ti, df::Matrix{Complex{Tv}}, iw::Ti, SR::Tv) where {Tv<:Float64, Ti<:Int64}
    nx = nrow + ncol - 1
    H = zeros(Complex{Float64}, nrow, ncol)
    idx = Int64[];
    for i = 1 : nx
        if rand() < SR
           push!(idx, i)
        end
    end
    rowIdx = Int64[]
    colIdx = Int64[]
    dt     = Complex128[]
    for i2 = 1 : ncol
        for i1 = 1 : nrow
            tmp= i1 + i2 - 1
            for k = 1 : length(idx)
                if tmp == idx[k]
                   push!(rowIdx, i1)
                   push!(colIdx, i2)
                   push!(dt, df[iw, tmp])
                end
            end
        end
    end
    return rowIdx, colIdx, dt
end

function L1MC(d::Vector{Tv}, rowIdx::Vector{Ti}, colIdx::Vector{Ti}, M::Ti, N::Ti,
              rk::Ti, mu::Tv, maxIter::Ti) where {Tv<:Float64, Ti<:Int64}
    # initialization
    w = randn(rk)
    U = randn(M, rk)
    V = randn(N, rk)
    for i = 1 : rk
        U[:,i] = U[:,i] / norm(U[:,i])
        V[:,i] = V[:,i] / norm(V[:,i])
    end
    X = zeros(M, N)
    Xr= zeros(M, N)
    Z = zeros(M, N)
    for i = 1 : length(d)
        X[rowIdx[i], colIdx[i]] = d[i]
    end
    # iterations
    for iter = 1 : maxIter
        Xr = copy(X)
        for i = 1 : rk
            if w[i] != 0
               # update column factors
               U[:,i] = (Xr * V[:,i]) / w[i]
               U[:,i] = U[:,i] / norm(U[:,i])
               # update row    factors
               V[:,i] = (Xr'* U[:,i]) / w[i]
               V[:,i] = V[:,i] / norm(V[:,i])
               # update weight by soft shrinking
               w[i]   = sum(Xr .* (U[:,i] * V[:,i]'))
               if w[i] > mu
                  w[i] = w[i] - mu
               elseif w[i] < -mu
                  w[i] = w[i] + mu
               else
                  w[i] = 0.0
               end
               Xr = Xr - (w[i]*U[:,i]*V[:,i]')
            end
        end
        # update X
        Z = X - Xr; res = 0.0
        X = copy(Z)
        for i = 1 : length(d)
            X[rowIdx[i], colIdx[i]] = d[i]
            res = res + (Z[rowIdx[i], colIdx[i]]-d[i])^2
        end
        res = sqrt(res) / norm(d)
        println("$iter relerr: $res")
    end
    # final trim of rank
    threshold = 1e-3 * length(d) / (M*N) * sum(abs.(w))
    count = 0;
    for i = 1 : rk
        if abs(w[i]) > threshold
           count = count + 1
        else
           w[i] = 0.0
        end
    end
    return count, w
end

# ==============================================================================
#                          Complex128 cases
# ==============================================================================
function L1MC(d::Vector{Tv}, rowIdx::Vector{Ti}, colIdx::Vector{Ti}, M::Ti, N::Ti,
              rk::Ti, mu::Float64, maxIter::Ti) where {Tv<:Complex128, Ti<:Int64}
    # initialization
    w = randn(  rk) + im*randn(  rk)
    U = randn(M,rk) + im*randn(M,rk)
    V = randn(N,rk) + im*randn(N,rk)
    # normalize factor vectors
    for i = 1 : rk
        U[:,i] = U[:,i] / norm(U[:,i])
        V[:,i] = V[:,i] / norm(V[:,i])
    end
    X = zeros(Complex128, M, N)
    Xr= zeros(Complex128, M, N)
    Z = zeros(Complex128, M, N)
    # insert observations
    for i = 1 : length(d)
        X[rowIdx[i], colIdx[i]] = d[i]
    end
    # iterations
    for iter = 1 : maxIter
        Xr = copy(X)
        for i = 1 : rk
            if norm(w[i]) > 0.0
               # update column factors
               U[:,i] = conj(w[i]) / (norm(w[i]))^2 * (Xr * V[:,i])
               U[:,i] = U[:,i] / norm(U[:,i])
               # update row    factors
               V[:,i] =      w[i]  / (norm(w[i]))^2 * (Xr'* U[:,i])
               V[:,i] = V[:,i] / norm(V[:,i])
               # update weight by soft shrinking
               w[i]   = sum(Xr .* conj(U[:,i] * V[:,i]'))
               if norm(w[i]) > mu
                  w[i] = (norm(w[i])-mu)/norm(w[i]) * w[i]
               else
                  w[i] = zero(Complex128)
               end
               Xr = Xr - (w[i]*U[:,i]*V[:,i]')
            end
        end
        # update X
        Z = X - Xr; res = 0.0
        X = copy(Z)
        for i = 1 : length(d)
            X[rowIdx[i], colIdx[i]] = d[i]
            res = res + (norm(Z[rowIdx[i], colIdx[i]]-d[i]))^2
        end
        res = sqrt(res) / norm(d)
        println("$iter, relative error: $res")
    end
    # final trim of rank
    threshold = 1e-3 * length(d) / (M*N) * sum(abs.(w))
    count = 0;
    for i = 1 : rk
        if abs(w[i]) > threshold
           count = count + 1
        else
           w[i] = zero(Complex128)
        end
    end
    return count, w
end

# ==============================================================================
#                         learn the properties of circulant matrix
# ==============================================================================
n = 256; c = rand(Complex128, n);
C = zeros(Complex128, n, n)
for k = 1 : n
    for i = k : n
        C[i,k] = c[i+1-k]
    end
    for i = 1 : k-1
        C[i,k] = c[n-k+i+1]
    end
end

# the multiplication of circulant matrix with a vector is
# C * v = ifft * (Fn * c) \cdot (Fn * v)
v = rand(Complex128, n)
b = C * v
b1= ifft(fft(c) .* fft(v))
(vecnorm(b1-b)) / (vecnorm(b))

# ==============================================================================
#                    test the properties of Toeplitze matrix
# ==============================================================================
L = 556; K = 392;
N = L + K - 1;

t = rand(N); T = zeros(eltype(t), L, K);
# form toeplitze matrix
for j = 1 : K
    istart = K - j + 1
    for i = 1 : L
        T[i,j] = t[istart+i-1]
    end
end

v = rand(eltype(t), K);
b = T * v;
c = zeros(eltype(t), N); c[1:L] = t[K:end]; c[L+1:end] = t[1:K-1];
v1 = vcat(v, zeros(eltype(t), L-1));
b1= ifft(fft(c) .* fft(v1))[1:L];
(vecnorm(b1-b)) / (vecnorm(b))


# ==============================================================================
#                    test the properties of Hankel matrix
# ==============================================================================
N = 1234;
L = floor(Int64, N/2) + 1; K = N - L + 1;
d = rand(Complex128, N)
H = zeros(eltype(d), L, K)
for j = 1 : K
    for i = 1 : L
        H[i,j] = d[i+j-1]
    end
end

 # compute H \times v
v = rand(eltype(d), K);
b = H * v
c = vcat(d[K:end], d[1:K-1]);
vhat = vcat(reverse(v), zeros(eltype(d), L-1))
b1 = ifft(fft(c) .* fft(vhat))[1:L]
(vecnorm(b1-b)) / (vecnorm(b))

# test compute the conjugate transpose of H times a vector
v = rand(L);
r = H' * v
c = conj(vcat(d[K:end], d[1:K-1])); # very important to use the conjugate of the vector
w = vcat(zeros(eltype(d),K-1), reverse(v))
r1 = ifft(fft(c) .* fft(w))[L:N]
(vecnorm(r1-r)) / (vecnorm(r))

# an equvialent algorithm to compute the conjugate transpose of H times a vector
c = conj(vcat(d[L:N], d[1:L-1]))
w = vcat(reverse(v), zeros(eltype(d), K-1))
r2 = ifft(fft(c) .* fft(w))[1:K]
(vecnorm(r2-r)) / (vecnorm(r))






#
