# form a new tensor
# N = 5; I = [23, 42, 13, 14, 17]; D = rand(I...)
# X = initTensor(N, I, D)
#
# # form a cptensor
# N = 5; I = [23, 42, 13, 14, 17];
# R = 10; lambda = rand(N);
# A = Array{Array{Float64,2}}(N)
# for m = 1 : N
#     A[m] = rand(I[m], R)
# end
# X = initCptensor(N, I, R, lambda, A);

# matricization and unmatricization
using SAIGtensor
N = 3; I = [111, 121, 321]; D = reshape(collect(1:prod(I)), I...)
X = initTensor(N, I, D)
tmp = 0.
for n = 1 : N
    Xn = matricization(X, n)
    X1 = unmatricization(Xn, n, I)
    tmp = tmp + vecnorm(X.D - X1.D)
end
@test tmp < 10e-6

# # test Khatri_Rao product
# # 1D case
# m1 = 23545; m2 = 345;
# a1 = rand(m1); a2 = rand(m2);
# @time r1 = kron(a2, a1);
# @time r = KRP(a1, a2);
# norm(r-r1)
# # 2D case
# m1 = 23545; m2 = 345; n = 3;
# a1 = rand(m1, n); a2 = rand(m2, n);
# function foo(a1::Array{Float64,2}, a2::Array{Float64,2})
#     m1 = size(a1,1);
#     m2 = size(a2,1);
#     m  = m1*m2
#     n = size(a1, 2)
#     r = Array{Float64}(m, n)
#     for j = 1 : n
#         r[:, j] = kron(a2[:,j], a1[:,j])
#     end
#     return r
# end
# @time r  = foo(a2, a1);
# @time r1 = KRP(a2, a1);
# vecnorm(r-r1)
#
# # test recursive Khatri-Rao product for matrix
# I = [21, 32, 43, 24]; R = 3;
# N = length(I);
# A = Array{Array{Float64,2}}(N)
# for m = 1 : N
#     A[m] = rand(I[m], R)
# end
# function foo(A, I, R)
#     r = zeros(prod(I), R)
#     for j = 1 : R
#         r[:,j] = kron(A[4][:,j], kron(A[3][:,j], kron(A[2][:,j], A[1][:,j])))
#     end
#     return r
# end
# @time r = foo(A, I, R);
# @time r1 = recursiveKRP(A);
# vecnorm(r-r1)
# # test recursive Khatri-Rao product for vector
# I = [21, 32, 43, 24];
# N = length(I);
# A = Array{Array{Float64,1}}(N)
# for m = 1 : N
#     A[m] = rand(I[m])
# end
# function foo(A)
#     r = kron(A[4], kron(A[3], kron(A[2], A[1])))
#     return r
# end
# @time r = foo(A);
# @time r1 = recursiveKRP(A);
# vecnorm(r-r1)
#
# # test on tensor times a set of vectors
# # the order can be exchanged
# I = [77,34,89,67];
# a1 = rand(I[1]); a2 = rand(I[2]);
# a3 = rand(I[3]); a4 = rand(I[4]);
# A = Array{Array{Float64,1}}(length(I)-1)
# A[1] = a1; A[2] = a3; A[3] = a4;
# X = rand(I...); X = tensor(length(I), I, X);
# X2 = matricization(X,2);
# t = recursiveKRP(A); r = X2 * t;
# # change the order of third and fourth dimension
# p = [1,2,4,3]; Inew = I[p];
# Xp = permutedims(X.D, (1,2,4,3)); Xp = tensor(length(Inew), Inew, Xp);
# Xp2 = matricization(Xp,2);
# B = Array{Array{Float64,1}}(length(I)-1)
# B[1] = a1; B[2] = a4; B[3] = a3;
# tp = recursiveKRP(B)
# rp = Xp2 * tp;
# norm(rp-r)
#
# # for better understanding the properties of permutedims
# I = [73,27,33,17]; N = length(I)
# X = rand(I...)
# od1 = (1,2,3,4);
# od2 = (3,4,1,2);
# X1 = permutedims(X, od1); X1m = reshape(X1, I[od1[1]]*I[od1[2]], I[od1[3]]*I[od1[4]]);
# X2 = permutedims(X, od2); X2m = reshape(X2, I[od2[1]]*I[od2[2]], I[od2[3]]*I[od2[4]]);
# vecnorm(X1m-X2m')
#
# od1 = (1,3,2,4);
# od2 = (2,4,1,3);
# X1 = permutedims(X, od1); X1m = reshape(X1, I[od1[1]]*I[od1[2]], I[od1[3]]*I[od1[4]]);
# X2 = permutedims(X, od2); X2m = reshape(X2, I[od2[1]]*I[od2[2]], I[od2[3]]*I[od2[4]]);
# vecnorm(X1m-X2m')
#
# # ===== bettering understanding of the properties of tensor times vector========
# # ===== test left side =====================
# I = [73,27,33,17]; N = length(I)
# X = rand(I...); X1 = copy(X);
# a = Array{Array{Float64,1}}(length(I))
# for m = 1 : length(I)
#     a[m] = rand(I[m])
# end
# n = 2
# for i = 1 : n
#     m = reshape(X, I[i], prod(I[i+1:N]))
#     X = m' * a[i]
# end
# m = reshape(X1, prod(I[1:n]), prod(I[n+1:N]))
# X2 = m'* kron(a[2], a[1])
# norm(X-X2)
#
# # ＝======== test right side =================
# I = [33,27,35,41,29]; N = length(I)
# X = rand(I...); X1 = copy(X);
# a = Array{Array{Float64,1}}(length(I))
# for m = 1 : length(I)
#     a[m] = rand(I[m])
# end
# n = 2
# t = kron(a[5], kron(a[4], a[3]))
# Xm = reshape(X, prod(I[1:n]), prod(I[n+1:end]))
# tmp = Xm * t;
#
# X1 = tensor(N, I, X1)
# for i = n+1 : N
#     Xi = matricization(X1, i)
#     tmp = Xi' * a[i]
#     I[i] = 1
#     X1  = unmatricization(tmp, i, I)
# end
# norm(vec(X1.D)-tmp)
#
# # # ========test change the order of tensor times vector========
# I = [33,27,35,41,29]; N = length(I)
# X = rand(I...); X1 = copy(X);
# a = Array{Array{Float64,1}}(length(I))
# for m = 1 : length(I)
#     a[m] = rand(I[m])
# end
# X = tensor(N, I, X)
# for i = 1 : N
#     Xi = matricization(X, i)
#     tmp = Xi' * a[i]
#     I[i] = 1
#     X  = unmatricization(tmp, i, I)
# end
# # ======== change the order========
# I = [33,27,35,41,29];
# X1 = tensor(N, I, X1)
# for i in [5,1,4]
#     Xi = matricization(X1, i)
#     tmp = Xi' * a[i]
#     I[i] = 1
#     X1 = unmatricization(tmp, i, I)
# end
# for i = [3,2]
#     Xi = matricization(X1, i)
#     tmp = Xi' * a[i]
#     I[i] = 1
#     X1 = unmatricization(tmp, i, I)
# end
# vec(X.D)
# vec(X1.D)
#
# # ===========test tensor times a set of tensor==================================
# I = [33,27,35,41,29]; N = length(I)
# X = tensor(N, I, rand(I...));
# dims = [5, 2, 4]
# v = Array{Array{Float64,1}}(length(dims))
# for m = 1 : length(dims)
#     v[m] = rand(I[dims[m]])
# end
# r = ttv(X, v, dims)
#
# for i = 1 : length(dims)
#     Xm = matricization(X, dims[i])
#     tmp = Xm' * v[i]
#     I[dims[i]] = 1
#     X = unmatricization(tmp, dims[i], I)
# end
# r1 = vec(X.D)
# norm(r1-r)
# # =========efficient way to compute the frobenius norm of cptensor=============
# N = 3; I = [23, 42, 13];
# R = 10; lambda = rand(R);
# A = Array{Array{Float64,2}}(N)
# for m = 1 : N
#     A[m] = rand(I[m], R)
# end
# X = cptensor(N, I, R, lambda, A)
# l = cpnorm(X)
# function cp2tensor(X::cptensor)
#     T = zeros(I...)
#     # for i5 = 1 : X.I[5]
#     #     for i4 = 1 : X.I[4]
#             for i3 = 1 : X.I[3]
#                 for i2 = 1 : X.I[2]
#                     for i1 = 1 : X.I[1]
#                         for r = 1 : R
#                             # T[i1,i2,i3,i4,i5] = T[i1,i2,i3,i4,i5] + lambda[r]*X.A[5][i5,r]*X.A[4][i4,r]*X.A[3][i3,r]*X.A[2][i2,r]*X.A[1][i1,r]
#                             T[i1,i2,i3] = T[i1,i2,i3] + lambda[r]*X.A[3][i3,r]*X.A[2][i2,r]*X.A[1][i1,r]
#                         end
#                     end
#                 end
#             end
#     #     end
#     # end
#     return T
# end
# T = cp2tensor(X)
# l1 = sqrt(sum(T.^2))
# l-l1 < 1e-6
#
# # ==effecient way to compute the fitness between full tensor and cptensor=======
# N = 3; I=[44, 11, 13]; R = 7;
# A = Array{Array{Float64,2}}(N); lambda = rand(R);
# for m = 1 : N
#     A[m] = rand(I[m],R)
# end
# X = cptensor(N, I, R, lambda, A)
# Y = rand(I...)
# function foo(Y::Array{Float64}, X::cptensor)
#     Y = vec(Y)
#     tmp = zeros(prod(I))
#     for r = 1 : X.R
#         tmp = tmp + kron(X.lambda[r]*X.A[3][:,r], kron(X.A[2][:,r], X.A[1][:,r]))
#     end
#     return sum(tmp .* Y)
# end
# Y1 = tensor(N, I, Y);
# @time tmp = foo(Y, X)
# @time ipt = innerprod(X, Y1);
# Xfull = cp2tensor(X)
# l = cpnorm(X); fit = sqrt(l^2 + sum(Y.^2) - 2*ipt);
# fit1 = sqrt(sum((Xfull.D-Y).^2))
# fit-fit1
#
#
# # ========test compute the right multiplication=================================
# N = 5; I = [21,32,43,44,45]; R = 7;
# ns = cptfirst(I);
# X = tensor(N, I, rand(I...));
# A = Array{Array{Float64,2}}(R);
# for m = 1 : N
#     A[m] = rand(I[m],R)
# end
# B = Array{Array{Float64}}(N-ns);
# for m = ns+1:N
#     B[m-ns] = A[m]
# end
# function test_speed(X::tensor, A::Array{Array{Float64,2}}, n::Int64)
#     Rnn = zeros(prod(I[1:n]), R) #keep the result
#     Jn = prod(I[1:N-1])
#     Kn = I[N]
#     T = reshape(X.D, Jn, Kn)
#     T = T * A[N]
#     for ir = 1 : R
#         tr = T[:, ir]
#         for m = N-1 : -1 : n+1
#             Jn = prod(I[1:m-1])
#             Kn = I[m]
#             tr = reshape(tr, Jn, Kn)
#             tr = tr * A[m][:,ir]
#         end
#         @inbounds for i = 1 : prod(I[1:n])
#             Rnn[i,ir] = tr[i]
#         end
#     end
#     return Rnn
# end
# function test_speed_left(X::tensor, A::Array{Array{Float64,2}}, n::Int64)
#     Jn = prod(X.I[1:n-1])
#     Kn = prod(X.I[n:end])
#     Xn = reshape(X.D, Jn, Kn)
#     tmp = recursiveKRP(A[1:n-1])
#     Lnn = Xn'*tmp
#     return Lnn
# end
#
# @time Rnn = cptRight(X, A, ns);
# @time Rnn1 = test_speed(X, A, ns);
# ns = 3
# @time Lnn = cptLeft(X, A, ns);
# @time Lnn1 = test_speed_left(X, A, ns);
#
# # =====cp_gradient=====================
# N = 5; I = [21,32,43,44,45]; R = 15;
# X = tensor(N, I, rand(I...));
# A = Array{Array{Float64,2}}(N);
# for m = 1 : N
#     A[m] = rand(I[m],R)
# end
# # test compute X × A[-n] from Right size
# tmp = 0.
# for n = 1 : N-1
#     G = ttf(X, A, n);
#     B = ttf_slow(X, A, n);
#     tmp = tmp + vecnorm(G-B)
# end
#
# tmp = 0.
# for ns = N-1:2
#     Rn = cptRight(X, A, ns);
#     m = ns-1
#     Rnsm1 = cptRight(Rn, A[m+1], X.I, m)
#     check = cptRight(X, A, m)
#     tmp = tmp + vecnorm(Rnsm1-check)
# end
#
# tmp = 0.
# for ns = 2:N-1
#     Ln = cptLeft(X, A, ns);
#     m = ns+1
#     Lnsp1 = cptLeft(Ln, A[m-1], X.I, m)
#     check = cptLeft(X, A, m)
#     tmp = tmp + vecnorm(Lnsp1-check)
# end
#
# # =====cp_gradient=====================
# N = 4; I = [73,41,39,29]; R = 20;
# lambda = randn(R)
# A = Array{Array{Float64,2}}(N);
# for m = 1 : N
#     A[m] = rand(I[m],R)
# end
# P = initCptensor(N, I, R, lambda, A);
# X = cp2tensor(P);
# (lambda, A1) = cpals(X, R);
# path = "/Users/wenyue/Desktop/d.bin"
# fid = open(path, "w")
# write(fid, convert(Array{Float64,1},X.I))
# write(fid, X.D)
# close(fid)
# (lambda, A) = cpals(X, R, pflag=true);
#
# # =========test sampling Z[n]===========
# N = 4; I = [42,35,42,23, 54]; R = 20;
# ns = ceil(Int64, 10*R*log(R))
# s = rand(collect(1:Int64(prod(I)/I[n])), ns);
# Su = sidx2sub(I, s);
# n = 1; In = [35, 42, 23, 54];
# S1 = ind2sub((In...), s);
# norm(S1[1]-Su[1][1])+norm(S1[2]-Su[1][2])+norm(S1[3]-Su[1][3])+norm(S1[4]-Su[1][4])
# n = 2; In = [42, 42, 23, 54];
# S1 = ind2sub((In...), s);
# norm(S1[1]-Su[2][1])+norm(S1[2]-Su[2][2])+norm(S1[3]-Su[2][3])+norm(S1[4]-Su[2][4])
# n = 3;  In = [42,35,23, 54];
# S1 = ind2sub((In...), s);
# norm(S1[1]-Su[3][1])+norm(S1[2]-Su[3][2])+norm(S1[3]-Su[3][3])+norm(S1[4]-Su[3][4])
# n = 4;  In = [42,35,42, 54];
# S1 = ind2sub((In...), s);
# norm(S1[1]-Su[4][1])+norm(S1[2]-Su[4][2])+norm(S1[3]-Su[4][3])+norm(S1[4]-Su[4][4])
# n = 5;  In = [42,35,42,23];
# S1 = ind2sub((In...), s);
# norm(S1[1]-Su[5][1])+norm(S1[2]-Su[5][2])+norm(S1[3]-Su[5][3])+norm(S1[4]-Su[5][4])
#
# # =======testing sampled KRP===================
# N = 4; I = [34,47,26,31]; R = 20;
# ns = ceil(Int64, 10*R*log(R))
# lambda = randn(R)
# A = Array{Array{Float64,2}}(N);
# Xst = Array{Array{Float64,2}}(N);
# for m = 1 : N
#     A[m] = rand(I[m],R)
#     Xst[m] = zeros(ns,I[m])
# end
# P = cptensor(N, I, R, lambda, A);
# T = cp2tensor(P);
# Zs = zeros(ns, R);
# n = 3
# @time idx = randSpl!(Zs, Xst, A, T, ns, n);
# It = (setdiff(I,[I[n]])...)
# rind = zeros(Int64, ns)
# for i = 1 : ns
#     rind[i] = sub2ind(It, idx[1][i], idx[2][i], idx[3][i])
# end
# dim = setdiff(collect(1:N), n);
# tmp = recursiveKRP(A[dim]);
# Z1 = tmp[rind, :];
# X = matricization(T, n);
# X1 = X[:,rind];
# vecnorm(Xst[n]-X1')
# vecnorm(Z1-Zs)
#
# # =======testing random cpals===================
N = 4; I = [61,62,63,64]; R = 20;
ns = ceil(Int64, 10*R*log(R));
lambda = randn(R);
A = Array{Array{Float64,2}}(N);
Xst = Array{Array{Float64,2}}(N);
for m = 1 : N
    A[m] = rand(I[m],R)
    Xst[m] = zeros(ns,I[m])
end
P = cptensor(N, I, R, lambda, A);
X = cp2tensor(P);
@time (lambda1, A1) = randCpAls(X, R, pflag=true);
@time (lambda2, A2) = cpals(X, R, pflag=true);
@time (lambda3, A3) = randCpAls_simplify(X, R);

normX = tnorm(X)
normresidual = sqrt(normX^2 + cpnorm(lambda,A)^2 - 2*innerprod(lambda, A, X))
fit = 1 - normresidual/normX


# =========measure efficiency========================
N = 3; I = [500,300,300]; R = 10;
ns = ceil(Int64, 10*R*log(R));
lambda = randn(R);
A = Array{Array{Float64,2}}(N);
Xst = Array{Array{Float64,2}}(N);
for m = 1 : N
    A[m] = rand(I[m],R)
    Xst[m] = zeros(ns,I[m])
end
P = cptensor(N, I, R, lambda, A);
X = cp2tensor(P);
noise = rand(I...);
r = vecnorm(noise)/(vecnorm(X.D *0.05)
noise = noise / r;
X.D = X.D + noise;
his = randCpAls_fit(X, R)
his = vcat(0.94, his);
t   = randCpAls_time(X, R)
t = cumsum(t)
t = vcat(0., t)
his1 = cpals_fit(X, R)
his1 = vcat(0.93, his1);
t1 = cpals_time(X, R)
t1 = cumsum(t1)
t1 = hcat(0., t1)
fig = figure(figsize=(6,6))
plot(t, his, linewidth=2, marker=".", "k", markersize=8)
plot(t1, his1, linewidth=2, marker="*", "k", markersize=8)
xlabel("Time (s)")
ylabel("relative error")




# ===========for 4D tensor=============
N = 4; I = [100,60,60,60,]; R = 10;
ns = ceil(Int64, 10*R*log(R));
lambda = randn(R);
A = Array{Array{Float64,2}}(N);
Xst = Array{Array{Float64,2}}(N);
for m = 1 : N
    A[m] = rand(I[m],R)
    Xst[m] = zeros(ns,I[m])
end
P = cptensor(N, I, R, lambda, A);
X = cp2tensor(P);
@time (lambda1, A1) = randCpAls(X, R, pflag=true);
@time (lambda2, A2) = cpals(X, R, pflag=true);
@time (lambda3, A3) = randCpAls_simplify(X, R);
