# ==============================================================================
#     compare the performance of CP decomposition against randomized version
# ==============================================================================
N = 3; I = [500,300,300]; R = 20;
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
raw = copy(X.D);
root= homedir(); path= join([root "/Desktop/tensor3D.bin"]);
fid = open(path, "w"); write(fid, vec(raw)); close(fid);

# add some random noise to pure data
noise = rand(I...); r = vecnorm(noise)/(vecnorm(X.D) *0.1)
noise = noise / r ; X.D = raw + noise;
# Compute relative fitness is more expensive than decomposition
his = randCpAls_fit(X, R); his = vcat(1,0, his);
t   = randCpAls_time(X, R); t = cumsum(t); t = vcat(0., t)
his1 = cpals_fit(X, R); his1 = vcat(1.0, his1);
t1 = cpals_time(X, R); t1 = cumsum(t1); t1 = vcat(0., t1)

# ===== 4D tensor ==============================================================
N = 4; I = [200,80,80,80]; R = 20;
ns = ceil(Int64, 10*R*log(R));
lambda = randn(R);
A = Array{Array{Float64,2}}(N);
Xst = Array{Array{Float64,2}}(N);
for m = 1 : N
    A[m] = rand(I[m],R)
    Xst[m] = zeros(ns,I[m])
end
P = cptensor(N, I, R, lambda, A);
X = cp2tensor(P); # it probably take more than 10 mins to convert a cptensor to dense tensor.
raw = reshape(raw, 200, 80, 80, 80);
noise = rand(I...);
r = vecnorm(noise)/(vecnorm(X.D) *0.1)
noise = noise / r;
X.D = X.D + noise;
his2 = randCpAls_fit(X, R); his2 = vcat(1.0, his2);
t2   = randCpAls_time(X, R);  t2 = cumsum(t2); t2 = vcat(0., t2)
his3 = cpals_fit(X, R); his3 = vcat(1.0, his3);
t3 = cpals_time(X, R) ; t3 = cumsum(t3); t3 = vcat(0., t3)

# plotting computation time VS relative fitness
fig = figure("1", figsize=(6,3))
tmp = 12
subplot(1,2,1);
plot(t[1:5:end], his[1:5:end], linewidth=2, marker=">", "k", markersize=8, label="RandALS")
plot(t1[1:5:end], his1[1:5:end], linewidth=2, marker="*", "k", markersize=8, label="ALS")
xlabel("Time (s)", fontsize=tmp); xticks(fontsize=tmp)
ylabel("relative error", fontsize=tmp); yticks(fontsize=tmp)
ax=gca(); ax[:set_xticks]([0, 4, 8, 12, 16])
ax[:set_yticks]([0., 0.2, 0.4, 0.6, 0.8])
ax[:text](-1.6, 1.01, "a)", fontsize=tmp, fontweight="bold")
legend(loc="upper right", fontsize=10)

subplot(1,2,2);
plot(t2[1:5:end], his2[1:5:end], linewidth=2, marker=">", "k", markersize=8, label="RandALS")
plot(t3[1:5:end], his3[1:5:end], linewidth=2, marker="*", "k", markersize=8, label="ALS")
ax=gca(); ax[:set_yticks]([]);
ax[:set_xticks]([0, 15, 30, 45, 60])
ax[:text](-6.0, 1.01, "b)", fontsize=tmp, fontweight="bold")
xlabel("Time (s)", fontsize=tmp); xticks(fontsize=tmp)
legend(loc="upper right", fontsize=10)
tight_layout()
