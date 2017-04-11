using PyPlot, Seismic, SAIGtensor

pin = "/Users/wenlei/Desktop/randCP/data_seis/john_lake"
(d,h,e) = SeisRead(pin);
n = length(h)
sx = zeros(Float32, n)
sy = zeros(Float32, n)
gx = zeros(Float32, n)
gy = zeros(Float32, n)
for i = 1 : n
    sx[i] = h[i].sx
    sy[i] = h[i].sy
    gx[i] = h[i].gx
    gy[i] = h[i].gy
end
tol = (maximum(sx) - minimum(sx)) / 2
p = []
tmp = 0
for i = 1 : n-1
    if abs(sx[i+1]-sx[i]) < tol
       tmp = tmp + 1
    else
       push!(p, tmp+1)
       tmp = 0
    end
    if i == n-1
       push!(p, tmp+1)
    end
end
n1 = size(d,1); n2 = p[1]; n3 = length(p);
d = reshape(d,n1,n2,n3);
d = d[:,1:203, 65:191];
# get a regular cubic from the data


# path to put data
pin = "/Users/wenlei/Desktop/randCP/rank50/"
# window size and window overlapping
it_wl = 60; it_wo=10;
x1_wl = 60; x1_wo=10;
x2_wl = 60; x2_wo=10;
# divided into patches
dir = patch(pin, d, it_wl, it_wo, x1_wl, x1_wo, x2_wl, x2_wo)
# loop over each patch
par = Array{Tuple{String, Int64}}(length(dir));
R   = 50;
for i = 1 : length(dir)
    par[i] = (dir[i], R)
end
pmap(WrapRandCpAls, par)

# taper each patch
par = Array{Tuple{String, Int64, Int64, Int64}}(length(dir));
for i = 1 : length(dir)
    par[i] = (dir[i], it_wo, x1_wo, x2_wo)
end
pmap(WrapTaper, par);
# unpatch


# ======check the results===============================
using PyPlot, Seismic, SAIGtensor
pin = "/Users/wenyue/Desktop/randCP/data_seis/john_lake"
(d,h,e) = SeisRead(pin);
n = length(h)
sx = zeros(Float32, n)
sy = zeros(Float32, n)
gx = zeros(Float32, n)
gy = zeros(Float32, n)
for i = 1 : n
    sx[i] = h[i].sx
    sy[i] = h[i].sy
    gx[i] = h[i].gx
    gy[i] = h[i].gy
end
tol = (maximum(sx) - minimum(sx)) / 2
p = []
tmp = 0
for i = 1 : n-1
    if abs(sx[i+1]-sx[i]) < tol
       tmp = tmp + 1
    else
       push!(p, tmp+1)
       tmp = 0
    end
    if i == n-1
       push!(p, tmp+1)
    end
end
n1 = size(d,1); n2 = p[1]; n3 = length(p);
d = reshape(d,n1,n2,n3);
d = d[:,1:203, 65:191];
ds = SeisAddNoise(d, 3, L=5);

# processing
pin = "/Users/wenyue/Desktop/randCP/rank60add/"
# window size and window overlapping
it_wl = 40; it_wo=10;
x1_wl = 40; x1_wo=10;
x2_wl = 40; x2_wo=10;
# divided into patches
dir = patch(pin, ds, it_wl, it_wo, x1_wl, x1_wo, x2_wl, x2_wo)
# loop over each patch
par = Array{Tuple{String, Int64}}(length(dir));
R   = 40;
for i = 1 : length(dir)
    par[i] = (dir[i], R)
end
pmap(WrapRandCpAls, par)

# taper each patch
par = Array{Tuple{String, Int64, Int64, Int64}}(length(dir));
for i = 1 : length(dir)
    par[i] = (dir[i], it_wo, x1_wo, x2_wo)
end
pmap(WrapTaper, par);
# end

pin = "/Users/wenyue/Desktop/randCP/rank40add/"
dir = getPatchDir(pin, 25, 7, 4)
d1 = UnPatch(dir);
dd = d1 - ds;

fig = figure("1", figsize=(8,4))
it = 205; v=quantile(vec(ds[it,:,:]), 0.90); tmp=15;
subplot(1,3,1); SeisPlot(ds[it,:,:], vmax=v,vmin=-v, fignum="1", cmap="gray", ylabel="CMP y number", ticksize=tmp)
ax=gca(); ax[:set_yticks]([50, 100, 150]); ax[:set_xticks]([50, 100]);
ax[:text](-14., -3, "a)", fontsize=tmp, fontweight="bold")
subplot(1,3,1); SeisPlot(d1[it,:,:], vmax=v,vmin=-v, fignum="1", cmap="gray", xlabel="CMP x number", ticksize=tmp)
ax=gca(); ax[:set_yticks]([]); ax[:set_xticks]([50, 100]);
ax[:text](-14., -3, "b)", fontsize=tmp, fontweight="bold")
subplot(1,3,3); SeisPlot(dd[it,:,:], vmax=v,vmin=-v, fignum="1", cmap="gray", ticksize=tmp)
ax=gca(); ax[:set_yticks]([]); ax[:set_xticks]([50, 100]);
ax[:text](-14, -3, "c)", fontsize=tmp, fontweight="bold")
tight_layout()
# savefig("/Users/wenyue/Desktop/fig5.pdf")

fig = figure("2", figsize=(8,8))
ix = 27; v=quantile(vec(ds[:,ix,:]), 0.95); tmp=14;
subplot(2,3,1); SeisPlot(ds[:,ix,:], vmax=v,vmin=-v, fignum="2", dy=0.004, cmap="gray", ylabel="Time (s)", ticksize=tmp)
ax=gca(); ax[:set_yticks]([0.5, 1.0, 1.5, 2.0, 2.5]); ax[:set_xticks]([]);
ax[:text](-15., -0.03, "a)", fontsize=tmp, fontweight="bold")

subplot(2,3,2); SeisPlot(d1[:,ix,:], vmax=v,vmin=-v, fignum="2", dy=0.004, cmap="gray", ylabel="Time (s)", ticksize=tmp)
ax=gca(); ax[:set_xticks]([]); ax[:set_yticks]([0.5, 1.0, 1.5, 2.0, 2.5]);
ax[:text](-15., -0.03, "b)", fontsize=tmp, fontweight="bold")

subplot(2,3,3); SeisPlot(dd[:,ix,:], vmax=v,vmin=-v, fignum="2", dy=0.004, cmap="gray", ylabel="Time (s)", xlabel="CMP x number", ticksize=tmp)
ax=gca(); ax[:set_xticks]([25, 50, 75]); ax[:set_yticks]([0.5, 1.0, 1.5, 2.0, 2.5]);
ax[:text](-15., -0.03, "c)", fontsize=tmp, fontweight="bold")

iy = 23; v=quantile(vec(ds[:,ix,:]), 0.95); tmp=14;

subplot(2,3,4); SeisPlot(ds[:,:,iy], vmax=v,vmin=-v, fignum="2", dy=0.004, cmap="gray", ticksize=tmp)
ax=gca(); ax[:set_yticks]([]); ax[:set_xticks]([]);
ax[:text](-20., -0.03, "d)", fontsize=tmp, fontweight="bold")

subplot(2,3,5); SeisPlot(d1[:,:,iy], vmax=v,vmin=-v, fignum="2", dy=0.004, cmap="gray", ticksize=tmp)
ax=gca(); ax[:set_yticks]([]); ax[:set_xticks]([]);
ax[:text](-20., -0.03, "e)", fontsize=tmp, fontweight="bold")

subplot(2,3,6); SeisPlot(dd[:,:,iy], vmax=v,vmin=-v, fignum="2", dy=0.004, cmap="gray",  xlabel="CMP y number", ticksize=tmp)
ax=gca(); ax[:set_yticks]([]); ax[:set_xticks]([50, 100, 150]);
ax[:text](-15., -0.03, "f)", fontsize=tmp, fontweight="bold")
tight_layout()
savefig("/Users/wenyue/Desktop/fig6.pdf")









# ============example with synthetic data==============
using PyPlot, Seismic, SAIGtensor
(d,e) = SeisLinearEvents(nt= 500, nx1=100, nx2=100, tau=[0.4, 0.8, 1.2, 1.6],
                         p1=[ 0.0003,  0.0004, -0.0004, -0.00027],
                         p2=[-0.0002, -0.0004,  0.0001,  0.0002 ],
                         p3=[0., 0., -0., -0.],
                         p4=[0., 0., -0., -0.],
                         amp=[1., 1.2, -0.5, -0.9])
SeisPlot(d[:,:,1], cmap="gray")
snr = 1.;
ds = SeisAddNoise(d, snr, L=5);
# =====decompose as a whole patch============
N = 3;
I = collect(size(ds));
X = tensor(N, I, ds);
R = 120;
(lambda2, A2) = cpals(X, R, pflag=true);
X = cptensor(N, I, R, lambda2, A2);
X = cp2tensor(X)
dp = X.D;
gap = zeros(250,1); r1=1:2:500; r2=1:10:100;
dd = ds - dp;

i3 = 1; d1 = hcat(ds[r1,r2,i3], gap, dp[r1,r2,i3], gap, dd[r1,r2,i3])
i3 = 25; d2 = hcat(ds[r1,r2,i3], gap, dp[r1,r2,i3], gap, dd[r1,r2,i3])
i3 = 50; d3 = hcat(ds[r1,r2,i3], gap, dp[r1,r2,i3], gap, dd[r1,r2,i3])


fig = figure("1", figsize=(6, 6)); tmp = 13;
subplot(3,1,1); SeisPlot(d1,fignum="1", style="wiggles", dy=0.004, dx=10)
ax=gca(); ax[:set_yticks]([0.3, 0.6, 0.9])
ax[:set_xticks]([])
ax[:text](-25., -0.05, "a)", fontsize=tmp, fontweight="bold")
subplot(3,1,2); SeisPlot(d2,fignum="1", style="wiggles", dy=0.004, dx=10, ylabel="Time (s)")
ax=gca(); ax[:set_yticks]([0.3, 0.6, 0.9])
ax[:set_xticks]([])
ax[:text](-25., -0.05, "b)", fontsize=tmp, fontweight="bold")
subplot(3,1,3); SeisPlot(d3,fignum="1", style="wiggles", dy=0.004, dx=10, xlabel="CMP x number")
ax=gca(); ax[:set_yticks]([0.3, 0.6, 0.9])
ax[:text](-25., -0.05, "c)", fontsize=tmp, fontweight="bold")







# ========for curve events=================
(d, e) = SeisParabEvents(nx2=100, nt=500, p2=[0.4, 0.4, 0.4], amp=[1.0, -1.0, 0.7]);
snr = 1.;
ds = SeisAddNoise(d, snr, L=5);
N = 3;
I = collect(size(ds));
X = tensor(N, I, ds);
R = 120;
(lambda2, A2) = cpals(X, R, pflag=true);
X = cptensor(N, I, R, lambda2, A2);
X = cp2tensor(X)

dp = X.D;
gap = zeros(250,1); r1=1:2:500; r2=10:10:90;
dd = ds - dp;

i3 = 1; d1 = hcat(ds[r1,r2,i3], gap, dp[r1,r2,i3], gap, dd[r1,r2,i3])
i3 = 40; d2 = hcat(ds[r1,r2,i3], gap, dp[r1,r2,i3], gap, dd[r1,r2,i3])
i3 = 80; d3 = hcat(ds[r1,r2,i3], gap, dp[r1,r2,i3], gap, dd[r1,r2,i3])

fig = figure("1", figsize=(6, 6)); tmp = 13;
subplot(3,1,1); SeisPlot(d1,fignum="1", style="wiggles", dy=0.004, dx=10)
ax=gca(); ax[:set_yticks]([0.3, 0.6, 0.9])
ax[:set_xticks]([])
ax[:text](-25., -0.05, "a)", fontsize=tmp, fontweight="bold")
subplot(3,1,2); SeisPlot(d2,fignum="1", style="wiggles", dy=0.004, dx=10, ylabel="Time (s)")
ax=gca(); ax[:set_yticks]([0.3, 0.6, 0.9])
ax[:set_xticks]([])
ax[:text](-25., -0.05, "b)", fontsize=tmp, fontweight="bold")
subplot(3,1,3); SeisPlot(d3,fignum="1", style="wiggles", dy=0.004, dx=10, xlabel="CMP x number")
ax=gca(); ax[:set_yticks]([0.3, 0.6, 0.9])
ax[:text](-25., -0.05, "c)", fontsize=tmp, fontweight="bold")
