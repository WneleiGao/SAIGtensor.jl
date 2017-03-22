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

pin = "/Users/wenlei/Desktop/randCP/rank50/"
dir = getPatchDir(pin, 15, 4, 3)
d1 = UnPatch(dir);
it = 200
SeisPlot(hcat(d[it,:,:], d1[it,:,:], d1[it,:,:]-d[it,:,:]), cmap="gray",
         pclip=95, name=join([pin, "slice200.pdf"]))
ix = 53
SeisPlot(hcat(d[:,ix,:], d1[:,ix,:], d1[:,ix,:]-d[:,ix,:]), cmap="gray",
         pclip=95, name=join([pin, "x53.pdf"]))
iy = 56
SeisPlot(hcat(d[:,:,iy], d1[:,:,iy], d1[:,:,iy]-d[:,:,iy]), cmap="gray",
         pclip=95, name=join([pin, "y56.pdf"]))


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
SeisPlot(ds[:,:,27], cmap="gray", style="wiggles", wiggle_trace_increment=3)

# =====decompose as a whole patch============
N = 3;
I = collect(size(ds));
X = tensor(N, I, ds);
R = 100;
(lambda2, A2) = cpals(X, R, pflag=true);
X = cptensor(N, I, R, lambda2, A2);
X = cp2tensor(X)
dp = X.D;
SeisPlot(ds[:,:,27], cmap="gray", style="wiggles", wiggle_trace_increment=3)

tmp = hcat(ds[:,1:10:end,27],zeros(500,1), dp[:,1:10:end,25], zeros(500,1), dp[:,1:10:end,27]-ds[:,1:10:end,27]);
SeisPlot(tmp, style="wiggles", wbox=6, hbox=4, dy=0.004, ylabel="Time (s)", xlabel="Trace number")

# ========for curve events=================
(d, e) = SeisParabEvents(nx2=100);
snr = 1.;
ds = SeisAddNoise(d, snr, L=5);
N = 3;
I = collect(size(ds));
X = tensor(N, I, ds);
R = 100;
(lambda2, A2) = cpals(X, R, pflag=true);
X = cptensor(N, I, R, lambda2, A2);
X = cp2tensor(X)
dp = X.D;
tmp = hcat(ds[:,1:10:end,27],zeros(500,1), dp[:,1:10:end,25], zeros(500,1), dp[:,1:10:end,27]-ds[:,1:10:end,27]);
SeisPlot(tmp, style="wiggles", wbox=6, hbox=4, dy=0.004, ylabel="Time (s)", xlabel="Trace number")
ax = gca()
ax[:xaxis][:tick_top]()
ax[:xaxis][:set_label_position]("top")
tight_layout()
savefig("/Users/wenlei/Desktop/fig2.pdf")
