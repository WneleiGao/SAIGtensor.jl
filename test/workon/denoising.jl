using PyPlot, Seismic, SAIGtensor

# read seismic data
path_out = "Desktop/Channel3D/john_lake"
(d, h, ext) = SeisRead(path_out)

# separate 2D data into 3D cube
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
d  = reshape(d,n1,n2,n3);
d = d[:,1:203, 65:191];
# get a regular cubic from the data

# add bandlimited random noise to the data.
dn = SeisAddNoise(d, 1.0, L=5);
# ===================================================

# check the result
# ===================================================
pin = "/Users/wenlei/Desktop/randCP/rank20/"
dir = getPatchDir(pin, 50, 14, 9)
# pin = "/Users/wenlei/Desktop/randCP/rank15/"
# dir = getPatchDir(pin, 75, 20, 13)
dp  = UnPatch(dir)
# examine
dd = dp - dn;

SeisPlot(d[:,:,123], cmap="gray", pclip=99.5)
SeisPlot(dn[:,:,123], cmap="gray", pclip=99.5)
SeisPlot(dp[:,:,123], cmap="gray", pclip=99.5)
tmp = hcat(dn[:,:,123], dp[:,:,123], dd[:,:,123])
SeisPlot(tmp, cmap="gray", pclip=99.5)

SeisPlot(d[206,:,:], cmap="gray", pclip=99.5)
SeisPlot(dn[206,:,:], cmap="gray", pclip=99.5)
SeisPlot(dp[206,:,:], cmap="gray", pclip=99.5)
tmp = hcat(d[206,:,:], dn[206,:,:], dp[206,:,:], dd[206,:,:])
SeisPlot(tmp, cmap="gray", pclip=99.5)

# path to put data
pin = "/Users/wenlei/Desktop/randCP/rank15_20/"
# window size and window overlapping
it_wl = 15; it_wo=5;
x1_wl = 15; x1_wo=5;
x2_wl = 15; x2_wo=5;
# divided into patches
dir = patch(pin, dn, it_wl, it_wo, x1_wl, x1_wo, x2_wl, x2_wo)

# loop over each patch
par = Array{Tuple{String, Int64}}(length(dir));
R   = 30;
for i = 1 : length(dir)
    par[i] = (dir[i], R)
end
pmap(WrapRandCpAls, par)

dn = convert(Array{Float32,3}, dn);
dp = convert(Array{Float32,3}, dp);
dd = convert(Array{Float32,3}, dd);

path = "/Users/wenlei/Desktop/dn.bin"
fid = open(path, "w"); write(fid, vec(dn)); close(fid);
path = "/Users/wenlei/Desktop/dp.bin"
fid = open(path, "w"); write(fid, vec(dp)); close(fid);
path = "/Users/wenlei/Desktop/dd.bin"
fid = open(path, "w"); write(fid, vec(dd)); close(fid);

# taper each patch
par = Array{Tuple{String, Int64, Int64, Int64}}(length(dir));
for i = 1 : length(dir)
    par[i] = (dir[i], it_wo, x1_wo, x2_wo)
end
pmap(WrapTaper, par);
# unpatch
dp = UnPatch(dir)

fig = figure("1", figsize=(12,15));
it = 206; ix = 101; iy = 77; tmp=12; pclip=0.998; vdpi=300;
v=quantile(abs(vec(dn[it,:,:])), pclip);
subplot(3,3,1); SeisPlot(dn[it,:,:], vmax=v,vmin=-v, fignum="1", cmap="gray", ylabel="CMP y number", ticksize=tmp, labelsize=tmp, dpi=vdpi)
ax=gca(); ax[:set_yticks](collect(30:30:180)); ax[:set_xticks](collect(25:25:100));
ax[:text](-8., -3, "a)", fontsize=tmp, fontweight="bold")
subplot(3,3,2); SeisPlot(dp[it,:,:], vmax=v,vmin=-v, fignum="1", cmap="gray", xlabel="CMP x number", ticksize=tmp, labelsize=tmp, dpi=vdpi)
ax=gca(); ax[:set_yticks]([]); ax[:set_xticks](collect(25:25:100));
ax[:text](-8., -3, "b)", fontsize=tmp, fontweight="bold")
subplot(3,3,3); SeisPlot(dd[it,:,:], vmax=v,vmin=-v, fignum="1", cmap="gray", ticksize=tmp, labelsize=tmp, dpi=vdpi)
ax=gca(); ax[:set_yticks]([]); ax[:set_xticks](collect(25:25:100));
ax[:text](-8., -3, "c)", fontsize=tmp, fontweight="bold")

v=quantile(abs(vec(dn[:,ix,:])), pclip);
subplot(3,3,4); SeisPlot(dn[:,ix,:], vmax=v,vmin=-v, fignum="1", dy=0.004, cmap="gray", ylabel="Time (s)", ticksize=tmp, labelsize=tmp, dpi=vdpi)
ax=gca(); ax[:set_xticks](collect(25:25:100)); ax[:set_yticks]([0.5, 1.0, 1.5, 2.0, 2.5]);
ax[:text](-8., -0.04, "d)", fontsize=tmp, fontweight="bold")
subplot(3,3,5); SeisPlot(dp[:,ix,:], vmax=v,vmin=-v, fignum="1", dy=0.004, cmap="gray", xlabel="CMP x number", ticksize=tmp, labelsize=tmp, dpi=vdpi)
ax=gca(); ax[:set_xticks](collect(25:25:100)); ax[:set_yticks]([]);
ax[:text](-8., -0.04, "e)", fontsize=tmp, fontweight="bold")
subplot(3,3,6); SeisPlot(dd[:,ix,:], vmax=v,vmin=-v, fignum="1", dy=0.004, cmap="gray", ticksize=tmp, labelsize=tmp)
ax=gca(); ax[:set_xticks](collect(25:25:100)); ax[:set_yticks]([]);
ax[:text](-8., -0.04, "f)", fontsize=tmp, fontweight="bold")

v=quantile(abs(vec(dn[:,ix,:])), pclip);
subplot(3,3,7); SeisPlot(dn[:,:,iy], vmax=v,vmin=-v, fignum="1", dy=0.004, cmap="gray", ylabel="Time (s)", ticksize=tmp, labelsize=tmp, dpi=vdpi)
ax=gca(); ax[:set_xticks](collect(30:30:180)); ax[:set_yticks]([0.5, 1.0, 1.5, 2.0, 2.5]);
ax[:text](-12., -0.04, "g)", fontsize=tmp, fontweight="bold")
subplot(3,3,8); SeisPlot(dp[:,:,iy], vmax=v,vmin=-v, fignum="1", dy=0.004, cmap="gray", xlabel="CMP y number", ticksize=tmp, labelsize=tmp, dpi=vdpi)
ax=gca(); ax[:set_xticks](collect(30:30:180)); ax[:set_yticks]([]);
ax[:text](-12., -0.04, "h)", fontsize=tmp, fontweight="bold")
subplot(3,3,9); SeisPlot(dd[:,:,iy], vmax=v,vmin=-v, fignum="1", dy=0.004, cmap="gray", ticksize=tmp, labelsize=tmp, dpi=vdpi)
ax=gca();  ax[:set_xticks](collect(30:30:180)); ax[:set_yticks]([]);
ax[:text](-12., -0.04, "k)", fontsize=tmp, fontweight="bold")
tight_layout()
savefig("/Users/wenlei/Desktop/fig1.pdf", dpi=300);


# ================================ second real data example================================
# ================================ second real data example================================
using PyPlot, Seismic, SAIGtensor
# real data example on F3dome
path = "/Users/wenlei/Desktop/F3dome/preprocessed.bin"
fid = open(path, "r"); n = read(fid, Int64, 3); d = read(fid, Float64, prod(n)); close(fid);
d   = reshape(d, n[1], n[2], n[3])
dn  = SeisAddNoise(d, 2.0, L=5)
# path to put data


pin = "/Users/wenlei/Desktop/randCP/F3rank30-30/"
# window size and window overlapping
# it_wl = 90; it_wo=30;
# x1_wl = 30; x1_wo=5;
# x2_wl = 30; x2_wo=5;
# R   = 30;
# more complex data
it_wl = 18; it_wo=6;
x1_wl = 18; x1_wo=6;
x2_wl = 18; x2_wo=6;
R   = [7, 7, 7];
# divided into patches
dir = patch(pin, dn, it_wl, it_wo, x1_wl, x1_wo, x2_wl, x2_wo)

# loop over each patch by random CP decomposition
# par = Array{Tuple{String, Int64}}(length(dir));
# for i = 1 : length(dir)
#     par[i] = (dir[i], R)
# end
# pmap(WrapRandCpAls, par)

# Tucker ALS decomposition
par = Array{Tuple{String, Vector{Int64}, Float64}}(length(dir));
tol = 1e-4;
for i = 1 : length(dir)
    par[i] = (dir[i], R, tol)
end
pmap(WrapTuckerAls, par)

# taper each patch
par = Array{Tuple{String, Int64, Int64, Int64}}(length(dir));
for i = 1 : length(dir)
    par[i] = (dir[i], it_wo, x1_wo, x2_wo)
end
pmap(WrapTaper, par);
dp = UnPatch(dir)

# unpatch
dn  = SeisAddNoise(d, 1.5, L=5)
dd = dn-dp; ddori = d - dp;
iy = 1;
tmp = hcat(d[:,:,iy], dn[:,:,iy], dp[:,:,iy], dd[:,:,iy], ddori[:,:,iy]);
tmp = hcat(dn[:,:,iy], dp[:,:,iy], dd[:,:,iy]);
v = quantile(abs.(vec(tmp)), 0.99)
SeisPlot(tmp, cmap="gray", vmax=v*1.2, vmin=-v*0.8, wbox=20, hbox=5);
itl = 200; ixl=200; dn = SeisAddNoise(d[:,:,iy], 1.7, L=5)
tmp1 = hcat(dp[itl:itl+100, ixl:ixl+100, iy], d[itl:itl+100, ixl:ixl+100, iy], dn[itl:itl+100, ixl:ixl+100, iy], dn[itl:itl+100, ixl:ixl+100, iy]-dp[itl:itl+100,ixl:ixl+100,iy])
SeisPlot(tmp1, cmap="gray", pclip=99, hbox=4, wbox=16);

# write data =====
root = "/Users/wenlei/Desktop/randCP/result"
path = joinpath(root, "F3_18_6_7.bin")
fid  = open(path, "w"); write(fid, collect(size(dp)));
write(fid, vec(dp)); close(fid);

root = "/Users/wenlei/Desktop/randCP/result"
path = joinpath(root, "F3noise.bin")
fid  = open(path, "w"); dn = convert(Array{Float32,3}, dn)
write(fid, vec(dn)); close(fid);


# write into vtk data form.
# ======check the results===============================
using WriteVTK
root = homedir()
path = joinpath(root, "Desktop/seis")
(nt, n1, n2) = size(dp)
dt = 0.004; dx = 10; dy = 10;
taxis = collect(0.0:dt:(nt-1)*dt)
xaxis = collect(0.0:dx:(n1-1)*dx)
yaxis = collect(0.0:dy:(n2-1)*dy)

vtkfile = vtk_grid(path, taxis, xaxis, yaxis)
vtk_point_data(vtkfile, dp, "Amplitude")
outfile = vtk_save(vtkfile)
# write into vtk data form.
# ======check the results===============================
# read data =====
root = "/Users/wenlei/Desktop/randCP/result"
path = joinpath(root, "F3_18_6_7.bin")
path = "/Users/wenlei/Desktop/F3dome/preprocessed.bin"
fid  = open(path, "r"); n = read(fid, Int64, 3);
dp   = read(fid, Float64, prod(n)); close(fid);
dp   = reshape(dp, n...)
using PyPlot, Seismic
SeisPlot3D(dp, wbox=3, hbox=3)
