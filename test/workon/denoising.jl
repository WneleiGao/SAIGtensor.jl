using PyPlot, Seismic, SAIGtensor

# ==============================================================================
#                          load data
# ==============================================================================
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
# get a regular cubic from the data
d = d[1:600,21:200, 66:185];
# add bandlimited random noise to the data.
snr=2.0; dn = SeisAddNoise(d, snr, L=5);

tmp="/Users/wenlei/Desktop/clean.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(d))); close(tmpfid);
tmp="/Users/wenlei/Desktop/input.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dn))); close(tmpfid);

pscube < /Users/wenlei/Desktop/clean.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=600 d1=0.004 d1num=0.4 f1num=0.4 label1="Time (s)" n2=180 d2num=40 f2num=40 label2="Crossline" n3=120 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 bclip=-14387.4 wclip=13873.6 > /Users/wenlei/Desktop/clean.eps
pscube < /Users/wenlei/Desktop/input.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=600 d1=0.004 d1num=0.4 f1num=0.4 label1="Time (s)" n2=180 d2num=40 f2num=40 label2="Crossline" n3=120 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 bclip=-14387.4 wclip=13873.6 > /Users/wenlei/Desktop/input.eps

# ==============================================================================
#                           CPALS
# ==============================================================================
pin = "/Users/wenlei/Desktop/randCP/"
it_wl = 30; it_wo=15;
x1_wl = 30; x1_wo=15;
x2_wl = 30; x2_wo=15;
R   = [10, 10, 10];
# divided into patches
dir = patch(pin, dn, it_wl, it_wo, x1_wl, x1_wo, x2_wl, x2_wo)

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
dd = dn - dp;

tmp="/Users/wenlei/Desktop/output.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dp))); close(tmpfid);
tmp="/Users/wenlei/Desktop/diff.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dd))); close(tmpfid);
it=207;
tmp="/Users/wenlei/Desktop/tscl.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec( d[it,:,:]))); close(tmpfid);
tmp="/Users/wenlei/Desktop/tspp.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dp[it,:,:]))); close(tmpfid);
tmp="/Users/wenlei/Desktop/tsin.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dn[it,:,:]))); close(tmpfid);
tmp="/Users/wenlei/Desktop/tsdd.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dd[it,:,:]))); close(tmpfid);


pscube < /Users/wenlei/Desktop/output.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=600 d1=0.004 d1num=0.4 f1num=0.4 label1="Time (s)" n2=180 d2num=40 f2num=40 label2="Crossline" n3=120 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 bclip=-14387.4 wclip=13873.6 > /Users/wenlei/Desktop/output.eps
pscube < /Users/wenlei/Desktop/diff.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=600 d1=0.004 d1num=0.4 f1num=0.4 label1="Time (s)" n2=180 d2num=40 f2num=40 label2="Crossline" n3=120 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 bclip=-14387.4 wclip=13873.6 > /Users/wenlei/Desktop/diff.eps
psimage < /Users/wenlei/Desktop/tscl.bin height=3.0 width=2.0 labelsize=12 n1=180 d1=1 d1num=40 f1num=40 label1="Crossline" n2=120 d2num=20 f2num=20 label2="Inline" xbox=0.5 ybox=0.5 bclip=-14387.4 wclip=13873.6 > /Users/wenlei/Desktop/tscl.eps
psimage < /Users/wenlei/Desktop/tspp.bin height=3.0 width=2.0 labelsize=12 n1=180 d1=1 d1num=40 f1num=40 label1="Crossline" n2=120 d2num=20 f2num=20 label2="Inline" xbox=0.5 ybox=0.5 bclip=-14387.4 wclip=13873.6 > /Users/wenlei/Desktop/tsdp.eps
psimage < /Users/wenlei/Desktop/tsin.bin height=3.0 width=2.0 labelsize=12 n1=180 d1=1 d1num=40 f1num=40 label1="Crossline" n2=120 d2num=20 f2num=20 label2="Inline" xbox=0.5 ybox=0.5 bclip=-14387.4 wclip=13873.6 > /Users/wenlei/Desktop/tsin.eps
psimage < /Users/wenlei/Desktop/tsdd.bin height=3.0 width=2.0 labelsize=12 n1=180 d1=1 d1num=40 f1num=40 label1="Crossline" n2=120 d2num=20 f2num=20 label2="Inline" xbox=0.5 ybox=0.5 bclip=-14387.4 wclip=13873.6 > /Users/wenlei/Desktop/tsdd.eps

# ==============================================================================
#                           MSSA
# ==============================================================================
pin1 = "/Users/wenlei/Desktop/MSSA/"
it_wl = 60; it_wo=30;
x1_wl = 30; x1_wo=15;
x2_wl = 30; x2_wo=15;
dir1 = patch(pin1, dn, it_wl, it_wo, x1_wl, x1_wo, x2_wl, x2_wo)

# decomposition
rk   = 10; dt=0.004f0; flow=1.0f0; fhigh=60.0f0;
par1 = Array{Tuple{String, Float32, Int64, Float32, Float32}}(length(dir1));
for i = 1 : length(dir1)
    par1[i] = (dir1[i], dt, rk, flow, fhigh)
end
pmap(WrapDenoiseOnePatchMSSA3D, par1)

# tapering and unpatching
par1 = Array{Tuple{String, Int64, Int64, Int64}}(length(dir1));
for i = 1 : length(dir1)
    par1[i] = (dir1[i], it_wo, x1_wo, x2_wo)
end
pmap(WrapTaper, par1);
dp1 = UnPatch(dir1)
dd1 = dn - dp1;

tmp="/Users/wenlei/Desktop/output1.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dp1))); close(tmpfid);
tmp="/Users/wenlei/Desktop/diff1.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dd1))); close(tmpfid);
it=207;
tmp="/Users/wenlei/Desktop/tspp1.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dp1[it,:,:]))); close(tmpfid);
tmp="/Users/wenlei/Desktop/tsdd1.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dd1[it,:,:]))); close(tmpfid);

pscube < /Users/wenlei/Desktop/output1.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=600 d1=0.004 d1num=0.4 f1num=0.4 label1="Time (s)" n2=180 d2num=40 f2num=40 label2="Inline" n3=120 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 bclip=-14387.4 wclip=13873.6 > /Users/wenlei/Desktop/output1.eps
pscube < /Users/wenlei/Desktop/diff1.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=600 d1=0.004 d1num=0.4 f1num=0.4 label1="Time (s)" n2=180 d2num=40 f2num=40 label2="Inline" n3=120 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 bclip=-14387.4 wclip=13873.6 > /Users/wenlei/Desktop/diff1.eps
psimage < /Users/wenlei/Desktop/tspp1.bin height=3.0 width=2.0 labelsize=12 n1=180 d1=1 d1num=40 f1num=40 label1="Crossline" n2=120 d2num=20 f2num=20 label2="Inline" xbox=0.5 ybox=0.5 bclip=-14387.4 wclip=13873.6 > /Users/wenlei/Desktop/tsdp1.eps
psimage < /Users/wenlei/Desktop/tsdd1.bin height=3.0 width=2.0 labelsize=12 n1=180 d1=1 d1num=40 f1num=40 label1="Crossline" n2=120 d2num=20 f2num=20 label2="Inline" xbox=0.5 ybox=0.5 bclip=-14387.4 wclip=13873.6 > /Users/wenlei/Desktop/tsdd1.eps

# ==============================================================================
#                         Julia plotting
# ==============================================================================
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

# ================================ second real data example================================
# ================================ second real data example================================
using PyPlot, Seismic, SAIGtensor
# real data example on F3dome
path = "/Users/wenlei/Desktop/F3dome/preprocessed.bin"
fid = open(path, "r"); n = read(fid, Int64, 3); d = read(fid, Float64, prod(n)); close(fid);
d   = reshape(d, n[1], n[2], n[3])
dn  = SeisAddNoise(d, 2.0, L=5)
# path to put data

pin = "/Users/wenlei/Desktop/randCP/"
# more complex data
it_wl = 18; it_wo=6;
x1_wl = 18; x1_wo=6;
x2_wl = 18; x2_wo=6;
R   = [7, 7, 7];
# divided into patches
dir = patch(pin, dn, it_wl, it_wo, x1_wl, x1_wo, x2_wl, x2_wo)

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
dd = dn-dp;

tmp="/Users/wenlei/Desktop/clean.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(d))); close(tmpfid);
tmp="/Users/wenlei/Desktop/input.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dn))); close(tmpfid);
tmp="/Users/wenlei/Desktop/output.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dp))); close(tmpfid);
tmp="/Users/wenlei/Desktop/diff.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dd))); close(tmpfid);
it=251;
tmp="/Users/wenlei/Desktop/tscl.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec( d[it,:,:]'))); close(tmpfid);
tmp="/Users/wenlei/Desktop/tspp.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dp[it,:,:]'))); close(tmpfid);
tmp="/Users/wenlei/Desktop/tsin.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dn[it,:,:]'))); close(tmpfid);
tmp="/Users/wenlei/Desktop/tsdd.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dd[it,:,:]'))); close(tmpfid);
ixl=221; ixu=370; iyl=201; iyu=300;
tmp="/Users/wenlei/Desktop/tspproom.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dp[it,ixl:ixu,iyl:iyu]'))); close(tmpfid);
tmp="/Users/wenlei/Desktop/tsinroom.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dn[it,ixl:ixu,iyl:iyu]'))); close(tmpfid);
tmp="/Users/wenlei/Desktop/tsddroom.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dd[it,ixl:ixu,iyl:iyu]'))); close(tmpfid);

iy=100; itl=251; itu=350; ixl=300; ixu=500;
tmp="/Users/wenlei/Desktop/ssroom.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec( d[itl:itu,ixl:ixu,iy]))); close(tmpfid);
tmp="/Users/wenlei/Desktop/inroom.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dn[itl:itu,ixl:ixu,iy]))); close(tmpfid);
tmp="/Users/wenlei/Desktop/pproom.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dp[itl:itu,ixl:ixu,iy]))); close(tmpfid);
tmp="/Users/wenlei/Desktop/ddroom.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dd[itl:itu,ixl:ixu,iy]))); close(tmpfid);



pscube < /Users/wenlei/Desktop/clean.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=400 d1=0.004 d1num=0.4 f1num=0.4 label1="Time (s)" n2=600 d2num=80 f2num=80 label2="Crossline" n3=400 d3num=80 f3num=80 label3="Inline" xbox=0.5 ybox=0.5 bclip=-2.36042 wclip=2.40392 > /Users/wenlei/Desktop/clean.eps
pscube < /Users/wenlei/Desktop/input.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=400 d1=0.004 d1num=0.4 f1num=0.4 label1="Time (s)" n2=600 d2num=80 f2num=80 label2="Crossline" n3=400 d3num=80 f3num=80 label3="Inline" xbox=0.5 ybox=0.5 bclip=-2.36042 wclip=2.40392 > /Users/wenlei/Desktop/input.eps
pscube < /Users/wenlei/Desktop/output.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=400 d1=0.004 d1num=0.4 f1num=0.4 label1="Time (s)" n2=600 d2num=80 f2num=80 label2="Crossline" n3=400 d3num=80 f3num=80 label3="Inline" xbox=0.5 ybox=0.5 bclip=-2.36042 wclip=2.40392 > /Users/wenlei/Desktop/output.eps
pscube < /Users/wenlei/Desktop/diff.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=400 d1=0.004 d1num=0.4 f1num=0.4 label1="Time (s)" n2=600 d2num=80 f2num=80 label2="Crossline" n3=400 d3num=80 f3num=80 label3="Inline" xbox=0.5 ybox=0.5 bclip=-2.36042 wclip=2.40392 > /Users/wenlei/Desktop/diff.eps
psimage < /Users/wenlei/Desktop/tscl.bin height=2.0 width=3.0 labelsize=12 n1=400 d1=1 d1num=100 f1num=100 label1="Crossline" n2=600 d2num=100 f2num=100 label2="Inline" xbox=0.5 ybox=0.5 bclip=-2.36042 wclip=2.40392 > /Users/wenlei/Desktop/tscl.eps
psimage < /Users/wenlei/Desktop/tspp.bin height=2.0 width=3.0 labelsize=12 n1=400 d1=1 d1num=100 f1num=100 label1="Crossline" n2=600 d2num=100 f2num=100 label2="Inline" xbox=0.5 ybox=0.5 bclip=-2.36042 wclip=2.40392 > /Users/wenlei/Desktop/tsdp.eps
psimage < /Users/wenlei/Desktop/tsin.bin height=2.0 width=3.0 labelsize=12 n1=400 d1=1 d1num=100 f1num=100 label1="Crossline" n2=600 d2num=100 f2num=100 label2="Inline" xbox=0.5 ybox=0.5 bclip=-2.36042 wclip=2.40392 > /Users/wenlei/Desktop/tsin.eps
psimage < /Users/wenlei/Desktop/tsdd.bin height=2.0 width=3.0 labelsize=12 n1=400 d1=1 d1num=100 f1num=100 label1="Crossline" n2=600 d2num=100 f2num=100 label2="Inline" xbox=0.5 ybox=0.5 bclip=-2.36042 wclip=2.40392 > /Users/wenlei/Desktop/tsdd.eps

psimage < /Users/wenlei/Desktop/tspproom.bin height=2.0 width=3.0 labelsize=12 n1=100 d1=1 f1=201 d1num=30 f1num=231 label1="Crossline" n2=150 f2=221 d2num=30 f2num=250 label2="Inline" xbox=0.5 ybox=0.5 bclip=-2.36042 wclip=2.40392 > /Users/wenlei/Desktop/tsdproom.eps
psimage < /Users/wenlei/Desktop/tsinroom.bin height=2.0 width=3.0 labelsize=12 n1=100 d1=1 f1=201 d1num=30 f1num=231 label1="Crossline" n2=150 f2=221 d2num=30 f2num=250 label2="Inline" xbox=0.5 ybox=0.5 bclip=-2.36042 wclip=2.40392 > /Users/wenlei/Desktop/tsinroom.eps
psimage < /Users/wenlei/Desktop/tsddroom.bin height=2.0 width=3.0 labelsize=12 n1=100 d1=1 f1=201 d1num=30 f1num=231 label1="Crossline" n2=150 f2=221 d2num=30 f2num=250 label2="Inline" xbox=0.5 ybox=0.5 bclip=-2.36042 wclip=2.40392 > /Users/wenlei/Desktop/tsddroom.eps

psimage < /Users/wenlei/Desktop/ssroom.bin height=3.0 width=2.0 labelsize=12 n1=100 f1=1.0 d1=0.004 d1num=0.1 f1num=1.1 label1="Time (s)" n2=200 f2=400 d2num=40 f2num=340 label2="Crossline" xbox=0.5 ybox=0.5 bclip=-2.36042 wclip=2.40392 > /Users/wenlei/Desktop/ssroom.eps
psimage < /Users/wenlei/Desktop/inroom.bin height=3.0 width=2.0 labelsize=12 n1=100 f1=1.0 d1=0.004 d1num=0.1 f1num=0.1 label1="Time (s)" n2=200 f2=400 d2num=40 f2num=40 label2="Crossline" xbox=0.5 ybox=0.5 bclip=-2.36042 wclip=2.40392 > /Users/wenlei/Desktop/inroom.eps
psimage < /Users/wenlei/Desktop/pproom.bin height=3.0 width=2.0 labelsize=12 n1=100 f1=1.0 d1=0.004 d1num=0.1 f1num=0.1 label1="Time (s)" n2=200 f2=400 d2num=40 f2num=40 label2="Crossline" xbox=0.5 ybox=0.5 bclip=-2.36042 wclip=2.40392 > /Users/wenlei/Desktop/pproom.eps
psimage < /Users/wenlei/Desktop/ddroom.bin height=3.0 width=2.0 labelsize=12 n1=100 f1=1.0 d1=0.004 d1num=0.1 f1num=0.1 label1="Time (s)" n2=200 f2=400 d2num=40 f2num=40 label2="Crossline" xbox=0.5 ybox=0.5 bclip=-2.36042 wclip=2.40392 > /Users/wenlei/Desktop/ddroom.eps



















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

# window size and window overlapping
# it_wl = 90; it_wo=30;
# x1_wl = 30; x1_wo=5;
# x2_wl = 30; x2_wo=5;
# R   = 30;
# loop over each patch by random CP decomposition
# par = Array{Tuple{String, Int64}}(length(dir));
# for i = 1 : length(dir)
#     par[i] = (dir[i], R)
# end
# pmap(WrapRandCpAls, par)
