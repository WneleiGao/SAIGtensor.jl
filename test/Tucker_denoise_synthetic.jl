# ==============================================================================
#                            test MSSA on Linear events
# ==============================================================================
using PyPlot, Seismic, SAIGtensor
(d,ext)=SeisLinearEvents(nt= 250, dt=0.008, nx1=200, nx2=200,
                         tau=[0.4, 0.8, 1.2, 1.6],
                         p1=[ 0.0003,  0.0004, -0.0004, -0.00027],
                         p2=[-0.0002, -0.0004,  0.0001,  0.0002 ],
                         p3=[0., 0., -0., -0.],
                         p4=[0., 0., -0., -0.],
                         amp=[1., 1.2, -0.5, -0.9])
snr = 1.0;
dn = SeisAddNoise(d, snr, L=5);

# patching
root = homedir(); pin = join([root "/Desktop/MSSA/"])
it_wl = 60; it_wo=20;
x1_wl = 30; x1_wo=15;
x2_wl = 30; x2_wo=15;
dir = patch(pin, dn, it_wl, it_wo, x1_wl, x1_wo, x2_wl, x2_wo)
# decomposition
rk   = 4; dt=0.008f0; flow=1.0f0; fhigh=40.0f0;
par = Array{Tuple{String, Float32, Int64, Float32, Float32}}(length(dir));
for i = 1 : length(dir)
    par[i] = (dir[i], dt, rk, flow, fhigh)
end
pmap(WrapDenoiseOnePatchMSSA3D, par)

# tapering and unpatching
par = Array{Tuple{String, Int64, Int64, Int64}}(length(dir));
for i = 1 : length(dir)
    par[i] = (dir[i], it_wo, x1_wo, x2_wo)
end
pmap(WrapTaper, par);
dp = UnPatch(dir)
dd = dn - dp;

tmp="/Users/wenlei/Desktop/input.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dn))); close(tmpfid);
tmp="/Users/wenlei/Desktop/output.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dp))); close(tmpfid);
tmp="/Users/wenlei/Desktop/diff.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dn-dp))); close(tmpfid);

# pscube < /Users/wenlei/Desktop/input.bin size1=4 size2=3 size3=1.5 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 perc=0.96 > /Users/wenlei/Desktop/input.eps
pscube < /Users/wenlei/Desktop/input.bin size1=4 size2=3 size3=1.5 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 bclip=0.96 wclip=-0.96 > /Users/wenlei/Desktop/input.eps
pscube < /Users/wenlei/Desktop/output.bin size1=4 size2=3 size3=1.5 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 bclip=0.96 wclip=-0.96 > /Users/wenlei/Desktop/output.eps
pscube < /Users/wenlei/Desktop/diff.bin size1=4 size2=3 size3=1.5 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 bclip=0.96 wclip=-0.96 > /Users/wenlei/Desktop/diff.eps

gap = zeros(size(dn,1),12)
# plotting
iy = 60
tmp = hcat(dn[:,:,iy], gap, dp[:,:,iy], gap, dd[:,:,iy])
SeisPlot(tmp, hbox=5, wbox=9, style="wiggles", xcur=1.6, wiggle_trace_increment=6, ticksize=15, ylabel="Time (s)", xlabel="CDP X", dy=0.008, dx=1)
ax=gca(); ax[:set_yticks](collect(0.2:0.4:2.0))

# ==============================================================================
#                           CP on linear events
# ==============================================================================
using PyPlot, Seismic, SAIGtensor
(d,ext)=SeisLinearEvents(nt= 250, dt=0.008, nx1=200, nx2=200,
                         tau=[0.4, 0.8, 1.2, 1.6],
                         p1=[ 0.0003,  0.0004, -0.0004, -0.00027],
                         p2=[-0.0002, -0.0004,  0.0001,  0.0002 ],
                         p3=[0., 0., -0., -0.],
                         p4=[0., 0., -0., -0.],
                         amp=[1., 1.2, -0.5, -0.9])
snr = 1.0;
dn = SeisAddNoise(d, snr, L=5);

# patching
root = homedir(); pin = join([root "/Desktop/randCP/"])
it_wl = 30; it_wo=15;
x1_wl = 30; x1_wo=15;
x2_wl = 30; x2_wo=15;
dir = patch(pin, dn, it_wl, it_wo, x1_wl, x1_wo, x2_wl, x2_wo)

# decomposition
R   = [10,10,10]; tol = 1e-4;
par = Array{Tuple{String, Vector{Int64}, Float64}}(length(dir));
for i = 1 : length(dir)
    par[i] = (dir[i], R, tol)
end
pmap(WrapTuckerAls, par)

# tapering and unpatching
par = Array{Tuple{String, Int64, Int64, Int64}}(length(dir));
for i = 1 : length(dir)
    par[i] = (dir[i], it_wo, x1_wo, x2_wo)
end
pmap(WrapTaper, par);
dp = UnPatch(dir)
dd = dn - dp;

tmp="/Users/wenlei/Desktop/input.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dn))); close(tmpfid);
tmp="/Users/wenlei/Desktop/output.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dp))); close(tmpfid);
tmp="/Users/wenlei/Desktop/diff.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dn-dp))); close(tmpfid);

# pscube < /Users/wenlei/Desktop/input.bin size1=4 size2=3 size3=1.5 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 perc=0.96 > /Users/wenlei/Desktop/input.eps
pscube < /Users/wenlei/Desktop/input.bin size1=4 size2=3 size3=1.5 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 bclip=0.96 wclip=-0.96 > /Users/wenlei/Desktop/input.eps
pscube < /Users/wenlei/Desktop/output.bin size1=4 size2=3 size3=1.5 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 bclip=0.96 wclip=-0.96 > /Users/wenlei/Desktop/output.eps
pscube < /Users/wenlei/Desktop/diff.bin size1=4 size2=3 size3=1.5 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 bclip=0.96 wclip=-0.96 > /Users/wenlei/Desktop/diff.eps

gap = zeros(size(dn,1),12)
# plotting
iy = 60
tmp = hcat(dn[:,:,iy], gap, dp[:,:,iy], gap, dd[:,:,iy])
SeisPlot(tmp, hbox=5, wbox=9, style="wiggles", xcur=1.6, wiggle_trace_increment=6, ticksize=15, ylabel="Time (s)", xlabel="CDP X", dy=0.008, dx=1)
ax=gca(); ax[:set_yticks](collect(0.2:0.4:2.0))

iy = 160
tmp = hcat(dn[:,:,iy], gap, dp[:,:,iy], gap, dd[:,:,iy])
SeisPlot(tmp, hbox=5, wbox=9, style="wiggles", xcur=1.6, wiggle_trace_increment=6, ticksize=15, ylabel="Time (s)", xlabel="CDP X", dy=0.008, dx=1)
ax=gca(); ax[:set_yticks](collect(0.2:0.4:2.0))

# ==============================================================================
#                            test MSSA on hyperbolic events
# ==============================================================================
using PyPlot, Seismic, SAIGtensor

(d, ext) = SeisParabEvents(nx1=200, nx2=200, dx1=5, dx2=5, nt=250, dt=0.008, p2=[0.4, -0.4, 0.6], amp=[1.0, -1.0, 0.7]);
snr = 1.0; dn = SeisAddNoise(d, snr, L=5);

# patching
root = homedir(); pin = join([root "/Desktop/MSSA/"])
it_wl = 60; it_wo=20;
x1_wl = 30; x1_wo=15;
x2_wl = 30; x2_wo=15;
dir = patch(pin, dn, it_wl, it_wo, x1_wl, x1_wo, x2_wl, x2_wo)
# decomposition
rk   = 3; dt=0.008f0; flow=1.0f0; fhigh=40.0f0;
par = Array{Tuple{String, Float32, Int64, Float32, Float32}}(length(dir));
for i = 1 : length(dir)
    par[i] = (dir[i], dt, rk, flow, fhigh)
end
pmap(WrapDenoiseOnePatchMSSA3D, par)

# tapering and unpatching
par = Array{Tuple{String, Int64, Int64, Int64}}(length(dir));
for i = 1 : length(dir)
    par[i] = (dir[i], it_wo, x1_wo, x2_wo)
end
pmap(WrapTaper, par);
dp = UnPatch(dir)
dd = dn - dp;

gap = zeros(size(dn,1),12); iy = 60;
tmp = hcat(dn[:,:,iy], gap, dp[:,:,iy], gap, dd[:,:,iy])
SeisPlot(tmp, hbox=5, wbox=9, style="wiggles", xcur=1.6, wiggle_trace_increment=6, ticksize=15, ylabel="Time (s)", xlabel="CDP X", dy=0.008, dx=1)
ax=gca(); ax[:set_yticks](collect(0.2:0.4:2.0))

tmp="/Users/wenlei/Desktop/input.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dn))); close(tmpfid);
tmp="/Users/wenlei/Desktop/output.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dp))); close(tmpfid);
tmp="/Users/wenlei/Desktop/diff.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dn-dp))); close(tmpfid);

# pscube < /Users/wenlei/Desktop/input.bin size1=4 size2=3 size3=1.5 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 perc=0.96 > /Users/wenlei/Desktop/input.eps
pscube < /Users/wenlei/Desktop/input.bin size1=4 size2=3 size3=1.5 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 bclip=0.96 wclip=-0.96 > /Users/wenlei/Desktop/input.eps
pscube < /Users/wenlei/Desktop/output.bin size1=4 size2=3 size3=1.5 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 bclip=0.96 wclip=-0.96 > /Users/wenlei/Desktop/output.eps
pscube < /Users/wenlei/Desktop/diff.bin size1=4 size2=3 size3=1.5 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 bclip=0.96 wclip=-0.96 > /Users/wenlei/Desktop/diff.eps


# ==============================================================================
#                            parabolic events
# ==============================================================================
(d, ext) = SeisParabEvents(nx1=200, nx2=200, dx1=5, dx2=5, nt=250, dt=0.008, p2=[0.4, -0.4, 0.6], amp=[1.0, -1.0, 0.7]);
snr = 1.0;
dn  = SeisAddNoise(d, snr, L=5);

root = homedir(); pin = join([root "/Desktop/randCP/"])
it_wl = 40; it_wo=20;
x1_wl = 40; x1_wo=20;
x2_wl = 40; x2_wo=20;
dir = patch(pin, dn, it_wl, it_wo, x1_wl, x1_wo, x2_wl, x2_wo)

R   = [15, 15, 15]; tol = 1e-4;
par = Array{Tuple{String, Vector{Int64}, Float64}}(length(dir));
for i = 1 : length(dir)
    par[i] = (dir[i], R, tol)
end
pmap(WrapTuckerAls, par)

par = Array{Tuple{String, Int64, Int64, Int64}}(length(dir));
for i = 1 : length(dir)
    par[i] = (dir[i], it_wo, x1_wo, x2_wo)
end
pmap(WrapTaper, par);
dp = UnPatch(dir)
dd = dn - dp;

tmp="/Users/wenlei/Desktop/input.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dn))); close(tmpfid);
tmp="/Users/wenlei/Desktop/output.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dp))); close(tmpfid);
tmp="/Users/wenlei/Desktop/diff.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dn-dp))); close(tmpfid);
pscube < /Users/wenlei/Desktop/input.bin size1=4 size2=3 size3=1.5 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 bclip=0.69 wclip=-0.69 > /Users/wenlei/Desktop/input.eps
pscube < /Users/wenlei/Desktop/output.bin size1=4 size2=3 size3=1.5 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 perc=96 bclip=0.69 wclip=-0.69 > /Users/wenlei/Desktop/output.eps
pscube < /Users/wenlei/Desktop/diff.bin size1=4 size2=3 size3=1.5 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 perc=96 bclip=0.69 wclip=-0.69 > /Users/wenlei/Desktop/diff.eps


gap = zeros(size(dn,1),12)
iy = 40
tmp = hcat(dn[:,:,iy], gap, dp[:,:,iy], gap, dd[:,:,iy])
SeisPlot(tmp, hbox=5, wbox=9, style="wiggles", xcur=1.4, wiggle_trace_increment=6, ticksize=15, ylabel="Time (s)", xlabel="CDP X", dy=0.008, dx=1)
ax=gca(); ax[:set_yticks](collect(0.2:0.4:2.0))
savefig("/Users/wenlei/Desktop/Ciy40.pdf")

iy = 100
tmp = hcat(dn[:,:,iy], gap, dp[:,:,iy], gap, dd[:,:,iy])
SeisPlot(tmp, hbox=5, wbox=9, style="wiggles", xcur=1.4, wiggle_trace_increment=6, ticksize=15, ylabel="Time (s)", xlabel="CDP X", dy=0.008, dx=1)
ax=gca(); ax[:set_yticks](collect(0.2:0.4:2.0))
savefig("/Users/wenlei/Desktop/Ciy40.pdf")

iy = 160
tmp = hcat(dn[:,:,iy], gap, dp[:,:,iy], gap, dd[:,:,iy])
SeisPlot(tmp, hbox=5, wbox=9, style="wiggles", xcur=1.4, wiggle_trace_increment=6, ticksize=15, ylabel="Time (s)", xlabel="CDP X", dy=0.008, dx=1)
ax=gca(); ax[:set_yticks](collect(0.2:0.4:2.0))
savefig("/Users/wenlei/Desktop/Ciy160.pdf")


# ==============================================================================
#                            test MSSA on irregular events
# ==============================================================================
using PyPlot, Seismic, SAIGtensor

(d, ext) = Seismic.SeisMixEvents3D(nt=250, nx1=200, nx2=200, dt=0.008, dx1=10, dx2=10,
                           v1=[6000. ,2000,-2700.,-2700.,6000. ],
                           v2=[14000.,4000, 4000., 4000.,14000.],
                           tau=[0.1  ,0.4 ,  0.9 ,  1.2 ,1.4   ],
                           amp=[1.0  ,0.9 ,  1.0 , -1.0 ,0.6   ],
                           apx1=[0.0 ,0.0 ,  0.0 ,  0.0 ,0.0   ],
                           apx2=[0.0 ,0.0 ,  0.0 ,  0.0 ,0.0   ],
                           eventType=['l', 'h', 'p', 'p', 'l'  ],
                           L=13, sigma=0.1)

snr = 1.0; dn = SeisAddNoise(d, snr, L=5);
10 * log10(vecnorm(d)^2 / vecnorm(d-dn)^2)
# patching
root = homedir(); pin = join([root "/Desktop/MSSA/"])
it_wl = 60; it_wo=20;
x1_wl = 30; x1_wo=15;
x2_wl = 30; x2_wo=15;
dir = patch(pin, dn, it_wl, it_wo, x1_wl, x1_wo, x2_wl, x2_wo)

# decomposition
rk   = 10; dt=0.008f0; flow=1.0f0; fhigh=40.0f0;
par = Array{Tuple{String, Float32, Int64, Float32, Float32}}(length(dir));
for i = 1 : length(dir)
    par[i] = (dir[i], dt, rk, flow, fhigh)
end
pmap(WrapDenoiseOnePatchMSSA3D, par)

# tapering and unpatching
par = Array{Tuple{String, Int64, Int64, Int64}}(length(dir));
for i = 1 : length(dir)
    par[i] = (dir[i], it_wo, x1_wo, x2_wo)
end
pmap(WrapTaper, par);
dp = UnPatch(dir)
dd = dn - dp;

gap = zeros(size(dn,1),12); iy = 60;
tmp = hcat(dn[:,:,iy], gap, dp[:,:,iy], gap, dd[:,:,iy])
SeisPlot(tmp, hbox=5, wbox=9, style="wiggles", xcur=1.6, wiggle_trace_increment=6, ticksize=15, ylabel="Time (s)", xlabel="CDP X", dy=0.008, dx=1)
ax=gca(); ax[:set_yticks](collect(0.2:0.4:2.0))

tmp="/Users/wenlei/Desktop/clean.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(d))); close(tmpfid);
tmp="/Users/wenlei/Desktop/minput.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dn))); close(tmpfid);
tmp="/Users/wenlei/Desktop/moutput.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dp))); close(tmpfid);
tmp="/Users/wenlei/Desktop/mdiff.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dd))); close(tmpfid);
iy = 100;
tmp="/Users/wenlei/Desktop/msection.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dp[:,:,iy]))); close(tmpfid);
tmp="/Users/wenlei/Desktop/msecdiff.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dd[:,:,iy]))); close(tmpfid);


# pscube < /Users/wenlei/Desktop/input.bin size1=4 size2=3 size3=1.5 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 perc=0.96 > /Users/wenlei/Desktop/input.eps
pscube < /Users/wenlei/Desktop/minput.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 bclip=0.96 wclip=-0.96 > /Users/wenlei/Desktop/minput.eps
pscube < /Users/wenlei/Desktop/moutput.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 bclip=0.96 wclip=-0.96 > /Users/wenlei/Desktop/moutput.eps
pscube < /Users/wenlei/Desktop/mdiff.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 bclip=0.96 wclip=-0.96 > /Users/wenlei/Desktop/mdiff.eps
pscube < /Users/wenlei/Desktop/clean.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 bclip=0.96 wclip=-0.96 > /Users/wenlei/Desktop/clean.eps


psimage < /Users/wenlei/Desktop/msection.bin height=3.0 width=2.0 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=60 d2num=40 f2num=40 label2="Crossline" xbox=0.5 ybox=0.5 bclip=0.96 wclip=-0.96 > /Users/wenlei/Desktop/msection.eps
psimage < /Users/wenlei/Desktop/msecdiff.bin height=3.0 width=2.0 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=60 d2num=40 f2num=40 label2="Crossline" xbox=0.5 ybox=0.5 bclip=0.96 wclip=-0.96 > /Users/wenlei/Desktop/msecdiff.eps

pswigb < /Users/wenlei/Desktop/msection.bin hbox=8.0 wbox=6.0 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=60 d2num=40 f2num=40 label2="Crossline" xbox=0.5 ybox=0.5 bclip=0.96 wclip=-0.96 > /Users/wenlei/Desktop/msection.eps
pswigb < /Users/wenlei/Desktop/msecdiff.bin hbox=8.0 wbox=6.0 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=60 d2num=40 f2num=40 label2="Crossline" xbox=0.5 ybox=0.5 bclip=0.96 wclip=-0.96 > /Users/wenlei/Desktop/msecdiff.eps

SeisPlot(dp[:,20:60,iy], hbox=5, wbox=4, style="wiggles", xcur=1.6, wiggle_trace_increment=1, ticksize=15, ylabel="Time (s)", xlabel="CDP X", dy=0.008, dx=1)

# ==============================================================================
#                           irregular events
# ==============================================================================
using PyPlot, Seismic, SAIGtensor
(d, ext) = Seismic.SeisMixEvents3D(nt=250, nx1=200, nx2=200, dt=0.008, dx1=10, dx2=10,
                           v1=[6000. ,2000,-2700.,-2700.,6000. ],
                           v2=[14000.,4000, 4000., 4000.,14000.],
                           tau=[0.1  ,0.4 ,  0.9 ,  1.2 ,1.4   ],
                           amp=[1.0  ,0.9 ,  1.0 , -1.0 ,0.6   ],
                           apx1=[0.0 ,0.0 ,  0.0 ,  0.0 ,0.0   ],
                           apx2=[0.0 ,0.0 ,  0.0 ,  0.0 ,0.0   ],
                           eventType=['l', 'h', 'p', 'p', 'l'  ],
                           L=13, sigma=0.1)

snr = 1.5; dn = SeisAddNoise(d, snr, L=5);
tmp = 10 * log10(vecnorm(d)^2 / vecnorm(dn-d)^2)

root = homedir(); pin = join([root "/Desktop/randCP/"])
it_wl = 30; it_wo=15;
x1_wl = 30; x1_wo=15;
x2_wl = 30; x2_wo=15;
dir1 = patch(pin, dn, it_wl, it_wo, x1_wl, x1_wo, x2_wl, x2_wo)

R   = [10, 10, 10]; tol = 1e-4;
par = Array{Tuple{String, Vector{Int64}, Float64}}(length(dir1));
for i = 1 : length(dir1)
    par[i] = (dir1[i], R, tol)
end
pmap(WrapTuckerAls, par)

par = Array{Tuple{String, Int64, Int64, Int64}}(length(dir1));
for i = 1 : length(dir1)
    par[i] = (dir1[i], it_wo, x1_wo, x2_wo)
end
pmap(WrapTaper, par);
dp1 = UnPatch(dir1)
dd1 = dn - dp1;

tmp="/Users/wenlei/Desktop/output.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dp1))); close(tmpfid);
tmp="/Users/wenlei/Desktop/diff.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dd1))); close(tmpfid);

# pscube < /Users/wenlei/Desktop/input.bin size1=4 size2=3 size3=1.5 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 perc=0.96 > /Users/wenlei/Desktop/input.eps
pscube < /Users/wenlei/Desktop/output.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 bclip=0.96 wclip=-0.96 > /Users/wenlei/Desktop/Routput.eps
pscube < /Users/wenlei/Desktop/diff.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="Crossline" n3=200 d3num=40 f3num=40 label3="Inline" xbox=0.5 ybox=0.5 bclip=1.5 wclip=-1.5 > /Users/wenlei/Desktop/Rdiff.eps

itl = 126; itu=200; ixl=5; ixu=35; iy=1;
SeisPlot(  d[itl:itu,ixl:ixu,iy], hbox=4, wbox=3, style="wiggles", xcur=1.2, wiggle_trace_increment=1, labelsize=12, ticksize=12, ylabel="Time (s)", xlabel="CDP X", oy=1.0, dy=0.008, dx=1); tight_layout(); ax=gca(); ax[:set_yticks](collect(1.1:0.1:1.5)); savefig("/Users/wenlei/Desktop/r1.pdf");
SeisPlot( dn[itl:itu,ixl:ixu,iy], hbox=4, wbox=3, style="wiggles", xcur=1.2, wiggle_trace_increment=1, labelsize=12, ticksize=12, ylabel="Time (s)", xlabel="CDP X", oy=1.0, dy=0.008, dx=1); tight_layout(); ax=gca(); ax[:set_yticks](collect(1.1:0.1:1.5)); savefig("/Users/wenlei/Desktop/r2.pdf");
SeisPlot( dp[itl:itu,ixl:ixu,iy], hbox=4, wbox=3, style="wiggles", xcur=1.2, wiggle_trace_increment=1, labelsize=12, ticksize=12, ylabel="Time (s)", xlabel="CDP X", oy=1.0, dy=0.008, dx=1); tight_layout(); ax=gca(); ax[:set_yticks](collect(1.1:0.1:1.5)); savefig("/Users/wenlei/Desktop/r3.pdf");
SeisPlot( dd[itl:itu,ixl:ixu,iy], hbox=4, wbox=3, style="wiggles", xcur=1.2, wiggle_trace_increment=1, labelsize=12, ticksize=12, ylabel="Time (s)", xlabel="CDP X", oy=1.0, dy=0.008, dx=1); tight_layout(); ax=gca(); ax[:set_yticks](collect(1.1:0.1:1.5)); savefig("/Users/wenlei/Desktop/r4.pdf");
SeisPlot(dp1[itl:itu,ixl:ixu,iy], hbox=4, wbox=3, style="wiggles", xcur=1.8, wiggle_trace_increment=1, labelsize=12, ticksize=12, ylabel="Time (s)", xlabel="CDP X", oy=1.0, dy=0.008, dx=1); tight_layout(); ax=gca(); ax[:set_yticks](collect(1.1:0.1:1.5)); savefig("/Users/wenlei/Desktop/r5.pdf");
SeisPlot(dd1[itl:itu,ixl:ixu,iy], hbox=4, wbox=3, style="wiggles", xcur=1.2, wiggle_trace_increment=1, labelsize=12, ticksize=12, ylabel="Time (s)", xlabel="CDP X", oy=1.0, dy=0.008, dx=1); tight_layout(); ax=gca(); ax[:set_yticks](collect(1.1:0.1:1.5)); savefig("/Users/wenlei/Desktop/r6.pdf");

gap = zeros(size(dn,1),12)
iy = 40
tmp = hcat(dn[:,:,iy], gap, dp[:,:,iy], gap, dd[:,:,iy])
SeisPlot(tmp, hbox=5, wbox=9, style="wiggles", xcur=1.4, wiggle_trace_increment=6, ticksize=15, ylabel="Time (s)", xlabel="CDP X", dy=0.008, dx=1)
ax=gca(); ax[:set_yticks](collect(0.2:0.4:2.0))
