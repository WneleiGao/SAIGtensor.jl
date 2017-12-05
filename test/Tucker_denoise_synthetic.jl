# ==============================================================================
#                            linear events
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
gap = zeros(size(dn,1),12)

# plotting
iy = 60
tmp = hcat(dn[:,:,iy], gap, dp[:,:,iy], gap, dd[:,:,iy])
SeisPlot(tmp, hbox=5, wbox=9, style="wiggles", xcur=1.4, wiggle_trace_increment=6, ticksize=15, ylabel="Time (s)", xlabel="CDP X", dy=0.008, dx=1)
ax=gca(); ax[:set_yticks](collect(0.2:0.4:2.0))

iy = 160
tmp = hcat(dn[:,:,iy], gap, dp[:,:,iy], gap, dd[:,:,iy])
SeisPlot(tmp, hbox=5, wbox=9, style="wiggles", xcur=1.4, wiggle_trace_increment=6, ticksize=15, ylabel="Time (s)", xlabel="CDP X", dy=0.008, dx=1)
ax=gca(); ax[:set_yticks](collect(0.2:0.4:2.0))

# ==============================================================================
#                            parabolic events
# ==============================================================================
(d, ext) = SeisParabEvents(nx1=200, nx2=200, dx1=5, dx2=5, nt=250, dt=0.008, p2=[0.4, 0.4, 0.4], amp=[1.0, -1.0, 0.7]);
snr = 0.5;
dn  = SeisAddNoise(d, snr, L=5);

root = homedir(); pin = join([root "/Desktop/randCP/"])
it_wl = 60; it_wo=30;
x1_wl = 60; x1_wo=30;
x2_wl = 60; x2_wo=30;
dir = patch(pin, dn, it_wl, it_wo, x1_wl, x1_wo, x2_wl, x2_wo)

R   = [20, 20, 20]; tol = 1e-4;
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
#                           irregular events
# ==============================================================================
using PyPlot, Seismic, SAIGtensor
(d, ext) = SeisMixEvents3D(nt=250, nx1=200, nx2=200, dt=0.008, dx1=10, dx2=10,
                           v1=[6000. ,2000,-2700.,-2700.,6000. ],
                           v2=[14000.,4000, 4000., 4000.,14000.],
                           tau=[0.1  ,0.4 ,  0.9 ,  1.2 ,1.4   ],
                           amp=[1.0  ,0.9 ,  1.0 , -1.0 ,0.6   ],
                           apx1=[0.0 ,0.0 ,  0.0 ,  0.0 ,0.0   ],
                           apx2=[0.0 ,0.0 ,  0.0 ,  0.0 ,0.0   ],
                           eventType=['l', 'h', 'p', 'p', 'l'  ],
                           L=13, sigma=0.1)

snr = 1.0;
dn = SeisAddNoise(d, snr, L=5);

root = homedir(); pin = join([root "/Desktop/Irregular/"])
it_wl = 60; it_wo=30;
x1_wl = 60; x1_wo=30;
x2_wl = 60; x2_wo=30;
dir = patch(pin, dn, it_wl, it_wo, x1_wl, x1_wo, x2_wl, x2_wo)

R   = [23, 23, 23]; tol = 1e-4;
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
gap = zeros(size(dn,1),12)

iy = 40
tmp = hcat(dn[:,:,iy], gap, dp[:,:,iy], gap, dd[:,:,iy])
SeisPlot(tmp, hbox=5, wbox=9, style="wiggles", xcur=1.4, wiggle_trace_increment=6, ticksize=15, ylabel="Time (s)", xlabel="CDP X", dy=0.008, dx=1)
ax=gca(); ax[:set_yticks](collect(0.2:0.4:2.0))

iy = 100
tmp = hcat(dn[:,:,iy], gap, dp[:,:,iy], gap, dd[:,:,iy])
SeisPlot(tmp, hbox=5, wbox=9, style="wiggles", xcur=1.4, wiggle_trace_increment=6, ticksize=15, ylabel="Time (s)", xlabel="CDP X", dy=0.008, dx=1)
ax=gca(); ax[:set_yticks](collect(0.2:0.4:2.0))


iy = 160
tmp = hcat(dn[:,:,iy], gap, dp[:,:,iy], gap, dd[:,:,iy])
SeisPlot(tmp, hbox=5, wbox=9, style="wiggles", xcur=1.4, wiggle_trace_increment=1, ticksize=15, ylabel="Time (s)", xlabel="CDP X", dy=0.008, dx=1)
ax=gca(); ax[:set_yticks](collect(0.2:0.4:2.0))
