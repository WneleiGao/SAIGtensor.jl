# ============example with synthetic data==============
using PyPlot, Seismic, SAIGtensor
(d,e) = SeisLinearEvents(nt= 250, dt=0.008, nx1=200, nx2=200,
                         tau=[0.4, 0.8, 1.2, 1.6],
                         p1=[ 0.0003,  0.0004, -0.0004, -0.00027],
                         p2=[-0.0002, -0.0004,  0.0001,  0.0002 ],
                         p3=[0., 0., -0., -0.],
                         p4=[0., 0., -0., -0.],
                         amp=[1., 1.2, -0.5, -0.9])

snr = 0.5;
dn = SeisAddNoise(d, snr, L=5);

pin = "/Users/wenlei/Desktop/randCP/synthetic/"
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
savefig("/Users/wenlei/Desktop/Liy40.pdf")

iy = 100
tmp = hcat(dn[:,:,iy], gap, dp[:,:,iy], gap, dd[:,:,iy])
SeisPlot(tmp, hbox=5, wbox=9, style="wiggles", xcur=1.4, wiggle_trace_increment=6, ticksize=15, ylabel="Time (s)", xlabel="CDP X", dy=0.008, dx=1)
ax=gca(); ax[:set_yticks](collect(0.2:0.4:2.0))
savefig("/Users/wenlei/Desktop/Liy40.pdf")

iy = 160
tmp = hcat(dn[:,:,iy], gap, dp[:,:,iy], gap, dd[:,:,iy])
SeisPlot(tmp, hbox=5, wbox=9, style="wiggles", xcur=1.4, wiggle_trace_increment=6, ticksize=15, ylabel="Time (s)", xlabel="CDP X", dy=0.008, dx=1)
ax=gca(); ax[:set_yticks](collect(0.2:0.4:2.0))
savefig("/Users/wenlei/Desktop/Liy160.pdf")


# # =====decompose as a whole patch============
# N = 3;
# I = collect(size(ds));
# X = tensor(N, I, ds);
# R = 120;
# (lambda2, A2) = cpals(X, R, pflag=true);
# X = cptensor(N, I, R, lambda2, A2);
# X = cp2tensor(X)
# dp = X.D;
# gap = zeros(250,1); r1=1:2:500; r2=1:10:100;
# dd = ds - dp;
#
# i3 = 1; d1 = hcat(ds[r1,r2,i3], gap, dp[r1,r2,i3], gap, dd[r1,r2,i3])
# i3 = 25; d2 = hcat(ds[r1,r2,i3], gap, dp[r1,r2,i3], gap, dd[r1,r2,i3])
# i3 = 50; d3 = hcat(ds[r1,r2,i3], gap, dp[r1,r2,i3], gap, dd[r1,r2,i3])
#
# fig = figure("1", figsize=(6, 6)); tmp = 13;
# subplot(3,1,1); SeisPlot(d1,fignum="1", style="wiggles", dy=0.004, dx=10)
# ax=gca(); ax[:set_yticks]([0.3, 0.6, 0.9])
# ax[:set_xticks]([])
# ax[:text](-25., -0.05, "a)", fontsize=tmp, fontweight="bold")
# subplot(3,1,2); SeisPlot(d2,fignum="1", style="wiggles", dy=0.004, dx=10, ylabel="Time (s)")
# ax=gca(); ax[:set_yticks]([0.3, 0.6, 0.9])
# ax[:set_xticks]([])
# ax[:text](-25., -0.05, "b)", fontsize=tmp, fontweight="bold")
# subplot(3,1,3); SeisPlot(d3,fignum="1", style="wiggles", dy=0.004, dx=10, xlabel="CMP x number")
# ax=gca(); ax[:set_yticks]([0.3, 0.6, 0.9])
# ax[:text](-25., -0.05, "c)", fontsize=tmp, fontweight="bold")

# ========for curve events=================
(d, e) = SeisParabEvents(nx1=200, nx2=200, dx1=5, dx2=5, nt=250, dt=0.008, p2=[0.4, 0.4, 0.4], amp=[1.0, -1.0, 0.7]);
snr = 0.5;
dn  = SeisAddNoise(d, snr, L=5);

pin = "/Users/wenlei/Desktop/randCP/synthetic/"
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


#
# N = 3;
# I = collect(size(ds));
# X = tensor(N, I, ds);
# R = 120;
# (lambda2, A2) = cpals(X, R, pflag=true);
# X = cptensor(N, I, R, lambda2, A2);
# X = cp2tensor(X)
#
# dp = X.D;
# gap = zeros(250,1); r1=1:2:500; r2=10:10:90;
# dd = ds - dp;
#
# i3 = 1; d1 = hcat(ds[r1,r2,i3], gap, dp[r1,r2,i3], gap, dd[r1,r2,i3])
# i3 = 40; d2 = hcat(ds[r1,r2,i3], gap, dp[r1,r2,i3], gap, dd[r1,r2,i3])
# i3 = 80; d3 = hcat(ds[r1,r2,i3], gap, dp[r1,r2,i3], gap, dd[r1,r2,i3])
#
# fig = figure("1", figsize=(6, 6)); tmp = 13;
# subplot(3,1,1); SeisPlot(d1,fignum="1", style="wiggles", dy=0.004, dx=10)
# ax=gca(); ax[:set_yticks]([0.3, 0.6, 0.9])
# ax[:set_xticks]([])
# ax[:text](-25., -0.05, "a)", fontsize=tmp, fontweight="bold")
# subplot(3,1,2); SeisPlot(d2,fignum="1", style="wiggles", dy=0.004, dx=10, ylabel="Time (s)")
# ax=gca(); ax[:set_yticks]([0.3, 0.6, 0.9])
# ax[:set_xticks]([])
# ax[:text](-25., -0.05, "b)", fontsize=tmp, fontweight="bold")
# subplot(3,1,3); SeisPlot(d3,fignum="1", style="wiggles", dy=0.004, dx=10, xlabel="CMP x number")
# ax=gca(); ax[:set_yticks]([0.3, 0.6, 0.9])
# ax[:text](-25., -0.05, "c)", fontsize=tmp, fontweight="bold")

# ===========================for data with irregular events=====================
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

snr = 0.5;
dn = SeisAddNoise(d, snr, L=1);

pin = "/Users/wenlei/Desktop/randCP/synthetic/"
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
savefig("/Users/wenlei/Desktop/Iiy40.pdf")

iy = 100
tmp = hcat(dn[:,:,iy], gap, dp[:,:,iy], gap, dd[:,:,iy])
SeisPlot(tmp, hbox=5, wbox=9, style="wiggles", xcur=1.4, wiggle_trace_increment=6, ticksize=15, ylabel="Time (s)", xlabel="CDP X", dy=0.008, dx=1)
ax=gca(); ax[:set_yticks](collect(0.2:0.4:2.0))
savefig("/Users/wenlei/Desktop/Iiy100.pdf")

iy = 160
tmp = hcat(dn[:,:,iy], gap, dp[:,:,iy], gap, dd[:,:,iy])
SeisPlot(tmp, hbox=5, wbox=9, style="wiggles", xcur=1.4, wiggle_trace_increment=6, ticksize=15, ylabel="Time (s)", xlabel="CDP X", dy=0.008, dx=1)
ax=gca(); ax[:set_yticks](collect(0.2:0.4:2.0))
savefig("/Users/wenlei/Desktop/Iiy160.pdf")
