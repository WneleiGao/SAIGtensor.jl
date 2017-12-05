using DSP

function Ricker(; dt::Real=0.002, f0::Real=20.0)
    nw = 2.0/(f0*dt)
    nc = floor(Int, nw/2)
    t = dt*collect(-nc:1:nc)
    b = (pi*f0*t).^2
    w = (1.-2.*b) .* exp.(-b)
end

function SeisMixEvents3D(; ot =0.0, dt =0.004, nt =500,
                           ox1=0.0, dx1=10.0 , nx1=200,
                           ox2=0.0, dx2=10.0 , nx2=200,
                           ox3=0.0, dx3=10.0 , nx3=1,
                           ox4=0.0, dx4=10.0 , nx4=1,
                           apx1=[0.,0.,0.], apx2=[0.,0.,0.],
                           v1  =[6000,2000,-2700], v2  =[14000,4000,4000],
                           tau =[0.1,0.4,0.9], amp=[1.0,-1.0,0.7], eventType=['l', 'p', 'h'], f0=20.0,
                           etp=Float64, pertube=true, sigma=0.075, L=13)

    mod(L,2) == 1 || error("L must be an odd number")
    w = Ricker(dt=dt,f0=f0)
    nf = nextpow2(nt)
    nw = length(w)
    t_delay = (nw-1)*dt/2     # can show the full wavelet
    w = vcat(w, zeros(nf-nw))
    W = fft(w)
    x1 = ox1 + collect(0:1:nx1-1)*dx1
    x2 = ox2 + collect(0:1:nx2-1)*dx2
    nevents = length(amp)
    D = zeros(Complex{etp}, nf, nx1, nx2)
    nfh = round(Int64, floor(nf/2)) + 1
    wrs = collect(0:1:nfh-1)*2*pi/(nf*dt)     # Frequency in rad/sec
    for ie = 1:nevents
        # generate random pertubation of traveltime
        if pertube
           f = hamming(L); f = f / sum(f);
           drop = floor(Int64, L/2)*2
           pt = conv2(f,f,randn(nx1+drop, nx2+drop))[drop+1:end-drop, drop+1:end-drop] * sigma
        else
           pt = zeros(nx1, nx2)
        end
        if eventType[ie] == 'l'    # linear
           p1 = 1/v1[ie]; p2 = 1/v2[ie];
           for ix2 = 1:nx2
               for ix1 = 1:nx1
                   for iw = 2:nfh-1
                       phase = wrs[iw] * (tau[ie] + p1*(x1[ix1]-apx1[ie])
                                                  + p2*(x2[ix2]-apx2[ie])
                                                  - t_delay + pt[ix1,ix2])
                       D[iw,ix1,ix2] += W[iw]*amp[ie]*exp(-im * phase)
                       D[nf-iw+2,ix1,ix2] = conj(D[iw,ix1,ix2])
                   end
               end
           end
        elseif eventType[ie] == 'p'  # parabolic
           p1 = sign(v1[ie]) / v1[ie]^2; p2 = sign(v2[ie]) / v2[ie]^2;
           for ix2 = 1:nx2
               for ix1 = 1:nx1
                   for iw = 2:nfh-1
                       phase = wrs[iw] * (tau[ie] + p1*(x1[ix1]-apx1[ie])^2
                                                  + p2*(x2[ix2]-apx2[ie])^2
                                                  - t_delay + pt[ix1,ix2])
                       D[iw,ix1,ix2] += W[iw]*amp[ie]*exp(-im * phase)
                       D[nf-iw+2,ix1,ix2] = conj(D[iw,ix1,ix2])
                   end
               end
           end
        elseif eventType[ie] == 'h'   # hyperbola
           p1 = sign(v1[ie]) / v1[ie]^2; p2 = sign(v2[ie]) / v2[ie]^2;
           for ix2 = 1:nx2
               for ix1 = 1:nx1
                   for iw = 2:nfh-1
                       phase = wrs[iw] * (sqrt(tau[ie]^2 + p1*(x1[ix1]-apx1[ie])^2
                                                         + p2*(x2[ix2]-apx2[ie])^2
                                                         - t_delay + pt[ix1,ix2]))
                       D[iw,ix1,ix2] += W[iw]*amp[ie]*exp(-im * phase)
                       D[nf-iw+2,ix1,ix2] = conj(D[iw,ix1,ix2])
                   end
               end
           end
        end
    end
    d = ifft(D,1)
    d = real(d[1:nt,:,:,:,:])

    if nx4 == 1 && nx3 == 1 && nx2 == 1
        sqdims = (3,4,5)
    elseif  nx4 == 1 && nx3 == 1
        sqdims = (4,5)
    end
    d = squeeze(d, sqdims)
    return d
end
