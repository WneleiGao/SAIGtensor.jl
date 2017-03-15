"""
d is a 3D array, it_wl is the window length along time axis, it_wo is the
overlapping of window.
"""

function patch{Ti<:AbstractFloat, Tn}(path::String, d::Array{Ti,3},
                                      it_wl::Tn, it_wo::Tn,
                                      x1_wl::Tn, x1_wo::Tn,
                                      x2_wl::Tn, x2_wo::Tn)
    it_ws = it_wl - it_wo
    x1_ws = x1_wl - x1_wo
    x2_ws = x2_wl - x2_wo
    (nt, n1, n2) = size(d)
    it_wl = it_wl > nt ? nt : it_wl
    x1_wl = x1_wl > n1 ? n1 : x1_wl
    x2_wl = x2_wl > n2 ? n2 : x2_wl
    it_nw = floor(Int64, nt/it_ws)
    x1_nw = floor(Int64, n1/x1_ws)
    x2_nw = floor(Int64, n2/x2_ws)
    # be carefull with the boundary
    if it_nw*it_ws + it_wo < nt
       it_nw = it_nw + 1
    end
    if x1_nw*x1_ws + x1_wo < n1
       x1_nw = x1_nw + 1
    end
    if x2_nw*x2_ws + x2_wo < n2
       x2_nw = x2_nw + 1
    end
    for i2 = 1 : x2_nw
        i2l = (i2-1)*x2_ws + 1
        i2u = i2l + x2_wl - 1
        if i2u > n2
           i2u = n2
        end
        for i1 = 1 : x1_nw
            i1l = (i1-1)*x1_ws + 1
            i1u = i1l + x1_wl - 1
            if i1u > n1
               i1u = n1
            end
            for it = 1 : it_nw
                itl = (it-1)*it_ws + 1
                itu = itl + it_wl - 1
                if itu > nt
                   itu = nt
                end
                pout = join([path "patch" "_" it "_" i1 "_" i2 ".bin"])
                fid = open(pout, "w")
                tmp = d[itl:itu, i1l:i1u, i2l:i2u]
                write(fid, Int32(nt) , Int32(n1) , Int32(n2))
                write(fid, Int32(itl), Int32(itu), Int32(i1l), Int32(i1u), Int32(i2l), Int32(i2u))
                write(fid, convert(Array{Float32}, vec(tmp)))
                close(fid)
            end
        end
    end
    dir = getPatchDir(path, it_nw, x1_nw, x2_nw)
    return dir
end

function getPatchDir(pin::String, it_nw::Int64, x1_nw::Int64, x2_nw::Int64)
    dir = Array{String}(it_nw*x1_nw*x2_nw)
    for i2 = 1 : x2_nw
        for i1 = 1 : x1_nw
            for it = 1 : it_nw
                idx= (i2-1)*x1_nw*it_nw + (i1-1)*it_nw + it
                dir[idx] = join([pin "patch" "_" it "_" i1 "_" i2 ".bin"])
            end
        end
    end
    return dir
end

function getOnePatch(pin::String; rflag=false)
    fid = open(pin, "r")
    nt  = read(fid, Int32); n1  = read(fid, Int32); n2 = read(fid, Int32);
    itl = read(fid, Int32); itu = read(fid, Int32); nt_w = itu - itl + 1;
    x1l = read(fid, Int32); x1u = read(fid, Int32); n1_w = x1u - x1l + 1;
    x2l = read(fid, Int32); x2u = read(fid, Int32); n2_w = x2u - x2l + 1;
    d   = reshape(read(fid,Float32,nt_w*n1_w*n2_w), nt_w, n1_w, n2_w)
    close(fid)
    if rflag
       return d, nt, n1, n2, itl, itu, x1l, x1u, x2l, x2u
    else
       return d
    end
end

function Taper{Ti<:Integer}(pin::Array{String,1}, it_wo::Ti, x1_wo::Ti, x2_wo::Ti)
    npatch = length(pin)
    for i = 1 : npatch
        Taper(pin[i], it_wo, x1_wo, x2_wo)
    end
    return nothing
end

function Taper{Ti<:Integer}(pin::String, it_wo::Ti, x1_wo::Ti, x2_wo::Ti)
    fid = open(pin, "r+")
    (d, nt, n1, n2, itl, itu, x1l, x1u, x2l, x2u) = getOnePatch(pin, rflag=true)
    nt_w = itu-itl+1;
    n1_w = x1u-x1l+1;
    n2_w = x2u-x2l+1;
    tapti = itl >  1  ? it_wo : 0
    taptf = itu <  nt ? it_wo : 0
    tap1i = x1l >  1  ? x1_wo : 0
    tap1f = x1u <  n1 ? x1_wo : 0
    tap2i = x2l >  1  ? x2_wo : 0
    tap2f = x2u <  n2 ? x2_wo : 0
    tt = 1.; t1 = 1.; t2 = 1.;
    for i2 = 1 : n2_w
        if i2 >= 1 && i2 <= tap2i
           t2 = 1. - cos(pi/2*((i2-1)/tap2i))
        end
        if i2 > tap2i && i2 <= n2_w-tap2f
           t2 = 1.
        end
        if i2 > n2_w-tap2f && i2 <= n2_w
           t2 = cos(pi/2*(i2-1-n2_w+tap2f)/tap2f)
        end
        for i1 = 1 : n1_w
            if i1 >= 1 && i1 <= tap1i
               t1 = 1. - cos(pi/2*((i1-1)/tap1i))
            end
            if i1 > tap1i && i1 <= n1_w-tap1f
               t1 = 1.
            end
            if i1 > n1_w-tap1f && i1 <= n1_w
               t1 = cos(pi/2*(i1-1-n1_w+tap1f)/tap1f)
            end
            for it = 1 : nt_w
                if it >= 1 && it <= tapti
                   tt = 1. - cos(pi/2*((it-1)/tapti))
                end
                if it > tapti && it <= nt_w-taptf
                   tt = 1.
                end
                if it > nt_w-taptf && it <= nt_w
                   tt = cos(pi/2*(it-1-nt_w+taptf)/taptf)
                end
                d[it,i1,i2] = d[it,i1,i2] * tt * t1 * t2
            end
        end
    end
    d = convert(Array{Float32,3}, d);
    pos = sizeof(Int32) * 9
    seek(fid, pos)
    write(fid, vec(d))
    close(fid)
    return nothing
end

function UnPatch(pin::Array{String,1})
    (d1, nt, n1, n2, itl, itu, x1l, x1u, x2l, x2u)= getOnePatch(pin[1], rflag=true);
    d = zeros(Float32, nt, n1, n2)
    d[itl:itu, x1l:x1u, x2l:x2u] = d1[:,:,:]
    for i = 2 : length(pin)
        (d1, nt, n1, n2, itl, itu, x1l, x1u, x2l, x2u)= getOnePatch(pin[i], rflag=true);
        d[itl:itu, x1l:x1u, x2l:x2u] = d[itl:itu, x1l:x1u, x2l:x2u] + d1[:,:,:]
    end
    return d
end

function WrapRandCpAls(par::Tuple{String, Int64})
    path = par[1];
    R    = par[2];
    dp = getOnePatch(path)
    X = tensor(3, collect(size(dp)), convert(Array{Float64,3}, dp));
    (lambda, A) = randCpAls_simplify(X, R);
    X = cptensor(3, collect(size(dp)), R, lambda, A);
    X = cp2tensor(X)
    dp = convert(Array{Float32,1}, vec(X.D))
    fid = open(path, "r+"); pos = sizeof(Int32)*9;
    seek(fid, pos); write(fid, dp); close(fid)
end

funtion WrapTaper(par::Tuple{String, Int64, Int64, Int64})
    dir = par[1]
    it_wo = par[2]
    x1_wo = par[3]
    x2_wo = par[4]
    Taper(dir, it_wo, x1_wo, x2_wo)
    return nothing
end
