"""
   tensor completion based on Parallel matrix factorization
"""
function TensorCP(dobs::Matrix{Tv}, dt::Tv, known::Vector{Ti}, I::Vector{Ti},
                  rank::Vector{Ti}, fmin::Tv, fmax::Tv; alpha=0.9f0, maxIter=100) where {Ti<:Int64, Tv<:Float32}
    N = length(I)
    nt= size(dobs,1); nf = nextpow2(nt);

    ilow = floor(Int64, fmin*dt*nf)+1
    ilow = ilow > 2 ? ilow : 2
    ihigh = floor(Int64, fmax*dt*nf)+1
    ihigh = ihigh < nf/2+1 ? ihigh : nf/2+1

    # padding data with zeros
    if nf > nt
       dobs = vcat(dobs, zeros(eltype(dobs), nf-nt, size(dobs,2)))
    end
    fobs = fft(dobs, 1)

    # pre-allocate memory for computation
    dim  = vcat(nf, collect(I)); etp = eltype(fobs);
    dest = zeros(eltype(fobs), dim...)
    CI = zeros(Int64, N); total = prod(I);
    X = Vector{Matrix{etp}}(N); Xt= Vector{Matrix{etp}}(N);
    Y = Vector{Matrix{etp}}(N); Yt= Vector{Matrix{etp}}(N);
    for i = 1 : N
        CI[i] = ceil(Int64, total/I[i])
        X[i]  = zeros(etp,I[i]   ,rank[i]); Xt[i]=zeros(etp,rank[i],I[i]   );
        Y[i]  = zeros(etp,rank[i],CI[i]  ); Yt[i]=zeros(etp,CI[i]  ,rank[i]);
    end
    Z = zeros(etp, I...)

    # do the reconstruction
    for i = ilow : ihigh
        slice = view(fobs, i, :)
        Z = TensorCP(slice, known, Z, X, Xt, Y, Yt, N, I, CI, rank, alpha, maxIter=maxIter)
        dest[i,:,:,:,:] = Z
        dest[nf-i+2,:,:,:,:] = conj(Z)
    end
    ifft!(dest, 1)
    return real(dest[1:nt,:,:,:,:])
end

"""
   tensor completion for just one frequency component
"""
function TensorCP(dobs::SubArray, known::Vector{Ti}, Z::Array{Tv},
                  X::Vector{Matrix{Tv}}, Xt::Vector{Matrix{Tv}},
                  Y::Vector{Matrix{Tv}}, Yt::Vector{Matrix{Tv}},
                  N::Ti, I::Vector{Ti}, CI::Vector{Ti}, rank::Union{Ti,Vector{Ti}},
                  alpha::Float32; maxIter=50, weight=0.25f0, convergence=false) where {Ti<:Int64, Tv<:Complex64}
    # checking input
    if length(rank) == 1
       rank = rank * ones(typeof(rank), N)
    end
    length(rank) == N || error("size of rank does not match")
    if length(weight) != N
       weight = Float32(1/N) * ones(Float32, N)
    end
    # normalize input data and initialization
    etp    = eltype(dobs)
    l2norm = vecnorm(dobs);
    eptl2norm = sqrt(l2norm^2 * prod(I) / length(known));
    for i = 1 : N
        rand!(X[i]); scaler=eptl2norm^(I[i] /(I[i]+CI[i]))/vecnorm(X[i]); X[i]=X[i]*scaler;
        rand!(Y[i]); scaler=eptl2norm^(CI[i]/(I[i]+CI[i]))/vecnorm(Y[i]); Y[i]=Y[i]*scaler;
        transpose!(Xt[i],X[i]); conj!(Xt[i]);
        transpose!(Yt[i],Y[i]); conj!(Yt[i]);
    end
    fill!(Z, zero(etp)); nobs = length(known);
    for i = 1 : nobs
        Z[known[i]] = dobs[i]
    end

    slope = (alpha-1.0f0)/(maxIter-1)
    for iter = 1 : maxIter
        # update X[n] and Y[n]
        for n = 1 : N
            Zn = matricization(Z,I,n); A_mul_B!(X[n],Zn,Yt[n])  ; transpose!(Xt[n],X[n]); conj!(Xt[n]); # update X[n]
            Xsq= Xt[n]*X[n]          ; Y[n]=pinv(Xsq)*(Xt[n]*Zn); transpose!(Yt[n],Y[n]); conj!(Yt[n]); # update Y[n]
        end
        # update Z
        fill!(Z, zero(etp))
        for n = 1 : N
            Zn = X[n]*Y[n];
            Zn = unmatricization(Zn, I, n)
            Z  = Z + weight[n]*Zn
        end
        if convergence==true
           relerr   = vecnorm(Z[known]-dobs) / vecnorm(dobs)
           println("relative error: $relerr")
        end
        # insert observations
        beta = slope*(iter-1) + 1
        for i = 1 : nobs
            Z[known[i]] =(1-beta)*Z[known[i]] + beta*dobs[i]
        end
    end
    return Z
end

"""
   parallel matrix factorization for tensor completion
"""
function PMF(dobs::Vector{Tv}, known::Vector{Ti},
             N::Ti, I::Vector{Ti}, rank::Union{Ti,Vector{Ti}},
             alpha::Float64; maxIter=100, weight=0.25, groundTrue=0.0) where {Ti<:Int64, Tv<:Number}

    CI = zeros(Int64, N)
    total = prod(I)
    for i = 1 : N
        CI[i] = round(Int64, total/I[i])
    end

    # checking input
    if length(rank) == 1
       rank = rank * ones(Int64, N)
    end
    length(rank) == N || error("size of rank does not match")
    if length(weight) != N
       weight = 1/N * ones(Float64, N)
    end

    # normalize input data and initialization
    etp    = eltype(dobs)
    l2norm = vecnorm(dobs)
    eptl2norm = sqrt(l2norm^2 * prod(I) / length(known));
    X = Vector{Matrix{Complex64}}(N); Xt = Vector{Matrix{Complex64}}(4);
    Y = Vector{Matrix{Complex64}}(N); Yt = Vector{Matrix{Complex64}}(4);
    for i = 1 : N
        X[i]=rand(etp,I[i]   ,rank[i]); scaler=eptl2norm^(I[i] /(I[i]+CI[i]))/vecnorm(X[i]); X[i]=X[i]*scaler;
        Y[i]=rand(etp,rank[i],CI[i]  ); scaler=eptl2norm^(CI[i]/(I[i]+CI[i]))/vecnorm(Y[i]); Y[i]=Y[i]*scaler;
        Xt[i]=zeros(etp,rank[i],I[i]   ); transpose!(Xt[i], X[i]); conj!(Xt[i]);
        Yt[i]=zeros(etp,CI[i]  ,rank[i]); transpose!(Yt[i], Y[i]); conj!(Yt[i]);
    end
    Z = zeros(etp, I...)
    Z[known] = dobs[:]

    for iter = 1 : maxIter
        # update X[n] and Y[n]
        for n = 1 : N
            Zn = matricization(Z,I,n); X[n]=Zn*Yt[n]; transpose!(Xt[n],X[n]); conj!(Xt[n]); # update X[n]
            Xsq= Xt[n]*X[n];  Y[n]= pinv(Xsq)*(Xt[n]*Zn); transpose!(Yt[n],Y[n]); conj!(Yt[n]); # update Y[n]
        end

        # update Z
        fill!(Z, zero(etp))
        for n = 1 : N
            Zn = X[n]*Y[n];
            Zn = unmatricization(Zn, I, n)
            Z  = Z + weight[n]*Zn
        end

        if groundTrue != 0
           relerr = vecnorm(Z - groundTrue) / vecnorm(groundTrue)
           println("relative error: $relerr")
        end

        # insert observations, alpha controal the contribution from observations
        Z[known] = (1-alpha)*Z[known] + alpha*dobs
    end
    return Z
end

# function TuckerCP(dobs::Matrix{Tv}, dt::Tv, known::Vector{Ti}, I::Vector{Ti},
#                   rank::Vector{Ti}, fmin::Tv, fmax::Tv; alpha=0.9f0, maxIter=100) where {Ti<:Int64, Tv<:Float32}
#     N = length(I)
#     nt= size(dobs,1); nf = nextpow2(nt);
#
#     ilow = floor(Int64, fmin*dt*nf)+1
#     ilow = ilow > 2 ? ilow : 2
#     ihigh = floor(Int64, fmax*dt*nf)+1
#     ihigh = ihigh < nf/2+1 ? ihigh : nf/2+1
#
#     # padding data with zeros
#     if nf > nt
#        dobs = vcat(dobs, zeros(eltype(dobs), nf-nt, size(dobs,2)))
#     end
#     fobs = fft(dobs, 1)
#
#     # pre-allocate memory for computation
#     dim  = vcat(nf, collect(I)); etp = eltype(fobs);
#     dest = zeros(eltype(fobs), dim...)
#     CI = zeros(Int64, N); total = prod(I);
#     X = Vector{Matrix{etp}}(N); Xt= Vector{Matrix{etp}}(N);
#     Y = Vector{Matrix{etp}}(N); Yt= Vector{Matrix{etp}}(N);
#     for i = 1 : N
#         CI[i] = ceil(Int64, total/I[i])
#         X[i]  = zeros(etp,I[i]   ,rank[i]); Xt[i]=zeros(etp,rank[i],I[i]   );
#         Y[i]  = zeros(etp,rank[i],CI[i]  ); Yt[i]=zeros(etp,CI[i]  ,rank[i]);
#     end
#     Z = zeros(etp, I...)
#
#     # do the reconstruction
#     for i = ilow : ihigh
#         slice = view(fobs, i, :)
#         Z = TensorCP(slice, known, Z, X, Xt, Y, Yt, N, I, CI, rank, alpha, maxIter=maxIter)
#         dest[i,:,:,:,:] = Z
#         dest[nf-i+2,:,:,:,:] = conj(Z)
#         println("$i")
#     end
#     ifft!(dest, 1)
#     return real(dest[1:nt,:,:,:,:])
# end
#
# """
#    tensor completion for just one frequency component
# """
# function TensorCP(dobs::SubArray, known::Vector{Ti}, Z::Array{Tv},
#                   X::Vector{Matrix{Tv}}, Xt::Vector{Matrix{Tv}},
#                   Y::Vector{Matrix{Tv}}, Yt::Vector{Matrix{Tv}},
#                   N::Ti, I::Vector{Ti}, CI::Vector{Ti}, rank::Union{Ti,Vector{Ti}},
#                   alpha::Float32; maxIter=50, weight=0.25f0, convergence=false) where {Ti<:Int64, Tv<:Complex64}
#     # checking input
#     if length(rank) == 1
#        rank = rank * ones(typeof(rank), N)
#     end
#     length(rank) == N || error("size of rank does not match")
#     if length(weight) != N
#        weight = Float32(1/N) * ones(Float32, N)
#     end
#     # normalize input data and initialization
#     etp    = eltype(dobs)
#     l2norm = vecnorm(dobs);
#     eptl2norm = sqrt(l2norm^2 * prod(I) / length(known));
#     for i = 1 : N
#         rand!(X[i]); scaler=eptl2norm^(I[i] /(I[i]+CI[i]))/vecnorm(X[i]); X[i]=X[i]*scaler;
#         rand!(Y[i]); scaler=eptl2norm^(CI[i]/(I[i]+CI[i]))/vecnorm(Y[i]); Y[i]=Y[i]*scaler;
#         transpose!(Xt[i],X[i]); conj!(Xt[i]);
#         transpose!(Yt[i],Y[i]); conj!(Yt[i]);
#     end
#     fill!(Z, zero(etp)); nobs = length(known);
#     for i = 1 : nobs
#         Z[known[i]] = dobs[i]
#     end
#
#     slope = (alpha-1.0f0)/(maxIter-1)
#     for iter = 1 : maxIter
#         # update X[n] and Y[n]
#         for n = 1 : N
#             Zn = matricization(Z,I,n); A_mul_B!(X[n],Zn,Yt[n])  ; transpose!(Xt[n],X[n]); conj!(Xt[n]); # update X[n]
#             Xsq= Xt[n]*X[n]          ; Y[n]=pinv(Xsq)*(Xt[n]*Zn); transpose!(Yt[n],Y[n]); conj!(Yt[n]); # update Y[n]
#         end
#         # update Z
#         fill!(Z, zero(etp))
#         for n = 1 : N
#             Zn = X[n]*Y[n];
#             Zn = unmatricization(Zn, I, n)
#             Z  = Z + weight[n]*Zn
#         end
#         if convergence==true
#            relerr   = vecnorm(Z[known]-dobs) / vecnorm(dobs)
#            println("relative error: $relerr")
#         end
#         # insert observations
#         beta = slope*(iter-1) + 1
#         for i = 1 : nobs
#             Z[known[i]] =(1-beta)*Z[known[i]] + beta*dobs[i]
#         end
#     end
#     return Z
# end



#  test the code on synthetic data set
# rank = [3,4,5,6]; I = [8,24,20,20]; N = length(I)
# CI = zeros(Int64, 4); etp = Float64;
# X = Vector{Matrix{etp}}(4)
# Y = Vector{Matrix{etp}}(4)
# weight = [0.25, 0.25, 0.25, 0.25]
# Z  = zeros(etp, I...)      # ground true result
# for i = 1 : N
#     CI[i] = round(Int64, prod(I)/I[i])
#     X[i]  = rand(etp, I[i], rank[i])
#     Y[i]  = rand(etp, rank[i], CI[i])
#     Z = Z + weight[i]*unmatricization(X[i]*Y[i], I, i)
# end
#
# total = prod(I); perc = 0.6; nobs = round(Int64, total*perc);
# knownIdx = rand(collect(1:total), nobs); knownIdx = sort(knownIdx);
# known = [knownIdx[1]]; tmp = knownIdx[1];
# for i = 2 : length(knownIdx)
#     if knownIdx[i] > tmp
#        tmp = knownIdx[i]
#        push!(known, tmp)
#     end
# end
# pmiss = 1 - length(known)/total
# print("missing rate is $pmiss")
# dobs = Z[known]  #observed data
# alpha = 1.0
# Zc = PMF(dobs, known, N, I, rank, alpha, groundTrue=Z)
