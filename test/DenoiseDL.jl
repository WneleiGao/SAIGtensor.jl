# tests on one cube
using PyPlot, Seismic, SAIGtensor

(d, e) = SeisParabEvents(nx2=100, nt=500, p2=[0.4, 0.4, 0.4], amp=[1.0, -1.0, 0.7]);
snr = 1.0;

d = dc[201:400, 1:200, 1:200]
ds = SeisAddNoise(d, snr, L=5);

X1 = reshape(ds, 200, 200*200)
X2 = permutedims(ds, (2,1,3));
X2 = reshape(X2, 200, 200*200)
X3 = permutedims(ds, (3,2,1));
X3 = reshape(X3, 200, 200*200);

K = 400; L=10;
maxIter = 20;
(D, coef) = AKSVD(X1, size(X1,1), L, maxIter)

dr = zeros(X1)
for i = 1 : length(coef)
    dr[:,i] = D * coef[i]
end
dr = reshape(dr, 200, 200, 200);
X1 = reshape(X1, 200, 200, 200);

i = 4;
SeisPlot(dr[:,:,i]); SeisPlot(X[:,:,i])
SeisPlot(dr[:,:,i], style="wiggles", wiggle_trace_increment=2);
SeisPlot( X[:,:,i], style="wiggles", wiggle_trace_increment=2);


(D2, coef2) = AKSVD(X2, size(X2,1), L, maxIter)
dr2 = similar(X2)
for i = 1: length(coef2)
    dr2[:,i] = D2 * coef2[i]
end
dr2 = reshape(dr2, 100, nt, 100);
dr2 = permutedims(dr2, (2,1,3));


(D3, coef3) = AKSVD(X3, size(X3,1), L, maxIter)
dr3 = similar(X3)
for i = 1: length(coef3)
    dr3[:,i] = D3 * coef3[i]
end
dr3 = reshape(dr3, 100, 100, nt);
dr3 = permutedims(dr3, (3,2,1));

result = 1/3 .* (dr .+ dr2 .+ dr3)
i = 1;
SeisPlot(result[:,:,i]); SeisPlot(X[:,:,i]);
SeisPlot(result[:,:,i], style="wiggles", wiggle_trace_increment=2);
SeisPlot(X[:,:,i], style="wiggles", wiggle_trace_increment=2);
tmp = hcat(X[:,:,i], result[:,:,i], result[:,:,i]-X[:,:,i]);
SeisPlot(tmp, style="wiggles", wiggle_trace_increment=4);

i = 2
tmp = hcat(dc[:,:,i], X[:,:,i], result[:,:,i], result[:,:,i]-X[:,:,i], result[:,:,i]-dc[:,:,i]);
SeisPlot(tmp, style="wiggles", wiggle_trace_increment=5);
dc = copy(d[istart:istart+nt-1,:,:]);

CalcSNRdB(dc, result)



# function DenoiseWithDL{Tv<:AbstractFloat, Ti<:Int64}(Img::Matrix{Tv}, sigma::Tv, K::Ti,
#                        maxIter::Ti; blockSize=8, maxBlockNum=900000, maxBlocksToTrain = 65000,
#                        slideRate=1, reduceDC=true)
#
#     etp = eltype(Img)
#     removeMean = 1
#     (nr, nc) = size(Img)
#     numIterKSVD = sigma > 5.0 ? 10 : 5
#     C = 1.15
#
#     numbr = nr - blockSize + 1
#     numbc = nc - blockSize + 1
#     totalBlocks = numbr * numbc
#
#     # re-organize blocks of image into training matrix
#     if totalBlocks > maxBlocksToTrain
#        X = zeros(etp, blockSize^2, maxBlocksToTrain)
#        randSelect = rand(1:totalBlocks, maxBlocksToTrain)
#        for i = 1 : maxBlocksToTrain
#            blockId = randSelect[i]
#            (ridx, cidx) = ind2sub((numbr, numbc), blockId)
#            rl = ridx; ru = rl + blockSize - 1;
#            cl = cidx; cu = cl + blockSize - 1;
#            X[:, i] = vec(view(Img, rl:ru, cl:cu))
#        end
#     else
#        X = zeros(etp, blockSize^2, totalBlocks)
#        for cidx = 1 : numbc
#            for ridx = 1 : numbr
#                rl = ridx; ru = rl + blockSize - 1;
#                cl = cidx; cu = cl + blockSize - 1;
#                blockId = (cidx-1)*numbr + ridx
#                X[:, blockId] = vec(view(Img, rl:ru, cl:cu))
#            end
#        end
#     end
#
#     # the error stoping criteria for OMP
#     tol = sigma * C
#     # initialize Dictionary as over-complete a cosine basis
#     sqrtK = ceil(Int64, sqrt(K))
#     dct   = zeros(blockSize, sqrtK)
#     increment = pi / sqrtK
#     for j = 1 : sqrtK
#         for i = 1 : blockSize
#             dct[i,j] = cos((i-1)*(j-1)*increment)
#         end
#         if j > 1
#            dct[:,j] = dct[:,j] - mean(dct[:,j])
#         end
#     end
#     D = kron(dct, dct)[:,1:K]
#
#     # remove the direct component of block matrix
#     if removeDC
#        vecOfMeans = mean(X,1)
#        for i = 1 : size(X,2)
#            for j = 1 : blockSize
#                X[j,i] = X[j,i] - vecOfMeans[i]
#            end
#        end
#     end
#
#     (learnedD, data) = AKSVD(X, K, L, maxIter, InitialMethod="GievnMatrix", GivenMatrix=D)
# end

# function PMF{Tv<:AbstractFloat, Ti<:Integer}(X::Array{Tv,4}, rk::Vector{Ti}, alpha::Vector{Tv}; maxIter=100)
#     length(rk) == 4 || error("number of rank is wrong")
#     length(alpha) == 4 || error("number of weight is wrong")
#     (n1, n2, n3, n4) = size(X)
#     X1 = copy(reshape(X, n1, n2*n3*n4))
#     X2 = permutedims(X, [2,1,3,4]); X2 = reshape(X2, n2, n1*n3*n4);
#     X3 = permutedims(X, [3,2,1,4]); X3 = reshape(X3, n3, n2*n1*n4);
#     X4 = permutedims(X, [4,2,3,1]); X4 = reshape(X4, n4, n2*n3*n1);
#     for iter = 1 : maxIter
#
#     end
#
#
# end
