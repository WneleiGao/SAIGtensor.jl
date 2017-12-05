export tensor,
       initTensor,
       matricization,
       unmatricization,
       KRP,
       recursiveKRP,
       tnorm,
       ttv,
       ttm,
       patch,
       getPatchDir,
       getOnePatch,
       WrapTaper,
       Taper,
       UnPatch,
       PatchProcess,
       SeisMixEvents3D

# include("basicOperations.jl")
include("tensor.jl")
include("patchAndUnpatch.jl")
include("SeisMixEvents.jl")
