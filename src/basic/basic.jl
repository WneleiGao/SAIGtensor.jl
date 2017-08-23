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
       Taper,
       UnPatch,
       PatchProcess,
       WrapRandCpAls,
       WrapTaper

# include("basicOperations.jl")
include("tensor.jl")
include("patchAndUnpatch.jl")
