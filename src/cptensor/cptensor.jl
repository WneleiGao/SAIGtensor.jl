export cptensor,
       initCptensor,
       cpnorm,
       innerprod,
       fitness,
       updateAn,
       cptfirst,
       cptRight,
       cptLeft,
       cp_gradient,
       cpals,
       cpals_fit,
       cpals_time,
       ttf,
       ttf_slow,
       cpnormalize!,
       cp2tensor,
       sidx2sub,
       randSpl!,
       formXst!,
       formZs!,
       updateAn!,
       randCpAls,
       WrapCpAls,
       WrapRandCpAls,
       randCpAls_fit,
       randCpAls_time

include("conventionalCP.jl")
include("randCpAls.jl")
