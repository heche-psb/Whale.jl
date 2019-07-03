module Whale
    using Distributed
    using PhyloTrees
    using Distributions
    using DistributedArrays
    using Optim
    using ForwardDiff
    using Random
    using MCMCChains
    using DataFrames
    import ProgressMeter: @showprogress
    import Distributions: @check_args, logpdf

    # I guess order matters for type dependencies
    include("slicedtree.jl")
    include("ccd.jl")
    include("core.jl")
    include("mle.jl")
    include("gbm.jl")
    include("mcmc.jl")

    export
        SlicedTree, WhaleModel, read_ale, logpdf, nwgd, nrates, nslices, ntaxa,
        CCD, gradient, asvector1, mle, set_constantrates!, describe, WhaleChain,
        GBMModel, IRModel, mcmc!
end
