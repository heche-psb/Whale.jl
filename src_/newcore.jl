#= It should be possible to use Mamba as MCMC engine, since we are actually not
doing a very special MCMC algorithm like e.g. in tree inference (where the
tree structure requires a lot of non-standard stuff in the algorithms). In the
(WH)ALE case however, the MCMC itself is more amenable to implementation in
external MCMC engines.

What we need would be the following:

- A struct analogous to a UnivariateDistribution but where the variate is a
CCD. This 'distribution' object should have a `logpdf` method which returns
the likelihood of a CCD given S, λ, μ, q, ... (which are the 'parameters' of
the distribution)

- Similar structs for the rate priors (GBM, IID).

That will be some work to implement but it's a nice opportunity to have alook
at the core code again and we might get a more efficient/reliable/maintainable
result in the end.

More tricky is how to recompute all those speedups like the partial recompute
etc.? By storing the parameter values associated with the last computation in
the CCD object? Perhaps a separate CCDArray type would work for thoe purposes

XXX PROBLEM: when Mamba is exectuted with julia in parallel, it runs chains in
parallel, but I want the likelihood to be computed in parallel.
=#

# XXX: OK here I start some rewrites with new structs etc. to make the code
# generally better and allow MCMC with Mamba (at least in principle).

#= Mamba sketch; I'm starting to like Turing.jl better...
model = Model(
    ν = Logical(() -> 0.1, false)
    η = Stochastic(() -> Beta(4, 2))
    θ = Stochastic(1,
        () -> MvNormal([log(0.2), log(0.2)], [0.5 0.25; 0.25 0.5]), false)
    λ0 = Logical((θ)->exp(θ[1]), false)
    μ0 = Logical((θ)->exp(θ[2]), false)
    λ = Stochastic(1, (λ0, ν, tree) -> GeometricBrownianMotion(tree, λ0, ν))
    μ = Stochastic(1, (μ0, ν, tree) -> GeometricBrownianMotion(tree, μ0, ν))
    q = Stochastic(1, (tree) -> UnivariateDistribution[
        Beta(1,1) for i=1:length(tree.qindex)])
    X = Stochastic(1, (λ, μ, q, η, tree) -> WhaleSamplingDist(λ, μ, q, η, tree,
        "oib"))
)
=#

#= Turing.jl
@model Whale(X, S) = begin
    # X is the observed data (a DArray of CCDs?)
    # S is the sliced species tree
    ν = 0.1
    η ~ Beta
    θ ~ MultivariateNormal
    l = θ[1]
    m = θ[2]
    λ ~ GeometricBrownianMotion(l, ν, S)
    μ ~ GeometricBrownianMotion(m, ν, S)
    q = tzeros(length(S.qindex))
    for i=1:length(S.qindex)
        q[i] ~ Beta(1., 1.)
    end
    X ~ WhaleSamplingDist(λ, μ, q, η, S, "oib")
end =#

# Every CCD object stores it's own rec. matrix
# Where should the extinction and propagation probabilities be stored? I guess
# at the WhaleModel level, because they result from applying the
# WhaleSamplingDist to the tree?

const PDict = Dict{Int64,Array{Float64,1}}
const DPMat = Dict{Int64,Array{Float64,2}}

mutable struct WhaleParams{T<:Real}
    λ::Array{T}
    μ::Array{T}
    q::Array{T}
    η::T

    function WhaleParams(λ::Array{T}, μ::Array{T}, q::Array{T}, η::T,) where T
        @check_args(WhaleParams, all([λ; μ] .> 0) &&
            all(1. .>= [[η]; q] .>= 0))
        return new{T}(λ, μ, q, η)
    end
end

struct WhaleModel
    S::SlicedTree
    M::WhaleParams
    ε::PDict
    ϕ::PDict
    cond::String

    function WhaleModel(S::SlicedTree, M::WhaleParams, cond="oib")
        ε = get_ε(S, M)
        ϕ = get_ϕ(S, M, ε)
        new(S, M, ε, ϕ, cond)
    end
end

# this should automatically result in partial recomputation when applicable...
# somewhere (maybe in the SlicedTree) there should be a `lastnode` field or so
# or if the update is in a fixed order, we might do it more clever?
logpdf(x::CCD, m::WhaleModel, node::Int64=-1) =
    x.Γ == -1 ? 0. : whale!(x, m.S, m.M, m.ε, m.ϕ, m.cond, node::Int64)

# main whale algorithm
function whale!(x::CCD, s::SlicedTree, m::WhaleParams, ε, ϕ, cond, node)
    branches = (node == -1) ? s.border[1:end-1] : get_branches(s, node)
    set_tmpmat!(x, s, branches)

    for e in branches  # skip the root branch
        qnode = haskey(s.qindex, e)
        sleaf = isleaf(s.tree, e)
        λe = m.λ[s.rindex[e]]
        μe = m.μ[s.rindex[e]]

        for γ in x.clades
            !(x.species[γ] ⊆ s.clades[e]) ? (continue) : nothing
            γleaf = haskey(x.m3, γ)

            for i in 1:length(s.slices[e])
                # speciation or WGD node
                if i == 1
                    if γleaf && x.m3[γ] == e
                        x.tmpmat[e][γ, i] = 1.0
                    elseif !(sleaf || qnode)
                        f, g = childnodes(s.tree, e)
                        if !(γleaf)
                            Π_speciation!(x, e, γ, f, g)
                        end
                        Π_loss!(x, e, γ, ε, f, g)
                    elseif qnode
                        qe = m.q[s.qindex[e]]
                        f = childnodes(s.tree, e)[1]
                        if !(γleaf)
                            Π_wgd_retention!(x, e, γ, qe, f)
                        end
                        Π_wgd_non_retention!(x, e, γ, qe, f)
                        Π_wgd_loss!(x, e, γ, qe, ε[f][end], f)
                    end

                # internal of branch (slice boundary)
                else
                    Δt = s.slices[e][i]
                    x.tmpmat[e][γ, i] += ϕ[e][i] * x.tmpmat[e][γ, i-1]
                    if !(γleaf)
                        Π_duplication!(x, e, i, γ, Δt, λe, μe)
                    end
                end
            end
        end
    end

    whale_root!(x, s, ε, m.η)
    cond_lhood(x, s, ε, m.η, cond)
end

function cond_lhood(x::CCD, s::SlicedTree, ε::PDict, η::Float64, cond::String)
    root = s.border[end]
    f, g = childnodes(s.tree, root)
    if cond == "oib"
        return oib(x, root, f, g, ε, η)
    else
        @error "Not implemented"
    end
end

function oib(x::CCD, e::Int64, f::Int64, g::Int64, ε::PDict, η::Float64)
    ε_root = geometric_extinctionp(ε[e][1], η)
    ε_left = geometric_extinctionp(ε[f][end], η)
    ε_rght = geometric_extinctionp(ε[g][end], η)
    nf = 1 - ε_left - ε_rght + ε_root
    x.tmpmat[e][x.Γ, 1] > 0. && nf > 0. ? log(x.tmpmat[e][x.Γ,1]/nf) : -Inf
end

geometric_extinctionp(ε::Float64, η::Float64) = η * ε / (1 - (1 - η)*ε)

function whale_root!(x::CCD, s::SlicedTree, ε::PDict, η::Float64)
    root = s.border[end]
    f, g = childnodes(s.tree, root)
    ε0 = ε[root][1]
    η_ = 1.0/(1. - (1. - η) * ε0)^2
    for γ in x.clades
        γleaf = haskey(x.m3, γ)
        x.tmpmat[root][γ, 1] = 0.
        Π_loss!(x, root, γ, ε, f, g)
        if !(γleaf)
            Π_speciation!(x, root, γ, f, g)
        end
        x.tmpmat[root][γ, 1] *= η_
        if !(γleaf)
            Π_root!(x, root, γ, η, ε0)
        end
    end
    x.tmpmat[root][x.Γ, 1] *= η
end

function Π_root!(x::CCD, root::Int64, γ::Int64, η::Float64, ε0::Float64)
    p = 0.
    for (γ1, γ2, count) in x.m2[γ]
        p += x.ccp[(γ, γ1, γ2)] * x.tmpmat[root][γ1, 1] * x.tmpmat[root][γ2, 1]
    end
    p *= (1. - η) * (1. - (1. - η) * ε0)
    x.tmpmat[root][γ, 1] += p
end

function Π_speciation!(x::CCD, e::Int64, γ::Int64, f::Int64, g::Int64)
    p = 0.
    for (γ1, γ2, count) in x.m2[γ]
        p += x.ccp[(γ, γ1, γ2)] * x.tmpmat[f][γ1, end] * x.tmpmat[g][γ2, end]
        p += x.ccp[(γ, γ1, γ2)] * x.tmpmat[g][γ1, end] * x.tmpmat[f][γ2, end]
    end
    x.tmpmat[e][γ, 1] += p
end

function Π_loss!(x::CCD, e::Int64, γ::Int64, ε::PDict, f::Int64, g::Int64)
    x.tmpmat[e][γ, 1] += x.tmpmat[f][γ,end] * ε[g][end] +
        x.tmpmat[g][γ,end] * ε[f][end]
end

function Π_wgd_retention!(x::CCD, e::Int64, γ::Int64, q::Float64, f::Int64)
    p = 0.
    for (γ1, γ2, count) in x.m2[γ]
        p += x.ccp[(γ, γ1, γ2)] * x.tmpmat[f][γ1, end] * x.tmpmat[f][γ2, end]
    end
    x.tmpmat[e][γ, 1] += p * q
end

function Π_wgd_non_retention!(x::CCD, e::Int64, γ::Int64, q::Float64, f::Int64)
    x.tmpmat[e][γ, 1] += (1-q) * x.tmpmat[f][γ, end]
end

function Π_wgd_loss!(x::CCD, e::Int64,γ::Int64,q::Float64, ε::Float64, f::Int64)
    x.tmpmat[e][γ, 1] += 2 * q * ε * x.tmpmat[f][γ, end]
end

function Π_duplication!(x::CCD, e::Int64, i::Int64, γ::Int64, Δt::Float64,
        λe::Float64, μe::Float64)
    p = 0.
    for (γ1, γ2, count) in x.m2[γ]
        p += x.ccp[(γ, γ1, γ2)] * x.tmpmat[e][γ1, i-1] * x.tmpmat[e][γ2, i-1]
    end
    # p12 = p_transition_kendall(2, Δt, λ, μ)
    # return p * p12
    #= NOTE: it would be more correct to use the BD transition probability to go
    from one lineage to two lineages. In practice it doesn't make a
    difference as long as the slice lengths (Δt) are  short enough (∼ [λ/10,
    λ/5]) and using the approximation seems to counteract some bias. =#
    x.tmpmat[e][γ, i] += p * λe * Δt #* 2
end

function set_εϕ!(w::WhaleModel)
    w.ε = get_ε(w.S, w.M)
    w.ϕ = get_ϕ(w.S, w.M, w.ε)
end

function get_ε(s::SlicedTree, w::WhaleParams)
    ε = get_pdict(s)
    for e in s.border
        if isleaf(s.tree, e)
            ε[e][1] = 0.
        elseif haskey(s.qindex, e)
            qe = m.q[s.qindex[e]]
            f = childnodes(s.tree, e)[1]
            ε_wgd = ε[f][length(s.slices[f])]
            ε[e][1] = qe * ε_wgd^2 + (1-qe) * ε_wgd
        else
            f, g = childnodes(s.tree, e)
            ε[e][1] = ε[f][length(s.slices[f])] * ε[g][length(s.slices[g])]
        end
        if isroot(s.tree, e)
            return ε
        end
        for i in 2:length(s.slices[e])
            λe = m.λ[s.rindex[e]]
            μe = m.μ[s.rindex[e]]
            ε[e][i] = ε_slice(λe, μe, s.slices[e][i], ε[e][i-1])
        end
    end
end

ε_slice(λ::Float64, μ::Float64, t::Float64, ε::Float64) =
    isapprox(λ, μ, atol=1e-5) ? 1. + (1. - ε)/(μ * (ε - 1.) * t - 1.) :
        (μ + (λ - μ)/(1. + exp((λ - μ)*t)*λ*(ε - 1.)/(μ - λ*ε)))/λ

function get_ϕ(s::SlicedTree, m::WhaleParams, ε::PDict)
    ϕ = get_pdict(s)
    for e in s.border[1:end-1]
        λe = m.λ[s.rindex[e]]
        μe = m.μ[s.rindex[e]]
        ϕ[e][1] = 1.
        for i in 2:length(s.slices[e])
            ϕ[e][i] = ϕ_slice(λe, μe, s.slices[e][i], ε[e][i-1])
        end
    end
    return ϕ
end

function ϕ_slice(λ::Float64, μ::Float64, t::Float64, ε::Float64)
    if isapprox(λ, μ, atol=1e-5)
        return 1. / (μ * (ε - 1.) * t - 1.)^2
    else
        x = exp((μ - λ)*t)
        a = x * (λ - μ)^2
        b = λ - (x * μ)
        c = (x - 1.) * λ * ε
        return a / (b + c)^2
    end
end

get_pdict(s::SlicedTree) = PDict(k=>zero(v) for (k, v) in s.slices)

function set_tmpmat!(x::CCD, s::SlicedTree, branches::Array{Int64})
    for b in [branches; s.border[end]]
        x.tmpmat[b] = zeros(length(x.clades), length(s.slices[b]))
    end
end
