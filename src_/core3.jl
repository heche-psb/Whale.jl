#= Reimplementation of the core likelihood routine with goal to support
automatic differentiation.

The main problem is that during AD, all structures should allow <:Real
types and be unrestricted. This makes the current tmpmat approach
infeasible, since this is modifying an array in place, whcih has already
a determined Float type. =#

# Arthur Zwaenepoel - 2019
# some helper types
const PDict{T} = Dict{Int64,Array{T,1}} where T<:Real

Base.getindex(d::PDict{T}, e::Int64, i::Int64) where T<:Real = d[e][i]

Base.setindex!(d::PDict{T}, x::T, e::Int64, i::Int64) where T<:Real =
    d[e][i] = x

"""
    $(TYPEDEF)

The full Whale model, containing both the sliced species tree, parameters and
fields for extinction and propagation probabilities.
"""
struct WhaleModel{T<:Real,CCD} <: DiscreteUnivariateDistribution
    S::SlicedTree
    λ::Array{T}
    μ::Array{T}
    q::Array{T}
    η::T
    ε::PDict{T}
    ϕ::PDict{T}
    cond::String

    function WhaleModel(S::SlicedTree, λ::Array{T}, μ::Array{T},
            q::Array{T}=T[], η::T=0.9, cond="oib") where {T<:Real}
        ε = get_ε(S, λ, μ, q, η)
        ϕ = get_ϕ(S, λ, μ, q, η, ε)
        new{T,CCD}(S, λ, μ, q, η, ε, ϕ, cond)
    end
end

eltype(w::WhaleModel{_,T}) where {_,T} = T

function Base.show(io::IO, w::WhaleModel)
    write(io, "Tree of $(ntaxa(w.S)) taxa with $(nrates(w.S)) rate classes")
    write(io, " and $(nwgd(w.S)) WGDs")
end

function get_ε(s::SlicedTree, λ, μ, q, η)
    ε = PDict{typeof(η)}(e => zeros(nslices(s, e)) for e in s.border)
    for e in s.border
        if isleaf(s.tree, e)
            ε[e][1] = 0.
        elseif haskey(s.qindex, e)
            qe = q[s[e, :q]]
            f = childnodes(s.tree, e)[1]
            ε_wgd = ε[f][nslices(s, f)]
            ε[e, 1] = qe * ε_wgd^2 + (1-qe) * ε_wgd
        else
            f, g = childnodes(s.tree, e)
            ε[e, 1] = ε[f, nslices(s, f)] * ε[g, nslices(s, g)]
        end
        if isroot(s.tree, e)
            return ε
        end
        for i in 2:nslices(s, e)
            λe = λ[s[e, :λ]]
            μe = μ[s[e, :μ]]
            ε[e, i] = ε_slice(λe, μe, s[e, i], ε[e, i-1])
        end
    end
end

ε_slice(λ, μ, t, ε) = isapprox(λ, μ, atol=1e-5) ?
    1. + (1. - ε)/(μ * (ε - 1.) * t - 1.) :
        (μ + (λ - μ)/(1. + exp((λ - μ)*t)*λ*(ε - 1.)/(μ - λ*ε)))/λ

function get_ϕ(s::SlicedTree, λ, μ, q, η, ε::PDict)
    ϕ = PDict{typeof(η)}(e => zeros(nslices(s, e)) for e in s.border)
    for e in s.border[1:end-1]
        λe = λ[s[e, :λ]]
        μe = μ[s[e, :μ]]
        ϕ[e][1] = 1.
        for i in 2:nslices(s, e)
            ϕ[e, i] = ϕ_slice(λe, μe, s[e, i], ε[e, i-1])
        end
    end
    return ϕ
end

function ϕ_slice(λ, μ, t, ε)
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

check_args(m::WhaleModel) =
    all([m.λ; m.μ] .> 0) && all(1. .>= [[m.η]; m.q] .>= 0)

function logpdf(m::WhaleModel, x::CCD, node::Int64=-1; matrix=true)
    if ~check_args(m) ; return -Inf; end
    if x.Γ == -1
        return 0.
    else
        l, M = whale!(x, m.S, m.λ, m.μ, m.q, m.η, m.ε, m.ϕ, m.cond, node::Int64)
        return matrix ? (l, M) : l
    end
end

# main whale algorithm
function whale!(x::CCD, s::SlicedTree, λ::Array{T}, μ::Array{T}, q::Array{T},
        η::T, ε::PDict{T}, ϕ::PDict{T}, cond::String, node=-1) where T<:Real
    branches = (node == -1) ? s.border : get_parentbranches(s, node)
    M = DPMat{typeof(η)}()
    init_matrix!(M, x, s, branches)

    for e in s.border[1:end-1]  # skip the root branch
        qnode = haskey(s.qindex, e)
        sleaf = isleaf(s, e)
        λe = λ[s[e, :λ]]
        μe = μ[s[e, :μ]]

        for γ in x.clades
            !(x.species[γ] ⊆ s.clades[e]) ? (continue) : nothing
            γleaf = isleaf(x, γ)

            for i in 1:nslices(s, e)
                # speciation or WGD node
                if i == 1
                    if γleaf && x.m3[γ] == e
                        M[e][γ, i] = 1.0
                    elseif !(sleaf || qnode)
                        f, g = childnodes(s.tree, e)
                        if !γleaf
                            Π_speciation!(M, x, e, γ, f, g)
                        end
                        Π_loss!(M, x, e, γ, ε, f, g)
                    elseif qnode
                        qe = q[s[e, :q]]
                        f = childnodes(s.tree, e)[1]
                        if !γleaf
                            Π_wgd_retention!(M, x, e, γ, qe, f)
                        end
                        Π_wgd_non_retention!(M, x, e, γ, qe, f)
                        Π_wgd_loss!(M, x, e, γ, qe, ε[f][end], f)
                    end

                # internal of branch (slice boundary)
                else
                    Δt = s[e, i]
                    M[e][γ, i] += ϕ[e, i] * M[e][γ, i-1]
                    if !γleaf
                        Π_duplication!(M, x, e, i, γ, Δt, λe, μe)
                    end
                end
            end
        end
    end
    whale_root!(M, x, s, ε, η)
    l = cond_lhood(M, x, s, ε, η, cond)
    return l, M
end

function cond_lhood(M, x::CCD, s::SlicedTree, ε, η, cond::String)
    root = s.border[end]
    f, g = childnodes(s.tree, root)
    #if cond == "oib"
    return oib(M, x, root, f, g, ε, η)
    #end
end

function oib(M, x::CCD, e::Int64, f::Int64, g::Int64, ε, η)
    ε_root = geometric_extinctionp(ε[e, 1], η)
    ε_left = geometric_extinctionp(ε[f][end], η)
    ε_rght = geometric_extinctionp(ε[g][end], η)
    nf = 1 - ε_left - ε_rght + ε_root
    M[e, x.Γ, 1] > 0. && nf > 0. ? log(M[e, x.Γ, 1]/nf) : -Inf
end

geometric_extinctionp(ε, η) = η * ε / (1 - (1 - η)*ε)

function whale_root!(M, x::CCD, s::SlicedTree, ε, η)
    root = s.border[end]
    f, g = childnodes(s.tree, root)
    ε0 = ε[root, 1]
    η_ = 1.0/(1. - (1. - η) * ε0)^2
    for γ in x.clades
        γleaf = isleaf(x, γ)
        M[root][γ, 1] = 0.
        Π_loss!(M, x, root, γ, ε, f, g)
        if ~γleaf
            Π_speciation!(M, x, root, γ, f, g)
        end
        M[root][γ, 1] *= η_
        if ~γleaf
            Π_root!(M, x, root, γ, η, ε0)
        end
    end
    M[root, x.Γ, 1] *= η
end

function Π_root!(M, x::CCD, root::Int64, γ::Int64, η, ε0)
    p = 0.
    for (γ1, γ2, count) in x.m2[γ]
        p += x.ccp[(γ, γ1, γ2)] * M[root, γ1, 1] * M[root, γ2, 1]
    end
    p *= (1. - η) * (1. - (1. - η) * ε0)
    M[root, γ, 1] += p
end

function Π_speciation!(M, x::CCD, e::Int64, γ::Int64, f::Int64, g::Int64)
    p = 0.
    for (γ1, γ2, count) in x.m2[γ]
        p += x.ccp[(γ, γ1, γ2)] * M[f][γ1, end] * M[g][γ2, end]
        p += x.ccp[(γ, γ1, γ2)] * M[g][γ1, end] * M[f][γ2, end]
    end
    M[e, γ, 1] += p
end

function Π_loss!(M, x::CCD, e::Int64, γ::Int64, ε, f::Int64, g::Int64)
    M[e, γ, 1] += M[f][γ, end] * ε[g][end] + M[g][γ, end] * ε[f][end]
end

function Π_wgd_retention!(M, x::CCD, e::Int64, γ::Int64, q, f::Int64)
    p = 0.
    for (γ1, γ2, count) in x.m2[γ]
        p += x.ccp[(γ, γ1, γ2)] * M[f][γ1,end] * M[f][γ2,end]
    end
    M[e, γ, 1] += p * q
end

function Π_wgd_non_retention!(M, x::CCD, e::Int64, γ::Int64, q, f::Int64)
    M[e, γ, 1] += (1-q) * M[f][γ, end]
end

function Π_wgd_loss!(M, x::CCD, e::Int64, γ::Int64, q, ε, f::Int64)
    M[e, γ, 1] += 2 * q * ε * M[f][γ, end]
end

function Π_duplication!(M, x::CCD, e::Int64, i::Int64, γ::Int64, Δt, λe, μe)
    p = 0.
    for (γ1, γ2, count) in x.m2[γ]
        p += x.ccp[(γ, γ1, γ2)] * M[e, γ1, i-1] * M[e, γ2, i-1]
    end
    # p12 = p_transition_kendall(2, Δt, λ, μ)
    # return p * p12
    #= NOTE: it would be more correct to use the BD transition probability to go
    from one lineage to two lineages. In practice it doesn't make a
    difference as long as the slice lengths (Δt) are  short enough (∼ [λ/10,
    λ/5]) and using the approximation seems to counteract some bias. =#
    M[e, γ, i] += p * λe * Δt
end

function init_matrix!(M::DPMat{T}, x::CCD, s::SlicedTree, bs) where T<:Real
    for b in s.border
        if b in bs
            M[b] = zeros(T, length(x.clades), length(s.slices[b]))
        else
            M[b] = x.tmpmat[b]
        end
    end
end
