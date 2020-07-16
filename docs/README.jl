# [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://arzwa.github.io/Whale.jl/dev/index.html)
# [![Build Status](https://travis-ci.com/arzwa/Whale.jl.svg?branch=master)](https://travis-ci.com/arzwa/Whale.jl)

# # Whale: Bayesian gene tree reconciliation and whole-genome duplication inference by amalgamated likelihood estimation

#```
#
#                             .-------------'```'----....,,__                        _,
#                            |                               `'`'`'`'-.,.__        .'(
#                            |                                             `'--._.'   )
#                            |                                                   `'-.<
#                            \               .-'`'-.                            -.    `\
#                             \               -.o_.     _                     _,-'`\    |
#                              ``````''--.._.-=-._    .'  \            _,,--'`      `-._(
#                                (^^^^^^^^`___    '-. |    \  __,,..--'                 `
#                                 `````````   `'--..___\    |`
#                                                       `-.,'
#```

# Whale.jl is a julia library implementing joint inference of gene tree topologies and their reconciliations to a species tree using the **amalgamation** method of Szollosi et al. (2014) to compute the marginalize the reconciliation likelihood over a distribution over tree topologies. Whale implements the duplication-loss (DL) model of gene family evolution as well as a duplication-loss and whole-genome duplication (DLWGD) model (Rabier et al. 2014, Zwaenepoel et al. 2019). The latter can be used for the inference of ancient whole-genome duplications (WGDs) from gene trees while taking into account gene tree and reconciliation uncertainty.

# The likelihood routines implemented in Whale support **automatic differentiation** using `ForwardDiff.jl`, allowing for efficient gradient-based Maximum-likelihood estimation and Hamiltonian Monte Carlo (HMC) based Bayesian inference. The library focuses on the Bayesian case, and implements relaxed clock priors to model the evolution of gene duplication and loss rates. Lastly, Whale allows to sample reconciled trees from the posterior distribution or a parameterized DL(+WGD) model using a stochastic backtracking agorithm (as in [ALE](https://github.com/ssolo/ALE)).

# Please have a look at the [docs](https://arzwa.github.io/Whale.jl/dev/index.html) for usage instructions and documentation. You might want to get some minimal familiarity with the Julia REPL and its package manager when using Whale, see [the julia docs](https://docs.julialang.org/en/v1/).

# Note that the scripts in the `scripts` directory might be helpful to prepare data for Whale analyses.

# ## Quickstart using Turing and a constant-rates model
using Whale, NewickTree, Distributions, Turing, DataFrames

# Get the tree
t = deepcopy(Whale.extree)
n = length(postwalk(t))  # number of internal nodes

# Now we add two WGD nodes to the tree. We do this by specifying
# the last common ancestor node for the lineages that share the
# hypothetical WGD. By default, the added node is halfway between
# the specified node and its parent.
insertnode!(getlca(t, "ATHA", "ATHA"), name="wgd")
insertnode!(getlca(t, "ATHA", "ATRI"), name="wgd")

# and we obtain a reference model object, here we will use a constant-rates
# model
θ = ConstantDLWGD(λ=0.1, μ=0.2, q=[0.2, 0.1], η=0.9)
r = Whale.RatesModel(θ, fixed=(:p,))
w = WhaleModel(r, t, .1)

# next we get the data (we need a model object for that)
ccd = read_ale(joinpath("example/example-1/ale"), w)

# Now we define the Turing model
@model constantrates(model, ccd) = begin
    r  ~ MvLogNormal(ones(2))
    η  ~ Beta(3,1)
    q1 ~ Beta()
    q2 ~ Beta()
    ccd ~ model((λ=r[1], μ=r[2], η=η, q=[q1, q2]))
end

model = constantrates(w, ccd)
chain = sample(model, NUTS(0.65), 100)
pdf = DataFrame(chain)

# We can sample reconciled trees from the posterior using a backtracking algorithm
fun = (m, x)-> Array(x) |> x->m((λ=x[3], μ=x[4], η=x[5], q=x[1:2]))
tt = TreeTracker(w, ccd[end-1:end], pdf, fun)
trees = track(tt)

# ## Citation

# If you use Whale, please cite:

# >[Zwaenepoel, A. and Van de Peer, Y., 2019. Inference of Ancient Whole-Genome Duplications and the Evolution of Gene Duplication and Loss Rates. *Molecular biology and evolution*, 36(7), pp.1384-1404.](https://academic.oup.com/mbe/article-abstract/36/7/1384/5475503)

using Literate  #src
Literate.markdown(joinpath(@__DIR__, "README.jl"), joinpath(@__DIR__, "../"), documenter=false, execute=true)  #src
