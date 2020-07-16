var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Index",
    "title": "Index",
    "category": "page",
    "text": ""
},

{
    "location": "#Introduction-1",
    "page": "Index",
    "title": "Introduction",
    "category": "section",
    "text": "warning: Warning\nThe latest Whale version is a thorough rewrite of the Whale library, and is still work in progress. For the version as used in Zwaenepoel & Van de Peer (2019), refer to this release (v0.2). Nevertheless, the current version should be safe to use and is more efficient and convenient (if you know a bit of julia).Whale provides tools for (genome-wide) amalgamated likelihood estimation (ALE) under a DL+WGD model, which is an approach to infer reconciled gene trees and parameters of a model of gene family evolution given a known species tree.The ALE approach takes into account uncertainty in the gene tree topology by marginalizing over all tree topologies that can be amalgamated from the so-called conditional clade distribution (CCD). This CCD can be constructed from a sample of the posterior distribution of tree topologies (which can be obtained using any standard software for Bayesian phylogenetics).More specifically, this library can be used toStatistically test the absence or presence of hypothetical whole-genome duplication (WGD) events in a species phylogeny\nInfer lineage-specific gene duplication and loss rates for a species phylogeny\nInfer high-quality (reconciled) gene trees given a known species tree cf. Szöllősi et al.\nAll of the above at oncenote: Note\nThis library implements the DL and DL+WGD models. It does not implement models of gene family evolution that take into account horizontal gene transfer, incomplete lineage sorting or gene conversion."
},

{
    "location": "#Installation-1",
    "page": "Index",
    "title": "Installation",
    "category": "section",
    "text": "You will need julia-1.x. Fire up a julia REPL by typing julia at the command line enter the package manager interface by typing ], and execute add Whale."
},

{
    "location": "#Data-preparation-1",
    "page": "Index",
    "title": "Data preparation",
    "category": "section",
    "text": "To perform analyses with Whale, you will need  An ultrametric species tree, with ideally branch lengths in geological time (since this allows straightforward interpretation of parameter estimates.)\nA bunch of ALE files, which summarize the conditional clade distributions (CCDs) for the same bunch of gene families. These can be obtained from a sample of the posterior distribution of gene trees using the ALEobserve tool. A pipeline to obtain these from a set of gene family protein fasta files is available at github.note: Note\nGene IDs should be prefixed by the name of the species to which the gene belongs as used in the species tree. For example if Arabidopsis thaliana is represented by ATHA in the species tree newick file, then the genes should be prefixed with ATHA_, e.g. ATHA_AT1G05000.note: Note\nAnalyzing CCDs (ALE files) with a very large number of clades or for very large families can be prohibitive computationally. It is therefore generally advisable that large orthogroups are filtered out based on some criterion (for example using the script orthofilter.py in the scripts directory of the Whale repository). To filter out families with very large numbers of clades in the CCD (which reflects that there is a lot of uncertainty in the gene tree), the scripts ccddata.py and ccdfilter.py can be used. This is a rather ad hoc filtering procedure, but can be useful to filter out families that trouble the analysis.warning: Warning\nMost analyses in Whale assume that for each family, there is at least one gene in both clades stemming from the root of the species tree. The likelihood in Whale is the conditional likelihood under this assumption. This is to rule out the possibility of de novo gain of a gene family along a branch of the species tree. The orthogroup data should therefore always be filtered to be in accordance with this criterion. This can also be done using the orthofilter.py script."
},

{
    "location": "#References-1",
    "page": "Index",
    "title": "References",
    "category": "section",
    "text": "Whale.jl is developed by Arthur Zwaenepoel at the VIB-UGent center for plant systems biology (bioinformatics & evolutionary genomics group). If you use Whale, please cite:Zwaenepoel, A. and Van de Peer, Y., 2019. Inference of Ancient Whole-Genome Duplications and the Evolution of Gene Duplication and Loss Rates. Molecular biology and evolution, 36(7), pp.1384-1404.The methods in Whale are heavily inspired by previous work done by other researchers. If you use Whale, consider citing the following two particularly important studies:[ALE] Szöllősi, G.J., Rosikiewicz, W., Boussau, B., Tannier, E. and Daubin, V., 2013. Efficient exploration of the space of reconciled gene trees. Systematic biology, 62(6), pp.901-912.[DL+WGD model] Rabier, C.E., Ta, T. and Ané, C., 2013. Detecting and locating whole genome duplications on a phylogeny: a probabilistic approach. Molecular biology and evolution, 31(3), pp.750-762."
},

{
    "location": "wgd-dhmc/#",
    "page": "Bayesian inference using NUTS with DynamicHMC.jl",
    "title": "Bayesian inference using NUTS with DynamicHMC.jl",
    "category": "page",
    "text": ""
},

{
    "location": "wgd-dhmc/#Bayesian-inference-using-NUTS-with-DynamicHMC.jl-1",
    "page": "Bayesian inference using NUTS with DynamicHMC.jl",
    "title": "Bayesian inference using NUTS with DynamicHMC.jl",
    "category": "section",
    "text": "using Whale, DynamicHMC, Random, NewickTree, Distributions, DataFrames\nusing DynamicHMC.Diagnostics\nRandom.seed!(562)Set up the model and the data, here I use a model with constant duplication and loss rates across the species tree. Note that the tree contains two WGD events.tree  = readnw(\"((MPOL:4.752,(PPAT:2.752)wgd_1:2.0):0.292,(SMOE:4.457,((((OSAT:1.555,(ATHA:0.5548,CPAP:0.5548):1.0002):0.738,ATRI:2.293):1.0)wgd_2:0.225,(GBIL:3.178,PABI:3.178):0.34):0.939):0.587);\")\nn = length(postwalk(tree))\nntaxa = (n+1)÷2\nrates = RatesModel(ConstantDLWGD(λ=0.1, μ=0.1, q=[0.2, 0.3], η=0.9))\nmodel = WhaleModel(rates, tree, 0.1)\ndata  = read_ale(joinpath(@__DIR__, \"../../example/example-1/ale\"), model, true)\nprior = Whale.CRPrior()\nproblem = WhaleProblem(data, model, prior)Run NUTS using DynamicHMC, (of course this is a ridicuously short run, and it\'s better to keep doubling_stages >= 3)results = mcmc_with_warmup(Random.GLOBAL_RNG, problem, 100,\n    warmup_stages=DynamicHMC.default_warmup_stages(doubling_stages=2))\nsummarize_tree_statistics(results.tree_statistics)Obtain the posterior distributionposterior = Whale.transform(problem, results.chain)\ndf = Whale.unpack(posterior)\ndescribe(df)Obtain reconciled trees sampled from the posteriortrees = track(problem, posterior)Consider the first gene familyfamily1 = trees[1].treesget the MAP tree as a newick stringnwstr(family1[1].tree)The support values are posterior probabilities for the associated reconciled split. Note that the tree does not contain branch lengths.The events field for each gene family contains a summary of the expected number of events for each branchtrees[1].eventsWe can get for every gene pair the posterior reconciliation probabilitypair_pps = Whale.getpairs(trees, model)\nfirst(pair_pps, 5)which sum to one, as we can verifymap(sum, eachrow(pair_pps[!,1:end-2]))Take for instance the following gene pair (second row)x = pair_pps[2,:]\nfor (n, v) in zip(names(x), Array(x))\n    (!(typeof(v)<:Number) || v > 0.) && println(n, \": \", v)\nendThe posterior probability (under the DL model) that this gene pair traces back to the speciation corresponding to node 17 (i.e. the root) is approximately 0.86, whereas the posterior probability that this gene pair traces back to an acnestral duplication event is 0.14.We can get a WGD-centric view as well. The following retrieves a table for each WGD with all gene tree nodes that have a non-zero psoterior probability of being reconciled to that particular WGD nodetables = Whale.getwgdtables(trees, data, model)\ntables[1]This page was generated using Literate.jl."
},

{
    "location": "wgd-turing/#",
    "page": "Bayesian inference using Turing.jl",
    "title": "Bayesian inference using Turing.jl",
    "category": "page",
    "text": ""
},

{
    "location": "wgd-turing/#Bayesian-inference-using-Turing.jl-1",
    "page": "Bayesian inference using Turing.jl",
    "title": "Bayesian inference using Turing.jl",
    "category": "section",
    "text": "using Whale, NewickTree, Distributions, Turing, DataFrames, LinearAlgebra"
},

{
    "location": "wgd-turing/#Using-a-constant-rates-model-1",
    "page": "Bayesian inference using Turing.jl",
    "title": "Using a constant-rates model",
    "category": "section",
    "text": "Get the treet = deepcopy(Whale.extree)\nn = length(postwalk(t))  # number of internal nodesNow we add two WGD nodes to the tree. We do this by specifying the last common ancestor node for the lineages that share the hypothetical WGD. By default, the added node is halfway between the specified node and its parent.insertnode!(getlca(t, \"PPAT\", \"PPAT\"), name=\"wgd_1\")\ninsertnode!(getlca(t, \"ATHA\", \"ATRI\"), name=\"wgd_2\")and we obtain a reference model object, here we will use a constant-rates modelθ = ConstantDLWGD(λ=0.1, μ=0.2, q=[0.2, 0.1], η=0.9)\nr = Whale.RatesModel(θ, fixed=(:p,))\nw = WhaleModel(r, t, .1)next we get the data (we need a model object for that)ccd = read_ale(joinpath(@__DIR__, \"../../example/example-1/ale\"), w)Now we define the Turing model@model constantrates(model, ccd) = begin\n    r  ~ MvLogNormal(ones(2))\n    η  ~ Beta(3,1)\n    q1 ~ Beta()\n    q2 ~ Beta()\n    ccd ~ model((λ=r[1], μ=r[2], η=η, q=[q1, q2]))\nend\n\nmodel = constantrates(w, ccd)\nchain = sample(model, NUTS(0.65), 100)warning: Warning\nOf course such a chain should be run much longer than in this example! Here a very short chain is presented to ensure reasonable build times for this documentation."
},

{
    "location": "wgd-turing/#Using-a-branch-specific-rates-model-1",
    "page": "Bayesian inference using Turing.jl",
    "title": "Using a branch-specific rates model",
    "category": "section",
    "text": "We\'ll use the same tree as above. The relevant model now is the DLWGD model:params = DLWGD(λ=randn(n), μ=randn(n), q=rand(2), η=rand())\nr = Whale.RatesModel(params, fixed=(:p,))\nw = WhaleModel(r, t, 0.5)\nccd = read_ale(joinpath(@__DIR__, \"../../example/example-1/ale\"), w)Note that the duplication and loss rates should here be specified on a log-scale for the DLWGD model.@model branchrates(model, ccd, ::Type{T}=Matrix{Float64}) where {T} = begin\n    η ~ Beta(3,1)\n    ρ ~ Uniform(-1, 1.)\n    τ ~ truncated(Cauchy(0, 1), 0, Inf)\n    S = [τ 0. ; 0. τ]\n    R = [1. ρ ; ρ 1.]\n    Σ = S*R*S\n    !isposdef(Σ) && return -Inf\n    r = T(undef, 2, n)\n    r[:,1] ~ MvNormal(zeros(2), ones(2))\n    for i=2:n\n        r[:,i] ~ MvNormal(r[:,1], Σ)\n    end\n    q1 ~ Beta()\n    q2 ~ Beta()\n    ccd ~ model((λ=r[1,:], μ=r[2,:], η=η, q=[q1, q2]))\nend\n\nmodel = branchrates(w, ccd)\nchain = sample(model, NUTS(0.65), 100)warning: Warning\nOf course such a chain should be run much longer than in this example! Here a very short chain is presented to ensure reasonable build times for this documentation.Let\'s obtain reconciled treespdf = DataFrame(chain)\nfun = (m, x)-> Array(x) |> x->m((λ=x[3:2:36], μ=x[4:2:36], η=x[end-2], q=x[1:2]))\ntt = TreeTracker(w, ccd[end-1:end], pdf, fun)\ntrees = track(tt)Let\'s have a look at the first familytrees[1].treesOr maybe all the gene pairsps = Whale.getpairs(trees, w);\nnothing #hideNow let\'s look at the gene pairs which have a non-zero posterior probability of being derived from WGD node 18 (the Arabidopsis WGD, execute @show w to check the model structure)p = filter(x->x[Symbol(\"18_wgd\")] > 0.0, ps)[!,:pair]The full (approximate) probability distribution over reconciliation events for this gene pair isrow = ps[ps[!,:pair] .== p[1],1:end-2]\nfor (k,v) in zip(names(row), Array(row))\n    v > 0. && println(k, \": \", v)\nendThe following can also be helpfultables = Whale.getwgdtables(trees, ccd, w)\ntables[1]This page was generated using Literate.jl."
},

{
    "location": "cytp450/#",
    "page": "Reconciled tree inference example",
    "title": "Reconciled tree inference example",
    "category": "page",
    "text": ""
},

{
    "location": "cytp450/#Reconciled-tree-inference-example-1",
    "page": "Reconciled tree inference example",
    "title": "Reconciled tree inference example",
    "category": "section",
    "text": "using Whale, DynamicHMC, DynamicHMC.Diagnostics, Random, NewickTree, DataFrames\nRandom.seed!(624)This is a use case I haven\'t been exploring before, namely large gene families. Here we consider a family of about 100 leaves.base  = joinpath(@__DIR__, \"../../example/example-4\")\ntree  = readnw(readline(joinpath(base, \"tree.nw\")))\nmodel = WhaleModel(RatesModel(ConstantDLWGD(λ=0.1, μ=0.2, η=0.9)), tree, .1)\ndata  = read_ale(joinpath(base, \"cytp450.ale\"), model, true)Reading in the single CCD is already a very heavy operation. The CCD has about 5000 unique clades.We will use the DynamicHMC interface, using a constant-rates model (i.e. a single duplication and loss rate for the entire tree).prior = Whale.CRPrior()\nproblem = WhaleProblem(data, model, prior)Now run the actual HMC samplerresults = mcmc_with_warmup(Random.GLOBAL_RNG, problem, 100,\n    warmup_stages=DynamicHMC.default_warmup_stages(doubling_stages=2))\nposterior = Whale.transform(problem, results.chain)\n@info summarize_tree_statistics(results.tree_statistics)We can get a data framedf = Whale.unpack(posterior)\ndescribe(df)Get reconciled trees from the posterior rectrees = track(problem, posterior[1:10])trees = track(problem, posterior)\ntrees[1].treesNote that there are many trees with similar posterior probability!Now we will plot the MAP (maximum a posteriori) reconciled tree, showing duplication and loss eventsusing PalmTree, Luxor\nimport Luxor: RGB\n\nrectree = trees[1].trees[1].tree\noutdir  = mkpath(joinpath(@__DIR__, \"../assets/\"))\noutpath = joinpath(outdir, \"cytp450-map.svg\")\n\ntl = TreeLayout(rectree, cladogram=true, dims=(400,800))\n@svg begin\n    Luxor.origin(Point(0,20))\n    Luxor.setline(1)\n    setfont(\"Noto sans italic\", 7)\n    colfun = n->n.data.label != \"loss\" ? RGB() : RGB(0.99,0.99,0.99)\n    drawtree(tl, color=colfun)\n    nodemap(tl, prewalk(rectree),\n        (n, p) -> !isleaf(n) ?\n            settext(\"  $(n.data.cred)\", p, valign=\"center\") :\n            settext(\"  $(n.data.name)\", p, valign=\"center\"))\n    nodemap(tl, prewalk(rectree),\n        (n, p) -> n.data.label == \"duplication\" && box(p, 4, 4, :fill))\nend 500 850 outpathSquares show duplication events and internal node labels show the posterior probability of observing the relevant node in the reconciled tree.This page was generated using Literate.jl."
},

]}
