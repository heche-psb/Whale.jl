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
    "text": "In this example case, the basic workflow for assessing WGD hypotheses using Whale will be illustrated. We will use the DynamicHMC library for Bayesian inference.note: Note\nInference in Whale with DynamicHMC.jl supports distributed computing. To use distributed parallelism, start up julia with -p <ncores> (or do using Distributed; addprocs(ncores). Instead of loading Whale with using Whale, use @everywhere using Whale and normally all log-likelihood computations should now run in parallel.using Whale, DynamicHMC, Random, NewickTree, Distributions, DataFrames\nusing DynamicHMC.Diagnostics\nRandom.seed!(562);\nnothing #hideSet up the model and the data, here I will use a model with constant duplication and loss rates across the species tree. Note that the tree contains two WGD events (internal nodes labeled with a name starting with wgd).tree  = readnw(\"((MPOL:4.752,(PPAT:2.752)wgd_1:2.0):0.292,(SMOE:4.457,((((OSAT:1.555,(ATHA:0.5548,CPAP:0.5548):1.0002):0.738,ATRI:2.293):1.0)wgd_2:0.225,(GBIL:3.178,PABI:3.178):0.34):0.939):0.587);\")\nn = length(postwalk(tree))\nntaxa = (n+1)÷2\nrates = RatesModel(ConstantDLWGD(λ=0.1, μ=0.1, q=[0.2, 0.3], η=0.9))\nmodel = WhaleModel(rates, tree, 0.1)\ndata  = read_ale(joinpath(@__DIR__, \"../../example/example-1/ale\"), model, true)note: Note\nTo use the DynamicHMC interface, the third argument of read_ale should be set to true.And next we set up the Bayesian inference \'problem\', using the default priors:prior = Whale.CRPrior()\nproblem = WhaleProblem(data, model, prior)Now we run NUTS (of course this is a ridicuously short run, and in reality you want to use something like a 1000 iterations. Also, it\'s better to keep doubling_stages >= 3).results = mcmc_with_warmup(Random.GLOBAL_RNG, problem, 200,\n    warmup_stages=DynamicHMC.default_warmup_stages(doubling_stages=2))\nsummarize_tree_statistics(results.tree_statistics)Now we obtain the posterior distribution in the form of a data frameposterior = Whale.transform(problem, results.chain)\ndf = Whale.unpack(posterior)\ndescribe(df, :mean, :q025=>x->quantile(x, 0.025), :q975=>x->quantile(x, 0.975))And we can visualize the marginal posterior distributions:using Plots\nkwargs = (bins=20, color=:white, grid=false, legend=false)\nps1 = [histogram(df[!,x], xlabel=x; kwargs...) for x in names(df)]\nkwargs = (color=:black, grid=false, legend=false)\nps2 = [plot(df[!,x], xlabel=x; kwargs...) for x in names(df)]\nplot(ps1..., ps2...)From these results (NB: which are based on a mere 12 gene families), we find little support for the second genome duplication (in the angiosperm branch), i.e. the retention rate q_2 is not markedly different from 0. The WGD on the P. patens tip branch however seems to gain some support, with a posterior mean retention rate (q_1) of about 0.4, which is quite high. However, this is definitely too small a data set to make substantial conclusions!We can obtain reconciled trees sampled from the posteriortrees = track(problem, posterior)Consider the first gene familyfamily1 = trees[1].treesNote that the freq field gives the approximate posterior probability of this tree (estimated by the sample frequency). We can get the MAP tree as a newick stringnwstr(family1[1].tree)Now we\'ll plot the MAP treeusing PalmTree, Luxor\nimport Luxor: RGB\n\noutdir  = mkpath(joinpath(@__DIR__, \"../assets/\"))\noutpath = joinpath(outdir, \"dhmc-fam1-map.svg\")\n\nrectree = family1[1].tree\ntl = TreeLayout(rectree, cladogram=true, dims=(300,300))\n@svg begin\n    Luxor.origin(Point(0,20))\n    Luxor.setline(2)\n    setfont(\"Noto sans italic\", 12)\n    colfun = n->n.data.label != \"loss\" ? RGB() : RGB(0.99,0.99,0.99)\n    drawtree(tl, color=colfun)\n    nodemap(tl, prewalk(rectree),\n        (n, p) -> !isleaf(n) ?\n            settext(\"  $(n.data.cred)\", p, valign=\"center\") :\n            settext(\"  $(split(n.data.name, \"_\")[1])\", p, valign=\"center\"))\n    nodemap(tl, prewalk(rectree),\n        (n, p) -> n.data.label == \"duplication\" && box(p, 8, 8, :fill))\n    nodemap(tl, prewalk(rectree),\n        (n, p) -> startswith(n.data.label, \"wgd\") && star(p,3,5,3,0.5,:fill))\nend 500 400 outpathThe support values are posterior probabilities for the associated reconciled split. Note that the tree does not contain branch lengths. Duplication events are marked by squares, whereas the WGDs are marked by stars.The events field for each gene family contains a summary of the expected number of events for each branch (where each branch is identified by the node to which the branch leads, as shown in the node column)trees[1].eventsWe can get for every gene pair the posterior reconciliation probability. The following data frame can therefore be used to probabilistically assess whether two homologous genes are orthologs, WGD-derived paralogs or non-WGD derived paralogs.pair_pps = Whale.getpairs(trees, model)\nfirst(pair_pps, 5)Every row of this data frame is a probability distribution over reconciliation events, so each row sums to one, as we can verify:map(sum, eachrow(pair_pps[!,1:end-2]))Take for instance the following gene pair (second row)x = pair_pps[2,1:end-2]\nfor (n, v) in zip(names(x), Array(x))\n    v > 0 && println(n, \": \", v)\nendThe posterior probability (under the DL model) that this gene pair traces back to the speciation corresponding to node 17 (i.e. the root) is approximately 0.86, whereas the posterior probability that this gene pair traces back to an ancestral duplication event is 0.14.We can get a WGD-centric view as well. The following retrieves a table for each WGD with all gene tree nodes that have a non-zero posterior probability of being reconciled to that particular WGD nodetables = Whale.getwgdtables(trees, data, model)\ntablesThis page was generated using Literate.jl."
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
    "text": "In this example we will use the probabilistic programming language implemented in Turing.jl with Whale to specify Bayesian hierarchical models for gene tree reconciliation in a flexible wayusing Whale, NewickTree, Distributions, Turing, DataFrames, LinearAlgebra, Random\nRandom.seed!(7137);\nnothing #hide"
},

{
    "location": "wgd-turing/#Using-a-constant-rates-model-1",
    "page": "Bayesian inference using Turing.jl",
    "title": "Using a constant-rates model",
    "category": "section",
    "text": "First we will do inference for a simple constant-rates model (i.e. assuming a single duplication and loss rate for the entire species tree). First we load the species tree (using the example tree available in the WHale library)t = deepcopy(Whale.extree)\nn = length(postwalk(t))  # number of internal nodesNow we add two WGD nodes to the tree. We do this by specifying the last common ancestor node for the lineages that share the hypothetical WGD. By default, the added node is halfway between the specified node and its parent.insertnode!(getlca(t, \"PPAT\", \"PPAT\"), name=\"wgd_1\")\ninsertnode!(getlca(t, \"ATHA\", \"ATRI\"), name=\"wgd_2\")and we obtain a reference model object, using the constant-rates model with two WGDsθ = ConstantDLWGD(λ=0.1, μ=0.2, q=[0.2, 0.1], η=0.9)\nr = Whale.RatesModel(θ, fixed=(:p,))\nw = WhaleModel(r, t, .1)next we get the data (we need a model object for that)ccd = read_ale(joinpath(@__DIR__, \"../../example/example-1/ale\"), w)Now we define the Turing model@model constantrates(model, ccd) = begin\n    r  ~ MvLogNormal(ones(2))  # prior on the duplication and loss rate\n    η  ~ Beta(3,1)  # hyperprior for the parameter of the geometric prior distribution on the number of genes at the root of the species tree\n    q1 ~ Beta()  # prior for the WGD retention rate of `wgd_1`\n    q2 ~ Beta()  # prior for the WGD retention rate of `wgd_2`\n    ccd ~ model((λ=r[1], μ=r[2], η=η, q=[q1, q2]))\nend\n\nmodel = constantrates(w, ccd)\nchain = sample(model, NUTS(0.65), 100)warning: Warning\nOf course such a chain should be run much longer than in this example! Here a very short chain is presented to ensure reasonable build times for this documentation."
},

{
    "location": "wgd-turing/#Using-a-branch-specific-rates-model-1",
    "page": "Bayesian inference using Turing.jl",
    "title": "Using a branch-specific rates model",
    "category": "section",
    "text": "Now we will consider a model with branch-specific duplication and loss rates, using a more complicated hierarchical model with an bivariate uncorrelated relaxed clock prior. We\'ll use the same tree as above. The relevant model now is the DLWGD model:params = DLWGD(λ=randn(n), μ=randn(n), q=rand(2), η=rand())\nr = Whale.RatesModel(params, fixed=(:p,))\nw = WhaleModel(r, t, 0.5)\nccd = read_ale(joinpath(@__DIR__, \"../../example/example-1/ale\"), w)Note that the duplication and loss rates should here be specified on a log-scale for the DLWGD model. We use an LKJ prior for the covariance matrix, specifying a prior for the correlation of duplication and loss rates (ρ) and a prior for the scale parameter τ, see e.g. the stan docs@model branchrates(model, ccd, ::Type{T}=Matrix{Float64}) where {T} = begin\n    η ~ Beta(3,1)\n    ρ ~ Uniform(-1, 1.)\n    τ ~ truncated(Cauchy(0, 1), 0, Inf)\n    S = [τ 0. ; 0. τ]\n    R = [1. ρ ; ρ 1.]\n    Σ = S*R*S\n    !isposdef(Σ) && return -Inf\n    r = T(undef, 2, n)\n    r[:,1] ~ MvNormal(zeros(2), ones(2))\n    for i=2:n\n        r[:,i] ~ MvNormal(r[:,1], Σ)\n    end\n    q1 ~ Beta()\n    q2 ~ Beta()\n    ccd ~ model((λ=r[1,:], μ=r[2,:], η=η, q=[q1, q2]))\nend\n\nmodel = branchrates(w, ccd)\nchain = sample(model, NUTS(0.65), 100)warning: Warning\nOf course such a chain should be run much longer than in this example! Here a very short chain is presented to ensure reasonable build times for this documentation.Now let\'s obtain reconciled treespdf = DataFrame(chain)\nfun = (m, x)-> Array(x) |> x->m((λ=x[3:2:36], μ=x[4:2:36], η=x[end-2], q=x[1:2]))\ntt = TreeTracker(w, ccd[end-1:end], pdf, fun)\ntrees = track(tt)Let\'s have a look at the first familytrees[1].treesOr maybe all the gene pairsps = Whale.getpairs(trees, w);\nnothing #hideNow let\'s look at the gene pairs which have a non-zero posterior probability of being derived from WGD node 18 (the Arabidopsis WGD, execute @show w to check the model structure)p = filter(x->x[Symbol(\"18_wgd\")] > 0.0, ps)[!,:pair]The full (approximate) probability distribution over reconciliation events for this gene pair isrow = ps[ps[!,:pair] .== p[1],1:end-2]\nfor (k,v) in zip(names(row), Array(row))\n    v > 0. && println(k, \": \", v)\nendThe following can also be helpfultables = Whale.getwgdtables(trees, ccd, w)\ntablesThis page was generated using Literate.jl."
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
    "text": "using Whale, DynamicHMC, DynamicHMC.Diagnostics, Random, NewickTree, DataFrames\nRandom.seed!(624);\nnothing #hideIn this case study, we will perform Bayesian gene tree reconciliation for a single (large) gene family. The data can be found in the example4 directory in the Whale git repository. We first load the data:base  = joinpath(@__DIR__, \"../../example/example-4\")\ntree  = readnw(readline(joinpath(base, \"tree.nw\")))\nmodel = WhaleModel(RatesModel(ConstantDLWGD(λ=0.1, μ=0.2, η=0.9)), tree, .1)\ndata  = read_ale(joinpath(base, \"cytp450.ale\"), model, true)Reading in the single CCD (in the read_ale step is already a rather heavy operation. The CCD has about 5000 unique clades.For Bayesian inference we will use the DynamicHMC interface, using a constant-rates model (i.e. assuming a single duplication and loss rate for the entire tree). The default prior should of course not be chosen lightly, although for our current purposes it is reasonable:prior = Whale.CRPrior()\nproblem = WhaleProblem(data, model, prior)Now we run the actual HMC sampler. Note that we perform a very short run here to reduce build times of the documentation, in reality you\'d rather use something like a 1000 iterations.results = mcmc_with_warmup(Random.GLOBAL_RNG, problem, 100,\n    warmup_stages=DynamicHMC.default_warmup_stages(doubling_stages=2))\nposterior = Whale.transform(problem, results.chain)\n@info summarize_tree_statistics(results.tree_statistics)A data frame may be easier to work with (and save to disk)df = Whale.unpack(posterior)\ndescribe(df, :mean, :q025=>x->quantile(x, 0.025), :q975=>x->quantile(x, 0.975))Now we will obtain reconciled trees from the posteriortrees = track(problem, posterior)\ntrees[1].treesNote that there are many trees with similar posterior probability, so in other words the maximum a posteriori (MAP) tree is not that meaningful in itself. We can however plot the MAP tree with posterior node probabilities to get an idea of the reconciled tree and the nodes with considerable posterior uncertainty. I will use Luxor.jl together with my small helper library for plotting trees:using PalmTree, Luxor\nimport Luxor: RGB\n\nrectree = trees[1].trees[1].tree\noutdir  = mkpath(joinpath(@__DIR__, \"../assets/\"))\noutpath = joinpath(outdir, \"cytp450-map.svg\")\n\ntl = TreeLayout(rectree, cladogram=true, dims=(400,800))\n@svg begin\n    Luxor.origin(Point(0,20))\n    Luxor.setline(1)\n    setfont(\"Noto sans italic\", 7)\n    colfun = n->n.data.label != \"loss\" ? RGB() : RGB(0.99,0.99,0.99)\n    drawtree(tl, color=colfun)\n    nodemap(tl, prewalk(rectree),\n        (n, p) -> !isleaf(n) ?\n            settext(\"  $(n.data.cred)\", p, valign=\"center\") :\n            settext(\"  $(n.data.name)\", p, valign=\"center\"))\n    nodemap(tl, prewalk(rectree),\n        (n, p) -> n.data.label == \"duplication\" && box(p, 4, 4, :fill))\nend 500 850 outpathSquares show duplication events and internal node labels show the posterior probability of observing the relevant node in the reconciled tree.This page was generated using Literate.jl."
},

]}
