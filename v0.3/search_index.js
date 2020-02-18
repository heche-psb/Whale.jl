var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "#Introduction-1",
    "page": "Introduction",
    "title": "Introduction",
    "category": "section",
    "text": "warning: Warning\nThe latest Whale version is a thorough rewrite of the Whale library, and is still work in progress. For the version as used in Zwaenepoel & Van de Peer (2019), refer to this release (v0.2). Nevertheless, the current version should be safe to use and is more efficient and convenient (if you know a bit of julia).Whale provides tools for genome-wide amalgamated likelihood estimation (ALE) under a DL+WGD model, which is an approach to infer reconciled gene trees and parameters of a model of gene family evolution given a known species tree.The ALE approach takes into account uncertainty in the gene tree topology by marginalizing over all tree topologies that can be amalgamated from the so-called conditional clade distribution (CCD). This CCD can be constructed from a sample of the posterior distribution of tree topologies (which can be obtained using any standard software for Bayesian phylogenetics).More specifically, this library can be used toStatistically test the absence or presence of hypothetical whole-genome duplication (WGD) events in a species phylogeny\nInfer lineage-specific gene duplication and loss rates for a species phylogeny\nInfer high-quality (reconciled) gene trees given a known species tree cf. Szöllősi et al.\nAll of the above at oncenote: Note\nThis library implements the DL and DL+WGD models. It does not implement models of gene family evolution that take into account horizontal gene transfer, incomplete lineage sorting or gene conversion."
},

{
    "location": "#Installation-1",
    "page": "Introduction",
    "title": "Installation",
    "category": "section",
    "text": "You will need julia-1.x. Fire up a julia REPL by typing julia at the command line enter the package manager interface by typing ], and run the command(v1.1) pkg> add https://github.com/arzwa/PhyloTrees.jl\n(v1.1) pkg> add https://github.com/arzwa/ConsensusTrees.jl\n(v1.1) pkg> add https://github.com/arzwa/Whale.jl"
},

{
    "location": "#Data-preparation-1",
    "page": "Introduction",
    "title": "Data preparation",
    "category": "section",
    "text": "To perform analyses with Whale, you will need  An ultrametric species tree, with ideally branch lengths in geological time (since this allows straightforward interpretation of parameter estimates.)\nA bunch of ALE files, which summarize the conditional clade distributions (CCDs) for the same bunch of gene families. These can be obtained from a sample of the posterior distribution of gene trees using the ALEobserve tool. A pipeline to obtain these from a set of gene family protein fasta files is available at github.note: Note\nGene IDs should be prefixed by the name of the species to which the gene belongs as used in the species tree. For example if Arabidopsis thaliana is represented by ATHA in the species tree newick file, then the genes should be prefixed with ATHA_, e.g. ATHA_AT1G05000.note: Note\nAnalyzing CCDs (ALE files) with a very large number of clades or for very large families can be prohibitive computationally. It is therefore generally advisable that large orthogroups are filtered out based on some criterion (for example using the script orthofilter.py in the scripts directory of the Whale repository). To filter out families with very large numbers of clades in the CCD (which reflects that there is a lot of uncertainty in the gene tree), the scripts ccddata.py and ccdfilter.py can be used. This is a rather ad hoc filtering procedure, but can be useful to filter out families that trouble the analysis.warning: Warning\nMost analyses in Whale assume that for each family, there is at least one gene in both clades stemming from the root of the species tree. The likelihood in Whale is the conditional likelihood under this assumption. This is to rule out the possibility of de novo gain of a gene family along a branch of the species tree. The orthogroup data should therefore always be filtered to be in accordance with this criterion. This can also be done using the orthofilter.py script."
},

{
    "location": "#Quick-start-1",
    "page": "Introduction",
    "title": "Quick start",
    "category": "section",
    "text": "If you\'re not familiar with julia, and you simply want to run analyses as performed for instance in Zwaenepoel & Van de Peer (2019) the following scripts will be helpful. If you want to get a more detailed view of the library, please consult the Manual."
},

{
    "location": "#Maximum-likelihood-estimation-1",
    "page": "Introduction",
    "title": "Maximum likelihood estimation",
    "category": "section",
    "text": "I assume that your tree is called tree_file.nw and your .ale files are in a directory ccd_dir. For maximum likelihood estimation with a constant rates model, the following script should work@everywhere using Whale\n\n# data and config\nst = SlicedTree(\"tree_file.nw\")\nccd = read_ale(\"ccd_dir\", st)\nconstant = true  # set to false for branch-wise rates\n\n# inference\nconstant ? set_constantrates!(st) : nothing\nw = WhaleModel(st)\nout = mle(w, ccd)Save this script to a file (say whale-mle.jl). To run the inference using 16 processors for instance, run julia -p 16 whale-mle.jl. Do not forget the @everywhere before the using statement, as this will load the Whale library on all available processors.To add WGDs, a wgd_conf argument should be provide to SLicedTree, please see the docs for SlicedTree. See the manual section on Rate indices on how to specify local-clock models.warning: Warning\nDepending on the time scale and data set at hand, you may need to tweak the initial values of the WhaleModel to prevent Numerical issues!warning: Warning\nML estimation with branch-wise rates may result in poor convergence for large species trees."
},

{
    "location": "#Bayesian-inference-with-MCMC-1",
    "page": "Introduction",
    "title": "Bayesian inference with MCMC",
    "category": "section",
    "text": "A similar script for Bayesian inference using MCMC looks like (here with a hypothetical WGD configuration)@everywhere using Whale\n\n# data and config\nwgd_conf = Dict(\n    \"wgd1\" => (\"taxon1,taxon2\", 1.2),  # (LCA, time BP)\n    \"wgd2\" => (\"taxon1,taxon5\", 3.2)\n)\nst = SlicedTree(\"tree_file.nw\", wgd_conf)\nccd = read_ale(\"ccd_dir\", st)\nmodel = IRModel(st)      # independent log-normal rates\n# model = GBMModel(st)   # autocorrelated rates\nn = 11000\n\nw = WhaleChain(st, model)\nchain = mcmc!(w, ccd, n, show_every=10)again, saving this as whale-bay.jl and running this with julia -p <nCPU> whale-bay.jl will start the MCMC with the likelihood evaluation performed on nCPU cores.note: Note\nYou might want to save the MCMC simulation to a file. In the above example one can either save the chain variable to a file using JLD, save the data frame in the field chain.df using CSV.write or alternatively one can just set the show_every argument to the desired thinning level and  redirect stdout to a file while running the MCMC.warning: Warning\nPlease take you time to understand the hierarchical model used in Whale and to modify the prior distributions to suit your data set! In particular, note that mixing can be very poor when the η and/or ν parameters are considere random variables and assigned hyperpriors. Fixing η and/or ν is therefore often necessary. Please consult the Bayesian inference section of the manual. "
},

{
    "location": "#References-1",
    "page": "Introduction",
    "title": "References",
    "category": "section",
    "text": "Whale.jl is developed by Arthur Zwaenepoel at the VIB-UGent center for plant systems biology (bioinformatics & evolutionary genomics group). If you use Whale, please cite:Zwaenepoel, A. and Van de Peer, Y., 2019. Inference of Ancient Whole-Genome Duplications and the Evolution of Gene Duplication and Loss Rates. Molecular biology and evolution, 36(7), pp.1384-1404.The methods in Whale are heavily inspired by previous work done by other researchers. If you use Whale, consider citing the following two particularly important studies:[ALE] Szöllősi, G.J., Rosikiewicz, W., Boussau, B., Tannier, E. and Daubin, V., 2013. Efficient exploration of the space of reconciled gene trees. Systematic biology, 62(6), pp.901-912.[DL+WGD model] Rabier, C.E., Ta, T. and Ané, C., 2013. Detecting and locating whole genome duplications on a phylogeny: a probabilistic approach. Molecular biology and evolution, 31(3), pp.750-762."
},

{
    "location": "manual/#",
    "page": "Manual",
    "title": "Manual",
    "category": "page",
    "text": ""
},

{
    "location": "manual/#Manual-1",
    "page": "Manual",
    "title": "Manual",
    "category": "section",
    "text": "Below the major components of the Whale library are discussed.DocTestSetup = quote\n    using Whale\n    using Random\n    Random.seed!(1234)\nend"
},

{
    "location": "manual/#Sliced-species-tree-1",
    "page": "Manual",
    "title": "Sliced species tree",
    "category": "section",
    "text": "The ALE approach to probabilistic gene tree - species tree reconciliation uses a discretization of the branches of the species tree into small time intervals. This \'sliced\' species tree defines the main structure of the model.julia> st = Whale.example_tree()\nSlicedTree(9, 17, 7)\n\njulia> st.tree\n\r[0mPhylogenetic tree with 24 nodes and 23 branches\n\njulia> st.leaves\nDict{Int64,String} with 9 entries:\n  4  => \"PPAT\"\n  13 => \"CPAP\"\n  10 => \"OSAT\"\n  14 => \"ATRI\"\n  3  => \"MPOL\"\n  16 => \"GBIL\"\n  17 => \"PABI\"\n  6  => \"SMOE\"\n  12 => \"ATHA\"\n\njulia> wgds(st)  # WGD ID → WGD node → q index\nPPAT → node 18 → q1\nCPAP → node 19 → q2\nBETA → node 20 → q3\nANGI → node 21 → q4\nSEED → node 22 → q5\nMONO → node 23 → q6\nALPH → node 24 → q7\n\njulia> st[3, 4]  # length of 4th slice in branch 3\n0.049499999999999995\n\njulia> nslices(st, 3)  # number of slices in branch 3\n97To get a tree in Newick format into a SlicedTree, one can simply use SlicedTree(tree_file).note: Note\nNote that the tree is assumed to be ultrametric and that you might need to change the default Δt value for your purposes. WGDs can be specified by using a configuration dictionary (see SlicedTree).For visualizing tree structures, the PalmTree library can be used. It is often useful for example to plot the tree with internal node labels for specifying models in Whalejulia> using PalmTree\njulia> drawtree(st, nodelabels=true)Here, nodes 18 to 24 are WGD \'nodes\', marking hypothetical WGDs along the sliced species tree.(Image: )"
},

{
    "location": "manual/#Rate-indices-1",
    "page": "Manual",
    "title": "Rate indices",
    "category": "section",
    "text": "The SlicedTree structure has two fields that store mappings from nodes/branches in the tree to indices in hypothetical parameter vectors. The qindex field is a mapping (Dict) from WGD nodes to indiced for a vector of retention rates, whereas the rindex serves as a mapping from species tree branches to indices for the duplication and loss rate vectors. The default rindex has a different index for each branch of the species tree, and with the same index for the part of a branch before and after a WGD (note that branches are identified by the index of there downstream (leafward) node).julia> st = Whale.example_tree();\n\njulia> st.qindex\nDict{Int64,Int64} with 7 entries:\n  20 => 3\n  23 => 6\n  24 => 7\n  19 => 2\n  21 => 4\n  22 => 5\n  18 => 1\n\njulia> st.rindex\nDict{Int64,Int64} with 24 entries:\n  18 => 4\n  2  => 2\n  16 => 16\n  11 => 11\n  21 => 8\n  7  => 7\n  9  => 9\n  10 => 10\n  19 => 13\n  17 => 17\n  8  => 8\n  22 => 7\n  6  => 6\n  24 => 12\n  4  => 4\n  3  => 3\n  5  => 5\n  20 => 12\n  23 => 10\n  ⋮  => ⋮In this example, branches 20, 24 and 12 (which are all part of the same species tree branch but refer to different segments marked by WGD nodes) all point to index 12, which means that they are associated with the same duplication and loss rates.The rindex can be modified to specify arbitrary rate models (for instance fixing a particular clade to a one shared duplication and loss rate). In order to specify a consjuliatant-rates model, one can dojulia> st = Whale.example_tree();\n\njulia> set_constantrates!(st)\n\njulia> st.rindex\nDict{Int64,Int64} with 24 entries:\n  18 => 1\n  2  => 1\n  16 => 1\n  11 => 1\n  21 => 1\n  7  => 1\n  9  => 1\n  10 => 1\n  19 => 1\n  17 => 1\n  8  => 1\n  22 => 1\n  6  => 1\n  24 => 1\n  4  => 1\n  3  => 1\n  5  => 1\n  20 => 1\n  23 => 1\n  ⋮  => ⋮"
},

{
    "location": "manual/#Conditional-clade-distribution(s)-1",
    "page": "Manual",
    "title": "Conditional clade distribution(s)",
    "category": "section",
    "text": "The conditional clade distributions (CCDs) for a set of gene families provide the main input data (observations) for Whale analyses. These can be read from .ale files generated by ALEobserve. The read_ale function accepts either a single .ale file, a text file with on each line the path to a .ale file or a directory of .ale files. When an empty file is provided, a dummy CCD object will be created (which is useful when one wants to run an MCMC chain without data to check the prior specification).julia> st = Whale.example_tree();\n\njulia> ccd = read_ale(\"../example/example-ale/\", st)\n[ Info:  .. read 12 ALE files\n12-element DistributedArrays.DArray{CCD,1,Array{CCD,1}}:\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 83 clades, 5001 samples)\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 55 clades, 5001 samples)\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 89 clades, 5001 samples)\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 131 clades, 5001 samples)\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 107 clades, 5001 samples)\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 59 clades, 5001 samples)\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 53 clades, 5001 samples)\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 83 clades, 5001 samples)\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 59 clades, 5001 samples)\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 95 clades, 5001 samples)\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 67 clades, 5001 samples)\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 65 clades, 5001 samples)By default, read_ale will distribute the resulting CCD array over all available processors."
},

{
    "location": "manual/#The-WhaleModel-1",
    "page": "Manual",
    "title": "The WhaleModel",
    "category": "section",
    "text": "The last object of importance to do inference with Whale is the WhaleModel type. This structure is used for computing the probability of observing the data conditional on the model and its parameters (logpdf).julia> st = Whale.example_tree();\n\njulia> ccd = read_ale(\"../example/example-ale/\", st);\n[ Info:  .. read 12 ALE files\n\njulia> w = WhaleModel(st)\nWhaleModel{Float64,CCD}(\nλ: [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]\nμ: [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3]\nq: [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]\nη: 0.9\n)\n\njulia> logpdf(w, ccd[1])  # single CCD\n-32.70481615232666\n\njulia> logpdf(w, ccd)     # multiple CCDs, distributed by default\n-298.98170493684256\n\njulia> w = WhaleModel(st, rand(17), rand(17), rand(7), 0.66)  # full constructor\nWhaleModel{Float64,CCD}(\nλ: [0.590845, 0.766797, 0.566237, 0.460085, 0.794026, 0.854147, 0.200586, 0.298614, 0.246837, 0.579672, 0.648882, 0.0109059, 0.066423, 0.956753, 0.646691, 0.112486, 0.276021]\nμ: [0.651664, 0.0566425, 0.842714, 0.950498, 0.96467, 0.945775, 0.789904, 0.82116, 0.0341601, 0.0945445, 0.314926, 0.12781, 0.374187, 0.931115, 0.438939, 0.246862, 0.0118196]\nq: [0.0460428, 0.496169, 0.732, 0.299058, 0.449182, 0.875096, 0.0462887]\nη: 0.66\n)\n\njulia> logpdf(w, ccd)\n-370.067694892539An informative description of the model can be printed using describejulia> st = Whale.example_tree();\n\njulia> w = WhaleModel(st);\n\njulia> describe(w)\nLeaves\n======\n4 	→ PPAT\n13 	→ CPAP\n10 	→ OSAT\n14 	→ ATRI\n3 	→ MPOL\n16 	→ GBIL\n17 	→ PABI\n6 	→ SMOE\n12 	→ ATHA\nRates (λ, μ)\n============\n3 	| λ, μ = 0.2,0.3	| (3)\n4 	| λ, μ = 0.2,0.3	| (4)\n18 	| λ, μ = 0.2,0.3	| (4)\n2 	| λ, μ = 0.2,0.3	| (4,3)\n6 	| λ, μ = 0.2,0.3	| (6)\n16 	| λ, μ = 0.2,0.3	| (16)\n17 	| λ, μ = 0.2,0.3	| (17)\n15 	| λ, μ = 0.2,0.3	| (16,17)\n13 	| λ, μ = 0.2,0.3	| (13)\n19 	| λ, μ = 0.2,0.3	| (13)\n12 	| λ, μ = 0.2,0.3	| (12)\n24 	| λ, μ = 0.2,0.3	| (12)\n20 	| λ, μ = 0.2,0.3	| (12)\n11 	| λ, μ = 0.2,0.3	| (13,12)\n10 	| λ, μ = 0.2,0.3	| (10)\n23 	| λ, μ = 0.2,0.3	| (10)\n9 	| λ, μ = 0.2,0.3	| (13,10,12)\n14 	| λ, μ = 0.2,0.3	| (14)\n8 	| λ, μ = 0.2,0.3	| (13,10,14,12)\n21 	| λ, μ = 0.2,0.3	| (13,10,14,12)\n7 	| λ, μ = 0.2,0.3	| (13,10,14,16,17,12)\n22 	| λ, μ = 0.2,0.3	| (13,10,14,16,17,12)\n5 	| λ, μ = 0.2,0.3	| (13,10,14,16,17,6,12)\n1 	| λ, μ = 0.2,0.3	| (4,13,10,14,3,16,17,6,12)\nWGDs (q)\n========\n20, q = 0.2\n23, q = 0.2\n24, q = 0.2\n19, q = 0.2\n21, q = 0.2\n22, q = 0.2\n18, q = 0.2\nOther\n=====\n   η = 0.9\ncond = oibnote: Note\nThe default initial rate values (~0.2) might not be appropriate for you data set and lead to numerical difficulties. Good initial values depend on the unit of time the branch lengths of the SlicedTree are expressed in."
},

{
    "location": "manual/#Maximum-likelihood-estimation-1",
    "page": "Manual",
    "title": "Maximum likelihood estimation",
    "category": "section",
    "text": "Maximum likelihood estimation is performed using Optim.jl with ForwardDiff.jl automatic differentiation. By default the LBFGS optimizer is used, but other Optimizers from Optim work as well.julia> st = Whale.example_tree();\n\njulia> set_constantrates!(st)\n\njulia> w = WhaleModel(st, 0.2, 0.3)\nWhaleModel{Float64,CCD}(\nλ: [0.2]\nμ: [0.3]\nq: [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]\nη: 0.9\n)\n\njulia> mle(w, ccd)\nFminbox\n-------\nInitial mu = 0.0038519\n\nFminbox iteration 1\n-------------------\nCalling inner optimizer with mu = 0.0038519\n\n(numbers below include barrier contribution)\nIter     Function value   Gradient norm\n     0     2.990420e+02     9.056730e+01\n    10     2.788972e+02     3.959786e+00\n  ... # a lot more\n, Results of Optimization Algorithm\n * Algorithm: Fminbox with L-BFGS\n * Starting Point: [0.2,0.3,0.2,0.2,0.2,0.2,0.2,0.2,0.2]\n * Minimizer: [0.08839734901061815,0.15058644779285166, ...]\n * Minimum: 2.787698e+02\n * Iterations: 6\n * Convergence: true\n   * |x - x\'| ≤ 0.0e+00: true\n     |x - x\'| = 0.00e+00\n   * |f(x) - f(x\')| ≤ 0.0e+00 |f(x)|: true\n     |f(x) - f(x\')| = 0.00e+00 |f(x)|\n   * |g(x)| ≤ 1.0e-08: false\n     |g(x)| = 1.38e+01\n   * Stopped by an increasing objective: true\n   * Reached Maximum Number of Iterations: false\n * Objective Calls: 1755\n * Gradient Calls: 1755)"
},

{
    "location": "manual/#Bayesian-inference-1",
    "page": "Manual",
    "title": "Bayesian inference",
    "category": "section",
    "text": "Currently, a model-specific MCMC algorithm (following an adaptive metropolis-within-Gibbs scheme) is used. Specifying arbitrary complex models in Turing.jl is possible, but currently does not support distributed likelihood evaluation and is therefore not yet possible for the kinds of problems tackled with Whale. This is a major goal for future developments."
},

{
    "location": "manual/#Independent-rates-model-1",
    "page": "Manual",
    "title": "Independent rates model",
    "category": "section",
    "text": "The default structure of the independent rates model is as follows.begineqnarray\nnu sim mathrmInverseGamma \neta sim mathrmBeta \nlambda_0 sim mathrmExponential \nmu_0 sim mathrmExponential \nlambda_i sim mathrmLogNormal(lambda_0 nu) \nmu_i sim mathrmLogNormal(mu_0 nu) \nq_i sim mathrmBeta(1 1)\nendeqnarrayBut other distributions from the Distributions.jl library can be used. It is also possible to set parameters to fixed values. This is often desirable for either η or ν to ensure proper mixing of the chain.julia> st = Whale.example_tree()\nSlicedTree(9, 17, 7)\n\njulia> w = WhaleChain(st, IRModel(st))\nERROR: MethodError: no method matching IRModel(::SlicedTree)\nClosest candidates are:\n  IRModel(::T, !Matched::U, !Matched::V, !Matched::W, !Matched::X) where {T, U, V, W, X} at /home/travis/.julia/packages/Parameters/l76EM/src/Parameters.jl:499\n  IRModel(; ν, η, λ, μ, q) at /home/travis/.julia/packages/Parameters/l76EM/src/Parameters.jl:518\n  IRModel(!Matched::IRModel; kws...) at /home/travis/.julia/packages/Parameters/l76EM/src/Parameters.jl:528\n  ...\nStacktrace:\n [1] top-level scope at none:0To run the MCMC simulation, use the mcmc! functionjulia> chain = mcmc!(w, D, 100, show_every=10)The resulting Chains object is fairly intuitive, see the docs for MCMCChains.jl."
},

{
    "location": "manual/#Autocorrelated-rates-model-(Geometric-Brownian-motion)-1",
    "page": "Manual",
    "title": "Autocorrelated rates model (Geometric Brownian motion)",
    "category": "section",
    "text": "The default structure is as above but withbegineqnarray\nnu sim mathrmExponential \neta sim mathrmBeta \nlambda_0 sim mathrmExponential \nmu_0 sim mathrmExponential lambda_i sim mathrmGeometricBrownianMotion(lambda_0 nu) \nmu_i sim mathrmGeometricBrownianMotion(mu_0 nu) \nq_i sim mathrmBeta(1 1)\nendeqnarrayjulia> st = Whale.example_tree()\njulia> w = WhaleChain(st, GBMModel(st))\njulia> chain = mcmc!(w, D, 100, show_every=10)"
},

{
    "location": "manual/#MCMC-mixing-issues-1",
    "page": "Manual",
    "title": "MCMC mixing issues",
    "category": "section",
    "text": "When WGD hypotheses and branch-wise rates across the tree are combined, MCMC in the Whale model can be quite sensitive to the (informative) priors used. For some data sets and prior settings the MCMC algorithm may have a hard time converging, or some parameter may wander of in an unrealistic area of parameter space. Since there is usually a lot of data, the influence of the prior is often very limited and the likelihood dominates the posterior (which is of course desirable), and it may be necessary to constrain some elements of the model to attain convergence."
},

{
    "location": "manual/#Fixing-η-and/or-ν-1",
    "page": "Manual",
    "title": "Fixing η and/or ν",
    "category": "section",
    "text": "To prevent poor mixing in the MCMC, it is often necessary to fix the η and/or ν parameters. For η this is usually not very problematic, since it embodies already a distributional assumption that allows for uncertainty in prior beliefs (since it is the parameter of the geometric prior distribution  on the number of lineages at the root). Choosing for instance η = 0.8, The probability of one lineage at the root is 0.8, two lineages 0.16, three lineages  0.032 etc.Fixing ν (which controls the variation in duplication and loss rates across lineages) can be more troublesome since it is hard to specify a cogent prior. In the context of WGD inference however, the \'true\' values of the duplication and loss rates might not matter too much, and we are mostly interested whether allowing more rate variation across the tree alters are posterior beliefs with regard to WGDs. When the goal is WGD inference, it is therefore advisable to run chains for different ν values and see whether this alters the posterior distributions for the retention rates. Often when the rates are constrained to be very similar across the tree (small ν values), some duplication/loss rate variation is captured by the retention rate, and in this case, for larger ν values, a previously significant non-zero retention rate might shift towards zero.Below the chain is fixed for the parameter values η=0.9 and ν=0.1.julia> st = Whale.example_tree()\njulia> w = WhaleChain(st, IRModel(st, 0.1, 0.9))\njulia> chain = mcmc!(w, D, 100, :ν, :η, show_every=10)"
},

{
    "location": "manual/#Constraining-rates-on-branches-stemming-from-the-root-1",
    "page": "Manual",
    "title": "Constraining rates on branches stemming from the root",
    "category": "section",
    "text": "Sometimes, in particular when there is a long outgroup branch (possibly with WGDs), it can help to constrain the branches stemming left and right from the root of the species tree to have the same duplication and loss rates. This can be done as follows:julia> st = Whale.example_tree();\n\njulia> set_equalrootrates!(st);\n\njulia> st.rindex\nDict{Int64,Int64} with 24 entries:\n  18 => 4\n  2  => 2\n  16 => 15\n  11 => 10\n  21 => 7\n  7  => 6\n  9  => 8\n  10 => 9\n  19 => 12\n  17 => 16\n  8  => 7\n  22 => 6\n  6  => 5\n  24 => 11\n  4  => 4\n  3  => 3\n  5  => 2\n  20 => 11\n  23 => 9\n  ⋮  => ⋮\n"
},

{
    "location": "manual/#Sampling-from-the-prior-1",
    "page": "Manual",
    "title": "Sampling from the prior",
    "category": "section",
    "text": "It is generally advisable to run a chain without data, to investigate the prior distributions one has assigned, and validate the correctness of the MCMC algorithm. This can be done as follows:julia> st = Whale.example_tree()\njulia> w = WhaleChain(st, IRModel(st, InverseGamma(10), 0.9, Exponential(1.), Exponential(1.)))\njulia> chain = mcmc!(w, 10000, :η)   # e.g. with fixed η"
},

{
    "location": "manual/#Backtracking-and-consensus-reconciled-trees-1",
    "page": "Manual",
    "title": "Backtracking and consensus reconciled trees",
    "category": "section",
    "text": "warning: Warning\nBacktracking functionalities are still quite experimental and prone to change.To backtrack reconciled trees for a particular CCD from a parameterized WhaleModel and compute majority-vote consensus trees one can use the following methodsjulia> x  # a single CCD\nCCD{Float64,PhyloTrees.RecTree}(13 taxa, 83 clades, 5001 samples)\n\njulia> w = WhaleModel(st, 0.2, 0.3)\nWhaleModel{Float64,CCD}(\nλ: [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]\nμ: [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3]\nq: [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]\nη: 0.9\n)\n\njulia> backtrack!(x, w, 100)\nCCD{Float64,PhyloTrees.RecTree}(13 taxa, 83 clades, 5001 samples)\n\njulia> drawtree(x.rectrs[1)(Image: )note: Note\nCurrently the branch lengths in backtracked trees are not yet meaningful, although they are related to the branch lengths in the CCD.Majority vote consensus trees can than be obtainedjulia> contree = consensus(x, st)\nPhylogenetic tree with 25 nodes and 24 branches, [...]\n\njulia> drawtree(contree)(Image: )Reconciled trees are stored in the rectrs field of the CCD object. By default, during MCMC, every iteration a tree is sampled from the posterior predictive distribution.julia> D = read_ale(\"example/example-ale/\", st);\n[ Info:  .. read 12 ALE files\n\njulia> w = WhaleChain(st, IRModel(st, 0.1, 0.9));\n\njulia> chain = mcmc!(w, D, 100, :ν, :η);\n\njulia> crts = consensus(D, st);(Consensus) reconciled trees can be written in a newick-like format, where support values x1-x2 with x1 the clade support, which is the frequency by which the clade was observed in the sample of reconciled trees, and x2 the reconciliation support, being the frequency of this nodes majority vote reconciliation (e.g. duplication or speciation). Note that loss events do not appear in consensus reconciliations (but they can be easily determined based on the reconciliations of other nodes).julia> write(stdout, crts[1])\n(((PABI_PAB00009793.1_PAB00009793,((((CPAP_Cpa.t.sc25.3_Cpa.g.sc25.3,((CPAP_Cpa.t.sc25.8_Cpa.g.sc25.8,CPAP_Cpa.t.sc25.5_Cpa.g.sc25.5)0.6086956521739131-0.33043478260869563,ATHA_AT5G48120.1_AT5G48120)0.5391304347826087-0.5043478260869565)0.9217391304347826-0.6173913043478261,(OSAT_LOC_Os07g08050.1_LOC_Os07g08050,CPAP_Cpa.t.sc25.4_Cpa.g.sc25.4)0.9043478260869565-0.9043478260869565)1.0-0.9304347826086956,ATRI_ATR0705G185.1_ATR0705G185)1.0-0.991304347826087,(PABI_PAB00012681.1_PAB00012681,GBIL_Gb_13638)1.0-0.9826086956521739)0.991304347826087-0.991304347826087)1.0-0.33043478260869563,SMOE_SMO111G0185.1_SMO111G0185)1.0-0.991304347826087,(MPOL_Mapoly0036s0119.1_Mapoly0036s0119,PPAT_Pp3c9_4950V3.1_Pp3c9_4950)1.0-1.0);(Consensus) reconciled trees can also be written in PhyloRecXML format:julia> write(stdout, crts[1], st, family=D[1].fname)\n<recGeneTree\n    xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n    xmlns=\"http://www.recgenetreexml.org\" \"txsi:schemaLocation=\"http://www.recgenetreexml.org ../../xsd/recGeneTreeXML.xsd\">\n    <phylogeny rooted=\"true\">\n        <id>example/example-ale/OG0004512.fasta.nex.treesample.ale</id>\n        <clade>\n            <name>1</name>\n            <eventsRec><speciation speciesLocation=\"1\"></speciation></eventsRec>\n            <clade>\n                <name>2</name>\n                <eventsRec><speciation speciesLocation=\"5\"></speciation></eventsRec>\n                <clade>\n                    <name>3</name>\n                    <eventsRec><duplication speciesLocation=\"7\"></duplication></eventsRec>\n                    <clade>\n                        <name>PABI_PAB00009793.1_PAB00009793</name>\n                        <eventsRec><leaf speciesLocation=\"PABI\" geneName=\"PABI_PAB00009793.1_PAB00009793\"></leaf></eventsRec>\n                    </clade>\n                    [...]  # a lot more\n                    </clade>\n               <clade>\n                   <name>PPAT_Pp3c9_4950V3.1_Pp3c9_4950</name>\n                   <eventsRec><leaf speciesLocation=\"PPAT\" geneName=\"PPAT_Pp3c9_4950V3.1_Pp3c9_4950\"></leaf></eventsRec>\n               </clade>\n           </clade>\n       </clade>\n   </phylogeny>\n</recGeneTree>"
},

{
    "location": "api/#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "api/#Whale.CCD",
    "page": "API",
    "title": "Whale.CCD",
    "category": "type",
    "text": "CCD{<:Real,RecTree}\n\nConditional clade distribution with many helper fields. See read_ale for details on IO.\n\njulia> x\nCCD{Float64,PhyloTrees.RecTree}(13 taxa, 83 clades, 5001 samples)\n\njulia> x.ccp\nDict{Tuple,Float64} with 203 entries:\n  (83, 65, 66) => 0.014797\n  (46, 6, 60)  => 0.0185963\n  (25, 16, 23) => 1.0\n  (43, 7, 38)  => 0.995601\n  ⋮            => ⋮\n\njulia> Whale.get_triples(x, 68)\n3-element Array{Tuple{Int64,Int64,Int64},1}:\n (24, 75, 1)\n (2, 74, 22)\n (6, 36, 103)\n\n\n\n\n\n"
},

{
    "location": "api/#Whale.GBMModel",
    "page": "API",
    "title": "Whale.GBMModel",
    "category": "type",
    "text": "GBMModel(st::SlicedTree, ν::T, η::T, λ::T, μ::T, q::T) where\n    T::Union{<:Distribution,Array{<:Distribution,1},<:Real}\n\nHierarchical model for Whale using a GBM prior.\n\nν: autocorrelation strength (variance of geometric Brownian motion)\nη: parameter of geometric prior on the number of genes at the root of the species tree S\nλ: duplication rate at the root of S\nμ: loss rate at the root of S\nq: retention rates\n\nExample: InverseGamma prior on ν, Beta prior on η, Exponential priors on λ and μ at the root and Beta(1,1) priors on the retention rates q.\n\njulia> m = GBMModel(st, InverseGamma(15), Beta(10, 1), Exponential(1), Exponential(1), [Beta(1,1) for i=1:nwgd(st)]);\n\njulia> rand(m, st)  # sample a state from the prior\nDict{Symbol,Union{Float64, Array{Float64,1}}} with 5 entries:\n  :ν => 0.0771695\n  :μ => [0.498764, 0.500195, 0.459161, 0.615959, 0.482646, 0.502002, 0.402976, 0.354524, 0.378785, …\n  :λ => [0.760383, 0.733029, 0.696493, 0.592979, 0.718899, 1.01307, 0.642244, 0.639109, 0.600577, 0…\n  :η => 0.915222\n  :q => [0.761553, 0.729741, 0.116805, 0.278579, 0.200221, 0.293774, 0.869646]\n\n\n\n\n\n"
},

{
    "location": "api/#Whale.IRModel",
    "page": "API",
    "title": "Whale.IRModel",
    "category": "type",
    "text": "IRModel(st::SlicedTree, ν::T, η::T, λ::T, μ::T, q::T) where\n    T::Union{<:Distribution,Array{<:Distribution,1},<:Real}\n\nHierarchical model for Whale using a independent rates prior. As GBMModel but ν is the variance of the log-Normal prior on the branch-wise rates.\n\nExample:\n\njulia> m = IRModel(st, InverseGamma(10), Beta(6,3));\n\njulia> m.ν\nInverseGamma{Float64}(\ninvd: Gamma{Float64}(α=10.0, θ=1.0)\nθ: 1.0\n)\n\njulia> m.λ\nExponential{Float64}(θ=1.0)\n\njulia> rand(m, st)\nDict{Symbol,Union{Float64, Array{Float64,1}}} with 5 entries:\n  :ν => 0.0955299\n  :μ => [0.456988, 0.43639, 0.626645, 0.464031, 0.540543, 0.629841, 0.574996, 0.498506, 0.428841, 0…\n  :λ => [0.498921, 0.69414, 0.54769, 0.533595, 0.521085, 0.515748, 0.534503, 0.594045, 0.583829, 0.…\n  :η => 0.486934\n  :q => [0.970779, 0.965934, 0.788127, 0.405944, 0.827755, 0.253859, 0.395734]\n\n\n\n\n\n"
},

{
    "location": "api/#Whale.SlicedTree",
    "page": "API",
    "title": "Whale.SlicedTree",
    "category": "type",
    "text": "SlicedTree(tree::Arboreal, wgdconf=Dict(), Δt=0.05, minn=5, maxn=Inf)\n\nSliced species tree with optional WGD nodes. WGDs are specified using a configuration dict such as for example:\n\nwgdconf = Dict(\n    \"PPAT\" => (\"PPAT\", 0.655),\n    \"CPAP\" => (\"CPAP\", 0.275),\n    \"BETA\" => (\"ATHA\", 0.55),\n    \"ANGI\" => (\"ATRI,ATHA\", 3.08),\n    \"SEED\" => (\"GBIL,ATHA\", 3.9),\n    \"MONO\" => (\"OSAT\", 0.91),\n    \"ALPH\" => (\"ATHA\", 0.501))\n\nwhere the keys are WGD IDs (names), and values are given as a tuple or array with as first entry a single taxon or pair of taxa (comma-separated) that specifies the last common ancestor node that is preceded by the WGD of interest and the second entry the estimate time (before present) of the WGD event. Δt is the desired slice length (time interval), minn the minimum number of slices a branch should have, and maxn the maximum number of slices a branch should have.\n\n\n\n\n\n"
},

{
    "location": "api/#Whale.WhaleChain",
    "page": "API",
    "title": "Whale.WhaleChain",
    "category": "type",
    "text": "WhaleChain(st::SlicedTree, π::Model)\n\nChain object for performing MCMC under various hierarchical models defined in Whale.\n\njulia> w = WhaleChain(st, IRModel(st))\nWhaleChain{IRModel}(SlicedTree(9, 17, 7))\n\njulia> w[:q]  # current state\'s q values\n7-element Array{Float64,1}:\n 0.4068351791361286\n 0.5893350551453437\n 0.21472854312451684\n 0.8068526003250731\n 0.05239284527812649\n 0.8432325466244709\n 0.9006557436550706\n\njulia> w[:μ, 3]  # current state of μ for branch 3\n0.6578486200316909\n\njulia> logpdf(w)  # prior logpdf\n59.40417835345352\n\njulia> logpdf(w, :ν=>0.1)  # prior logpdf, but with ν = 0.1\n58.49660562613573\n\njulia> logpdf(w, :η=>0.2)  # prior logpdf, but with η = 0.2\n45.22412219063766\n\njulia> WhaleModel(w)  # WhaleModel corresponding to current state of the chain\nWhaleModel{Float64,CCD}(\nλ: [0.440648, 0.486472, 0.412741, 0.432659, 0.392268, 0.485483, 0.485181, 0.422687, 0.431686, 0.477667, 0.405392, 0.414238, 0.442341, 0.390893, 0.429505, 0.41379, 0.435038]\nμ: [0.635544, 0.575146, 0.657849, 0.657189, 0.626361, 0.668117, 0.619245, 0.600366, 0.682475, 0.570257, 0.58954, 0.576526, 0.6966, 0.729826, 0.561087, 0.614903, 0.668041]\nq: [0.406835, 0.589335, 0.214729, 0.806853, 0.0523928, 0.843233, 0.900656]\nη: 0.966691254252368\n\n\n\n\n\n"
},

{
    "location": "api/#Whale.WhaleModel",
    "page": "API",
    "title": "Whale.WhaleModel",
    "category": "type",
    "text": "WhaleModel{T<:Real,CCD}(S::SlicedTree, λ, μ, q, η=0.9, cond=\"oib\")\n\nThe Whale probabilistic model, containing both the sliced species tree, parameters and fields for extinction and propagation probabilities.\n\njulia> st = Whale.example_tree();\n\njulia> ccd = read_ale(\"example/example-ale/\", st);\n[ Info:  .. read 12 ALE files\n\njulia> w = WhaleModel(st, 0.2, 0.3, 0.1, 0.9)\nWhaleModel{Float64,CCD}(\nλ: [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]\nμ: [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3]\nq: [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]\nη: 0.9\n)\n\njulia> logpdf(w, ccd)\n-296.67809475568924\n\njulia> ccd[1].tmpmat  # this did not store the dynamic programming (DP) matrix\nDict{Int64,Array{Float64,2}} with 0 entries\n\njulia> logpdf(w, ccd, matrix=true)  # this will store the DP matrix\n-296.67809475568924\n\njulia> ccd[1].tmpmat  # this is the DP matrix, required for backtracking\nDict{Int64,Array{Float64,2}} with 24 entries:\n  1  => [0.000160506; 0.00276326; … ; 6.63525e-24; 4.18032e-15]\n  2  => [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]\n  16 => [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]\n  ...\n\n\n\n\n\n"
},

{
    "location": "api/#Whale.backtrack!",
    "page": "API",
    "title": "Whale.backtrack!",
    "category": "function",
    "text": "backtrack!(ccd::CCD, w::WhaleModel, n=1)\nbacktrack!(D::CCDArray, w::WhaleModel, n=1)\nbacktrack!(D::CCDArray, w::WhaleChain, n=1000)\nbacktrack!(D::CCDArray, df::DataFrame, st::SlicedTree, n=1000)\n\nSample reconciled trees (RecTree) from the Whale model conditional on the conditional clade distribution (CCD) by stochastic backtracking along the dynamic programming matrix as in Szollosi 2013. If the recmat field of the ccd argument is non-empty, this matrix will be used, otherwise the reconciliation matrix is computed from scratch. Reconciled trees are stored in the CCD object for convenience.\n\nnote: Note\nIf the second argument is a WhaleChain object (or data frame), than post hoc backtracking is performed, that is, reconciled trees are simulated from the posterior predictive distribution after running an MCMC chain.\n\n\n\n\n\n"
},

{
    "location": "api/#Whale.consensus",
    "page": "API",
    "title": "Whale.consensus",
    "category": "function",
    "text": "consensus(x::CCD, S::SlicedTree)\nconsensus(D::DArray, S::SlicedTree)\n\nGet consensus trees from the backtracked RecTrees in the rectrs field of the conditional clade distribution object(s) (CCD).\n\nwarning: Warning\nBranch lengths are not yet correct (they are meaningfully related to the true branch length estimates nevertheless).\n\nwarning: Warning\nA consensus reconciled tree does not include loss events. This is because different reconciled trees for the same family can have different numbers of loss events.\n\n\n\n\n\n"
},

{
    "location": "api/#Whale.contreetable-Tuple{IO,PhyloTrees.ConRecTree,SlicedTree}",
    "page": "API",
    "title": "Whale.contreetable",
    "category": "method",
    "text": "contreetable([io::IO,] contree::ConRecTree, S::SlicedTree)\n\nWrite a table representation of a consensus reconciled tree.\n\n\n\n\n\n"
},

{
    "location": "api/#Whale.describe-Tuple{WhaleModel}",
    "page": "API",
    "title": "Whale.describe",
    "category": "method",
    "text": "describe(w::WhaleModel)\n\nGet a detailed description of a particular WhaleModel.\n\n\n\n\n\n"
},

{
    "location": "api/#Whale.mcmc!-Tuple{WhaleChain,DistributedArrays.DArray{CCD,1,Array{CCD,1}},Int64,Vararg{Any,N} where N}",
    "page": "API",
    "title": "Whale.mcmc!",
    "category": "method",
    "text": "mcmc!(w::WhaleChain, D::CCDArray, n::Int64, args...; kwargs...)\n\nPerform n generations of MCMC sampling for a WhaleChain given a bunch of observed CCDs.\n\n\n\n\n\n"
},

{
    "location": "api/#Whale.mle",
    "page": "API",
    "title": "Whale.mle",
    "category": "function",
    "text": "mle(w::WhaleModel, ccd::CCDArray, optimizer=LBFGS(); kwargs...)\n\n\n\n\n\n"
},

{
    "location": "api/#Whale.read_ale-Tuple{String,SlicedTree}",
    "page": "API",
    "title": "Whale.read_ale",
    "category": "method",
    "text": "read_ale(fname::String, s::SlicedTree)\n\nRead in a bunch of conditional clade distributions (CCD) from ALEobserve (.ale) files. Either provide\n\na file with a filename on each line\na directory with .ale files\na single .ale file\nan empty file, for running MCMC under the prior alone\n\njulia> st = Whale.example_tree()\nSlicedTree(9, 17, 7)\n\njulia> ccd = read_ale(\"example/example-ale/\", st)\n[ Info:  .. read 12 ALE files\n12-element DistributedArrays.DArray{CCD,1,Array{CCD,1}}:\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 83 clades, 5001 samples)\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 55 clades, 5001 samples)\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 89 clades, 5001 samples)\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 131 clades, 5001 samples)\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 107 clades, 5001 samples)\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 59 clades, 5001 samples)\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 53 clades, 5001 samples)\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 83 clades, 5001 samples)\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 59 clades, 5001 samples)\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 95 clades, 5001 samples)\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 67 clades, 5001 samples)\n CCD{Float64,PhyloTrees.RecTree}(13 taxa, 65 clades, 5001 samples)\n\n\n\n\n\n"
},

{
    "location": "api/#Whale.ConstantDistribution",
    "page": "API",
    "title": "Whale.ConstantDistribution",
    "category": "type",
    "text": "ConstantDistribution(x)\n\nA \'constant\' distribution (Dirac mass), sometimes useful.\n\n\n\n\n\n"
},

{
    "location": "api/#Whale.GeometricBrownianMotion",
    "page": "API",
    "title": "Whale.GeometricBrownianMotion",
    "category": "type",
    "text": "GeometricBrownianMotion{T<:Real}\nGBM(t::SlicedTree, r::T, v::T) where {T<:Real}\n\nDistribution induced by a Geometric Brownian Motion (GBM) over the SlicedTree. r is the value at the root of the tree, ν is the variance (autocorrelation strength).\n\nThe log density for the GBM distribution is computed with an implementation of the GBM prior on rates based on Ziheng Yang\'s MCMCTree, described in Rannala & Yang (2007). This uses the approach whereby rates are defined for midpoints of branches, and where a correction is performed to ensure that the correlation is proper (in contrast with Thorne et al. 1998).\n\n\n\n\n\n"
},

{
    "location": "api/#API-1",
    "page": "API",
    "title": "API",
    "category": "section",
    "text": "Modules = [Whale]"
},

]}
