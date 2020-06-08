using Pkg; Pkg.activate(@__DIR__)
using Whale, NewickTree
using Test
using Random

@testset "Whale tests" begin
    @testset "likelihood" begin
        data = joinpath(@__DIR__, "../example/example-1/ale")
        t = deepcopy(Whale.extree)
        n = length(postwalk(t))
        l = length(getleaves(t))
        insertnode!(getlca(t, "ATHA", "ATHA"), name="wgd_1")
        insertnode!(getlca(t, "ATHA", "ATRI"), name="wgd_2")
        r = Whale.RatesModel(
            DLWGD(λ=ones(n), μ=ones(n), q=[0.2, 0.1], η=0.9, p=zeros(l)),
            fixed=(:p,))
        w = WhaleModel(r, t, 0.05)
        ccd = read_ale(data, w)
        # @test logpdf!(w, ccd) ≈ -588.0647568294178
        # @test logpdf(w, ccd) ≈ -588.0647568294178
        @test logpdf!(w, ccd) ≈ -570.9667405899105
        @test logpdf(w, ccd) ≈ -570.9667405899105
        ccd = read_ale(data, w, true)
        # @test logpdf!(w, ccd) ≈ -588.0647568294178
        # @test logpdf(w, ccd) ≈ -588.0647568294178
        @test logpdf!(w, ccd) ≈ -570.9667405899105
        @test logpdf(w, ccd) ≈ -570.9667405899105

        # root rates do not have an effect
        r.params.λ[id(getroot(w))] = NaN
        r.params.μ[id(getroot(w))] = NaN
        w = WhaleModel(r, t, 0.05)
        ccd = read_ale(data, w)
        @test logpdf!(w, ccd) ≈ -570.9667405899105
    end

    include(joinpath(@__DIR__, "mle.jl"))
end
