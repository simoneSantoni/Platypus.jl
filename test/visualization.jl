@testset "Plot labels, overlays and numerical color keys" begin
    result = OutcomeResult([1.5],[:intercept],[1.0,2.0],[1.5,1.5],[-0.5,0.5],
                           0.0,[1,2],ones(2,1),(;))
    spec = plot_outcome_fit(result)
    recipe = RecipesBase.apply_recipe(Dict{Symbol,Any}(),spec)
    @test length(recipe) == 2
    @test recipe[1].plotattributes[:xlabel] == "Observed"
    @test recipe[1].plotattributes[:ylabel] == "Fitted"
    @test recipe[1].plotattributes[:title] == "Descriptive outcome fit"
    # Plots applies subplot defaults per series: even omitted overlay labels
    # can erase the first series' labels. Both series carry the same labels.
    @test recipe[2].plotattributes[:xlabel] == "Observed"
    @test recipe[2].plotattributes[:ylabel] == "Fitted"
    @test recipe[2].plotattributes[:title] == "Descriptive outcome fit"
    custom = RecipesBase.apply_recipe(Dict{Symbol,Any}(:xlabel=>"Measured"),spec)
    @test all(r.plotattributes[:xlabel] == "Measured" for r in custom)
    result.observed[1] = 99
    @test recipe[1].args[1] == [1.0,2.0]

    es = EntitySet([0.0 2.0;1.0 0.0;2.0 1.0],[:x,:y])
    fit = fit_typicality(build_pairs(es))
    for heatmap in (plot_correlation_heatmap(es),plot_pair_heatmap(fit))
        recipe = RecipesBase.apply_recipe(Dict{Symbol,Any}(),heatmap)
        @test only(recipe).plotattributes[:colorbar] === true
        hidden = RecipesBase.apply_recipe(Dict{Symbol,Any}(:colorbar=>false),heatmap)
        @test only(hidden).plotattributes[:colorbar] === false
    end
end
