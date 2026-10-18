# Optional integration validation: Plots is deliberately not a package dependency.
# Run from an environment containing Platypus and Plots. The optional first
# argument is an output directory; otherwise a retained temporary one is created.
using Platypus
using Plots
using Random
using Test

function validate_rendering(output_dir)
    mkpath(output_dir)
    gr()
    rng = MersenneTwister(9)
    es = EntitySet(rand(rng,30,4),[:a,:b,:c,:d])
    ps = build_pairs(es)
    fitted = fit_typicality(ps)
    outcome = fit_outcome(fitted,rand(rng,30))
    graph = fit_graph_typicality(es)
    specifications = [
        ("features",plot_feature_distributions(es)),
        ("cloud",plot_pair_cloud(ps,1)),
        ("correlations",plot_correlation_heatmap(es)),
        ("redundancy",plot_pair_redundancy(ps)),
        ("profile",plot_typicality_profile(fitted,1)),
        ("atypicality",plot_atypicality_distribution(fitted)),
        ("graph",plot_atypicality_distribution(graph)),
        ("novelty",plot_novelty_conventionality(fitted)),
        ("outcome",plot_outcome_fit(outcome)),
        ("heatmap",plot_pair_heatmap(fitted)),
    ]
    @testset "Plots/GR backend rendering" begin
        for (name,specification) in specifications
            figure = plot(specification)
            if name == "outcome"
                @test figure[1][:xaxis][:guide] == "Observed"
                @test figure[1][:yaxis][:guide] == "Fitted"
                @test figure[1][:title] == "Descriptive outcome fit"
            elseif name in ("correlations","heatmap")
                @test !(figure[1][:colorbar] in (false,:none))
            elseif name == "features"
                for (i,feature) in enumerate(feature_names(es))
                    @test figure[i][:xaxis][:guide] == string(feature)
                    @test figure[i][:yaxis][:guide] == "Count"
                end
            end
            path = joinpath(output_dir,name * ".svg")
            savefig(figure,path)
            @test filesize(path) > 0
            println("Rendered ",path)
        end
    end
    println("Rendered all nine plotting helpers and the graph-result overload.")
    return output_dir
end

output_dir = isempty(ARGS) ? mktempdir(;cleanup=false) : abspath(only(ARGS))
validate_rendering(output_dir)
