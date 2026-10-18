@testset "Research workflow: feature selection, provenance and cohorts" begin
    rng = MersenneTwister(412)
    df = DataFrame(id=1:24, year=repeat(2020:2022;inner=8),
        f1=rand(rng,24), f2=rand(rng,24), outcome=rand(rng,24))
    selected = [:f1,:f2]
    es = load_entities(df;feature_cols=selected)
    @test feature_names(es) == selected
    @test n_features(es) == 2
    # Changing an outcome cannot change explicitly selected measurement data.
    changed = copy(df)
    changed.outcome .= 1000
    @test load_entities(changed;feature_cols=selected).matrix == es.matrix
    # Automatic numeric discovery is convenient, but not a leakage safeguard.
    @test :outcome in feature_names(load_entities(df))
    @test :year in feature_names(load_entities(df))

    provenance = (dataset="fixture-v1",cohort="year <= 2021",features=copy(selected),
        preprocessing="complete-case",cutoff=2021)
    past = df[df.year .<= 2021,:]
    fitted = fit_typicality(build_pairs(load_entities(past;feature_cols=selected));
        method=:kde,provenance=provenance)
    metadata = analysis_metadata(fitted)
    @test metadata.provenance.cutoff == 2021
    @test metadata.provenance.dataset == "fixture-v1"
    @test metadata.transform === :rank
    @test metadata.feature_names == selected
    provenance.features[1] = :mutated
    metadata.provenance.features[2] = :mutated
    @test analysis_metadata(fitted).provenance.features == selected
    # Future records are excluded by cohort construction, including marginal
    # preprocessing. This checks the protocol, not a historical prediction API.
    future_changed = copy(df)
    future_changed[future_changed.year .> 2021,:f1] .= 1e6
    revised_past = future_changed[future_changed.year .<= 2021,:]
    refit = fit_typicality(build_pairs(load_entities(revised_past;feature_cols=selected));
        method=:kde,provenance=(dataset="fixture-v1",cohort="year <= 2021"))
    @test refit.pair_scores == fitted.pair_scores
    @test fitted.entity_ids == past.id

    incomplete = copy(df)
    allowmissing!(incomplete,:f1)
    incomplete[1,:f1] = missing
    @test_throws ArgumentError load_entities(incomplete;feature_cols=selected)
    complete = dropmissing(incomplete,selected)
    @test n_entities(load_entities(complete;feature_cols=selected)) == 23
    @test !(1 in load_entities(complete;feature_cols=selected).entity_ids)
end
