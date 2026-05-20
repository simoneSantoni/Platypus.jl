"""
Entity loading and parsing utilities.

Entities are accepted from CSV files or in-memory `DataFrame`s. Feature
columns may be named explicitly or detected automatically (every numeric
column that is not the id or label column). Entity identifiers and an
optional category label are carried through unchanged.
"""

"""
    load_entities(filepath::String; kwargs...) -> EntitySet

Load entities from a CSV file. The CSV is read with
`CSV.read(filepath, DataFrame)` and then dispatched to the `DataFrame`
method below; see that method for the keyword arguments.
"""
function load_entities(filepath::AbstractString;
                       feature_cols::Union{AbstractVector, Nothing}=nothing,
                       id_col::Union{Symbol, Nothing}=:id,
                       label_col::Union{Symbol, Nothing}=nothing)
    df = CSV.read(filepath, DataFrame)
    load_entities(df; feature_cols, id_col, label_col)
end

"""
    load_entities(df::DataFrame; kwargs...) -> EntitySet

Build an [`EntitySet`](@ref) from a `DataFrame` of entity records.

# Keyword Arguments
- `feature_cols::Union{AbstractVector, Nothing}=nothing`: Columns holding
  feature values. If `nothing`, every numeric column other than `id_col`
  and `label_col` is used, in their order of appearance.
- `id_col::Union{Symbol, Nothing}=:id`: Column with integer entity IDs.
  If the column is absent (or `nothing`), IDs default to `1:nrow(df)`.
- `label_col::Union{Symbol, Nothing}=nothing`: Optional column with a
  category label per entity (stored as `String`, never used to define
  the reference set).
"""
function load_entities(df::DataFrame;
                       feature_cols::Union{AbstractVector, Nothing}=nothing,
                       id_col::Union{Symbol, Nothing}=:id,
                       label_col::Union{Symbol, Nothing}=nothing)
    cols = propertynames(df)

    fcols = if feature_cols !== nothing
        Symbol.(feature_cols)
    else
        excl = Symbol[]
        id_col !== nothing && id_col in cols && push!(excl, id_col)
        label_col !== nothing && label_col in cols && push!(excl, label_col)
        [c for c in cols if c ∉ excl && eltype(df[!, c]) <: Union{Missing, Real}]
    end
    isempty(fcols) &&
        throw(ArgumentError("no feature columns were found or specified"))

    n = nrow(df)
    ids = (id_col !== nothing && id_col in cols) ? Int.(df[!, id_col]) : collect(1:n)
    labels = (label_col !== nothing && label_col in cols) ?
        df[!, label_col] : fill(nothing, n)

    entities = Entity[]
    sizehint!(entities, n)
    for r in 1:n
        feats = Float64[parse_feature(df[r, c]) for c in fcols]
        lab = labels[r]
        labstr = (lab === nothing || lab isa Missing) ? nothing : String(string(lab))
        push!(entities, Entity(ids[r], feats, labstr))
    end
    return EntitySet(entities, fcols)
end

"""
    parse_feature(val) -> Float64

Coerce a single feature cell to `Float64`. Numbers are converted
directly; strings are parsed. `missing` is rejected — impute or drop
incomplete rows before loading, because the density estimators in
[`fit_typicality`](@ref) have no missing-value semantics.
"""
parse_feature(val::Real) = Float64(val)
parse_feature(val::AbstractString) = parse(Float64, val)
parse_feature(::Missing) = throw(ArgumentError(
    "missing feature values are not supported; impute or drop rows before loading"))
