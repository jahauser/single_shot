# jobs/collate.jl
#
# Usage: julia jobs/collate.jl <output_dir>
#
# Collates all .jld2 files in <output_dir> into a single collated_results.jld2.
# Files sharing the same (tag, algo, L, T, p, q, c) key are merged:
#   - (1) M1s, M2s, E1s, E2s: weighted-averaged time series
#   - (2) Ms_snapshots, Es_snapshots: concatenated across files (rows = samples)
#   - (3) instant_failures: weighted-averaged time series
#   - (4) cumulative_failures: weighted-averaged time series

using JLD2
using Glob

if length(ARGS) != 1
    println("Usage: julia jobs/collate.jl <output_dir>")
    exit(1)
end

outdir = ARGS[1]
files  = glob("*.jld2", outdir)
filter!(f -> !endswith(f, "collated_results.jld2"), files)

if isempty(files)
    println("No .jld2 files found in ", outdir)
    exit(0)
end

results = Dict{NTuple{7,Any}, NamedTuple}()

for f in files
    @load f algo L T p q c samples tag M1s M2s E1s E2s instant_failures cumulative_failures Ms_snapshots Es_snapshots dt
    key = (tag, algo, L, T, p, q, c)

    if !haskey(results, key)
        results[key] = (
            sum_M1              = samples .* copy(M1s),
            sum_M2              = samples .* copy(M2s),
            sum_E1              = samples .* copy(E1s),
            sum_E2              = samples .* copy(E2s),
            sum_instant_fail    = samples .* instant_failures,
            sum_cumulative_fail = samples .* cumulative_failures,
            Ms_snapshots        = copy(Ms_snapshots),
            Es_snapshots        = copy(Es_snapshots),
            total_samples       = samples,
            total_time          = float(dt),
        )
    else
        r = results[key]
        results[key] = (
            sum_M1              = r.sum_M1 .+ samples .* M1s,
            sum_M2              = r.sum_M2 .+ samples .* M2s,
            sum_E1              = r.sum_E1 .+ samples .* E1s,
            sum_E2              = r.sum_E2 .+ samples .* E2s,
            sum_instant_fail    = r.sum_instant_fail .+ samples .* instant_failures,
            sum_cumulative_fail = r.sum_cumulative_fail .+ samples .* cumulative_failures,
            Ms_snapshots        = vcat(r.Ms_snapshots, Ms_snapshots),
            Es_snapshots        = vcat(r.Es_snapshots, Es_snapshots),
            total_samples       = r.total_samples + samples,
            total_time          = r.total_time + dt,
        )
    end
end

# Normalise weighted sums to weighted averages
for (key, r) in results
    results[key] = (
        M1s                 = r.sum_M1 ./ r.total_samples,
        M2s                 = r.sum_M2 ./ r.total_samples,
        E1s                 = r.sum_E1 ./ r.total_samples,
        E2s                 = r.sum_E2 ./ r.total_samples,
        instant_failures    = r.sum_instant_fail ./ r.total_samples,
        cumulative_failures = r.sum_cumulative_fail ./ r.total_samples,
        Ms_snapshots        = r.Ms_snapshots,
        Es_snapshots        = r.Es_snapshots,
        total_samples       = r.total_samples,
        total_time          = r.total_time,
    )
end

outfile = joinpath(outdir, "collated_results.jld2")
@save outfile results
println("Collated $(length(results)) parameter sets from $(length(files)) files -> ", outfile)