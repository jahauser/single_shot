# jobs/collate.jl
using JLD2
using Glob

outdir = joinpath(@__DIR__, "..", "output")
files = glob("*.jld2", outdir)

if isempty(files)
    println("No .jld2 files found in ", outdir)
    exit(0)
end

# results[(L,T,lambda,delta,q)] = (sum_E, sum_E2, total_samples, obs, total_time)
results = Dict{NTuple{6,Any}, Tuple{Vector, Vector, Int, Float64}}()

for f in files
    @load f algo L T p q peff samples M1s M2s dt
    key = (algo, L, T, p, q, peff)

    # Convert per-file averages to weighted sums
    # sum_E   += samples * E[x]
    # sum_E2  += samples * E[x^2]
    if !haskey(results, key)
        sum_E  = samples .* copy(M1s)
        sum_E2 = samples .* copy(M2s)
        results[key] = (sum_E, sum_E2, samples, float(dt))
    else
        sum_E, sum_E2, acc_samples, acc_dt = results[key]
        sum_E  .+= samples .* M1s
        sum_E2 .+= samples .* M2s
        results[key] = (sum_E, sum_E2, acc_samples + samples, acc_dt + dt)
    end
end

# Normalize to weighted averages: E[x] and E[x^2]
for (key, (sum_E, sum_E2, total_samples, total_time)) in results
    sum_E  ./= total_samples
    sum_E2 ./= total_samples
    # Overwrite with normalized values
    results[key] = (sum_E, sum_E2, total_samples, total_time)
end

outfile = joinpath(outdir, "collated_results.jld2")
@save outfile results
println("Collated $(length(results)) parameter sets from $(length(files)) files -> ", outfile)
