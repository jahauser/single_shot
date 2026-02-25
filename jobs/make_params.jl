#!/usr/bin/env julia
# scripts/gen_params.jl
#
# Generates scripts/params.txt where each line is a CLI for scripts/run.jl

const PARAMS_FILE = joinpath(@__DIR__, "params.txt")

"Format one CLI line for run.jl."
function format_line(algo::AbstractString, L::Int, T::Int, p::Float64, q::Float64; peff::Float64=NaN, samples::Int=1)
    base = "--algo $algo --L $L --T $T --p $p --q $q --samples $samples"
    return algo == "adv_MWPM" ? "$base --peff $peff" : base
end

"Write or append parameter lines to scripts/params.txt."
function write_params(lines::Vector{String}; append::Bool=false)
    open(PARAMS_FILE, append ? "a" : "w") do io
        for ln in lines
            println(io, ln)
        end
    end
    action = append ? "Appended" : "Wrote"
    println("$action $(length(lines)) lines to $(PARAMS_FILE)")
end

# -------------------------------
# Define your sweeps here
# -------------------------------

# algos    = ["MV", "basic_MWPM", "adv_MWPM"]

# Ls       = collect(4:4:16)                  # example: 8,16,24,32,40
# T_of_L   = L -> 10L                      # example mapping; edit as needed

# ps       = 0.0:0.05:0.5
# qs       = 0.0:0.05:0.3
# peffs(q) = collect(q:0.05:0.5)

# # per-L sampling and repeats (edit as needed)
# samples  = Dict(
#     4  => 100000,
#     8 => 10000,
#     12 => 1000,
#     16 => 100,
# )
# repeats  = Dict(
#     4 => 1,
#     8  => 1,
#     12 => 1,
#     16 => 10,
# )

algos    = ["adv_MWPM"]

Ls       = 8:8:24           # example: 8,16,24,32,40
T_of_L   = L -> 10L                      # example mapping; edit as needed

ps       = 0.0:0.02:0.5
qs       = 0.0:0.05:0.3
peffs(p,q) = [q, q+2p*(1-p)*(1-2q), 0.4, 0.5]

# per-L sampling and repeats (edit as needed)
samples  = Dict(
    # 4  => 100000,
    8 => 10000,
    # 12 => 1000,
    16 => 10000,
    # 20 => 100,
    24 => 1000,
)
repeats  = Dict(
    # 4 => 1,
    8  => 1,
    # 12 => 2,
    16 => 1,
    # 20 => 20,
    24 => 10,
)

# algos    = ["MV"]

# Ls       = 8:8:104           # example: 8,16,24,32,40
# T_of_L   = L -> L                      # example mapping; edit as needed

# ps       = 0.0:0.01:0.3
# qs       = 0.0:0.01:0.3
# # peffs(p,q) = [q+2p*(1-p)]
# samples = 10000

append = false  # set true to append to existing params.txt

# -------------------------------
# Build lines
# -------------------------------
lines = String[]
for L in Ls
    T = T_of_L(L)
    S = samples[L]
    R = repeats[L]
    for _ in 1:R, p in ps, q in qs, algo in algos
        if algo == "adv_MWPM"
            for peff in peffs(p,q)
                push!(lines, format_line(algo, L, T, p, q; peff=peff, samples=S))
            end
        else
            push!(lines, format_line(algo, L, T, p, q; samples=S))
        end
    end
end

write_params(lines; append=append)
