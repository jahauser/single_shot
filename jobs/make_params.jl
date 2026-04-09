#!/usr/bin/env julia
# jobs/make_params.jl
#
# Generates jobs/params.txt where each line is a CLI for jobs/run.jl

const PARAMS_FILE = joinpath(@__DIR__, "params.txt")

"Format one CLI line for run.jl."
function format_line(algo::AbstractString, L::Int, T::Int, p::Float64, q::Float64; c::Float64=NaN, samples::Int=1, tag::AbstractString="")
    base = "--algo $algo --L $L --T $T --p $p --q $q --samples $samples"
    base = algo == "IsingML" ? "$base --c $c" : base
    return isempty(tag) ? base : "$base --tag $tag"
end

"Write or append parameter lines to jobs/params.txt."
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

algos   = ["IsingML"]
tag     = ""           # label for this batch of jobs, e.g. "sweep_Mar25"

Ls      = 8:8:24
T_of_L  = L -> 10L

ps      = 0.0:0.02:0.5
qs      = 0.01:0.05:0.3

# c = (peff - q) / (0.5 - q), so the four sweep points simplify to:
#   peff = q                  -> c = 0
#   peff = q + 2p(1-p)(1-2q) -> c = 4p(1-p)   [the (0.5-q) cancels]
#   peff = 0.4                -> c = (0.4-q)/(0.5-q)
#   peff = 0.5                -> c = 1
cs(p, q) = [0.0, 4p*(1-p), (0.4-q)/(0.5-q), 1.0]

# per-L sampling and repeats
samples = Dict(
    8  => 10000,
    16 => 10000,
    24 => 1000,
)
repeats = Dict(
    8  => 1,
    16 => 1,
    24 => 10,
)

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
        if algo == "IsingML"
            for c in cs(p, q)
                push!(lines, format_line(algo, L, T, p, q; c=c, samples=S, tag=tag))
            end
        else
            push!(lines, format_line(algo, L, T, p, q; samples=S, tag=tag))
        end
    end
end

write_params(lines; append=append)