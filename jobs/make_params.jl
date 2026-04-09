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

algos   = ["MV"]
tag     = "MV_L80_sweep"

Ls      = [80]
T_of_L  = L -> 20L

ps      = 0.0:0.01:0.5
qs      = 0.0:0.01:0.5

samples = Dict(80 => 100000)
repeats = Dict(80 => 1)

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
        push!(lines, format_line(algo, L, T, p, q; samples=S, tag=tag))
    end
end

write_params(lines; append=append)