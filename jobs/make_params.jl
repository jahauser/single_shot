#!/usr/bin/env julia
# jobs/make_params.jl
#
# Generates jobs/params.txt where each line is a CLI for jobs/run.jl
const PARAMS_FILE = joinpath(@__DIR__, "params.txt")

"Format one CLI line for run.jl."
function format_line(algo::AbstractString, L::Int, T::Int, p::Float64, q::Float64;
                     c::Float64=NaN, samples::Int=1, tag::AbstractString="")
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

function main()
    algo   = "IsingML"
    tag    = "IsingML_fineq_sweep"
    T_of_L = L -> 20L
    append = false

    samples = Dict(8 => 1000, 12 => 1000, 16 => 100, 20 => 100, 24 => 100, 28 => 10, 32 => 10)
    repeats = Dict(8 => 10,    12 => 10,    16 => 10,   20 => 10,   24 => 10,   28 => 100,  32 => 100)

    # Each block: fixed p, fixed c, specific q values to probe
    blocks = [
        (p=0.05, c=0.05, qs=[0.186, 0.188, 0.192, 0.194, 0.196, 0.198]),
        (p=0.30, c=0.80, qs=[0.094, 0.096, 0.098, 0.102, 0.104, 0.106]),
    ]

    lines = String[]

    for L in [8, 12, 16, 20, 24, 28, 32]
        T = T_of_L(L)
        S = samples[L]
        R = repeats[L]
        for blk in blocks
            for _ in 1:R, q in blk.qs
                push!(lines, format_line(algo, L, T, Float64(blk.p), Float64(q);
                                         c=Float64(blk.c), samples=S, tag=tag))
            end
        end
    end

    println("Total jobs: $(length(lines))")
    write_params(lines; append=append)
end

main()