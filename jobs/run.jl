#!/usr/bin/env julia
using ArgParse
using Dates
using JLD2
using Printf

using StatsBase
using Graphs
using SimpleWeightedGraphs
using GraphsMatching

const ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(ROOT, "src", "circuits.jl"))
include(joinpath(ROOT, "src", "local.jl"))
include(joinpath(ROOT, "src", "main.jl"))
include(joinpath(ROOT, "src", "mwpm.jl"))

function build_parser()
    s = ArgParseSettings(; description = "Run one simulation and save .jld2")
    @add_arg_table s begin
        "--algo"
            help     = "MV | TOOM | UV | IsingML"
            arg_type = String
            required = true
        "--L"
            arg_type = Int; required = true
        "--T"
            arg_type = Int; required = true
        "--p"
            arg_type = Float64; required = true
        "--q"
            arg_type = Float64; required = true
        "--c"
            help     = "IsingML only: correlation parameter; peff = q + c*(0.5 - q)"
            arg_type = Float64
            default  = NaN
        "--samples"
            arg_type = Int
            default  = 1
        "--tag"
            help     = "optional label to group/categorise related jobs"
            arg_type = String
            default  = ""
    end
    return s
end

f3(x) = replace(@sprintf("%.3f", x), "." => "p")

function make_decoder(algo, q, c)
    if algo == "MV"
        return MajorityVote()
    elseif algo == "TOOM"
        return Toom()
    elseif algo == "UV"
        return UnanimousVote()
    elseif algo == "IsingML"
        isnan(c) && error("IsingML requires --c")
        return IsingML(q, c)
    else
        error("Unknown --algo=$(algo)")
    end
end

function main(args)
    o       = parse_args(args, build_parser())
    algo    = o["algo"]
    L, T    = o["L"], o["T"]
    p, q, c = o["p"], o["q"], o["c"]
    samples = o["samples"]
    tag     = o["tag"]

    outdir = joinpath(ROOT, "output", isempty(tag) ? "untagged" : tag)
    mkpath(outdir)

    decoder = make_decoder(algo, q, c)

    # Warmup to avoid timing JIT compilation
    sample(2, 2, 0.0, 0.0, decoder, 1)

    t0 = time()
    M1s, M2s, E1s, E2s, instant_failures, cumulative_failures, Ms_hist, Es_hist = sample(L, T, p, q, decoder, samples)
    dt = time() - t0

    peff = algo == "IsingML" ? q + c * (0.5 - q) : NaN

    fname_tag  = "algo=$(algo)_L=$(L)_T=$(T)_p=$(f3(p))_q=$(f3(q))"
    fname_tag *= algo == "IsingML" ? "_c=$(f3(c))" : ""
    fname_tag *= "_samples=$(samples)"
    ts    = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    fname = joinpath(outdir, "sample_$(fname_tag)_$(ts)_$(rand(UInt32)).jld2")

    @save fname algo L T p q c peff samples tag M1s M2s E1s E2s instant_failures cumulative_failures Ms_hist Es_hist dt
    println("Wrote: ", fname)
end

main(ARGS)