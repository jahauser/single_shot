#!/usr/bin/env julia
using ArgParse
using Dates
using JLD2
using Printf

using StatsBase
using Graphs
using SimpleWeightedGraphs
using BlossomV
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
            help = "MV | TOOM | basic_MWPM | adv_MWPM"
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
        "--peff"
            help = "adv_MWPM only"
            arg_type = Float64
            default = NaN
        "--samples"
            arg_type = Int
            default = 1
    end
    return s
end

f3(x) = replace(@sprintf("%.3f", x), "." => "p")

function main(args)
    o = parse_args(args, build_parser())
    algo, L, T = o["algo"], o["L"], o["T"]
    p, q, peff, samples = o["p"], o["q"], o["peff"], o["samples"]

    # Hardcoded output dir
    outdir = joinpath(ROOT, "output")
    mkpath(outdir)

    t0 = time()

    if algo == "MV"
        _, _ = MV_sample(2, 2, 0.0, 0.0, 1) 
        M1s, M2s = MV_sample(L, T, p, q, samples)
    elseif algo == "TOOM"
        _, _ = TOOM_sample(2, 2, 0.0, 0.0, 1)
        M1s, M2s = TOOM_sample(L, T, p, q, samples)
    elseif algo == "basic_MWPM"
        _, _ = basic_MWPM_sample(2, 2, 0.0, 0.1, 1)
        M1s, M2s = basic_MWPM_sample(L, T, p, q, samples)
    elseif algo == "adv_MWPM"
        if isnan(peff)
            error("adv_MWPM requires --peff")
        end
        _, _ = adv_MWPM_sample(2, 2, 0.0, 0.5, 0.1, 1)
        M1s, M2s = adv_MWPM_sample(L, T, p, peff, q, samples)
    else
        error("Unknown --algo=$(algo)")
    end

    dt = time() - t0


    tag = "algo=$(algo)_L=$(L)_T=$(T)_p=$(f3(p))_q=$(f3(q))"
    if algo == "adv_MWPM"; tag *= "_peff=$(f3(peff))"; end
    tag *= "_samples=$(samples)"
    ts = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    fname = joinpath(outdir, "sample_$(tag)_$(ts)_$(rand(UInt32)).jld2")

    @save fname algo L T p q peff samples M1s M2s dt
    println("Wrote: ", fname)
end

main(ARGS)
