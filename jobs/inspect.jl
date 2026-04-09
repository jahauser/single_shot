#!/usr/bin/env julia
# jobs/inspect.jl
#
# Usage: julia jobs/inspect.jl <path/to/file.jld2>

using JLD2
using Printf
using Statistics

if isempty(ARGS)
    println("Usage: julia jobs/inspect.jl <file.jld2> [file2.jld2 ...]")
    println("       julia jobs/inspect.jl output/test/*.jld2")
    exit(1)
end

for fname in ARGS
    @load fname algo L T p q c peff samples tag M1s M2s E1s E2s instant_failures cumulative_failures Ms_hist Es_hist dt

    println("File:     ", fname)
    println("Tag:      ", isempty(tag) ? "(none)" : tag)
    println()
    println("Algo:     ", algo)
    println("L=$(L)  T=$(T)  p=$(p)  q=$(q)  c=$(c)  peff=$(@sprintf("%.4f", peff))")
    println("Samples:  ", samples)
    println("Run time: $(@sprintf("%.2f", dt))s")
    println()
    println("Final values (t=T):")
    println("  M1  (mean magnetization):    $(@sprintf("%.6f", M1s[end]))")
    println("  M2  (mean mag²):             $(@sprintf("%.6f", M2s[end]))")
    println("  E1  (mean energy):           $(@sprintf("%.6f", E1s[end]))")
    println("  E2  (mean energy²):          $(@sprintf("%.6f", E2s[end]))")
    println("  instant failure rate:        $(@sprintf("%.6f", instant_failures[end]))")
    println("  cumulative failure rate:     $(@sprintf("%.6f", cumulative_failures[end]))")
    println()
    # Reconstruct bin edges from L and report mode of final snapshot histogram
    M_edges = (0:L^2) ./ L^2
    E_edges = (-2*L*(L-1):2:2*L*(L-1))
    M_mode  = M_edges[argmax(Ms_hist[:, end])]
    E_mode  = E_edges[argmax(Es_hist[:, end])]
    println("Snapshot histograms (t=T, last snapshot):")
    println("  M  mode: $(@sprintf("%.4f", M_mode))  ($(sum(Ms_hist[:,end])) counts)")
    println("  E  mode: $(@sprintf("%.1f", Float64(E_mode)))  ($(sum(Es_hist[:,end])) counts)")
    println("-"^60)
end