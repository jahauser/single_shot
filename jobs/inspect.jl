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
    @load fname algo L T p q c peff samples tag M1s M2s E1s E2s instant_failures cumulative_failures Ms_snapshots Es_snapshots dt

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
    println("Snapshot distributions (t=T, i.e. last snapshot):")
    println("  M  min/median/max:  $(@sprintf("%.4f", minimum(Ms_snapshots[:,end]))) / $(@sprintf("%.4f", median(Ms_snapshots[:,end]))) / $(@sprintf("%.4f", maximum(Ms_snapshots[:,end])))")
    println("  E  min/median/max:  $(@sprintf("%.4f", minimum(Es_snapshots[:,end]))) / $(@sprintf("%.4f", median(Es_snapshots[:,end]))) / $(@sprintf("%.4f", maximum(Es_snapshots[:,end])))")
    println("-"^60)
end