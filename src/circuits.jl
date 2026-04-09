abstract type AbstractDecoder end

struct MajorityVote <: AbstractDecoder end
struct Toom <: AbstractDecoder end
struct UnanimousVote <: AbstractDecoder end

struct IsingML <: AbstractDecoder
    q::Float64
    c::Float64
end

function correct(ρ::AbstractMatrix, checks::Tuple, ::MajorityVote)
    return majority_vote(ρ, checks)
end

function correct(ρ::AbstractMatrix, checks::Tuple, ::Toom)
    return toom(ρ, checks)
end

function correct(ρ::AbstractMatrix, checks::Tuple, ::UnanimousVote)
    return unanimous_vote(ρ, checks)
end

function correct(ρ::AbstractMatrix, checks::Tuple, decoder::IsingML)
    q, c = decoder.q, decoder.c
    peff = q + c * (1/2 - q)
    horizontal_checks, vertical_checks = heal(checks, peff, q)
    domain = track_domains((horizontal_checks, vertical_checks))
    if magnetization(domain) > 0.5
        domain = .!domain
    end
    
    return ρ .⊻ domain
end

function sample(L::Int, T::Int, p::Float64, q::Float64, decoder::AbstractDecoder)
    ρ = initialstate(L)
    Ms = zeros(Float64, T+1)
    Es = zeros(Float64, T+1)


    Ms[1] = magnetization(ρ)
    Es[1] = energy(ρ)

    for t in 1:T
        ρ = noiselayer(ρ, p)
        checks = measure(ρ, q)
        ρ = correct(ρ, checks, decoder)
        Ms[t+1] = magnetization(ρ)
        Es[t+1] = energy(ρ)
    end
    return Ms, Es
end

function sample(L::Int, T::Int, p::Float64, q::Float64, decoder::AbstractDecoder, samples::Int)
    # snapshot indices: t = 0, L, 2L, ..., T  (i.e. every L steps)
    snapshot_ts  = 0:L:T
    n_snapshots  = length(snapshot_ts)

    # histogram bin counts:
    #   M = k/L² for k in 0:L²  =>  n_M_bins = L²+1
    #   E is even integer in [-2L(L-1), 2L(L-1)]  =>  n_E_bins = 2L(L-1)+1
    # bin index helpers (1-based):
    #   M_bin(M) = round(Int, M * L²) + 1
    #   E_bin(E) = round(Int, E) ÷ 2 + L*(L-1) + 1
    n_M_bins = L^2 + 1
    n_E_bins = 2*L*(L-1) + 1

    # (1) averaged time series
    M1s = zeros(Float64, T+1)
    M2s = zeros(Float64, T+1)
    E1s = zeros(Float64, T+1)
    E2s = zeros(Float64, T+1)

    # (2) histograms at each snapshot: shape (n_bins × n_snapshots)
    Ms_hist = zeros(Int32, n_M_bins, n_snapshots)
    Es_hist = zeros(Int32, n_E_bins, n_snapshots)

    # (3) averaged instant failure rate at each timestep
    instant_failures    = zeros(Float64, T+1)

    # (4) averaged cumulative failure rate at each timestep
    cumulative_failures = zeros(Float64, T+1)

    for i in 1:samples
        Ms, Es = sample(L, T, p, q, decoder)

        # (1)
        M1s .+= Ms
        M2s .+= Ms.^2
        E1s .+= Es
        E2s .+= Es.^2

        # (2)
        for (si, t) in enumerate(snapshot_ts)
            Ms_hist[round(Int, Ms[t+1] * L^2) + 1, si] += one(Int32)
            Es_hist[round(Int, Es[t+1]) ÷ 2 + L*(L-1) + 1, si] += one(Int32)
        end

        # (3)
        instant_failures .+= Ms .> 0.5

        # (4)
        ever_failed = false
        for t in 1:T+1
            ever_failed |= Ms[t] > 0.5
            cumulative_failures[t] += ever_failed
        end
    end

    M1s                ./= samples
    M2s                ./= samples
    E1s                ./= samples
    E2s                ./= samples
    instant_failures   ./= samples
    cumulative_failures ./= samples

    return M1s, M2s, E1s, E2s, instant_failures, cumulative_failures, Ms_hist, Es_hist
end