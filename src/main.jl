function initialstate(L::Int)
    ρ = zeros(Bool,L,L)
    return ρ
end

function noiselayer(ρ::AbstractMatrix, p::Float64)
    L = size(ρ)[1]
    noise = rand(L,L) .< p
    return ρ .⊻ noise
end

function magnetization(ρ::AbstractMatrix)
    return mean(ρ)
end

function measure(ρ::AbstractMatrix, q::Float64)
    horizontal_checks = noiselayer(ρ .⊻ circshift(ρ,(-1,0)),q)
    vertical_checks = noiselayer(ρ .⊻ circshift(ρ,(0,-1)),q)
    return horizontal_checks, vertical_checks
end

function energy(ρ::AbstractMatrix)
    horizontal_checks, vertical_checks = measure(ρ, 0.0)
    horizontal_ss = sum(1 .- 2*horizontal_checks[1:end-1,:])
    vertical_ss = sum(1 .- 2*vertical_checks[:,1:end-1])

    return horizontal_ss + vertical_ss
end