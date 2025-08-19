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

# WARNING: now "horizontal" and "vertical" are on the dual lattice; also now open boundary conditions
function measure(ρ::AbstractMatrix, q::Float64)
    horizontal_checks = noiselayer(ρ .⊻ circshift(ρ,(-1,0)),q)
    vertical_checks = noiselayer(ρ .⊻ circshift(ρ,(0,-1)),q)
    return horizontal_checks, vertical_checks
end