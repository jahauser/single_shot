function majority_vote(ρ::AbstractMatrix, checks::Tuple)
    horizontal_checks, vertical_checks = checks
    total = horizontal_checks + vertical_checks + circshift(horizontal_checks,(1,0)) + circshift(vertical_checks,(0,1))
    ρ = ρ .⊻ (total .> 2)
    return ρ
end

function toom(ρ::AbstractMatrix, checks::Tuple)
    horizontal_checks, vertical_checks = checks
    ρ = ρ .⊻ (vertical_checks .& horizontal_checks)
    return ρ
end

function unanimous_vote(ρ::AbstractMatrix, checks::Tuple)
    horizontal_checks, vertical_checks = checks
    total = horizontal_checks + vertical_checks + circshift(horizontal_checks,(1,0)) + circshift(vertical_checks,(0,1))
    ρ = ρ .⊻ (total .> 3)
    return ρ
end

# function glauber(ρ::AbstractMatrix, checks::Tuple, β::Float64)
#     horizontal_checks, vertical_checks = checks
#     total = horizontal_checks + vertical_checks + circshift(horizontal_checks,(1,0)) + circshift(vertical_checks,(0,1))
#     ΔE = 2 .- total
#     return ρ
# end