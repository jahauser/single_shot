function MV(ρ::AbstractMatrix, checks::Tuple)
    vertical_checks, horizontal_checks = checks
    total = vertical_checks + horizontal_checks + circshift(vertical_checks,(1,0)) + circshift(horizontal_checks,(0,1))
    ρ = ρ .⊻ (total .> 2)
    return ρ
end

function toom(ρ::AbstractMatrix, checks::Tuple)
    vertical_checks, horizontal_checks = checks
    ρ = ρ .⊻ (vertical_checks .& horizontal_checks)
    return ρ
end