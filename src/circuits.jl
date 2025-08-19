function adv_MWPM_sample(L::Int, T::Int, p::Float64, peff::Float64, q::Float64)
    ρ = initialstate(L)
    Ms = zeros(Float64, T+1)

    for t in 1:T
        ρ = noiselayer(ρ, p)
        ρ = correct(ρ, peff, q)
        Ms[t+1] = magnetization(ρ)
    end

    return Ms
end

function adv_MWPM_sample(L::Int, T::Int, p::Float64, peff::Float64, q::Float64, samples::Int)
    M1s = zeros(Float64, T+1)
    M2s = zeros(Float64, T+1)
    for i in 1:samples
        Ms = adv_MWPM_sample(L, T, p, peff, q)
        M1s .+= Ms
        M2s .+= Ms.^2
    end
    M1s ./= samples
    M2s ./= samples
    return M1s, M2s
end

function basic_MWPM_sample(L::Int, T::Int, p::Float64, q::Float64)
    ρ = initialstate(L)
    Ms = zeros(Float64, T+1)

    for t in 1:T
        ρ = noiselayer(ρ, p)
        ρ  = basic_correct(ρ, q)
        Ms[t+1] = magnetization(ρ)
    end

    return Ms
end

function basic_MWPM_sample(L::Int, T::Int, p::Float64, q::Float64, samples::Int)
    M1s = zeros(Float64, T+1)
    M2s = zeros(Float64, T+1)
    for i in 1:samples
        Ms = basic_MWPM_circuit(L, T, p, q)
        M1s .+= Ms
        M2s .+= Ms.^2
    end
    M1s ./= samples
    M2s ./= samples
    return M1s, M2s
end

function MV_sample(L::Int, T::Int, p::Float64, q::Float64)
    ρ = initialstate(L)
    Ms = zeros(Float64, T+1)

    Ms[1] = magnetization(ρ)

    for t in 1:T
        ρ = noiselayer(ρ, p)
        checks = measure(ρ, q)
        ρ = MV(ρ, checks)
        Ms[t+1] = magnetization(ρ)
    end
    return Ms
end

function MV_sample(L::Int, T::Int, p::Float64, q::Float64, samples::Int)
    M1s = zeros(Float64, T+1)
    M2s = zeros(Float64, T+1)
    for _ in 1:samples
        Ms = MV_sample(L, T, p, q)
        M1s .+= Ms
        M2s .+= Ms.^2
    end
    M1s ./= samples
    M2s ./= samples
    return M1s, M2s
end

function TOOM_sample(L::Int, T::Int, p::Float64, q::Float64)
    ρ = initialstate(L)
    Ms = zeros(Float64, T+1)

    Ms[1] = magnetization(ρ)

    for t in 1:T
        ρ = noiselayer(ρ, p)
        checks = measure(ρ, q)
        ρ = toom(ρ, checks)
        Ms[t+1] = magnetization(ρ)
    end
    return Ms
end

function TOOM_sample(L::Int, T::Int, p::Float64, q::Float64, samples::Int)
    M1s = zeros(Float64, T+1)
    M2s = zeros(Float64, T+1)
    for _ in 1:samples
        Ms = TOOM_sample(L, T, p, q)
        M1s .+= Ms
        M2s .+= Ms.^2
    end
    M1s ./= samples
    M2s ./= samples
    return M1s, M2s
end