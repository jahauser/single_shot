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

function boundary_filter(checks::Tuple)
    horizontal_checks, vertical_checks = deepcopy(checks)
    L = size(horizontal_checks)[1]
    
    horizontal_checks[end,:] .= false
    vertical_checks[:,end] .= false
    
    return horizontal_checks, vertical_checks
end

function boundary_pad(checks::Tuple)
    horizontal_checks, vertical_checks = checks
    L = size(vertical_checks)[1]
    vcs = zeros(eltype(vertical_checks), L+1, L+1)
    vcs[1:L, 2:L] = vertical_checks[:,1:L-1]

    hcs = zeros(eltype(horizontal_checks), L+1, L+1)
    hcs[2:L, 1:L] = horizontal_checks[1:L-1,:]

    return hcs, vcs
end

function measure(ρ::AbstractMatrix, q::Float64)
    horizontal_checks = noiselayer(ρ .⊻ circshift(ρ,(-1,0)),q)
    vertical_checks = noiselayer(ρ .⊻ circshift(ρ,(0,-1)),q)
    checks = horizontal_checks, vertical_checks
    return boundary_pad(boundary_filter(checks))
end



function detect_charges(hcs::AbstractMatrix, vcs::AbstractMatrix)
    L = size(vcs)[1] - 1

    sites = hcs .⊻ circshift(hcs,(0,1)) .⊻ vcs .⊻ circshift(vcs,(1,0))
    return [(j,i) for i in 0:L, j in 0:L][Bool.(sites)]
end

prob_weight(x::Float64) = -log(x/(1-x))
weight(peff::Float64, q::Float64, m::Bool) = q > 0 ? ((-1)^m * prob_weight(peff) + prob_weight(q)) : (m == true ? 0.0 : Inf)

site(L,x,y) = (L+1)*mod(y,L+1)+mod(x,L+1)+1
unsite(L,s) = (mod1(s,L+1)-1, div(s-1,L+1))

function build_matching_graph(checks::Tuple, peff::Float64, q::Float64)
    horizontal_checks, vertical_checks = checks
    L = size(horizontal_checks)[1] - 1
    
    horizontal_sources = [site(L,x,y) for y in 0:L for x in 0:L-1]
    vertical_sources = [site(L,x,y) for y in 0:L-1 for x in 0:L]
    horizontal_destinations = [site(L,x+1,y) for y in 0:L for x in 0:L-1]
    vertical_destinations = [site(L,x,y+1) for y in 0:L-1 for x in 0:L]
    
    horizontal_weights = zeros(Float64, L+1, L)
    horizontal_weights[2:L, :] = weight.(peff, q, Bool.(horizontal_checks)[2:L,1:L])
    horizontal_weights = reshape(horizontal_weights', L*(L+1), 1)[:,1]

    vertical_weights = zeros(Float64, L, L+1)
    vertical_weights[:,2:L] = weight.(peff, q, Bool.(vertical_checks)[1:L,2:L])
    vertical_weights = reshape(vertical_weights', L*(L+1), 1)[:,1]
    
    return SimpleWeightedGraph([horizontal_sources; vertical_sources], [horizontal_destinations; vertical_destinations], [horizontal_weights; vertical_weights])
end

function build_matching_graph_syndrome_only(checks::Tuple)
    hcs, vcs = checks
    L = size(hcs, 1) - 1

    srcs  = Int[]
    dsts  = Int[]
    wts   = Float64[]

    # horizontal edges
    for y in 1:L-1, x in 0:L-1
        if Bool(hcs[y+1,x+1])
            push!(srcs, site(L,x,y))
            push!(dsts, site(L,x+1,y))
            push!(wts, 1.0)
        end
    end
    for x in 0:L-1
        push!(srcs, site(L, x, 0))
        push!(dsts, site(L, x+1, 0))
        push!(wts, 0.0)

        push!(srcs, site(L, x, L))
        push!(dsts, site(L, x+1, L))
        push!(wts, 0.0)
    end

    # vertical edges
    for y in 0:L-1, x in 1:L-1
        if Bool(vcs[y+1,x+1])
            push!(srcs, site(L, x, y))
            push!(dsts, site(L, x, y+1))
            push!(wts, 1.0)
        end
    end
    for y in 0:L-1
        push!(srcs, site(L, 0, y))
        push!(dsts, site(L, 0, y+1))
        push!(wts, 0.0)

        push!(srcs, site(L, L, y))
        push!(dsts, site(L, L, y+1))
        push!(wts, 0.0)
    end

    return SimpleWeightedGraph(srcs, dsts, wts)
end

function match_charges(fw::Graphs.FloydWarshallState, charges::Vector, L::Int)
    subgraph = complete_graph(length(charges))

    weights = Dict{Edge,Float64}()
    for i in 1:length(charges)-1
        for j in i+1:length(charges)
            if fw.dists[site(L,charges[i]...),site(L,charges[j]...)] < Inf
                weights[Edge(i, j)] = fw.dists[site(L,charges[i]...),site(L,charges[j]...)]
            else
                weights[Edge(i, j)] = 1000.0
            end
        end
    end

    match = minimum_weight_perfect_matching(subgraph, weights)
    return match
end


function heal(checks::Tuple, peff::Float64, q::Float64)
    hcs, vcs = checks
    L = size(hcs)[1]-1

    charges = detect_charges(hcs, vcs)

    if length(charges) == 0
        return hcs, vcs
    end
    
    if peff > q
        g = build_matching_graph((hcs, vcs), peff, q)
    else
        g = build_matching_graph_syndrome_only((hcs, vcs))
    end
    fw = floyd_warshall_shortest_paths(g)
    match = match_charges(fw, charges, L)

    hcs, vcs = apply_paths((hcs, vcs), fw, match, charges, L)

    return hcs, vcs
end

function apply_paths(checks::Tuple, fw::Graphs.FloydWarshallState, match::MatchingResult, charges::Vector, L::Int)
    hcs, vcs = checks
    mated = Int[]
    
    backup_checks = deepcopy(checks)

    for i in 1:length(charges)
        j = match.mate[i]
        if j in mated
            continue
        end
        
        s0 = site(L, charges[i]...)
        s2 = site(L, charges[j]...)
        while s0 != s2
            s1 = fw.parents[s0, s2]
            x1, y1 = unsite(L, s1)
            x2, y2 = unsite(L, s2)
            if x1 == x2
                vcs[max(y1, y2), x1+1] ⊻= true
            elseif y1 == y2
                hcs[y1+1, max(x1, x2)] ⊻= true
            end

            s2 = s1
        end

        push!(mated, i)
        push!(mated, j)
    end
    return hcs, vcs
end

function track_domains(checks::Tuple)
    hcs, vcs = checks
    hcs = Bool.(hcs)
    vcs = Bool.(vcs)
    L = size(hcs)[1]-1

    domain = zeros(Bool, L, L)
    for y in 1:L
        for x in 1:L
            if x == 1 && y == 1
                continue
            elseif x == 1
                domain[y, x] = domain[y-1, x] ⊻ hcs[y, x]
            else
                domain[y, x] = domain[y, x-1] ⊻ vcs[y, x]
            end
        end
    end
    return domain
end

function correct(ρ::AbstractMatrix, peff::Float64, q::Float64)
    checks = measure(ρ, q)
    
    horizontal_checks, vertical_checks = heal(checks, peff, q)
    domain = track_domains((horizontal_checks, vertical_checks))
    if magnetization(domain) > 0.5
        domain = .!domain
    end
    
    return ρ .⊻ domain
end