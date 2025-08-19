function detect_charges(horizontal_checks::AbstractMatrix, vertical_checks::AbstractMatrix)
    L = size(vertical_checks)[1]
    sites = vertical_checks .⊻ horizontal_checks .⊻ circshift(vertical_checks,(-1,0)) .⊻ circshift(horizontal_checks,(0,-1))
    return [(j,i) for i in 1:L, j in 1:L][Bool.(sites)]
end

prob_weight(x::Float64) = -log(x/(1-x))
weight(peff::Float64, q::Float64, m::Bool) = (-1)^m * prob_weight(peff) + prob_weight(q)

site(L,x,y) = L*(mod1(y,L)-1)+mod1(x,L)
unsite(L,s) = (mod1(s,L), div(s-1,L)+1)

function build_matching_graph(checks::Tuple, peff::Float64, q::Float64)
    horizontal_checks, vertical_checks = checks
    L = size(horizontal_checks)[1]
    
    sources = [site(L,x,y) for y in 1:L for x in 1:L]
    horizontal_destinations = [site(L,x+1,y) for y in 1:L for x in 1:L]
    vertical_destinations = [site(L,x,y+1) for y in 1:L for x in 1:L]
    
    horizontal_weights = weight.(peff, q, Bool.(horizontal_checks))
    horizontal_weights[end,:] .= 0
    horizontal_weights = reshape(circshift(horizontal_weights,(0,-1))', L^2, 1)[:,1]

    vertical_weights = weight.(peff, q, Bool.(vertical_checks))
    vertical_weights[:,end] .= 0
    vertical_weights = reshape(circshift(vertical_weights,(-1,0))', L^2, 1)[:,1]
    
    return SimpleWeightedGraph([sources; sources], [horizontal_destinations; vertical_destinations], [horizontal_weights; vertical_weights])
end

function match_charges(fw::Graphs.FloydWarshallState, charges::Vector, L::Int; depth=50)
    subgraph = complete_graph(length(charges))
    weights = Dict{Edge,Float64}()
    for i in 1:length(charges)-1
        for j in i+1:length(charges)
            # if fw.dists[charges[i],charges[j]] < Inf
                weights[Edge(i, j)] = fw.dists[site(L,charges[i]...),site(L,charges[j]...)]
            # end
        end
    end

    match = minimum_weight_perfect_matching(subgraph, weights, depth)
    return match
end

function boundary_filter(checks::Tuple)
    horizontal_checks, vertical_checks = deepcopy(checks)
    L = size(horizontal_checks)[1]
    
    horizontal_checks[end,:] .= false
    vertical_checks[:,end] .= false
    
    return horizontal_checks, vertical_checks
end

cyclicmax(L, a, b) = mod(a-b-1, L) == 0 ? a : b 


function heal(checks::Tuple, peff::Float64, q::Float64)
    horizontal_checks, vertical_checks = deepcopy(boundary_filter(checks))
    L = size(horizontal_checks)[1]
    charges = detect_charges(horizontal_checks, vertical_checks)

    # println(charges)

    if length(charges) == 0
        return horizontal_checks, vertical_checks
    end
    
    g = build_matching_graph((horizontal_checks, vertical_checks), peff, q)
    fw = floyd_warshall_shortest_paths(g)
    match = match_charges(fw, charges, L)

    mated = Int[]


    for i in 1:length(charges)
        j = match.mate[i]
        if j in mated
            continue
        end

        # println(charges[i], charges[j])

        path = enumerate_paths(fw)[site(L,charges[i]...)][site(L,charges[j]...)]

        sites = unsite.(L, path)
        steps = [((x1,y1),(x2,y2)) for ((x1,y1),(x2,y2)) in zip(sites[1:end-1], sites[2:end])]
        # println(steps)
        for ((x1, y1), (x2, y2)) in steps
            if x1 == x2
                vertical_checks[cyclicmax(L, y1, y2), x1] ⊻= true
            elseif y1 == y2
                horizontal_checks[y1, cyclicmax(L, x1, x2)] ⊻= true
            end
        end

        push!(mated, i)
        push!(mated, j)
    end

    return horizontal_checks, vertical_checks
end

function linear_balanced_path(dx::Int, dy::Int; X=:X, Y=:Y)
    q, r = divrem(dx, dy + 1)                # distribute dx X's into dy+1 bins
    steps = Symbol[]
    for i in 1:dy+1
        append!(steps, fill(X, q + (i <= r ? 1 : 0)))
        if i <= dy
            push!(steps, Y)
        end
    end
    steps
end

function basic_heal(checks::Tuple)
    horizontal_checks, vertical_checks = deepcopy(boundary_filter(checks))
    L = size(horizontal_checks)[1]
    charges = detect_charges(horizontal_checks, vertical_checks)

    # println(charges)

    if length(charges) == 0
        return horizontal_checks, vertical_checks
    end
    
    g = build_matching_graph((horizontal_checks, vertical_checks), 0.5, 0.1)
    fw = floyd_warshall_shortest_paths(g)
    pairings = match_charges(fw, charges, L)

    paired = Int[]
    
    for i in 1:length(charges)
        if i in paired
            continue
        end
        j = pairings.mate[i]
        push!(paired, j)

        path = enumerate_paths(fw)[site(L,charges[i]...)][site(L,charges[j]...)]
        sites = unsite.(L, path)
        if L in [max(pos...) for pos in sites]
            # println("whoops")
            steps = [((x1,y1),(x2,y2)) for ((x1,y1),(x2,y2)) in zip(sites[1:end-1], sites[2:end])]
            # println(steps)
            for ((x1, y1), (x2, y2)) in steps
                if x1 == x2
                    vertical_checks[cyclicmax(L, y1, y2), x1] ⊻= true
                elseif y1 == y2
                    horizontal_checks[y1, cyclicmax(L, x1, x2)] ⊻= true
                end
            end 
        else
            if sites[1][1] < sites[end][1]
                xi, yi = sites[1]
                xf, yf = sites[end]
            else
                xi, yi = sites[end]
                xf, yf = sites[1]
            end
            
            dx = abs(xi - xf)
            dy = abs(yi - yf)
            if dx > dy
                seq = linear_balanced_path(dx, dy, X=:X, Y=:Y)
            else
                seq = linear_balanced_path(dy, dx, X=:Y, Y=:X)
            end

            x1 = xi
            y1 = yi

            # println(seq)

            steps = Tuple{Tuple{Int,Int},Tuple{Int,Int}}[]


            for step in seq
                if step == :X
                    steps = push!(steps, ((x1, y1), (x1+1, y1)))
                    x1 = x1+1
                else
                    steps = push!(steps, ((x1, y1), (x1, y1+sign(yf-yi))))
                    y1 = y1+sign(yf-yi)
                end
            end
            # println(steps)
            for ((x1, y1), (x2, y2)) in steps
                if x1 == x2
                    vertical_checks[cyclicmax(L, y1, y2), x1] ⊻= true
                    # println("y1: $y1, y2: $y2, modmax: $(cyclicmax(L, y1, y2))")
                elseif y1 == y2
                    horizontal_checks[y1, cyclicmax(L, x1, x2)] ⊻= true
                    # println("x1: $x1, x2: $x2, modmax: $(cyclicmax(L, x1, x2))")
                end
            end
        end
    end
    return horizontal_checks, vertical_checks
end

function track_domains(checks::Tuple)
    horizontal_checks, vertical_checks = checks
    horizontal_checks = Bool.(horizontal_checks)
    vertical_checks = Bool.(vertical_checks)
    L = size(horizontal_checks)[1]
    
    domain = zeros(Bool, L, L)
    for y in 1:L
        for x in 1:L
            if x == 1 && y == 1
                continue
            elseif x == 1
                domain[y, x] = domain[y-1, x] ⊻ horizontal_checks[y-1, x]
            else
                domain[y, x] = domain[y, x-1] ⊻ vertical_checks[y, x-1]
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


function basic_correct(ρ::AbstractMatrix, q::Float64)
    checks = measure(ρ, q)
    
    horizontal_checks, vertical_checks = basic_heal(checks)
    domain = track_domains((horizontal_checks, vertical_checks))
    if magnetization(domain) > 0.5
        domain = .!domain
    end
    
    return ρ .⊻ domain
end