using WGLMakie
using Bonito
using Printf
using StatsBase
using Graphs
using SimpleWeightedGraphs
using GraphsMatching

# include("src/main.jl")
# include("src/mwpm.jl")
include("src/web_src.jl")

# =============================================================
# Simulation  (abort-aware version)
# =============================================================

# Thrown when the user clicks Abort mid-simulation
struct AbortSimulation <: Exception end

function adv_MWPM_sample_hotstart(ρ0::Array{Bool,2}, L::Int, T::Int,
                                   p::Float64, peff::Float64, q::Float64;
                                   abort_flag::Threads.Atomic{Bool},
                                   progress_obs::Observable{Float64})
    ρ = deepcopy(ρ0)
    Ms                = zeros(Bool, 2T+2, L, L)
    # checks are now (L+1, L+1) after boundary_pad
    horizontal_checks = zeros(Bool, 2T+2, L+1, L+1)
    vertical_checks   = zeros(Bool, 2T+2, L+1, L+1)
    for t in 1:T
        abort_flag[] && throw(AbortSimulation())
        ρ = noiselayer(ρ, p)
        Ms[2t+1, :, :] = ρ
        ρ, orig, new_c = correct(ρ, peff, q)
        horizontal_checks[2t+1, :, :] = orig[1]
        vertical_checks[2t+1,   :, :] = orig[2]
        horizontal_checks[2t+2, :, :] = new_c[1]
        vertical_checks[2t+2,   :, :] = new_c[2]
        Ms[2t+2, :, :] = ρ
        progress_obs[] = t / T   # 0.0 → 1.0
    end
    return Ms, horizontal_checks, vertical_checks
end

function correct(ρ::AbstractMatrix, peff::Float64, q::Float64)
    checks = measure(ρ, q)
    original_checks = deepcopy(checks)
    hc, vc = heal(checks, peff, q)
    domain = track_domains((hc, vc))
    magnetization(domain) > 0.5 && (domain = .!domain)
    return ρ .⊻ domain, original_checks, (hc, vc)
end

function make_initial_state(L::Int, mode::String)::Matrix{Bool}
    ρ = zeros(Bool, L, L)
    if mode == "Square patch"
        ℓ = max(2, L ÷ 4)
        r1 = clamp(L÷2 - ℓ÷2, 1, L)
        r2 = clamp(L÷2 + ℓ÷2, 1, L)
        ρ[r1:r2, r1:r2] .= true
    elseif mode == "Random"
        ρ .= rand(Bool, L, L)
    end  # "All zeros" → already zero
    return ρ
end

function run_simulation(L::Int, T::Int, p::Float64, peff::Float64, q::Float64,
                        init_mode::String="Square patch";
                        abort_flag::Threads.Atomic{Bool}  = Threads.Atomic{Bool}(false),
                        progress_obs::Observable{Float64} = Observable(0.0))
    ρ = make_initial_state(L, init_mode)
    M1s, hc, vc = adv_MWPM_sample_hotstart(ρ, L, T, p, peff, q;
                                            abort_flag=abort_flag,
                                            progress_obs=progress_obs)
    M1s[1, :, :] = ρ
    M1s[2, :, :] = ρ
    Nt_old, Lx, Ly = size(M1s)
    Nt = 2Nt_old
    K  = Nt ÷ 4
    M = similar(M1s, Nt, Lx, Ly)
    # H and V are (L+1, L+1) per frame
    H = similar(hc,  Nt, L+1, L+1)
    V = similar(vc,  Nt, L+1, L+1)
    for k in 1:K
        t0 = 4k - 3
        M[t0:t0+2, :, :] .= reshape(M1s[2k-1, :, :], 1, Lx, Ly)
        M[t0+3,    :, :]  .= M1s[2k, :, :]
    end
    fill!(H, false); fill!(V, false)
    for k in 1:K
        t0 = 4k - 3
        H[t0+1, :, :]      .= hc[2k-1, :, :]
        V[t0+1, :, :]      .= vc[2k-1, :, :]
        H[t0+2:t0+3, :, :] .= reshape(hc[2k, :, :], 1, L+1, L+1)
        V[t0+2:t0+3, :, :] .= reshape(vc[2k, :, :], 1, L+1, L+1)
    end
    raw_H = Dict{Int,Matrix{Bool}}()
    raw_V = Dict{Int,Matrix{Bool}}()
    for k in 2:K
        sidx = 4k - 5
        raw_H[sidx] = Bool.(hc[2k-1, :, :])
        raw_V[sidx] = Bool.(vc[2k-1, :, :])
    end
    return M[4:end,:,:], H[4:end,:,:], V[4:end,:,:], raw_H, raw_V
end

# =============================================================
# JIT warmup
# =============================================================

@info "Warming up JIT..."
@async begin
    try
        run_simulation(8,  4, 0.1, 0.2, 0.1, "Square patch")
        run_simulation(16, 8, 0.1, 0.2, 0.1, "Square patch")
        @info "JIT warmup complete."
    catch e
        @warn "Warmup failed (non-fatal)" exception=e
    end
end

# =============================================================
# Input parsers
# =============================================================

function parse_T_field(s::AbstractString, L::Int)::Union{Int,Nothing}
    s = strip(s)
    v = tryparse(Int, s)
    v !== nothing && return max(1, v)
    m = match(r"^([0-9]*\.?[0-9]+)\s*\*?\s*[Ll]$", s)
    if m !== nothing
        c = tryparse(Float64, m.captures[1])
        c !== nothing && return max(1, round(Int, c * L))
    end
    return nothing
end

function parse_prob(s::AbstractString, lo::Float64=0.0, hi::Float64=0.5)::Union{Float64,Nothing}
    v = tryparse(Float64, strip(s))
    (v === nothing || v < lo || v > hi) && return nothing
    return v
end

# =============================================================
# Edge segments
# =============================================================

function edge_segments(hc::AbstractMatrix, vc::AbstractMatrix)
    # hc and vc are (L+1, L+1) padded check matrices from boundary_pad.
    # L = size - 1; domain cells live at 1..L.
    #
    # WGLMakie heatmap maps first matrix index → x-axis, second → y-axis,
    # so matrix [row, col] appears at screen position (row, col).
    # "row" (y in track_domains) is therefore horizontal on screen.
    #
    # hc[y, x]: bond between rows y-1 and y at col x
    #   → domain wall is a vertical line on screen at screen-x = y-0.5,
    #     spanning screen-y from x-0.5 to x+0.5
    #   → Point2f(y-0.5, x±0.5)
    #
    # vc[y, x]: bond between cols x-1 and x at row y
    #   → domain wall is a horizontal line on screen at screen-y = x-0.5,
    #     spanning screen-x from y-0.5 to y+0.5
    #   → Point2f(y±0.5, x-0.5)
    #
    # Boundary edges (those that would fall outside the 1..L cell range on
    # either axis) are suppressed — they live on the padding rows/cols and
    # have no corresponding heatmap cell to border.

    L = size(hc, 1) - 1   # hc is (L+1, L+1)

    adj  = Dict{Tuple{Int,Int}, Vector{Tuple{Int,Int}}}()
    add_edge! = function(a, b)
        push!(get!(adj, a, Tuple{Int,Int}[]), b)
        push!(get!(adj, b, Tuple{Int,Int}[]), a)
    end

    key(px, py) = (round(Int, 2*px), round(Int, 2*py))

    # hc[y,x] → vertical wall on screen at (y-0.5, x±0.5)
    # Only interior bonds: y in 2..L (walls between rows 1..L-1 and 2..L),
    # x in 1..L (wall spans screen-y 0.5..1.5 .. L-0.5..L+0.5, all inside grid)
    for y in 2:L, x in 1:L
        hc[y, x] || continue
        yw = Float32(y) - 0.5f0
        add_edge!(key(yw, Float32(x) - 0.5f0), key(yw, Float32(x) + 0.5f0))
    end

    # vc[y,x] → horizontal wall on screen at (y±0.5, x-0.5)
    # Only interior bonds: y in 1..L, x in 2..L
    for y in 1:L, x in 2:L
        vc[y, x] || continue
        xw = Float32(x) - 0.5f0
        add_edge!(key(Float32(y) - 0.5f0, xw), key(Float32(y) + 0.5f0, xw))
    end

    isempty(adj) && return [Point2f(-10, -10), Point2f(-10, -10)]

    # Greedily trace paths: prefer to continue in the same direction (avoids
    # unnecessary zig-zags through degree-2 nodes), fall back to any unvisited neighbour.
    used   = Set{Tuple{Tuple{Int,Int},Tuple{Int,Int}}}()
    result = Point2f[]

    pt(k) = Point2f(k[1] * 0.5f0, k[2] * 0.5f0)

    for start in keys(adj)
        while any(n -> !((start,n) in used), adj[start])
            chain = Tuple{Int,Int}[start]
            prev  = start
            cur = first(n for n in adj[prev] if !((prev,n) in used))
            push!(used, (prev, cur)); push!(used, (cur, prev))
            push!(chain, cur)

            while true
                dx = chain[end][1] - chain[end-1][1]
                dy = chain[end][2] - chain[end-1][2]
                node = chain[end]
                straight = (node[1]+dx, node[2]+dy)
                nxt = nothing
                if straight in keys(adj) && !((node, straight) in used) &&
                        straight in adj[node]
                    nxt = straight
                else
                    for n in adj[node]
                        (node, n) in used && continue
                        nxt = n; break
                    end
                end
                nxt === nothing && break
                push!(used, (node, nxt)); push!(used, (nxt, node))
                push!(chain, nxt)
            end

            isempty(result) || push!(result, Point2f(NaN, NaN))
            append!(result, pt.(chain))
        end
    end

    return result
end

# =============================================================
# App
# =============================================================

app = App() do session

    # ── Text inputs ──────────────────────────────────────────────────────────
    tf_L = Bonito.TextField("16")
    tf_T = Bonito.TextField("2L")
    tf_p = Bonito.TextField("0.05")
    tf_q = Bonito.TextField("0.10")

    # ── c slider (syndrome confidence, 0.0→1.0) ──────────────────────────────
    # peff is derived: peff = q + (0.5 - q) * c
    c_vals = collect(0.00:0.01:1.00)
    sl_c   = Bonito.Slider(c_vals; value=1.00)
    c_disp = Observable("1.00")
    # Shown below the slider when explore mode is active
    explore_status_obs = Observable("")
    # True while a heal() call is running — drives the spinner
    healing_obs = Observable(false)

    # ── Render worker ─────────────────────────────────────────────────────────
    render_ch = Channel{Tuple{Union{Matrix{Bool},Nothing}, Matrix{Bool}, Matrix{Bool}}}(1)

    @async begin
        for (heat_mat, hc_mat, vc_mat) in render_ch
            try
                new_edges = edge_segments(hc_mat, vc_mat)
                if heat_mat !== nothing
                    heat_obs[] = Float32.(heat_mat)
                end
                edges_obs[] = new_edges
            catch e
                @warn "render worker error" exception=e
            end
        end
    end

    # Helper: drain-and-replace put into a capacity-1 channel
    function latest_put!(ch, val)
        isopen(ch) || return
        isready(ch) && take!(ch)
        put!(ch, val)
    end

    # ── Debounced heal worker ─────────────────────────────────────────────────
    heal_ch = Channel{Tuple{Float64,Int,Matrix{Bool},Matrix{Bool},Matrix{Bool},Float64}}(1)

    @async begin
        for (v, sidx, rh, rv, rho_noisy, q_cur) in heal_ch
            healing_obs[] = true
            try
                peff = q_cur + (0.5 - q_cur) * v
                hc_new, vc_new = heal((deepcopy(rh), deepcopy(rv)), peff, q_cur)

                domain = track_domains((hc_new, vc_new))
                magnetization(domain) > 0.5 && (domain = .!domain)
                rho_corrected = rho_noisy .⊻ domain

                explore_hc[]     = hc_new
                explore_vc[]     = vc_new
                explore_rho[]    = rho_corrected
                explore_active[] = true
                explore_sidx[]   = sidx
                explore_status_obs[] = @sprintf(
                    "Explore mode — c = %.2f  (peff = %.3f) · data beyond this step hidden",
                    v, peff)

                cur_t = time_ref[]
                if cur_t == sidx
                    latest_put!(render_ch, (rho_noisy, rh, rv))
                elseif cur_t == sidx + 1
                    latest_put!(render_ch, (rho_noisy, hc_new, vc_new))
                elseif cur_t == sidx + 2
                    latest_put!(render_ch, (rho_corrected, hc_new, vc_new))
                end
            catch e
                @warn "heal worker error" exception=e
            finally
                healing_obs[] = false
            end
        end
    end

    on(sl_c.value) do v
        c_disp[] = @sprintf("%.2f", v)
        ready[] || return
        M, H, V, raw_H, raw_V = sim_ref[]
        cur_t = time_ref[]
        sidx = 0
        for candidate in cur_t:-1:1
            if haskey(raw_H, candidate)
                sidx = candidate
                break
            end
        end
        sidx == 0 && return
        rho_noisy = Bool.(M[sidx, :, :])
        latest_put!(heal_ch, (v, sidx, raw_H[sidx], raw_V[sidx], rho_noisy, q_ref[]))
    end

    # ── Time scrubber ─────────────────────────────────────────────────────────
    sl_time    = Bonito.Slider(1:1:1; value=1)
    time_label = Observable("—")

    # ── Initial state selector (graphical SVG radio buttons) ─────────────────
    init_mode_obs = Observable("Square patch")

    init_click_obs = Observable("")
    on(init_click_obs) do mode
        isempty(mode) && return
        init_mode_obs[] = mode
        init_click_obs[] = ""
    end

    _thumb_base   = "display:block;width:48px;height:48px;border-radius:2px;cursor:pointer;" *
                    "box-sizing:border-box;flex-shrink:0;"
    _sel_border   = "border:2px solid #111;"
    _unsel_border = "border:1px solid #ccc;"

    thumb_zeros  = DOM.div(;
        id="ibtn-zeros", class="init-thumb",
        style=_thumb_base * _unsel_border *
              "background:#000;")
    thumb_square = DOM.div(
        DOM.div(;
            style="width:22px;height:22px;background:#fff;border-radius:1px;");
        id="ibtn-square", class="init-thumb",
        style=_thumb_base * _sel_border *
              "background:#000;display:flex;align-items:center;justify-content:center;")
    thumb_random = DOM.div(;
        id="ibtn-random", class="init-thumb",
        style=_thumb_base * _unsel_border *
              "background-color:#000;" *
              "background-image:" *
                "linear-gradient(45deg,#fff 25%,transparent 25%)," *
                "linear-gradient(-45deg,#fff 25%,transparent 25%)," *
                "linear-gradient(45deg,transparent 75%,#fff 75%)," *
                "linear-gradient(-45deg,transparent 75%,#fff 75%);" *
              "background-size:16px 16px;" *
              "background-position:0 0,0 8px,8px -8px,-8px 0;")

    onjs(session, init_mode_obs, js"""(mode) => {
        const map = {"All zeros":"ibtn-zeros","Square patch":"ibtn-square","Random":"ibtn-random"};
        ["ibtn-zeros","ibtn-square","ibtn-random"].forEach(id => {
            const el = document.getElementById(id);
            if (el) { el.style.borderColor = "#ccc"; el.style.borderWidth = "1px"; }
        });
        const sel = document.getElementById(map[mode]);
        if (sel) { sel.style.borderColor = "#111"; sel.style.borderWidth = "2px"; }
    }""")
    evaljs(session, js"""
        ["ibtn-zeros","ibtn-square","ibtn-random"].forEach((id, i) => {
            const modes = ["All zeros","Square patch","Random"];
            const el = document.getElementById(id);
            if (el) el.addEventListener("click", () => $(init_click_obs).notify(modes[i]));
        });
    """)

    # ── Progress bar (0.0–1.0) ────────────────────────────────────────────────
    progress_obs = Observable(0.0)

    # ── Run params banner (frozen when simulation completes, shown above heatmap) ──
    run_params_obs = Observable("")

    # ── Status ────────────────────────────────────────────────────────────────
    status_obs = Observable("No sample yet — set parameters and click Generate.")

    # ── Buttons ──────────────────────────────────────────────────────────────
    btn_generate = Bonito.Button("▶  Generate")
    btn_abort    = Bonito.Button("■  Abort")

    # ── Abort flag & task ref ─────────────────────────────────────────────────
    abort_flag  = Ref(Threads.Atomic{Bool}(false))
    running_ref = Ref{Union{Task,Nothing}}(nothing)

    # ── Simulation state ─────────────────────────────────────────────────────
    L_init  = 4
    sim_ref = Ref{Tuple}((zeros(Bool,1,L_init,L_init),
                          zeros(Bool,1,L_init+1,L_init+1),
                          zeros(Bool,1,L_init+1,L_init+1),
                          Dict{Int,Matrix{Bool}}(),
                          Dict{Int,Matrix{Bool}}()))
    Nt_ref   = Ref{Int}(1)
    time_ref = Ref{Int}(1)
    ready    = Ref{Bool}(false)

    # ── Explore-mode state ────────────────────────────────────────────────────
    explore_active = Ref{Bool}(false)
    explore_sidx   = Ref{Int}(0)
    explore_hc     = Ref{Matrix{Bool}}(zeros(Bool, L_init+1, L_init+1))
    explore_vc     = Ref{Matrix{Bool}}(zeros(Bool, L_init+1, L_init+1))
    explore_rho    = Ref{Matrix{Bool}}(zeros(Bool, L_init, L_init))
    q_ref          = Ref{Float64}(0.1)

    # ── WGLMakie figure ──────────────────────────────────────────────────────
    heat_obs  = Observable(zeros(Float32, L_init, L_init))
    edges_obs = Observable([Point2f(-10,-10), Point2f(-10,-10)])

    fig = Figure(size=(460, 460), backgroundcolor=RGBf(0.18, 0.18, 0.18), figure_padding=8)
    ax  = Axis(fig[1, 1];
               aspect=DataAspect(), backgroundcolor=:black,
               xgridvisible=false, ygridvisible=false,
               bottomspinevisible=false, topspinevisible=false,
               leftspinevisible=false, rightspinevisible=false,
               xautolimitmargin=(0f0, 0f0), yautolimitmargin=(0f0, 0f0))
    hidedecorations!(ax)
    tightlimits!(ax)
    deregister_interaction!(ax, :scrollzoom)
    deregister_interaction!(ax, :rectanglezoom)
    deregister_interaction!(ax, :limitreset)
    heatmap!(ax, heat_obs; colormap=[RGBf(0,0,0), RGBf(1,1,1)], colorrange=(0f0, 1f0))
    lines!(ax, edges_obs; color=RGBf(0.85, 0.08, 0.08), linewidth=4.0,
           joinstyle=:round, linecap=:round)

    # ── Arrow key scrubbing via JS keydown listener ───────────────────────────
    key_obs = Observable(0)
    on(key_obs) do dir
        dir == 0 && return
        ready[] || return
        max_t = explore_active[] ? explore_sidx[] + 2 : Nt_ref[]
        new_t = clamp(time_ref[] + dir, 1, max_t)
        sl_time.index[] = new_t
        key_obs[] = 0
    end
    evaljs(session, js"""
        document.addEventListener('keydown', function(e) {
            if (e.key === 'ArrowLeft')  { $(key_obs).notify(-1); e.preventDefault(); }
            if (e.key === 'ArrowRight') { $(key_obs).notify(1);  e.preventDefault(); }
        });
    """)

    # ── Frame renderer ───────────────────────────────────────────────────────
    function show_frame(t::Int)
        ready[] || return
        M, H, V, raw_H, raw_V = sim_ref[]

        if explore_active[]
            sidx = explore_sidx[]
            tt = clamp(t, 1, sidx + 2)
            if tt < sidx
                explore_active[]     = false
                explore_status_obs[] = ""
                tt = clamp(t, 1, Nt_ref[])
            end
        else
            tt = clamp(t, 1, Nt_ref[])
        end

        time_ref[] = tt

        phys_t = (tt - 1) / 4.0
        phys_T = (Nt_ref[] - 1) / 4.0
        t_str = @sprintf("%.2f", phys_t)
        T_str = isinteger(phys_T) ? "$(round(Int, phys_T))" : "$(round(phys_T, digits=1))"
        time_label[] = "t = $t_str / $T_str"

        if explore_active[]
            sidx = explore_sidx[]
            if tt == sidx
                latest_put!(render_ch, (Bool.(M[sidx,:,:]), raw_H[sidx], raw_V[sidx]))
            elseif tt == sidx + 1
                latest_put!(render_ch, (Bool.(M[sidx,:,:]), explore_hc[], explore_vc[]))
            elseif tt == sidx + 2
                latest_put!(render_ch, (explore_rho[], explore_hc[], explore_vc[]))
            else
                latest_put!(render_ch, (Bool.(M[tt,:,:]), H[tt,:,:], V[tt,:,:]))
            end
        else
            latest_put!(render_ch, (Bool.(M[tt,:,:]), H[tt,:,:], V[tt,:,:]))
        end
    end

    on(sl_time.value) do v
        t = Int(v)
        if ready[] && explore_active[] && t < explore_sidx[]
            explore_active[]     = false
            explore_sidx[]       = 0
            explore_status_obs[] = ""
        end
        show_frame(t)
    end

    # ── Abort handler ─────────────────────────────────────────────────────────
    on(btn_abort.value) do _
        abort_flag[][] = true
    end

    # ── Generate handler ──────────────────────────────────────────────────────
    on(btn_generate.value) do _
        if !isnothing(running_ref[]) && !istaskdone(running_ref[])
            abort_flag[][] = true
            sleep(0.05)
        end

        L_val = tryparse(Int, strip(tf_L.value[]))
        if L_val === nothing || L_val < 2
            status_obs[] = "⚠  L must be an integer ≥ 2."; return
        end

        T_val = parse_T_field(tf_T.value[], L_val)
        T_val === nothing && (status_obs[] = "⚠  T: integer or e.g. '10L'."; return)

        p_val = parse_prob(tf_p.value[], 0.0, 0.5)
        p_val === nothing && (status_obs[] = "⚠  p must be in [0, 0.5]."; return)

        q_val = parse_prob(tf_q.value[], 0.0, 0.5)
        q_val === nothing && (status_obs[] = "⚠  q must be in [0, 0.5]."; return)

        c_val    = sl_c.value[]
        peff_val = q_val + (0.5 - q_val) * c_val

        flag = Threads.Atomic{Bool}(false)
        abort_flag[] = flag

        progress_obs[] = 0.0
        status_obs[]   = "⏳  Running simulation…"

        @info "Generating: L=$L_val T=$T_val p=$p_val q=$q_val peff=$peff_val"

        t_start = time()
        t = Threads.@spawn begin
            try
                M, H, V, raw_H, raw_V = run_simulation(L_val, T_val, p_val, peff_val, q_val,
                                              init_mode_obs[];
                                              abort_flag=flag,
                                              progress_obs=progress_obs)
                sim_ref[]  = (M, H, V, raw_H, raw_V)
                q_ref[]    = q_val
                explore_active[]     = false
                explore_sidx[]       = 0
                explore_status_obs[] = ""
                Nt_ref[]   = size(M, 1)
                ready[]    = true
                time_ref[] = 1

                sl_time.values[] = collect(1:Nt_ref[])
                sl_time.index[]  = 1

                heat_obs[]   = Float32.(M[1, :, :])
                latest_put!(render_ch, (Bool.(M[1,:,:]), H[1,:,:], V[1,:,:]))

                phys_T = (Nt_ref[] - 1) / 4.0
                T_str  = isinteger(phys_T) ? "$(round(Int, phys_T))" : "$(round(phys_T, digits=1))"
                time_label[] = "t = 0.00 / $T_str"
                L_cur = size(M, 2)
                limits!(ax, 0.5, L_cur + 0.5, 0.5, L_cur + 0.5)

                elapsed = round(time() - t_start; digits=1)
                progress_obs[] = 1.0
                status_obs[]   = "✓  Done — $T_str time steps  ($(elapsed)s)"
                run_params_obs[] = @sprintf("L = %d · T = %s · p = %.3f · q = %.3f · c = %.2f",
                                            L_val, T_str, p_val, q_val, c_val)
            catch e
                if e isa AbortSimulation
                    status_obs[] = "⚠  Aborted."
                    @info "Simulation aborted by user."
                else
                    status_obs[] = "⚠  Error: $(sprint(showerror, e))"
                    @error "Simulation failed" exception=(e, catch_backtrace())
                end
            end
        end
        running_ref[] = t
    end

    # ── CSS ──────────────────────────────────────────────────────────────────
    style_tag = DOM.style("""
        *, *::before, *::after { box-sizing: border-box; }

        body {
            margin: 0;
            padding: 20px 16px;
            min-height: 100vh;
            background: #f0f0f0;
            font-family: system-ui, sans-serif;
            font-size: 12px;
            color: #111;
            display: flex;
            justify-content: center;
            align-items: flex-start;
        }

        /* ── Card ── */
        .app-card {
            background: #fff;
            border: 1px solid #ddd;
            border-radius: 4px;
            padding: 24px;
            width: 100%;
            max-width: 1200px;
        }

        /* ── Title ── */
        .app-title {
            margin-bottom: 18px;
            padding-bottom: 10px;
            border-bottom: 1px solid #ddd;
            font-size: 13px;
            font-weight: 600;
            color: #111;
        }
        .app-title .subtitle {
            font-weight: 400;
            color: #888;
        }

        /* ── Two-column layout ── */
        .app-layout {
            display: grid;
            grid-template-columns: 280px 1fr;
            gap: 0 24px;
            align-items: start;
        }

        /* ── Settings panel ── */
        .app-settings {
            background: #fff;
            border: 1px solid #ddd;
            border-radius: 3px;
            padding: 12px 14px;
        }
        .param-row {
            display: flex;
            align-items: center;
            gap: 8px;
            padding: 6px 0;
            border-bottom: 1px solid #eee;
        }
        .param-row:last-of-type { border-bottom: none; }
        .param-label {
            width: 108px;
            flex-shrink: 0;
            color: #555;
            font-weight: 500;
        }
        .param-input input {
            width: 72px !important;
            background: #fff !important;
            color: #111 !important;
            border: 1px solid #ccc !important;
            border-radius: 3px !important;
            padding: 3px 6px !important;
            font-family: monospace !important;
            font-size: 12px !important;
            outline: none !important;
            transition: border-color 0.15s !important;
        }
        .param-input input:focus { border-color: #888 !important; }

        /* c-slider row */
        .c-row {
            padding: 6px 0;
            border-bottom: 1px solid #eee;
        }
        .c-row-header {
            display: flex;
            align-items: baseline;
            gap: 6px;
            margin-bottom: 4px;
        }
        .c-row-header .param-label { width: auto; flex: 1; }
        .c-val { font-family: monospace; color: #111; flex-shrink: 0; }
        #heal-spinner {
            flex-shrink: 0;
            display: none;
            animation: spin 0.8s linear infinite;
            transform-origin: center;
        }
        .c-slider-wrap { width: 100%; }
        .c-slider-wrap input[type=range] {
            width: 100% !important;
            accent-color: #111 !important;
            cursor: pointer !important;
        }

        /* Initial state thumbs */
        .init-section { padding: 6px 0; border-bottom: 1px solid #eee; }
        .init-label { color: #555; font-weight: 500; margin-bottom: 6px; display: block; }
        .init-thumbs { display: flex; gap: 10px; }
        .init-thumb-wrap { display: flex; flex-direction: column; align-items: center; gap: 3px; }
        .init-thumb {
            display: block;
            width: 44px; height: 44px;
            border-radius: 2px;
            cursor: pointer;
            border: 1px solid #ccc;
            flex-shrink: 0;
        }
        .init-thumb:hover { border-color: #555 !important; }
        .init-thumb-lbl { color: #888; font-size: 10px; }

        /* Buttons */
        .btn-row { margin-top: 14px; display: flex; flex-direction: column; gap: 6px; }
        .generate-btn button, .abort-btn button {
            width: 100% !important;
            border-radius: 3px !important;
            font-size: 12px !important;
            cursor: pointer !important;
            padding: 6px 0 !important;
        }
        .generate-btn button {
            background: #111 !important;
            color: #fff !important;
            border: none !important;
        }
        .generate-btn button:hover { background: #333 !important; }
        .abort-btn button {
            background: #fff !important;
            color: #666 !important;
            border: 1px solid #ccc !important;
        }
        .abort-btn button:hover { color: #c00 !important; border-color: #c00 !important; }

        /* Progress + status */
        .progress-track {
            width: 100%; height: 2px;
            background: #eee;
            margin-top: 12px;
            overflow: hidden;
        }
        .progress-fill { height: 100%; width: 0%; background: #111; transition: width 0.18s ease; }
        .status-text {
            margin-top: 8px;
            min-height: 15px;
            color: #888;
            line-height: 1.4;
        }

        /* ── Viz column ── */
        .app-viz {
            display: flex;
            flex-direction: column;
            align-items: stretch;
            min-width: 0;
        }

        /* Canvas: square, capped so both banners + time controls fit in viewport */
        .app-viz canvas {
            width: 100% !important;
            height: auto !important;
            aspect-ratio: 1 / 1;
            display: block;
            min-width: 200px;
            max-width: calc(100vh - 220px);
            max-height: calc(100vh - 220px);
            image-rendering: pixelated;
        }
        /* Makie wraps canvas in several divs — make them all fill width */
        .makie-wrap, .makie-wrap > div, .makie-wrap > div > div { width: 100% !important; }

        /* ── Banners ── */
        .banner {
            padding: 4px 8px;
            border-radius: 3px;
            width: 100%;
            font-family: monospace;
            font-size: 11px;
            display: none;
        }
        .banner-params {
            background: #f0f7f0;
            border: 1px solid #9fc89f;
            color: #1a3a1a;
            margin-bottom: 5px;
        }
        .banner-explore {
            background: #fff8e6;
            border: 1px solid #e8c96b;
            color: #7a4f00;
            margin-top: 5px;
        }

        /* ── Time controls ── */
        .time-controls {
            display: flex;
            align-items: center;
            gap: 6px;
            margin-top: 8px;
            width: 100%;
        }
        .time-nav {
            background: #fff;
            color: #555;
            border: 1px solid #ccc;
            border-radius: 3px;
            padding: 3px 8px;
            font-size: 11px;
            cursor: pointer;
            line-height: 1;
            flex-shrink: 0;
        }
        .time-nav:hover { color: #111; border-color: #888; }
        .time-slider {
            flex: 1;
            min-width: 0;
            display: flex;
            align-items: center;
        }
        .time-slider, .time-slider > div, .time-slider > div > div {
            width: 100% !important;
            flex: 1 !important;
            min-width: 0 !important;
        }
        .time-slider input[type=range] {
            width: 100% !important;
            accent-color: #111 !important;
            cursor: pointer !important;
        }
        .time-label {
            font-family: monospace;
            font-size: 11px;
            color: #888;
            white-space: nowrap;
            text-align: right;
            width: 17ch;
            flex-shrink: 0;
        }
        .explore-active .time-slider input[type=range] {
            accent-color: #c8a000 !important;
            opacity: 0.65;
        }
        .explore-active #nav-next { color: #bbb !important; border-color: #ddd !important; cursor: default !important; }

        /* ── Mobile ── */
        @media (max-width: 680px) {
            body { padding: 12px; }
            .app-card { padding: 16px; }
            .app-layout { grid-template-columns: 1fr; }
            .app-settings { margin-bottom: 16px; }
            .app-viz canvas { max-width: 100%; max-height: none; }
        }

        @keyframes spin { from { transform: rotate(0deg); } to { transform: rotate(360deg); } }
    """)

    lbl(text) = DOM.span(text; class="param-label")

    progress_fill_node  = DOM.div(; class="progress-fill")
    progress_track_node = DOM.div(progress_fill_node; class="progress-track")
    onjs(session, progress_obs, js"""(v) => {
        const el = $(progress_fill_node);
        if (el) el.style.width = Math.round(Math.min(v, 1.0) * 100) + '%';
    }""")
    onjs(session, healing_obs, js"""(busy) => {
        const el = document.getElementById('heal-spinner');
        if (el) el.style.display = busy ? 'inline' : 'none';
    }""")

    param_panel = DOM.div(
        DOM.div(lbl("L (grid size)"),       DOM.div(tf_L; class="param-input"); class="param-row"),
        DOM.div(lbl("T (timesteps)"),        DOM.div(tf_T; class="param-input"); class="param-row"),
        DOM.div(lbl("p  (physical noise)"),  DOM.div(tf_p; class="param-input"); class="param-row"),
        DOM.div(lbl("q  (meas. noise)"),     DOM.div(tf_q; class="param-input"); class="param-row"),
        DOM.div(
            DOM.div(
                DOM.span("c  (syndrome confidence)"; class="param-label", style="width:auto;flex:1;"),
                DOM.span(c_disp; class="c-val"),
                DOM.span("\u27f3"; id="heal-spinner");
                class="c-row-header"
            ),
            DOM.div(sl_c; class="c-slider-wrap");
            class="c-row"
        ),
        DOM.div(
            DOM.span("Initial state"; class="init-label"),
            DOM.div(
                DOM.div(thumb_zeros,  DOM.span("all zeros"; class="init-thumb-lbl"); class="init-thumb-wrap"),
                DOM.div(thumb_square, DOM.span("square";    class="init-thumb-lbl"); class="init-thumb-wrap"),
                DOM.div(thumb_random, DOM.span("random";    class="init-thumb-lbl"); class="init-thumb-wrap");
                class="init-thumbs"
            );
            class="init-section"
        ),
        DOM.div(
            DOM.div(btn_generate; class="generate-btn"),
            DOM.div(btn_abort;    class="abort-btn");
            class="btn-row"
        ),
        progress_track_node,
        DOM.div(status_obs; class="status-text");
        class="app-settings"
    )

    nav_prev = DOM.button("\u25c0"; id="nav-prev", class="time-nav")
    nav_next = DOM.button("\u25b6"; id="nav-next", class="time-nav")

    evaljs(session, js"""
        document.getElementById('nav-prev').addEventListener('click', () => $(key_obs).notify(-1));
        document.getElementById('nav-next').addEventListener('click', () => $(key_obs).notify(1));
    """)

    time_controls = DOM.div(
        nav_prev,
        DOM.div(sl_time; class="time-slider"),
        nav_next,
        DOM.span(time_label; class="time-label");
        id="time-controls",
        class="time-controls"
    )

    evaljs(session, js"""
        document.title = 'ML Single-Shot';
        setTimeout(() => { document.title = 'ML Single-Shot'; }, 500);
        // Ensure proper scaling on mobile — without this, mobile browsers
        // shrink the viewport which makes the WebGL canvas render incorrectly
        if (!document.querySelector('meta[name=viewport]')) {
            const m = document.createElement('meta');
            m.name = 'viewport';
            m.content = 'width=device-width, initial-scale=1';
            document.head.appendChild(m);
        }
    """)

    explore_banner = DOM.div(
        DOM.span(explore_status_obs);
        id="explore-banner", class="banner banner-explore"
    )
    onjs(session, explore_status_obs, js"""(msg) => {
        const banner = document.getElementById('explore-banner');
        const tc     = document.getElementById('time-controls');
        if (msg && msg.length > 0) {
            if (banner) banner.style.display = 'block';
            if (tc)     tc.classList.add('explore-active');
        } else {
            if (banner) banner.style.display = 'none';
            if (tc)     tc.classList.remove('explore-active');
        }
    }""")

    params_banner = DOM.div(
        DOM.span(run_params_obs);
        id="params-banner", class="banner banner-params"
    )
    onjs(session, run_params_obs, js"""(msg) => {
        const banner = document.getElementById('params-banner');
        if (banner) banner.style.display = (msg && msg.length > 0) ? 'block' : 'none';
    }""")

    DOM.div(
        style_tag,
        DOM.div(
            DOM.div(
                DOM.span("2D repetition code"),
                DOM.span(" \u00b7 approximate maximum-likelihood single-shot decoder"; class="subtitle");
                class="app-title"),
            DOM.div(
                DOM.div(param_panel; class="app-settings-col"),
                DOM.div(
                    params_banner,
                    DOM.div(fig; class="makie-wrap"),
                    time_controls,
                    explore_banner;
                    class="app-viz");
                class="app-layout"
            );
            class="app-card"
        )
    )
end

# =============================================================
# Start server
# =============================================================

# server = Bonito.Server(app, "0.0.0.0", 8000)
server = Bonito.Server(app, "0.0.0.0", 8000;
    proxy_url="https://ml-single-shot.mythb.ee")
@info "▶  Server running →  http://localhost:8000"
wait(server)