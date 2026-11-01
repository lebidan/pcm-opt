# Algorithm 3 — Progressive Selection of Parity-Check Rows (PRG).
#
# Grow H_a one row at a time from H_s (up to M_a rows). At each step every remaining candidate row is
# appended to H_a and 4-cycle-removed (Algorithm 1) to score it; the winner is chosen by the paper's
# criteria and committed. The best 4-cycle-free H_b over all steps i ≥ N−K (largest η̄_min) is
# returned — a single M_a-length run sweeps the row-size axis of Fig. 3 and picks its peak.
#
# Criteria: (1) while i ≤ N−K a candidate must raise the rank of H_a; (2) keep only candidates whose
# removal added equal ACNs and AVNs (M̃_b − i == Ñ_b − N); then (3a) minimize M̃_b, (3b) maximize the
# variance of H_a's column weights, (3c) break ties at random.

# Forward-reduction basis for rank tracking: reduce h against the echelon basis; the residual is
# nonzero iff h is independent of it.
function residual(basis::Vector{Tuple{Int,BitVector}}, h::BitVector)
    r = copy(h)
    for (pc, pv) in basis
        r[pc] && xor_bits!(r, pv)
    end
    r
end

is_independent(basis::Vector{Tuple{Int,BitVector}}, h::BitVector) =
    findfirst(residual(basis, h)) !== nothing

function basis_add!(basis::Vector{Tuple{Int,BitVector}}, h::BitVector)
    r = residual(basis, h)
    pc = findfirst(r)
    pc === nothing && return false
    push!(basis, (pc, r))
    true
end

# Variance of H_a's column weights if candidate row `hk` were appended (cw = current column weights).
function colweight_var(cw::Vector{Int}, hk::BitVector)
    N = length(cw)
    s = 0
    s2 = 0
    @inbounds for j in 1:N
        v = cw[j] + (hk[j] ? 1 : 0)
        s += v
        s2 += v * v
    end
    s2 / N - (s / N)^2
end

struct PRGResult
    Hb::GF2Matrix
    i::Int              # rows of H_a at the selected step
    eta_min::Int
end

# Pick the winner by column-weight variance (3b) then random (3c) among a candidate subset.
function pick_3bc(R::Vector{Int}, cand::Vector{Int}, Hs::GF2Matrix, cw::Vector{Int}, rng)
    vars = [colweight_var(cw, Hs.rows[cand[c]]) for c in R]
    vmax = maximum(vars)
    R = [R[t] for t in eachindex(R) if vars[t] == vmax]
    R[rand(rng, 1:length(R))]
end

function extract_pcm(Hs::GF2Matrix, N::Int, K::Int; M_a::Int = nrows(Hs), seed::Integer = 0,
                     rapid::Bool = false)
    target = N - K
    Ma = min(M_a, nrows(Hs))
    rng = Xoshiro(hash((seed, :prg)))

    Ha = GF2Matrix([copy(Hs.rows[1])], N)
    Hb_prev = rapid ? remove_4cycles(Ha) : Ha
    avail = collect(2:nrows(Hs))
    basis = Tuple{Int,BitVector}[]
    basis_add!(basis, Hs.rows[1])
    cw = zeros(Int, N)
    cw .+= Hs.rows[1]

    best_Hb, best_i, best_eta = Ha, 0, typemin(Int)
    i = 1
    while i < Ma && !isempty(avail)
        i += 1
        cand = i <= target ? filter(k -> is_independent(basis, Hs.rows[k]), avail) : avail
        isempty(cand) && break

        chosen = 0
        Hb_new = nothing

        # Algorithm 4: try to bypass full 4-cycle removal (add a row without a new AVN).
        if rapid
            bypass = Vector{Union{Nothing,GF2Matrix}}(nothing, length(cand))
            Threads.@threads for c in eachindex(cand)
                bypass[c] = rapid_bypass(Hb_prev, Hs.rows[cand[c]], N)
            end
            Rby = [c for c in eachindex(cand) if bypass[c] !== nothing]
            if !isempty(Rby)
                c = pick_3bc(Rby, cand, Hs, cw, rng)
                chosen, Hb_new = cand[c], bypass[c]
            end
        end

        # Full Algorithm 1 per candidate (exact mode, or rapid with no bypass available).
        if Hb_new === nothing
            Mb = Vector{Int}(undef, length(cand))
            Nb = Vector{Int}(undef, length(cand))
            Threads.@threads for c in eachindex(cand)
                Hbk = remove_4cycles(GF2Matrix(vcat(Ha.rows, [Hs.rows[cand[c]]]), N))
                Mb[c] = nrows(Hbk)
                Nb[c] = ncols(Hbk)
            end
            R = [c for c in eachindex(cand) if Mb[c] - i == Nb[c] - N]   # crit 2
            isempty(R) && break
            minMb = minimum(Mb[c] for c in R)                            # crit 3a
            R = [c for c in R if Mb[c] == minMb]
            c = pick_3bc(R, cand, Hs, cw, rng)                           # crit 3b, 3c
            chosen = cand[c]
            Hb_new = remove_4cycles(GF2Matrix(vcat(Ha.rows, [Hs.rows[chosen]]), N))
        end

        addrow!(Ha, copy(Hs.rows[chosen]))
        basis_add!(basis, Hs.rows[chosen])
        cw .+= Hs.rows[chosen]
        deleteat!(avail, findfirst(==(chosen), avail))
        Hb_prev = Hb_new

        if i >= target
            # η = −1 (no AVN on a 6-cycle) ranks below any real η̄_min here, while scripts/pcm_i.jl
            # ranks it first, as if η̄_min were infinite.
            η = eta_min(Hb_new, N)
            if η > best_eta
                best_Hb, best_i, best_eta = Hb_new, i, η
            end
        end
    end
    best_eta == typemin(Int) &&
        error("extract_pcm reached no full-rank PCM (needs M_a ≥ N−K = $target and a spanning H_s)")
    PRGResult(best_Hb, best_i, best_eta)
end
