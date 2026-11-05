# Structured sparse PCM (ssPCM) built from parity-check row blocks (PCRBs).
#
# Y. Shen, Z. Li, E. Boutillon, W. Song, Y. Ren, C. Zhang, X. You, A. Burg, "Belief propagation
# decoding for short codes on structured sparse parity-check matrices," IEEE ISIT 2025, Sec. III,
# doi:10.1109/ISIT63088.2025.11195685.
#
# A PCRB is a set of s ≥ 2 dual codewords that share a column set T and whose residuals h \ T are
# pairwise disjoint. One AVN v' replaces T: the block has the rows (h \ T) + v' and the ACN T + v'.
# Inside a block every coded column has weight ≤ 1, so the block is cycle-free, and AVNs of
# different blocks share no CN. Rows that pairwise share ≤ 1 column form an AVN-free block.
#
# The paper's greedy picks the largest block first. Two changes make it work beyond BCH codes:
#  - rank first: until the PCM has full rank, only rank-increasing blocks are considered;
#  - Γ cap: Γmax, the largest number of CNs shared by two coded columns, stays ≤ a cap that starts
#    at 1 and is raised by one only when no block fits. This follows Theorem 1 of Y. Shen,
#    H. S. Ünal, A. Burg, "HDPC codes with LDPC matrices: construction based on social golfer
#    problem," IEEE ISIT 2026, doi:10.1109/ISIT62367.2026.11653623: minimize Γmax, then maximize
#    the VN degrees, i.e. add as many blocks as possible.
# Ties: largest block, then fewest new 4-cycles, then smallest variance of the coded column weights,
# then at random. Breaking ties by index order instead picks structured blocks that can fail to
# reach full rank (RM(128,64)), so the seed matters: try a few.

struct PCRB
    rows::Vector{Int}                 # search-space rows used
    T::Union{BitVector,Nothing}       # columns replaced by the AVN (nothing: AVN-free block)
end

# Coded parts of the block's rows (the ACN T comes last).
coded_rows(S, b::PCRB) = b.T === nothing ? S[b.rows] : [[S[k] .& .!b.T for k in b.rows]; [b.T]]

# True if adding row r keeps every coded column pair within `cap` shared CNs. O[a,b] counts them.
function fits(r::BitVector, O::Matrix{Int}, cap::Int)
    idx = findall(r)
    @inbounds for a in eachindex(idx), b in a+1:length(idx)
        O[idx[a], idx[b]] >= cap && return false
    end
    true
end

# 4-cycles that row r would add: Σ over its column pairs of the CNs already sharing that pair.
function new_4cycles(r::BitVector, O::Matrix{Int})
    idx, q = findall(r), 0
    @inbounds for a in eachindex(idx), b in a+1:length(idx)
        q += O[idx[a], idx[b]]
    end
    q
end

# Insertion-ordered elimination basis: v reduces to zero against it iff v is in the span.
function reduce_basis!(v::BitVector, basis)
    for b in basis
        v[findfirst(b)] && xor_bits!(v, b)
    end
    v
end
function rank_gain(basis, rows)
    loc = BitVector[]
    for r in rows
        v = reduce_basis!(reduce_basis!(copy(r), basis), loc)
        any(v) && push!(loc, v)
    end
    loc
end

# Candidate anchors of one weight class: every distinct overlap T (|T| ≥ tmin) of two rows, with the
# class rows that strictly contain it. Sorted, so the construction is deterministic.
function anchors(S, cls::Vector{Int}, N::Int, tmin::Int)
    nt = Threads.nthreads()
    found = [Set{BitVector}() for _ in 1:nt]
    Threads.@threads :static for t in 1:nt
        for a in t:nt:length(cls), b in a+1:length(cls)
            and_count(S[cls[a]], S[cls[b]]) >= tmin && push!(found[t], S[cls[a]] .& S[cls[b]])
        end
    end
    Ts = sort!(collect(union(found...)); by = T -> T.chunks)
    post = [BitVector([S[k][j] for k in cls]) for j in 1:N]      # post[j][a]: row cls[a] has column j
    sup = Vector{Vector{Int}}(undef, length(Ts))
    Threads.@threads for i in eachindex(Ts)
        m = trues(length(cls))
        for j in findall(Ts[i])
            m .&= post[j]
        end
        sup[i] = [cls[a] for a in findall(m) if S[cls[a]] != Ts[i]]
    end
    Ts, sup
end

# The block grown from anchor T (or from row i when T === nothing) with the unused rows that keep the
# cap; `nothing` if fewer than two rows fit.
function grow(S, T, cand, used, O, cap)
    if T === nothing
        rows = Int[]
        for k in cand
            used[k] && continue
            all(m -> and_count(S[m], S[k]) <= 1, rows) && fits(S[k], O, cap) && push!(rows, k)
        end
        return isempty(rows) ? nothing : PCRB(rows, nothing)
    end
    fits(T, O, cap) || return nothing
    rows, R = Int[], falses(length(T))
    for k in cand
        used[k] && continue
        r = S[k] .& .!T
        and_count(r, R) == 0 && fits(r, O, cap) && (push!(rows, k); R .|= r)
    end
    length(rows) >= 2 ? PCRB(rows, T) : nothing
end

# Best block of one class under the cap, or nothing. need_rank: only rank-increasing blocks.
function best_block(S, cls, Ts, sup, used, O, colw, basis, cap, need_rank, avnfree, rng)
    free = avnfree ? [k for k in cls if !used[k]] : Int[]
    tasks = [[(Ts[i], sup[i]) for i in eachindex(Ts)]; [(nothing, @view free[a:end]) for a in eachindex(free)]]
    prio = rand(rng, length(tasks))           # random tie-break, drawn up front: thread-independent
    nt = Threads.nthreads()
    best = Vector{Any}(nothing, nt)
    Threads.@threads :static for t in 1:nt
        for x in t:nt:length(tasks)
            T, cand = tasks[x]
            b = grow(S, T, cand, used, O, cap)
            b === nothing && continue
            rows = coded_rows(S, b)
            need_rank && isempty(rank_gain(basis, S[b.rows])) && continue
            w = copy(colw)
            foreach(r -> (w .+= r), rows)
            m = sum(w) / length(w)
            key = (length(rows), -sum(r -> new_4cycles(r, O), rows), -sum(abs2, w .- m), prio[x])
            (best[t] === nothing || key > best[t][1]) && (best[t] = (key, b))
        end
    end
    filter!(!isnothing, best)
    isempty(best) ? nothing : last(argmax(first, best))
end

function assemble(S, blocks, N)
    navn = count(b -> b.T !== nothing, blocks)
    A, a = GF2Matrix(0, N + navn), N
    for b in blocks
        b.T === nothing || (a += 1)
        for r in coded_rows(S, b)
            row = [r; falses(navn)]
            b.T === nothing || (row[a] = true)
            addrow!(A, row)
        end
    end
    A
end

"""
    sspcm(Hs, N, K; tmin = 2, avnfree = true, max_edges = 4000, max_blocks, seed = 1, onstep)

Build ssPCMs from the search space `Hs` (low-weight dual codewords). Block by block, the matrix
reaches full rank (ssPCM-I), then keeps growing (ssPCM-II) until it exceeds `max_edges`, has
`max_blocks` blocks, or uses up the search space. `onstep(A, Γmax, B)` is called on every
full-rank matrix A with B blocks, so the caller sees the whole nested sequence. `tmin` is the
smallest anchor size. `avnfree` allows AVN-free blocks, at a cost of O(|class|²) per step.
"""
function sspcm(Hs::GF2Matrix, N::Int, K::Int; tmin::Int = 2, avnfree::Bool = true,
               max_edges::Int = 4000, max_blocks::Int = typemax(Int), seed::Integer = 1,
               onstep = (A, Γ, B) -> nothing)
    S = Hs.rows
    wt = count.(S)
    classes = [findall(==(w), wt) for w in sort(unique(wt))]
    cand = [anchors(S, cls, N, tmin) for cls in classes]
    used, blocks, basis = falses(length(S)), PCRB[], BitVector[]
    O, colw = zeros(Int, N, N), zeros(Int, N)
    edges, cap = 0, 1
    rng = Xoshiro(seed)

    function add!(b)
        push!(blocks, b)
        used[b.rows] .= true
        for r in coded_rows(S, b)
            colw .+= r
            edges += count(r) + (b.T !== nothing)
            idx = findall(r)
            for x in idx, y in idx
                x != y && (O[x, y] += 1)
            end
        end
        append!(basis, rank_gain(basis, S[b.rows]))
        length(basis) == N - K && onstep(assemble(S, blocks, N), maximum(O), length(blocks))
    end
    # Add a block from the lightest class that has one under the cap, raising the cap when none
    # has. Once the cap exceeds Γmax every block fits, so failing then means none is left.
    function step!(ci_range, need_rank)
        while true
            for ci in ci_range
                b = best_block(S, classes[ci], cand[ci]..., used, O, colw, basis, cap, need_rank, avnfree, rng)
                b === nothing || (add!(b); return true)
            end
            cap > maximum(O) && return false
            cap += 1
        end
    end

    for ci in eachindex(classes)          # phase I: rank, class by class
        lim = rank_gf2(S[reduce(vcat, classes[1:ci])], N)
        while length(basis) < lim
            step!(ci:ci, true) || break
        end
    end
    length(basis) == N - K || error("the search space spans rank $(length(basis)) < $(N - K)")
    while edges <= max_edges && length(blocks) < max_blocks && !all(used)  # phase II
        step!(eachindex(classes), false) || break
    end
    assemble(S, blocks, N)
end

"""
    represents_code(A, G, N, K)

True if the PCM `A` (coded columns 1..N, then AVNs, each row with at most one AVN) defines exactly
the code generated by `G`: rank(A) − #AVN = N − K, and every row of G, extended with the AVN values
it implies, satisfies every row of A.
"""
function represents_code(A::GF2Matrix, G::AbstractMatrix{Bool}, N::Int, K::Int)
    navn = ncols(A) - N
    rank_gf2(A) - navn == N - K || return false
    defining = [findfirst(r -> r[N+a], A.rows) for a in 1:navn]
    for g in eachrow(G)
        x = [BitVector(g); falses(navn)]
        for a in 1:navn
            x[N+a] = isodd(and_count(A.rows[defining[a]], x))   # x[N+a] is still 0 here
        end
        all(r -> iseven(and_count(r, x)), A.rows) || return false
    end
    true
end
