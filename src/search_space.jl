# Algorithm 2 — Generation of the Search Space PCM H_s.
#
# H_s is a redundant, low-weight PCM of rank N−K, built by collecting low-weight codewords of the
# dual code (rows of any valid PCM are dual codewords). For increasing target weight d_c, a
# meet-in-the-middle search over random systematic bases finds weight-d_c dual codewords: reduce a
# randomly column-permuted H_std to systematic form [I|M], split its rows into two halves, and look
# for a sum of up to p rows from each half (at most 2p ones in the identity part) whose total weight is d_c.
#
# The paper additionally hashes the two halves on σ random columns to bucket the O(|P|·|Q|)
# comparison; here N−K is small enough that we compare all pairs directly (finding, per basis, every
# such codeword), so σ/the hash-join is left as a scaling optimization. Randomizing the basis over
# T_max trials accumulates the full set of low-weight dual codewords.

# Reduce `rows` (m×n, full row rank) so columns 1..m are the identity. Returns false if the first m
# columns are not an information set (caller resamples the permutation).
function try_systematic!(rows::Vector{BitVector}, m::Int)
    @inbounds for i in 1:m
        p = 0
        for r in i:m
            rows[r][i] && (p = r; break)
        end
        p == 0 && return false
        rows[i], rows[p] = rows[p], rows[i]
        for r in 1:m
            r != i && rows[r][i] && xor_bits!(rows[r], rows[i])
        end
    end
    true
end

# Sums of 0..p rows drawn from `idxs` (size 0 = the zero row), so a candidate codeword can be split
# flexibly across the two halves rather than requiring exactly p rows from each.
function subset_sums(rows::Vector{BitVector}, idxs, p::Int, n::Int)
    out = [falses(n)]
    ids = collect(idxs)
    for q in 1:p, comb in combinations(ids, q)
        s = copy(rows[comb[1]])
        for k in 2:q
            xor_bits!(s, rows[comb[k]])
        end
        push!(out, s)
    end
    out
end

# One random-basis trial: collect every found codeword (sum of ≤2p rows) into per-weight buckets,
# keyed in natural column order.
function search_trial!(buckets::Dict{Int,Set{BitVector}}, Hstd::GF2Matrix, m::Int, p::Int,
                       max_weight::Int, rng)
    N = Hstd.ncols
    π = randperm(rng, N)
    prows = [Hstd.rows[i][π] for i in 1:m]       # column j of prows[i] is original column π[j]
    try_systematic!(prows, m) || return
    half = m ÷ 2
    P = subset_sums(prows, 1:half, p, N)
    Q = subset_sums(prows, half+1:m, p, N)
    for pr in P, qr in Q
        w = xor_count(pr, qr)                     # weight without allocating on rejected pairs
        (1 <= w <= max_weight) || continue
        hn = falses(N)
        hn[π] = pr .⊻ qr                          # unpermute back to natural column order
        push!(get!(() -> Set{BitVector}(), buckets, w), hn)
    end
end

# All low-weight dual codewords found over T_max random-basis trials, bucketed by weight.
function pcm_gen_weight(Hstd::GF2Matrix; T_max::Int = 4000, p::Int = 2, seed::Integer = 0,
                        max_weight::Int = Hstd.ncols)
    m = nrows(Hstd)
    nt = Threads.nthreads()
    locals = [Dict{Int,Set{BitVector}}() for _ in 1:nt]
    Threads.@threads for t in 1:nt
        rng = Xoshiro(hash((seed, t)))
        for _ in 1:cld(T_max, nt)
            search_trial!(locals[t], Hstd, m, p, max_weight, rng)
        end
    end
    merged = Dict{Int,Set{BitVector}}()
    for d in locals, (w, s) in d
        union!(get!(() -> Set{BitVector}(), merged, w), s)
    end
    merged
end

# Assemble H_s: add codewords by increasing weight until the collected rows reach rank N−K. The
# minimum dual weight is at most the lightest H_std row, so we cap collection there (avoiding bloated
# high-weight buckets) and only lift the cap if the low-weight rows fail to span. `wcap` overrides
# that initial cap; codes with a sparse low-weight dual spectrum need many increments to reach a
# spanning weight, and each increment repeats the whole search.
function gen_search_space(Hstd::GF2Matrix; T_max::Int = 4000, p::Int = 2, seed::Integer = 0,
                          wcap::Union{Nothing,Int} = nothing)
    N = Hstd.ncols
    target = nrows(Hstd)
    rank_gf2(Hstd) == target || throw(ArgumentError("gen_search_space: H must have full row rank"))
    wcap = something(wcap, minimum(count(r) for r in Hstd.rows))
    1 <= wcap <= N || throw(ArgumentError("gen_search_space: wcap=$wcap outside 1..$N"))
    while true
        buckets = pcm_gen_weight(Hstd; T_max, p, seed, max_weight = wcap)
        rows = BitVector[]
        for w in sort!(collect(keys(buckets)))
            append!(rows, sort!(collect(buckets[w]); by = bv -> bv.chunks))  # deterministic order
            # Every codeword up to the cap enters H_s, so a cap set above the weight that first
            # spans buys extra redundancy. When the cap is the default (or has climbed to the
            # spanning weight), that weight is the last bucket and this is the plain rank test.
            w >= wcap && rank_gf2(rows, N) >= target && return GF2Matrix(rows, N)
        end
        wcap >= N && error("gen_search_space: dual codewords up to weight $N span only rank " *
                           "$(rank_gf2(rows, N)) < $target; increase T_max")
        wcap += 1
    end
end

