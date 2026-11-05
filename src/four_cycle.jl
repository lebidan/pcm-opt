# Algorithm 1 — Optimized 4-Cycle Removal Strategy.
#
# Greedy: while 4-cycles remain, pick the all-ones block (rows S_c × cols S_v) that removes the most
# 4-cycles, collapse its shared columns by adding one AVN + one ACN (sub-matrix (9)), then run the
# degree-1 and degree-2 CN avoidance cleanups over the whole matrix. Blocks are ranked by
# (δ = C(|S_v|,2)·C(|S_c|,2), |S_c|, |S_v|) with a ≥ tie-break (last in scan order wins). On RM(32,16)
# this gives Table I's PCM-I size and edges with 39 6-cycles instead of 38; the randomized tie-break
# with seed 1815 gives 38.
#
# Cleanups (Alg.1 lines 17-23): a degree-1 CN on an AVN forces that AVN to 0 — drop the AVN column
# and the CN. An AVN-related degree-2 CN makes its two VNs equal — merge them (keeping a real VN
# when present), dropping the AVN column and the now-trivial CN. Both shrink the matrix and remove
# the un-informative degree-1/2 CNs that seed small trapping sets. They are repeated until neither
# applies, so the result has no degree-1 CN (BP decoders cannot process one).
#
# The original N columns are never deleted or moved; AVN columns are appended at the end and only
# AVN columns are deleted, so "column j is an AVN" ⟺ j > N throughout.

# The all-ones block seeded by row pair (i1,i2): shared columns S_v (with |S_v| = t) and the rows S_c
# that contain them. Returns nothing if the pair overlaps in <2 columns (no 4-cycle to remove).
@inline function pair_block(A::GF2Matrix, i1::Int, i2::Int)
    t = overlap(A, i1, i2)
    t > 1 || return nothing
    sv = A.rows[i1] .& A.rows[i2]
    sc = Int[]
    @inbounds for r in 1:nrows(A)
        subset_of(sv, A.rows[r]) && push!(sc, r)
    end
    (sv, sc, t)
end

function select_block(A::GF2Matrix)
    m = nrows(A)
    bestkey = (-1, -1, -1)
    Sv = falses(A.ncols)
    Sc = Int[]
    for i1 in 1:m-1, i2 in i1+1:m
        b = pair_block(A, i1, i2)
        b === nothing && continue
        sv, sc, t = b
        s = length(sc)
        key = (binom2(t) * binom2(s), s, t)
        key >= bestkey && ((bestkey, Sv, Sc) = (key, sv, sc))
    end
    Sv, Sc
end

# Remove any degree-1 CN sitting on an AVN (column > N), which forces that AVN to 0.
function cleanup_degree1!(A::GF2Matrix, N::Int)
    changed = true
    while changed
        changed = false
        for r in 1:nrows(A)
            count(A.rows[r]) == 1 || continue
            j = findfirst(A.rows[r])
            j > N || continue
            deleterow!(A, r)
            deletecol!(A, j)
            changed = true
            break
        end
    end
    A
end

drop_zero_rows!(A::GF2Matrix) = (filter!(any, A.rows); A)

# Merge the two VNs of any AVN-related degree-2 CN. Columns are in ascending order and AVNs come
# last, so the AVN is v_j2 and the kept VN v_j1 is a real one whenever the CN has one.
function cleanup_degree2!(A::GF2Matrix, N::Int)
    changed = true
    while changed
        changed = false
        for r in 1:nrows(A)
            count(A.rows[r]) == 2 || continue
            j1, j2 = findall(A.rows[r])
            j2 > N || continue
            @inbounds for row in A.rows        # substitute v_j2 := v_j1
                row[j2] && (row[j1] = !row[j1])
            end
            deleterow!(A, r)
            deletecol!(A, j2)
            changed = true
            break
        end
    end
    A
end

# Remove one all-ones block: add its AVN + ACN, XOR the ACN into the block rows, then clean up.
function apply_block!(A::GF2Matrix, Sv::BitVector, Sc::Vector{Int}, N::Int)
    addcol!(A)
    addrow!(A, push!(copy(Sv), true))   # ACN: the block's columns plus the new AVN
    mrow = nrows(A)
    @inbounds for i in Sc
        xorrow!(A, i, mrow)
    end
    cleanup!(A, N)
end

# Run both cleanups until neither applies. Each can create work for the other: a degree-2 merge can
# leave a degree-1 CN on an AVN, and dropping that AVN can leave a new degree-2 CN.
function cleanup!(A::GF2Matrix, N::Int)
    while true
        sz = (nrows(A), ncols(A))
        cleanup_degree1!(A, N)
        cleanup_degree2!(A, N)
        drop_zero_rows!(A)      # merge identical rows (a duplicate collapses to an all-zero row)
        (nrows(A), ncols(A)) == sz && return A
    end
end

# All maximal-δ blocks (δ = 4-cycles removed = C(|S_v|,2)·C(|S_c|,2)), for randomized tie-breaking.
function max_delta_blocks(A::GF2Matrix)
    m = nrows(A)
    bestδ = -1
    pool = Tuple{BitVector,Vector{Int}}[]
    for i1 in 1:m-1, i2 in i1+1:m
        b = pair_block(A, i1, i2)
        b === nothing && continue
        sv, sc, t = b
        δ = binom2(t) * binom2(length(sc))
        if δ > bestδ
            bestδ = δ
            empty!(pool)
            push!(pool, (sv, sc))
        elseif δ == bestδ
            push!(pool, (sv, sc))
        end
    end
    pool
end

# Operates on a copy of A0 (assumed AVN-free: ncols == N original columns) and returns the
# 4-cycle-free result, its AVN columns at indices N+1…Ñ.
function remove_4cycles(A0::GF2Matrix; max_iter::Int = 100_000)
    A = copy(A0)
    N = ncols(A)
    for _ in 1:max_iter
        Sv, Sc = select_block(A)            # empty S_c ⟺ no 4-cycles remain
        isempty(Sc) && return A
        apply_block!(A, Sv, Sc, N)
    end
    error("remove_4cycles did not converge in $max_iter iterations")
end

# As remove_4cycles, but breaks ties randomly among the maximal-δ blocks. Different rngs explore
# different greedy paths, which can reach sparser 6-cycle structures than the deterministic rule.
function remove_4cycles_rand(A0::GF2Matrix, rng::AbstractRNG; max_iter::Int = 100_000)
    A = copy(A0)
    N = ncols(A)
    for _ in 1:max_iter
        pool = max_delta_blocks(A)          # empty ⟺ no 4-cycles remain
        isempty(pool) && return A
        Sv, Sc = rand(rng, pool)
        apply_block!(A, Sv, Sc, N)
    end
    error("remove_4cycles_rand did not converge in $max_iter iterations")
end

# Algorithm 4 — Rapid Elimination of 4-Cycles for one candidate row.
# Predict whether candidate `hk` (length N) can be appended to the already-4-cycle-free `Hb_prev`
# and made 4-cycle-free again *without* adding an AVN, by preventing degree-2 CNs: shrink the padded
# candidate by XORing in any existing row it forms a 4-cycle with and that differs from the overlap
# by a single VN, then XOR the reduced candidate into any existing row it now almost-fits into.
# Returns the new matrix (M̃_prev+1 rows, same columns) on success, or nothing (caller runs Alg.1).
function rapid_bypass(Hb_prev::GF2Matrix, hk::BitVector, N::Int)
    A = copy(Hb_prev)
    addrow!(A, [hk; falses(A.ncols - N)])
    k = nrows(A)
    count_4cycles(A) == 0 && return A
    hp = A.rows[k]
    @inbounds for i in 1:k-1
        sh = and_count(A.rows[i], hp)
        (sh >= 2 && count(A.rows[i]) - sh == 1) && xor_bits!(hp, A.rows[i])
    end
    whp = count(hp)
    @inbounds for i in 1:k-1
        sh = and_count(A.rows[i], hp)
        (sh >= 2 && whp - sh == 1) && xor_bits!(A.rows[i], hp)
    end
    (count_4cycles(A) == 0 && all(any, A.rows)) ? A : nothing
end
