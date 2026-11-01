# Cycle enumeration and ACE (approximate cycle EMD) on the bipartite Tanner graph of a GF2Matrix
# (rows = check nodes, columns = variable nodes). Only 4- and 6-cycles are tracked; a girth of 8
# or more is reported as 8.

# 6-cycle count in closed form. A 6-cycle is 3 CNs {a,b,c} + 3 distinct VNs, one shared by each
# CN pair. For an unordered CN triple the number of such (v_ab, v_bc, v_ca) with the three VNs
# distinct is, by inclusion-exclusion (all pairwise VN-set intersections equal N(a)∩N(b)∩N(c)):
#     p_ab·p_bc·p_ca − t·(p_ab + p_bc + p_ca) + 2t,   t = |N(a)∩N(b)∩N(c)|.
function count_6cycles(A::GF2Matrix)
    m = nrows(A)
    m < 3 && return 0
    P = zeros(Int, m, m)
    @inbounds for i in 1:m-1, j in i+1:m
        P[i, j] = P[j, i] = overlap(A, i, j)
    end
    total = 0
    R = A.rows
    @inbounds for a in 1:m-2, b in a+1:m-1
        P[a, b] == 0 && continue
        for c in b+1:m
            pab, pbc, pca = P[a, b], P[b, c], P[c, a]
            (pbc == 0 || pca == 0) && continue
            t = and3_count(R[a], R[b], R[c])
            total += pab * pbc * pca - t * (pab + pbc + pca) + 2t
        end
    end
    total
end

# Invoke f(v_ab, v_bc, v_ca) once per 6-cycle (the three VNs). Meant for 4-cycle-free matrices,
# where each CN pair shares at most one VN, so the inner products are tiny.
function each_6cycle(f, A::GF2Matrix)
    m = nrows(A)
    R = A.rows
    @inbounds for a in 1:m-2, b in a+1:m-1
        sab = R[a] .& R[b]
        any(sab) || continue
        Sab = findall(sab)
        for c in b+1:m
            sbc = R[b] .& R[c]
            any(sbc) || continue
            sca = R[c] .& R[a]
            any(sca) || continue
            Sbc, Sca = findall(sbc), findall(sca)
            for vab in Sab, vbc in Sbc, vca in Sca
                (vab == vbc || vbc == vca || vab == vca) && continue
                f(vab, vbc, vca)
            end
        end
    end
end

function girth(A::GF2Matrix)
    count_4cycles(A) > 0 && return 4
    count_6cycles(A) > 0 && return 6
    return 8
end

# For each AVN (columns N+1..ncols), the minimum ACE (Σ_{v in cycle}(d_v − 2)) over the 6-cycles
# touching it. Unreferenced AVNs (no 6-cycle) keep typemax and are dropped by callers via eta_min.
# Requires a 4-cycle-free matrix: each_6cycle only enumerates genuine 6-cycles when each CN pair
# shares ≤1 VN.
function ace_in_avn(A::GF2Matrix, N::Int)
    count_4cycles(A) == 0 ||
        throw(ArgumentError("ace_in_avn/eta_min require a 4-cycle-free (girth ≥ 6) matrix"))
    navn = A.ncols - N
    navn <= 0 && return Int[]
    dv = colweights(A)
    best = fill(typemax(Int), navn)
    each_6cycle(A) do vab, vbc, vca
        ace = (dv[vab] - 2) + (dv[vbc] - 2) + (dv[vca] - 2)
        for v in (vab, vbc, vca)
            v > N && (best[v-N] = min(best[v-N], ace))
        end
    end
    best
end

# η̄_min: the worst (minimum) per-AVN ACE. Returns -1 when no AVN participates in any 6-cycle.
function eta_min(A::GF2Matrix, N::Int)
    best = ace_in_avn(A, N)
    seen = filter(!=(typemax(Int)), best)
    isempty(seen) ? -1 : minimum(seen)
end

# Aggregate structural metrics of a PCM whose first N columns are the original (channel) VNs.
function pcm_metrics(A::GF2Matrix, N::Int)
    m, n, e = nrows(A), ncols(A), num_edges(A)
    q4, sixcyc = count_4cycles(A), count_6cycles(A)
    (rows = m, cols = n, edges = e, density = e / (m * n), avn = n - N,
     rowwt = e / m, colwt = e / n, q4 = q4, girth = q4 > 0 ? 4 : sixcyc > 0 ? 6 : 8,
     sixcyc = sixcyc, eta_min = eta_min(A, N))
end
