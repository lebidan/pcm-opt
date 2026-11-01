# Trapping-set distribution 𝒯(a,b) on the bipartite Tanner graph of a GF2Matrix (rows = check nodes,
# columns = variable nodes). An (a,b) trapping set (Def. 1 of the JSAC paper / [50]) is a set of `a`
# VNs; `b` is the number of odd-degree (unsatisfied) CNs in the subgraph it induces.
#
# We count the dominant trapping sets 𝒯 as defined in Karimi & Banihashemi [50, §II]: the induced
# subgraph is **connected**, **contains a cycle** (their Lemma 6), and **every VN is connected to ≥2
# satisfied (even-degree) CNs**, with degree-2 VNs exempt (their relaxation for irregular codes). This
# excludes tree-like sets (e.g. two VNs sharing one CN) — matching the paper's Table I, where girth-6
# PCMs have no 𝒯(2,·) and the smallest trapping set is a 6-cycle 𝒯(3,·). The enumeration below is
# exhaustive. On the RM(32,16) PCM-II it reproduces the paper's counts exactly (𝒯(3,12)=896,
# 𝒯(4,12)=1088). The search method of [50] is not exhaustive on dense girth-4 graphs (e.g. a standard
# PCM), so the paper's counts for those may differ.
#
# Connected subsets are enumerated once each by the ESU reverse-search (Wernicke): each set is rooted
# at its minimum-index VN and grown only by exclusive neighbours of larger index. Enumeration is
# bounded by `a_max` (default 4, matching the paper's tables).

# Smallest trapping set present (minimize a, then b); `nothing` if the matrix has none.
min_trapping_set(dist::AbstractDict) = isempty(dist) ? nothing : minimum(keys(dist))

function trapping_distribution(A::GF2Matrix; a_max::Int = 4)
    a_max >= 2 || throw(ArgumentError("a_max must be ≥ 2"))
    m, n = nrows(A), ncols(A)

    # per-VN incident CNs, and VN adjacency (two VNs adjacent iff they share a CN)
    vn_cns = [Int[] for _ in 1:n]
    @inbounds for i in 1:m
        r = A.rows[i]
        for j in 1:n
            r[j] && push!(vn_cns[j], i)
        end
    end
    adj = [Int[] for _ in 1:n]
    adjset = [Set{Int}() for _ in 1:n]
    @inbounds for i in 1:m
        vs = findall(A.rows[i])
        for a in 1:length(vs)-1, b in a+1:length(vs)
            x, y = vs[a], vs[b]
            (y in adjset[x]) || (push!(adjset[x], y); push!(adj[x], y))
            (x in adjset[y]) || (push!(adjset[y], x); push!(adj[y], x))
        end
    end

    dist = Dict{Tuple{Int,Int},Int}()
    inc = zeros(Int, m)   # per-CN incidence with the current subset
    bcount = Ref(0)       # #CNs with odd incidence = b
    ccount = Ref(0)       # #CNs incident to the subset (induced check nodes)
    edeg = Ref(0)         # Σ VN degrees over the subset = #edges of the induced subgraph
    Vsub = Int[]

    inc_cn!(i) = (inc[i] += 1; inc[i] == 1 && (ccount[] += 1); bcount[] += isodd(inc[i]) ? 1 : -1)
    dec_cn!(i) = (bcount[] -= isodd(inc[i]) ? 1 : -1; inc[i] -= 1; inc[i] == 0 && (ccount[] -= 1))
    push_vn!(w) = (edeg[] += length(vn_cns[w]); @inbounds for i in vn_cns[w]; inc_cn!(i); end)
    pop_vn!(w) = (edeg[] -= length(vn_cns[w]); @inbounds for i in vn_cns[w]; dec_cn!(i); end)

    # `u` is exclusive to the current Vsub when it is neither in Vsub nor adjacent to any of its VNs.
    function is_exclusive(u)
        @inbounds for s in Vsub
            (u == s || u in adjset[s]) && return false
        end
        true
    end

    # [50] dominance: every VN (except degree-2 ones) touches ≥2 satisfied (even-degree) CNs.
    function is_dominant()
        @inbounds for v in Vsub
            length(vn_cns[v]) == 2 && continue
            sat = 0
            for cn in vn_cns[v]
                iseven(inc[cn]) && (sat += 1)
                sat >= 2 && break
            end
            sat >= 2 || return false
        end
        true
    end

    # Induced subgraph contains a cycle: #edges ≥ #vertices (|Vsub| VNs + ccount CNs). [50]'s Lemma 6
    # guarantees every dominant trapping set has a cycle; the degree-2 exemption in is_dominant must
    # not admit tree-like sets, so this is enforced explicitly.
    has_cycle() = edeg[] >= length(Vsub) + ccount[]

    function extend!(ext::Vector{Int}, root::Int)
        length(Vsub) >= 2 && has_cycle() && is_dominant() &&
            (dist[(length(Vsub), bcount[])] = get(dist, (length(Vsub), bcount[]), 0) + 1)
        length(Vsub) == a_max && return
        while !isempty(ext)
            w = pop!(ext)
            child = copy(ext)   # siblings no longer see w; child also gets w's exclusive neighbours
            @inbounds for u in adj[w]
                (u > root && is_exclusive(u) && !(u in child)) && push!(child, u)
            end
            push_vn!(w)
            push!(Vsub, w)
            extend!(child, root)
            pop!(Vsub)
            pop_vn!(w)
        end
    end

    for v in 1:n
        isempty(vn_cns[v]) && continue      # isolated VN (all-zero column)
        push_vn!(v)
        push!(Vsub, v)
        ext = Int[]
        @inbounds for u in adj[v]
            u > v && push!(ext, u)
        end
        extend!(ext, v)
        pop!(Vsub)
        pop_vn!(v)
    end

    dist
end
