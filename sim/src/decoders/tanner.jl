# Store Tanner graph using CSR (compressed sparse row) representation, for
# contiguous, cache-friendly access during BP decoding.
# Edge labeling uses Int32 to further optimize memory/cache access.
# Int32 is sufficient to handle up to 2 billions edges.
#
# Usage:
# - edges connected to a variable v: vn_edge[vn_ptr[v]:vn_ptr[v+1]-1]
# - edges connected to a check c: cn_edge[cn_ptr[c]:cn_ptr[c+1]-1]
# - variable nodes connected to a check c: cn_vn[cn_ptr[c]:cn_ptr[c+1]-1]
#
# Concrete example. Take a tiny H (3 checks × 4 variables):
#
#         v1 v2 v3 v4
#   c1:    1  1  0  0
#   c2:    0  1  1  1
#   c3:    1  0  1  0
#
# Building the graph gives n_edges=7, and:
#
#   vn_ptr  = [1, 3, 5, 7, 8]           # v1 owns slots 1:2, v2 owns 3:4, v3 owns 5:6, v4 owns 7:7
#   vn_edge = [1, 2, 3, 4, 5, 6, 7]
#   cn_ptr  = [1, 3, 6, 8]              # c1 owns slots 1:2, c2 owns 3:5, c3 owns 6:7
#   cn_edge = [1, 3, 4, 5, 7, 2, 6]     # edge id at each slot
#   cn_vn   = [1, 2, 2, 3, 4, 1, 3]     # variable node at each slot

struct TannerGraph
    n_vn::Int
    n_cn::Int
    n_edges::Int
    dc_max::Int
    cn_ptr::Vector{Int}     # n_cn+1 offsets (1-based) into cn_edge/cn_vn
    cn_edge::Vector{Int32}  # E: edge ids grouped by check
    cn_vn::Vector{Int32}    # E: the variable node of each edge from a check perspective
    vn_ptr::Vector{Int}     # n_vn+1 offsets into vn_edge
    vn_edge::Vector{Int32}  # E: edge ids grouped by variable
end

# Edges are numbered column by column (variable node by variable node), so the edges of each
# variable node are consecutive and vn_edge is simply 1:E.
function TannerGraph(H::AbstractMatrix{Bool})
    n_cn, n_vn = size(H)
    dv, dc = vec(count(H; dims = 1)), vec(count(H; dims = 2))
    # The check-node update computes extrinsic messages from the other edges, so it needs degree ≥ 2.
    c = findfirst(<(2), dc)
    c === nothing || throw(ArgumentError("check $c has degree $(dc[c]); BP needs every check degree ≥ 2"))
    vn_ptr, cn_ptr = cumsum([1; dv]), cumsum([1; dc])
    E = vn_ptr[end] - 1
    cn_edge, cn_vn = Vector{Int32}(undef, E), Vector{Int32}(undef, E)
    cpos = copy(cn_ptr)   # next free slot of each check
    e = 0
    for v in 1:n_vn, c in 1:n_cn
        @inbounds H[c, v] || continue
        e += 1
        cn_edge[cpos[c]] = e
        cn_vn[cpos[c]] = v
        cpos[c] += 1
    end
    return TannerGraph(n_vn, n_cn, E, maximum(dc), cn_ptr, cn_edge, cn_vn, vn_ptr, collect(Int32, 1:E))
end
