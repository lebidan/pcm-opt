@enum LdpcSchedule Layered Flooding
@enum CNRule       SPAExact FastSPA MinSum NormMinSum   # check-node "boxplus" rule

struct BPConfig
    schedule::LdpcSchedule
    max_iter::Int
    cn_rule::CNRule
    norm_factor::LLRFloat        # normalized-min-sum scaling α (used only by NormMinSum)
    llr_clip::LLRFloat
    function BPConfig(schedule, max_iter, cn_rule, norm_factor, llr_clip)
        norm, clip = LLRFloat(norm_factor), LLRFloat(llr_clip)
        isfinite(norm) && norm > 0 || throw(ArgumentError("norm_factor must be finite and positive in LLRFloat"))
        isfinite(clip) && clip > 0 || throw(ArgumentError("llr_clip must be finite and positive in LLRFloat"))
        new(schedule, max_iter, cn_rule, norm, clip)
    end
end
# Convenience constructor: α defaults to the paper's 0.7 (irrelevant for non-NMS rules).
BPConfig(schedule, max_iter, cn_rule) = BPConfig(schedule, max_iter, cn_rule, 0.7f0, 128f0)
BPConfig(schedule, max_iter, cn_rule, norm_factor) = BPConfig(schedule, max_iter, cn_rule, norm_factor, 128f0)

# One decoder per thread: the graph is shared read-only, the message buffers are its own.
struct BPDecoder
    graph::TannerGraph
    cfg::BPConfig
    llr_out::Vector{LLRFloat}
    msg_c2v::Vector{LLRFloat}
    msg_v2c::Vector{LLRFloat}
    α::Vector{LLRFloat}
    β::Vector{LLRFloat}
end

BPDecoder(g::TannerGraph, cfg::BPConfig) = BPDecoder(g, cfg, zeros(LLRFloat, g.n_vn),
    zeros(LLRFloat, g.n_edges), zeros(LLRFloat, g.n_edges), zeros(LLRFloat, g.dc_max), zeros(LLRFloat, g.dc_max))

decoder_name(cfg::BPConfig) = cfg.cn_rule == SPAExact ? "SPA" : cfg.cn_rule == FastSPA ? "FastSPA" :
    cfg.cn_rule == NormMinSum ? "NMS" : "MS"
schedule_name(s::LdpcSchedule) = lowercase(string(s))

# Human (label => value) pairs for the console info header.
function decoder_params(cfg::BPConfig)
    params = [
        "Algorithm"           => decoder_name(cfg),
        "Schedule"            => schedule_name(cfg.schedule),
        "Max. iterations (i)" => string(cfg.max_iter),
        "LLR clip"            => string(cfg.llr_clip),
    ]
    cfg.cn_rule == NormMinSum && push!(params, "Norm. factor (α)" => string(cfg.norm_factor))
    return params
end

# Variable-node update: sum of the channel LLR and incoming messages, extrinsic outputs.
function process_vn!(v::Int, vn_edge::Vector{Int32}, lo::Int, hi::Int, llr_in::Vector{LLRFloat}, msg_c2v::Vector{LLRFloat}, dec_bits::Vector{Bool}, llr_out::Vector{LLRFloat}, msg_v2c::Vector{LLRFloat}, clip::LLRFloat)
    acc = llr_in[v]
    @inbounds @simd for idx in lo:hi
        acc += msg_c2v[vn_edge[idx]]
    end
    isfinite(acc) || throw(OverflowError("BP variable-node sum overflowed; reduce llr_clip"))
    llr_out[v] = acc
    dec_bits[v] = acc < 0
    @inbounds @simd for idx in lo:hi
        e = vn_edge[idx]
        msg_v2c[e] = clamp(acc - msg_c2v[e], -clip, clip)
    end
end

@inline my_sign(x::LLRFloat) = x < 0 ? -1 : +1
@inline my_min(x::LLRFloat, y::LLRFloat) = x < y ? x : y
@inline boxplus_ms(x::LLRFloat, y::LLRFloat) = my_sign(x) * my_sign(y) * my_min(abs(x), abs(y))

# Exact CN "boxplus" operator
@inline function boxplus_spa_exact(x::LLRFloat, y::LLRFloat)
    absx, absy = abs(x), abs(y)
    sum_abs, diff_abs = absx + absy, abs(absx - absy)
    x1, x2 = log(1 + exp(-sum_abs)), log(1 + exp(-diff_abs))
    return my_sign(x) * my_sign(y) * (my_min(absx, absy) + x1 - x2)
end

# Richter's single-line approximation of the CN "boxplus" operator
@inline function boxplus_spa_fast(x::LLRFloat, y::LLRFloat)
    absx, absy = abs(x), abs(y)
    sum_abs, diff_abs = absx + absy, abs(absx - absy)
    x1 = sum_abs > 2.5f0 ? 0f0 : 0.6f0 - 0.24f0 * sum_abs
    x2 = diff_abs > 2.5f0 ? 0f0 : 0.6f0 - 0.24f0 * diff_abs
    return my_sign(x) * my_sign(y) * (my_min(absx, absy) + x1 - x2)
end

# Forward-backward check-node update over cn_edge[lo:hi]; parametrized on the boxplus operator
# (the where {F} clause) so the compiler specializes and statically dispatches it in the inner loop.
# Each check-to-variable message is scaled once by `norm` — the normalized-min-sum factor α (norm=1
# for the other rules, an exact no-op in Float32). Scaling here, not inside boxplus, keeps α from
# compounding to the check-degree power.
function process_cn!(cn_edge::Vector{Int32}, lo::Int, hi::Int, msg_v2c::Vector{LLRFloat}, α::Vector{LLRFloat}, β::Vector{LLRFloat}, msg_c2v::Vector{LLRFloat}, boxplus::F, norm::LLRFloat, clip::LLRFloat) where {F}
    dc = hi - lo + 1
    @inbounds begin
        α[1] = msg_v2c[cn_edge[lo]]
        for i in 2:dc-1
            α[i] = boxplus(α[i-1], msg_v2c[cn_edge[lo+i-1]])
        end
        β[dc] = msg_v2c[cn_edge[hi]]
        for i in dc-1:-1:2
            β[i] = boxplus(β[i+1], msg_v2c[cn_edge[lo+i-1]])
        end
        msg_c2v[cn_edge[lo]] = clamp(norm * β[2], -clip, clip)
        for i in 2:dc-1
            msg_c2v[cn_edge[lo+i-1]] = clamp(norm * boxplus(α[i-1], β[i+1]), -clip, clip)
        end
        msg_c2v[cn_edge[hi]] = clamp(norm * α[dc-1], -clip, clip)
    end
end

function is_codeword(code_bits::Vector{Bool}, g::TannerGraph)
    @inbounds for c in 1:g.n_cn
        checksum = false
        for i in g.cn_ptr[c]:g.cn_ptr[c+1]-1
            checksum ⊻= code_bits[g.cn_vn[i]]
        end
        checksum && return false
    end
    return true
end

# Reset the messages and run the initial variable-node pass from the channel LLRs alone
# (msg_c2v ≡ 0); returns true if the raw hard decision is already a valid codeword.
function init_pass!(dec::BPDecoder, llr_in::Vector{LLRFloat}, dec_bits::Vector{Bool})
    g = dec.graph
    msg_c2v, msg_v2c, llr_out = dec.msg_c2v, dec.msg_v2c, dec.llr_out
    fill!(msg_c2v, 0)
    @inbounds for v in 1:g.n_vn
        process_vn!(v, g.vn_edge, g.vn_ptr[v], g.vn_ptr[v+1]-1, llr_in, msg_c2v, dec_bits, llr_out, msg_v2c, dec.cfg.llr_clip)
    end
    return is_codeword(dec_bits, g)
end

# Returns (result, iterations): 0 iterations if the channel hard decision is already a codeword,
# max_iter on a Failure.
function decode!(dec::BPDecoder, llr_in::Vector{LLRFloat}, dec_bits::Vector{Bool})
    all(isfinite, llr_in) || throw(ArgumentError("channel LLRs must be finite"))
    cfg = dec.cfg
    boxplus = cfg.cn_rule == SPAExact ? boxplus_spa_exact :
              cfg.cn_rule == FastSPA ? boxplus_spa_fast : boxplus_ms   # MinSum & NormMinSum
    norm = cfg.cn_rule == NormMinSum ? cfg.norm_factor : one(LLRFloat)
    return cfg.schedule == Flooding ?
        flooding_bp!(boxplus, norm, dec, llr_in, dec_bits, cfg.max_iter) :
        layered_bp!(boxplus, norm, dec, llr_in, dec_bits, cfg.max_iter)
end

function layered_bp!(boxplus::FBP, norm::LLRFloat, dec::BPDecoder, llr_in::Vector{LLRFloat}, dec_bits::Vector{Bool}, max_iter) where {FBP}
    g = dec.graph
    msg_c2v, msg_v2c, llr_out = dec.msg_c2v, dec.msg_v2c, dec.llr_out
    α, β = dec.α, dec.β
    init_pass!(dec, llr_in, dec_bits) && return (Codeword, 0)
    for it in 1:max_iter
        @inbounds for c in 1:g.n_cn
            lo, hi = g.cn_ptr[c], g.cn_ptr[c+1]-1
            process_cn!(g.cn_edge, lo, hi, msg_v2c, α, β, msg_c2v, boxplus, norm, dec.cfg.llr_clip)
            for i in lo:hi
                v = Int(g.cn_vn[i])
                process_vn!(v, g.vn_edge, g.vn_ptr[v], g.vn_ptr[v+1]-1, llr_in, msg_c2v, dec_bits, llr_out, msg_v2c, dec.cfg.llr_clip)
            end
        end
        is_codeword(dec_bits, g) && return (Codeword, it)
    end
    return (Failure, max_iter)
end

function flooding_bp!(boxplus::FBP, norm::LLRFloat, dec::BPDecoder, llr_in::Vector{LLRFloat}, dec_bits::Vector{Bool}, max_iter) where {FBP}
    g = dec.graph
    msg_c2v, msg_v2c, llr_out = dec.msg_c2v, dec.msg_v2c, dec.llr_out
    α, β = dec.α, dec.β
    init_pass!(dec, llr_in, dec_bits) && return (Codeword, 0)
    for it in 1:max_iter
        @inbounds for c in 1:g.n_cn
            process_cn!(g.cn_edge, g.cn_ptr[c], g.cn_ptr[c+1]-1, msg_v2c, α, β, msg_c2v, boxplus, norm, dec.cfg.llr_clip)
        end
        @inbounds for v in 1:g.n_vn
            process_vn!(v, g.vn_edge, g.vn_ptr[v], g.vn_ptr[v+1]-1, llr_in, msg_c2v, dec_bits, llr_out, msg_v2c, dec.cfg.llr_clip)
        end
        is_codeword(dec_bits, g) && return (Codeword, it)
    end
    return (Failure, max_iter)
end

# Machine (key => value) pairs for the YAML sidecar's `decoder:` block.
decoder_spec(cfg::BPConfig) = [
    "algorithm" => decoder_name(cfg),
    "schedule"  => schedule_name(cfg.schedule),
    "max_iter"  => cfg.max_iter,
    "norm_factor" => cfg.norm_factor,
    "llr_clip" => cfg.llr_clip,
]
