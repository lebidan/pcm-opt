mutable struct Counters
    n_cw::Int
    n_bit_err::Int          # information-bit errors (over K, via Ginv)
    n_code_bit_err::Int     # code-bit errors (over N)
    n_cw_err::Int
    n_mld_err::Int
    n_fail::Int             # decoder returned Failure (no valid codeword)
    w_cw_min::Int
    iters_sum::Int
    iters_min::Int
    iters_max::Int
end
Counters(n::Int) = Counters(0, 0, 0, 0, 0, 0, n, 0, typemax(Int), 0)

function reduce_counters(cs::Vector{Counters})
    r = Counters(typemax(Int))
    for c in cs
        r.n_cw += c.n_cw
        r.n_bit_err += c.n_bit_err
        r.n_code_bit_err += c.n_code_bit_err
        r.n_cw_err += c.n_cw_err
        r.n_mld_err += c.n_mld_err
        r.n_fail += c.n_fail
        r.w_cw_min = min(r.w_cw_min, c.w_cw_min)
        r.iters_sum += c.iters_sum
        r.iters_min = min(r.iters_min, c.iters_min)
        r.iters_max = max(r.iters_max, c.iters_max)
    end
    return r
end

safe_ratio(a, b) = b == 0 ? 0.0 : a / b

struct PointReport                     # one Eb/N0 point's aggregate counters at one instant
    ebn0::Float64
    snr::Float64
    c::Counters
    elapsed::Float64
end

# Info-bit errors via the reverse encoder (msg_hat vs true message) and code-bit errors on the
# n channel bits.
function count_errors!(msg_hat::Vector{Bool}, Ginv::Matrix{Bool}, data::Vector{Bool}, tx_cw::Vector{Bool}, rx_cw::Vector{Bool})
    extract_message!(msg_hat, Ginv, rx_cw)
    n_err_data = 0
    @inbounds @simd for l in eachindex(data)
        n_err_data += msg_hat[l] ⊻ data[l]
    end
    n_err_total = 0
    @inbounds @simd for i in eachindex(tx_cw)
        n_err_total += tx_cw[i] ⊻ rx_cw[i]
    end
    return n_err_total, n_err_data
end

# MLD correlation metric, accumulated in Float64. Auxiliary VNs have zero LLR and add nothing.
function corr_metric(llr::Vector{LLRFloat}, cw::Vector{Bool})
    s = 0.0
    @inbounds @simd for i in eachindex(cw)
        s += Float64(llr[i]) * (1 - 2*cw[i])
    end
    return s
end
