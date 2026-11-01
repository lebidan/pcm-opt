# BPSK over BI-AWGN: sample y = (1-2·bit) + σ·w, w ~ N(0,1). With `scale`, feed the true LLR
# 2y/σ²; without it, feed the raw sample y (unscaled) — magnitude matters for SPA, not min-sum.
# Noise drawn in Float64, stored at LLR precision. Only the first n entries of `rx_llr` are written.
function transmit!(rng::AbstractRNG, tx_bits::Vector{Bool}, rx_llr::Vector{LLRFloat}, σ::Float64, scale::Bool=true)
    g = scale ? 2. / (σ*σ) : 1.0
    @inbounds for i in eachindex(tx_bits)
        rx_llr[i] = LLRFloat(g * (1 - 2*tx_bits[i] + σ*randn(rng)))
    end
end
