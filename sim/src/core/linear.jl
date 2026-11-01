# The decoding matrix is the file's `H`: the optimized PCM (M × Ñ, with Ñ = n + n_avn) in a produced
# PCM file, the plain PCM (m × n, n_avn = 0) in a code file. Columns 1..n are the channel VNs; the
# trailing `n_avn` columns (listed in `puncturing_set`) are auxiliary VNs (redundant GF(2) parities
# of the channel bits) that carry no channel observation — zero channel LLR at decode time.
# Produced files also keep the un-optimized PCM as `H_base`, which `use_base_pcm=true` selects.
struct LinearCode
    name::String
    n::Int                              # channel/code length (columns of G, transmitted bits)
    k::Int
    n_avn::Int                          # auxiliary VN columns in H (0 for a base PCM)
    m::Int                              # rows (check nodes) of the decoding matrix H
    rate::Float64
    dmin::Union{Int,Nothing}
    G::Matrix{Bool}                     # k×n generator
    H::Matrix{Bool}                     # decoding matrix: M × (n + n_avn)
    Ginv::Matrix{Bool}                  # n×k reverse-encoder: (Ginv' * c) mod 2 recovers the message
end

# Coerce a MAT variable (BitMatrix / Bool / numeric 0-1, dense or sparse) to a dense Bool matrix.
to_bool_dense(A) = A isa AbstractMatrix{Bool} ? Matrix{Bool}(A) : Matrix{Bool}(A .!= 0)

# Load a code and pick the decoding matrix. With the optimized-PCM files under results/, `H` is the
# optimized PCM; `use_base_pcm=true` falls back to `H_base` with no auxiliary VNs — the baseline
# comparison curve. Plain code files have no `H_base`, so both modes decode their `H`.
function load_code(fname::AbstractString; use_base_pcm::Bool=false)
    endswith(fname, ".mat") || (fname *= ".mat")
    data = matread(fname)
    haskey(data, "G") || error("Matrix G not found in $fname")
    (haskey(data, "n") && haskey(data, "k")) || error("Parameters n and k not found in $fname")

    G = convert(Matrix{Bool}, data["G"])
    n, k = Int(data["n"]), Int(data["k"])
    (k, n) == size(G) || error("size(G)=$(size(G)) inconsistent with (k,n)=($k,$n)")

    key = (use_base_pcm && haskey(data, "H_base")) ? "H_base" : "H"
    haskey(data, key) || error("Matrix $key not found in $fname")
    H = to_bool_dense(data[key])
    n_avn = size(H, 2) - n
    n_avn >= 0 || error("$key has fewer than n=$n columns in $fname")
    if key == "H" && haskey(data, "puncturing_set")
        ps = Int.(vec(collect(data["puncturing_set"])))   # a 1-element set round-trips as a scalar
        ps == collect(n+1:n+n_avn) ||
            error("puncturing_set $(ps) is not the trailing $(n_avn) columns of H in $fname")
    end

    m = size(H, 1)
    dmin = haskey(data, "dmin") ? Int(data["dmin"]) : nothing
    name = haskey(data, "name") ? String(data["name"]) : "Linear"

    Ginv = build_ginv(G, k, n, get(data, "Ginv", nothing))
    return LinearCode(name, n, k, n_avn, m, k / n, dmin, G, H, Ginv)
end

function build_ginv(G::Matrix{Bool}, k::Int, n::Int, stored)
    Ginv = zeros(Bool, n, k)
    Ik = Matrix{Bool}(I, k, k)
    if stored !== nothing
        Ginv = convert(Matrix{Bool}, stored)
        size(Ginv) == (n, k) || error("stored Ginv has size $(size(Ginv)), expected ($n,$k)")
    elseif G[:, 1:k] == I
        Ginv[1:k, :] = Ik
    elseif G[:, n-k+1:n] == I
        Ginv[n-k+1:n, :] = Ik
    else
        error("can't find the identity matrix at the beginning or end of G")
    end
    (Int.(G) * Int.(Ginv)) .% 2 == I || error("can't find or build a valid reverse-encoding matrix Ginv")
    return Ginv
end

function Base.show(io::IO, c::LinearCode)
    dm = c.dmin === nothing ? "" : ",$(c.dmin)"
    aux = c.n_avn == 0 ? "" : " +$(c.n_avn) AVN"
    print(io, "$(c.name) ($(c.n),$(c.k)$dm)$aux")
end

# Random message into `data` (k bits), codeword cw = (Gᵀ·data) mod 2 into `cw` (n bits).
function encode!(rng::AbstractRNG, data::Vector{Bool}, cw::Vector{Bool}, G::Matrix{Bool})
    k, n = size(G)
    rand!(rng, data)
    @inbounds for j in 1:n
        acc = false
        @simd for i in 1:k
            acc ⊻= G[i, j] & data[i]
        end
        cw[j] = acc
    end
    return data
end

# Recover the message from a codeword via the reverse encoder: msg_hat = (Ginv'·cw) mod 2.
# Reads only the first n entries of `cw`, so an Ñ-length decoded vector (channel bits + AVNs)
# is handled correctly — the auxiliary tail is ignored.
function extract_message!(msg_hat::Vector{Bool}, Ginv::Matrix{Bool}, cw::Vector{Bool})
    n, k = size(Ginv)
    @inbounds for l in 1:k
        acc = false
        @simd for j in 1:n
            acc ⊻= Ginv[j, l] & cw[j]
        end
        msg_hat[l] = acc
    end
    return msg_hat
end
