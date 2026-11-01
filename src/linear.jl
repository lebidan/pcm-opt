# A code read from a .mat file: parity-check matrix H (m×n), generator G (k×n), and optionally its
# minimum distance and a display name. A sparse H is stored dense: the codes here are short.
struct LinearCode
    name::String
    n::Int
    k::Int
    m::Int
    dmin::Union{Int,Nothing}
    G::Matrix{Bool}
    H::Matrix{Bool}
end

function load_code(fname::AbstractString)
    endswith(fname, ".mat") || (fname *= ".mat")
    data = matread(fname)
    haskey(data, "G") || error("Matrix G not found in $fname")
    haskey(data, "H") || error("Matrix H not found in $fname")
    (haskey(data, "n") && haskey(data, "k")) || error("Parameters n and k not found in $fname")

    G = Matrix{Bool}(data["G"])
    H = Matrix{Bool}(data["H"])
    n, k = Int(data["n"]), Int(data["k"])
    (k, n) == size(G) || error("size(G)=$(size(G)) inconsistent with (k,n)=($k,$n)")
    dmin = haskey(data, "dmin") ? Int(data["dmin"]) : nothing
    name = haskey(data, "name") ? String(data["name"]) : "Linear"
    return LinearCode(name, n, k, size(H, 1), dmin, G, H)
end

function Base.show(io::IO, c::LinearCode)
    dm = c.dmin === nothing ? "" : ",$(c.dmin)"
    print(io, "$(c.name) ($(c.n),$(c.k)$dm)")
end
