# Growable binary matrix over GF(2), stored row-major as bit-packed rows (BitVector).
# Algorithm 1 grows the matrix by rows (ACNs) and columns (AVNs) and merges/deletes them,
# so rows are independently resizable and column ops touch every row.
mutable struct GF2Matrix
    rows::Vector{BitVector}
    ncols::Int
end

GF2Matrix(m::Integer, n::Integer) = GF2Matrix([falses(n) for _ in 1:m], Int(n))

nrows(A::GF2Matrix) = length(A.rows)
ncols(A::GF2Matrix) = A.ncols

function from_bool(H::AbstractMatrix{Bool})
    m, n = size(H)
    GF2Matrix([BitVector(@view H[i, :]) for i in 1:m], n)
end

function to_bool(A::GF2Matrix)
    H = falses(nrows(A), A.ncols)
    @inbounds for i in 1:nrows(A), j in 1:A.ncols
        A.rows[i][j] && (H[i, j] = true)
    end
    H
end

Base.copy(A::GF2Matrix) = GF2Matrix([copy(r) for r in A.rows], A.ncols)

# Count of set bits in a & b — the row-pair overlap / column-support intersection size.
@inline function and_count(a::BitVector, b::BitVector)
    s = 0
    ca, cb = a.chunks, b.chunks
    @inbounds @simd for k in eachindex(ca)
        s += count_ones(ca[k] & cb[k])
    end
    s
end

# Count of set bits in a & b & c — size of a three-way column-support intersection.
@inline function and3_count(a::BitVector, b::BitVector, c::BitVector)
    s = 0
    ca, cb, cc = a.chunks, b.chunks, c.chunks
    @inbounds @simd for k in eachindex(ca)
        s += count_ones(ca[k] & cb[k] & cc[k])
    end
    s
end

# Count of set bits in a ⊻ b — the weight of the XOR (number of differing bits).
@inline function xor_count(a::BitVector, b::BitVector)
    s = 0
    ca, cb = a.chunks, b.chunks
    @inbounds @simd for k in eachindex(ca)
        s += count_ones(ca[k] ⊻ cb[k])
    end
    s
end

# a ← a ⊻ b in place (unused high bits stay zero since both operands have zeros there).
@inline function xor_bits!(a::BitVector, b::BitVector)
    ca, cb = a.chunks, b.chunks
    @inbounds @simd for k in eachindex(ca)
        ca[k] ⊻= cb[k]
    end
    a
end

# a ⊆ b : every set bit of a is set in b.
@inline function subset_of(a::BitVector, b::BitVector)
    ca, cb = a.chunks, b.chunks
    @inbounds for k in eachindex(ca)
        (ca[k] & ~cb[k]) != 0 && return false
    end
    true
end

overlap(A::GF2Matrix, i::Int, j::Int) = and_count(A.rows[i], A.rows[j])
xorrow!(A::GF2Matrix, dst::Int, src::Int) = xor_bits!(A.rows[dst], A.rows[src])
rowweight(A::GF2Matrix, i::Int) = count(A.rows[i])
num_edges(A::GF2Matrix) = sum(count, A.rows; init = 0)

function colweights(A::GF2Matrix)
    w = zeros(Int, A.ncols)
    @inbounds for r in A.rows, j in 1:A.ncols
        r[j] && (w[j] += 1)
    end
    w
end

function addrow!(A::GF2Matrix, r::BitVector)
    length(r) == A.ncols || error("row length $(length(r)) != ncols $(A.ncols)")
    push!(A.rows, r)
    A
end

# Append one all-zero column; return the new column index.
function addcol!(A::GF2Matrix)
    @inbounds for r in A.rows
        push!(r, false)
    end
    A.ncols += 1
    A.ncols
end

function deletecol!(A::GF2Matrix, j::Int)
    @inbounds for r in A.rows
        deleteat!(r, j)
    end
    A.ncols -= 1
    A
end

deleterow!(A::GF2Matrix, i::Int) = (deleteat!(A.rows, i); A)

binom2(x::Integer) = x * (x - 1) ÷ 2

# Number of 4-cycles: every row pair sharing t≥2 columns contributes C(t,2).
function count_4cycles(A::GF2Matrix)
    q = 0
    m = nrows(A)
    @inbounds for i in 1:m-1, j in i+1:m
        q += binom2(overlap(A, i, j))
    end
    q
end

# GF(2) rank via bit-packed Gaussian elimination on a copy.
function rank_gf2(rows_in::Vector{BitVector}, n::Int)
    isempty(rows_in) && return 0
    rows = [copy(r) for r in rows_in]
    m = length(rows)
    rank = 0
    @inbounds for c in 1:n
        p = 0
        for i in rank+1:m
            if rows[i][c]
                p = i
                break
            end
        end
        p == 0 && continue
        rank += 1
        rows[rank], rows[p] = rows[p], rows[rank]
        for i in 1:m
            i != rank && rows[i][c] && xor_bits!(rows[i], rows[rank])
        end
        rank == m && break
    end
    rank
end
rank_gf2(A::GF2Matrix) = rank_gf2(A.rows, A.ncols)
