# Build the three CRC24C polar mother codes used in the reproduction.
# Run from the repository root: julia --project src/build_polar_codes.jl [payload=20|40|72]
# With no argument, write all three files in codes/. No rate matching or CRC mask.
# Minimum distances are not certified, so the dmin field is omitted.
# Tables: 3GPP TS 38.212 V17.7.0, Tables 5.3.1.1-1 and 5.3.1.2-1.
# https://www.etsi.org/deliver/etsi_ts/138200_138299/138212/17.07.00_60/ts_138212v170700p.pdf
# Extracted from Sionna commit 2f0a3820269b53fd022780a6f2f5f720f7ceb9b4:
# src/sionna/phy/fec/polar/codes/polar_5G.csv and polar/encoding.py.
# Reliability entries >=128 are omitted; all table indices are zero-based.
using MAT, LinearAlgebra, Random
const ROOT = normpath(joinpath(@__DIR__, ".."))
const N, L = 128, 24
const POLYS = Dict("C" => [24,23,21,20,17,15,13,12,8,4,2,1,0])
const Q = [
    0, 1, 2, 4, 8, 16, 32, 3, 5, 64, 9, 6, 17, 10, 18, 12,
    33, 65, 20, 34, 24, 36, 7, 66, 11, 40, 68, 19, 13, 48, 14, 72,
    21, 35, 26, 80, 37, 25, 22, 38, 96, 67, 41, 28, 69, 42, 49, 74,
    70, 44, 81, 50, 73, 15, 52, 23, 76, 82, 56, 27, 97, 39, 84, 29,
    43, 98, 88, 30, 71, 45, 100, 51, 46, 75, 104, 53, 77, 54, 83, 57,
    112, 78, 85, 58, 99, 86, 60, 89, 101, 31, 90, 102, 105, 92, 47, 106,
    55, 113, 79, 108, 59, 114, 87, 116, 61, 91, 120, 62, 103, 93, 107, 94,
    109, 115, 110, 117, 118, 121, 122, 63, 124, 95, 111, 119, 123, 125, 126, 127,
]
const IL = [
    0, 2, 4, 7, 9, 14, 19, 20, 24, 25, 26, 28, 31, 34, 42, 45,
    49, 50, 51, 53, 54, 56, 58, 59, 61, 62, 65, 66, 67, 69, 70, 71,
    72, 76, 77, 81, 82, 83, 87, 88, 89, 91, 93, 95, 98, 101, 104, 106,
    108, 110, 111, 113, 115, 118, 119, 120, 122, 123, 126, 127, 129, 132, 134, 138,
    139, 140, 1, 3, 5, 8, 10, 15, 21, 27, 29, 32, 35, 43, 46, 52,
    55, 57, 60, 63, 68, 73, 78, 84, 90, 92, 94, 96, 99, 102, 105, 107,
    109, 112, 114, 116, 121, 124, 128, 130, 133, 135, 141, 6, 11, 16, 22, 30,
    33, 36, 44, 47, 64, 74, 79, 85, 97, 100, 103, 117, 125, 131, 136, 142,
    12, 17, 23, 37, 48, 75, 80, 86, 137, 143, 13, 18, 38, 144, 39, 145,
    40, 146, 41, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
    160, 161, 162, 163,
]
@assert sort(Q) == 0:127 && sort(IL) == 0:163
positions(A) = (K = A + L; info = sort(Q[end-K+1:end]) .+ 1; (info, setdiff(1:N, info)))

function crc(msg, poly)
    A = length(msg)
    work = vcat(copy(msg), falses(L))
    for i in 1:A
        if work[i]
            for exponent in POLYS[poly]
                work[i + L - exponent] ⊻= true
            end
        end
    end
    work[A+1:end]
end

function polar!(x)
    stride = 1
    while stride < N
        for base in 1:2stride:N, j in 0:stride-1
            x[base+j] ⊻= x[base+j+stride]
        end
        stride *= 2
    end
    x
end

function build(; payload=20)
    poly = "C"
    payload in (20, 40, 72) || throw(ArgumentError("payload must be 20, 40, or 72"))
    A = payload; K = A + L
    INFO, FROZEN = positions(A)
    perm = [p-(164-K)+1 for p in IL if p >= 164-K]
    @assert sort(perm) == 1:K
    # Each row encodes one message basis vector, using polynomial long division.
    C = falses(A,L)
    for i in 1:A
        b = falses(A); b[i] = true
        C[i,:] = crc(b,poly)
    end
    F = falses(N,N)
    for i in 1:N
        b = falses(N); b[i] = true; F[i,:] = polar!(b)
    end
    G = falses(A,N)
    for i in 1:A
        b = falses(A); b[i] = true
        u = falses(N); u[INFO] = vcat(b, C[i,:])[perm]
        G[i,:] = polar!(u)
    end
    # x F = u since F is self inverse. CRC checks are p_j + sum_i C_ij a_i.
    P = falses(K,N)
    for j in 1:K
        P[perm[j],:] = F[:,INFO[j]]
    end
    H = falses(N-A,N)
    H[1:length(FROZEN),:] = transpose(F[:,FROZEN])
    for j in 1:L
        H[length(FROZEN)+j,:] = P[A+j,:]
        for i in 1:A
            C[i,j] && (H[length(FROZEN)+j,:] .⊻= P[i,:])
        end
    end
    Ginv = permutedims(P[1:A,:])
    G,H,Ginv,perm
end

gfprod(X,Y) = isodd.(Int.(X)*Int.(Y))
function rank2(M)
    B=copy(M); row=1
    for col in axes(B,2)
        p=findfirst(i->B[i,col], row:size(B,1))
        isnothing(p) && continue
        pivot=row+p-1
        B[row,:],B[pivot,:]=copy(B[pivot,:]),copy(B[row,:])
        for i in row+1:size(B,1)
            B[i,col] && (B[i,:] .⊻= B[row,:])
        end
        row+=1; row>size(B,1) && break
    end
    row-1
end

function check(G,H,Ginv,perm,poly)
    A = size(G, 1)
    INFO, _ = positions(A)
    @assert rank2(G)==A && rank2(H)==N-A
    @assert !any(gfprod(G,transpose(H)))
    @assert gfprod(G,Ginv)==Matrix{Bool}(I,A,A)
    # Independent shift-register CRC and recursive butterfly; compare random and basis inputs.
    mask = sum(1 << e for e in POLYS[poly] if e != 24)
    rng=MersenneTwister(38212)
    for t in 1:120
        msg = t<=A ? [i==t for i in 1:A] : rand(rng,Bool,A)
        reg=0
        for bit in msg
            top = ((reg >> 23) & 1) ⊻ Int(bit)
            reg = (reg << 1) & 0xffffff
            top==1 && (reg ⊻= mask)
        end
        parity = [Bool((reg >> i)&1) for i in 23:-1:0]
        u=falses(N); u[INFO]=vcat(msg,parity)[perm]
        # Recursive F tensor power, deliberately distinct from matrix multiplication.
        function rec(v)
            length(v)==1 && return v
            h=length(v)÷2; a=rec(v[1:h]); b=rec(v[h+1:end])
            vcat(a .⊻ b,b)
        end
        cw=rec(u)
        @assert cw == vec(gfprod(reshape(msg,1,:),G))
        @assert !any(gfprod(H,reshape(cw,:,1)))
    end
end

function main()
    payloads = if isempty(ARGS)
        (20, 40, 72)
    elseif length(ARGS) == 1 && ARGS[1] in ("payload=20", "payload=40", "payload=72")
        (parse(Int, split(ARGS[1], "=")[2]),)
    else
        error("usage: julia --project src/build_polar_codes.jl [payload=20|40|72]")
    end
    for A in payloads
        G, H, Ginv, perm = build(; payload=A)
        check(G, H, Ginv, perm, "C")
        info, frozen = positions(A)
        path = joinpath(ROOT, "codes", "polar.128.$A.crc24c.mat")
        matwrite(path, Dict(
            "name" => "5G NR polar CRC24C", "n" => N, "k" => A,
            "G" => Matrix(G), "H" => Matrix(H), "Ginv" => Matrix(Ginv),
            "crc" => "CRC24C", "input_interleaving" => 1,
            "info_positions_0based" => info .- 1, "frozen_positions_0based" => frozen .- 1,
            "input_permutation_0based" => perm .- 1, "crc_mask" => 0,
            "standard" => "3GPP TS 38.212 V17.7.0 sections 5.1, 5.3.1.1, 5.3.1.2",
            "scope" => "N=128 mother code before rate matching; no CRC initialization or RNTI mask"))
        println("$path: G $(size(G)), H $(size(H)), rank(H) $(rank2(H)), edges $(sum(H))")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
