using Test, Random, MAT
using PCMOpt

const CODES = joinpath(@__DIR__, "..", "codes")
loadA(name) = (c = load_code(joinpath(CODES, name)); (c, from_bool(c.H)))

# RM(32,16) search space and Table I PCM-II, shared by several test sets.
const RM32 = load_code(joinpath(CODES, "rm.32.16"))
const RM32_HS = gen_search_space(from_bool(RM32.H); T_max = 8000, seed = 1)
const RM32_PCM2 = extract_pcm(RM32_HS, RM32.n, RM32.k; M_a = 48, seed = 6)

@testset "PCMOpt" begin

    @testset "code loading" begin
        for (name, m, n) in (("rm.8.4", 4, 8), ("rm.32.16", 16, 32),
                             ("bch.63.45", 18, 63), ("ebch.64.24", 40, 64))
            c = load_code(joinpath(CODES, name))
            @test (c.m, c.n) == (m, n)
        end
    end

    @testset "GF(2) + cycle utilities" begin
        # Hand-built single hexagon c1-v1-c2-v2-c3-v3-c1, 4th CN lifts v3 degree to 3.
        A = from_bool(Bool[1 0 1; 1 1 0; 0 1 1; 0 0 1])
        @test count_4cycles(A) == 0
        @test count_6cycles(A) == 1
        @test girth(A) == 6
        @test eta_min(A, 2) == 1          # ACE of the lone 6-cycle: (2-2)+(2-2)+(3-2)

        # RM(32,16) standard PCM reproduces Table I's Standard row exactly.
        _, H = loadA("rm.32.16")
        @test (nrows(H), ncols(H)) == (16, 32)
        @test num_edges(H) == 192
        @test rank_gf2(H) == 16
        @test count_6cycles(H) == 91_580
        @test girth(H) == 4               # standard PCM has 4-cycles

        # GF(2) primitive spot checks.
    end

    @testset "Algorithm 1 (4-cycle removal)" begin
        # RM(32,16) PCM-I vs Table I: size 32×48, 122 edges, 16 AVNs, Q4=0, girth 6, η̄_min=2.
        # (6-cycles: this greedy path yields 39; the paper reports 38 — a one-cycle tie-break
        # difference at identical size/edge count, within greedy 4-cycle-removal sensitivity.)
        c, H = loadA("rm.32.16")
        A = remove_4cycles(H)
        @test (nrows(A), ncols(A)) == (32, 48)
        @test num_edges(A) == 122
        @test ncols(A) - c.n == 16              # AVNs
        @test count_4cycles(A) == 0
        @test girth(A) == 6
        @test eta_min(A, c.n) == 2
        @test count_6cycles(A) == 39

        # 4-cycle removal produces a valid girth-6, 4-cycle-free PCM on every target code.
        for name in ("rm.8.4", "bch.63.45", "ebch.64.24")
            c, H = loadA(name)
            A = remove_4cycles(H)
            @test count_4cycles(A) == 0
            @test girth(A) == 6
            @test ncols(A) >= c.n               # only AVN columns appended
        end

        # Randomized tie-break also yields valid girth-6 PCMs; seed 1042 gives a PCM-I that
        # dominates the deterministic one (fewer edges & 6-cycles, higher η̄_min).
        c, H = loadA("rm.32.16")
        for s in (1, 7, 1042)
            A = remove_4cycles_rand(H, Xoshiro(s))
            @test count_4cycles(A) == 0 && girth(A) == 6
        end
        m = pcm_metrics(remove_4cycles_rand(H, Xoshiro(1042)), c.n)
        @test (m.rows, m.cols, m.edges, m.sixcyc, m.eta_min) == (32, 48, 121, 24, 3)
    end

    @testset "Algorithm 1 cleanups repeat until none applies" begin
        # N = 3 real VNs, AVNs 4–6. Merging the degree-2 CN {4,5} turns {4,5,6} into the degree-1
        # CN {6}, whose AVN must then go too, leaving {1,2} and {1,3,4}. A single pass of each
        # cleanup (the old behaviour) left that degree-1 CN in place.
        A = from_bool(Bool[0 0 0 1 1 0; 0 0 0 1 1 1; 1 1 0 0 0 1; 1 0 1 1 0 0])
        PCMOpt.cleanup!(A, 3)
        @test PCMOpt.to_bool(A) == Bool[1 1 0 0; 1 0 1 1]
    end

    @testset "Algorithm 2 (search space)" begin
        is_dual(G, r, k, n) = all(a -> !reduce(⊻, (G[a, j] & r[j] for j in 1:n)), 1:k)
        # Complete low-weight dual-codeword sets (deterministic given enough trials); sizes and
        # uniform weights match the paper: RM(8,4)→14×wt4, RM(32,16)→620×wt8 (its 620×32 search
        # space), BCH(63,45)→189×wt16 (the paper's 189×63 H_s), eBCH(64,24)→2520×wt8.
        for (name, nrow, wt) in (("rm.8.4", 14, 4), ("rm.32.16", 620, 8),
                                 ("bch.63.45", 189, 16), ("ebch.64.24", 2520, 8))
            c, _ = loadA(name)
            Hs = gen_search_space(from_bool(c.H); T_max = 8000, seed = 1)
            @test nrows(Hs) == nrow
            @test rank_gf2(Hs) == c.n - c.k
            @test all(i -> rowweight(Hs, i) == wt, 1:nrows(Hs))
            @test all(i -> is_dual(c.G, Hs.rows[i], c.k, c.n), 1:nrows(Hs))
        end
    end

    @testset "Algorithm 3 (PRG)" begin
        # RM(8,4): the full pipeline reproduces eq.(16) — H_b^{(8)} is 10×10, 32 edges,
        # column weight 3.2, girth 6, 2 AVNs, η̄_min=4.
        c, _ = loadA("rm.8.4")
        Hs = gen_search_space(from_bool(c.H); T_max = 4000, seed = 1)
        r = extract_pcm(Hs, c.n, c.k; M_a = 14, seed = 1)
        @test (nrows(r.Hb), ncols(r.Hb)) == (10, 10)
        @test num_edges(r.Hb) == 32
        @test ncols(r.Hb) - c.n == 2
        @test count_4cycles(r.Hb) == 0 && girth(r.Hb) == 6
        @test r.eta_min == 4

        # RM(32,16): reproduces Table I's PCM-II exactly — peak at i=42: 48×38, 240 edges,
        # 2336 6-cycles, 6 AVNs, η̄_min=14 (deterministic given the fixed search-space/PRG seeds).
        c, r = RM32, RM32_PCM2
        @test r.i == 42
        @test (nrows(r.Hb), ncols(r.Hb)) == (48, 38)
        @test num_edges(r.Hb) == 240
        @test ncols(r.Hb) - c.n == 6
        @test count_4cycles(r.Hb) == 0 && girth(r.Hb) == 6
        @test count_6cycles(r.Hb) == 2336
        @test r.eta_min == 14
    end

    @testset "BCH(63,45) PCM-II + pcm_metrics" begin
        # BCH(63,45): reproduces the paper text — PCM-II is 63×72, 567 edges, uniform row weight 9.
        c, _ = loadA("bch.63.45")
        Hs = gen_search_space(from_bool(c.H); T_max = 8000, seed = 1)
        r = extract_pcm(Hs, c.n, c.k; M_a = 70, seed = 1)
        m = pcm_metrics(r.Hb, c.n)
        @test (m.rows, m.cols) == (63, 72)
        @test m.edges == 567
        @test m.rowwt == 9.0            # uniform row weight
        @test m.avn == 9
        @test m.q4 == 0 && m.girth == 6

        # pcm_metrics agrees with the standalone functions.
        @test m.edges == num_edges(r.Hb)
        @test m.sixcyc == count_6cycles(r.Hb)
        @test m.eta_min == eta_min(r.Hb, c.n)
    end

    @testset "Algorithm 4 (rapid path)" begin
        c, Hs = RM32, RM32_HS
        # Rapid mode yields a valid girth-6, 4-cycle-free PCM with good connectivity (peak η̄_min ≥ 12).
        r = extract_pcm(Hs, c.n, c.k; M_a = 55, seed = 3, rapid = true)
        @test count_4cycles(r.Hb) == 0
        @test girth(r.Hb) == 6
        @test ncols(r.Hb) >= c.n
        @test r.eta_min >= 12
        # rapid_bypass, when it succeeds, adds exactly one row and no AVN column.
        Hb0 = remove_4cycles(GF2Matrix([copy(Hs.rows[1])], c.n))
        by = rapid_bypass(Hb0, Hs.rows[2], c.n)
        if by !== nothing
            @test ncols(by) == ncols(Hb0)
            @test nrows(by) == nrows(Hb0) + 1
            @test count_4cycles(by) == 0
        end
    end

    @testset "trapping-set distribution 𝒯(a,b)" begin
        # Single hexagon (one 6-cycle): the triple {v1,v2,v3} is a 𝒯(3,1) dominant TS (CN4 is the lone
        # odd CN); the 2-VN subsets are trees (no cycle) and are excluded.
        hex = from_bool(Bool[1 0 1; 1 1 0; 0 1 1; 0 0 1])
        @test trapping_distribution(hex; a_max = 4) == Dict((3, 1) => 1)
        @test min_trapping_set(trapping_distribution(hex)) == (3, 1)

        # A 4-cycle (v1,v2 share CN1,CN2) with v1 degree 3: a 𝒯(2,1) — v2 is degree-2 (exempt), v1 has
        # its two shared CNs satisfied. Exercises the 𝒯(2,·)/degree-2-exemption path.
        c4 = from_bool(Bool[1 1; 1 1; 1 0])
        @test count_4cycles(c4) == 1
        @test trapping_distribution(c4; a_max = 4) == Dict((2, 1) => 1)

        # RM(32,16) PCM-II reproduces Table I's dominant trapping sets exactly (same matrix as the Algorithm 3 test).
        d = trapping_distribution(RM32_PCM2.Hb; a_max = 4)
        @test d[(3, 12)] == 896
        @test d[(4, 12)] == 1088
        @test min_trapping_set(d) == (3, 12)               # paper's headline: smallest TS is 𝒯(3,12)
        @test !any(k -> k[1] == 2, keys(d))                # girth 6 ⇒ no 2-VN dominant TS
        @test !any(k -> k[1] == 3 && k[2] < 12, keys(d))   # "nothing lower or in-between" below 𝒯(3,12)

        # PCM-I (girth 6): the deterministic matrix differs from the paper's, but its signature holds —
        # no 2-VN sets and the smallest trapping set is a 6-cycle 𝒯(3,·).
        A1 = remove_4cycles(from_bool(RM32.H))
        d1 = trapping_distribution(A1; a_max = 4)
        @test !any(k -> k[1] == 2, keys(d1))
        @test min_trapping_set(d1)[1] == 3
    end

    @testset "save_optpcm output format" begin
        c, H = loadA("rm.8.4")
        A = remove_4cycles(H)
        tmp = tempname() * ".mat"
        save_optpcm(tmp, A, joinpath(CODES, "rm.8.4"))
        d = matread(tmp)
        # source variables preserved with original types
        @test d["G"] isa Matrix{Bool} && d["n"] isa UInt32 && d["name"] == "RM"
        @test Int(d["n"]) == c.n && Int(d["k"]) == c.k
        # optimized PCM stored as H (logical dense), source PCM kept as H_base
        @test d["H"] isa Matrix{Bool}
        @test size(d["H"]) == (nrows(A), ncols(A)) && count(d["H"]) == num_edges(A)
        @test size(d["H_base"], 2) == c.n
        # puncturing set: 1-based auxiliary columns; eta_min as Int32
        @test eltype(d["puncturing_set"]) == Int32       # a 1-element set round-trips as a scalar
        @test vec(collect(d["puncturing_set"])) == Int32[c.n+1:ncols(A);]
        @test d["eta_min"] isa Int32
        rm(tmp; force = true)
    end

    @testset "guards" begin
        c, H = loadA("rm.32.16")
        # eta_min rejects matrices that still contain 4-cycles (the standard PCM).
        @test count_4cycles(H) > 0
        @test_throws ArgumentError eta_min(H, c.n)
        # extract_pcm errors rather than returning a bogus 1-row result when no full-rank PCM is
        # reached (here M_a < N−K).
        @test_throws ErrorException extract_pcm(RM32_HS, c.n, c.k; M_a = 5, seed = 1)
    end

end
