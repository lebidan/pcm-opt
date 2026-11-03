using BPsim
using BPsim: encode!, transmit!, count_errors!, run_snr_point!, parse_algorithm, corr_metric,
    parse_decoder
using Test, Random, LinearAlgebra

# Resolve files against the repo root (one level above the package dir) so tests work
# regardless of the working directory.
repofile(rel) = normpath(joinpath(pkgdir(BPsim), "..", rel))

# Optimized-PCM artifacts and their expected auxiliary-VN counts (as in the README).
const OPT_FILES = [
    "results/rm.8.4.pcm2.mat"    => 2,
    "results/rm.32.16.pcm2.mat"  => 6,
    "results/rm.32.16.pcm1.mat"  => 12,
    "results/bch.63.45.pcm2.mat" => 9,
    "results/ebch.64.24.pcm2.mat" => 17,
    "results/ebch.64.24.pcm2.eta52.mat" => 21,
    "results/rm.128.64.pcm2.mat" => 33,
]

bp(code, cn, sched=Flooding, max_iter=50) =
    BPDecoder(TannerGraph(code.H), BPConfig(sched, max_iter, cn))

# GF(2) consistency of A·a = b (does a solution `a` exist?) — Gaussian elimination over GF(2).
function gf2_consistent(A::AbstractMatrix{Bool}, b::AbstractVector{Bool})
    M, N = size(A)
    aug = hcat(Matrix{Bool}(A), reshape(collect(Bool, b), M, 1))
    row = 1
    for col in 1:N
        row > M && break
        piv = findnext(@view(aug[:, col]), row)
        piv === nothing && continue
        aug[[row, piv], :] = aug[[piv, row], :]
        for r in 1:M
            (r != row && aug[r, col]) && (aug[r, :] .⊻= aug[row, :])
        end
        row += 1
    end
    return !any(r -> !any(@view aug[r, 1:N]) && aug[r, N+1], 1:M)
end

@testset "BPsim" begin

    @testset "load_code: optimized-H dims + base fallback ($f)" for (f, want_avn) in OPT_FILES
        code = load_code(repofile(f))
        @test code.n_avn == want_avn
        @test size(code.H, 2) == code.n + code.n_avn
        @test size(code.G) == (code.k, code.n)
        @test size(code.Ginv) == (code.n, code.k)
        @test (Int.(code.G) * Int.(code.Ginv)) .% 2 == Matrix{Int}(I, code.k, code.k)
        g = TannerGraph(code.H)
        @test g.n_vn == code.n + code.n_avn
        # base-PCM fallback: no auxiliary VNs, decoding matrix is the m×n H_base
        base = load_code(repofile(f); use_base_pcm=true)
        @test base.n_avn == 0
        @test size(base.H, 2) == base.n
        @test base.n == code.n && base.k == code.k
    end

    @testset "plain code file has no H_base → n_avn=0" begin
        code = load_code(repofile("codes/rm.32.16.mat"))
        @test code.n_avn == 0
        @test size(code.H, 2) == code.n
    end

    @testset "algorithm names" begin
        @test parse_algorithm("MS") == MinSum
        @test parse_algorithm("FastSPA") == FastSPA
        @test parse_algorithm("spa") == SPAExact
        @test_throws ErrorException parse_algorithm("XYZ")
    end

    @testset "LLR scaling toggle (off = raw samples)" begin
        n, σ = 64, 0.7
        tx = rand(Xoshiro(9), Bool, n)
        scaled, raw, default = zeros(LLRFloat, n), zeros(LLRFloat, n), zeros(LLRFloat, n)
        transmit!(Xoshiro(42), tx, scaled, σ, true)
        transmit!(Xoshiro(42), tx, raw, σ, false)
        transmit!(Xoshiro(42), tx, default, σ)
        @test scaled == default
        @test scaled ./ LLRFloat(2/σ^2) ≈ raw rtol=1e-3
    end

    @testset "auxiliary-VN LLR tail stays zero after transmit!" begin
        code = load_code(repofile("results/rm.32.16.pcm2.mat"))
        Ñ = code.n + code.n_avn
        tx = zeros(Bool, code.n)
        llr = zeros(LLRFloat, Ñ)
        transmit!(Xoshiro(3), tx, llr, 0.8)
        @test all(llr[code.n+1:Ñ] .== 0)      # the auxiliary tail is never written
        @test any(llr[1:code.n] .!= 0)         # the channel positions are
    end

    @testset "error counting & MLD metric ignore the auxiliary tail" begin
        code = load_code(repofile("results/rm.32.16.pcm2.mat"))
        Ñ = code.n + code.n_avn
        rng = Xoshiro(21)
        data, msg_hat = zeros(Bool, code.k), zeros(Bool, code.k)
        tx = zeros(Bool, code.n)
        rx = zeros(Bool, Ñ)
        llr = zeros(LLRFloat, Ñ)
        encode!(rng, data, tx, code.G)
        transmit!(rng, tx, llr, 0.9)
        rand!(rng, rx)                                    # arbitrary decoded word (incl. aux tail)
        tot1, dat1 = count_errors!(msg_hat, code.Ginv, data, tx, rx)
        cm1 = corr_metric(llr, rx)
        rx[code.n+1:Ñ] .= .!rx[code.n+1:Ñ]                # flip only the auxiliary tail
        tot2, dat2 = count_errors!(msg_hat, code.Ginv, data, tx, rx)
        cm2 = corr_metric(llr, rx)
        @test (tot1, dat1) == (tot2, dat2)                # counts depend only on the first n bits
        @test cm1 == cm2                                  # AVN LLRs are 0 → no metric contribution
        @test tot1 == count(rx[1:code.n] .⊻ tx)
    end

    @testset "every base codeword extends into the augmented code ($f)" for (f, _) in OPT_FILES
        code = load_code(repofile(f))
        n, Ñ = code.n, code.n + code.n_avn
        Hc, Ha = code.H[:, 1:n], code.H[:, n+1:Ñ]        # channel block | auxiliary block
        rng = Xoshiro(99)
        data, c = zeros(Bool, code.k), zeros(Bool, n)
        for _ in 1:50
            encode!(rng, data, c, code.G)                # base codeword (H·c = 0)
            b = Bool.((Int.(Hc) * Int.(c)) .% 2)         # need Ha·a = b for some auxiliary a
            @test gf2_consistent(Ha, b)
        end
    end

    @testset "all algorithms × schedules recover a clean codeword through the optimized H" begin
        code = load_code(repofile("results/rm.32.16.pcm2.mat"))
        Ñ = code.n + code.n_avn
        rng = Xoshiro(7)
        data, tx, rx = zeros(Bool, code.k), zeros(Bool, code.n), zeros(Bool, Ñ)
        llr = zeros(LLRFloat, Ñ)
        encode!(rng, data, tx, code.G)
        transmit!(rng, tx, llr, 0.2)                     # very low noise → must converge
        g = TannerGraph(code.H)
        for cn in instances(CNRule), sched in instances(LdpcSchedule)
            dec = bp(code, cn, sched, 50)
            @test first(decode!(dec, llr, rx)) == Codeword
            @test rx[1:code.n] == tx                      # channel bits recovered
            @test is_codeword(rx, g)                      # full Ñ-vector satisfies the optimized H
        end
    end

    @testset "high-SNR FER → 0 on optimized H ($f)" for (f, _) in OPT_FILES
        code = load_code(repofile(f))
        cfg = load_config(nothing, Dict{String,Any}("code_file" => repofile(f),
            "ebn0_min" => 7.0, "ebn0_max" => 7.0, "min_cw" => 3000, "min_cw_err" => 0,
            "batch_size" => 500, "decoder" => Dict("algorithm" => "MS", "schedule" => "flooding", "max_iter" => 50)))
        graph = TannerGraph(code.H)
        σ = sqrt(10^(-7.0/10) / (2*code.rate))
        c, _ = run_snr_point!(code, graph, cfg, 1, σ, (a...) -> nothing)
        @test c.n_cw_err / c.n_cw < 0.01                  # deep in the waterfall at 7 dB
    end

    @testset "deterministic at fixed thread count (optimized H)" begin
        cfg = load_config(nothing, Dict{String,Any}(
            "code_file" => repofile("results/rm.32.16.pcm2.mat"),
            "ebn0_min" => 3.0, "ebn0_max" => 3.0, "min_cw" => 4000, "min_cw_err" => 40, "batch_size" => 500))
        code = load_code(cfg.code_file)
        graph = TannerGraph(code.H)
        σ = sqrt(10^(-3.0/10) / (2*code.rate))
        noop = (args...) -> nothing
        c1, _ = run_snr_point!(code, graph, cfg, 1, σ, noop)
        c2, _ = run_snr_point!(code, graph, cfg, 1, σ, noop)
        @test (c1.n_cw, c1.n_cw_err, c1.n_bit_err, c1.n_code_bit_err) ==
              (c2.n_cw, c2.n_cw_err, c2.n_bit_err, c2.n_code_bit_err)
    end

    @testset "all-zero mode matches random-codeword FER (optimized H)" begin
        base = Dict{String,Any}("code_file" => repofile("results/rm.32.16.pcm2.mat"),
            "ebn0_min" => 2.0, "ebn0_max" => 2.0, "min_cw" => 6000, "min_cw_err" => 200, "batch_size" => 500)
        code = load_code(base["code_file"])
        graph = TannerGraph(code.H)
        σ = sqrt(10^(-2.0/10) / (2*code.rate))
        noop = (args...) -> nothing
        c_rand, _ = run_snr_point!(code, graph, load_config(nothing, base), 1, σ, noop)
        c_zero, _ = run_snr_point!(code, graph, load_config(nothing, merge(base, Dict("all_zero" => true))), 1, σ, noop)
        fer_rand = c_rand.n_cw_err / c_rand.n_cw
        fer_zero = c_zero.n_cw_err / c_zero.n_cw
        @test 0.5 < fer_zero / fer_rand < 2.0
        @test c_zero.n_code_bit_err ≥ c_zero.n_cw_err
        @test c_zero.n_code_bit_err > 0
    end

    @testset "normalized min-sum" begin
        @test parse_algorithm("NMS") == NormMinSum
        @test parse_decoder(Dict("algorithm"=>"NMS")).cn_rule == NormMinSum
        @test parse_decoder(Dict("algorithm"=>"NMS")).norm_factor == 0.7f0   # default α
        @test parse_decoder(Dict("algorithm"=>"NMS", "norm_factor"=>0.8)).norm_factor == 0.8f0
        @test_throws ErrorException parse_decoder(Dict("norm_factor"=>0))
        @test_throws ErrorException parse_decoder(Dict("norm_factor"=>-0.5))
        @test decoder_name(BPConfig(Flooding, 10, NormMinSum)) == "NMS"

        # α = 1 is exactly plain min-sum (1f0·x is exact in Float32); α = 0.7 scales the messages.
        code = load_code(repofile("results/rm.32.16.pcm2.mat"))
        Ñ = code.n + code.n_avn
        rng = Xoshiro(2024)
        data, tx = zeros(Bool, code.k), zeros(Bool, code.n)
        rx_ms, rx_nms1 = zeros(Bool, Ñ), zeros(Bool, Ñ)
        llr = zeros(LLRFloat, Ñ)
        dec_ms   = BPDecoder(TannerGraph(code.H), BPConfig(Flooding, 50, MinSum))
        dec_nms1 = BPDecoder(TannerGraph(code.H), BPConfig(Flooding, 50, NormMinSum, 1f0))
        dec_nms07 = BPDecoder(TannerGraph(code.H), BPConfig(Flooding, 50, NormMinSum, 0.7f0))
        encode!(rng, data, tx, code.G)
        transmit!(rng, tx, llr, 0.9)
        decode!(dec_ms, llr, rx_ms)
        decode!(dec_nms1, llr, rx_nms1)
        @test rx_nms1 == rx_ms                       # α = 1 reproduces min-sum exactly ...
        @test dec_nms1.llr_out == dec_ms.llr_out     # ... down to the soft output
        # α = 0.7 must alter the soft output on some frame where BP actually iterates.
        differed = false
        for _ in 1:20
            encode!(rng, data, tx, code.G)
            transmit!(rng, tx, llr, 1.1)
            decode!(dec_ms, llr, rx_ms)
            decode!(dec_nms07, llr, rx_nms1)
            dec_nms07.llr_out != dec_ms.llr_out && (differed = true; break)
        end
        @test differed
    end

    @testset "LLR saturation" begin
        @test BPConfig(Flooding, 5, MinSum).llr_clip == 128f0
        @test BPConfig(Flooding, 5, NormMinSum, 0.8).llr_clip == 128f0
        @test parse_decoder(Dict()).llr_clip == 128f0
        @test parse_decoder(Dict("llr_clip" => 16)).llr_clip == 16f0
        @test Dict(BPsim.decoder_spec(parse_decoder(Dict("llr_clip" => 16))))["llr_clip"] == 16f0
        for value in (0, -1, Inf, -Inf, NaN, 1e100, 1e-100)
            @test_throws Exception parse_decoder(Dict("llr_clip" => value))
            @test_throws Exception parse_decoder(Dict("norm_factor" => value))
        end

        # Below saturation, retain the exact single-check extrinsic updates.
        graph = TannerGraph(Bool[1 1 1])
        for schedule in (Flooding, Layered)
            dec = BPDecoder(graph, BPConfig(schedule, 1, MinSum))
            bits = zeros(Bool, 3)
            @test decode!(dec, LLRFloat[-2, 3, 4], bits) == (Codeword, 1)
            @test bits == [false, false, false]
            @test dec.llr_out == LLRFloat[1, 1, 2]
            for value in (NaN, Inf, -Inf)
                @test_throws ArgumentError decode!(dec, LLRFloat[value, 1, 1], bits)
            end
        end

        # Unbounded NMS overflows on this graph and falsely accepts all zeros.
        graph = TannerGraph(trues(4, 4))
        for schedule in (Flooding, Layered), rule in (SPAExact, FastSPA, MinSum, NormMinSum), clip in (16f0, 128f0)
            dec = BPDecoder(graph, BPConfig(schedule, 200, rule, 0.7, clip))
            bits = zeros(Bool, 4)
            result, iterations = decode!(dec, LLRFloat[-1, 1, 1, 1], bits)
            @test all(isfinite, dec.llr_out)
            for buffer in (dec.msg_v2c, dec.msg_c2v, dec.α, dec.β)
                @test all(isfinite, buffer) && maximum(abs, buffer) <= clip
            end
            if schedule == Flooding && rule == NormMinSum
                @test (result, iterations) == (BPsim.Failure, 200)
                @test bits == [true, false, false, false]
            end
        end

        mktempdir() do dir
            path = repofile("codes/rm.32.16.mat")
            cfg = load_config(nothing, Dict("code_file" => path, "decoder" => Dict("llr_clip" => 16)))
            BPsim.write_sidecar(joinpath(dir, "run.yaml"), cfg, load_code(path), "test")
            @test BPsim.YAML.load_file(joinpath(dir, "run.yaml"))["decoder"]["llr_clip"] == 16
        end
    end

    @testset "config validation" begin
        @test_throws ErrorException parse_decoder(Dict("max_iter" => 0))
        @test parse_decoder(Dict("max_iter" => 10)).max_iter == 10
        @test_throws ErrorException parse_decoder(Dict("type" => "bp"))    # unknown key
    end
end
