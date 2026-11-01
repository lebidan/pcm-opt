# Extract an optimized sparse girth-6 PCM for a code and save it to a .mat file.
#
#   julia --project -t auto scripts/pcm_ii.jl codes/rm.32.16 [M_a=rows of H_s] [T_max=8000] [wcap=]
#                          [seed_s=1] [seed_p=6] [rapid=0] [out=results/rm.32.16.pcm2.mat]
#
# rapid=1 uses Algorithm 4 (faster, approximate). The reported PCM is the best over the M_a sweep.
# wcap sets Algorithm 2's weight cap (default: the lightest H row, raised until the collected
# codewords span rank n-k). Every dual codeword up to the cap enters H_s, so raising it past the
# spanning weight is how a code whose lightest dual codewords already span gets a redundant H_s.
#
# The output file carries every variable of the source code file (original types), with the source
# PCM renamed H_base and H holding the optimized PCM: its first n columns are the original VNs, the
# remaining ones are auxiliary (zero-LLR under BP) and listed in puncturing_set. Also: eta_min.

using Printf
using PCMOpt
include(joinpath(@__DIR__, "cli_helpers.jl"))

function main()
    code_path, o = parse_args(ARGS, ("M_a", "seed_p", "seed_s", "T_max", "wcap", "rapid", "out"))
    c = load_code(code_path)
    N, K = c.n, c.k
    T_max = geti(o, "T_max", 8000)
    wcap = geti(o, "wcap")
    seed_s = geti(o, "seed_s", 1)
    seed_p = geti(o, "seed_p", 6)

    @printf("Code: %s  (n=%d, k=%d, n-k=%d)\n", string(c), N, K, N - K)
    println("Threads: ", Threads.nthreads())

    t0 = time()
    Hs = gen_search_space(from_bool(c.H); T_max = T_max, seed = seed_s, wcap = wcap)
    wmin, wmax = extrema(rowweight(Hs, i) for i in 1:nrows(Hs))
    @printf("Search space H_s: %d×%d (weight %s), %.1fs\n", nrows(Hs), ncols(Hs),
            wmin == wmax ? "$wmin" : "$(wmin)–$(wmax)", time() - t0)

    M_a = geti(o, "M_a", nrows(Hs))
    rapid = geti(o, "rapid", 0)
    rapid in (0, 1) || error("rapid must be 0 or 1 (got $rapid)")
    t0 = time()
    res = extract_pcm(Hs, N, K; M_a = M_a, seed = seed_p, rapid = rapid == 1)
    el = time() - t0
    m = pcm_metrics(res.Hb, N)

    println("\nExtracted PCM-II (best over sweep):")
    @printf("  best at i   = %d rows of H_a (M_a = %d)\n", res.i, M_a)
    @printf("  size        = %d×%d\n", m.rows, m.cols)
    @printf("  edges       = %d  (density %.2f%%)\n", m.edges, 100 * m.density)
    @printf("  row weight  = %.2f   col weight = %.2f\n", m.rowwt, m.colwt)
    @printf("  AVNs        = %d\n", m.avn)
    @printf("  4-cycles    = %d   girth = %d   6-cycles = %d\n", m.q4, m.girth, m.sixcyc)
    @printf("  η̄_min       = %d\n", m.eta_min)
    @printf("  extract time= %.1fs\n", el)

    out = resolve_out(code_path, o, "pcm2.mat")
    save_optpcm(out, res.Hb, code_path)
    println("\nSaved: ", out)
end

main()
