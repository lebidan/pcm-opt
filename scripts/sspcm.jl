# Build structured sparse PCMs (ssPCMs) for a code and save them to .mat files.
#
#   julia --project -t auto scripts/sspcm.jl codes/<code> [wcap=…] [T_max=8000] [seed_s=1]
#                          [tmin=2] [avnfree=1] [max_edges=4000] [seed=1] [B=…] [out=…]
#
# Without B, saves every full-rank matrix of the nested sequence to sweep/<code>/B<blocks>.mat and
# prints its metrics: pick one by simulation with sim/scripts/fer_screen.jl. With B, stops at B blocks
# and saves that matrix to results/<code>.sspcm.mat (max_edges is then ignored).
# wcap, T_max and seed_s set the search space (Algorithm 2, as in pcm_ii.jl). tmin is the smallest
# anchor, avnfree=0 forbids AVN-free blocks, and seed breaks ties between equally good blocks.

using Printf
using PCMOpt
include(joinpath(@__DIR__, "cli_helpers.jl"))

function main()
    code_path, o = parse_args(ARGS, ("wcap", "T_max", "seed_s", "tmin", "avnfree", "max_edges", "seed", "B", "out"))
    c = load_code(code_path)
    N, K = c.n, c.k
    B = geti(o, "B")
    stem = replace(basename(code_path), r"\.mat$" => "")
    @printf("Code: %s  (n=%d, k=%d), %d threads\n", string(c), N, K, Threads.nthreads())

    t0 = time()
    Hs = gen_search_space(from_bool(c.H); T_max = geti(o, "T_max", 8000), seed = geti(o, "seed_s", 1),
                          wcap = geti(o, "wcap"))
    wts = [count(r) for r in Hs.rows]
    @printf("Search space: %d rows (%s), %.0fs\n", nrows(Hs),
            join(["$(count(==(w), wts))×w$w" for w in sort(unique(wts))], ", "), time() - t0)

    dir = get(o, "out", joinpath(@__DIR__, "..", "sweep", stem))
    B === nothing && mkpath(dir)
    function onstep(A, Γ, nb)
        @assert represents_code(A, c.G, N, K)
        @printf("B=%3d  %3d×%3d  AVNs=%3d  edges=%4d  4-cycles=%6d  Γmax=%d\n", nb, nrows(A),
                ncols(A), ncols(A) - N, num_edges(A), count_4cycles(A), Γ)
        if B === nothing
            save_optpcm(joinpath(dir, "B$nb.mat"), A, code_path)
        elseif nb == B
            out = get(o, "out", joinpath(@__DIR__, "..", "results", "$stem.sspcm.mat"))
            save_optpcm(out, A, code_path)
            println("Saved: ", out)
        end
    end
    sspcm(Hs, N, K; tmin = geti(o, "tmin", 2), avnfree = geti(o, "avnfree", 1) == 1,
          max_edges = B === nothing ? geti(o, "max_edges", 4000) : typemax(Int),
          max_blocks = something(B, typemax(Int)), seed = geti(o, "seed", 1), onstep)
    B === nothing && println("Saved the sequence in ", dir)
end

main()
