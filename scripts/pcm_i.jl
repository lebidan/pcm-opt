# Extract a 4-cycle-free PCM-I (Algorithm 1) from a code's PCM and save it to a .mat file.
#
#   julia --project [-t auto] scripts/pcm_i.jl <code> [seed=X | n_seed=M] [out=results/<code>.pcm1.mat]
#
# Three modes (seed and n_seed are mutually exclusive):
#   (default)  deterministic Algorithm 1 (largest block, then last found, among equal δ).
#   seed=X     one randomized-tie-break run seeded with X (reproducible; may differ from default).
#   n_seed=M   sweep seeds 1..M, keep the best PCM-I found. The deterministic result is folded in as
#              the baseline candidate, so the saved matrix is never worse than the default.
#
# "Best" is the lexicographic minimum of (AVNs, −η̄_min, 6-cycles, edges): fewest auxiliary VNs, then
# largest connectivity η̄_min, then fewest 6-cycles, then fewest edges. In n_seed mode each new best
# is logged with its seed so any listed candidate can be regenerated with seed=<that seed>.

using Printf
using PCMOpt
using Random
include(joinpath(@__DIR__, "cli_helpers.jl"))

# η̄_min = −1 means no AVN sits on any 6-cycle (best case) → sort ahead of every non-negative value.
eta_cost(eta) = eta == -1 ? typemin(Int) : -eta
sortkey(r) = (r.avn, eta_cost(r.eta_min), r.sixcyc, r.edges)   # r: (seed, avn, eta_min, sixcyc, edges)

metrics(A, N) = (avn = ncols(A) - N, eta_min = eta_min(A, N),
                 sixcyc = count_6cycles(A), edges = num_edges(A))

function report(label, r)
    @printf("%-14s AVNs=%d  η̄_min=%d  6-cycles=%d  edges=%d\n",
            label, r.avn, r.eta_min, r.sixcyc, r.edges)
end

function main()
    code_path, opts = parse_args(ARGS, ("seed", "n_seed", "out"))
    haskey(opts, "seed") && haskey(opts, "n_seed") &&
        error("seed and n_seed are mutually exclusive")
    c = load_code(code_path)
    H = from_bool(c.H)

    @printf("Code: %s - Input H has shape %d×%d\n", string(c), size(c.H, 1), size(c.H, 2))
    @printf("Using %d thread(s)\n", Threads.nthreads())
    println("Starting greedy 4-cycles elimination process (Alg. 1)")

    if haskey(opts, "seed")
        seed = geti(opts, "seed", 0)
        A = remove_4cycles_rand(H, Xoshiro(seed))
        report("seed $seed:", metrics(A, c.n))
    elseif haskey(opts, "n_seed")
        n_seed = geti(opts, "n_seed", 0)
        n_seed >= 1 || error("n_seed must be ≥ 1")
        A = search_best(H, c.n, n_seed)
    else
        A = remove_4cycles(H)
        report("deterministic:", metrics(A, c.n))
    end

    m = pcm_metrics(A, c.n)
    @assert m.q4 == 0 && m.girth == 6
    out = resolve_out(code_path, opts, "pcm1.mat")
    save_optpcm(out, A, code_path)
    println("Saved: ", out)
end

# Sweep seeds 1..n_seed for the best randomized PCM-I, with the deterministic result as the baseline.
# Threads collect immutable records only; the winner is recomputed serially from its seed (matrices
# are never shared across threads). The improvement log is folded over seed order so it — and the
# winner among key-ties — is independent of thread timing.
function search_best(H, N, n_seed)
    det = merge((seed = 0,), metrics(remove_4cycles(H), N))   # seed 0 labels the deterministic run
    nt = Threads.nthreads()
    recs = [typeof(det)[] for _ in 1:nt]
    Threads.@threads for t in 1:nt
        for s in t:nt:n_seed
            push!(recs[t], merge((seed = s,), metrics(remove_4cycles_rand(H, Xoshiro(s)), N)))
        end
    end
    ordered = sort(reduce(vcat, recs), by = r -> r.seed)

    best = det
    report("deterministic:", det)
    for r in ordered
        if sortkey(r) < sortkey(best)
            best = r
            report("seed $(r.seed):", r)
        end
    end
    @printf("Best of %d seeds + deterministic: seed %s\n",
            n_seed, best.seed == 0 ? "deterministic" : string(best.seed))
    best.seed == 0 ? remove_4cycles(H) : remove_4cycles_rand(H, Xoshiro(best.seed))
end

main()
