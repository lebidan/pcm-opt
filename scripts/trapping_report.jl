# Report the dominant trapping-set distribution 𝒯(a,b) of a PCM.
#
#   julia --project scripts/trapping_report.jl <file.mat> [a_max=4]
#
# <file.mat> may be a produced PCM or a code file; both are read from their `H` variable.
# 𝒯(a,b) = number of size-a variable-node subsets whose induced subgraph is a *dominant* trapping set
# per Karimi & Banihashemi [50, §II]: connected, contains a cycle, and every VN touches ≥2 satisfied
# (even-degree) check nodes (degree-2 VNs exempt); b = number of unsatisfied (odd-degree) check nodes.
# See src/trapping.jl for the full definition. The minimum trapping set — smallest a, then smallest b
# — is the harmfulness metric the paper reports (Tables I, II).

using Printf
using PCMOpt
include(joinpath(@__DIR__, "cli_helpers.jl"))

function main()
    path, opts = parse_args(ARGS, ("a_max",))
    a_max = geti(opts, "a_max", 4)
    A, N = load_optpcm(path)

    @printf("PCM: %s\n", basename(path))
    @printf("  size  = %d×%d  (N=%d VNs, %d AVNs)   girth = %d\n",
            nrows(A), ncols(A), N, ncols(A) - N, girth(A))

    dist = trapping_distribution(A; a_max = a_max)
    println("\nDominant trapping sets 𝒯(a,b)  (a ≤ ", a_max, "):")
    if isempty(dist)
        println("  (none)")
    else
        for k in sort!(collect(keys(dist)))
            @printf("  𝒯(%d,%2d) = %d\n", k[1], k[2], dist[k])
        end
        mts = min_trapping_set(dist)
        @printf("\nMinimum trapping set: 𝒯(%d,%d)\n", mts[1], mts[2])
    end
end

main()
