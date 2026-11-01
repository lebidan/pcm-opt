module PCMOpt

using MAT, Random, Printf
using Combinatorics: combinations

include("linear.jl")
include("gf2.jl")
include("cycles.jl")
include("trapping.jl")
include("io.jl")
include("four_cycle.jl")
include("search_space.jl")
include("prg.jl")

export LinearCode, load_code,
    GF2Matrix, nrows, ncols, rowweight, rank_gf2, from_bool,
    count_4cycles, count_6cycles, girth, eta_min, num_edges, pcm_metrics,
    trapping_distribution, min_trapping_set,
    save_optpcm, load_optpcm,
    remove_4cycles, remove_4cycles_rand, rapid_bypass,
    gen_search_space,
    extract_pcm

end
