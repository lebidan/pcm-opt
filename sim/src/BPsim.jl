module BPsim

using MAT, LinearAlgebra, Random, Printf, Dates
using DataFrames, CSV, YAML

# Single message/LLR precision: Float32 halves memory and widens SIMD vs Float64
const LLRFloat = Float32

@enum DecResult Codeword Failure

include("decoders/tanner.jl")
include("core/linear.jl")
include("core/channel.jl")
include("decoders/bp.jl")
include("core/counters.jl")
include("io/config.jl")
include("io/console.jl")
include("io/csv.jl")
include("core/engine.jl")

export LLRFloat, Codeword, decode!, decoder_name, TannerGraph, is_codeword,
    CNRule, SPAExact, FastSPA, MinSum, NormMinSum, LdpcSchedule, Layered, Flooding, BPConfig, BPDecoder,
    LinearCode, load_code, load_config, parse_commandline, run_simulation

end
