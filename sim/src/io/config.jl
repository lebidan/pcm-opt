struct SimConfig
    code_file::String
    decoder::BPConfig
    use_base_pcm::Bool      # decode with the base PCM H_base (n_avn=0) instead of H — baseline curve
    all_zero::Bool
    llr_scaling::Bool
    ebn0_min::Float64
    ebn0_max::Float64
    ebn0_step::Float64
    min_cw::Int
    min_cw_err::Int
    seed::Int
    batch_size::Int
    checkpoint_seconds::Float64
    output_dir::String
end

ebn0_points(cfg::SimConfig) = cfg.ebn0_min:cfg.ebn0_step:cfg.ebn0_max

const SIM_DEFAULTS = Dict{String,Any}(
    "use_base_pcm" => false,
    "all_zero" => false,
    "llr_scaling" => true,
    "ebn0_min" => 0.0,
    "ebn0_max" => 4.0,
    "ebn0_step" => 1.0,
    "min_cw" => 100000,
    "min_cw_err" => 100,
    "seed" => 8731,
    "batch_size" => 1000,
    "checkpoint_seconds" => 5.0,
    "output_dir" => "logs",
)

# --- parsing helpers ---

parse_enum(names, s, what) = get(() -> error("unknown $what '$s'"), names, uppercase(String(s)))
function check_keys(d, known, ctx)
    for k in keys(d)
        String(k) in known ||
            error("unknown $ctx key '$k' (known: $(join(sort!(collect(known)), ", ")))")
    end
end

# --- the `decoder:` settings ---

const SCHEDULE_NAMES = Dict("LAYERED"=>Layered, "FLOODING"=>Flooding)
const ALGORITHM_NAMES = Dict("SPA"=>SPAExact, "FASTSPA"=>FastSPA, "MS"=>MinSum, "NMS"=>NormMinSum)

parse_algorithm(s) = parse_enum(ALGORITHM_NAMES, s, "algorithm")

const DECODER_KEYS = Set(["algorithm", "schedule", "max_iter", "norm_factor"])
function parse_decoder(d)
    check_keys(d, DECODER_KEYS, "decoder")
    max_iter = Int(get(d, "max_iter", 32))
    max_iter >= 1 || error("decoder max_iter must be ≥ 1 (got $max_iter)")
    norm_factor = Float64(get(d, "norm_factor", 0.7))
    norm_factor > 0 || error("decoder norm_factor must be > 0 (got $norm_factor)")
    return BPConfig(
        parse_enum(SCHEDULE_NAMES, get(d, "schedule", "layered"), "schedule"),
        max_iter,
        parse_algorithm(get(d, "algorithm", "FastSPA")),
        LLRFloat(norm_factor),
    )
end

# --- config assembly: defaults < YAML < CLI overrides, all as nested dicts ---

function merge_config!(dst, src)
    for (k, v) in src
        if v isa AbstractDict && get(dst, k, nothing) isa AbstractDict
            merge_config!(dst[k], v)
        else
            dst[k] = v
        end
    end
    return dst
end

const KNOWN_TOP = Set(String.(fieldnames(SimConfig)))
function load_config(yaml_path::Union{Nothing,AbstractString}, overrides::AbstractDict=Dict{String,Any}())
    cfg = deepcopy(SIM_DEFAULTS)
    yaml_path === nothing || merge_config!(cfg, YAML.load_file(yaml_path))
    merge_config!(cfg, deepcopy(overrides))
    haskey(cfg, "code_file") || error("`code_file` must be set (config file or code_file=…)")
    check_keys(cfg, KNOWN_TOP, "config")
    ebn0_min, ebn0_max, ebn0_step = Float64(cfg["ebn0_min"]), Float64(cfg["ebn0_max"]), Float64(cfg["ebn0_step"])
    ebn0_step > 0 || error("`ebn0_step` must be > 0 (got $ebn0_step)")
    ebn0_max >= ebn0_min || error("`ebn0_max` ($ebn0_max) must be ≥ `ebn0_min` ($ebn0_min)")
    return SimConfig(
        String(cfg["code_file"]),
        parse_decoder(get(cfg, "decoder", Dict{String,Any}())),
        Bool(cfg["use_base_pcm"]),
        Bool(cfg["all_zero"]),
        Bool(cfg["llr_scaling"]),
        ebn0_min, ebn0_max, ebn0_step,
        Int(cfg["min_cw"]),
        Int(cfg["min_cw_err"]),
        Int(cfg["seed"]),
        Int(cfg["batch_size"]),
        Float64(cfg["checkpoint_seconds"]),
        String(cfg["output_dir"]),
    )
end

# --- CLI: optional YAML path followed by Hydra-style dotted overrides (decoder.max_iter=30) ---

const USAGE = """
Usage: simulate.jl [config.yaml] [key=value ...]

Dotted keys set nested fields; the config file is optional. Examples:
  simulate.jl config/rm.32.16.pcm2.yaml decoder.max_iter=30 min_cw=5000
  simulate.jl config/rm.32.16.pcm2.yaml use_base_pcm=true            # baseline (base PCM) curve
  simulate.jl code_file=../results/bch.63.45.pcm2.mat ebn0_min=1 ebn0_max=5 decoder.algorithm=MS
"""

function parse_value(s::AbstractString)
    s in ("true", "false") && return s == "true"
    for T in (Int, Float64)
        v = tryparse(T, s)
        v === nothing || return v
    end
    return String(s)
end

function set_nested!(d, path, value)
    for k in @view path[1:end-1]
        d = get!(() -> Dict{String,Any}(), d, String(k))
    end
    d[String(path[end])] = value
end

function parse_overrides(tokens)
    overrides = Dict{String,Any}()
    for t in tokens
        occursin('=', t) || error("expected key=value, got '$t'")
        key, value = split(t, '=', limit=2)
        set_nested!(overrides, split(key, '.'), parse_value(value))
    end
    return overrides
end

function parse_commandline(argv)
    args = collect(String, argv)
    (isempty(args) || args[1] in ("-h", "--help")) && (print(USAGE); exit(0))
    config = occursin('=', args[1]) ? nothing : popfirst!(args)
    return config, parse_overrides(args)
end
