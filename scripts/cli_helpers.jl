# CLI helpers shared by the scripts: parse `<code> [key=value ...]` and resolve the output path.

function parse_args(args, allowed)
    isempty(args) && error("usage: <code> [key=value ...]")
    opts = Dict{String,String}()
    for a in args[2:end]
        kv = split(a, '=', limit = 2)
        length(kv) == 2 || error("bad option '$a' (expected key=value)")
        kv[1] in allowed || error("unknown option '$(kv[1])' (expected one of: $(join(allowed, ", ")))")
        opts[kv[1]] = kv[2]
    end
    first(args), opts
end

geti(opts, key, default) = haskey(opts, key) ? parse(Int, opts[key]) : default

# For options whose "unset" is `nothing` rather than a number (the callee picks the default).
geti(opts, key) = haskey(opts, key) ? parse(Int, opts[key]) : nothing

# results/<code stem>.<suffix>, unless overridden by opts["out"]. Only a trailing .mat is stripped
# (the code names themselves contain dots, e.g. rm.32.16).
function resolve_out(code_path, opts, suffix)
    base = basename(code_path)
    endswith(base, ".mat") && (base = base[1:end-4])
    get(opts, "out", joinpath(@__DIR__, "..", "results", "$base.$suffix"))
end
