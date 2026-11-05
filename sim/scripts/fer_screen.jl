# Compare candidate PCMs and normalization factors by a short simulation, and print a FER table.
#
#   julia --project=. -t auto scripts/fer_screen.jl config.yaml <file.mat | dir>... [alpha=0.5,0.7] [key=value …]
#
# Each .mat file (or every .mat file of a directory) is simulated with the settings of config.yaml,
# once per value of alpha (decoder.norm_factor), with the key=value overrides of simulate.jl. Keep
# the run short: two Eb/N0 points and min_cw=50000 are usually enough to rank candidates.
# The simulation records go to logs/screen/.

using BPsim, CSV, DataFrames, Printf

function main()
    yaml = ARGS[1]
    files = String[]
    for a in ARGS[2:end]
        contains(a, '=') && continue
        isdir(a) ? append!(files, sort(filter(endswith(".mat"), readdir(a; join = true)); by = natural)) :
                   push!(files, a)
    end
    kv = BPsim.parse_overrides(filter(a -> contains(a, '=') && !startswith(a, "alpha="), ARGS[2:end]))
    alphas = [parse(Float64, s) for a in ARGS if startswith(a, "alpha=") for s in split(a[7:end], ',')]
    isempty(alphas) && (alphas = [0.7])

    for f in files, α in alphas
        o = merge(kv, Dict("code_file" => f, "output_dir" => "logs/screen"))
        o["decoder"] = merge(get(kv, "decoder", Dict{String,Any}()), Dict("norm_factor" => α))
        csv = redirect_stdout(() -> BPsim.run_simulation(BPsim.load_config(yaml, o)), devnull)
        df = CSV.read(csv, DataFrame)
        @printf("%-32s α=%.2f  %s\n", basename(f), α,
                join([@sprintf("%.2f dB: %.3e", e, p) for (e, p) in zip(df[!, "Eb/N0"], df.FER)], "   "))
    end
end

# B9.mat before B10.mat
natural(path) = [isnothing(tryparse(Int, s)) ? s : lpad(s, 8, '0') for s in split(basename(path), r"(?<=\D)(?=\d)|(?<=\d)(?=\D)")]

main()
