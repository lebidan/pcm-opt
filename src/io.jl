# Save an optimized PCM. The output carries every variable of the source code file with its original
# type (copied verbatim from the file), except that the source PCM is renamed `H_base` and `H` now
# holds the optimized PCM (logical dense). Columns 1..n of `H` are the original VNs; the remaining
# ones are auxiliary (punctured) VNs, listed 1-based in `puncturing_set` (Int32). The connectivity
# metric `eta_min` (Int32) is also stored when `H` is 4-cycle-free; sizes and AVN counts are
# inferable from `H`.
function save_optpcm(out_path::AbstractString, A::GF2Matrix, src_path::AbstractString)
    src = endswith(src_path, ".mat") ? src_path : src_path * ".mat"
    d = matread(src)
    n = Int(d["n"])
    d["H_base"] = d["H"]
    d["H"] = to_bool(A)                          # logical, dense
    d["puncturing_set"] = Int32[n+1:ncols(A);]
    count_4cycles(A) == 0 && (d["eta_min"] = Int32(eta_min(A, n)))
    mkpath(dirname(out_path))
    matwrite(out_path, d)
    out_path
end

# Load a PCM from a .mat file for analysis: `H` is the decoding matrix in both cases — the optimized
# PCM in a produced file, the standard PCM in a plain code file. Returns (A::GF2Matrix, N) with
# N = #original VNs.
function load_optpcm(path::AbstractString)
    p = endswith(path, ".mat") ? path : path * ".mat"
    d = matread(p)
    from_bool(Matrix{Bool}(d["H"])), Int(d["n"])
end
