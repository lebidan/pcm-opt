# Crash-resistant CSV logging: periodic checkpoint writes plus a self-describing sidecar.

const CSV_COLS = ["Eb/N0", "SNR", "Frames", "Frame Errors", "Info Bit Errors",
    "Code Bit Errors", "FER", "MLD FE", "MLD LB", "FAILS", "Info BER", "Code BER", "w_min",
    "Iter Min", "Iter Max", "Iter Avg"]
const CSV_INT_COLS = ("Frames", "Frame Errors", "Info Bit Errors", "Code Bit Errors", "MLD FE", "w_min",
    "Iter Min", "Iter Max")

# One CSV row, in CSV_COLS order; FAILS/MLD LB are fractions in [0,1]
csv_row(pt::PointReport, code::LinearCode) = (
    pt.ebn0, pt.snr, pt.c.n_cw, pt.c.n_cw_err, pt.c.n_bit_err, pt.c.n_code_bit_err,
    safe_ratio(pt.c.n_cw_err, pt.c.n_cw), pt.c.n_mld_err, safe_ratio(pt.c.n_mld_err, pt.c.n_cw),
    safe_ratio(pt.c.n_fail, pt.c.n_cw), safe_ratio(pt.c.n_bit_err, pt.c.n_cw * code.k),
    safe_ratio(pt.c.n_code_bit_err, pt.c.n_cw * code.n), pt.c.w_cw_min,
    pt.c.iters_min, pt.c.iters_max, safe_ratio(pt.c.iters_sum, pt.c.n_cw))

codeslug(path) = first(splitext(basename(path)))

function write_sidecar(path, cfg::SimConfig, code::LinearCode, ts)
    lines = [
        "# BP optimized-PCM simulation — parameters and run metadata",
        "# results are reproducible at a fixed nthreads (deterministic per-thread seeding",
        "# + round-synchronized batching); the thread count affects the total frame count.",
        "code_file: $(cfg.code_file)",
        "code_name: $(code.name)",
        "n: $(code.n)",
        "k: $(code.k)",
        "n_avn: $(code.n_avn)",
        "m: $(code.m)",
        "rate: $(code.rate)",
        "dmin: $(code.dmin === nothing ? "null" : code.dmin)",
        "use_base_pcm: $(cfg.use_base_pcm)",
        "decoder:",
        ("  $k: $v" for (k, v) in decoder_spec(cfg.decoder))...,
        "all_zero: $(cfg.all_zero)",
        "llr_scaling: $(cfg.llr_scaling)",
        "ebn0_min: $(cfg.ebn0_min)",
        "ebn0_max: $(cfg.ebn0_max)",
        "ebn0_step: $(cfg.ebn0_step)",
        "min_cw: $(cfg.min_cw)",
        "min_cw_err: $(cfg.min_cw_err)",
        "seed: $(cfg.seed)",
        "batch_size: $(cfg.batch_size)",
        "checkpoint_seconds: $(cfg.checkpoint_seconds)",
        "nthreads: $(Threads.nthreads())",
        "julia_version: \"$(VERSION)\"",
        "started: $ts",
    ]
    write(path, join(lines, "\n") * "\n")
end

function write_csv(path, df, extra_row)
    out = extra_row === nothing ? df : push!(copy(df), extra_row)
    CSV.write(path, out)
end

# --- CSV output: rewrites the in-progress row and appends committed rows, so an interrupted run
#     still leaves an up-to-date file ---

struct CsvReporter
    csv_path::String
    code::LinearCode
    df::DataFrame
end
CsvReporter(csv_path, code) =
    CsvReporter(csv_path, code, DataFrame([c => (c in CSV_INT_COLS ? Int[] : Float64[]) for c in CSV_COLS]))

report_update(r::CsvReporter, pt::PointReport) =
    write_csv(r.csv_path, r.df, csv_row(pt, r.code))

function report_point(r::CsvReporter, pt::PointReport)
    push!(r.df, csv_row(pt, r.code))
    write_csv(r.csv_path, r.df, nothing)
end
