# Console output: a grouped parameters header + a live, in-place box-drawing results table.

# Results table columns as (title, unit, width). The same layout drives the header rows and
# every data row so they cannot drift.
const COLS = [("SNR", "(dB)", 6), ("Eb/N0", "(dB)", 6),
    ("FRA", "", 9), ("FE", "", 7), ("MLD FE", "", 7), ("iBE", "", 7), ("cBE", "", 7),
    ("FER", "", 9), ("FAILS", "", 8), ("MLD LB", "", 9), ("iBER", "", 9), ("cBER", "", 9),
    ("ITER", "(avg)", 6), ("C_THR", "(Mb/s)", 9), ("RT/ET", "(hhmmss)", 9)]

# Each cell gets one space of padding inside its `│` bars, so a column spans `w + 2` chars
# between junctions; that keeps the border segments below aligned with the data cells.
render_row(cells)      = "  │ " * join(cells, " │ ") * " │"
render_border(l, m, r) = "  " * l * join(["─"^(w + 2) for (_, _, w) in COLS], m) * r

fmt_hms(s::Real) = isfinite(s) ?
    (t = round(Int, s); @sprintf("%02dh%02d'%02d", t ÷ 3600, (t % 3600) ÷ 60, t % 60)) : "--h--'--"

# Estimate remaining time for the running point: it stops once BOTH targets (min_cw frames and
# min_cw_err frame errors) are met, so the less-satisfied one is binding. Extrapolate linearly.
function estimate_remaining(c::Counters, elapsed, cfg::SimConfig)
    frac = min(safe_ratio(c.n_cw, cfg.min_cw), safe_ratio(c.n_cw_err, cfg.min_cw_err))
    frac ≤ 0 ? Inf : elapsed * (1 / frac - 1)
end
rate_str(k, n) = (g = gcd(k, n); @sprintf("%.3f (%d/%d)", k / n, k ÷ g, n ÷ g))
ebn0_range_str(cfg::SimConfig) = (r = ebn0_points(cfg);
    @sprintf("%g:%g:%g (%d points)", first(r), step(r), last(r), length(r)))

const SECTION_WIDTH = 45
section(title) = "  " * title * " " * "─"^(SECTION_WIDTH - length(title) - 1)
param_line(label, value) = "    " * rpad(label, 21) * string(value)

function print_sim_info(cfg::SimConfig, code::LinearCode)
    rule = "  " * "─"^SECTION_WIDTH
    println(rule)
    println("  BP Decoder — Optimized-PCM Performance Benchmark")
    println(rule)
    println()
    println(section("Simulation"))
    for (l, v) in ["Seed" => cfg.seed, "Multi-threading (t)" => "$(Threads.nthreads()) thread(s)",
        "Target frame errors" => cfg.min_cw_err, "Min frames per point" => cfg.min_cw,
        "Batch size" => cfg.batch_size, "Refresh (s)" => cfg.checkpoint_seconds,
        "LLR precision" => string(LLRFloat), "Eb/N0 range (dB)" => ebn0_range_str(cfg)]
        println(param_line(l, v))
    end
    println()
    println(section("Source"))
    println(param_line("Type", cfg.all_zero ? "all-zero (encoding skipped)" : "random codewords"))
    println()
    println(section("Code"))
    for (l, v) in ["Name" => code.name, "File" => cfg.code_file, "Info. bits (K)" => code.k,
        "Codeword size (N)" => code.n, "Parity checks (M)" => code.m,
        "Aux. VNs (n_avn)" => code.n_avn, "Decoding VNs (Ñ)" => code.n + code.n_avn,
        "Code rate (R)" => rate_str(code.k, code.n),
        "Min. distance (dmin)" => code.dmin === nothing ? "n/a" : code.dmin]
        println(param_line(l, v))
    end
    println()
    println(section("Channel"))
    println(param_line("Type", "BI-AWGN"))
    println(param_line("Modulation", "BPSK"))
    println(param_line("LLR scaling", cfg.llr_scaling ? "on  (2/σ²)" : "off (raw samples)"))
    println()
    println(section("Decoder"))
    for (l, v) in decoder_params(cfg.decoder)
        println(param_line(l, v))
    end
end

# Top border, title row, units row, then the header/body separator. The bottom border is
# printed by report_finish, since data rows stream in one at a time below it.
function print_table_header()
    println(render_border("┌", "┬", "┐"))
    println(render_row([lpad(t, w) for (t, _, w) in COLS]))
    println(render_row([lpad(u, w) for (_, u, w) in COLS]))
    println(render_border("├", "┼", "┤"))
end

# While a point runs, the last column shows estimated remaining time (RT); on the committed
# row `remaining === nothing`, so it shows the point's elapsed time (ET) instead — hence "RT/ET".
function format_row(snr, ebn0, c::Counters, elapsed, code::LinearCode; remaining=nothing)
    cthr = elapsed > 0 ? code.n * c.n_cw / elapsed / 1e6 : 0.0
    values = (@sprintf("%.2f", snr), @sprintf("%.2f", ebn0),
        string(c.n_cw), string(c.n_cw_err), string(c.n_mld_err),
        string(c.n_bit_err), string(c.n_code_bit_err),
        @sprintf("%.2e", safe_ratio(c.n_cw_err, c.n_cw)),
        @sprintf("%.2e", safe_ratio(c.n_fail, c.n_cw)),
        @sprintf("%.2e", safe_ratio(c.n_mld_err, c.n_cw)),
        @sprintf("%.2e", safe_ratio(c.n_bit_err, c.n_cw * code.k)),
        @sprintf("%.2e", safe_ratio(c.n_code_bit_err, c.n_cw * code.n)),
        @sprintf("%.2f", safe_ratio(c.iters_sum, c.n_cw)), @sprintf("%.2f", cthr), fmt_hms(remaining === nothing ? elapsed : remaining))
    return render_row([lpad(values[i], COLS[i][3]) for i in eachindex(values)])
end

# --- The live table: each checkpoint rewrites the current row in place, the end of a point commits it ---

struct ConsoleReporter
    cfg::SimConfig
    code::LinearCode
end

function report_update(r::ConsoleReporter, pt::PointReport)
    print("\r", format_row(pt.snr, pt.ebn0, pt.c, pt.elapsed, r.code;
        remaining=estimate_remaining(pt.c, pt.elapsed, r.cfg)))
    flush(stdout)
end

report_point(r::ConsoleReporter, pt::PointReport) =
    println("\r", format_row(pt.snr, pt.ebn0, pt.c, pt.elapsed, r.code))

function report_finish(r::ConsoleReporter, wall_clock, w_min)
    println(render_border("└", "┴", "┘"))
    wmin_str = w_min < r.code.n ? "$w_min   (upper bound on the true code dmin)" : "none found"
    println()
    println("End of the simulation.")
    println()
    println("  Total wall-clock  = $(fmt_hms(wall_clock))")
    println("  Min. weight found = $wmin_str")
    println()
end
