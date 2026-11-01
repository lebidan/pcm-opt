# Persistent per-worker state: own decoder, RNG (deterministically seeded) and buffers.
struct WorkerState
    dec::BPDecoder
    rng::Xoshiro
    data::Vector{Bool}
    msg_hat::Vector{Bool}
    tx_cw::Vector{Bool}
    rx_cw::Vector{Bool}
    llr::Vector{LLRFloat}
end
# The decoder-facing buffers span all Ñ = n + n_avn variable nodes; `tx_cw`/`data` stay on the
# n channel bits. `llr` is zero-initialized, and since `transmit!` writes only its first n
# entries and BP treats `llr` as read-only, the auxiliary tail llr[n+1:Ñ] stays 0 for the whole
# run — exactly the zero channel LLR the auxiliary VNs require.
function WorkerState(code::LinearCode, graph::TannerGraph, cfg::SimConfig, snr_index::Int, tid::Int)
    n_vn = code.n + code.n_avn
    WorkerState(BPDecoder(graph, cfg.decoder), Xoshiro(hash((cfg.seed, snr_index, tid))),
        zeros(Bool, code.k), zeros(Bool, code.k), zeros(Bool, code.n), zeros(Bool, n_vn), zeros(LLRFloat, n_vn))
end

function run_batch!(ws::WorkerState, c::Counters, code::LinearCode, cfg::SimConfig, σ::Float64)
    G, Ginv = code.G, code.Ginv
    for _ in 1:cfg.batch_size
        # all-zero mode: leave data/tx_cw at their zero-initialized state (no encoding)
        cfg.all_zero || encode!(ws.rng, ws.data, ws.tx_cw, G)
        transmit!(ws.rng, ws.tx_cw, ws.llr, σ, cfg.llr_scaling)
        res, iters = decode!(ws.dec, ws.llr, ws.rx_cw)
        n_err_total, n_err_data = count_errors!(ws.msg_hat, Ginv, ws.data, ws.tx_cw, ws.rx_cw)
        c.n_cw += 1
        c.n_bit_err += n_err_data
        c.n_code_bit_err += n_err_total
        res == Failure && (c.n_fail += 1)
        c.iters_sum += iters
        c.iters_min = min(c.iters_min, iters)
        c.iters_max = max(c.iters_max, iters)
        if n_err_total > 0
            c.n_cw_err += 1
            if res == Codeword
                n_err_total < c.w_cw_min && (c.w_cw_min = n_err_total)
                corr_metric(ws.llr, ws.rx_cw) > corr_metric(ws.llr, ws.tx_cw) && (c.n_mld_err += 1)
            end
        end
    end
end

# Round-synchronized batching: every round all workers run exactly one batch, then we
# aggregate and test the stopping criterion. Each worker runs the same number of batches
# in every run, so — with per-worker deterministic seeding — the result is reproducible
# at a fixed thread count (independent of scheduling jitter).
function run_snr_point!(code::LinearCode, graph::TannerGraph, cfg::SimConfig, snr_index::Int, σ::Float64, on_checkpoint)
    nt = Threads.nthreads()
    counters = [Counters(code.n) for _ in 1:nt]
    workers = [WorkerState(code, graph, cfg, snr_index, t) for t in 1:nt]
    t_point = time()
    last_ckpt = t_point
    tot_cw = tot_cwe = 0
    while !(tot_cw ≥ cfg.min_cw && tot_cwe ≥ cfg.min_cw_err)
        @sync for t in 1:nt
            Threads.@spawn run_batch!(workers[t], counters[t], code, cfg, σ)
        end
        tot_cw  = sum(c.n_cw for c in counters)
        tot_cwe = sum(c.n_cw_err for c in counters)
        if time() - last_ckpt ≥ cfg.checkpoint_seconds
            on_checkpoint(reduce_counters(counters), time() - t_point)
            last_ckpt = time()
        end
    end
    return reduce_counters(counters), time() - t_point
end

# Loop over the Eb/N0 points; the console and the CSV file are updated at every checkpoint and
# at the end of each point.
function run_simulation(cfg::SimConfig)
    code = load_code(cfg.code_file; use_base_pcm=cfg.use_base_pcm)
    graph = TannerGraph(code.H)

    mkpath(cfg.output_dir)
    ts = Dates.format(Dates.now(), "yyyy-mm-dd-HH-MM-SS")
    base = joinpath(cfg.output_dir, "sim-$(codeslug(cfg.code_file))-$(decoder_name(cfg.decoder))-$ts")
    csv_path = base * ".csv"
    console, csv = ConsoleReporter(cfg, code), CsvReporter(csv_path, code)
    write_sidecar(base * ".yaml", cfg, code, ts)
    print_sim_info(cfg, code)
    println()
    println("Output file: $csv_path")
    println()
    println("The simulation is running...")
    println()
    print_table_header()

    t_start = time()
    w_min_global = typemax(Int)
    for (si, ebn0) in enumerate(ebn0_points(cfg))
        σ = sqrt(10^(-ebn0/10) / (2 * code.rate))
        snr = ebn0 + 10*log10(2 * code.rate)   # BPSK: SNR = 2·Es/N0 = 1/σ² (Es/N0 = R·Eb/N0)

        on_ckpt = (c::Counters, elapsed) -> begin
            pt = PointReport(ebn0, snr, c, elapsed)
            report_update(console, pt)
            report_update(csv, pt)
        end

        c, elapsed = run_snr_point!(code, graph, cfg, si, σ, on_ckpt)
        pt = PointReport(ebn0, snr, c, elapsed)
        report_point(console, pt)
        report_point(csv, pt)
        w_min_global = min(w_min_global, c.w_cw_min)
    end

    report_finish(console, time() - t_start, w_min_global)
    return csv_path
end
