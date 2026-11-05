# BPsim: BP simulator for optimized PCMs

BPsim measures the frame and bit error rates of belief propagation (BP) decoding by Monte-Carlo
simulation, with BPSK over the BI-AWGN channel. It decodes the PCMs produced by the scripts of this
repository, including those with auxiliary variable nodes. It is multithreaded. It is adapted from
my own BP simulator, reduced to what this repository needs.

## Installation

From the repository root:

```
julia --project=sim -e 'using Pkg; Pkg.instantiate()'
```

## Auxiliary variable nodes

An optimized PCM `H` has `n + n_avn` columns: the `n` code bits, followed by `n_avn` auxiliary VNs
(listed in `puncturing_set`). Auxiliary VNs are not transmitted, so the simulator gives them a zero
channel LLR. BP runs on the whole Tanner graph, and the decoded message is recovered from the first
`n` bits.

## Usage

From the `sim/` directory:

```
julia --project=. -t auto scripts/simulate.jl [config.yaml] [key=value …]
```

Settings come from the YAML file, and `key=value` arguments override them. Nested settings use a dot,
e.g. `decoder.max_iter=50`. The config file is optional if `code_file` is given on the command line:

```
julia --project=. -t auto scripts/simulate.jl config/rm.32.16.pcm2.yaml                       # optimized PCM
julia --project=. -t auto scripts/simulate.jl config/rm.32.16.pcm2.yaml use_base_pcm=true     # original PCM
julia --project=. -t auto scripts/simulate.jl code_file=../codes/bch.63.45.mat ebn0_max=6 decoder.algorithm=MS
```

To compare several PCMs, or several values of α, by a short simulation, `scripts/fer_screen.jl`
runs `simulate.jl` on each `.mat` file given (or each one in a directory), once per value of α, and
prints one line of FER per run:

```
julia --project=. -t auto scripts/fer_screen.jl config/rm.32.16.sspcm.yaml ../sweep/rm.32.16 alpha=0.5,0.7 min_cw=50000
```

`sim/config/` contains one config for each optimized PCM in `results/`, except RM(8,4). They use the decoder
settings of the paper. The [reference records](logs/README.md) include CSV results and run metadata;
use the recorded thread count for exact replay. See the [top-level README](../README.md) for the results.

### Settings

| Key | Default | Meaning |
|---|---|---|
| `code_file` | (required) | Code or optimized PCM (`.mat`), relative to the current directory |
| `use_base_pcm` | `false` | Decode with the original PCM `H_base` instead of the optimized `H` |
| `ebn0_min`, `ebn0_max`, `ebn0_step` | `0`, `4`, `1` | Eb/N0 points, in dB |
| `min_cw` | `100000` | Minimum number of frames per point |
| `min_cw_err` | `100` | Minimum number of frame errors per point |
| `all_zero` | `false` | Transmit the all-zero codeword instead of random codewords |
| `llr_scaling` | `true` | Feed the decoder with LLRs `2y/σ²` (`false`: raw channel samples `y`) |
| `seed` | `8731` | Random seed |
| `batch_size` | `1000` | Frames simulated by each thread between two checks of the stopping rule |
| `checkpoint_seconds` | `5` | Interval between updates of the console and the CSV file, in seconds |
| `output_dir` | `logs` | Output directory |
| `decoder.algorithm` | `FastSPA` | See below |
| `decoder.schedule` | `layered` | `flooding` or `layered` |
| `decoder.max_iter` | `32` | Maximum number of iterations |
| `decoder.norm_factor` | `0.7` | Normalization factor α of NMS |
| `decoder.llr_clip` | `128` | LLR saturation limit |

A point stops when both `min_cw` and `min_cw_err` are reached. The stopping rule is checked after
each round of `batch_size` frames on every thread, so the number of frames of a point is a multiple
of `batch_size` × number of threads, possibly well above `min_cw`. Each thread draws its own noise, so a simulation is reproducible for a given `seed`,
`batch_size` and number of threads. With another number of threads, the error rates differ within
the usual Monte-Carlo accuracy. The decoder enforces LLR saturation because messages can overflow in
`Float32` arithmetic, even with normalized min-sum, when using many iterations with high-degree variable nodes.

### Decoders

| `decoder.algorithm` | Check-node update |
|---|---|
| `SPA` | Sum-product (exact box-plus) |
| `FastSPA` | Sum-product, with Richter's approximation of the box-plus |
| `MS` | Min-sum |
| `NMS` | Normalized min-sum: min-sum scaled by α (`decoder.norm_factor`) |

Decoding stops as soon as the hard decision is a codeword.

### Output

Each run writes to `output_dir` a CSV file with one line per Eb/N0 point, updated every few seconds
so that an interrupted run keeps its results, and a YAML file with the full settings and code
parameters. The main columns are:

| Column | Meaning |
|---|---|
| `FER`, `Info BER`, `Code BER` | Frame error rate, and bit error rates on the message and on the codeword |
| `FAILS` | Fraction of frames where BP did not converge to a codeword |
| `MLD FE` | Frame errors where BP returned a codeword more likely than the one transmitted, so that an ML decoder would fail as well |
| `MLD LB` | `MLD FE` / frames: a lower bound on the ML frame error rate |
| `w_min` | Smallest Hamming distance observed between a wrong decoded codeword and the transmitted one (`n` if there was none) |
| `Iter Min`, `Iter Max`, `Iter Avg` | Number of BP iterations |

## Tests

```
julia --project=. -t 4 test/runtests.jl
```
