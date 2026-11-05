# PCM optimization for BP decoding of short block codes

Belief propagation (BP) is excellent for decoding long LDPC codes, but performs poorly on short
non-LDPC codes such as RM, BCH, eBCH or polar codes. Their standard parity-check matrices (PCMs) are
often dense, which leads to Tanner graphs riddled with 4-cycles and small trapping sets. In 2025,
Yifei Shen and co-authors proposed a general method to extract from the PCM of a short binary block
code an optimized sparse PCM with enhanced short-cycle connectivity, on which BP decoding performs
remarkably well:

> Y. Shen, Z. Li, Y. Ren, E. Boutillon, A. Balatsoukas-Stimming, C. Zhang, X. You, A. Burg,
> "Toward Universal Belief Propagation Decoding for Short Binary Block Codes,"
> *IEEE Journal on Selected Areas in Communications*, vol. 43, no. 4, pp. 1135–1152, April 2025.
> [IEEE Xplore](https://ieeexplore.ieee.org/document/10887536)

The authors then refined the method, building the sparse PCM from structured blocks of parity checks:

> Y. Shen, Z. Li, E. Boutillon, W. Song, Y. Ren, C. Zhang, X. You, A. Burg,
> "Belief Propagation Decoding for Short Codes on Structured Sparse Parity-Check Matrices,"
> *IEEE International Symposium on Information Theory (ISIT)*, 2025.
> [IEEE Xplore](https://doi.org/10.1109/ISIT63088.2025.11195685)

and turned it into design rules for constructing new codes with sparse PCMs:

> Y. Shen, H. S. Ünal, A. Burg,
> "HDPC Codes with LDPC Matrices: Construction Based on Social Golfer Problem,"
> *IEEE International Symposium on Information Theory (ISIT)*, 2026.
> [IEEE Xplore](https://doi.org/10.1109/ISIT62367.2026.11653623)

I like this approach because it can turn a dense PCM into a sparser one suited to BP decoding.
As with short LDPC codes, the main challenge then becomes bringing BP closer to maximum-likelihood
decoding (MLD) with as little computational overhead as possible, which is another story.

Unfortunately, to the best of my knowledge, neither the optimized PCMs of these papers nor the code
to construct them have been made publicly available. This repository started as my attempt to fill
this gap, with a Julia implementation of the JSAC method and a BP simulator to measure the resulting
frame error rates. It now provides a more general tool that searches, for a given short code, the
PCM augmented with auxiliary variable nodes that gives the best FER under BP. It comes with the best
PCMs found for the codes of the JSAC paper, and the recipes to rebuild them.

This is my own implementation of the papers, plus a few ideas of my own, developed with the help of
Opus 5.5 and GPT-6.1-Sol. It is not the authors' code, so any error or misinterpretation of the
algorithms is mine.

**Current status (November 2026):** for every code of the JSAC paper, one tool builds a PCM whose FER
under normalized min-sum decoding is better than the published one: only marginally for some codes,
by a factor of about 2 for others, and 30 for RM(128,64). [Have a look at this plot](media/benchmark.png)
to compare them with the paper's results. The reproduction of the PCM-I and PCM-II matrices of the
JSAC paper is documented [here](docs/jsac_reproduction.md).

**Lessons learned:**

- Choosing the parity checks of the starting PCM deserves as much attention as optimizing the
  resulting graph.
- Graph metrics are useful guides when searching large design spaces, but they are only approximate
  proxies for FER: the ACE-based selection rule of the JSAC paper sometimes preferred matrices that
  decoded worse than others found along the way. The final choice is best made by a short
  simulation.
- Removing every 4-cycle is not always the best strategy: in a structured PCM, a few 4-cycles
  between blocks do much less harm than the dozens or hundreds of auxiliary VNs needed to remove
  them. What matters is to keep low the number of checks shared by any two code bits.
- Tuning the normalization factor α of min-sum decoding can lower the FER by more than an order of
  magnitude compared with α = 0.7, the value used in the paper, especially on structured PCMs.

## Installation

Requires Julia 1.12 or later.

```
git clone <this repository>
cd <repository>
julia --project -e 'using Pkg; Pkg.instantiate()'       # the PCM construction
julia --project=sim -e 'using Pkg; Pkg.instantiate()'   # the BP simulator
```

## Best PCMs

All are structured sparse PCMs (ssPCMs) built by `scripts/sspcm.jl`, and decoded with flooding NMS
and the number of iterations of the paper. FER is given at one Eb/N0 point, against the best
published curve, always that of the JSAC paper (the ISIT 2025 paper gives the same PCM for
BCH(63,45)).

| Code | Size | AVNs | Edges | 4-cycles | α | Iterations | FER here / JSAC |
|---|---|---:|---:|---:|---:|---:|---|
| RM(32,16) | 72×41 | 9 | 360 | 0 | 0.5 | 6 | 2.7e-4 / 4.5e-4 at 5 dB |
| BCH(63,45) | 63×72 | 9 | 567 | 0 | 0.5 | 20 | 1.6e-4 / 2.0e-4 at 5 dB |
| eBCH(64,24) | 336×85 | 21 | 1680 | 0 | 0.35 | 20 | 2.9e-4 / 5.2e-4 at 4 dB |
| RM(128,64) | 560×163 | 35 | 5040 | 7936 | 0.3 | 20 | 2.6e-4 / 8.7e-3 at 4 dB |
| Polar(128,20) + CRC24 | 288×150 | 22 | 1334 | 275 | 0.55 | 200 | 3.0e-5 / 7.9e-5 at 2 dB |
| Polar(128,40) + CRC24 | 308×154 | 26 | 2635 | 6493 | 0.45 | 200 | 2.6e-4 / 5.1e-4 at 3 dB |
| Polar(128,72) + CRC24 | 231×194 | 66 | 2773 | 14129 | 0.5 | 200 | 1.9e-4 / 2.2e-4 at 4 dB |

Each PCM is in `results/<code>.sspcm.mat`. To rebuild them:

```
julia --project -t auto scripts/sspcm.jl codes/rm.32.16 B=9
julia --project -t auto scripts/sspcm.jl codes/bch.63.45 B=9
julia --project -t auto scripts/sspcm.jl codes/ebch.64.24 avnfree=0 B=21
julia --project -t auto scripts/sspcm.jl codes/rm.128.64 tmin=8 avnfree=0 B=35
julia --project -t auto scripts/sspcm.jl codes/polar.128.20.crc24c wcap=5 B=30
julia --project -t auto scripts/sspcm.jl codes/polar.128.40.crc24c B=42
julia --project -t auto scripts/sspcm.jl codes/polar.128.72.crc24c wcap=23 T_max=30000 B=72
```

Each takes less than a minute on 96 threads. To simulate one, at its best α, then at the paper's:

```
cd sim
julia --project=. -t auto scripts/simulate.jl config/<code>.sspcm.yaml
julia --project=. -t auto scripts/simulate.jl config/<code>.sspcm.yaml decoder.norm_factor=0.7
```

A few notes:

- BCH(63,45) and eBCH(64,24): up to the order of rows and auxiliary VNs, the PCMs are those
  already found with the JSAC method (the 63×72 PCM-II, and the 336×85 PCM-II with η̄_min = 52, see
  the [reproduction notes](docs/jsac_reproduction.md)). Only α changes. In both, the blocks cover every
  pair of code bits at most once, and there are as many of them as this allows.
- RM(128,64): FER keeps decreasing slowly as blocks are added, at the cost of edges. With `B=17`,
  the PCM is 272×145 with 2448 edges and no 4-cycle, for about twice the FER.
- Polar codes: Eb/N0 counts the 24 CRC bits as information bits, as inferred for Fig. 8 of the
  paper (see the [polar notes](docs/polar_codes.md)). `wcap=23` adds heavier checks to the search
  space, which lets more blocks in before code bits start to share many checks.

## Building a PCM for another code

1. Build the sequence of ssPCMs, from the first full-rank one to an edge budget:

   ```
   julia --project -t auto scripts/sspcm.jl codes/<code> [max_edges=4000]
   ```

   Every matrix is saved to `sweep/<code>/B<blocks>.mat`, and its size, edges, 4-cycles and Γmax
   (the largest number of checks shared by two code bits) are printed.
2. Compare them by a short simulation, for a few values of α:

   ```
   cd sim
   julia --project=. -t auto scripts/fer_screen.jl config/<a similar code>.yaml ../sweep/<code> alpha=0.4,0.5,0.7 ebn0_min=3 ebn0_max=4 min_cw=50000
   ```

3. Rebuild the best one with `B=<blocks>`, which saves it to `results/<code>.sspcm.mat`, and refine
   α around the best value.

If the result is poor, try `avnfree=0` (it helped eBCH and RM codes), a few values of `seed`, and a
larger `wcap`. For a very large search space, such as the 94 488 checks of RM(128,64), use
`avnfree=0` and set `tmin` to half the weight of the lightest checks.

## Choosing α

α is a decoder setting: it does not change the PCM, so it is tuned after the construction, with
`fer_screen.jl`. On the codes here, ssPCMs need a smaller α than the JSAC PCMs, the more so as the
code bits get more checks: 0.5 for RM(32,16) and BCH(63,45) (8 or 9 checks per code bit), 0.35 for
eBCH(64,24) (21) and 0.3 for RM(128,64) (35). The optimum is flat within about ±0.05. A small α
slows down BP, so it pays to check the iteration budget at the same time.

The ranking of the PCMs of one sequence hardly depends on α, so `fer_screen.jl` can pick the number of
blocks at a fixed α. The ranking of different kinds of PCMs does depend on it: an ssPCM decoded with
α = 0.7 can lose to a JSAC PCM that it beats at its own best α.

## Contents

### Algorithms

| Paper | Function | File |
|---|---|---|
| JSAC, Algorithm 1: greedy 4-cycle removal | `remove_4cycles` | `src/four_cycle.jl` |
| JSAC, Algorithm 2: search space of low-weight dual codewords | `gen_search_space` | `src/search_space.jl` |
| JSAC, Algorithm 3: PCM extraction by row growth (PRG) | `extract_pcm` | `src/prg.jl` |
| JSAC, Algorithm 4: rapid bypass | `rapid_bypass` | `src/four_cycle.jl` |
| ISIT 2025: ssPCM from parity-check row blocks | `sspcm` | `src/sspcm.jl` |
| Check that a PCM defines exactly the code | `represents_code` | `src/sspcm.jl` |
| η̄_min, girth, cycle counts | `pcm_metrics` | `src/cycles.jl` |
| Trapping-set distribution 𝒯(a,b) | `trapping_distribution` | `src/trapping.jl` |

An ssPCM is made of parity-check row blocks. A block takes s low-weight dual codewords that share
the same t code bits and are disjoint elsewhere, and replaces those t bits by a single auxiliary VN:
s rows without the shared bits, plus one row that defines the auxiliary VN. No cycle lies inside a
block, and auxiliary VNs of different blocks share no check. `sspcm` adds blocks greedily. It
departs from the ISIT 2025 paper in two ways, without which it fails on polar codes: until the PCM
has full rank, only blocks that increase the rank are considered; and no block may let two code
bits share more than Γ checks, where Γ starts at 1 and is raised only when no block fits (Theorem 1
of the ISIT 2026 paper). Ties go to the largest block, then the one adding the fewest 4-cycles,
then the most even column weights, then at random.

η̄_min is the cycle-connectivity metric of the JSAC paper: the smallest ACE (the sum of d_v − 2 over
the VNs of the cycle) among the 6-cycles that go through an auxiliary VN. The larger, the better.
The scripts report −1 when no auxiliary VN lies on a 6-cycle.

### Scripts

All commands are run from the repository root. Options are given as `key=value`, and an unknown
option is an error.

The most expensive computations are multithreaded. Julia uses a single thread by default: use
`julia -t X` for X threads, or `julia -t auto` for one thread per logical core. For the codes of
this repository, the PCMs do not depend on the number of threads.

The three original CRC24C polar codes can be regenerated with
`julia --project scripts/build_polar_codes.jl`. Add `payload=20`, `payload=40` or `payload=72`
to generate only one of them.

**ssPCM:**

```
julia --project -t auto scripts/sspcm.jl codes/<code> [wcap=…] [T_max=8000] [seed_s=1] [tmin=2] [avnfree=1] [max_edges=4000] [seed=1] [B=…] [out=…]
```

- `wcap`, `T_max`, `seed_s`: search space of low-weight dual codewords, as for PCM-II below.
- `tmin`: smallest number of code bits shared by the checks of a block. Default: 2.
- `avnfree=0`: forbid blocks without an auxiliary VN (rows that pairwise share at most one bit).
- `max_edges`: stop the sequence once the PCM has more edges. Default: 4000.
- `seed`: seed of the random tie-breaks. Default: 1.
- `B`: stop at B blocks and save that PCM to `results/<code>.sspcm.mat` (or `out`). Without it,
  the whole sequence is saved to `sweep/<code>/` (or the directory `out`).

**PCM-I (JSAC Algorithm 1 on the code's PCM):**

Removes the 4-cycles of the code's PCM and saves the resulting girth-6 PCM to `results/<code>.pcm1.mat`.

```
julia --project scripts/pcm_i.jl codes/<code> [seed=X | n_seed=M] [out=file.mat]
```

At each step, Algorithm 1 picks the all-ones block of the PCM whose removal eliminates the most
4-cycles. Without options, ties between blocks are broken deterministically (largest block, then
last one found). `seed=X` breaks them at random with seed `X`;
`n_seed=M` tries seeds 1 to M and keeps the PCM with the fewest auxiliary VNs, then the largest
η̄_min, then the fewest 6-cycles, then the fewest edges.

**PCM-II (JSAC Algorithms 2 and 3):**

Builds the search space of low-weight dual codewords, extracts from it the girth-6 PCM with the
largest η̄_min, and saves it to `results/<code>.pcm2.mat`.

```
julia --project -t auto scripts/pcm_ii.jl codes/<code> [M_a=…] [seed_p=…] [seed_s=…] [T_max=…] [wcap=…] [rapid=1] [out=file.mat]
```

- `M_a`: maximum number of rows of the auxiliary matrix `H_a` (the x-axis of Fig. 5 in the paper).
  Default: the number of codewords in the search space.
- `seed_p`: seed of the random tie-breaks of Algorithm 3. Default: 6.
- `seed_s`: seed of the random column permutations of Algorithm 2. Default: 1.
- `T_max`: number of random trials of Algorithm 2. Default: 8000.
- `wcap`: initial weight cap of the codewords collected by Algorithm 2. Default: the weight of the
  lightest row of the code's PCM. The cap is raised automatically until the codewords span the dual
  code. A larger value adds heavier codewords to the search space. For some codes the default
  collects far too many codewords: the polar (128,20) code needs `wcap=5`.
- `rapid=1`: use Algorithm 4.

**Trapping-set distribution of a PCM:**

```
julia --project scripts/trapping_report.jl <file.mat> [a_max=4]
```

**PCM performance under BP over the BI-AWGN channel:**

```
cd sim
julia --project=. -t auto scripts/simulate.jl config/<config>.yaml [key=value …]
julia --project=. -t auto scripts/fer_screen.jl config/<config>.yaml <file.mat | dir>… [alpha=…] [key=value …]
```

`simulate.jl` estimates the frame and bit error rates over a range of Eb/N0. Any setting can be
overridden on the command line (e.g. `decoder.max_iter=50`), and `use_base_pcm=true` decodes with
the original PCM instead, for comparison. `fer_screen.jl` runs short simulations of several PCMs, for
each value of α, and prints one line of FER per run. Results are written to `sim/logs/`. See
[`sim/README.md`](sim/README.md) for the available decoders, and the
[simulation records](sim/logs/README.md) for the reference runs and their replay settings.
`uv run scripts/plot_benchmark.py` redraws [the benchmark plot](media/benchmark.png) from the
reference records.

### Files

Codes are stored in `codes/<code>.mat` (`H`, `G`, `n`, `k`, optional `dmin`, and `Ginv` for
non-systematic encoders). An optimized PCM is saved in the same format, with the optimized PCM in
`H`, the original one in `H_base`, the indices of the auxiliary VNs in `puncturing_set`, and its
`eta_min` when it has no 4-cycle. The best PCMs are `results/<code>.sspcm.mat`; the other files in
`results/` come from the [reproduction of the JSAC paper](docs/jsac_reproduction.md).

The test suite (`julia --project test/runtests.jl`, about 1 minute) checks the JSAC algorithms on
RM(8,4), RM(32,16), BCH(63,45) and eBCH(64,24), and the ssPCMs of RM(32,16) and BCH(63,45).

## License

MIT (see [LICENSE](LICENSE)). This code is provided for educational and research purposes.
