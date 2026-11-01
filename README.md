# PCM optimization for BP decoding of short block codes

Belief propagation (BP) is excellent for decoding long LDPC codes, but performs poorly on short
non-LDPC codes such as RM, BCH, eBCH or polar codes. Their standard parity-check matrices (PCMs) are
dense, and their Tanner graphs are riddled with 4-cycles and small trapping sets. In 2025, Yifei Shen
and co-authors proposed a general method to extract from the PCM of a short binary block code an
optimized sparse PCM with enhanced short-cycle connectivity, on which BP decoding performs
remarkably well:

> Y. Shen, Z. Li, Y. Ren, E. Boutillon, A. Balatsoukas-Stimming, C. Zhang, X. You, A. Burg,
> "Toward Universal Belief Propagation Decoding for Short Binary Block Codes,"
> *IEEE Journal on Selected Areas in Communications*, vol. 43, no. 4, pp. 1135–1149, April 2025.
> [IEEE Xplore](https://ieeexplore.ieee.org/document/10887536)

Unfortunately, and to the best of my knowledge, neither the optimized PCMs of the paper nor the code
to construct them have been made publicly available. This repository is my attempt to fill this
gap. It is a Julia implementation of the method, which builds sparse PCMs of girth 6, often redundant,
possibly with auxiliary (punctured) variable nodes, together with a BP simulator to measure the
resulting frame error rates. Its main purpose is to reproduce the results of the paper.

This is my own implementation of the paper, written with the help of Claude Code. It is not the
authors' code, so any error or misinterpretation of the algorithms is mine.

## Installation

Requires Julia 1.12 or later.

```
git clone <this repository>
cd <repository>
julia --project -e 'using Pkg; Pkg.instantiate()'       # the paper's algorithms
julia --project=sim -e 'using Pkg; Pkg.instantiate()'   # the BP simulator
```

## Contents

### Algorithms

| Paper | Function | File |
|---|---|---|
| Algorithm 1: greedy 4-cycle removal | `remove_4cycles` | `src/four_cycle.jl` |
| Algorithm 2: search space of low-weight dual codewords | `gen_search_space` | `src/search_space.jl` |
| Algorithm 3: PCM extraction by row growth (PRG) | `extract_pcm` | `src/prg.jl` |
| Algorithm 4: rapid bypass | `rapid_bypass` | `src/four_cycle.jl` |
| η̄_min, girth, cycle counts | `pcm_metrics` | `src/cycles.jl` |
| Trapping-set distribution 𝒯(a,b) | `trapping_distribution` | `src/trapping.jl` |

η̄_min is the cycle-connectivity metric of the paper: the smallest ACE (the sum of d_v − 2 over the
VNs of the cycle) among the 6-cycles that go through an auxiliary VN. The larger, the better. The
scripts report −1 when no auxiliary VN lies on a 6-cycle.

### Scripts

All commands are run from the repository root. Options are given as `key=value`, and an unknown
option is an error.

The most expensive computations are multithreaded: the random trials of Algorithm 2, the evaluation
of the candidate rows in Algorithm 3, the seed search of `pcm_i.jl` (`n_seed=M`), and the BP
simulations. Julia uses a single thread by default. Use `julia -t X` to run calculations on X
threads, and `julia -t auto` to use one thread per logical core of the machine. For the codes of
this repository, the extracted PCMs do not depend on the number of threads.

**PCM-I (Algorithm 1 on the code's PCM):**

Removes the 4-cycles of the code's PCM and saves the resulting girth-6 PCM to `results/<code>.pcm1.mat`.

```
julia --project scripts/pcm_i.jl codes/<code> [seed=X | n_seed=M] [out=file.mat]
```

At each step, Algorithm 1 picks the all-ones block of the PCM whose removal eliminates the most
4-cycles. Without options, ties between blocks are broken deterministically (largest block, then
last one found). `seed=X` breaks them at random with seed `X`;
`n_seed=M` tries seeds 1 to M and keeps the PCM with the fewest auxiliary VNs, then the largest
η̄_min, then the fewest 6-cycles, then the fewest edges.

**PCM-II (Algorithms 2 and 3):**

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
  code. A larger value adds heavier codewords to the search space.
- `rapid=1`: use Algorithm 4.

**Trapping-set distribution of a PCM:**

Calculates the trapping-set distribution 𝒯(a,b) up to `a = a_max` for the specified PCM, and its
smallest trapping set.

```
julia --project scripts/trapping.jl <file.mat> [a_max=4]
```

**PCM performance under BP over the BI-AWGN channel:**

Estimates the frame and bit error rates of BP decoding with the specified PCM, over a range of Eb/N0.

```
cd sim
julia --project=. -t auto scripts/simulate.jl config/<config>.yaml [key=value …]
```

The configs in `sim/config/` use normalized min-sum (α = 0.7) with flooding, as in the paper. Any
setting can be overridden on the command line (e.g. `decoder.max_iter=50`), and `use_base_pcm=true`
decodes with the original PCM instead, for comparison. Results are written to `sim/logs/`. See
[`sim/README.md`](sim/README.md) for the available decoders.

### Files

Codes are stored in `codes/<code>.mat` (`H`, `G`, `n`, `k`, `dmin`, and `Ginv` for non-systematic
encoders). An optimized PCM is saved in the same format, with the optimized PCM in `H`, the original
one in `H_base`, the indices of the auxiliary VNs in `puncturing_set`, and its `eta_min`. The
optimized PCMs of the next section are in `results/`.

The test suite (`julia --project test/runtests.jl`, about 1 minute) checks the RM(8,4) and RM(32,16)
PCM-II of the next section, the size of the BCH(63,45) PCM-II, the search spaces of RM(8,4),
RM(32,16), BCH(63,45) and eBCH(64,24), and the trapping sets of the RM(32,16) PCM-II.

## Reproduction of the paper's results

### RM(8,4): example of eq. (16)

```
julia --project scripts/pcm_ii.jl codes/rm.8.4 M_a=14 seed_p=1
```

Reproduces the 10×10 PCM of eq. (16): 32 edges, 2 auxiliary VNs, η̄_min = 4.
Result: `results/rm.8.4.pcm2.mat`.

### RM(32,16): Table I and Fig. 4

PCM-II. Reproduces Table I exactly: search space of 620 codewords of weight 8, PCM-II of size 48×38
with 240 edges, 6 auxiliary VNs, 2336 6-cycles, η̄_min = 14, and the same trapping-set distribution.
Result: `results/rm.32.16.pcm2.mat`.

```
julia --project scripts/pcm_ii.jl codes/rm.32.16 M_a=48 seed_p=6
(cd sim && julia --project=. -t auto scripts/simulate.jl config/rm.32.16.pcm2.yaml)
```

PCM-I. Algorithm 1 with the default tie-break gives a 32×48 PCM with 16 auxiliary VNs, 122 edges and
η̄_min = 2, but 39 6-cycles instead of 38. Random tie-breaks with `seed=1815` give a PCM-I with the
size, number of edges, number of 6-cycles and η̄_min of Table I:

```
julia --project scripts/pcm_i.jl codes/rm.32.16 seed=1815 out=results/rm.32.16.pcm1.seed1815.mat
```

However, this PCM decodes 2 to 3 times worse than the PCM-I curve of Fig. 4 (NMS, 6 iterations), and so
do all the other PCM-Is with 16 auxiliary VNs I tried. I have not found the cause. Searching over 20 000
seeds finds a smaller PCM-I, 28×44 with 12 auxiliary VNs, whose FER matches Fig. 4.
Result: `results/rm.32.16.pcm1.mat`.

```
julia --project scripts/pcm_i.jl codes/rm.32.16 seed=11474
(cd sim && julia --project=. -t auto scripts/simulate.jl config/rm.32.16.pcm1.yaml)
```

| PCM-I | FER at 5 dB | 6 dB | 7 dB |
|---|---|---|---|
| paper (Fig. 4) | 2.3e-2 | 4.0e-3 | 3.5e-4 |
| 32×48, 16 auxiliary VNs (seed 1815) | 5.2e-2 | 1.1e-2 | 1.2e-3 |
| 28×44, 12 auxiliary VNs | 2.9e-2 | 4.6e-3 | 3.9e-4 |

Trapping sets. The 𝒯(a,b) distributions of Table I are obtained with:

```
julia --project scripts/trapping.jl results/rm.32.16.pcm2.mat                  # PCM-II
julia --project scripts/trapping.jl results/rm.32.16.pcm1.seed1815.mat         # PCM-I (seed 1815)
julia --project scripts/trapping.jl codes/rm.32.16.mat                         # original PCM
```

The script counts, by exhaustive enumeration, every trapping set that meets the definition of
Karimi and Banihashemi (ref. [50] of the paper). For PCM-II, I obtain the distribution of Table I
exactly. For PCM-I, I do not. In a girth-6 PCM, the trapping sets with 3 VNs are exactly the
6-cycles, so the 𝒯(3,·) counts must add up to the number of 6-cycles. Mine do (38 for PCM-I, 2336
for PCM-II). But the PCM-I column of Table I lists 3 + 120 + 473 + 1208 + 2447 = 4251 trapping sets
with 3 VNs, for a PCM-I with 38 6-cycles. So either there is an error in Table I, or the paper
counted the trapping sets differently. I do not obtain the distribution of the original PCM either.
It has girth 4, and the search method of [50], which the paper presumably used, is not exhaustive
on such dense matrices.

### BCH(63,45)

Reproduces the PCM-II given in the text: search space of 189 codewords of weight 16, PCM-II of size
63×72 with 567 edges, row weight 9, η̄_min = 17. Result: `results/bch.63.45.pcm2.mat`.

```
julia --project scripts/pcm_ii.jl codes/bch.63.45 M_a=70 seed_p=1
(cd sim && julia --project=. -t auto scripts/simulate.jl config/bch.63.45.pcm2.yaml)
```

Larger values of `M_a` give the same PCM: Algorithm 3 stops by itself at 54 rows.

### eBCH(64,24)

Reproduces the PCM-II of the paper: 272×81, 1360 edges, 17 auxiliary VNs, η̄_min = 44 (about 30 min
on 48 threads). Result: `results/ebch.64.24.pcm2.mat`.

```
julia --project -t auto scripts/pcm_ii.jl codes/ebch.64.24 M_a=255 seed_p=1
(cd sim && julia --project=. -t auto scripts/simulate.jl config/ebch.64.24.pcm2.yaml)
```

The search space contains 2520 codewords of weight 8, not 2560 as written in the paper. I checked
this number by enumerating all 2²⁴ codewords and applying the MacWilliams identity.

Improvement. η̄_min = 44 is not the maximum. Past 255 rows, η̄_min keeps increasing by 2 every 15
rows, with drops in between (the saw-tooth curve of Fig. 5). With `M_a=345`, the best PCM is found at
315 rows: 336×85, 1680 edges, 21 auxiliary VNs, η̄_min = 52 (about 1.5 h on 48 threads).
Result: `results/ebch.64.24.pcm2.eta52.mat`.

```
julia --project -t auto scripts/pcm_ii.jl codes/ebch.64.24 M_a=345 seed_p=1 out=results/ebch.64.24.pcm2.eta52.mat
(cd sim && julia --project=. -t auto scripts/simulate.jl config/ebch.64.24.pcm2.eta52.yaml)
```

It gains about 0.1 dB over the paper's PCM, for 24% more edges (NMS, 20 iterations):

| PCM-II | η̄_min | FER at 2 dB | 3 dB | 4 dB | 5 dB |
|---|---|---|---|---|---|
| 272×81 (paper's PCM) | 44 | 5.60e-2 | 8.28e-3 | 6.04e-4 | 1.19e-5 |
| 336×85 | 52 | 4.88e-2 | 6.77e-3 | 4.00e-4 | 8.68e-6 |

The row count at which η̄_min stops increasing depends on `seed_p`, so other seeds with a large `M_a`
may do better.

### RM(128,64): not reproduced

The search space matches the paper (94 488 codewords of weight 16), but I could not obtain the
141×141 PCM-II with 13 auxiliary VNs of the paper. The best PCM-II I found is 147×161 with 1078
edges, 33 auxiliary VNs and η̄_min = 9 (about 4.5 h on 96 threads). Other values of `seed_p` gave the same
number of auxiliary VNs. Result: `results/rm.128.64.pcm2.mat`.

```
julia --project -t auto scripts/pcm_ii.jl codes/rm.128.64 M_a=150 seed_p=6
(cd sim && julia --project=. -t auto scripts/simulate.jl config/rm.128.64.pcm2.yaml)
```

Its FER is nevertheless close to the paper's (NMS, 20 iterations):

| PCM-II | FER at 3 dB | 4 dB |
|---|---|---|
| 141×141 (paper, Fig. 9) | ≈ 1.2e-1 | ≈ 8e-3 |
| 147×161 | 1.29e-1 | 1.03e-2 |

If you have been able to reproduce the paper's PCM-II for this code, I would love to hear from you:
please open an issue or contact me directly.

### Not covered yet

I have not run the (128, 20), (128, 40) and (128, 72) polar codes of the paper yet.

## License

MIT (see [LICENSE](LICENSE)). This code is provided for educational and research purposes.
