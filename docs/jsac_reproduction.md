# Reproduction of Shen et al., JSAC 2025

This page records how the results of the paper were reproduced with the JSAC algorithms
(`scripts/pcm_i.jl` and `scripts/pcm_ii.jl`): what matches, what does not, and the commands to
regenerate each PCM. The best PCMs found since then are listed in the [README](../README.md).

## RM(8,4): example of eq. (16)

```
julia --project scripts/pcm_ii.jl codes/rm.8.4 M_a=14 seed_p=1
```

Reproduces the 10×10 PCM of eq. (16): 32 edges, 2 auxiliary VNs, η̄_min = 4.
Result: `results/rm.8.4.pcm2.mat`.

## RM(32,16): Table I and Fig. 4

PCM-II. Reproduces Table I exactly: search space of 620 codewords of weight 8, PCM-II of size 48×38
with 240 edges, 6 auxiliary VNs, 2336 6-cycles, η̄_min = 14, and the same trapping-set distribution.
Result: `results/rm.32.16.pcm2.mat`.

```
julia --project scripts/pcm_ii.jl codes/rm.32.16 M_a=48 seed_p=6
(cd sim && julia --project=. -t 16 scripts/simulate.jl config/rm.32.16.pcm2.yaml)
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
(cd sim && julia --project=. -t 16 scripts/simulate.jl config/rm.32.16.pcm1.yaml)
```

| PCM-I | FER at 5 dB | 6 dB | 7 dB |
|---|---|---|---|
| paper (Fig. 4) | 2.3e-2 | 4.0e-3 | 3.5e-4 |
| 32×48, 16 auxiliary VNs (seed 1815) | 5.2e-2 | 1.1e-2 | 1.2e-3 |
| 28×44, 12 auxiliary VNs | 2.9e-2 | 4.6e-3 | 3.9e-4 |

Trapping sets. The 𝒯(a,b) distributions of Table I are obtained with:

```
julia --project scripts/trapping_report.jl results/rm.32.16.pcm2.mat                  # PCM-II
julia --project scripts/trapping_report.jl results/rm.32.16.pcm1.seed1815.mat         # PCM-I (seed 1815)
julia --project scripts/trapping_report.jl codes/rm.32.16.mat                         # original PCM
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

## BCH(63,45)

Reproduces the PCM-II given in the text: search space of 189 codewords of weight 16, PCM-II of size
63×72 with 567 edges, row weight 9, η̄_min = 17. Result: `results/bch.63.45.pcm2.mat`.

```
julia --project scripts/pcm_ii.jl codes/bch.63.45 M_a=70 seed_p=1
(cd sim && julia --project=. -t 16 scripts/simulate.jl config/bch.63.45.pcm2.yaml)
```

Larger values of `M_a` give the same PCM: Algorithm 3 stops by itself at 54 rows.

## eBCH(64,24)

Reproduces the PCM-II of the paper: 272×81, 1360 edges, 17 auxiliary VNs, η̄_min = 44 (about 30 min
on 48 threads). Result: `results/ebch.64.24.pcm2.mat`.

```
julia --project -t auto scripts/pcm_ii.jl codes/ebch.64.24 M_a=255 seed_p=1
(cd sim && julia --project=. -t 16 scripts/simulate.jl config/ebch.64.24.pcm2.yaml)
```

The search space contains 2520 codewords of weight 8, not 2560 as written in the paper. I checked
this number by enumerating all 2²⁴ codewords and applying the MacWilliams identity.

Improvement. η̄_min = 44 is not the maximum. Past 255 rows, η̄_min keeps increasing by 2 every 15
rows, with drops in between (the saw-tooth curve of Fig. 5). With `M_a=345`, the best PCM is found at
315 rows: 336×85, 1680 edges, 21 auxiliary VNs, η̄_min = 52 (about 1.5 h on 48 threads).
Result: `results/ebch.64.24.pcm2.eta52.mat`.

```
julia --project -t auto scripts/pcm_ii.jl codes/ebch.64.24 M_a=345 seed_p=1 out=results/ebch.64.24.pcm2.eta52.mat
(cd sim && julia --project=. -t 16 scripts/simulate.jl config/ebch.64.24.pcm2.eta52.yaml)
```

It gains about 0.1 dB over the paper's PCM, for 24% more edges (NMS, 20 iterations):

| PCM-II | η̄_min | FER at 2 dB | 3 dB | 4 dB | 5 dB |
|---|---|---|---|---|---|
| 272×81 (reproduction) | 44 | 5.60e-2 | 8.47e-3 | 6.05e-4 | 1.06e-5 |
| 336×85 | 52 | 4.83e-2 | 6.82e-3 | 4.30e-4 | 8.15e-6 |

The row count at which η̄_min stops increasing depends on `seed_p`, so other seeds with a large `M_a`
may do better.

## RM(128,64)

The search space matches the paper (94 488 codewords of weight 16), but I could not obtain the
141×141 PCM-II with 13 auxiliary VNs of the paper. The best PCM-II I found is 147×161 with 1078
edges, 33 auxiliary VNs and η̄_min = 9 (about 4.5 h on 96 threads). Other values of `seed_p` gave the same
number of auxiliary VNs. Result: `results/rm.128.64.pcm2.mat`.

```
julia --project -t auto scripts/pcm_ii.jl codes/rm.128.64 M_a=150 seed_p=6
(cd sim && julia --project=. -t 16 scripts/simulate.jl config/rm.128.64.pcm2.yaml)
```

However, GPT-6.1-Sol found an optimized PCM with *fewer* auxiliary VNs and *better FER* than the one reported in the paper by
starting from a different, redundant PCM constructed using finite geometry to have many overlapping checks.
The new base PCM, `codes/rm.128.64.geom.mat`, has 135 weight-16 checks.
Passing it to Algorithm 1 gives a 144×137 PCM with 9 auxiliary VNs, 1296 edges and η̄_min = 28:

```
julia --project scripts/pcm_i.jl codes/rm.128.64.geom
(cd sim && julia --project=. -t 16 scripts/simulate.jl config/rm.128.64.geom.pcm1.yaml)
```

See the [construction details and recipe](rm_128_64.md).

Both candidates have FER close to the paper's (NMS, 20 iterations), with the PCM-I obtained from the geometric base
matrix slightly outperforming the paper's PCM-II:

| PCM | FER at 3 dB | 4 dB |
|---|---|---|
| 141×141 (paper, Fig. 9) | ≈ 1.2e-1 | ≈ 8e-3 |
| 147×161 | 1.31e-1 | 1.03e-2 |
| 144×137 (geometric) | 1.13e-1 | 7.54e-3 |

## Polar codes: partially reproduced

The original PCMs for the three length-128 polar codes, using the 5G NR downlink (DCI)
construction, can be regenerated with:

```
julia --project scripts/build_polar_codes.jl
```

The files are saved as `codes/polar.128.<k>.crc24c.mat`. Add `payload=20`, `payload=40` or
`payload=72` to generate only one code. Each code uses CRC24C, the NR input interleaver and a
length-128 polar transform, before rate matching. The dimension counts payload bits only;
the unfrozen positions also include the 24 CRC bits.

The reliability and input-interleaving tables were taken from NVIDIA's
[Sionna](https://github.com/nvlabs/sionna) library.

| Payload | Unfrozen bits | Frozen checks + CRC checks | Original PCM | Edges |
|---|---|---|---|---|
| 20 | 44 | 84 + 24 | 108×128 | 2401 |
| 40 | 64 | 64 + 24 | 88×128 | 2481 |
| 72 | 96 | 32 + 24 | 56×128 | 2425 |

The PCM dimensions and edge counts match those reported in Fig. 8 of the paper.

The optimized PCMs are only partially reproduced. The candidates below were built using low-weight
checks from Algorithm 2 and row growth from Algorithm 3, with Algorithm 1 removing 4-cycles.
They are the optimized matrices obtained after selecting 236, 174 and 145 checks, respectively,
into the redundant input PCM, rather than the final matrices returned by Algorithm 3's
ACE-based selection rule.

| Code | PCM-II | Size | Auxiliary VNs | Edges |
|---|---|---|---|---|
| (128,20) | Paper | 323×215 | 87 | 1224 |
| (128,20) | Here | 323×215 | 87 | 1225 |
| (128,40) | Paper | 329×283 | 155 | 1517 |
| (128,40) | Here | 329×283 | 155 | 1538 |
| (128,72) | Paper | 409×420 | 292 | 1972 |
| (128,72) | Here | 472×458 | 330 | 2197 |

Results: [20-bit candidate](../results/polar.128.20.crc24c.pcm2.mat),
[40-bit candidate](../results/polar.128.40.crc24c.pcm2.mat) and
[72-bit candidate](../results/polar.128.72.crc24c.pcm2.mat).

```
(cd sim && julia --project=. -t 16 scripts/simulate.jl config/polar.128.20.crc24c.pcm2.yaml)
(cd sim && julia --project=. -t 16 scripts/simulate.jl config/polar.128.40.crc24c.pcm2.yaml)
(cd sim && julia --project=. -t 16 scripts/simulate.jl config/polar.128.72.crc24c.pcm2.yaml)
```

The first two have FER close to the paper's; the third still has about three times the paper's
FER at 4 dB. These comparisons use an inferred energy normalization that includes CRC bits.
The exact published PCMs have not been reproduced. See the [construction notes](polar_codes.md)
for the selected inputs, replay recipe, graph statistics and FER measurements.
