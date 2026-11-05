# Reference simulations

See this [plot](../../media/benchmark.png) to compare the best PCMs with the paper's.

Each supplied optimized PCM except RM(8,4) has a CSV performance record and a
YAML sidecar recording the decoder settings, code parameters, random seed,
batch size, thread count and Julia version. The archived files use
`.reference` filenames; new simulations write timestamped files.

## Replay

From `sim/`, run:

```
julia --project=. -t 90 scripts/simulate.jl config/<code>.sspcm.yaml                       # best α
julia --project=. -t 90 scripts/simulate.jl config/<code>.sspcm.yaml decoder.norm_factor=0.7
julia --project=. -t 16 scripts/simulate.jl config/<PCM>.yaml                              # JSAC PCMs
```

The ssPCM records were made with 90 threads, the others with 16. Use the same count for an exact
replay. The [main README](../../README.md) lists the command for each code.

### Best PCMs (ssPCM)

Each has two records: `<code>.sspcm` at the best α of its configuration, and
`<code>.sspcm.alpha0.7` at α = 0.7.

| PCM | Configuration | Best α | CSV | Eb/N0 range, dB |
|---|---|---:|---|---|
| RM(32,16) | [YAML](../config/rm.32.16.sspcm.yaml) | 0.5 | [best α](rm.32.16.sspcm.reference.csv), [0.7](rm.32.16.sspcm.alpha0.7.reference.csv) | 0–7 |
| BCH(63,45) | [YAML](../config/bch.63.45.sspcm.yaml) | 0.5 | [best α](bch.63.45.sspcm.reference.csv), [0.7](bch.63.45.sspcm.alpha0.7.reference.csv) | 1–6 |
| eBCH(64,24) | [YAML](../config/ebch.64.24.sspcm.yaml) | 0.35 | [best α](ebch.64.24.sspcm.reference.csv), [0.7](ebch.64.24.sspcm.alpha0.7.reference.csv) | 0–5.5 |
| RM(128,64) | [YAML](../config/rm.128.64.sspcm.yaml) | 0.3 | [best α](rm.128.64.sspcm.reference.csv), [0.7](rm.128.64.sspcm.alpha0.7.reference.csv) | 0–5.5 |
| Polar(128,20) | [YAML](../config/polar.128.20.crc24c.sspcm.yaml) | 0.55 | [best α](polar.128.20.crc24c.sspcm.reference.csv), [0.7](polar.128.20.crc24c.sspcm.alpha0.7.reference.csv) | 0–2.5, CRC-inclusive |
| Polar(128,40) | [YAML](../config/polar.128.40.crc24c.sspcm.yaml) | 0.45 | [best α](polar.128.40.crc24c.sspcm.reference.csv), [0.7](polar.128.40.crc24c.sspcm.alpha0.7.reference.csv) | 0–3.5, CRC-inclusive |
| Polar(128,72) | [YAML](../config/polar.128.72.crc24c.sspcm.yaml) | 0.5 | [best α](polar.128.72.crc24c.sspcm.reference.csv), [0.7](polar.128.72.crc24c.sspcm.alpha0.7.reference.csv) | 1–4, CRC-inclusive |

### JSAC reproduction

| PCM | Configuration | CSV | Eb/N0 range, dB |
|---|---|---|---|
| RM(32,16), PCM-I | [YAML](../config/rm.32.16.pcm1.yaml) | [CSV](rm.32.16.pcm1.reference.csv) | 0–7 |
| RM(32,16), PCM-II | [YAML](../config/rm.32.16.pcm2.yaml) | [CSV](rm.32.16.pcm2.reference.csv) | 0–7 |
| BCH(63,45), PCM-II | [YAML](../config/bch.63.45.pcm2.yaml) | [CSV](bch.63.45.pcm2.reference.csv) | 1–6.5 |
| eBCH(64,24), PCM-II | [YAML](../config/ebch.64.24.pcm2.yaml) | [CSV](ebch.64.24.pcm2.reference.csv) | 0–6 |
| eBCH(64,24), η̄_min = 52 | [YAML](../config/ebch.64.24.pcm2.eta52.yaml) | [CSV](ebch.64.24.pcm2.eta52.reference.csv) | 0–6 |
| RM(128,64), PCM-II | [YAML](../config/rm.128.64.pcm2.yaml) | [CSV](rm.128.64.pcm2.reference.csv) | 0–6 |
| RM(128,64), geometric PCM-I | [YAML](../config/rm.128.64.geom.pcm1.yaml) | [CSV](rm.128.64.geom.pcm1.reference.csv) | 0–6 |
| Polar(128,20), PCM-II candidate | [YAML](../config/polar.128.20.crc24c.pcm2.yaml) | [CSV](polar.128.20.crc24c.pcm2.reference.csv) | 0–3, CRC-inclusive |
| Polar(128,40), PCM-II candidate | [YAML](../config/polar.128.40.crc24c.pcm2.yaml) | [CSV](polar.128.40.crc24c.pcm2.reference.csv) | 0–4, CRC-inclusive |
| Polar(128,72), PCM-II candidate | [YAML](../config/polar.128.72.crc24c.pcm2.yaml) | [CSV](polar.128.72.crc24c.pcm2.reference.csv) | 1–4, CRC-inclusive |

All grids use 0.5 dB spacing. The upper endpoints are chosen to retain at least
50 observed frame errors, rather than extending the curves with poorly
sampled or zero-error points.

The runs use random codewords, flooding NMS (factor 0.7 for the JSAC PCMs), Float32
messages and LLR saturation at 128. RM(32,16) uses six iterations, as in Fig. 4; the
other non-polar codes use 20, as in Fig. 9; the polar codes use 200, as in
Fig. 8. The BCH record includes BER, but uses 20 iterations rather than the
five iterations of Fig. 7.

## Stopping and precision

Each point runs until both 200000 frames and 50 frame errors have been reached.
There is no frame limit. With 16 threads and batches of 250 frames per thread,
the stopping rule is checked every 4000 frames (22500 with 90 threads), so the recorded count can
exceed either target. Every retained point, including the final one, therefore
contains at least 50 errors. At 50 errors, the approximate relative sampling
standard error of a small FER is 14%.

## Polar Eb/N0 convention

The simulator uses payload-normalized Eb/N0. Comparison with Fig. 8 uses the
inferred convention that includes the 24 CRC bits in the information-bit
count. The polar configurations shift the simulated Eb/N0 grid by
`10 log10((k+24)/k)`. Both the YAML settings and CSV `Eb/N0` column contain
the shifted, payload-normalized values. To plot them on Fig. 8's inferred
scale, subtract the shift (the CRC-inclusive ranges are in the tables above):

| Payload k | Shift, dB |
|---|---:|
| 20 | 3.424226808222 |
| 40 | 2.041199826559 |
| 72 | 1.249387366083 |

This changes only the noise level associated with an Eb/N0 label, not the
encoder dimension. The paper does not explicitly specify this convention;
see the [polar construction notes](../../docs/polar_codes.md).
