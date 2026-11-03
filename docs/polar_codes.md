# Polar code reproduction

The original length-128 CRC24C polar codes match the standard PCM dimensions
and edge counts in Fig. 8 of the paper. The optimized matrices reproduce only
part of the reported results: the 20- and 40-bit candidates follow the plotted
FER closely, while the 72-bit candidate remains worse. None is the exact
published PCM-II.

## Supplied matrices

All paths below are relative to the repository root. The optimized files use
the existing format: `H` is the decoding PCM, `H_base` is the original PCM,
and `puncturing_set` lists the trailing auxiliary columns. `G`, `Ginv`, `n`
and `k` are retained from the original code. Minimum distances are not certified,
so `dmin` is omitted.

| Payload | Optimized file in `results/` | Size | Edges | Auxiliary VNs | η̄_min |
|---|---|---|---:|---:|---:|
| 20 | `polar.128.20.crc24c.pcm2.mat` | 323×215 | 1225 | 87 | 8 |
| 40 | `polar.128.40.crc24c.pcm2.mat` | 329×283 | 1538 | 155 | 5 |
| 72 | `polar.128.72.crc24c.pcm2.mat` | 482×465 | 2232 | 337 | 3 |

All three have girth 6. Their auxiliary-related six-cycle counts are 179,
1313 and 1553, respectively, versus 193, 1105 and 1429 reported for the paper's
PCMs. The last paper count should be read alongside its different matrix size.

Each optimized file has a companion ending in `.input.mat` instead of `.mat`.
Its `H` is the selected redundant PCM on the original 128 coordinates, before
Algorithm 1. It retains the original encoder and includes the number of selected
checks (the `checkpoint` field), seed, first-row index, removal policy and selected search-space indices.
These files preserve the exact selected rows and their order, independently
of the randomized search that found them.

## Construction

Algorithm 2 searches for low-weight words of the dual code. Algorithm 3 grows
a redundant matrix from those checks, applying Algorithm 1 to remove 4-cycles
when evaluating and accepting each candidate. The retained graphs use full
Algorithm 1, without Algorithm 4's rapid bypass.

| Payload | Search rows | Weights | Initial row | PRG seed | Selected checks | Algorithm 1 tie rule |
|---|---:|---|---:|---:|---:|---|
| 20 | 266 | 4–5 | 1 | 44 | 236 | Existing `(δ,s,t)` preference |
| 40 | 619 | 8–10 | 80 | 1 | 174 | Maximum δ, last row pair wins |
| 72 | 474 | 14–21 | 13 | 1 | 145 | Maximum δ, last row pair wins |

The first-row indices refer to the stored search-space ordering, not coordinate
indices. The 40-bit initial check has support `{7,17,39,49,71,81,103,113}`
in one-based coordinates. Both tie policies use stable row deletion and the
existing degree-one/two cleanup rules.

For payload 20, exhaustive enumeration certifies 18 weight-4 and 248 weight-5
dual words, together spanning rank 108. For payload 40, it certifies 310
weight-8, 12 weight-9 and 297 weight-10 dual words, together spanning rank 88;
there are no lower-weight nonzero dual words. The corresponding randomized
searches recovered those complete sets.

For payload 72, the 474-row search space is a sample, not a completeness
certificate. It contains 19 weight-14, 152 weight-16, one weight-17, 134
weight-18, six weight-19, 25 weight-20 and 137 weight-21 checks, spanning rank
56. It was obtained by removing the single weight-two check from a stable
475-row sample. Thus this search-space choice is a further departure from
using every discovered check through the first spanning weight.

### Selection of the retained matrices

Algorithm 3 adds checks one at a time to a redundant input PCM, applying
Algorithm 1 to obtain an optimized PCM at each step. We saved the optimized
matrices obtained when the input contained 236, 174 and 145 checks for the
20-, 40- and 72-bit codes, respectively. These are not Algorithm 3's final
ACE-selected outputs. The first two were chosen to match the paper's matrix
dimensions; the larger 72-bit matrix was retained for its improved measured FER.

The final ACE rule can prefer a smaller graph with worse decoding: for example,
for the 40-bit code it prefers the matrix obtained from 116 selected checks
with η̄_min = 6 to the one obtained from 174 checks with η̄_min = 5. For the
72-bit code, it selects the matrix obtained from 56 checks rather than 145.

This is an unresolved reproduction issue, not a replacement selection
algorithm or evidence of a general improvement over the paper's method.
Changing the initial row or row order can also change the greedy trajectory.
Seeds and selected-check counts alone do not reconstruct it without the exact
search-space ordering and tie policy.

## Replay from the selected inputs

The supplied inputs allow the final Algorithm 1 step to be repeated without
regenerating the search spaces or importing the research C++ backend. Save
this snippet as a Julia script and run it from the repository root with
`julia --project <script>`:

```julia
using PCMOpt, MAT

# Algorithm 1 with delta-only ties; retain the last pair in scan order.
function remove_delta_last(H)
    A = copy(H)
    N = ncols(H)
    while true
        best = -1
        Sv, Sc = falses(ncols(A)), Int[]
        for i in 1:nrows(A)-1, j in i+1:nrows(A)
            block = PCMOpt.pair_block(A, i, j)
            block === nothing && continue
            sv, sc, t = block
            delta = binomial(t, 2) * binomial(length(sc), 2)
            delta >= best && ((best, Sv, Sc) = (delta, sv, sc))
        end
        isempty(Sc) && return A
        PCMOpt.apply_block!(A, Sv, Sc, N)
    end
end

for k in (20, 40, 72)
    stem = "results/polar.128.$k.crc24c.pcm2"
    H = from_bool(Matrix{Bool}(matread(stem * ".input.mat")["H"]))
    A = k == 20 ? remove_4cycles(H) : remove_delta_last(H)
    @assert PCMOpt.to_bool(A) == matread(stem * ".mat")["H"]
    println("payload $k: exact replay verified")
end
```

This verifies the final factorization from the saved redundant input. It does
not rerun Algorithms 2 and 3. The complete search-space generation and PRG
exploration remain on the local research branches; the selected inputs are
the portable record needed to replay these particular final graphs.

## FER and energy normalization

The measurements below are from the construction investigation. The comparisons
use random codewords, flooding NMS with factor 0.7 and at
most 200 iterations. Auxiliary VNs have zero channel LLR. The 40- and 72-bit
curves use LLR saturation at 128. The 20-bit curve predates saturation;
a later saturated control was statistically consistent with it.

Agreement with Fig. 8 uses an inferred Eb/N0 convention counting payload plus
CRC bits: the normalization rates are 44/128, 64/128 and 96/128. The actual
code dimensions remain 20, 40 and 72. The paper does not explicitly state this
energy convention, so it remains an inference. With the current simulator's
payload-rate convention, obtain the same noise levels by adding
`10 log10((k+24)/k)` to the tabulated Eb/N0: approximately 3.42423, 2.04120 and
1.24939 dB, respectively. The simulator then labels its output with the
shifted, payload-normalized values.

| Payload | Eb/N0 (CRC-inclusive), dB | Frames | Errors | Candidate FER | Paper FER, approximate |
|---|---:|---:|---:|---:|---:|
| 20 | 1 | 100000 | 317 | 0.003170 | ≈0.003 |
| 20 | 2 | 1692000 | 150 | 0.00008865 | ≈0.0001 |
| 40 | 2 | 200000 | 2490 | 0.01245 | ≈0.01195 |
| 40 | 2.5 | 200000 | 453 | 0.002265 | ≈0.002761 |
| 40 | 3 | 534000 | 200 | 0.0003745 | ≈0.0005110 |
| 40 | 3.5 | 5738000 | 200 | 0.00003486 | ≈0.00005300 |
| 72 | 3 | 200000 | 4250 | 0.02125 | ≈0.01379 |
| 72 | 3.5 | 200000 | 934 | 0.004670 | ≈0.002088 |
| 72 | 4 | 284000 | 200 | 0.0007042 | ≈0.0002177 |

The 20-bit curve used seed 8731; the other curves used seed 984731. All used
eight threads and batches of 250 frames per thread. Paper values are readings
from Fig. 8, not the authors' raw samples; the 20-bit values are coarser visual
estimates. The exact matrices and selection trajectory remain unreproduced,
even where FER follows the paper closely.
