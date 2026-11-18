# LinearKnot

LinearKnot computes the pseudoknotted partition function and base-pair probabilities of an RNA sequence under the Dirks and Pierce (2003) model (DP03, the NUPACK pseudoknot grammar) in time and memory that scale linearly with sequence length. It replaces the O(n^5) exact DP03 recursions with three beam-pruned left-to-right passes in the style of LinearPartition.

- Input: RNA sequences, one per line, optional `>name` headers.
- Output: log partition function `logZ` and base-pair probabilities, including pseudoknotted (crossing) pairs.
- Energy model: NUPACK 3 `rna1999` (Turner 1999) nested-loop parameters with DP03 pseudoknot penalties.

## Requirements

- A C++17 compiler (tested with g++ 13) and `make`.
- Linux or macOS. No other dependencies (one header-only library is vendored in `cpp/third_party/`).

## Build

```bash
git clone https://github.com/gautammilan/LinearKnot.git
cd LinearKnot
make            # builds build/dp03
make check      # optional self-test, about 1.5 min and 0.6 GB
```

## Quick start

```bash
printf ">example\nGGGAAACCCAAAGGGAAACCC\n" > example.txt
build/dp03 < example.txt > example.out 2> example.log
```

`example.out`:

```
>example
logZ 13.222636047479279
sec 0.001885
pairs 31
1 8 0.79930398516105061
2 7 0.81022478379273788
...
```

With no flags, LinearKnot runs with the beams used in the paper. Keep `build/dp03` inside the repository so it finds its parameter file (`dp03/params/rna1999.dG`).

## Input format

- One sequence per line. Lines starting with `>` are names and are copied to the output.
- Do not wrap long sequences over several lines; each line is read as a separate sequence.
- Upper or lower case is fine; `T` is read as `U`.

## Output format

Per sequence, on stdout:

| line | meaning |
|---|---|
| `>name` | the header, if given |
| `logZ x` | natural log of the partition function |
| `sec t` | wall time for this sequence, in seconds |
| `pairs N` | number of pairs listed below |
| `i j p` | pair (i, j), 1-based, with probability p >= 1e-5 |

Diagnostics go to stderr; redirect them with `2> file`.

## Beams

| flag | meaning | default |
|---|---|---|
| `--b1` | beam of the nested pre-pass (pass I) | 200 |
| `--b2` | beam of the two-index states (pass III) | 100 |
| `--b4` | beam of the four-index pseudoknot states (pass III) | 3000 |
| `--bV` | beam of V (pass II) | 256 |
| `--bG` | beam of G (pass II) | 32 |

Any beam not given keeps its default. `--b4` controls the main accuracy and speed trade-off, for example `build/dp03 --b4 1000 < seqs.txt`. See `docs/ARCHITECTURE.md` for what each beam prunes.

Rough speed: a 150 nt sequence takes about 25 s on one core. The binary is single-threaded; to run many sequences, split the input and run several processes in parallel.

## Reproducing the paper

The LinearKnot rows of the main accuracy table were run with:

```bash
build/dp03 < sequences.txt > predictions.txt            # b4 = 3000 row
build/dp03 --b4 1000 < sequences.txt > predictions.txt  # b4 = 1000 row
```

Structures are taken from the pair probabilities with ThreshKnot (threshold 0.3).

## Repository layout

| path | contents |
|---|---|
| `cpp/main.cpp` | entry point; a unity build that includes all parts |
| `cpp/pass1/`, `cpp/pass2/`, `cpp/pass3/` | the three beamed passes |
| `cpp/energy/` | energy parameter loaders |
| `cpp/exact/` | exact four-index DP03 tables |
| `dp03/params/rna1999.dG` | NUPACK rna1999 parameters |
| `tests/` | `make check` script and test sequence |
| `docs/ARCHITECTURE.md` | file-by-file map of the algorithm and beams |

## Citation

The paper describing LinearKnot is under review. Citation details will be added here.

## Contact

Milan Gautam, Oregon State University. Please open a GitHub issue for questions or bugs.
