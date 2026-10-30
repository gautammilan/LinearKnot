# 2026-10-30 LinearPass2: split for readability, values unchanged

Tag `paper-v1` = 7c4148d, the commit of the binary used for the paper runs of 2026-10-30
(`/tmp/new_exactness/bin/dp03_default` and `dp03_g256`, md5 48d4b172c657d3825eef86d7a31b5793,
`g++ -O3 -std=c++17 -ffp-contract=off -Wall -Wextra`). Branch `LinearPass2` starts there.

Note: at paper-v1, `cpp/gate_index.hpp` and `cpp/stem_est.hpp` (included by
`dp03_l2r_merged.cpp`) were not tracked, so the tag alone did not build. The split commit adds
them as `cpp/common/gate_index.hpp` and `cpp/common/stem_est.hpp` (same content).

## Step 1: split (6f04eae)

`cpp/dp03_l2r_merged.cpp` (23,820 lines) moved into 50 parts under `cpp/{common,energy,pass1,
pass2,pass3,exact,diag,cli}/`, included by `cpp/main.cpp` in the original order (unity build).
Largest part 1,277 lines. Only comments were added (a header per part).

## Step 2: same values as paper-v1

1. Code: `g++ -E -P` of the split code equals that of paper-v1, byte for byte (rechecked after
   the comment edit in step 4).
2. Binary: built through the shim `cpp/dp03_l2r_merged.cpp`, the split binary has the md5 of
   the paper binary, 48d4b172c657d3825eef86d7a31b5793. Built from `cpp/main.cpp` (`make`), it
   differs only in the source-file name in the symbol table.
3. Runs: the `main.cpp` build at the paper preset (`--preset accurate --material rna1999
   --dangles some --param-file rna1999.dG`), compared with the paper-v1 outputs of
   `/tmp/new_exactness/growth_read/step2/`. Every output line is identical except `sec`
   (CPU time).

| sequence | logZ split | logZ paper-v1 | output lines | identical |
|---|---:|---:|---:|---|
| Ther | 327.39766236390301 | 327.39766236390301 | 1193 | yes |
| Deha | 233.58001232046206 | 233.58001232046206 | 2969 | yes |
| Pseu | 215.53327699794767 | 215.53327699794767 | 5307 | yes |
| RP150 (RNase P D. radiodurans 1-150) | 116.28132370133829 | 116.28132370133829 | 1572 | yes |
| RP486 (RNase P D. radiodurans) | 440.33933303544205 | 440.33933303544205 | 3550 | yes |

4. Identity (RP150, `--bV 10000000 --bG 10000000 --out-id`): V pairs 2948 max_rel 7.816e-13,
   G-states pairs 2948 max_rel 7.390e-13, the same numbers as paper-v1. Now `make check`.

Exact reference: not rerun. In paper-v1 the plain exact command fails validation (default b_V 256
needs `--stem-est-principled on`), and the preset with `--pk-l2r off --b1 0 --b2 0 --b4 inf
--bV 0 --bG 0` segfaults; the paper-v1 binary segfaults the same way, so it predates the split.
The stored exact logZ values come from the older `cpp/dp03_cpp_stem --beam inf`.

## Step 3: unused code (5bd3a10)

`docs/UNUSED_CODE.md`, from a coverage build on the paper preset (RP150 and M. maripaludus):
365 never-called functions, 9,940 source lines, plus 2,486 never-executed lines inside used
functions; 181 of the parsed flags are not used by the preset or the exact command (5 of them
already only print "removed"). Nothing deleted; waiting for approval.

## Steps 4 and 5 (f4887c2, 408dabc)

`docs/ARCHITECTURE.md`, `.cursor/rules/architecture.mdc` (always applied); `make check` runs
`tests/identity_check.sh` (75 s, 0.6 GB; `PARAM_FILE=` if `dp03/params/rna1999.dG` is not in
the checkout, since it is not tracked).

## Step 6: push

Not done: the repository has no remote and `gh` is not installed.
