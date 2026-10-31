# Dead-code pass: what was deleted and how it was checked

Branch `LinearPass2`. Starting point: commit 6652db2 (code of tag `paper-v1-build`, which builds the
paper-v1 binary byte for byte). Code that no run below executes was deleted; values are unchanged.
50 files under `cpp/` changed, about 10,950 lines removed.

## Runs

Coverage build of the starting tree (`g++ -O1 --coverage`, with a SIGSEGV / SIGABRT handler that
flushes the counters). Paper preset = `--preset accurate --material rna1999 --dangles some
--param-file dp03/params/rna1999.dG`.

| run | input | command | rc |
|---|---|---|---|
| Ther | tmRNA Ther.ther. TRW-300852 1-349 | paper preset | 0 |
| Deha | tmRNA Deha.ethe. CP000027 1-352 | paper preset | 0 |
| Pseu | tmRNA Pseu.aeru. AE004091 1-353 | paper preset | 0 |
| RP486 | RNase P D. radiodurans (486 nt) | paper preset | 0 |
| HXB2_3000 | HIV-1 HXB2 1-3000 | paper preset | see "HXB2" below |
| exact100_nobpp | RNase P D. radiodurans 1-100 | paper preset + `--pk-l2r off --b1 0 --b2 0 --b4 inf --bV 0 --bG 0 --pr-v-key xsys --grow-r-from-v off --no-bpp` (the exact command) | 0 |
| exact100 | same | the exact command with pair probabilities | 139 (segfault in outside, also in paper-v1; not fixed here) |
| exact100_plain | same | the exact command without the preset | 2 (rejected by the flag checks, as in paper-v1) |
| RP150_id | RNase P D. radiodurans 1-150 | paper preset + `--bV 10000000 --bG 10000000 --out-id --no-bpp` (identity check) | 0 |
| RP150_gtrace | same | paper preset + `--g-trace 20,120 --out-id --no-bpp` (pair traces) | 0 |

A line counts as executed if any run executed it (union of the gcov counts).

## What was deleted

- Functions, branches and members no run executes: direct-G cuts (topq / collapse / entries,
  `--g-eq`), the old lazy growth and principled-B variants, Beam-DP03 ranking variants not used
  by the exact command, arm score / access file / gate, scan variants (spawn modes other than
  inc, ctx, kx, ninf, arm fallback, `--l2r-p1`, split, sub, tails), V / arm / Q^z variants not
  used (ratio, sample, persistent arm index, sep), the sampling profiler and the memory log.
- Diagnostics and dumps, except the identity check (`--out-id`, pass2/sweep.hpp) and the pair
  traces (`--g-trace`, pass2/g_states.hpp), which are kept. Deleted files: `diag/dumps.hpp`,
  `diag/g_eq.hpp`, `diag/stem_diag.hpp`, `diag/sweep_checks.hpp`, `diag/v_diag.hpp` (its ArmQ
  struct moved to `pass2/arms.hpp`), `common/gate_index.hpp`.
- Options: all of the 181 options not used by the preset or the exact command, except
  `--g-trace` (pair traces) and `--help`, plus three preset options whose code never runs
  (`--stem-est-chain`, `--stem-est-profile`, `--grow-rec`). 182 removed; 57 remain. The preset
  expansion no longer emits the three preset options. The removed list is below.
- Members or globals that no run writes were folded to their constant value. Report lines keep
  their text: where a removed variable was printed, the constant is printed instead (for
  example `diag 0 U_nobml 0` in `[growth]`, `g_cut topq`), so report fields stay the same.
  On RP150 stderr differs from the stored run only in the `[preset] expanded:` line, the
  `[deviation]` lines of removed features (gate, arm score, access, max_gap_scan, max_arm_len,
  v-floor, v-cut-fb) and timings / RSS.

Kept on purpose although they look optional: numeric tunables read on executed paths (`pr_loop`,
`chain_floor`, `g_ml_L`, `g_ml_R`, `g_tr_eps`, `l2r_spawn_k`), the `[mem-nt]` and `[nest]`
reports (`mem_items`, `stem_mem_items`, `log_nest`) and the tables they read (`prQbIn`, `prKept`),
and `--l2r-spawn`.

## Check that only unexecuted code went

1. **Forward check.** Every original line that was deleted or replaced was mapped to its
   coverage count (difflib per file). 10,991 lines went; 919 of them had a nonzero count. All
   919 fall into: option-parsing branches of removed options (the `if (a == "--x")` test itself
   runs for every argument), assignments of folded constants and their declarations, `if`
   lines unwrapped because the condition is now constant (the body stays, re-indented; e.g. the
   `g_U_nobml || g_grow_rec` block in `pr_ub_init`), printf argument lists rewritten to literals,
   and code moved between files. Largest counts: cli/ (512 lines), rank_options (62),
   scan_prune and scan_scores (37 each).
2. **Reverse check.** The new tree was built with coverage and run on the small commands
   (exact100_nobpp, exact100_plain, RP150_id, RP150_gtrace, RP150 paper preset). Every executed
   new line was mapped back; an executed line whose original line never ran would mean a removed
   guard. This found one: the `dp.outside()` call of the deleted `DUMP_P2_ALL` block had been
   left in `run_report.hpp` and segfaulted the exact command (the same outside() crash as
   exact100 above). It is removed; the reverse check then found none.
3. **Values.** Release build (`-O3 -std=c++17 -ffp-contract=off`) of the final tree against the
   stored outputs: logZ and all non-timing stdout lines identical on Ther, Deha, Pseu, RP150 and
   RP486 (`/tmp/new_exactness/growth_read/step2/*/out.txt`); identity check identical (V 2948
   pairs max_rel 7.816e-13, G states 7.390e-13); exact100_nobpp, RP150_id, RP150_gtrace stdout
   identical to the coverage runs, and the 41 `[g-trace]` lines and the `[out-id]` lines
   identical. HXB2 1-3000: see below.

## HXB2

HXB2_3000 coverage: pending at the time of writing; see the task file
`tasks/2026-10-30-dead-code-and-refs.md`.

## Removed options

`--access-file` `--arm-score` `--beam-stats` `--beam-table` `--c-trace` `--dump-exact` `--dump-ext`
`--dump-kept` `--dump-min-w` `--dump-pruned` `--dump-pruned-min-w` `--dump-qb` `--dump-scores`
`--dump-states` `--g-cap` `--g-cut` `--g-cut-at` `--g-eq` `--g-fallback` `--g-key` `--g-mass`
`--g-mode` `--g-pairs` `--gate` `--gate-delta` `--gate-diag` `--gate-diag-taus` `--gate-long-span`
`--gate-reserve` `--gate-role` `--gate-span` `--gate-stats` `--gate-tables` `--gate-tau`
`--grow-ml-L` `--grow-ml-R` `--grow-ml-diag` `--grow-node-trunc` `--grow-rec` `--grow-rec-check`
`--grow-trunc-eps` `--k-bound` `--k-exact` `--k-exact-check` `--k-exact-ub` `--k-memo`
`--l2r-arm-fb` `--l2r-beam` `--l2r-ctx` `--l2r-ctx-rounds` `--l2r-ctx-tau` `--l2r-diag`
`--l2r-inc-check` `--l2r-keep-ninf` `--l2r-kx` `--l2r-lean` `--l2r-ninf-mass` `--l2r-oracle`
`--l2r-oracle-dump` `--l2r-p-aopen` `--l2r-p-hk` `--l2r-p1` `--l2r-p1-beam` `--l2r-p1-qb`
`--l2r-p3-split` `--l2r-pbeam` `--l2r-perf` `--l2r-qz-audit` `--l2r-spawn-k` `--l2r-spawn-maxpop`
`--l2r-spawn-pool` `--l2r-spawn-w` `--l2r-split` `--l2r-sub` `--l2r-tbeam` `--l2r-trace`
`--l2r-trace-key` `--lazy-ml-cap` `--loop-cap-p1` `--loop-cap-p2` `--loop-cap-scan`
`--max-arm-len` `--max-gap-scan` `--mem-log` `--miss-sample` `--mu-dump` `--nest-beam`
`--next-pair` `--no-gu` `--no-pk` `--no-single-stems` `--pair-bonus` `--params` `--pass1-only`
`--pbeam-ratio` `--pk-dbeta1` `--pk-dbeta1m` `--pk-dbeta1p` `--pk-l2r-checkR` `--pr-U-nobml`
`--pr-abl-pk` `--pr-aq-max` `--pr-arm-gaponly` `--pr-arm-persist` `--pr-arm-z` `--pr-arms`
`--pr-arms-eps` `--pr-arms-k` `--pr-arms-max` `--pr-arms-rank` `--pr-cut-by` `--pr-cut-roots`
`--pr-floor` `--pr-grow-B` `--pr-grow-check` `--pr-ml-close` `--pr-oracle` `--pr-oracle-set`
`--pr-partner-check` `--pr-partner-maxpairs` `--pr-partner-ml` `--pr-partner-ml-check`
`--pr-partner-span` `--pr-score-dump` `--pr-sep-check` `--pr-span` `--pr-static-only`
`--pr-test-qz` `--pr-v-cap` `--pr-v-diag` `--pr-v-mlsteps` `--pr-v-ratio` `--pr-v-ratio2`
`--pr-v-sample` `--pr-v2-dump` `--pre-beam` `--pre-key` `--profile` `--qz-beam` `--qz-key`
`--qzbar-old` `--remove-state` `--restrict-pairs` `--rlt` `--rlt-beam` `--sample-prof` `--score`
`--stem-bg-check` `--stem-bg-reorder-check` `--stem-diag` `--stem-diag-loops`
`--stem-est-bg-reorder` `--stem-est-chain` `--stem-est-chain-cut` `--stem-est-chain-floor`
`--stem-est-chain-loop` `--stem-est-chain-ml` `--stem-est-chain-span` `--stem-est-consts`
`--stem-est-ext-max` `--stem-est-grow-cut` `--stem-est-grow-nodes` `--stem-est-hr`
`--stem-est-hr-gamma` `--stem-est-ml-span` `--stem-est-ml-tail` `--stem-est-preset`
`--stem-est-principled-nested` `--stem-est-profile` `--stem-est-profile-eps`
`--stem-est-qg-combine` `--stem-ext-check` `--stem-hr-check` `--stem-pr-check`
`--stem-profile-check` `--sweep-check` `--v-cap` `--v-cut-fb` `--v-floor` `--v-pairs` `--v-query`
`--v-trim`
