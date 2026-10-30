# Code and flags not used by the paper preset or the exact reference

Branch `LinearPass2`, code at tag `paper-v1` (moved, not changed). Nothing here is deleted; deletion waits for approval.

**Method.** A coverage build of the split code (`g++ -O1 --coverage`, same source) ran the paper preset (`--preset accurate --material rna1999 --dangles some`) on RNase P D. radiodurans 1-150 and RNase P M. maripaludus (233 nt). A function counts as unused if it never executed in either run. Line counts are source lines inside those functions that never executed (union over functions, so nested lambdas and overloads are counted once). Two sequences cannot reach every data-dependent branch, so treat a small never-called function as a candidate to check, not a proof.

**Exact reference.** In paper-v1 the exact command cannot be run as before: without the preset, the default b_V 256 trips `--qz-beam / --v-pairs / --v-cap ... need --stem-est-principled on`; with the preset, `--pk-l2r off --b1 0 --b2 0 --b4 inf --bV 0 --bG 0 --pr-v-key xsys --grow-r-from-v off` passes the checks and segfaults (the paper-v1 binary itself segfaults the same way, so this is not from the split). The stored exact references came from the older binary `cpp/dp03_cpp_stem` (`--beam inf`). The group "exact-reference path" below lists the code that path is expected to use; it was not exercised.

## Unused code by feature

| feature | functions | source lines never executed |
|---|---:|---:|
| Direct G (`--g-sweep direct`) with its cuts topq / collapse / entries (`--g-cut`), and its diagnostic `--g-eq` | 2 | 566 |
| Old lazy growth and principled-B reads (`--sweep off`: pr_steps / rec_node; W, T, K_L, S_L / S_R, U reads; grow recursion) | 64 | 1565 |
| Stem-estimate ranking variants and Beam-DP03 4-index ranking (`--pk-l2r off`; under `--pk-l2r on` these scores never run) | 96 | 1969 |
| Arm score, access file and gate (`--arm-score`, `--access-file`, `--gate*`) | 22 | 354 |
| L2R scan variants not in the preset (spawn modes other than inc, tails, `--grammar simple`, checkR, ninf, kx, ctx, arm fallback, `--l2r-p1`, split) | 47 | 1131 |
| Pass-2 V, arm and Q^z variants not in the preset (layered xsys key, V ratio / sample, persistent arm index, arm selection, sep, `pass2_apply` flags) | 23 | 663 |
| Energy models not used (compact dp03 table, vienna2004) | 12 | 188 |
| Exact-reference path (`--pk-l2r off`, all beams infinite): expected users, not exercised (see note) | 21 | 462 |
| Diagnostics, dumps, traces, profilers, `--help` | 75 | 3038 |
| Other | 3 | 11 |
| **total** | | **9947** |

A few lines sit in two groups (a lambda inside a function of another group), so this total is slightly above the per-file total below. Lines never executed inside used functions (untaken branches) are not in this table; see the last column of the per-file table.

## Per file

Executable lines are what gcov counts as code (no comments, braces or declarations).

| file | source lines | executable | executed | never-called functions | their never-executed lines | other never-executed lines (branches inside used functions) |
|---|---:|---:|---:|---:|---:|---:|
| `cli/defaults_validate.hpp` | 547 | 299 | 201 | 0 | 0 | 98 |
| `cli/options.hpp` | 143 | 38 | 38 | 0 | 0 | 0 |
| `cli/parse_args.hpp` | 1099 | 665 | 332 | 1 | 4 | 330 |
| `cli/run_passes.hpp` | 282 | 164 | 93 | 0 | 0 | 71 |
| `cli/run_report.hpp` | 934 | 481 | 239 | 1 | 17 | 228 |
| `cli/run_setup.hpp` | 328 | 249 | 178 | 1 | 23 | 51 |
| `common/beam_tables.hpp` | 142 | 65 | 4 | 7 | 80 | 6 |
| `common/dp03_state.hpp` | 517 | 260 | 114 | 12 | 198 | 9 |
| `common/gate_index.hpp` | 167 | 83 | 0 | 9 | 103 | 2 |
| `common/globals.hpp` | 147 | 19 | 5 | 2 | 5 | 10 |
| `common/logspace.hpp` | 43 | 21 | 5 | 2 | 21 | 0 |
| `common/memlog.hpp` | 161 | 96 | 32 | 8 | 72 | 12 |
| `common/stem_est.hpp` | 182 | 91 | 64 | 4 | 31 | 1 |
| `common/tables.hpp` | 389 | 195 | 180 | 3 | 11 | 6 |
| `common/timers.hpp` | 25 | 12 | 12 | 0 | 0 | 0 |
| `diag/dumps.hpp` | 511 | 348 | 0 | 7 | 458 | 0 |
| `diag/g_eq.hpp` | 305 | 193 | 0 | 1 | 296 | 0 |
| `diag/memory.hpp` | 278 | 209 | 165 | 2 | 30 | 21 |
| `diag/scan_trace.hpp` | 358 | 219 | 46 | 13 | 208 | 20 |
| `diag/stem_diag.hpp` | 1069 | 733 | 0 | 9 | 1007 | 0 |
| `diag/sweep_checks.hpp` | 369 | 245 | 0 | 2 | 355 | 0 |
| `diag/v_diag.hpp` | 324 | 204 | 2 | 5 | 249 | 13 |
| `energy/loop_energy.hpp` | 416 | 272 | 161 | 9 | 112 | 20 |
| `energy/parameters.hpp` | 347 | 149 | 98 | 3 | 76 | 2 |
| `exact/pk_tables.hpp` | 207 | 142 | 0 | 6 | 194 | 0 |
| `main.cpp` | 407 | 78 | 44 | 3 | 56 | 7 |
| `pass1/inside.hpp` | 609 | 438 | 336 | 4 | 67 | 52 |
| `pass1/prepass_products.hpp` | 379 | 256 | 166 | 1 | 11 | 81 |
| `pass2/arms.hpp` | 551 | 373 | 39 | 15 | 437 | 7 |
| `pass2/g_accum.hpp` | 622 | 428 | 136 | 16 | 310 | 61 |
| `pass2/g_direct.hpp` | 279 | 188 | 0 | 1 | 270 | 0 |
| `pass2/g_states.hpp` | 941 | 634 | 354 | 2 | 109 | 206 |
| `pass2/growth_lazy.hpp` | 615 | 350 | 170 | 14 | 207 | 49 |
| `pass2/qz.hpp` | 239 | 167 | 109 | 2 | 48 | 24 |
| `pass2/reads_eval.hpp` | 615 | 433 | 41 | 11 | 473 | 31 |
| `pass2/reads_k.hpp` | 810 | 503 | 35 | 20 | 595 | 20 |
| `pass2/reads_v.hpp` | 386 | 277 | 86 | 7 | 134 | 93 |
| `pass2/sweep.hpp` | 365 | 227 | 45 | 1 | 21 | 166 |
| `pass2/v_build.hpp` | 655 | 430 | 123 | 7 | 172 | 178 |
| `pass2/v_static.hpp` | 292 | 191 | 79 | 1 | 35 | 86 |
| `pass3/arm_gate.hpp` | 261 | 119 | 4 | 10 | 165 | 3 |
| `pass3/outside.hpp` | 613 | 436 | 243 | 2 | 41 | 157 |
| `pass3/rank_options.hpp` | 604 | 249 | 75 | 8 | 208 | 17 |
| `pass3/rank_prune.hpp` | 1178 | 806 | 74 | 33 | 871 | 61 |
| `pass3/rank_score4.hpp` | 269 | 170 | 9 | 13 | 210 | 1 |
| `pass3/rank_stem_variants.hpp` | 880 | 511 | 32 | 32 | 635 | 13 |
| `pass3/scan_driver.hpp` | 230 | 134 | 95 | 1 | 10 | 31 |
| `pass3/scan_helpers.hpp` | 303 | 197 | 59 | 10 | 117 | 43 |
| `pass3/scan_inside.hpp` | 1277 | 821 | 332 | 21 | 583 | 45 |
| `pass3/scan_prune.hpp` | 498 | 328 | 92 | 3 | 270 | 55 |
| `pass3/scan_scores.hpp` | 662 | 397 | 243 | 9 | 108 | 78 |
| `pass3/scan_steps_outside.hpp` | 553 | 393 | 210 | 5 | 202 | 18 |
| `pass3/stem_est_tables.hpp` | 152 | 46 | 24 | 6 | 25 | 3 |
| **total** | **24535** | **15032** | **5224** | **365** | **9940** | **2486** |

## Paper-preset flags whose scoring code never ran

The preset passes these, but under `--pk-l2r on` the functions they configure never executed (they still feed flag validation and parts of the pass-1 build): `--stem-est-qg-ext loop`, `--stem-est-grow-loop 2` (stem_grow_nodes, stem_ext_score), `--stem-est-helpers grow` (helper_*), `--stem-est-profile on` (prof_*, pair_prof), `--stem-est-chain on` (chain_*), `--stem-est rank` scoring (stem_B, prune4_stem), `--grow-rec on` (rec_node, rec_gbuild), `--pr-grow-select bound` (pr_gselect), `--pr-armq-sep on` (pr_armq).

## Flags not used by the paper preset or the exact command (181 of the flags parsed)

Used = in the expanded `--preset accurate` line, or in the exact command (`--loop-cap --max-loop --wc-ends --rt --material --dangles --param-file --gap-beam --cutoff --b4`, plus `--no-bpp` for long runs). Flags that `--preset` sets are used even if a script never types them.

`--out-id` is used by `make check` (V and G identity) and stays.

- **Already removed in paper-v1** (5): the parser only prints "... is removed; use ..." and exits 2; only the parse branch is left: `--grow-node-trunc`, `--pr-v-cap`, `--qz-beam`, `--qz-key`, `--v-cap`
- **Older names of beams the preset sets another way** (4): `--nest-beam` writes b1 (same as `--b1`), `--l2r-pbeam` writes pbeam, `--l2r-beam` / `--rlt-beam` write the scan beams that `--b4` seeds: `--nest-beam`, `--l2r-beam`, `--l2r-pbeam`, `--rlt-beam`
- **G build and cut variants** (9): `--g-cap`, `--g-cut`, `--g-cut-at`, `--g-eq`, `--g-fallback`, `--g-key`, `--g-mass`, `--g-mode`, `--g-pairs`
- **Arm score, access file, gate** (15): `--access-file`, `--arm-score`, `--gate`, `--gate-delta`, `--gate-diag`, `--gate-diag-taus`, `--gate-long-span`, `--gate-reserve`, `--gate-role`, `--gate-span`, `--gate-stats`, `--gate-tables`, `--gate-tau`, `--max-arm-len`, `--max-gap-scan`
- **Stem-estimate variants and their checks** (25): `--stem-bg-check`, `--stem-bg-reorder-check`, `--stem-diag`, `--stem-diag-loops`, `--stem-est-bg-reorder`, `--stem-est-chain-cut`, `--stem-est-chain-floor`, `--stem-est-chain-loop`, `--stem-est-chain-ml`, `--stem-est-chain-span`, `--stem-est-consts`, `--stem-est-ext-max`, `--stem-est-grow-cut`, `--stem-est-grow-nodes`, `--stem-est-hr`, `--stem-est-hr-gamma`, `--stem-est-ml-span`, `--stem-est-ml-tail`, `--stem-est-principled-nested`, `--stem-est-profile-eps`, `--stem-est-qg-combine`, `--stem-ext-check`, `--stem-hr-check`, `--stem-pr-check`, `--stem-profile-check`
- **L2R scan variants and diagnostics** (36): `--k-bound`, `--k-exact`, `--k-exact-check`, `--k-exact-ub`, `--k-memo`, `--l2r-arm-fb`, `--l2r-ctx`, `--l2r-ctx-rounds`, `--l2r-ctx-tau`, `--l2r-diag`, `--l2r-inc-check`, `--l2r-keep-ninf`, `--l2r-kx`, `--l2r-lean`, `--l2r-ninf-mass`, `--l2r-oracle`, `--l2r-oracle-dump`, `--l2r-p1`, `--l2r-p1-beam`, `--l2r-p1-qb`, `--l2r-p3-split`, `--l2r-p-aopen`, `--l2r-perf`, `--l2r-p-hk`, `--l2r-qz-audit`, `--l2r-spawn-k`, `--l2r-spawn-maxpop`, `--l2r-spawn-pool`, `--l2r-spawn-w`, `--l2r-split`, `--l2r-sub`, `--l2r-tbeam`, `--l2r-trace`, `--l2r-trace-key`, `--pk-l2r-checkR`, `--rlt`
- **Principled-B / pass-2 options** (50): `--grow-ml-diag`, `--grow-ml-L`, `--grow-ml-R`, `--grow-rec-check`, `--grow-trunc-eps`, `--lazy-ml-cap`, `--next-pair`, `--pr-abl-pk`, `--pr-aq-max`, `--pr-arm-gaponly`, `--pr-arm-persist`, `--pr-arms`, `--pr-arms-eps`, `--pr-arms-k`, `--pr-arms-max`, `--pr-arms-rank`, `--pr-arm-z`, `--pr-cut-by`, `--pr-cut-roots`, `--pre-beam`, `--pre-key`, `--pr-floor`, `--pr-grow-B`, `--pr-grow-check`, `--pr-ml-close`, `--pr-oracle`, `--pr-oracle-set`, `--pr-partner-check`, `--pr-partner-maxpairs`, `--pr-partner-ml`, `--pr-partner-ml-check`, `--pr-partner-span`, `--pr-score-dump`, `--pr-sep-check`, `--pr-span`, `--pr-static-only`, `--pr-test-qz`, `--pr-U-nobml`, `--pr-v2-dump`, `--pr-v-diag`, `--pr-v-mlsteps`, `--pr-v-ratio`, `--pr-v-ratio2`, `--pr-v-sample`, `--qzbar-old`, `--v-cut-fb`, `--v-floor`, `--v-pairs`, `--v-query`, `--v-trim`
- **Legacy beams and caps (superseded by b1, b2, b4, b_V, b_G)** (6): `--beam-table`, `--loop-cap-p1`, `--loop-cap-p2`, `--loop-cap-scan`, `--pass1-only`, `--pbeam-ratio`
- **Model variants** (10): `--no-gu`, `--no-pk`, `--no-single-stems`, `--pair-bonus`, `--params`, `--pk-dbeta1`, `--pk-dbeta1m`, `--pk-dbeta1p`, `--restrict-pairs`, `--score`
- **Diagnostics, dumps, traces, profiling** (20): `--beam-stats`, `--c-trace`, `--dump-exact`, `--dump-ext`, `--dump-kept`, `--dump-min-w`, `--dump-pruned`, `--dump-pruned-min-w`, `--dump-qb`, `--dump-scores`, `--dump-states`, `--g-trace`, `--help`, `--mem-log`, `--miss-sample`, `--mu-dump`, `--profile`, `--remove-state`, `--sample-prof`, `--sweep-check`

## Every never-called function

`file:start-end name` (source lines of the definition, including lines that are not code).

### Direct G (`--g-sweep direct`) with its cuts topq / collapse / entries (`--g-cut`), and its diagnostic `--g-eq`

```
diag/g_eq.hpp:9-304 g_eq_write
diag/g_eq.hpp:15-19 g_eq_write
diag/g_eq.hpp:20-29 g_eq_write
diag/g_eq.hpp:30-46 g_eq_write
diag/g_eq.hpp:80-85 g_eq_write
diag/g_eq.hpp:177-182 g_eq_write
diag/g_eq.hpp:199-213 g_eq_write
pass2/g_direct.hpp:9-278 sw_fill
pass2/g_direct.hpp:11-23 sw_fill
pass2/g_direct.hpp:24-29 sw_fill
pass2/g_direct.hpp:32-139 sw_fill
pass2/g_direct.hpp:38-44 sw_fill
pass2/g_direct.hpp:153-156 sw_fill
pass2/g_direct.hpp:164-173 sw_fill
pass2/g_direct.hpp:174-183 sw_fill
pass2/g_direct.hpp:213-215 sw_fill
pass2/g_direct.hpp:245-253 sw_fill
```

### Old lazy growth and principled-B reads (`--sweep off`: pr_steps / rec_node; W, T, K_L, S_L / S_R, U reads; grow recursion)

```
pass2/g_accum.hpp:69-73 g_pref
pass2/g_accum.hpp:74-80 racc_self
pass2/g_accum.hpp:103-119 racc_add_self
pass2/g_accum.hpp:105-111 racc_add_self
pass2/g_accum.hpp:276-281 racc_finish_trunc
pass2/g_accum.hpp:312-338 rv_finish
pass2/g_accum.hpp:339-351 rec_hg_trim
pass2/g_accum.hpp:353-361 rec_topc
pass2/g_accum.hpp:363-389 rec_helper_terms
pass2/g_accum.hpp:369-369 rec_helper_terms
pass2/g_accum.hpp:390-394 g_in_set
pass2/g_accum.hpp:396-425 g_one_step
pass2/g_accum.hpp:432-444 rec_helper_ub
pass2/g_accum.hpp:446-474 rec_helper
pass2/g_accum.hpp:477-520 rec_ml_ub
pass2/g_accum.hpp:492-492 rec_ml_ub
pass2/g_accum.hpp:523-553 rec_node
pass2/g_accum.hpp:570-587 rec_stat
pass2/g_accum.hpp:591-620 rec_gbuild
pass2/growth_lazy.hpp:15-31 gz_distinct
pass2/growth_lazy.hpp:34-39 pr_qz
pass2/growth_lazy.hpp:55-58 pbytes
pass2/growth_lazy.hpp:80-85 pr_abl_qz
pass2/growth_lazy.hpp:132-160 pr_finish_fut_d
pass2/growth_lazy.hpp:132-160 pr_finish_out_d
pass2/growth_lazy.hpp:167-171 pr_finish_out_d
pass2/growth_lazy.hpp:170-170 pr_finish_out_d
pass2/growth_lazy.hpp:175-177 pr_finish_fut_d
pass2/growth_lazy.hpp:184-189 PAcc::add
pass2/growth_lazy.hpp:195-215 PAcc::finish
pass2/growth_lazy.hpp:195-215 pr_finish_fut
pass2/growth_lazy.hpp:195-215 pr_finish_out
pass2/growth_lazy.hpp:223-227 pr_finish_out
pass2/growth_lazy.hpp:226-226 pr_finish_out
pass2/growth_lazy.hpp:228-234 pr_finish_fut
pass2/growth_lazy.hpp:231-233 pr_finish_fut
pass2/growth_lazy.hpp:284-284 pr_steps
pass2/growth_lazy.hpp:368-418 pr_gchild
pass2/growth_lazy.hpp:368-418 rfv_node
pass2/growth_lazy.hpp:379-379 pr_gchild
pass2/growth_lazy.hpp:379-379 rfv_node
pass2/growth_lazy.hpp:425-457 pr_gchild
pass2/growth_lazy.hpp:425-457 pr_static_ats
pass2/growth_lazy.hpp:425-457 rfv_node
pass2/growth_lazy.hpp:436-436 pr_gchild
pass2/growth_lazy.hpp:436-436 pr_static_ats
pass2/growth_lazy.hpp:436-436 rfv_node
pass2/growth_lazy.hpp:463-556 pr_gchild
pass2/growth_lazy.hpp:463-556 rfv_node
pass2/reads_eval.hpp:189-193 tr_uL
pass2/reads_eval.hpp:194-198 tr_uR
pass2/reads_eval.hpp:199-286 rec_trunc
pass2/reads_eval.hpp:242-267 rec_trunc
pass2/reads_eval.hpp:301-310 rec_tau
pass2/reads_eval.hpp:311-431 pr_eval_rec
pass2/reads_eval.hpp:435-465 rec_check
pass2/reads_eval.hpp:448-462 rec_check
pass2/reads_eval.hpp:476-533 pr_eval_sep
pass2/reads_eval.hpp:535-553 pr_eval
pass2/reads_eval.hpp:557-591 pr_gcheck_ep
pass2/reads_eval.hpp:595-614 pr_actx
pass2/reads_k.hpp:7-34 pr_W
pass2/reads_k.hpp:40-47 acc_begin
pass2/reads_k.hpp:48-55 acc_add
pass2/reads_k.hpp:68-122 pr_T
pass2/reads_k.hpp:135-178 pr_K
pass2/reads_k.hpp:180-183 pr_ctx
pass2/reads_k.hpp:195-220 pr_SLR_abl
pass2/reads_k.hpp:222-271 pr_SLR
pass2/reads_k.hpp:315-319 pr_iqz
pass2/reads_k.hpp:355-368 pr_ub_sync
pass2/reads_k.hpp:370-424 pr_ul
pass2/reads_k.hpp:426-442 pr_ur
pass2/reads_k.hpp:445-482 pr_U
pass2/reads_k.hpp:487-554 pr_gselect
pass2/reads_k.hpp:501-506 pr_gselect
pass2/reads_k.hpp:587-593 pr_hg_trim
pass2/reads_k.hpp:595-595 pr_gchild
pass2/reads_k.hpp:598-660 pr_gbuild
pass2/reads_k.hpp:662-696 pr_root
pass2/reads_k.hpp:700-712 pr_root_eval
pass2/reads_k.hpp:733-788 pr_key
pass2/reads_v.hpp:64-69 rfv_get
pass2/reads_v.hpp:70-74 rfv_at
pass2/reads_v.hpp:75-84 rfv_prof
pass2/reads_v.hpp:85-124 rfv_node
pass2/reads_v.hpp:95-107 rfv_node
```

### Stem-estimate ranking variants and Beam-DP03 4-index ranking (`--pk-l2r off`; under `--pk-l2r on` these scores never run)

```
common/dp03_state.hpp:131-143 l2r_simple_close
common/dp03_state.hpp:131-143 tagged_outside
common/stem_est.hpp:98-103 stemidx::PersistentSum::range_logsum
common/stem_est.hpp:141-147 stemidx::PersistentSum::query
common/stem_est.hpp:158-163 stemidx::MinTree::set_min
common/stem_est.hpp:164-175 stemidx::MinTree::query
pass1/prepass_products.hpp:89-99 build_nested_beta
pass3/rank_options.hpp:9-51 bml_check
pass3/rank_options.hpp:58-83 build_hr_av
pass3/rank_options.hpp:90-121 build_hr_phi
pass3/rank_options.hpp:125-149 nested_lphi
pass3/rank_options.hpp:547-599 stem_bg_check
pass3/rank_prune.hpp:9-35 prof_trunc
pass3/rank_prune.hpp:37-42 prof_finish_key
pass3/rank_prune.hpp:45-54 qg_prof
pass3/rank_prune.hpp:57-66 prof_closings
pass3/rank_prune.hpp:70-86 ml_tail_sum
pass3/rank_prune.hpp:89-103 ls_prof
pass3/rank_prune.hpp:106-118 rs_prof
pass3/rank_prune.hpp:123-174 prof_eval
pass3/rank_prune.hpp:180-219 prof_mrow
pass3/rank_prune.hpp:223-236 prof_R_reorder
pass3/rank_prune.hpp:240-246 ro_cmp
pass3/rank_prune.hpp:250-258 ro_flips
pass3/rank_prune.hpp:263-270 prune4_helper
pass3/rank_prune.hpp:263-270 prune4_helper
pass3/rank_prune.hpp:263-270 prune4_stem
pass3/rank_prune.hpp:263-270 prune4_stem
pass3/rank_prune.hpp:273-277 prof_cmp
pass3/rank_prune.hpp:281-407 prune4_helper
pass3/rank_prune.hpp:317-319 prune4_helper
pass3/rank_prune.hpp:332-334 prune4_helper
pass3/rank_prune.hpp:410-415 order_desc
pass3/rank_prune.hpp:420-442 split_mask
pass3/rank_prune.hpp:444-450 stem_A
pass3/rank_prune.hpp:454-480 stem_add_A
pass3/rank_prune.hpp:494-527 StemEval::StemEval
pass3/rank_prune.hpp:529-536 StemEval::H
pass3/rank_prune.hpp:539-547 StemEval::Hs
pass3/rank_prune.hpp:549-555 StemEval::B
pass3/rank_prune.hpp:557-562 StemEval::s1h
pass3/rank_prune.hpp:565-573 StemEval::sigL
pass3/rank_prune.hpp:576-603 StemEval::sigR
pass3/rank_prune.hpp:605-621 StemEval::score
pass3/rank_prune.hpp:626-801 prune4_stem
pass3/rank_prune.hpp:669-671 prune4_stem
pass3/rank_prune.hpp:685-687 prune4_stem
pass3/rank_prune.hpp:809-853 prune4_gated
pass3/rank_prune.hpp:809-853 prune4_gated
pass3/rank_prune.hpp:855-886 prune4_gated
pass3/rank_prune.hpp:863-866 prune4_gated
pass3/rank_prune.hpp:880-884 prune4_gated
pass3/rank_prune.hpp:924-959 prune4_step
pass3/rank_prune.hpp:962-968 score_ranks
pass3/rank_prune.hpp:972-988 note_drop
pass3/rank_prune.hpp:991-996 note_l2r_drop
pass3/rank_prune.hpp:1000-1010 note_kept
pass3/rank_score4.hpp:8-60 pr_score4
pass3/rank_score4.hpp:62-64 pr_tab
pass3/rank_score4.hpp:99-104 lse
pass3/rank_score4.hpp:107-137 prn_step
pass3/rank_score4.hpp:139-148 prn_paths
pass3/rank_score4.hpp:150-157 prn_G
pass3/rank_score4.hpp:159-166 prn_V
pass3/rank_score4.hpp:168-175 prn_ctx
pass3/rank_score4.hpp:180-198 prn_future
pass3/rank_score4.hpp:202-217 prn_partner
pass3/rank_score4.hpp:219-229 prn_outqg
pass3/rank_score4.hpp:231-262 pr_naive4
pass3/rank_score4.hpp:264-268 prc_cmp
pass3/rank_stem_variants.hpp:9-19 stem_B
pass3/rank_stem_variants.hpp:54-61 stem_B_cached
pass3/rank_stem_variants.hpp:67-84 stem_Bg_raw
pass3/rank_stem_variants.hpp:86-95 stem_Bg
pass3/rank_stem_variants.hpp:101-101 hr_h
pass3/rank_stem_variants.hpp:102-102 hr_c
pass3/rank_stem_variants.hpp:107-131 hr_naive
pass3/rank_stem_variants.hpp:137-169 hr_naive_phi
pass3/rank_stem_variants.hpp:139-143 hr_naive_phi
pass3/rank_stem_variants.hpp:171-178 stem_H_cached
pass3/rank_stem_variants.hpp:191-301 stem_grow_nodes
pass3/rank_stem_variants.hpp:239-241 stem_grow_nodes
pass3/rank_stem_variants.hpp:306-326 stem_ext_score
pass3/rank_stem_variants.hpp:330-349 ext_paths
pass3/rank_stem_variants.hpp:351-375 stem_ext_naive
pass3/rank_stem_variants.hpp:381-396 helper_stop
pass3/rank_stem_variants.hpp:400-411 helper_O
pass3/rank_stem_variants.hpp:422-435 closings_ls
pass3/rank_stem_variants.hpp:465-469 helper_sigma
pass3/rank_stem_variants.hpp:473-502 helper_naive
pass3/rank_stem_variants.hpp:507-524 helper_nested_sum
pass3/rank_stem_variants.hpp:527-533 helper_sigma_nested
pass3/rank_stem_variants.hpp:536-552 helper_nested_naive
pass3/rank_stem_variants.hpp:560-565 ProfAcc::add
pass3/rank_stem_variants.hpp:568-582 ProfAcc::finish
pass3/rank_stem_variants.hpp:596-619 pair_prof
pass3/rank_stem_variants.hpp:603-612 pair_prof
pass3/rank_stem_variants.hpp:643-701 chain_nbrs
pass3/rank_stem_variants.hpp:667-667 chain_nbrs
pass3/rank_stem_variants.hpp:704-713 chain_stop
pass3/rank_stem_variants.hpp:716-723 chain_merge
pass3/rank_stem_variants.hpp:726-740 chain_trunc
pass3/rank_stem_variants.hpp:744-772 chain_build
pass3/rank_stem_variants.hpp:776-794 chain_root
pass3/rank_stem_variants.hpp:798-804 chain_paths
pass3/rank_stem_variants.hpp:810-841 chain_naive
pass3/stem_est_tables.hpp:60-63 StemEst::QbIn
pass3/stem_est_tables.hpp:68-72 StemEst::lean_abort
pass3/stem_est_tables.hpp:113-116 StemEst::H
pass3/stem_est_tables.hpp:118-121 StemEst::H2
pass3/stem_est_tables.hpp:122-125 StemEst::B
pass3/stem_est_tables.hpp:126-129 StemEst::Mbar
```

### Arm score, access file and gate (`--arm-score`, `--access-file`, `--gate*`)

```
common/dp03_state.hpp:356-402 arm_front
common/dp03_state.hpp:366-381 arm_front
common/dp03_state.hpp:404-407 softplus
common/dp03_state.hpp:409-443 arm_score_of
common/gate_index.hpp:28-53 gate::PartnerSums::PartnerSums
common/gate_index.hpp:59-68 gate::PartnerSums::r
common/gate_index.hpp:80-80 gate::PersistentGate::PersistentGate
common/gate_index.hpp:81-99 gate::PersistentGate::PersistentGate
common/gate_index.hpp:102-107 gate::PersistentGate::range_max
common/gate_index.hpp:119-124 gate::PersistentGate::node
common/gate_index.hpp:126-132 gate::PersistentGate::build
common/gate_index.hpp:134-143 gate::PersistentGate::update
common/gate_index.hpp:145-151 gate::PersistentGate::query
common/gate_index.hpp:155-165 gate::mirror_pairs
pass3/arm_gate.hpp:19-24 acc_log_add
pass3/arm_gate.hpp:37-40 AccessRec::log_acc
pass3/arm_gate.hpp:42-45 AccessRec::pair_prob
pass3/arm_gate.hpp:47-57 AccessRec::pairs
pass3/arm_gate.hpp:67-75 read_line
pass3/arm_gate.hpp:77-111 index_access_file
pass3/arm_gate.hpp:113-187 load_access
pass3/arm_gate.hpp:224-229 Gate::Gate
pass3/arm_gate.hpp:231-240 Gate::roles
pass3/arm_gate.hpp:242-246 Gate::g
```

### L2R scan variants not in the preset (spawn modes other than inc, tails, `--grammar simple`, checkR, ninf, kx, ctx, arm fallback, `--l2r-p1`, split)

```
pass3/scan_helpers.hpp:7-38 l2r_pk_suffix
pass3/scan_helpers.hpp:11-24 l2r_pk_suffix
pass3/scan_helpers.hpp:41-51 l2r_items
pass3/scan_helpers.hpp:52-58 l2r_sfx_pass
pass3/scan_helpers.hpp:220-225 l2r_cap_beam
pass3/scan_helpers.hpp:227-239 l2r_yw_rebuild
pass3/scan_helpers.hpp:255-257 l2r_smap
pass3/scan_helpers.hpp:258-262 l2r_tab
pass3/scan_helpers.hpp:263-273 l2r_reset_b
pass3/scan_helpers.hpp:274-288 l2r_collect_b
pass3/scan_helpers.hpp:289-302 l2r_split_run
pass3/scan_inside.hpp:224-229 tail_extend
pass3/scan_inside.hpp:243-264 l2r_prune_prefix
pass3/scan_inside.hpp:266-293 l2r_prune_tail
pass3/scan_inside.hpp:295-321 l2r_tail
pass3/scan_inside.hpp:323-339 l2r_tail_back
pass3/scan_inside.hpp:342-374 l2r_close_tails
pass3/scan_inside.hpp:345-351 l2r_close_tails
pass3/scan_inside.hpp:529-541 l2r_spawn_back
pass3/scan_inside.hpp:565-591 l2r_spawn_step_cols
pass3/scan_inside.hpp:592-630 l2r_spawn_cube_step
pass3/scan_inside.hpp:632-653 ltr_pool
pass3/scan_inside.hpp:654-697 ltr_rank
pass3/scan_inside.hpp:698-729 l2r_spawn_rank
pass3/scan_inside.hpp:736-751 l2r_spawn_cols
pass3/scan_inside.hpp:752-769 l2r_spawn_row
pass3/scan_inside.hpp:770-829 l2r_spawn_cube
pass3/scan_inside.hpp:775-781 l2r_spawn_cube
pass3/scan_inside.hpp:791-812 l2r_spawn_cube
pass3/scan_inside.hpp:832-886 l2r_spawn
pass3/scan_inside.hpp:925-937 lc_cmp
pass3/scan_inside.hpp:964-1000 l2r_checkR
pass3/scan_inside.hpp:1168-1194 l2r_inc_check_step
pass3/scan_inside.hpp:1226-1237 ninf_note
pass3/scan_inside.hpp:1240-1274 l2r_ninf_report
pass3/scan_prune.hpp:81-116 keep_split
pass3/scan_prune.hpp:207-291 l2r_pruneA_exact
pass3/scan_prune.hpp:233-247 l2r_pruneA_exact
pass3/scan_scores.hpp:27-41 growth_eager
pass3/scan_scores.hpp:67-78 lstack_best
pass3/scan_scores.hpp:79-98 larm_fb
pass3/scan_scores.hpp:145-156 lKx_sum
pass3/scan_scores.hpp:157-168 lKx
pass3/scan_scores.hpp:170-181 lKx_memo
pass3/scan_scores.hpp:210-217 l2r_ctx_init
pass3/scan_scores.hpp:226-237 lkbar
pass3/scan_scores.hpp:568-572 lpf_lap
pass3/scan_steps_outside.hpp:41-73 l2r_simple_m4
pass3/scan_steps_outside.hpp:42-63 l2r_simple_m4
pass3/scan_steps_outside.hpp:76-98 l2r_simple_m5
pass3/scan_steps_outside.hpp:86-94 l2r_simple_m5
pass3/scan_steps_outside.hpp:102-169 l2r_simple_close
pass3/scan_steps_outside.hpp:154-162 l2r_simple_close
pass3/scan_steps_outside.hpp:171-215 l2r_simple_pk
pass3/scan_steps_outside.hpp:217-249 l2r_simple_end
```

### Pass-2 V, arm and Q^z variants not in the preset (layered xsys key, V ratio / sample, persistent arm index, arm selection, sep, `pass2_apply` flags)

```
pass2/arms.hpp:8-21 pr_gap_begin
pass2/arms.hpp:22-25 pr_gap_covered
pass2/arms.hpp:145-159 pidx_pull
pass2/arms.hpp:160-182 pidx_set
pass2/arms.hpp:183-197 pidx_q
pass2/arms.hpp:198-232 pidx_build
pass2/arms.hpp:237-243 aq_better
pass2/arms.hpp:244-248 aq_tree
pass2/arms.hpp:249-257 aq_query
pass2/arms.hpp:259-342 pr_armq
pass2/arms.hpp:347-440 pr_select
pass2/arms.hpp:364-364 pr_select
pass2/arms.hpp:399-422 pr_select
pass2/arms.hpp:443-473 pr_sep_extend
pass2/arms.hpp:479-501 pr_sep_ext
pass2/g_states.hpp:120-126 sw_hyper_push
pass2/g_states.hpp:560-564 sw_hyper_push
pass2/g_states.hpp:565-616 sw_hyper_push
pass2/g_states.hpp:660-664 sw_hyper_push
pass2/g_states.hpp:721-725 sw_hyper_push
pass2/g_states.hpp:764-791 pass2_apply
pass2/g_states.hpp:849-856 pass2_apply
pass2/sweep.hpp:50-54 sw_as_rn
pass2/sweep.hpp:55-70 sw_as_rn
pass2/v_build.hpp:10-50 pr_finish_layers
pass2/v_build.hpp:17-19 pr_finish_layers
pass2/v_build.hpp:260-265 pr_vrand
pass2/v_build.hpp:270-292 pr_vsuf_build
pass2/v_build.hpp:295-333 pr_vb_group
pass2/v_build.hpp:336-347 pr_vr_group
pass2/v_build.hpp:431-443 pr_static_ats
pass2/v_build.hpp:529-534 pr_static_ats
pass2/v_static.hpp:40-44 pr_build_static
pass2/v_static.hpp:101-115 pr_build_static
pass2/v_static.hpp:161-175 pr_build_static
```

### Energy models not used (compact dp03 table, vienna2004)

```
energy/loop_energy.hpp:287-292 Energy::stack
energy/loop_energy.hpp:294-303 Energy::bulge
energy/loop_energy.hpp:305-312 Energy::internal
energy/loop_energy.hpp:318-327 Energy::vtype
energy/loop_energy.hpp:329-329 Energy::vext
energy/loop_energy.hpp:331-348 Energy::hairpin04
energy/loop_energy.hpp:351-390 Energy::interior04
energy/loop_energy.hpp:394-406 Energy::branch04
energy/loop_energy.hpp:409-414 Energy::closing04
energy/parameters.hpp:104-107 size_continuation
energy/parameters.hpp:128-131 param_error
energy/parameters.hpp:272-339 load_v04
energy/parameters.hpp:277-285 load_v04
```

### Exact-reference path (`--pk-l2r off`, all beams infinite): expected users, not exercised (see note)

```
common/dp03_state.hpp:131-143 fill_Qg
common/dp03_state.hpp:131-143 outside_from_future
common/dp03_state.hpp:244-251 run_weight
common/dp03_state.hpp:253-267 save_runs
common/dp03_state.hpp:254-262 save_runs
common/dp03_state.hpp:335-341 qz_log
common/dp03_state.hpp:348-352 four_score
common/dp03_state.hpp:447-516 prune4
exact/pk_tables.hpp:8-98 fill_Qg
exact/pk_tables.hpp:81-88 fill_Qg
exact/pk_tables.hpp:100-128 fill_Qgls
exact/pk_tables.hpp:130-146 fill_Qgl
exact/pk_tables.hpp:148-162 fill_Qgr
exact/pk_tables.hpp:164-184 fill_Qp
exact/pk_tables.hpp:186-206 fill_Qgrs
pass1/inside.hpp:422-426 build_qlead
pass1/inside.hpp:544-553 set_Qmp
pass1/inside.hpp:555-572 add_branch_p
pass1/inside.hpp:575-608 fill_Qmp
pass2/qz.hpp:7-52 fill_Qz
pass2/qz.hpp:211-215 add2b
pass3/outside.hpp:197-209 outside_from_future
pass3/outside.hpp:349-376 outside_same_j
```

### Diagnostics, dumps, traces, profilers, `--help`

```
cli/parse_args.hpp:29-32 lambda
cli/run_report.hpp:184-187 lambda
cli/run_report.hpp:188-191 lambda
cli/run_report.hpp:804-812 lambda
cli/run_setup.hpp:125-148 lambda
common/beam_tables.hpp:20-27 BeamStats::open
common/beam_tables.hpp:29-43 BeamStats::mass_share
common/beam_tables.hpp:44-47 BeamStats::log_Q
common/beam_tables.hpp:48-57 BeamStats::log_prune
common/beam_tables.hpp:59-81 BeamStats::log_prune_mask
common/beam_tables.hpp:82-96 BeamStats::flush_summary
common/beam_tables.hpp:112-116 beam_table_id
common/globals.hpp:19-19 on_cpu_stop
common/globals.hpp:81-84 miss_urand
common/logspace.hpp:19-26 g_mass_note
common/logspace.hpp:27-39 g_mass_flush
common/memlog.hpp:20-22 msz::vec
common/memlog.hpp:20-22 msz::vec
common/memlog.hpp:20-22 msz::vec
common/memlog.hpp:20-22 msz::vec
common/memlog.hpp:20-22 msz::vec
common/memlog.hpp:20-22 msz::vec
common/memlog.hpp:24-28 msz::vec2
common/memlog.hpp:30-33 msz::umap
common/memlog.hpp:30-33 msz::umap
common/memlog.hpp:40-47 mem_rss_mb
common/memlog.hpp:48-51 mem_heap_mb
common/memlog.hpp:106-126 sample_on_sig
common/memlog.hpp:127-151 sample_start
common/memlog.hpp:152-160 sample_stop
diag/dumps.hpp:7-11 json_num
diag/dumps.hpp:14-24 std::__cxx11::basic_string
diag/dumps.hpp:14-24 std::__cxx11::basic_string
diag/dumps.hpp:26-33 json_str
diag/dumps.hpp:41-126 gate_diag
diag/dumps.hpp:133-184 dump_states
diag/dumps.hpp:139-150 dump_states
diag/dumps.hpp:151-159 dump_states
diag/dumps.hpp:202-379 dump_scores
diag/dumps.hpp:207-212 dump_scores
diag/dumps.hpp:238-238 dump_scores
diag/dumps.hpp:239-268 dump_scores
diag/dumps.hpp:285-297 dump_scores
diag/dumps.hpp:393-510 dump_ext
diag/dumps.hpp:401-413 dump_ext
diag/memory.hpp:23-44 log_pr_side
diag/memory.hpp:115-122 vmap_sz
diag/scan_trace.hpp:9-75 l2r_oracle_dump
diag/scan_trace.hpp:76-109 l2r_oracle_load
diag/scan_trace.hpp:113-116 lors
diag/scan_trace.hpp:118-123 lS
diag/scan_trace.hpp:200-222 lsc_R
diag/scan_trace.hpp:229-232 ldiag_put
diag/scan_trace.hpp:244-246 ltr_match4
diag/scan_trace.hpp:247-249 ltr_match3
diag/scan_trace.hpp:250-264 ltr_line
diag/scan_trace.hpp:265-273 ltr_step4
diag/scan_trace.hpp:274-302 ltr_final
diag/scan_trace.hpp:336-340 vi_lap
diag/scan_trace.hpp:350-355 lpf_prune_stat
diag/stem_diag.hpp:11-150 tagged_outside
diag/stem_diag.hpp:64-71 tagged_outside
diag/stem_diag.hpp:155-168 stemdiag::avg_ranks
diag/stem_diag.hpp:171-185 stemdiag::pearson
diag/stem_diag.hpp:188-201 stemdiag::top_b
diag/stem_diag.hpp:203-208 stemdiag::num
diag/stem_diag.hpp:215-307 pr_build_hbar
diag/stem_diag.hpp:317-530 pr_hybrid_outside
diag/stem_diag.hpp:327-327 pr_hybrid_outside
diag/stem_diag.hpp:379-384 pr_hybrid_outside
diag/stem_diag.hpp:404-416 pr_hybrid_outside
diag/stem_diag.hpp:543-829 pr_oracle_diag
diag/stem_diag.hpp:577-579 pr_oracle_diag
diag/stem_diag.hpp:607-623 pr_oracle_diag
diag/stem_diag.hpp:845-1068 stem_diag
diag/stem_diag.hpp:852-857 stem_diag
diag/stem_diag.hpp:1044-1051 stem_diag
diag/sweep_checks.hpp:11-175 sw_verify
diag/sweep_checks.hpp:21-30 sw_verify
diag/sweep_checks.hpp:31-36 sw_verify
diag/sweep_checks.hpp:50-53 sw_verify
diag/sweep_checks.hpp:57-64 sw_verify
diag/sweep_checks.hpp:124-128 sw_verify
diag/sweep_checks.hpp:179-368 sw_growth_check
diag/sweep_checks.hpp:185-194 sw_growth_check
diag/sweep_checks.hpp:195-200 sw_growth_check
diag/sweep_checks.hpp:201-206 sw_growth_check
diag/sweep_checks.hpp:214-365 sw_growth_check
diag/sweep_checks.hpp:225-228 sw_growth_check
diag/sweep_checks.hpp:229-232 sw_growth_check
diag/sweep_checks.hpp:236-243 sw_growth_check
diag/sweep_checks.hpp:298-302 sw_growth_check
diag/v_diag.hpp:8-95 pr_vdiag_static
diag/v_diag.hpp:47-54 pr_vdiag_static
diag/v_diag.hpp:102-239 pr_vdiag_T
diag/v_diag.hpp:164-190 pr_vdiag_T
diag/v_diag.hpp:231-231 pr_vdiag_T
diag/v_diag.hpp:276-281 CTree::better
diag/v_diag.hpp:282-288 CTree::build
diag/v_diag.hpp:290-299 CTree::query
main.cpp:257-277 usage
main.cpp:280-303 load_states
main.cpp:305-315 parse_list
pass2/arms.hpp:66-97 k_bound_check
pass2/arms.hpp:505-550 pr_sep_diag
pass2/reads_eval.hpp:99-179 dump_p2_all
pass2/reads_eval.hpp:106-127 dump_p2_all
pass2/reads_eval.hpp:128-139 dump_p2_all
pass2/reads_v.hpp:16-23 v_cut_fb_mu
pass2/reads_v.hpp:149-196 trace_v
pass2/reads_v.hpp:151-158 trace_v
pass2/reads_v.hpp:200-216 vquery_dump
pass2/v_build.hpp:210-242 pr_v2_dump
pass3/rank_options.hpp:296-315 pr_sep_report
pass3/rank_options.hpp:323-327 psc_bin
pass3/rank_options.hpp:486-489 PrPartnerTimer::~PrPartnerTimer
pass3/scan_driver.hpp:189-198 dump_qb
pass3/scan_prune.hpp:349-497 l2r_diag_report
```

### Other

```
common/tables.hpp:209-217 SpMat::dense
common/tables.hpp:230-230 SpMat::Proxy::operator=(SpMat::Proxy const&)
common/tables.hpp:233-233 SpMat::operator[](unsigned long) const
```

