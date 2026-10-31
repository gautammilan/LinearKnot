# Architecture of the DP03 binary

One binary computes the DP03 partition function with pseudoknots (Dirks and Pierce 2003) and
its pair probabilities, either exactly (all beams infinite, 4-index tables) or in three beamed
linear passes (the paper setting). This file maps every source file to what it computes, the
beams, and the reference text.

## Build

```
make            # build/dp03, flags of paper-v1: -O3 -std=c++17 -ffp-contract=off -Wall -Wextra
make check      # V, G (states) and G (lazy) outside identity on RNase P D. radiodurans 1-150 (see below)
```

`cpp/main.cpp` is a unity build: it `#include`s the parts in the order of the original single
file (tag `paper-v1`, `cpp/dp03_l2r_merged.cpp` at 7c4148d), so the compiler sees the same code.
Parts under `cpp/pass*/`, `cpp/exact/`, `cpp/diag/` and `cpp/common/dp03_state.hpp` are
fragments of `class DP03`; parts under `cpp/cli/` are fragments of `main()`. None of them is a
standalone header. `cpp/dp03_l2r_merged.cpp` is a shim that includes `main.cpp`. At tag
`paper-v1-build` the shim build is byte-identical to the paper-v1 binary (md5
48d4b172c657d3825eef86d7a31b5793). After the dead-code pass (`docs/UNUSED_CODE.md`) the binary
differs; logZ on the paper preset and the identity check are unchanged.

Never add `-ffast-math` or FMA contraction: results are compared bit for bit.

## Passes and beams

| pass | what it computes | beam | flag | code variable |
|---|---|---|---|---|
| 1 | nested (pseudoknot-free) pre-pass, LinearPartition style: Q^b_pre, Q^m_pre, Q^m2_pre inside, nested outside, S̄_pre, C_pre | b1 | `--b1` | `nest_beam` |
| 2 | Q^z_pre (gap partition function) kept starts; nested tables Qb, Qm, Qm2, Qmp, Qp of the scan | b2 | `--b2` | `table_beam[BT_QB..BT_QP]`, `b2_width` |
| 2 | V^pre: future partner-stem profile, entries x_s.y_s per pair, cut by V S̄_pre[y_s+1] | b_V | `--bV` | `pr_vcap`, `p2_v_cap` |
| 2 | G (growth table G_L, Ĝ): entries (p', q') per state, cut by C_pre(p'-1) G S̄_pre(q'+1) | b_G | `--bG` | `g_ntr_K` |
| 2 / 3 | lazy G (`--g-sweep lazy`): pass 2 records the same hypergraph; the scan pulls H_s(q') = Σ_p' Ĉ(p'-1) G_s(p', q') per node, cut by H u_L(q') (`--bG-node-key ul`, default; `sbar`: H S̄_pre(q'+1)); the read keeps K_read right ends. `[g-read]` prints `g_mode node_K node_key read_K` | K (default 32), K_read (default K) | `--g-sweep lazy`, `--bG-node`, `--bG-node-key sbar\|ul`, `--bG-read` | `g_lz_K`, `g_lz_key`, `g_ntr_K` |
| 3 | L2R scan (LR-DP03 inside, `--pk-l2r on`): 4-index pseudoknot states ranked by the pass-2 reads | b4 | `--b4` | `beam`, `l2r_beam`, `rlt_beam` |
| 3 | scan partial states P1, P2, P3 | (b2, b4, b4) | derived | `l2r_pbeam`, `l2r_tb[PT_P1..]` |

Paper preset (`--preset accurate`): b1 200, b2 1000, b4 3000, b_V 256, b_G 256, P beams
(b2, b4, b4), G through multiloop states (`--g-sweep states`), linear pass-2 sweep
(`--sweep on --v-sweep on`). The energy model is not part of the preset:
`--material rna1999 --dangles some --param-file dp03/params/rna1999.dG`.
Width 0 (or `inf`) keeps every state.

## Files

References: **M** = `docs/LinearNupack_math.tex` (section / equation labels), **P2** =
`LINEAR_PASS2.md`, **R** = `docs/LinearNupack_report.tex` (M and R are in the main checkout,
not tracked yet). Paper section and equation numbers: to fill when the manuscript is in the
repository. `docs/UNUSED_CODE.md` records the dead-code pass: the coverage runs, what was
deleted (code no run executes, 181 options) and what was kept.

| file | computes | reference |
|---|---|---|
| `main.cpp` | entry, CLI usage text, preset expansion, per-sequence loop shell, output record | |
| `cli/options.hpp` | option variables and defaults | |
| `cli/parse_args.hpp` | the argument loop | |
| `cli/defaults_validate.hpp` | derived defaults (b1, b2, pbeam, b_V, sweep on), flag checks, beam audit | |
| `cli/run_setup.hpp` | per sequence: build DP03, pass 1 and the stem estimate | P2 Phase 0 |
| `cli/run_passes.hpp` | pass-2 setup (`stamp_p2`) on the pre-pass and the scan, pass 3 | P2 Phase 0 |
| `cli/run_report.hpp` | reports, checks, bpp, output | |
| `common/globals.hpp` | includes; global switches (g_sweep, g_ntr_K = b_G, traces) | |
| `common/logspace.hpp` | log_add, NEG_INF, RT | M sec:notation |
| `common/tables.hpp` | sparse tables, insertion-ordered maps, 4-index Tab4 | M sec:inside |
| `common/beam_tables.hpp` | table ids BT_*, state keys, beam statistics, drop log | M sec:beam |
| `common/dp03_state.hpp` | DP03 data members, run states, alloc, pass-3 Q^z filler, prune4 | M sec:beam, sec:two |
| `common/memlog.hpp`, `common/timers.hpp` | memory accounting (MemItem / msz) for the [mem-nt] line; CPU and peak-RSS timers | |
| `common/stem_est.hpp` | persistent log-sum index of the stem estimate (PersistentSum) | |
| `energy/parameters.hpp` | PkParams (beta1, beta2, beta3, ...), rna1999 / vienna2004 loaders | M sec:notation |
| `energy/loop_energy.hpp` | Energy: loop Boltzmann weights, dangles | M sec:notation |
| `pass1/inside.hpp` | nested inside at b1 (Q^b, Q^m, Q^m2, Qmp) | M sec:inside (two-index tables) |
| `pass1/prepass_products.hpp` | S, C0, nested outside, prS, prQz, PrMbar (S̄_pre, C_pre) | M sec:pre, sec:prpre, eq:prS, eq:prBml, eq:prctx |
| `pass2/qz.hpp` | Q^z, kept starts (top b2 by C(k-1) Q^z(k, j)), exterior C | P2 Phase 0 b), Phase 2 |
| `pass2/v_static.hpp`, `pass2/v_build.hpp` | V^pre build and the b_V cut (pr_vtrim) | M sec:prV, eq:prV, eq:prH; P2 Phase 0 c) |
| `pass2/arms.hpp` | arms: arm_raw / arm_phi give K^pre_L(q'; b, c) = Σ_{f>q'} Q^z_pre(q'+1, f-1) V, S̄_pk; argmax index, arm queries (pr_armq), arm selection (pr_select) | M eq:pra, eq:prK; P2 Phase 0 e) |
| `pass2/reads_k.hpp` | principled-B reads W, T, K_L, S_L / S_R, U_L / U_R and growth selection; the paper path runs only pr_ub_init and pr_ls_closings (it reads K_L through the scan's lK) | M sec:prstop, eq:prW, eq:prT, eq:prSL, eq:prK, eq:prSR, sec:prbound |
| `pass2/reads_v.hpp` | growth roots and keys read from V (`--grow-r-from-v`: rfv_root, rfv_key) | M sec:prG |
| `pass2/g_states.hpp` | G through the multiloop states (paper default); the b_G cut by C_pre(p'-1) G S̄_pre(q'+1) is cap_ga here. Under `--g-sweep lazy` the same loop records the edges (child, parent, w_e) and the stop terms into StemEst::lz* instead of pushing G maps | M sec:prG, eq:prG, eq:prR; P2 Phase 0 d) |
| `pass2/g_lazy.hpp` | lazy G: lz_get pulls H_s(q') = stop'_s(q') + Σ_{e: t→s} w_e H_t(q') on demand (explicit stack, memoized, cut to K by H S̄_pre(q'+1) or H u_L(q')); lz_out_id is the `--out-id` identity with Ĉ = C_pre | M sec:prG |
| `pass2/g_accum.hpp` | per-state G accumulators (racc_*, ntr_init) and the rec_* growth recursion | M eq:prKls, eq:prKrs |
| `pass2/sweep.hpp` | the linear V / G sweep (`--sweep on`), `--out-id` identity | |
| `pass2/reads_eval.hpp` | Ĝ reads for pass 3 (rec_root, rec_key), truncation, pr_eval* | M sec:prscore |
| `pass2/g_direct.hpp` | direct G, pair to pair through PrMbar (sw_fill) | |
| `pass2/growth_lazy.hpp` | stem steps (pr_steps, used by the pre-pass), distinct-key accumulators, profile finishing | |
| `pass3/scan_driver.hpp` | DP03 constructor, pfunc() (inside over j), bpp() | M sec:inside (order of a step) |
| `pass3/scan_inside.hpp` | L2R scan inside: PKM, growth, spawns, closings, P1-P3 | R Phase 1 |
| `pass3/scan_scores.hpp` | ranking scores read during the scan (arms, lK, S, lp_score) | M sec:prscore, eq:prO, eq:prQg, eq:prQgl |
| `pass3/scan_prune.hpp` | scan pruning (top b4 by score) | M sec:beam |
| `pass3/scan_steps_outside.hpp` | pseudoknot part of step j; scan outside (`l2r_*_back`) | R Phase 2 |
| `pass3/scan_helpers.hpp` | S̄_pk suffix, items, profile counters, beam caps | |
| `pass3/rank_score4.hpp`, `pass3/rank_prune.hpp` | 4-index state score and pruning (pr_score4, prune4_*, prune1 / prune2) | M sec:four, sec:prscore |
| `pass3/rank_options.hpp` | ranking option setters | |
| `pass3/rank_stem_variants.hpp`, `pass3/stem_est_tables.hpp` | stem-estimate scores σ_L, σ_R, σ_ext, helpers | M sec:sigL, sec:sigR, sec:grow, sec:helpgrow, sec:helpnest |
| `pass3/arm_gate.hpp` | per-table candidate counters of the Qg / Qgl rankings (GateStats) | |
| `pass3/outside.hpp` | outside and pair probabilities | M "Outside pass and pair probabilities"; R Phase 2 |
| `exact/pk_tables.hpp` | exact 4-index tables Qg, Qgls, Qgl, Qgr, Qp, Qgrs | M sec:inside (four-index tables) |
| `diag/memory.hpp`, `diag/scan_trace.hpp` | [nest] and [mem-nt] memory accounting; scan scores lsc_L / lsc_R. The identity check is `--out-id` (pass2/sweep.hpp), the pair traces `--g-trace` (pass2/g_states.hpp) | |

## The identity check

`--out-id` (in `pass2/sweep.hpp`) reconstructs the outside weight of every kept pair from the
pass-2 tables and compares it with pass 1's Q^b_b:
Σ_{p', q'} C_pre(p'-1) G(p', q') S̄_pre(q'+1) = Q^b_b(x, y), and the same for V. With b_V and
b_G uncapped on RNase P D. radiodurans 1-150 both hold to max relative error 7.8e-13 (V) and
7.4e-13 (G states). Under `--g-sweep lazy --bG-node 0` (in `pass2/g_lazy.hpp`) the check is
Σ_{q'} H_{x·y}(q') S̄_pre(q'+1) = Q^b_b(x, y) with C_pre as the prefix: 5.7e-13. `make check`
runs all three and fails above 1e-9.
