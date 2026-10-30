# 2026-10-30 growth-read: G on V's sweep

Branch `perf/growth-read`. Default remains `--g-sweep direct` and `--g-cut topq`.

## 1. Direct multiloop step is not beam-bounded

`sw_fill` (paper `--sweep on`) takes the `qm_ix` path: for a pair ending at `y`, the right-gap positions are the kept `PrMbar` right ends `j` with left end `y+1` and `j <= y + S`, `S = --pr-ml-span` (100000). The other path (every `yo` in the span) is not used.

`S` does not bind on these lengths. The count is how many columns kept that left end. Nested beam `b1 = 200` caps starts per column, not this list, so one pair can see more than 200 right ends.

| seq | pairs | yo mean | yo max | path | S |
|---|---:|---:|---:|---|---:|
| Ther n=349 | 17266 | 89.26 | 333 | qm_ix | 100000 |
| RNase P 1-150 | 2948 | 40.30 | 136 | qm_ix | 100000 |
| HXB2 1-3000 | running | | | | |

Ther direct entries `--bG 128`: logZ `327.39737671903737` (same as the earlier entries run).

## E3 multi / E4 multi-ext

Those edges are gated on `l2r_mode == 0`. The scan is `--pk-l2r on`, but `sw_hyper_push` runs on the nested pre-pass, which is constructed with the default `l2r_mode == 0` and `allow_pk == 0`. RP150 states log: `[sweep-kept] l2r_mode 0 allow_pk 0`. The pushes are that nested Multi path, not the scan.

## 2. `--g-sweep states|direct`

- `direct` (default): `sw_fill` as before.
- `states`: `do_g = true` in `sw_hyper_push` on the same Qm/Qm2/pair edges as V. Stop at pair `x.y` is `G[x,y] += E.at(x,y)` (no prefix). Every state keeps the top `bG` `(p',q')` by `C_pre(p'-1) * G * S_pre(q'+1)`. The scan still adds `Chat` at read (`sw_as_rn`).

## 3. Uncut check — stopped

RNase P 1-150, paper preset, `--bG 10000000` (pre max = post max on both, so the cap did not drop entries).

| | logZ |
|---|---|
| direct | 116.28133400511338 |
| states | 116.28132371946509 |
| states − direct | −1.029e-5 |

Threshold was 1e-9. Steps 4 and 5 were not run.

## Correction: logZ was the wrong test — outside identity instead

`--out-id` (RNase P 1-150, paper preset, `--bV 10000000 --bG 10000000`). For each kept end_ok pair,
compare the reconstructed outside to pass 1's `Qb_b`. Stop term is `E.at(x,y)`, the exterior-loop
pair factor pass 1 adds into `Qb_b`.

| table | pairs | max rel | median rel |
|---|---:|---:|---:|
| V | 2948 | 7.8e-13 | 5.7e-14 |
| G states | 2948 | 7.4e-13 | 5.7e-14 |
| G direct | 4045 | 1e300 (pairs with no Qb_b) | 0.715 |

G direct does not satisfy the identity: it walks `pr_steps` (pseudoknot-stem multiloop
parameters, `PrMbar` = nested Qm, pair_bonus), not the nested grammar.

States vs direct, uncapped, per entry: 1,925,690 in both, 48,426 only in states, 0 only in direct,
max relative difference 2.8e9.

V does not push one multiloop twice. In the nested pre-pass (`l2r_mode == 0`) the pair→Qm2 close
(E9, `l2r_mode != 0`) is off; the close goes pair→Multi (also counted as E9), then Multi→Qm2 (E3 multi)
and Multi→Multi (E4 multi-ext). E5–E8 are branch edges Qm/Qm2→pair. V's identity at 7.8e-13 confirms.

## Paper matrix (bG as listed; signed dlogZ vs topq bG 32)

| seq | direct entries 128 | states 128 | G pushes topq / direct / states | G fill s topq / direct / states | scan CPU s topq / direct / states | peak RSS MB topq / direct / states |
|---|---:|---:|---|---|---|---|
| Ther | -0.0004 | **-12.615** | 22.5M / 22.6M / 9.7M | 1617 / 254 / 90 | 1738 / 393 / 137 | 3212 / 1489 / 1471 |
| Deha | +0.0012 | +0.0004 | 21.0M / 21.1M / 9.7M | 1260 / 230 / 85 | 1387 / 378 / 145 | 3036 / 1516 / 1496 |
| Pseu | -0.0005 | -0.0074 | 20.2M / 20.2M / 9.3M | 1194 / 249 / 103 | 1371 / 486 / 239 | 2892 / 1517 / 1543 |
| RP150 | -0.0000 | -0.0000 | 2.35M / 2.38M / 1.27M | 33 / 15 / 7.5 | 54 / 37 / 23 | 509 / 441 / 438 |
| RP486 | +0.1055 | +0.2743 | running / 49.1M / 20.7M | running / 594 / 236 | running / 1061 / 481 | running / 2496 / 2732 |

RP486 topq logZ from the earlier gcut matrix (440.02910384486495); the rerun is still going.
NUPACK refs: Ther 327.5433 (topq -0.146, direct -0.146, states -12.761), Deha 235.1604, Pseu 215.7890.
