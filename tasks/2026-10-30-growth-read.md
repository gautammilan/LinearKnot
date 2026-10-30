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

## Ther diagnosis (signed dlogZ vs NUPACK exact 327.54328498495346)

Binary `dp03_diag7` (adds `--g-cut-at all|pair`, `--g-trace P,Q[,X,Y]`, and the `--out-id` breakdowns
below). Ther, paper preset, `--g-sweep states` unless noted.

| run | logZ | dlogZ vs exact | G fill s | peak RSS MB |
|---|---:|---:|---:|---:|
| topq bG 32 (default) | 327.39777562 | -0.145509 | 1617 | 3212 |
| direct entries bG 128 | 327.39737672 | -0.145908 | 254 | 1489 |
| states bG 128 | 314.78248426 | -12.760801 | 90 | 1471 |
| states no cut (bG 1e8) | 327.39782478 | **-0.145460** | 679 | 8280 |
| states bG 128, cut at pair states only | 314.78248430 | -12.760801 | 163 | 1748 |
| states bG 256 | 327.39766236 | -0.145623 | 165 | 1829 |
| states bG 512 | 327.39764313 | -0.145642 | 319 | 2624 |

No cut is above topq by 4.9e-5, so the loss is the cut, not the nested growth model.

Which cut: Ther's main crossing stems from the bpp (1-based) are A 31-203 .. 40-194 (min p 0.95)
and B 41-296 .. 48-289 (0.88). 0-based trace at bG 128:

- A, entry (30,202), root (39,193): the entry exists only in the map of pair (30,202) itself,
  where it is the stop-term seed, and is dropped there: rank 158 of 498, score 127.72 vs the 128th
  131.63. Seen at 1 pair state, 0 Qm, 0 Qm2. G(39,193) ends with 128 entries, none (30,202).
- B, entry (40,295), root (47,288): same, dropped at pair (40,295), rank 182 of 434.

The cut scores the seed as an exterior-loop pair, C_pre(p'-1) * E.at * S_pre(q'+1). In the nested
model (30,202) sits inside the 12-314 helix (p 0.99), so its exterior score is low next to the
entries pushed from that helix. For a pseudoknot stem the seed's real outside is the pseudoknot,
not the nested exterior. Cutting only at pair states gives the same loss (Qm/Qm2 cuts are not
involved); bG 256 keeps the seed and recovers.

## Item 4: the 1,107 pairs and the 48,426 entries (RP150, uncapped)

- The earlier "1,107 pairs with no pass-1 outside" was wrong. All 4,045 pairs in both tables have
  finite Qb_b. The 1,097 pairs outside the 2,948 checked are GU/UG pairs (`end_ok` false, cannot be
  terminal), which the check skipped. Counting them, G states matches Qb_b on all 4,045 (GU too).
  The infinite max for direct came from 40 pairs where direct's G is empty (see next item).
- The 48,426 states-only entries are direct's topq cut. It runs even when the map is under K and
  drops every q' with ntrUL[q'] = -inf (e.g. q' = n-1; the largest are stop terms x.149, share 1.0
  of their pair). With `--g-cut entries` on the direct side the entry sets are equal (0 / 0).
- Values still differ. Direct is below Qb_b on every pair that differs (none above): median
  log(direct/Qb_b) -1.13 (end_ok) / -1.11 (GU), min -9.8. The nearest mismatches are interior
  loops with 6 unpaired bases (-0.0011); at 5% the nearest is a 7-base left gap (-1.94). The
  multiloop constants are the same values under rna1999/DP03 (a_p = alpha1 = 3.4, b_p = alpha2 =
  0.4, c_p = alpha3 = 0) and pair_bonus is 0, so they are not the cause; the remaining gap is
  not yet localized (dangles none can't test it: the identity fails for V and G states there too).

## New default: G through multiloop states, bG 256 (Step 2)

`--g-sweep states` is the default; `--preset accurate` writes `--bG 256` (was 32). Cut key unchanged
(top-bG raw (p',q') by C_pre*G*S_pre). Paper preset: b1 200, b2 1000, b4 3000, bV 256, bG 256,
pbeam (b2,b4,b4), linear pass-2 sweep, other paper flags as before. Binary `/tmp/new_exactness/bin/dp03_g256`.
Accept rule: signed dlogZ vs exact >= old cut's - 0.01 (RNase P has no exact; compared on logZ).

| seq | n | logZ | dlogZ vs exact | old cut (topq bG 32) | max dp | G fill s | scan CPU s | total CPU s | peak RSS MB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Ther | 349 | 327.397662 | -0.1456 | -0.1455 | 0.049 | 165 | 162 | 425 | 1829 |
| Deha | 352 | 233.580012 | -1.5804 | -1.5817 | 0.627 | 144 | 173 | 407 | 1805 |
| Pseu | 353 | 215.533277 | -0.2557 | -0.2534 | 0.100 | 151 | 273 | 513 | 1978 |
| RP150 | 150 | 116.281324 | logZ -1.0e-5 vs old | 116.281334 | - | 10 | 23 | 43 | 453 |
| RP486 | 486 | 440.339333 | logZ +0.3102 vs old | 440.029104 | - | 375 | 584 | 1145 | 3733 |

All pass. Old cut for comparison (same table format; its G fill is lazy inside the scan):

| seq | G fill s | scan CPU s | total CPU s | peak RSS MB |
|---|---:|---:|---:|---:|
| Ther | 1552 | 1741 | 1803 | 3212 |
| Deha | 1203 | 1388 | 1453 | 3036 |
| Pseu | 1137 | 1372 | 1436 | 2892 |
| RP150 | 31 | 54 | 64 | 509 |
| RP486 | 5062 | 5595 | 5722 | 7809 |

New G fill runs in pass 2 (pass 2 CPU: Ther 214, Deha 187, Pseu 194, RP150 12, RP486 472 s).
