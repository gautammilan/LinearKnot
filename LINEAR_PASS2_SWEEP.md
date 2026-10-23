# Linear pass 2, one outside sweep

Flags off stays the current preset (`--preset default --material rna1999 --dangles some`).
Pass 2 only ranks. Nothing that accumulates `Z` in the scan changes.
Rows go to `stemest/pass2/sweep.csv`.

Bit-identity targets, flags off:

- 5S 119 nt `logZ 96.343443758865405`
- RNase P 486 nt `logZ 440.529136884630`

This note is appended after each step. Step 0 is read-only except the profile line.

## Step 0. What the engine does

Source: `cpp/dp03_l2r_merged.cpp`. Preset is LR-DP03, rna1999, dangles some, loop cap 30. Nested pre-pass is this same grammar with pseudoknots off (`allow_pk` false), not a separate LinearPartition binary.

### a) Pre-pass states, beam, C_pre, S̄_pre

Inside states the nested run actually fills, per right end `j`:

| state | role |
|---|---|
| `Qb(i,j)` | a pair. Hairpin, interior loop of span ≤ 30 onto a kept inner pair, or a multiloop close |
| `Qm(i,j)` | ≥ 1 multiloop branch, then unpaired bases |
| `Qm2(i,j)` | the subset of `Qm` with ≥ 2 branches. This is what a pair closes |
| `Qz(i,j)` | nested-in-a-gap chain. Not a multiloop state |
| `C(j)` | exterior inside on `0..j` |

There is no `M1` state and no stored `multi` state. `Qmp` exists only when `stem_ml` is on. Dirks–Pierce 2003, which rna1999 uses, has `stem_ml` false, so the nested multiloop table is `Qm`, not `Qmp`.

`Qm` is not built one branch at a time as a stored chain. `fill_Qm` does two things:

- extend every surviving `Qm` / `Qm2` at `j-1` by one unpaired base (`alpha3`, and with dangles the 3-state run);
- `add_branch` from each pair `d.j`: a closed form over the unpaired bases before that branch, `i` from `d - ml_left` to `d`, into `Qm` only; plus, for every kept `Qm(i, d-1)`, one new branch into both `Qm` and `Qm2`.

`ml_left` is 2^30 for rna1999, so that unpaired run is the whole prefix. A pair closes with `Qm2(i+1, j-1) + alpha1 + alpha2 + at_close` (`fill_Qb`). Dangles some makes the unpaired extension a P/U1/U2 run, not a single factor `u`.

`--nest-beam b` (and the automatic beam 100 for `n > 1000`) constructs that nested DP with beam `b` and runs it. `prune2` then keeps, independently, the top `b` starts of `Qb`, of `Qm`, and of `Qm2` at each right end. The key is `C(i-1) + alpha(i,j)`, and `C` at that moment is the exterior inside of this same run, so it is the prefix key. A dropped start is set to log 0 and is never read again. It does not record the hyperedges that built a state, and there is no inside-only key. The scan's own beams are separate (`--beam`, the P beams). `--nest-beam` does not touch them.

`C_pre` is that nested `C`. `S̄_pre` is not the outside of the exterior state. `build_principled` recomputes it from the inside `Qb`, right to left: suffix from `y` to `n`, with `prS[n] = 0`. On an exact pre-pass that equals the exterior outside. On a beam it is the suffix of the kept inside, not `C_b`.

### b) Multiloop weights in V and in growth

Both use `pr_steps`. A step from the current pair `x.y` to an enclosing pair `xo.yo` is one of:

| kind | weight |
|---|---|
| interior | `stem_interior`, both gaps ≤ 30 |
| multiloop, `Qm` on the left, unpaired on the right | `Qm(xo+1, x-1) + (w_bp + at(x,y)) + w_cp*(yo-y-1) + run + close` |
| unpaired on the left, `Qm` on the right | the symmetric term, left gap ≤ 30 |
| `Qm` on both sides | both interiors, no unpaired-run factor |

`close = w_ap + w_bp + at(xo,yo) + pair_bonus`. For these parameters `w_ap, w_bp, w_cp` are 3.4, 0.4, 0.0, the same numbers as `alpha1, alpha2, alpha3`, but they are the pseudoknot-stem multiloop parameters, not the nested ones. The left `Qm` during the scan is the scan's own `Qm` (`qmp_at`). Off the scan, and for the right side always, it is `PrMbar`, the nested pre-pass `Qm`. The right-hand span is `--pr-ml-span`, which the preset sets to 100,000, so it is the whole sequence. That is the direct pair-to-pair jump. At 486 nt it was 996e6 of the 1.02e9 growth steps.

`Qm2` is not what these steps read. They read `Qm`, which still contains the one-branch states.

### c) Where V and growth are read

`V^pre` is built before the scan in `pr_build_static`, one map per pair, by the same `pr_steps` (interior plus the direct multiloop jump). The preset then stores the trimmed map as flat arrays (`prVkX`, `prVkW`). U-type beam keys read it through `rfv_root` / `rfv_key` (`l2r_pruneB`), because `--grow-r-from-v` is on and `G_R` is not built. The L side of an R key does not read `V`. `lsc_L` dots the growth vector with `lK`, and `lK` / `lYwS` are arm sums built from those stored `V` entries times `S̄_pre`.

Growth is read from `rec_root` in `l2r_pruneA` (R-type beam key) and in `lh_get` (the H list inside P3). Spawn would also call it; the preset's spawn width is 0.

### d) Profile timers

`[profile]` is unchanged and still inclusive: `R` and `P3` contain the growth built while they run. `[profile-x]` is the exclusive split of the scan. `R` and `P3` there have `grow_R` and `grow_H` removed, up to the size of the bucket. `rest` is scan minus Qz, growth, the exclusive tables, and close. `over` is growth that the R or P3 bucket did not actually contain, so it could not be removed. `V` and init stay in the before-scan time and are not in this line. Printing it does not change `Z`.

## Step 1. Stopped

`--pre-beam` / `--pre-key` would need a beamed pre-pass that, for every kept state, stores the inside hyperedges, and an `M1` / `M` / `multi` factorization whose outside push is one branch at a time. `--nest-beam` does not do that.

What it does instead: one beam width, the prefix key only, applied separately to `Qb`, `Qm`, and `Qm2`, with no hyperedge list. The one-branch object is not a state. It is the closed-form loop inside `add_branch` (`i` from `d - ml_left` to `d`), and the ≥2-branch object is a second copy (`Qm2`) rather than `multi(x,y) = sum_k M * M1`. `S̄_pre` is recomputed from inside `Qb`, not from the exterior outside.

The sweep in Step 2 pushes along those hyperedges. Implementing it on `Qm` / `Qm2` would not be the grammar in the request, and implementing the requested grammar beside the real one would be a different partition function for the ranking tables. Stopped here, before any of that was coded.

Flags off was not changed except the extra `[profile-x]` line. 5S flags off after the sweep-check build: `logZ 96.343443758865405`.

## Step 1, continued. The sweep follows `pr_steps`, not a new grammar

`--nest-beam` still does not record hyperedges, and `Qm` / `Qm2` are not `M1` / `M` / `multi`. The sweep does not replace those tables. It rebuilds the pass-0 future-stem profile by walking the same steps `pr_steps` already uses, outermost pair first, with the multiloop jump replaced by a chain:

- A closing pair `xo.yo` (`end_ok`) deposits its profile times `close = w_ap + w_bp + at(xo,yo) + pair_bonus`.
- The right side is three buckets per left end `xo`. Gap 0 is the pair at `yo = y+1`. Gap 1 is `yo = y+2`, times `s1[y+1] + w_cp`. Longer gaps sit in one summed profile per `xo`, stored with `d3[yo-1] + w_cp*(yo-1)`, and the common `d5[y+1] - w_cp*y` is applied when the inner pair at right end `y` is reached. That is `E.run` for dangles some; `w_cp` is 0 for these parameters, and the factoring still holds when it is not.
- Left of the inner pair `x.y`: each finite `PrMbar(i, x-1)` (the nested `Qm`, span at most `--pr-ml-span`) multiplies that right-side profile by `br = w_bp + at(x,y)`. That is `pr_steps` case A, and the same left `Qm` times the right-`Qm` bucket is case C.
- Case B is the right-`Qm` bucket times an unpaired left gap of at most the loop cap, weight `w_cp*(x-xo-1) + run`.
- Interior steps stay inside the loop cap: `stem_interior + pair_bonus`, pushed from the outer pair onto the inner one.
- The stop at a pair `x.y` is `at(x,y) + beta1 + Sbar[y+1]` on key `x` (`--pr-v-key xs`).

`--sweep-check` builds this beside the current `V` and compares every pair and every key. The current trim is off (`--pr-v-cap 0`, preset `--pr-v-eps 0`). Max relative difference on the linear weights:

| n | pairs | entries | max rel |
|---|---|---|---|
| 40 (5S prefix) | 214 | 1832 | 1.78e-14 |
| 150 (RNase P prefix) | 4045 | 170414 | 1.01e-12 |
| 300 | 15738 | 1475442 | 5.40e-12 |

The check does not feed the scan, so `Z` is whatever the preset produces with the cap removed. Flags off, cap left at 512, is the 5S value above.

## Step 4b. Stopped

Caps off (`--pr-v-cap 0 --grow-node-trunc 0 --grow-ml-c 0 --grow-r-from-v off`), 5S prefix, 40 nt. The chain's growth, summed onto `q'` with the scan's final `C`, was compared with every `rec_root` layer vector (98 roots, 507 finite entries). Indices are 0-based.

Using the scan's `Qm` on the left, the same table `rec_root` reads: max relative difference `8.88e-15`.

Using the pre-pass `Qm` (`PrMbar`), which is all a sweep before the scan can see: the first root that moves is pair `15.28`, `q' = 30`, sweep `-1.42641505`, `rec_root` `-1.42641504`. The largest is pair `23.32`, `q' = 34`, sweep `2.26543508`, `rec_root` `2.26544029`, relative `5.22e-6`.

That pair's left `Qm` already disagrees. The first cell is `Qm(1, 14)`: pre-pass `-0.826539125`, scan `-0.82653883`, relative `2.95e-7`. The scan fills `Qm` with pseudoknot branches; the nested pre-pass does not. Step 4c has not been run. `V` already reads `PrMbar`, and Step 4a stands.

## `--sweep` wired for growth

`--sweep` (default off) builds the gamma chain before the scan from the pre-pass `Qm`, twice: right gap unlimited, and right gap at most the loop cap. `rec_root` then reads that table. The value at `q'` is `gamma(p', q') + C(p'-1)` at the moment the root is first built, kept to the current node cap by the same `g + u_L` key. The pair-to-pair walk is not run. Flags off do not enter this path.

5S flags off: `logZ 96.343443758865405`.

| run | logZ | gamma entries | chain | scan steps |
|---|---|---|---|---|
| 5S `--sweep` | `96.343443758865632` | 1,390,280 | 7.66 s | |
| 40 nt off and `--sweep` | `13.890361241973714` both | 12,398 | 0.022 s | |
| 150 off | `116.26879396844541` | | | 8,992,446 (7.5e6 multiloop) |
| 150 `--sweep` | `116.26878463324961` | 3,948,232 | 32.2 s | 0 |

The 150 nt scan no longer walks. `grow_s` 0.582 → 0.243, roots 6,721, child nodes 0. `logZ` moved by `-9.3e-6`, the pre-pass `Qm` gap above, on top of the preset's 16-entry cap. The 32 s is the uncapped chain (no `--pre-beam`), and it is counted inside `v_s` / `init` because the fill sits in `pr_build_static`. Stored gamma in that run was 3.9e6 entries at 150 nt.

## Gamma cap during the sweep

With `--sweep` and the preset's `--grow-node-trunc 16`, each pair's map is cut to 16 entries before that map is read by a lower pair. The rank is `C_pre(p'-1) + gamma + u_L(q')`. `u_L` is the same arm weight the scan uses, so the fill runs after the `V` arms exist. `--g-key suffix` ranks by `Sbar(q'+1)` instead. A column is marked finished only after that cut, and a later push into it, or a read of it before the cut, throws. Flags off do not copy `C_pre` and do not fill the table.

| run | logZ | gamma entries | chain | scan steps |
|---|---|---|---|---|
| 5S off | `96.343443758865405` | | | |
| 5S `--sweep` | `96.343443758864623` | 72,996 | 1.19 s | 0 |
| 40 nt off and `--sweep` | `13.890361241973714` both | 5,398 | 0.014 s | |
| 150 off | `116.26879396844541` | | | 8,992,446 |
| 150 `--sweep` | `116.26879680558258` | 124,564 | 2.92 s | 0 |

`logZ` moved by `-7.8e-13` at 119 nt and `+2.8e-6` at 150 nt. Stored entries are about 15 per pair per right-gap table, so the cap is binding. The chain is still superlinear in `n`: every pair is kept, and each column still folds every higher pair.

## `--pre-beam` / `--pre-key`

`--pre-beam b_N` (default off) runs the nested pre-pass with beam `b_N` on `Qb`, `Qm`, and `Qm2`. `--pre-key prefix` ranks by `C(i-1) + inside`. `--pre-key inside` ranks by the inside weight. A dropped start is erased and set to log 0, so later reads, including `PrMbar`, miss it. No hyperedge list is stored: the sweep reads the kept `Qm` value, which already sums the inside edges. This is not an `M1` / `M` / `multi` push. When `--pre-beam` is set it chooses the nested beam width (`--nest-beam` is not also applied).

The full-span right-gap bucket is now accumulated once per pair as the sweep moves left. `Qm` cells are visited from the kept-cell index. Flags off do not build that index. With the preset cap, `--sweep` and no pre-beam still matches the previous capped run: 5S `96.343443758864623`, 150 nt `116.26879680558258`. 5S flags off is `96.343443758865405`.

The beam does not bind at 150 nt. The largest column is `Qm` 146, `Qb` 68, both under 200, so `--pre-beam 200` repeats the unbounded pre-pass.

| run | logZ | Qm cells | max column | pairs | chain |
|---|---|---|---|---|---|
| 150 `--sweep` | `116.26879680558258` | 10,495 | 146 | 4,045 | 2.67 s |
| 150 `--pre-beam 64` | `116.26884098455179` | 7,210 | 64 | 4,041 | 2.06 s |
| 150 `--pre-beam 32` | `115.98521961527207` | 4,146 | 32 | 3,177 | 1.00 s |
| 300 `--sweep` | `247.78897043042861` | 43,486 | 293 | 15,738 | 39.5 s |
| 300 `--pre-beam 200` | `247.78898134099288` | 38,959 | 200 | 15,703 | 34.6 s |

At 300 nt, `Qm` and `Qm2` hit 200. `Qb` peaks at 170, so the pair list barely changes. `logZ` moves by `+1.1e-5` from that beam, on top of the sweep. `--pre-beam 32` at 150 nt drops pairs (3,177 vs 4,045) and moves `logZ` by `-0.284`. The chain is still far from linear: 2.7 s at 150 nt, 39.5 s at 300 nt. The right-gap aggregate is summed from every higher pair before the 16-cap, and that aggregate is what each inner pair copies.

## Right-gap aggregate cut to 16

Before an inner pair copies `gap0`, `gap1`, the long-gap bucket, or the right-`Qm` bucket, each of those buckets is cut to 16 with the same rank as a pair map. The long-gap bucket stays cut, so the next column adds into a map of size at most 16. Flags off are unchanged.

| run | logZ | chain |
|---|---|---|
| 5S off | `96.343443758865405` | |
| 5S `--sweep` | `96.343443758864623` | 0.84 s |
| 40 nt off and `--sweep` | `13.890361241973714` both | 0.014 s |
| 150 `--sweep` | `116.26879680569769` | 1.83 s |
| 150 `--pre-beam 200` | `116.26879680569769` | 1.83 s |
| 300 `--sweep` | `247.78897191520156` | 16.4 s |
| 300 `--pre-beam 200` | `247.78898069882658` | 17.0 s |

Against the previous sweep, `logZ` is unchanged at 119 nt, moves by `+1.2e-10` at 150 nt, and by `+1.5e-6` at 300 nt. The 300 nt chain dropped from 39.5 s to 16.4 s. `--pre-beam 200` still does not bind at 150 nt. At 300 nt it still only cuts `Qm` / `Qm2`, and the chain stays at 17 s.

## Hyperedge push

The right-gap aggregate is no longer what `--sweep` reads. One right-to-left pass over the nested pre-pass pushes a map along every `Qb` / `Qm` / `Qm2` hyperedge the outside uses. A map sits on a kept `Qm` or `Qm2` state and on a kept pair. For an edge that builds `s` from a child, the child gets `edge weight * inside of the other children * map_s`. Columns run `n-1` down to `0`. At each column the `Qm` and `Qm2` maps are cut to `n_G` (`--grow-node-trunc`, preset 16) before they push, then freed. Pair maps are cut the same way before the close pushes into `Qm2(i+1, j-1)`. Splits iterate kept states at that right end (`Qb_at`, `Qm_at`), so a beam of width `b_N` gives each state at most `b_N` splits. With the beam off, a column is as wide as the exact pre-pass (`Qm` max 146 at 150 nt).

`n_V` is `--v-cap`. It cuts the `V` map by `C_pre(b-1) + mu` before that map pushes. `--v-floor` returns the smallest kept `mu` when the scan asks for a `b` that was cut. `--v-sweep` makes `rfv_root` and `rfv_key` read this map. Arms and `lK` are still built from the existing `V` pass. There is one gamma table: the nested grammar has no separate inner right-gap. Both `rec_root` flags read it. Flags off do not enter the push. 5S flags off on this binary is `96.343443758865405`.

The pre-pass `Qm` differs from the scan `Qm` because the scan's multiloop admits pseudoknot branches and this pass does not. That gap is by design. Pass 2 is pseudoknot-free. At 40 nt the first cell is `Qm(1,14)`, pre `-0.826539125` vs scan `-0.82653883`, relative `2.95e-7`; the largest growth root that uses the scan `Qm` is pair `23.32` at `q'=34`, relative `5.22e-6`. At 150 nt those roots reach relative `3.939e-4` (pair `65.113`, `q'=119`). The push is not adjusted to close that gap.

Step 4c, caps off (`sum_b C_pre(b-1) + mu_beta[b]` against the pre-pass outside of `x.y`):

| n | pairs | max relative | at |
|---|---|---|---|
| 40 | 172 | `2.842e-14` | `18.24` |
| 150 | 2948 | `6.395e-13` | `110.119` |

The same push against the scan's pass-0 `V` (caps off, `--pr-v-cap 0`):

| n | entries | one-sided | max relative | at | mu | pass-0 V |
|---|---|---|---|---|---|---|
| 40 | 1832 | 0 | `2.132e-14` | `11.15` b `4` | | |
| 150 | 170414 | 0 | `5.644e-05` | `64.70` b `30` | `51.9087648` | `51.9087084` |

40 nt is inside `1e-12`. 150 nt is not. The first 150 nt entry past `1e-12` is pair `32.42` at `b=0` (`79.7064292` both, difference past the printed digits). 4c matches, so the map is the nested outside. The miss against pass-0 `V` is the pseudoknot-free pre-pass, the same design gap as the `Qm` comparison. `mu` is not replaced by the `pr_steps` chain.

States that carry a map, and pushes, with the beam off:

| n | gamma states | V states | pushes | Qm max | chain | logZ |
|---|---|---|---|---|---|---|
| 40, caps off | 762 | 762 | 10,489 | 33 | 0.007 s | `13.890314599382474` |
| 40, flags off, same trunc | | | | | | `13.890361241973714` |
| 5S, preset cap 16, `--sweep` | 13,133 | 0 | 545,065 | 115 | 0.296 s | `96.343443730096979` |
| 5S, preset cap 16, `--sweep --v-sweep` | 13,133 | 13,133 | 545,065 | 115 | 0.614 s | `96.343443730096979` |
| 150, caps off | 20,937 | 20,937 | 1,091,530 | 146 | 10.5 s | `106.2904257426754` |

The 40 nt flags-off row uses `--grow-node-trunc 0 --grow-ml-c 0 --grow-r-from-v off`, so its `logZ` is the right baseline for that sweep row. The difference is `-4.66e-5`. The 5S difference from flags off `96.343443758865405` is `-2.88e-8`. Adding `--v-sweep` does not move the 5S `logZ` further. The 150 nt caps-off `logZ` is not comparable to the preset flags-off value, because the node cap is off.

