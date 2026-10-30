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
