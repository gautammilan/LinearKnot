# Performance pass, 2026-10-31 (branch `perf/merged`)

Five fixes were profiled on branch `LinearPass2` at `764e8ce` (the version 1.7 code, default
G mode lazy). Four are in and one was dropped after measurement. Nothing here changes a
result: every step was required to reproduce the base binary's stdout bit for bit before it
was committed, and the merged binary passes the same check.

## Protocol

- Base binary: `bin/dp03_base`, built from `cpp/main.cpp` at `764e8ce` with
  `g++ -O3 -std=c++17 -ffp-contract=off`, md5 `b004bc3c6611fc44186db8c3c8f34933` (the same
  md5 as the version 1.7 `CURRENT` binary in the notebook).
- Four cases, run by `/tmp/perf5/runcases.sh BIN TAG`:
  - b1: RP150 1-150, `--g-sweep lazy --bG-node 32 --bG-node-key ul`;
  - b2: RP150 1-150, `--g-sweep states`;
  - b3: Ther 1-349, `--g-sweep lazy --bG-node 32 --bG-node-key ul`;
  - b4: `tests/identity_check.sh` (uncapped `bV`/`bG`; V pairs, G states, G lazy, all
    `<= 1e-9` at the paper preset).
  Common flags: `--preset accurate --material rna1999 --dangles some --param-file
  dp03/params/rna1999.dG`.
- Comparison: `/tmp/perf5/cmp.py base TAG CASE` compares every stdout line except `sec`
  (and the `[fine]`/`[stage]` stderr lines are not part of the comparison). A step is
  committed only on `RESULT PASS` for b1 to b3 plus three `ok (<= 1e-9)` lines for b4.
- CPU is `User time` from `/usr/bin/time -v`, on the interactive node cn-v-1, one process at
  a time (the serial runs below), same flags, same sequences.

## Result

CPU is `User time` in seconds, serial on cn-v-1 (interactive node), same flags, same sequences.
Base and merged were run interleaved, two runs each for RP150 (below: the better of the two),
one run for Ther.

| case | base | merged | delta | peak RSS |
| --- | --- | --- | --- | --- |
| RP150 1-150, lazy | 27.46 | 24.00 | -12.6 % | 452 -> 411 MB |
| RP150 1-150, states | 38.12 | 34.95 | -8.3 % | 463 -> 436 MB |
| Ther 1-349, lazy | 168.95 | 152.05 | -10.0 % | 1566 -> 1464 MB |

Merged binary: `bin/dp03_merged`, md5 `6222c11e886f3b037854192354045eef`, from `perf/merged`
at `6e8e3dc`. All four harness cases bit-identical (b4: V pairs, G states, G lazy all
`ok (<= 1e-9)`).

Which step carries it, each against the same base binary, one run each, nothing else on the
node:

| case | base | 2 reads | 3 maps | merged |
| --- | --- | --- | --- | --- |
| RP150 lazy | 27.46 | 27.50 | 24.77 | 24.00 |
| RP150 states | 38.12 | 38.12 | 35.44 | 34.95 |
| Ther lazy | 168.95 | 169.30 | (not re-run) | 152.05 |

The maps carry the result on their own (within a couple of percent of the merged binary) and the
reads step measures zero: `mu_at` did one pair lookup and one entry lookup where the old code
did two of each, on maps that are small and hot (4045 pairs, 170,414 entries at n=150), and the
saving does not show up. It stays because it is shorter and exact, not because it pays. Whether
GCC had already commoned the repeated inlined lookups is not verified from the assembly.

Steps 1 and 4 were measured in the earlier single-run series, against its own base
(29.59 / 42.55 / 185.20 s):

| case | timers | log-space skip |
| --- | --- | --- |
| RP150 lazy | 28.41 | 29.35 |
| RP150 states | 41.73 | 42.69 |
| Ther lazy | 180.39 | 185.90 |

Both are neutral: the timers step is instrumentation (its value is the `[fine]` split), and the
skip is an exact simplification that does not fire often enough to matter.

An earlier series also measured the merged binary at 23.98 / 35.05 / 152.72 s against the same
29.59 / 42.55 / 185.20 s base, i.e. -19.0 / -17.6 / -17.5 %. That series ran with a compiler
running in the background, which inflates base and step alike; the percentages land close to the
interleaved ones above but the absolute seconds are not comparable.




## 1. Timers: attribute pass 3 before changing anything (`3827ccd`)

- `[stage]` keeps its pass-1/pass-2/scan CPU split and gains a fourth line
  `[stage] outside` (the pass-3 `outside()` plus `bpp()`, otherwise unattributed), timed by
  `sec_outside`/`sec_bpp` in `scan_driver.hpp`.
- Every sub-stage clock (`SecTimer`, `PrGrowTimer`, `PrPartnerTimer`, the K-read charge, the
  `rec_key` charge) is now behind `-DLK_FINE_TIMERS`: compiled out by default, so the default
  binary keeps zeros in those fields and the `[fine]` line is only printed by an instrumented
  build. Build the profiling binary with
  `make BIN=bin/dp03_fine CXXFLAGS='-O3 -std=c++17 -ffp-contract=off -DLK_FINE_TIMERS'`.
- Bug fixed on the way: `PrPartnerTimer` had no destructor, so its seconds field was 0 in
  every report.
- Measured: neutral to slightly positive (RP150 lazy 29.59 -> 28.41 s, RP150 states 42.55 ->
  41.73 s, Ther lazy 185.20 -> 180.39 s), which is inside the noise of a single run; its
  value is the attribution, not the microseconds.
- What the instrumented binary says for RP150 lazy (`[fine]` line, seconds of CPU at n=150):
  pass 2 growth 1.39, pass 2 select 0.00, pass 2 static 0.26, pass 3 nested tables (qb 0.13,
  qm 0.06, branch 0.04, qz 0.11, c 0.06, init 0.28), **reads K 5.06**, reads key 0.67,
  outside 0.00, bpp 0.03. Pass 2 is 1.7 s of 29.6 s and the named pass-3 sub-stages are
  0.7 s, so the two places to look next were the read paths and the pass-3 state tables,
  which is what steps 2 and 3 did.
- Open: the `[stage]` lines account for 21.9 s of the 29.6 s user time at RP150 lazy. About
  7.7 s is not attributed to any stage yet (not the sweep, not `outside()`, not `bpp()`).

## 2. Reads: one `prMu` lookup per read instead of three (`e80753f`)

- `mu_at(x, y, xs)` returned `prMu[key2(x,y)][key2(xs,0)]` and first called `v_map_kind`,
  which itself looked up `prMu[key2(x, y)]` and then `[key2(xs, 0)]` to decide the kind. That is
  four hash lookups on a hit (two of them in `v_map_kind`, one for the outer map again, one
  for the inner), on the hottest reads path: `rfv_key` over the closings of a state calls it
  once per closing.
- Measured call count: **65,135,618 `mu_at` calls for one RP150 1-150 run** (`prMu` holds
  4045 pairs and 170,414 entries at that size, 15.6 closings per `rfv_key` call over 4.17
  million U-key lookups).
- Fix: `mu_at` does the single `prMu.find` itself and derives the kind from that iterator
  (`kind 3` when the pair is absent, `0` when the entry exists and is `> HALF_NEG`, else `1`
  when the inner map is at the `--bV` cap, `2` otherwise). The `vm_i`/`vm_ii`/`vm_iii`
  counters see exactly the same increments on the same calls, because the same conditions are
  tested in the same order; `v_map_kind` stays for its other callers (`k_map_kind`).
- Effect: **neutral** (27.50, 38.12, 169.30 against the base 27.46, 38.12, 168.95). The maps it
  reads are small and hot (4045 pairs, 170,414 entries at n = 150), so the two saved lookups per
  call do not show up in the total; whether GCC had already commoned the repeated inlined
  lookups was not checked in the assembly. The change stays because it is shorter and exact
  (same counters, same stdout), not because it pays. Bit-identical on b1 to b4.
- Note: this path is on by default, because the `accurate` preset sets `--v-sweep on`.

## 3. Hot index maps: `ankerl::unordered_dense` behind `FastMap` (`e6a8e2e`)

- The pass-3 state tables look up integers in chained maps: `Tab4::Grp::idx`
  (`unordered_map<uint64_t,int>`, one per vertex), `Tab4::vidx` (`unordered_map<int,int>`),
  `SMap::idx` (`scan_inside.hpp`: key -> slot in the keys/states vectors) and `lArm`
  (`unordered_map<uint64_t,int>`, `larm()` is called per seed lookup).
- Vendored `ankerl::unordered_dense` 4.5.0 (MIT, single self-contained header) as
  `cpp/third_party/ankerl_unordered_dense.h`, with `cpp/third_party/README.md` recording the
  URL, version, license and sha256; `cpp/common/fast_map.hpp` includes it once and defines
  `template <class K, class V> using FastMap`. The swap is a declaration change per map, plus
  the local map in `filter_group` (`rank_prune.hpp`) whose type must match the one it is
  swapped into.
- Why it cannot change a result: all four maps are lookup-only index maps. Their iteration
  order never reaches a floating-point sum (`SMap` keeps its own `keys`/`v` vectors, which do
  keep insertion order).
- Effect: the whole of it; see the table (RP150 lazy 27.46 -> 24.77, states 38.12 -> 35.44 with
  nothing else on the node; the Ther maps-only run was not repeated after the first series).
  Bit-identical on b1 to b4.

## 4. Log-space skip in `log_add` (`67e46db`)

- `log_add(a, b)` with `a >= b` needs no `exp`/`log1p` when `b` is far below `a`: for
  `d = b - a < (e_a - 55) ln 2`, with `a` in `[2^e, 2^(e+1))`, `exp(d) < 2^(e-55) = ulp(a)/8`
  and `log1p(x) < 1.01 x` puts the added term below `ulp(a)/4`, inside the half-spacing at
  which `a + t` would round away from `a`, so `log(exp(a) + exp(b)) = a + log1p(exp(d))` returns
  the same double as before. The exponent `e_a` is read from the bits of `a`; zero, subnormal,
  infinite and NaN inputs and `NEG_INF` stay on the old path.
- Bit-identical on b1 to b4 (b4 exact).
- Measured: **neutral** (RP150 lazy 29.35 against 29.59 s, RP150 states 42.69 against 42.55 s,
  Ther lazy 185.90 against 185.20 s, all inside the noise of a single run). The terms that
  reach `log_add` are mostly within the band where the skip cannot fire, so the win is not
  there to take; it is kept because it is exact and free, not because it was measured to pay.

## 5. Batched `lse` behind `--lse-batch`: not applied

The fifth item on the list was to merge the per-term `log_add` chains of the hot accumulators
(`Tab4::acc`, the ledge/lacc loop, the `sw_hyper_push` V-map merge) into one `logsumexp` over a
drained buffer, behind a `--lse-batch` flag with a `1e-9` acceptance bound and a
max-minus-`54 ln 2` band drop. It was dropped, for two reasons that the measurements above
give:

- The step that removes exponentials from `log_add` (4) is neutral on all three timed cases, so
  the `exp`/`log` pairs that batching would remove are not a measurable share of the run. A
  batch rewrite removes at most a superset of the same work.
- Batching changes the summation order, so it cannot be bit-identical: it needs a tolerance
  (`<= 1e-9` relative) and therefore a second, weaker acceptance test for every future change,
  which is a bad trade for a fix whose upside the profile does not show.

Where the time actually sits, from the attribution build and the step measurements: the pass-3
state tables (step 3, the only measured win), the ~7.7 s per RP150-lazy run that no stage claims
yet (open item in section 1), and the reads, whose 65 million `mu_at` calls per RP150 run are
cheap individually (step 2 measured neutral, so the lookups themselves are not a lever). If a
transcendental-heavy kernel is one day found, batching is the way to attack it; nothing measured
here points at one, so the first move is to attribute that 7.7 s.


## 6. Reproduce

```
# base and merged binaries
cd LinearNupack-LinearPass2            # at 764e8ce, clean
make BIN=bin/dp03_base
git checkout -B perf/merged 764e8ce     # merges 3827ccd, e80753f, e6a8e2e, 67e46db
make BIN=bin/dp03_merged

# bit-identity: four cases per binary, stdout compared line by line except `sec`
/tmp/perf5/runcases.sh bin/dp03_base   base
/tmp/perf5/runcases.sh bin/dp03_merged merged
python3 /tmp/perf5/cmp.py base merged b1     # also b2, b3; b4 is tests/identity_check.sh

# attribution build (the `[fine]` line)
make BIN=bin/dp03_fine CXXFLAGS='-O3 -std=c++17 -ffp-contract=off -DLK_FINE_TIMERS'
```

The harness scripts, sequences and per-run outputs of this pass are in `/tmp/perf5/`
(`runcases.sh`, `cmp.py`, `seq/`, `out/`, `timing.txt`).

