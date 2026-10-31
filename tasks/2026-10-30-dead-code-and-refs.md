# 2026-10-30 Dead code, exact-reference provenance, paper build tag, NL4-3

## 1. Dead code (commit on LinearPass2)

Deleted all code that none of the coverage runs executes; kept the identity check (`--out-id`)
and the pair traces (`--g-trace`); deleted the other dumps and diagnostics; removed the unused
options and the three preset options whose code never runs (`--stem-est-chain`,
`--stem-est-profile`, `--grow-rec`). Record of what went, the runs, and the checks:
`docs/UNUSED_CODE.md`. 182 options removed (the 181 listed except `--g-trace` and `--help`, plus
the three preset options); 57 remain. 50 files under `cpp/`, about 10,950 lines removed; six files
deleted (`diag/dumps.hpp`, `diag/g_eq.hpp`, `diag/stem_diag.hpp`, `diag/sweep_checks.hpp`,
`diag/v_diag.hpp`, `common/gate_index.hpp`).

Coverage runs (gcov, union over runs): Ther, Deha, Pseu, RP486, HXB2 1-3000 at the paper preset;
the exact command (RNase P D. radiodurans 1-100, with and without bpp; with bpp it segfaults in
outside, as paper-v1 does; not fixed here); RP150 with `--out-id` (uncapped b_V, b_G) and with
`--g-trace`.

Checks on the final tree (release build `-O3 -std=c++17 -ffp-contract=off`, binary
`/tmp/lp2/chk/dp03_final3`, md5 018486b58c2469133a1b53d12d23f355; reference outputs
`/tmp/new_exactness/growth_read/step2/`):

| sequence | logZ (final tree) | stdout except `sec` vs paper-v1 |
|---|---:|---|
| Ther | 327.39766236390301 | identical (1193 lines) |
| Deha | 233.58001232046206 | identical (2969 lines) |
| Pseu | 215.53327699794767 | identical (5307 lines) |
| RP150 | 116.28132370133829 | identical |
| RP486 | 440.33933303544205 | identical (3550 lines) |

- Identity check (`make check` / `tests/identity_check.sh`): V 2948 pairs max_rel 7.816e-13, G
  states 7.390e-13, as before.
- Exact command (no bpp), RP150 `--out-id`, RP150 `--g-trace`: stdout identical to the coverage
  runs of the original tree; `[out-id]` lines and the 41 `[g-trace]` lines identical.
- Forward check: 10,991 original lines deleted or replaced, 919 of them executed; all are
  option-parser branches, folded constants, unwrapped `if` lines, printf literals or moved code.
- Reverse check: coverage build of the final tree on the small runs; no executed line comes from
  an original line that never ran. It caught one bug before the commit: the `dp.outside()` of
  the deleted `DUMP_P2_ALL` block had been left in `run_report.hpp` and crashed the exact
  command; removed.

**HXB2 1-3000 coverage.** Run `/tmp/lp2/cov2/runs/HXB2_3000` (coverage build of the original
tree, started 19:18:54, guarded by `hxb2_guard.sh`, which kills it if cgroup headroom falls
below 12 GB). It had not finished at commit time. `/tmp/lp2/cov2/hxb2_post.sh` waits for it,
adds it to `union.json`, reruns the forward check and writes `/tmp/lp2/cov2/hxb2_post.txt`:
the deleted lines that only HXB2 executed. If that list is not empty, those lines go back in a
follow-up commit.

**Incident: hxb2_3000 (step 3) killed.** driver3 killed its own job `c/hxb2_3000` at 19:18:23
(cgroup headroom 15.0 GB, its floor 15 GB), after pass 2 had finished (6585 s, peak 40.8 GB). It
does not requeue. My intermediate RP486 check run (3.8 GB, 19:14 to 19:32) probably pushed the
headroom under the floor. I did not relaunch it.

## 2. Exact references: provenance

Full table with evidence: `/tmp/lp2/exactref/provenance.md`. Summary:

| reference | binary | commit | full command | ML leading-run cap 30 |
|---|---|---|---|---|
| ArchiveII, 16 sequences (RNaseP M. maripaludus + 15 tmRNA, incl. Ther, Deha, Pseu) | `Algorithm/nupack3.0.6-dp03cap/bin/pairs`, md5 ceddfd2cf8c3c2c2a186e41ab1fbfeb1 (NUPACK 3.0.6 + local patch `-DDP03_MAXLOOP=30`) | none: the NUPACK tree is not under git | `NUPACKHOME=Algorithm/nupack3.0.6-dp03cap pairs -material rna1999 -dangles some -pseudo -T 37 -validate <id>`; logZ = -dG / 0.6163207755 | off (the patch caps interior / bulge loops only; the ML-lead builds are from 10-25, after the runs of 10-02 to 10-08) |
| RNase P D. radiodurans 1-150 | `LinearNupack/cpp/dp03_cpp_stem`, md5 1096a4999fee4280596500a0b5104e41 | none: built from untracked `cpp/dp03.cpp` before the git history starts (10-23) | `dp03_cpp_stem --loop-cap 30 --max-loop 30 --wc-ends --rt 0.6163207755 --material rna1999 --dangles some --gap-beam step --cutoff 1e-5 --beam inf ... < rnasep150.fa` | off (no such option in that binary; a rerun with the cap off matches to 8.6e-11, with the cap on differs by 3.3e-8) |
| RNase P D. radiodurans full length (486 nt) | `LinearNupack/cpp/dp03_cpp` as built on 10-01 (since overwritten; md5 unknown); not `dp03_cpp_stem` | none (same reason) | via `archiveII/run_step.py --beams 0`: `dp03_cpp --beam 0 --gap-beam step --loop-cap 30 --max-loop 30 --rt 0.6163207755 --cutoff 1e-05 --wc-ends --material rna1999 --dangles some` | off (no ML-lead code in `dp03.cpp` or any build of that period) |

Notes:
- The `--beam inf` command stored with the 486-nt reference was typed by hand when the file was
  normalised; the actual run used `dp03_cpp --beam 0`.
- The ArchiveII comparisons ran LinearNupack with `--ml-lead-cap 30` (the accurate preset)
  against references without the cap.

**Tag `exact-ref`: not made.** None of these binaries maps to a commit. The NUPACK build is not
under git, and the DP03 exact engine (`cpp/dp03.cpp`) was never tracked; the LinearNupack history
starts on 10-23, after all reference runs (10-01 to 10-08). Tagging any existing commit would
claim a provenance that does not exist. To make the tag meaningful, the source of the reference
binary would have to be committed (for example `cpp/dp03.cpp` as it is now, with a note that it
is not proven identical to the 10-01 build), or a new exact reference produced from a tracked
commit. The current exact mode was not fixed, as asked.

## 3. Paper build tag

`cpp/gate_index.hpp` and `cpp/stem_est.hpp` committed on top of 7c4148d (b0987a9, nothing else);
a clean checkout builds `cpp/dp03_l2r_merged.cpp` with md5 48d4b172c657d3825eef86d7a31b5793;
tag `paper-v1-build` = b0987a9. `dp03/params/rna1999.dG` tracked on LinearPass2 (6652db2).

## 4. NL4-3

Not started: `b/nl43_b4_1000` and `b/nl43_b4_3000` are queued in step-3 driver3 (since 15:06);
the 15:04 launch was killed at 15:05 during the OOM cleanup, before pass 1 finished. No NL4-3
process now: 0 CPU, no RSS. Not relaunched (N2 also not relaunched).
