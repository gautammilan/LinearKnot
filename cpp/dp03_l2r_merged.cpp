// C++ port of dp03/engine.py + dp03/energy.py (beamed DP03 inside/outside).
//
// Mirrors the Python engine method-for-method.  Python dicts iterate in
// insertion order and the beam uses a stable sort, so every "at" table and
// 4-index group here keeps insertion order too; that keeps the log-space sums
// in the same order and the results equal to the Python ones up to rounding.
//
//   dp03_cpp [--beam B] [--loop-cap L] [--max-loop M] [--no-gu] [--rt RT]
//            [--no-single-stems] [--wc-ends]
//            [--score qz|prefix] [--no-pk] [--cutoff C] [--no-bpp] [--gap-beam group|step]
//            [--params dp03|dp09|dp09st|dp09ml|dp09full]
//            [--material dp03|rna1999|vienna2004 [--param-file PATH] [--dangles none|some|all|d2]]
//            [--arm-score on|off --access-file PATH [--max-gap-scan N] [--max-arm-len W]]
//            [--pair-bonus X] [--restrict-pairs]
//            [--gate off|on [--gate-tau T] [--gate-span M] [--gate-delta D]
//             [--gate-role stem1|both] [--gate-tables qg|qgl|qg,qgl] [--gate-reserve K]]
//            [--gate-stats] [--gate-diag PATH [--gate-diag-taus T1,...] [--gate-long-span L]]
//            [--beam-table T=B ...] [--dump-states PATH [--dump-min-w W]]
//            [--dump-pruned PATH [--dump-exact PATH [--dump-pruned-min-w X]]]
//            [--remove-state T,a,b,c,d ...]
//            < sequences
//
// --params: multiloop and pseudoknot penalties (PkParams below); the nested-loop
//   energies are the same in all sets.
// --material: nested-loop energies; dp03 = the compact table below (default), rna1999 =
//   NUPACK 3's Turner 1999 file (default path ../dp03/params/rna1999.dG next to this
//   binary) with terminal AU/GU penalties; --dangles some|all adds NUPACK's dangles
//   (ExplDangle: d5 + d3 per unpaired run, "some" = the smaller one for a single base);
//   vienna2004 = LinearPartition-V's Turner 2004 model (../dp03/params/vienna2004.par),
//   --loop-cap 70 by default, --dangles none|d2. Pseudoknot penalties still come from
//   --params (dp03 or dp09); the nested loops use the Turner 2004 tables.
// --beam inf (or 0) = infinite beam: nothing is pruned.
// --arm-score on: Qg/Qgl pruning score gets the LinearAccess downstream-arm term
//   (engine.py docstring); per-sequence stats go to stderr as "[arm] ...".
//   --max-arm-len W: longest arm window, 0 (default) = sequence length.
// --gap-beam group|step: 4-index tables keep `beam` states per (v, j) group (default)
//   or per step j over all groups.
//   Build with -DDP03_DEBUG (make debug) to assert kappa >= 1 and G >= F(gap).
// --pair-bonus X: add X (log units) for every pair formed (structure-weight evaluation only;
//   not in the Python engine).
// --restrict-pairs: each sequence line is followed by a line of 1-based pairs "i j i j ...";
//   only those pairs may form.  With a large --pair-bonus, logZ - X * #pairs is the log
//   Boltzmann weight of that structure (archiveII/ref_prob.py).
// --gate on: eligibility gate on the Qg/Qgl candidates (see GateCfg below), from the
//   P lines of --access-file; independent of --arm-score (best used with it off).
//   Defaults: tau 0.01, span 0 (= n, no cap), delta 1e-3, role both, tables qg,qgl,
//   reserve 0.  --gate on --gate-tau 0 --gate-reserve 0 gives the same output as
//   --gate off.  With tau > 0 states are removed even at infinite beam, so
//   infinite-beam exactness holds only for tau = 0.  "[gate] ..." stats go to stderr
//   with --gate on or --gate-stats (candidate counts also with the gate off).
// --gate-diag PATH: after outside, one JSON line per sequence with every Qg/Qgl state's
//   gate values against its posterior weight (gate_recall.py; run at --beam inf).
// --beam-table T=B (repeatable): beam B (or inf) for table T only, T in Qb, Qm, Qm2, Qmp,
//   PKM, Qp, Qg, Qgls, Qgl, Qgr, Qgrs; the others keep --beam.  Qm2 is dropped together
//   with Qm while Qm2's beam is finite, then held to its own beam.
// --b2 N: beam of the scan's two-index columns Qb, Qm, Qm2, Qp, Qmp, and of Q^z / Q^z_pre
//   (writes those --beam-table entries; a later --beam-table wins for that table). Unset: 200.
// --b4 N: four-index scan beam (was --beam). Unset: 1000 under --pk-l2r on defaults, else 0.
// --bV N: V left-end cap per map (was --v-cap / --pr-v-cap). Sets both p2_v_cap and pr_vcap.
// --bG N: growth node trunc K (was --grow-node-trunc). Unset: 16 when grow-rec defaults apply.
// --g-cut topq|collapse|entries: sw_fill node cap (default topq). topq = all p' for top-bG q';
//   collapse = Ghat of size bG; entries = top-bG raw (p',q') by C_pre*G*S_pre.
// --g-sweep direct|states (default direct): how G is built for growth reads.
//   direct = sw_fill (pair→pair via PrMbar); states = sw_hyper_push do_g on Qm/Qm2/pair (E1–E9),
//   cut every state to bG entries by C_pre*G*S_pre; scan still applies Chat at read.
// --qzbar-old: build StemEst Qzbar as the Oct-8/incperf path does — dense full nested Qz→prQz,
//   no beamed sparse Qz_pre (qz_on). Used by H_ex, K_L, pr_sf_build, and every other Qzbar read.
// --b1 N: --nest-beam N, the nested pre-pass and the stem-estimate run at every length.
//   Unset: 200, independent of --b4. 0 keeps every state.
// --pbeam-ratio a:b:c: P1,P2,P3 beams = (a,b,c) times --b4. Unset: pbeam = (b2, b4, b4).
//   An explicit --l2r-pbeam or --pbeam-ratio on the command line wins.
// Removed (error if passed): --beam, --v-cap, --pr-v-cap, --grow-node-trunc.
// --dump-states PATH: after outside, every stored state's posterior weight (format at
//   DP03::dump_states); run at --beam inf for exact weights.  --dump-min-w W: write only
//   states with w >= W (the rest are summed per table).
// --dump-pruned PATH: "D id table src j v a b c d rank cands score w" for every candidate
//   a pruning event drops (src = the table whose ranking dropped it; v = group in group
//   mode, else -1; rank 0 = best of the event's cands), then per table "P id table events E
//   dropped D matched M w W pos W0 ... W9".  --dump-exact PATH (a --dump-states file):
//   write only the dropped states found there, with their exact w (else w = nan); M, W
//   count them and Wk is their w with floor(10 j / n) = k.  --dump-pruned-min-w X: D lines
//   only for w >= X (the P sums keep all).  Not with --gate on.  With --stem-est rank the
//   stem ranking's drops are logged (rank among the feasible states); the exact removal of
//   infeasible states is not.  --dump-kept PATH (needs --dump-pruned): "K ..." lines (fields
//   of D; w nan if not in the exact file) for the kept states of every stem / helper ranking.
// --dump-scores PATH: exact step-mode run (every table at inf): every candidate's rank under
//   the default score and under the --stem-est rank score with role all and loop
//   --stem-est-loop, with its exact w (>= --dump-min-w) and stem role class (format at
//   DP03::dump_scores).
// --stem-est-qg-r mirror|partner (default mirror; partner needs role all): Qg's stem-2 score
//   from the mirrored index H_L, or from the partner tables: log Bg(d, u, v) + alpha +
//   log S[q+1] (DP03::stem_Bg).  --stem-est-consts as-is|consistent (default as-is): the
//   pseudoknot factors each role adds ("[stem-consts]" stderr line when either flag is set).
//   Only the ranking changes; the removal of infeasible states does not.
// --stem-bg-check: with --stem-est rank at --beam inf, compare stem_Bg with a brute-force sum
//   ("[bg-check]" stderr line).
// --stem-est-qg-ext off|stack|bulge1 (default off; needs role all and qg-r partner): Qg's
//   score also sums sigma_ext, the stem-1 and stem-2 scores of the state grown outward by
//   case-(b) wraps (DP03::stem_ext_score); at most --stem-est-ext-max pairs (default 30).
//   --stem-est-qg-combine sum|split (default sum): split keeps up to ceil(b/3) states each by
//   sigma_L, sigma_R and sigma_ext, then fills to b by the sum (DP03::split_mask).
//   "[stem-ext]" stderr line.  --stem-ext-check: compare sigma_ext with a path-by-path
//   recomputation at every Qg ranking ("[ext-check]" stderr line).
//   --stem-est-qg-ext loop: growth steps with up to --stem-est-grow-loop Lg (default 3)
//   unpaired bases per side (case-(b) wraps); per layer, nodes with log G < log(c) + the layer
//   max are dropped (--stem-est-grow-cut c, default 1e-8; 0 = no cut), and at most
//   --stem-est-grow-nodes N (default 2000; 0 = no cap) nodes are kept per outer pair
//   (DP03::stem_grow_nodes).  --stem-est-helpers off|grow|nested (default off; needs
//   --stem-est-qg-ext): Qgls / Qgrs are ranked by their closings into a Qg (cases f / g / h;
//   f over j <= q + 1 + --stem-est-ml-span, default 30) times, for grow, the stop-or-grow
//   outside O of the closed stem (DP03::helper_sigma) or, for nested, the outside weight of
//   the closing pair in the nested-only run (DP03::helper_sigma_nested; one lookup per state).
// --stem-est-profile off|on (default off; needs --stem-est-qg-ext loop, --stem-est-helpers grow
//   and --stem-est-qg-combine sum): the same Qg / helper scores, evaluated from per-key
//   profiles over the final right ends y' (L part) and left ends x' (R part) of the growth
//   nodes (DP03::pair_prof, prof_eval).  --stem-est-profile-eps e (default 0 = off): per key,
//   keep the smallest set of entries holding 1 - e of each profile's weight.
//   --stem-est-ml-tail off|nested (default off): Qgls also adds its closings j > q+1+L_ml with
//   the nested outside of Result 20 (gap-free tail times Qz).  --stem-profile-check: compare
//   every profile score with the Result-19 score ("[profile-check]" stderr line).
// --stem-est-chain off|on (default off; needs --stem-est-profile on, no --stem-est-ml-tail):
//   pair profiles over arm chains of kept pairs (DP03::chain_node): pairs with nested pair
//   probability >= --stem-est-chain-floor f (default 1e-7; the b1 nested run), joined by
//   interior steps (<= --stem-est-chain-loop Lc
//   unpaired per side, default 30, capped by --loop-cap) and, with --stem-est-chain-ml on
//   (default), multiloop steps (f: Qm on the left, <= ml-span unpaired on the right; g: <= cap
//   unpaired on the left, nested Qm on the right; h: both; Qm sides <= --stem-est-chain-span
//   S, default 150).  Per chain node, profile entries below --stem-est-chain-cut c (default
//   1e-8; 0 = none) times its max are dropped.  With --stem-profile-check, scores are
//   compared with every chain path summed separately.
// --stem-est-preset v1: expanded in place (later flags override it) to variant A of Result 27:
//   --stem-est rank --stem-est-role all --stem-est-loop 0 --stem-est-qg-r partner
//   --stem-est-qg-ext loop --stem-est-grow-loop 2 --stem-est-helpers grow --stem-est-profile on
//   --stem-est-chain on --stem-est-chain-loop 2 --stem-est-chain-ml off.  Defaults unchanged.
// --stem-est-bg-reorder off|on (default off; needs --stem-est-profile on): the R part
//   sum_x' P_R(x') Bg(x', u, v) of a profile score as beta2 + at(u, v) + log sum_e
//   Qz[u+1][e-1] M(e, v), with M(e, v) = sum_x' P_R(x') B(x', e, v) built once per key
//   profile and v (DP03::prof_mrow); equal up to rounding.  --stem-bg-reorder-check also
//   computes the original R part of every scored state and the keep decisions it gives
//   ("[bg-reorder-check]" stderr line).
// --stem-est-hr orig|av|phi (default orig): in every L-stem pruning score, replace the original
//   future-partner H_R by the availability-weighted H_R^av built from the nested pair and
//   unpaired probabilities (av), or by H_R with phi^G multiplied into the original stem
//   table once per stem pair and nothing else changed (phi; G = 0 gives H_R bit for bit).
//   --stem-est-hr-gamma G (default 1; > 0 for av, >= 0 for phi) raises each pair availability
//   to G.  Only pruning scores change; exact-removal feasibility still uses the original H_R.
//   --stem-hr-check explicitly enumerates every partner stem and checks H_R^av / H_R'.
// --stem-est-principled off|on (default off; needs --stem-est rank, --stem-est-role all):
//   every pseudoknot table (Qg, Qgl, Qgr, Qgls, Qgrs, Qp) ranked by alpha + log Out_hat,
//   Out_hat following DP03's edges outward from the state: exact scan tables at positions
//   <= q, a nested-only pre-pass beyond q (exact for n <= 1000, else the nested beam,
//   every pair span), outside the whole pseudoknot Ctx(i, j) = e^beta1 C[i-1]
//   Sbar[j+1] + e^beta1m Bml[i][j] (Bml: the pseudoknot as one multiloop branch, from the
//   nested outside; a pseudoknot inside another one's gap is not represented).  Stems grow
//   and future stems are built by DP03's stem steps over the nested pairs with probability
//   >= --pr-floor f (default 1e-7): interior steps with --stem-est-chain-loop's cap, multiloop
//   steps if --stem-est-chain-ml on, multiloop parts of span <= --pr-ml-span S (default 150).
//   Both apply to the growth and to the future stems; --pr-partner-ml off (default on)
//   removes the multiloop steps from the future partner stems (V) only.
//   Qgls closings up to --pr-ml-close L unpaired (default 30; the engine's cap for seeds),
//   profiles cut at --pr-cut c (default 1e-8) times their max, --pr-arms-k k arms per gap by
//   their gap-free future factor (default 64; 0 = all).  Exact removal is unchanged.
//   --stem-pr-check: every ranked state's Out_hat against an enumeration of every DP03 path
//   in long double, and the Bml table cell by cell (n <= 200) ("[pr-check]" stderr line).
//   --pr-grow-select floor|bound (default floor): bound grows over every pair (implies
//   --pr-floor 0) and keeps the endpoints p'.q' of every growth node, root and key profile by
//   iota = w (U_L + U_R), U the gap-free bound of the step on the endpoint's L and R parts
//   relative to the gap's Qz[u+1][v-1] (DP03::pr_U), largest first until 1 - --pr-grow-eps
//   (default 1e-6) of the total, at most --pr-grow-B (default 512; 0 = no cap); nodes when
//   built, roots and keys at the step that evaluates them; replaces --pr-cut on the growth
//   profiles.  --pr-grow-check: every evaluated endpoint's L (all arms) and R against
//   Qz[u+1][v-1] U ("[pr-grow-check]" stderr line).  --pr-score-dump PATH: every principled
//   score as "table q p u v sigma".
// --pr-abl LIST (comma list; offline ablation of sigma_pr's couplings, with --pr-oracle-set
//   score): bml Ctx without its multiloop term; q the partner of Stop_L beyond the state's q,
//   T(q; u, v) / Qz[q+1][q'] for T(q'; u, v); p Stop_R with K(p; u, v) of the state's p for
//   K(p'; u, v); xs every arm of the gap with Qz[u+1][v-1] / Qz[xs][x'] for its two fillers;
//   ks the same bound for K's fillers, Qz[u+1][v-1] / Qz[e][f]; sep (implies bml, q) the
//   separable evaluation at DP03::pr_eval_sep, Stop_R over the top --pr-abl-pk K left ends p'
//   (default 0 = all); sepv (implies bml) the same with Stop_L over the right ends q' of the
//   profile (exact under bml); xb the xs bound on the selected arms only; qr Qz[0][q'] / Qz[0][q] for
//   q's Qz[q+1][q'].
//   --pr-oracle-set score: "A id tab q N W kept100 orc100 kept300 orc300 kept1000 orc1000 rho
//   nmiss wmiss".
// --stem-est-principled-nested off|on (default off): Qb, Qm, Qm2 ranked by alpha + their
//   outside in the nested pre-pass.
// --dump-ext PATH --dump-exact STATES: inside-only step-mode run (every table at inf, use
//   --no-bpp): every Qg candidate's ranks under the Result-18 scores, joined with the exact
//   weights of STATES (format at DP03::dump_ext).
// --remove-state T,a,b,c,d (repeatable): erase that state at its table's pruning point
//   (any beam); log Z then drops by log(1 - w) for an unambiguous grammar (tests the dump).
// stdin: one sequence per line; an optional ">id" line names the next one.
// stdout per sequence: [>id] / logZ / sec / pairs K / K lines "i j p" (1-based).

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <csignal>

// SIGXCPU (soft RLIMIT_CPU) / SIGUSR1: stop at the next step, print the counters so far.
static volatile std::sig_atomic_t g_cpu_stop = 0;
static void on_cpu_stop(int) { g_cpu_stop = 1; }
#include <string>
#include <vector>
static std::string g_v2dump;  // --pr-v2-dump
static bool g_static_only = false;  // --pr-static-only: exit after the V / arm build
static bool g_v2ratio = false;  // --pr-v-ratio2 on: Vr2 in place of V
// --grow-rec on: g_L / g_R from layered marginal vectors (DP03::RN) instead of growth profiles;
// --grow-ml-L / -R: highest multiloop layer m used by the L / R stop (0 none, 1 depth1, 2 all);
// --grow-ml-c C: top C choices per multiloop half-step (0: all); --grow-ml-diag: every layer
// built and each scored state's mass kept by the five (L, R) configurations logged;
// --grow-rec-check: every recursion score against the profile score of the same state.
static bool g_grow_rec = false, g_ml_diag = false, g_rec_check = false;
static int g_ml_L = 2, g_ml_R = 2, g_ml_c = 0;
// --grow-ml-ub on|off (default on): the --grow-ml-c cut of a node's multiloop terms builds only the
// helpers whose tot bound can reach the top C (same terms kept, same order)
static bool g_ml_ub = true;
// --pr-U-nobml: pr_U without the Bml term; g_skip_bml: build_principled leaves Bml empty.
static bool g_U_nobml = false, g_skip_bml = false;
// --grow-trunc K (with --grow-rec on): at root / key evaluation, the g_L entries ranked by
// g_L[q'] u_L(q') and the g_R entries by g_R[p'] u_R(p') (step-q pr_U halves without Bml),
// largest first until 1 - --grow-trunc-eps of the total, at most K; --grow-trunc-mode side:
// per side (K per side), pool: one ranking over both sides (K in total).
static int g_tr_K = 0, g_tr_pool = 0;
static double g_tr_eps = 1e-6;
// --grow-node-trunc K (with --grow-rec on, ml all/all): every growth node / root / helper / key
// keeps one layer-summed list per side, its top K entries by --grow-trunc's ranking at build
// time, built from its terms' lists; the full layer vectors are never stored.
static int g_ntr_K = 0;
static FILE* g_mass_f = nullptr;
// --g-cut topq|collapse|entries: how sw_fill caps each growth node (default topq).
//   topq:     ca8cd5e — top-bG q' by Ghat+ntrUL (or Sbar); retain all (p',q') for those q'.
//   collapse: Ghat[q'] = sum_p' G[p',q'] * C_pre(p'-1), bG entries as key2(0,q').
//   entries:  pass2_push_edges — top-bG raw (p',q') by C_pre(p'-1)*G*S_pre(q'+1).
enum { G_CUT_TOPQ = 0, G_CUT_COLLAPSE = 1, G_CUT_ENTRIES = 2 };
static int g_cut = G_CUT_TOPQ;
enum { G_SWEEP_DIRECT = 0, G_SWEEP_STATES = 1 };
static int g_sweep = G_SWEEP_DIRECT;
// --grow-node-trunc-r KR: K of the R side (0: K; -1: R not truncated, stored dense in R[0], the
// layers summed)
static int g_ntr_KR = 0;
// --grow-r-from-v on (with --grow-node-trunc K --grow-node-trunc-r -1, --pr-v-key xs): the
// U-type scores read V (kept after the build; seed roots through the arm entries) in place of
// the R side of the growth vectors, which is then not built (only its total, the --grow-ml-c rank).
static bool g_rfv = false;
// --miss-sample N PATH: reservoir of up to N lookups per class where the V map exists
// and the entry is absent (kind ii). Classes: V (every mu_at), K and H on P3 and S.
struct MissRec {
  int st, x, y, b;
};
static int g_miss_cap = 0;
static std::string g_miss_path;
static std::vector<MissRec> g_miss[3];
static long long g_miss_seen[3] = {};
static uint32_t g_miss_rng = 1u;
static uint32_t miss_urand() {
  g_miss_rng = g_miss_rng * 1664525u + 1013904223u;
  return g_miss_rng;
}
static void miss_take(int kind, int st, int x, int y, int b) {
  if (g_miss_cap <= 0 || kind < 0 || kind > 2) return;
  long long& seen = g_miss_seen[kind];
  auto& buf = g_miss[kind];
  ++seen;
  const MissRec rec{st, x, y, b};
  if (static_cast<int>(buf.size()) < g_miss_cap) {
    buf.push_back(rec);
    return;
  }
  if (seen <= 0) return;
  const uint64_t j = (static_cast<uint64_t>(miss_urand()) << 32 | miss_urand()) %
                     static_cast<uint64_t>(seen);
  if (j < static_cast<uint64_t>(g_miss_cap)) buf[static_cast<size_t>(j)] = rec;
}
// --v-query IN OUT: after the lazy V build, print each sample's log weight and the
// largest left-end in that map. H lines (b < 0) scan partner pairs (x, y>j).
struct VQ {
  int tag, st, x, y, b;
};
static std::string g_vquery_out;
static std::vector<VQ> g_vquery;
static std::string g_vtrace_path;
static std::string g_mu_dump;
static std::string g_eq_path;
// --grow-node-trunc-z bar: the node ranking's arm weight 1/Qz(xs, x') from the pre-pass Qzbar
// (scan: the scan's Qz, which is -inf for x' past the scan position at the first ranking call)
static bool g_ntr_zbar = false;
// --malloc-trim on: malloc_trim(0) after the large frees
static bool g_mtrim = false;
struct CpuStop {};
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <queue>
#include <sstream>
#include <string>
#include <tuple>
#include <map>
#include <functional>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <unistd.h>
#include <malloc.h>
#include <execinfo.h>
#include <fcntl.h>
#include <sys/time.h>
#include <sys/resource.h>

static void mtrim() {
  if (g_mtrim) malloc_trim(0);
}

#include "gate_index.hpp"
#include "stem_est.hpp"

namespace {

constexpr double NEG_INF = -1e300;
constexpr double HALF_NEG = NEG_INF / 2;

inline double log_add(double a, double b) {
  if (a < b) std::swap(a, b);
  if (a < NEG_INF / 2) return NEG_INF;
  if (a - b > 50.0) return a;
  return a + std::log1p(std::exp(b - a));
}
// --g-mass: log-sum of each R-key growth score, per (kind, seed, p, j). Arm is summed in.
static std::unordered_map<uint64_t, double> g_mass_acc;
static void g_mass_note(int kind, int seed, int p, int j, double s) {
  if (!g_mass_f || s < NEG_INF / 2) return;
  const uint64_t k = (static_cast<uint64_t>(kind) << 56) | (static_cast<uint64_t>(seed) << 48) |
                      (static_cast<uint64_t>(p) << 24) | static_cast<uint64_t>(j & 0xffffff);
  auto it = g_mass_acc.find(k);
  if (it == g_mass_acc.end()) g_mass_acc.emplace(k, s);
  else it->second = log_add(it->second, s);
}
static void g_mass_flush() {
  if (!g_mass_f) return;
  for (const auto& kv : g_mass_acc) {
    const int rec[5] = {static_cast<int>(kv.first >> 56), static_cast<int>((kv.first >> 48) & 0xff),
                         static_cast<int>((kv.first >> 24) & 0xffffff),
                         static_cast<int>(kv.first & 0xffffff), -1};
    const double s = kv.second;
    std::fwrite(rec, sizeof rec, 1, g_mass_f);
    std::fwrite(&s, sizeof s, 1, g_mass_f);
  }
  std::fflush(g_mass_f);
  g_mass_acc.clear();
}

// ------------------------------------------------------------------ energy
double RT = 0.616330318;  // --rt overrides (e.g. NUPACK's kB*T at 37 C)

// Multiloop and pseudoknot penalties (kcal/mol), selected with --params.
//   stem:      per pseudoknot stem (Qgl); NUPACK's beta2, HotKnots' Pb (per band)
//   nest:      per nested pair inside a gap (Qz); beta2 / HotKnots' Pps
//   nest_pk_k: a pseudoknot nested in a gap pays nest_pk_k * nest (+ beta1p)
//   ml_pk_k:   a pseudoknot as a multiloop branch pays ml_pk_k * alpha2 (+ beta1m)
//   st_mult, int_mult: scale the energy of a stack / interior loop inside a pseudoknot
//              stem (Qg (a), (b)); HotKnots' stP when the left side has no unpaired base,
//              else intP
//   a_p, b_p, c_p: multiloop inside a pseudoknot stem (Qg (f)-(h), Qgls, Qgrs);
//              stem_ml gives its branches their own Qm table (Qmp) so b_p, c_p apply there
struct PkParams {
  double alpha1, alpha2, alpha3, beta1, beta1m, beta1p, stem, nest, beta3;
  double nest_pk_k, ml_pk_k;
  double st_mult, int_mult, a_p, b_p, c_p;
  bool stem_ml;
};
// Dirks & Pierce 2003, as in NUPACK 3.
constexpr PkParams DP03_PARAMS{3.4, 0.4, 0.0, 9.6, 15.0, 15.0, 0.1, 0.1, 0.1, 2.0, 2.0,
                               1.0, 1.0, 3.4, 0.4, 0.0, false};
// Andronescu et al. 2010 (HotKnots 2.0 params/parameters_DP09.txt, pkmodelDP a, b, c,
// Ps, Psm, Psp, Pb, Pps, Pup, stP, intP, a_p, b_p, c_p); HotKnots counts a nested
// pseudoknot as one closed region and a pseudoknot branch as one multiloop branch.
// dp09 = penalties only; dp09st adds stP/intP; dp09ml adds a_p/b_p/c_p; dp09full both.
constexpr PkParams DP09_PARAMS{3.39, 0.03, 0.02, -1.38, 10.07, 15.00, 2.46, 0.96, 0.06, 1.0, 1.0,
                               1.0, 1.0, 3.39, 0.03, 0.02, false};
constexpr PkParams DP09ST_PARAMS{3.39, 0.03, 0.02, -1.38, 10.07, 15.00, 2.46, 0.96, 0.06, 1.0, 1.0,
                                 0.89, 0.74, 3.39, 0.03, 0.02, false};
constexpr PkParams DP09ML_PARAMS{3.39, 0.03, 0.02, -1.38, 10.07, 15.00, 2.46, 0.96, 0.06, 1.0, 1.0,
                                 1.0, 1.0, 3.41, 0.56, 0.12, true};
constexpr PkParams DP09FULL_PARAMS{3.39, 0.03, 0.02, -1.38, 10.07, 15.00, 2.46, 0.96, 0.06, 1.0,
                                   1.0, 0.89, 0.74, 3.41, 0.56, 0.12, true};
PkParams PK = DP03_PARAMS;

constexpr int MIN_HAIRPIN = 3;
constexpr int MIN_PAIR = MIN_HAIRPIN + 1;
constexpr int MAX_LOOP = 30;
// Interior-loop size cap (l1 + l2). Python's energy.py hard-codes MAX_LOOP;
// raise it only to compare against uncapped NUPACK.
int max_loop = MAX_LOOP;
bool allow_gu = true;
// One-pair pseudoknot stems (Qg_{i,i,j,j} = 1), stored as Qg[j][j][(i,i)].
bool single_stems = true;
// NUPACK 3: pairs bordering exterior / multi / Qz regions must be Watson-Crick.
bool wc_ends = false;
double pair_bonus = 0.0;
// Multiloop leading unpaired run, every pass: closing left p and first branch i
// satisfy i - p - 1 <= ml_lead_cap - 1 (29 when the flag is 30, p >= i - 30).
// 0 turns the cap off. Interior loops, the trailing run, between-branch runs,
// the exterior, and Q^z stay on their own rules.
int ml_lead_cap = 30;
const std::vector<std::pair<int, int>>* allowed_pairs = nullptr;  // 0-based, nullptr = all

// pair types: AU CG GC GU UA UG
inline int base_code(char c) {
  switch (c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'U': return 3;
    default: return -1;
  }
}

inline int pair_type(int a, int b) {
  static const int T[4][4] = {
      {-1, -1, -1, 0},
      {-1, -1, 1, -1},
      {-1, 2, -1, 3},
      {4, -1, 5, -1},
  };
  if (a < 0 || b < 0) return -1;
  return T[a][b];
}

const double STACK[6][6] = {
    {-0.9, -2.2, -2.1, -0.6, -1.1, -1.4},
    {-2.1, -3.3, -2.4, -1.4, -2.1, -2.1},
    {-2.4, -3.4, -3.3, -1.5, -2.2, -2.5},
    {-1.3, -2.5, -2.1, -0.5, -1.4, 1.3},
    {-1.3, -2.4, -2.1, -1.0, -0.9, -1.3},
    {-1.0, -1.5, -1.4, -0.3, -0.6, -0.5},
};

const double HAIRPIN_INIT[21] = {0,   0,   0,   5.4, 5.6, 5.7, 5.4,
                                 6.0, 5.5, 6.4, 6.5, 6.6, 6.7, 6.8,
                                 6.9, 6.9, 7.0, 7.1, 7.1, 7.2, 7.2};
constexpr int HAIRPIN_K = 20;
const double BULGE_INIT[7] = {0, 3.8, 2.8, 3.2, 3.6, 4.0, 4.4};
constexpr int BULGE_K = 6;
const double INTERIOR_INIT[21] = {0,   0,   0,   0,   1.1, 2.0, 2.0,
                                  2.1, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8,
                                  2.9, 2.9, 3.0, 3.1, 3.1, 3.2, 3.3};
constexpr int INTERIOR_K = 20;

inline double log_weight(double g) { return -g / RT; }

inline double size_continuation(int n, const double* table, int k) {
  if (n <= k) return table[n];
  return table[k] + 1.75 * RT * std::log(static_cast<double>(n) / k);
}

// --material rna1999: NUPACK 3.0.6's Turner 1999 loop energies (dp03/rna1999.py) read from
// --param-file, in NUPACK's table layout (init.c LoadEnergies), plus its terminal AU/GU
// penalty on pairs bordering an exterior loop, multiloop or pseudoknot gap.  Bases in the
// NUPACK formulas are coded A=1, C=2, G=3, U=4.
bool mat_rna1999 = false;
std::string param_file;
// --dangles: 0 none; rna1999: 1 some, 2 all; vienna2004: 3 d2
int dangles_mode = 0;

struct Rna1999Tables {
  double stack[36], loop37[90], asym[4], max_asym;
  double tri[1024] = {}, tetra[4096] = {};  // base-4 index of the loop, A=0
  double mm_hp[96], mm_il[96], dangle[48];
  double alpha[3], at_penalty;
  double sint2[576], sint4[9216], asint1x2[2304];
  double polyc3, polyc_slope, polyc_int, beta[5];
};
Rna1999Tables T99;

[[noreturn]] void param_error(const std::string& msg) {
  std::fprintf(stderr, "%s: %s\n", param_file.c_str(), msg.c_str());
  std::exit(2);
}

std::vector<int> line_ints(const std::string& line) {
  std::vector<int> out;
  size_t k = 0;
  while (k < line.size()) {
    while (k < line.size() && std::isspace(static_cast<unsigned char>(line[k]))) ++k;
    size_t b = k;
    while (k < line.size() && !std::isspace(static_cast<unsigned char>(line[k]))) ++k;
    if (b == k) break;
    std::string t = line.substr(b, k - b);
    size_t s = (t[0] == '-' || t[0] == '+') ? 1 : 0;
    if (s < t.size() && t.find_first_not_of("0123456789", s) == std::string::npos)
      out.push_back(std::atoi(t.c_str()));
  }
  return out;
}

int loop_index(const std::string& loop) {
  int idx = 0;
  for (char c : loop) {
    int b = base_code(c == 'T' ? 'U' : c);
    if (b < 0) param_error("bad loop sequence " + loop);
    idx = idx * 4 + b;
  }
  return idx;
}

void load_rna1999(const std::string& path) {
  FILE* f = std::fopen(path.c_str(), "r");
  if (!f) param_error("cannot open (use --param-file)");
  std::vector<std::vector<std::string>> secs;
  std::string line;
  bool open = false;
  int c;
  while (true) {
    line.clear();
    while ((c = std::fgetc(f)) != EOF && c != '\n') line.push_back(static_cast<char>(c));
    if (c == EOF && line.empty()) break;
    if (!line.empty() && line[0] == '>') {
      open = false;
    } else {
      if (!open) secs.emplace_back();
      open = true;
      secs.back().push_back(line);
    }
    if (c == EOF) break;
  }
  std::fclose(f);
  size_t si = 0;
  auto next = [&]() -> const std::vector<std::string>& {
    if (si >= secs.size()) param_error("too few sections");
    return secs[si++];
  };
  auto rows = [&](const std::vector<std::string>& sec, size_t k, size_t width, double* out) {
    if (sec.size() != k) param_error("bad section");
    for (size_t r = 0; r < k; ++r) {
      std::vector<int> v = line_ints(sec[r]);
      if (v.size() != width) param_error("bad row " + sec[r]);
      for (size_t q = 0; q < width; ++q) out[r * width + q] = v[q] / 100.0;
    }
  };
  Rna1999Tables& T = T99;
  rows(next(), 6, 6, T.stack);
  for (int i = 0; i < 3; ++i) {
    std::vector<int> v = line_ints(next()[0]);
    if (v.size() != 30) param_error("loop rows must list sizes 1..30");
    for (int j = 0; j < 30; ++j) T.loop37[30 * (2 - i) + j] = v[j] / 100.0;
  }
  {
    std::vector<int> v = line_ints(next()[0]);
    for (int k = 0; k < 4; ++k) T.asym[k] = v[k] / 100.0;
    T.max_asym = v[4] / 100.0;
  }
  for (int which = 0; which < 2; ++which) {
    double* bonus = which == 0 ? T.tri : T.tetra;
    size_t len = which == 0 ? 5 : 6;
    for (const std::string& s : next()) {
      size_t b = s.find_first_not_of(" \t");
      size_t e = b == std::string::npos ? b : s.find_first_of(" \t", b);
      std::vector<int> v = e == std::string::npos ? std::vector<int>() : line_ints(s.substr(e));
      if (e == std::string::npos || e - b != len || v.size() != 1)
        param_error("bad loop bonus " + s);
      bonus[loop_index(s.substr(b, e - b))] = v[0] / 100.0;
    }
  }
  rows(next(), 16, 6, T.mm_hp);
  rows(next(), 16, 6, T.mm_il);
  rows(next(), 6, 4, T.dangle);
  rows(next(), 6, 4, T.dangle + 24);
  {
    std::vector<int> v = line_ints(next()[0]);
    for (int k = 0; k < 3; ++k) T.alpha[k] = v[k] / 100.0;
    T.at_penalty = line_ints(next()[0])[0] / 100.0;
  }
  auto blocks = [&](size_t k, auto&& put) {
    const std::vector<std::string>& sec = next();
    if (sec.size() != 5 * k) param_error("bad loop table");
    for (size_t b = 0; b < k; ++b)
      for (int j = 0; j < 4; ++j) {
        std::vector<int> v = line_ints(sec[5 * b + 1 + j]);
        if (v.size() != 4) param_error("bad table row " + sec[5 * b + 1 + j]);
        for (int q = 0; q < 4; ++q) put(static_cast<int>(b), j, q, v[q] / 100.0);
      }
  };
  blocks(36, [&](int i, int j, int k, double x) { T.sint2[i * 16 + j * 4 + k] = x; });
  blocks(36 * 16, [&](int i, int j, int k, double x) {
    T.sint4[1536 * (i / 96) + 256 * ((i % 96) / 16) + 64 * ((i % 16) / 4) + 4 * (i % 4) +
            k * 16 + j] = x;
  });
  blocks(144, [&](int i, int j, int k, double x) {
    T.asint1x2[384 * (i / 24) + 4 * ((i % 24) / 4) + 24 * (i % 4) + 96 * j + k] = x;
  });
  {
    std::vector<int> v = line_ints(next()[0]);
    T.polyc3 = v[0] / 100.0;
    T.polyc_slope = v[1] / 100.0;
    T.polyc_int = v[2] / 100.0;
    v = line_ints(next()[0]);
    for (int k = 0; k < 5; ++k) T.beta[k] = v[k] / 100.0;
  }
}

// --material vienna2004: LinearPartition-V's Turner 2004 model (dp03/vienna2004.py), nested only,
// read from --param-file (dp03/params/vienna2004.par, written by dp03/make_vienna2004.py);
// energies in dcal/mol.  LinearPartition-V's loop limits: interior loops of at most 70 unpaired,
// 69 on the left; at most 69 unpaired before a multiloop's first branch.
bool mat_v04 = false;
constexpr int V04_MAX_LOOP = 70;
constexpr int V04_MAX_LEFT = 69;
constexpr int NO_LIMIT = 1 << 30;

struct V04Tables {
  double lxc = 0.0;
  int ml_intern = 0, ml_closing = 0, ml_base = 0, max_ninio = 0, ninio = 0, au = 0;
  std::vector<int> stack, hairpin, bulge, interior, mmI, mmH, mmM, mm1n, mm23, mmExt, d5, d3,
      int11, int21, int22;
  std::unordered_map<std::string, int> tri, tetra, hexa;
};
V04Tables V04;

void load_v04(const std::string& path) {
  FILE* f = std::fopen(path.c_str(), "r");
  if (!f) param_error("cannot open (use --param-file)");
  int c;
  while ((c = std::fgetc(f)) != EOF && c != '\n') {}
  auto tok = [&]() -> std::string {
    std::string t;
    while ((c = std::fgetc(f)) != EOF && std::isspace(c)) {}
    while (c != EOF && !std::isspace(c)) {
      t.push_back(static_cast<char>(c));
      c = std::fgetc(f);
    }
    return t;
  };
  V04Tables& T = V04;
  std::unordered_map<std::string, int*> scalars = {
      {"ML_intern37", &T.ml_intern}, {"ML_closing37", &T.ml_closing}, {"ML_BASE37", &T.ml_base},
      {"MAX_NINIO", &T.max_ninio},   {"ninio37", &T.ninio},           {"TerminalAU37", &T.au}};
  std::unordered_map<std::string, std::vector<int>*> arrays = {
      {"stack37", &T.stack},       {"hairpin37", &T.hairpin},    {"bulge37", &T.bulge},
      {"internal_loop37", &T.interior}, {"mismatchI37", &T.mmI}, {"mismatchH37", &T.mmH},
      {"mismatchM37", &T.mmM},     {"mismatch1nI37", &T.mm1n},   {"mismatch23I37", &T.mm23},
      {"mismatchExt37", &T.mmExt}, {"dangle5_37", &T.d5},        {"dangle3_37", &T.d3},
      {"int11_37", &T.int11},      {"int21_37", &T.int21},       {"int22_37", &T.int22}};
  bool has_lxc = false;
  size_t n_scalars = 0;
  while (true) {
    std::string kind = tok();
    if (kind.empty()) break;
    std::string name = tok();
    if (kind == "scalar") {
      std::string v = tok();
      if (name == "lxc37") {
        T.lxc = std::strtod(v.c_str(), nullptr);
        has_lxc = true;
      } else {
        auto it = scalars.find(name);
        if (it == scalars.end()) param_error("unknown scalar " + name);
        *it->second = std::atoi(v.c_str());
        ++n_scalars;
      }
    } else if (kind == "array") {
      auto it = arrays.find(name);
      if (it == arrays.end()) param_error("unknown array " + name);
      int nd = std::atoi(tok().c_str());
      size_t size = 1;
      for (int k = 0; k < nd; ++k) size *= static_cast<size_t>(std::atoi(tok().c_str()));
      it->second->resize(size);
      for (size_t q = 0; q < size; ++q) (*it->second)[q] = std::atoi(tok().c_str());
    } else if (kind == "loops") {
      std::unordered_map<std::string, int>& m =
          name == "Triloops" ? T.tri : name == "Tetraloops" ? T.tetra : T.hexa;
      if (name != "Triloops" && name != "Tetraloops" && name != "Hexaloops")
        param_error("unknown loops " + name);
      int cnt = std::atoi(tok().c_str());
      for (int k = 0; k < cnt; ++k) {
        std::string s = tok();
        m[s] = std::atoi(tok().c_str());
      }
    } else {
      param_error("bad token " + kind);
    }
  }
  std::fclose(f);
  if (!has_lxc || n_scalars != scalars.size()) param_error("missing scalars");
  for (const auto& a : arrays)
    if (a.second->empty()) param_error("missing array " + a.first);
}

inline int mshift(int b1, int b2) {  // GetMismatchShift: AU CG GC UA GU UG
  if (b1 == 1) return 0;
  if (b1 == 2) return 1;
  if (b1 == 3) return b2 == 2 ? 2 : 4;
  return b2 == 1 ? 3 : 5;
}

struct Energy {
  std::string seq;
  int n = 0;
  std::vector<int> code;
  std::vector<char> cp;  // n*n, can_pair
  double w_alpha1, w_alpha2, w_alpha3, w_beta1, w_beta1m, w_beta1p, w_stem, w_nest,
      w_beta3, wu, nest_pk_k, ml_pk_k, st_mult, int_mult, w_ap, w_bp, w_cp;
  bool stem_mult, stem_ml;

  explicit Energy(const std::string& s) : seq(s), n(static_cast<int>(s.size())) {
    code.resize(n);
    for (int i = 0; i < n; ++i) code[i] = base_code(seq[i]);
    cp.assign(static_cast<size_t>(n) * n, 0);
    for (int i = 0; i < n; ++i)
      for (int j = i + MIN_PAIR; j < n; ++j) {
        int t = pair_type(code[i], code[j]);
        if (t < 0 || (!allow_gu && (t == 3 || t == 5))) continue;
        cp[static_cast<size_t>(i) * n + j] = 1;
      }
    if (allowed_pairs) {
      std::vector<char> keep(cp.size(), 0);
      for (const auto& p : *allowed_pairs) {
        size_t k = static_cast<size_t>(p.first) * n + p.second;
        if (0 <= p.first && p.first < p.second && p.second < n) keep[k] = cp[k];
      }
      cp.swap(keep);
    }
    w_alpha1 = log_weight(PK.alpha1);
    w_alpha2 = log_weight(PK.alpha2);
    w_alpha3 = log_weight(PK.alpha3);
    w_beta1 = log_weight(PK.beta1);
    w_beta1m = log_weight(PK.beta1m);
    w_beta1p = log_weight(PK.beta1p);
    w_stem = log_weight(PK.stem);
    w_nest = log_weight(PK.nest);
    w_beta3 = log_weight(PK.beta3);
    wu = w_beta3;
    nest_pk_k = PK.nest_pk_k;
    ml_pk_k = PK.ml_pk_k;
    st_mult = PK.st_mult;
    int_mult = PK.int_mult;
    stem_mult = st_mult != 1.0 || int_mult != 1.0;
    w_ap = log_weight(PK.a_p);
    w_bp = log_weight(PK.b_p);
    w_cp = log_weight(PK.c_p);
    stem_ml = PK.stem_ml;
    non_c.resize(n);
    crun.assign(n + 1, 0);
    for (int i = 0; i < n; ++i) {
      non_c[i] = code[i] != 1;
      crun[i + 1] = crun[i] + (code[i] == 1);
    }
    if (mat_rna1999) w_at_pen = log_weight(T99.at_penalty);
    if (mat_v04) {
      w_alpha1 = log_weight(V04.ml_closing / 100.0);
      w_alpha2 = log_weight(V04.ml_intern / 100.0);
      w_alpha3 = log_weight(V04.ml_base / 100.0);
      max_left = V04_MAX_LEFT;
      ml_left = V04_MAX_LEFT;
    }
    dangle = mat_rna1999 && dangles_mode != 0;
    if (dangle) {
      d5w.resize(n);
      d3w.resize(n);
      s1w.resize(n);
      for (int k = 0; k < n; ++k) {
        // Non-ACGU (N): dangle / mismatch terms that read it are 0 (LinearPartition-style).
        double d5 = 0.0, d3 = 0.0;
        if (k > 0 && code[k - 1] >= 0 && code[k] >= 0)
          d5 = T99.dangle[code[k - 1] * 4 + code[k]];
        if (k < n - 1 && code[k + 1] >= 0 && code[k] >= 0)
          d3 = T99.dangle[24 + (3 - code[k + 1]) * 4 + code[k]];
        d5w[k] = log_weight(d5);
        d3w[k] = log_weight(d3);
        s1w[k] = log_weight(dangles_mode == 1 && 0 < k && k < n - 1 ? std::min(d3, d5)
                                                                    : d3 + d5);
      }
    }
  }

  std::vector<char> non_c;
  std::vector<int> crun;  // crun[k] = number of C in seq[0, k)
  double w_at_pen = 0.0;

  int max_left = NO_LIMIT;  // interior loop: most unpaired on the left
  int ml_left = NO_LIMIT;   // multiloop: most unpaired before the first branch

  // Log weight of the terminal AU/GU penalty of pair i·j (0.0 for dp03); for vienna2004 also
  // its -d2 mismatch / dangle.
  double at(int i, int j) const {
    if (mat_v04) return log_weight(branch04(i, j) / 100.0);
    return non_c[i] && non_c[j] ? w_at_pen : 0.0;
  }

  // at() for pair i·j closing a multiloop.
  double at_close(int i, int j) const {
    if (mat_v04) return log_weight(closing04(i, j) / 100.0);
    return at(i, j);
  }

  // Dangles (ene.c ExplDangleRaw): an unpaired run s..t of an exterior loop, multiloop
  // or pseudoknot gap gets d5(s) + d3(t); a single base k gets s1w[k].
  bool dangle = false;
  std::vector<double> d5w, d3w, s1w;
  double run(int s, int t) const {
    if (t < s) return 0.0;
    if (s == t) return s1w[s];
    return d5w[s] + d3w[t];
  }

  bool pairs(int i, int j) const {
    return 0 <= i && i < j && j < n && cp[static_cast<size_t>(i) * n + j];
  }

  bool end_ok(int i, int j) const {
    if (!wc_ends) return true;
    int t = pair_type(code[i], code[j]);
    return t >= 0 && t != 3 && t != 5;
  }

  double hairpin(int i, int j) const {
    if (!pairs(i, j)) return NEG_INF;
    int u = j - i - 1;
    if (u < MIN_HAIRPIN) return NEG_INF;
    if (mat_rna1999) return log_weight(hairpin99(i, j));
    if (mat_v04) return log_weight(hairpin04(i, j) / 100.0);
    return log_weight(size_continuation(u, HAIRPIN_INIT, HAIRPIN_K));
  }

  double interior(int i, int d, int e, int j) const {
    if (!(pairs(i, j) && pairs(d, e))) return NEG_INF;
    if (!(i < d && d < e && e < j)) return NEG_INF;
    int l1 = d - i - 1;
    int l2 = j - e - 1;
    if (l1 + l2 > (mat_v04 ? V04_MAX_LOOP : max_loop) || l1 > max_left) return NEG_INF;
    if (mat_rna1999) return log_weight(interior99(i, d, e, j));
    if (mat_v04) return log_weight(interior04(i, j, d, e) / 100.0);
    if (l1 == 0 && l2 == 0) return stack(i, j);
    if (l1 == 0 || l2 == 0) return bulge(i, d, e, j, l1 + l2);
    return internal(l1, l2);
  }

  // interior() for consecutive pairs of a pseudoknot stem, scaled by st_mult / int_mult.
  double stem_interior(int i, int d, int e, int j) const {
    double w = interior(i, d, e, j);
    if (stem_mult && w > HALF_NEG) w *= (d == i + 1) ? st_mult : int_mult;
    return w;
  }

 private:
  // ---- rna1999 (kcal/mol), ene.c HairpinEnergy / InteriorEnergyFull / asymmetryEfn
  // code[k] is -1 for any non-ACGU base (N). Those bases cannot pair; loop lengths still
  // count them; every dangle / mismatch / special-hairpin term that would read an N is 0.
  bool okb(int k) const { return code[k] >= 0; }
  bool ok_range(int lo, int hi) const {
    for (int k = lo; k <= hi; ++k)
      if (code[k] < 0) return false;
    return true;
  }
  int b(int k) const { return code[k] + 1; }
  static double size_log(int size) { return 1.75 * RT * std::log(size / 30.0); }
  static double helix99(int a, int bb, int c, int d) {
    return T99.stack[mshift(a, bb) * 6 + mshift(c, d)];
  }
  static double interior_mm(int a, int bb, int x, int y) {
    return T99.mm_il[(4 * (x - 1) + (y - 1)) * 6 + mshift(a, bb)];
  }

  double hairpin99(int i, int j) const {
    const Rna1999Tables& T = T99;
    int size = j - i - 1;
    double g;
    if (size <= 30) {
      g = T.loop37[60 + size - 1];
    } else {
      g = T.loop37[60 + 30 - 1];
      g += size_log(size);
    }
    bool poly_c = crun[j] - crun[i + 1] == size;
    if (size == 3) {
      if (b(i) != 2 && b(j) != 2) g += T.at_penalty;
      if (ok_range(i, j)) {
        int idx = 0;
        for (int k = i; k <= j; ++k) idx = idx * 4 + code[k];
        g += T.tri[idx];
      }
      if (poly_c) g += T.polyc3;
    } else if (size == 4) {
      if (ok_range(i, j)) {
        int idx = 0;
        for (int k = i; k <= j; ++k) idx = idx * 4 + code[k];
        g += T.tetra[idx];
      }
      if (okb(i + 1) && okb(j - 1))
        g += T.mm_hp[(4 * (b(i + 1) - 1) + (b(j - 1) - 1)) * 6 + mshift(b(i), b(j))];
      if (poly_c) g += T.polyc_slope * size + T.polyc_int;
    } else {
      if (okb(i + 1) && okb(j - 1))
        g += T.mm_hp[(4 * (b(i + 1) - 1) + (b(j - 1) - 1)) * 6 + mshift(b(i), b(j))];
      if (poly_c) g += T.polyc_slope * size + T.polyc_int;
    }
    return g;
  }

  static double asymmetry99(int l1, int l2, int size) {
    const Rna1999Tables& T = T99;
    int asym = std::abs(l1 - l2);
    double g;
    if (size <= 30) {
      g = T.loop37[size - 1];
    } else {
      g = T.loop37[30 - 1];
      g += size_log(size);
    }
    int k = std::min(4, std::min(l1, l2));
    if (asym * T.asym[k - 1] < T.max_asym)
      g += asym * T.asym[k - 1];
    else
      g += T.max_asym;
    return g;
  }

  double interior99(int i, int d, int e, int j) const {
    const Rna1999Tables& T = T99;
    int l1 = d - i - 1;
    int l2 = j - e - 1;
    int size = l1 + l2;
    if (size == 0) return helix99(b(i), b(j), b(d), b(e));
    if (l1 == 0 || l2 == 0) {
      double g;
      if (size <= 30) {
        g = T.loop37[30 + size - 1];
      } else {
        g = T.loop37[30 + 30 - 1];
        g += size_log(size);
      }
      if (size == 1) {
        g += helix99(b(i), b(j), b(d), b(e));
      } else {
        if (b(i) != 2 && b(j) != 2) g += T.at_penalty;
        if (b(d) != 2 && b(e) != 2) g += T.at_penalty;
      }
      return g;
    }
    int asym = std::abs(l1 - l2);
    if (asym > 1 || size > 4) {
      double g = asymmetry99(l1, l2, size);
      if (l1 > 1 && l2 > 1) {
        if (okb(e + 1) && okb(d - 1)) g += interior_mm(b(e), b(d), b(e + 1), b(d - 1));
        if (okb(i + 1) && okb(j - 1)) g += interior_mm(b(i), b(j), b(i + 1), b(j - 1));
      } else {
        g += interior_mm(b(e), b(d), 1, 1);
        g += interior_mm(b(i), b(j), 1, 1);
      }
      return g;
    }
    if (asym == 0 && size == 2) {
      if (!(okb(i + 1) && okb(j - 1))) return 0.0;
      int cpi = mshift(b(i), b(j)), ipi = mshift(b(d), b(e));
      return 0.0 + T.sint2[96 * cpi + 16 * ipi + 4 * (b(i + 1) - 1) + (b(j - 1) - 1)];
    }
    if (asym == 0 && size == 4) {
      if (!(okb(i + 1) && okb(i + 2) && okb(j - 1) && okb(j - 2))) return 0.0;
      int cpi = mshift(b(i), b(j)), ipi = mshift(b(d), b(e));
      return 0.0 + T.sint4[cpi * 1536 + ipi * 256 + (4 * (b(i + 1) - 1) + b(j - 1) - 1) * 16 +
                           (4 * (b(i + 2) - 1) + b(j - 2) - 1)];
    }
    if (l1 == 1) {
      if (!(okb(i + 1) && okb(j - 1) && okb(j - 2))) return 0.0;
      int cpi = mshift(b(i), b(j)), ipi = mshift(b(d), b(e));
      return 0.0 + T.asint1x2[cpi * 384 + (b(j - 2) - 1) * 96 + (b(i + 1) - 1) * 24 + 4 * ipi +
                              ((b(j - 1) - 1) % 4)];
    }
    if (!(okb(i + 1) && okb(i + 2) && okb(j - 1))) return 0.0;
    int cpi = mshift(b(j), b(i)), ipi = mshift(b(e), b(d));
    return 0.0 + T.asint1x2[ipi * 384 + (b(i + 1) - 1) * 96 + (b(j - 1) - 1) * 24 + 4 * cpi +
                            ((b(i + 2) - 1) % 4)];
  }

  double stack(int i, int j) const {
    int a = pair_type(code[i], code[j]);
    int b = pair_type(code[i + 1], code[j - 1]);
    if (a < 0 || b < 0) return NEG_INF;
    return log_weight(STACK[a][b]);
  }

  double bulge(int i, int d, int e, int j, int nn) const {
    if (nn < 1) return NEG_INF;
    double g = size_continuation(nn, BULGE_INIT, BULGE_K);
    if (nn == 1) {
      int a = pair_type(code[i], code[j]);
      int b = pair_type(code[d], code[e]);
      g += (a >= 0 && b >= 0) ? STACK[a][b] : 0.0;
    }
    return log_weight(g);
  }

  double internal(int l1, int l2) const {
    int nn = l1 + l2;
    if (nn < 4) return NEG_INF;
    double g = size_continuation(nn, INTERIOR_INIT, INTERIOR_K);
    int asym = std::abs(l1 - l2);
    if (asym) g += std::min(3.0, 0.5 * asym);
    return log_weight(g);
  }

  // ---- vienna2004 (dcal/mol), LinearPartition-V's utility_v.h; bases A=1..U=4,
  // pair types CG=1 GC=2 GU=3 UG=4 AU=5 UA=6
  int vn(int k) const { return code[k] + 1; }

  int vtype(int i, int j) const {
    int a = code[i], c = code[j];
    if (a == 1 && c == 2) return 1;
    if (a == 2 && c == 1) return 2;
    if (a == 2 && c == 3) return 3;
    if (a == 3 && c == 2) return 4;
    if (a == 0 && c == 3) return 5;
    if (a == 3 && c == 0) return 6;
    return 0;
  }

  int vext(int size) const { return static_cast<int>(V04.lxc * std::log(size / 30.0)); }

  int hairpin04(int i, int j) const {
    const V04Tables& T = V04;
    int size = j - i - 1;
    int t = vtype(i, j);
    int e = size <= 30 ? T.hairpin[size] : T.hairpin[30] + vext(size);
    if (size == 4) {
      auto it = T.tetra.find(seq.substr(i, 6));
      if (it != T.tetra.end()) return it->second;
    } else if (size == 6) {
      auto it = T.hexa.find(seq.substr(i, 8));
      if (it != T.hexa.end()) return it->second;
    } else if (size == 3) {
      auto it = T.tri.find(seq.substr(i, 5));
      if (it != T.tri.end()) return it->second;
      return e + (t > 2 ? T.au : 0);
    }
    return e + T.mmH[t * 25 + vn(i + 1) * 5 + vn(j - 1)];
  }

  // Outer pair i·j, inner pair p·q (v_score_single).
  int interior04(int i, int j, int p, int q) const {
    const V04Tables& T = V04;
    int si1 = vn(i + 1), sj1 = vn(j - 1), sp1 = vn(p - 1), sq1 = vn(q + 1);
    int t = vtype(i, j);
    int t2 = vtype(q, p);
    int n1 = p - i - 1, n2 = j - q - 1;
    int nl = n1 > n2 ? n1 : n2;
    int ns = n1 > n2 ? n2 : n1;
    if (nl == 0) return T.stack[t * 8 + t2];
    if (ns == 0) {
      int e = nl <= 30 ? T.bulge[nl] : T.bulge[30] + vext(nl);
      if (nl == 1) {
        e += T.stack[t * 8 + t2];
      } else {
        if (t > 2) e += T.au;
        if (t2 > 2) e += T.au;
      }
      return e;
    }
    if (ns == 1) {
      if (nl == 1) return T.int11[((t * 8 + t2) * 5 + si1) * 5 + sj1];
      if (nl == 2) {
        if (n1 == 1) return T.int21[(((t * 8 + t2) * 5 + si1) * 5 + sq1) * 5 + sj1];
        return T.int21[(((t2 * 8 + t) * 5 + sq1) * 5 + si1) * 5 + sp1];
      }
      int e = nl + 1 <= 30 ? T.interior[nl + 1] : T.interior[30] + vext(nl + 1);
      e += std::min(T.max_ninio, (nl - ns) * T.ninio);
      return e + T.mm1n[t * 25 + si1 * 5 + sj1] + T.mm1n[t2 * 25 + sq1 * 5 + sp1];
    }
    if (ns == 2) {
      if (nl == 2) return T.int22[((((t * 8 + t2) * 5 + si1) * 5 + sp1) * 5 + sq1) * 5 + sj1];
      if (nl == 3)
        return T.interior[5] + T.ninio + T.mm23[t * 25 + si1 * 5 + sj1] +
               T.mm23[t2 * 25 + sq1 * 5 + sp1];
    }
    int u = nl + ns;
    int e = u <= 30 ? T.interior[u] : T.interior[30] + vext(u);
    e += std::min(T.max_ninio, (nl - ns) * T.ninio);
    return e + T.mmI[t * 25 + si1 * 5 + sj1] + T.mmI[t2 * 25 + sq1 * 5 + sp1];
  }

  // Pair i·j as an exterior-loop or multiloop branch, without ML_intern
  // (v_score_external_paired; E_MLstem uses the same values, mismatchM = mismatchExt).
  int branch04(int i, int j) const {
    const V04Tables& T = V04;
    int t = vtype(i, j);
    int e = t > 2 ? T.au : 0;
    if (dangles_mode == 3) {
      int si1 = i > 0 ? vn(i - 1) : -1;
      int sj1 = j < n - 1 ? vn(j + 1) : -1;
      if (si1 >= 0 && sj1 >= 0) e += T.mmExt[t * 25 + si1 * 5 + sj1];
      else if (si1 >= 0) e += T.d5[t * 5 + si1];
      else if (sj1 >= 0) e += T.d3[t * 5 + sj1];
    }
    return e;
  }

  // Pair i·j closing a multiloop, without ML_closing and ML_intern (v_score_multi).
  int closing04(int i, int j) const {
    int t = vtype(j, i);
    int e = t > 2 ? V04.au : 0;
    if (dangles_mode == 3) e += V04.mmM[t * 25 + vn(j - 1) * 5 + vn(i + 1)];
    return e;
  }
};

// ---------------------------------------------------- insertion-ordered maps
// int key in [0, n) -> double; deleted keys leave a tombstone (key = -1) so
// re-insertion goes to the end, like a Python dict.
// Open-addressed key -> slot. Keys are >= 0. Empty is -1, tombstone is -2.
struct Ix {
  std::vector<int> k, v;
  uint32_t mask = 0;
  int live = 0, tombs = 0;

  static uint32_t hash(int x) { return static_cast<uint32_t>(x) * 0x9e3779b1u; }
  void clear() { k.clear(); v.clear(); mask = 0; live = tombs = 0; }
  double bytes() const {
    return (k.capacity() ? static_cast<double>(k.capacity() * sizeof(int) + 16) : 0.0) +
           (v.capacity() ? static_cast<double>(v.capacity() * sizeof(int) + 16) : 0.0);
  }
  void grow() {
    size_t cap = k.empty() ? 8 : k.size() * 2;
    std::vector<int> nk(cap, -1), nv(cap);
    uint32_t m = static_cast<uint32_t>(cap - 1);
    for (size_t i = 0; i < k.size(); ++i) {
      if (k[i] < 0) continue;
      uint32_t s = hash(k[i]) & m;
      while (nk[s] >= 0) s = (s + 1) & m;
      nk[s] = k[i];
      nv[s] = v[i];
    }
    k.swap(nk);
    v.swap(nv);
    mask = m;
    tombs = 0;
  }
  int slot_of(int key) const {
    if (k.empty() || key < 0) return -1;
    uint32_t s = hash(key) & mask;
    for (;;) {
      int c = k[s];
      if (c == key) return v[s];
      if (c == -1) return -1;
      s = (s + 1) & mask;
    }
  }
  void put(int key, int slot) {
    if (k.empty() || (live + tombs + 1) * 2 > static_cast<int>(k.size())) grow();
    uint32_t s = hash(key) & mask;
    int tomb = -1;
    for (;;) {
      int c = k[s];
      if (c == key) {
        v[s] = slot;
        return;
      }
      if (c == -2 && tomb < 0) tomb = static_cast<int>(s);
      if (c == -1) {
        if (tomb >= 0) {
          k[tomb] = key;
          v[tomb] = slot;
          --tombs;
        } else {
          k[s] = key;
          v[s] = slot;
        }
        ++live;
        return;
      }
      s = (s + 1) & mask;
    }
  }
  void del(int key) {
    if (k.empty() || key < 0) return;
    uint32_t s = hash(key) & mask;
    for (;;) {
      int c = k[s];
      if (c == key) {
        k[s] = -2;
        --live;
        ++tombs;
        return;
      }
      if (c == -1) return;
      s = (s + 1) & mask;
    }
  }
};

struct OMap {
  std::vector<int> keys;
  std::vector<double> vals;
  Ix pos;
  int live = 0;

  void init(int) {
    keys.clear();
    vals.clear();
    pos.clear();
    live = 0;
  }
  bool has(int k) const { return pos.slot_of(k) >= 0; }
  double get(int k) const {
    int p = pos.slot_of(k);
    return p >= 0 ? vals[p] : NEG_INF;
  }
  void set(int k, double v) {
    int p = pos.slot_of(k);
    if (p >= 0) {
      vals[p] = v;
      return;
    }
    pos.put(k, static_cast<int>(keys.size()));
    keys.push_back(k);
    vals.push_back(v);
    ++live;
  }
  void erase(int k) {
    int p = pos.slot_of(k);
    if (p < 0) return;
    keys[p] = -1;
    pos.del(k);
    --live;
  }
  // Drops tombstones, keeping live keys in insertion order.  Must not be
  // called while iterating over slots().
  void compact() {
    size_t w = 0;
    for (size_t r = 0; r < keys.size(); ++r) {
      int k = keys[r];
      if (k < 0) continue;
      keys[w] = k;
      vals[w] = vals[r];
      ++w;
    }
    keys.resize(w);
    vals.resize(w);
    pos.clear();
    for (size_t r = 0; r < w; ++r) pos.put(keys[r], static_cast<int>(r));
    if (keys.capacity() > 2 * w + 16) {
      keys.shrink_to_fit();
      vals.shrink_to_fit();
    }
  }
  bool sparse() const { return keys.size() > 2 * static_cast<size_t>(live) + 16; }
  size_t slots() const { return keys.size(); }
};

// Nested table stored as the kept keys of one OMap column. ix(i, j) = i * n + j.
// A read of a missing key is NEG_INF. Assigning NEG_INF drops the key.
struct AtView {
  std::vector<OMap>* cols = nullptr;
  int n = 0;
  void bind(std::vector<OMap>* c, int n_) { cols = c; n = n_; }
  bool empty() const { return n == 0 || cols == nullptr; }
  std::vector<double> dense() const {
    std::vector<double> o(static_cast<size_t>(n) * n, NEG_INF);
    if (!cols) return o;
    for (int j = 0; j < n; ++j)
      for (size_t r = 0; r < (*cols)[j].slots(); ++r) {
        int i = (*cols)[j].keys[r];
        if (i >= 0) o[static_cast<size_t>(i) * n + j] = (*cols)[j].vals[r];
      }
    return o;
  }
  struct Proxy {
    OMap* col;
    int i;
    operator double() const { return col->get(i); }
    Proxy& operator=(double v) {
      if (v < HALF_NEG) {
        if (col->has(i)) col->erase(i);
      } else {
        col->set(i, v);
      }
      return *this;
    }
    Proxy& operator=(const Proxy& o) { return *this = static_cast<double>(o); }
  };
  Proxy operator[](size_t k) { return Proxy{&(*cols)[k % static_cast<size_t>(n)], static_cast<int>(k / n)}; }
  double operator[](size_t k) const { return (*cols)[k % static_cast<size_t>(n)].get(static_cast<int>(k / n)); }
};

// Outside (and other) values that are finite only on kept keys. Same ix layout.
struct SpMat {
  int n = 0;
  std::vector<OMap> col;
  void init(int n_) {
    n = n_;
    col.assign(n_, {});
  }
  bool empty() const { return n == 0; }
  void clear() { col.clear(); n = 0; }
  long long finite() const {
    long long e = 0;
    for (const OMap& m : col) e += m.live;
    return e;
  }
  double bytes() const {
    double b = col.capacity() ? static_cast<double>(col.capacity() * sizeof(OMap) + 16) : 0.0;
    for (const OMap& m : col)
      b += (m.keys.capacity() ? static_cast<double>(m.keys.capacity() * sizeof(int) + 16) : 0.0) +
           (m.vals.capacity() ? static_cast<double>(m.vals.capacity() * sizeof(double) + 16) : 0.0) +
           m.pos.bytes();
    return b;
  }
  std::vector<double> dense() const {
    std::vector<double> o(static_cast<size_t>(n) * n, NEG_INF);
    for (int j = 0; j < n; ++j)
      for (size_t r = 0; r < col[j].slots(); ++r) {
        int i = col[j].keys[r];
        if (i >= 0) o[static_cast<size_t>(i) * n + j] = col[j].vals[r];
      }
    return o;
  }
  struct Proxy {
    OMap* col;
    int i;
    operator double() const { return col->get(i); }
    Proxy& operator=(double v) {
      if (v < HALF_NEG) {
        if (col->has(i)) col->erase(i);
      } else {
        col->set(i, v);
      }
      return *this;
    }
    Proxy& operator=(const Proxy& o) { return *this = static_cast<double>(o); }
  };
  Proxy operator[](size_t k) { return Proxy{&col[k % static_cast<size_t>(n)], static_cast<int>(k / n)}; }
  double operator[](size_t k) const { return col[k % static_cast<size_t>(n)].get(static_cast<int>(k / n)); }
};

// Qz is finite on every i <= j, so the column is a packed triangle rather than a hash.
struct TriMat {
  int n = 0;
  std::vector<std::vector<double>> col;
  void init(int n_) {
    n = n_;
    col.resize(n_);
    for (int j = 0; j < n_; ++j) col[j].assign(j + 1, NEG_INF);
  }
  bool empty() const { return n == 0; }
  long long finite() const {
    long long e = 0;
    for (const auto& c : col)
      for (double x : c) e += x > HALF_NEG;
    return e;
  }
  double bytes() const {
    double b = col.capacity() ? static_cast<double>(col.capacity() * sizeof(std::vector<double>) + 16) : 0.0;
    for (const auto& c : col)
      if (c.capacity()) b += static_cast<double>(c.capacity() * sizeof(double) + 16);
    return b;
  }
  struct Proxy {
    TriMat* m;
    int i, j;
    operator double() const {
      if (j < 0 || j >= m->n || i < 0 || i >= static_cast<int>(m->col[j].size())) return NEG_INF;
      return m->col[j][i];
    }
    Proxy& operator=(double v) {
      if (j < 0 || j >= m->n || i < 0) return *this;
      if (i >= static_cast<int>(m->col[j].size())) m->col[j].resize(i + 1, NEG_INF);
      m->col[j][i] = v;
      return *this;
    }
    Proxy& operator=(const Proxy& o) { return *this = static_cast<double>(o); }
  };
  Proxy operator[](size_t k) {
    return Proxy{this, static_cast<int>(k / n), static_cast<int>(k % static_cast<size_t>(n))};
  }
  double operator[](size_t k) const {
    int i = static_cast<int>(k / n), j = static_cast<int>(k % static_cast<size_t>(n));
    if (j < 0 || j >= n || i < 0 || i >= static_cast<int>(col[j].size())) return NEG_INF;
    return col[j][i];
  }
};

// Values over [0, n) with NEG_INF for absent keys; keys() in first-insertion order.
struct Col {
  std::vector<double> v;
  std::vector<char> in;
  std::vector<int> keys;
  void init(int n) {
    v.assign(n, NEG_INF);
    in.assign(n, 0);
    keys.clear();
  }
  void clear() {
    for (int k : keys) {
      v[k] = NEG_INF;
      in[k] = 0;
    }
    keys.clear();
  }
  double get(int k) const { return v[k]; }
  void put(int k, double x) {
    if (!in[k]) {
      in[k] = 1;
      keys.push_back(k);
    }
    v[k] = x;
  }
  void add(int k, double x) { put(k, log_add(in[k] ? v[k] : NEG_INF, x)); }
};

// Trailing-run states of one column (see DP03::run_ext); outside keeps U1, U2 only.
struct RunSt {
  Col P, U1, U2;
  void init(int n) {
    P.init(n);
    U1.init(n);
    U2.init(n);
  }
  void clear() {
    P.clear();
    U1.clear();
    U2.clear();
  }
};

struct Ent {
  int p, u;
  double a, b;  // inside, outside
};

struct Grp {
  int v;
  std::vector<Ent> ents;
  std::unordered_map<uint64_t, int> idx;
};

inline uint64_t key2(int p, int u) {
  return (static_cast<uint64_t>(static_cast<uint32_t>(p)) << 32) |
         static_cast<uint32_t>(u);
}

// groups[v][(p, u)] with insertion order at both levels.
struct Tab4 {
  std::vector<Grp> grps;
  std::unordered_map<int, int> vidx;

  Grp* get(int v) {
    auto it = vidx.find(v);
    return it == vidx.end() ? nullptr : &grps[it->second];
  }
  Ent* find(int v, int p, int u) {
    Grp* g = get(v);
    if (!g) return nullptr;
    auto it = g->idx.find(key2(p, u));
    return it == g->idx.end() ? nullptr : &g->ents[it->second];
  }
  void acc(int v, int p, int u, double val) {
    if (val < HALF_NEG) return;
    Grp* g;
    auto it = vidx.find(v);
    if (it == vidx.end()) {
      vidx.emplace(v, static_cast<int>(grps.size()));
      grps.push_back(Grp{v, {}, {}});
      g = &grps.back();
    } else {
      g = &grps[it->second];
    }
    uint64_t k = key2(p, u);
    auto jt = g->idx.find(k);
    if (jt == g->idx.end()) {
      g->idx.emplace(k, static_cast<int>(g->ents.size()));
      g->ents.push_back(Ent{p, u, log_add(NEG_INF, val), NEG_INF});
    } else {
      Ent& e = g->ents[jt->second];
      e.a = log_add(e.a, val);
    }
  }
};

inline double g4(Tab4& t, int v, int p, int u) {
  Ent* e = t.find(v, p, u);
  return e ? e->b : NEG_INF;
}

inline void add4(Ent& e, double val) {
  if (val < HALF_NEG) return;
  e.b = log_add(e.b, val);
}

// ------------------------------------------ LinearAccess pre-pass (arm score, gate)
// Mirrors dp03/access.py.  Format 3 (natural log, 1-based): header line
// "#dp03-access-format 3 ...", then per sequence "S seq", "A a b logacc",
// "D a b y logQ", "P i j p" (every pair with p >= 1e-5).  Stored 0-based; a missing
// A line means logacc = ACC_LOG_FLOOR, missing D lines mean log D = -inf.
const char* const ACCESS_FORMAT_HEADER = "#dp03-access-format 3";
constexpr double ACC_LOG_FLOOR = -200.0;
constexpr int LK_HIST_BINS = 9;
constexpr double LK_HIST_EDGES[LK_HIST_BINS - 1] = {0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0};

// Same operations as access.log_add (true -inf, not NEG_INF).
inline double acc_log_add(double a, double b) {
  if (a == -INFINITY) return b;
  if (b == -INFINITY) return a;
  if (a > b) return a + std::log1p(std::exp(b - a));
  return b + std::log1p(std::exp(a - b));
}

struct AccWin {
  int t;
  std::vector<int> ys;       // innermost partners, ascending
  std::vector<double> lsuf;  // lsuf[k] = log sum_{m >= k} Q(ys[m])
};

struct AccessRec {
  int n = 0;
  std::unordered_map<uint64_t, double> logacc_tab, bpp;
  std::vector<std::vector<AccWin>> by_start;  // sorted by t

  double log_acc(int a, int b) const {
    auto it = logacc_tab.find(key2(a, b));
    return it == logacc_tab.end() ? ACC_LOG_FLOOR : it->second;
  }
  bool has_acc(int a, int b) const { return logacc_tab.count(key2(a, b)) != 0; }
  double pair_prob(int i, int j) const {
    auto it = bpp.find(key2(i, j));
    return it == bpp.end() ? 0.0 : it->second;
  }
  // The P lines as (i, j, p), 0-based, sorted by (i, j), one entry per pair.
  std::vector<gate::Pair> pairs() const {
    std::vector<gate::Pair> out;
    out.reserve(bpp.size());
    for (const auto& kv : bpp)
      out.push_back({static_cast<int>(kv.first >> 32), static_cast<int>(kv.first & 0xffffffffu),
                     kv.second});
    std::sort(out.begin(), out.end(), [](const gate::Pair& a, const gate::Pair& b) {
      return a.i != b.i ? a.i < b.i : a.j < b.j;
    });
    return out;
  }
};

// sequence -> byte span [start, end) of its record; the file stays open for seeks.
struct AccessIndex {
  FILE* f = nullptr;
  std::string path;
  std::unordered_map<std::string, std::pair<long, long>> span;
};

bool read_line(FILE* f, std::string& line) {
  line.clear();
  int c;
  while ((c = std::fgetc(f)) != EOF) {
    line.push_back(static_cast<char>(c));
    if (c == '\n') break;
  }
  return !line.empty();
}

AccessIndex index_access_file(const std::string& path) {
  AccessIndex idx;
  idx.path = path;
  idx.f = std::fopen(path.c_str(), "rb");
  if (!idx.f) {
    std::fprintf(stderr, "cannot open access file %s\n", path.c_str());
    std::exit(2);
  }
  std::string line;
  read_line(idx.f, line);
  if (line.compare(0, std::strlen(ACCESS_FORMAT_HEADER), ACCESS_FORMAT_HEADER) != 0) {
    std::fprintf(stderr,
                 "%s: not an access file of format 3 (older format?); "
                 "rebuild it with access/access_prepass.py\n",
                 path.c_str());
    std::exit(2);
  }
  std::string cur;
  bool have = false;
  long start = 0, pos = std::ftell(idx.f);
  while (read_line(idx.f, line)) {
    if (line.size() >= 2 && line[0] == 'S' && line[1] == ' ') {
      if (have) idx.span.emplace(cur, std::make_pair(start, pos));
      size_t b = 2, e = line.size();
      while (b < e && std::isspace(static_cast<unsigned char>(line[b]))) ++b;
      while (e > b && std::isspace(static_cast<unsigned char>(line[e - 1]))) --e;
      cur = line.substr(b, e - b);
      have = true;
      start = pos;
    }
    pos += static_cast<long>(line.size());
  }
  if (have) idx.span.emplace(cur, std::make_pair(start, pos));
  return idx;
}

AccessRec load_access(const AccessIndex& idx, const std::string& seq) {
  AccessRec rec;
  const int n = static_cast<int>(seq.size());
  rec.n = n;
  rec.by_start.assign(n, {});
  auto sp = idx.span.find(seq);
  if (sp == idx.span.end()) {
    std::fprintf(stderr, "[arm] warning: no access record for sequence of length %d in %s; D = 0\n", n,
                 idx.path.c_str());
    return rec;
  }
  std::fseek(idx.f, sp->second.first, SEEK_SET);
  std::string block(static_cast<size_t>(sp->second.second - sp->second.first), '\0');
  if (std::fread(&block[0], 1, block.size(), idx.f) != block.size()) {
    std::fprintf(stderr, "short read in access file %s\n", idx.path.c_str());
    std::exit(2);
  }
  std::unordered_map<uint64_t, std::vector<std::pair<int, double>>> qs;
  std::vector<uint64_t> order;
  size_t p0 = 0;
  while (p0 < block.size()) {
    size_t p1 = block.find('\n', p0);
    if (p1 == std::string::npos) p1 = block.size();
    const std::string line = block.substr(p0, p1 - p0);
    const char* ln = line.c_str();
    char tag[8], x1[64], x2[64], x3[64], x4[64], x5[8];
    int k = std::sscanf(ln, "%7s %63s %63s %63s %63s %7s", tag, x1, x2, x3, x4, x5);
    if (k == 4 && std::strcmp(tag, "A") == 0) {
      int a = std::atoi(x1) - 1, b = std::atoi(x2) - 1;
      if (0 <= a && a <= b && b < n) rec.logacc_tab[key2(a, b)] = std::strtod(x3, nullptr);
    } else if (k == 4 && std::strcmp(tag, "P") == 0) {
      int i = std::atoi(x1) - 1, j = std::atoi(x2) - 1;
      if (0 <= i && i < j && j < n) rec.bpp[key2(i, j)] = std::strtod(x3, nullptr);
    } else if (k == 5 && std::strcmp(tag, "D") == 0) {
      int a0 = std::atoi(x1) - 1, b0 = std::atoi(x2) - 1, y0 = std::atoi(x3) - 1;
      double lq = std::strtod(x4, nullptr);
      if (0 <= a0 && a0 <= b0 && b0 < y0 && y0 < n) {
        uint64_t key = key2(a0, b0);
        auto it = qs.find(key);
        if (it == qs.end()) {
          order.push_back(key);
          it = qs.emplace(key, std::vector<std::pair<int, double>>{}).first;
        }
        bool found = false;
        for (auto& e : it->second)
          if (e.first == y0) {
            e.second = acc_log_add(e.second, lq);
            found = true;
            break;
          }
        if (!found) it->second.emplace_back(y0, acc_log_add(-INFINITY, lq));
      }
    }
    p0 = p1 + 1;
  }
  for (uint64_t key : order) {
    auto slot = qs[key];
    std::sort(slot.begin(), slot.end());
    AccWin w;
    w.t = static_cast<int>(key & 0xffffffffu);
    w.ys.resize(slot.size());
    w.lsuf.assign(slot.size(), -INFINITY);
    double run = -INFINITY;
    for (int m = static_cast<int>(slot.size()) - 1; m >= 0; --m) {
      run = acc_log_add(run, slot[m].second);
      w.ys[m] = slot[m].first;
      w.lsuf[m] = run;
    }
    rec.by_start[static_cast<int>(key >> 32)].push_back(std::move(w));
  }
  for (auto& row : rec.by_start)
    std::stable_sort(row.begin(), row.end(),
                     [](const AccWin& x, const AccWin& y) { return x.t < y.t; });
  return rec;
}

struct ArmStats {
  long long scored = 0, arm_wins = 0, windows = 0, floor_windows = 0, chosen_floor = 0;
  double sum_log_kappa = 0.0, max_log_kappa = 0.0;
  long long lk_hist[LK_HIST_BINS] = {};
};

struct ArmCand {
  double a, z;
  bool floor;
};

// ------------------------------------------------------- eligibility gate (Qg / Qgl)
// A four-index state with outer pair p.q and gap [u+1, v-1] gets
//   g1 = min(1, max_{x in gap} r(x, q) / max(delta, 1 - P(p, q)))       (stem-1 role)
//   g2 = the same on the mirrored pairs, partners left of p              (stem-2 role)
//   g  = g1 (role stem1) or max(g1, g2) (role both); an empty gap gives g = 0,
// with r(x, q) = sum of P(x, y) over q < y <= q + span (gate_index.hpp).  At every
// pruning event of a gated table, states with g < tau are dropped before ranking; the
// rest are ranked by the unchanged pruning score, and `reserve` more eligible states
// with the highest g may be added.  With tau > 0 states are removed even at infinite
// beam, so infinite-beam exactness holds only for tau = 0.
enum GateTable { GATE_QG = 0, GATE_QGL = 1 };

struct GateCfg {
  bool on = false;
  double tau = 0.01;
  int span = 0;
  double delta = 1e-3;
  bool both = true;
  bool tables[2] = {true, true};
  int reserve = 0;
};

class Gate {
 public:
  Gate(int n, const AccessRec& acc, int span, double delta, bool both)
      : n_(n), acc_(acc), delta_(delta), both_(both) {
    std::vector<gate::Pair> pairs = acc.pairs();
    fwd_ = gate::PersistentGate(n, pairs, span);
    mir_ = gate::PersistentGate(n, gate::mirror_pairs(n, pairs), span);
  }

  void roles(int p, int u, int v, int q, double& g1, double& g2) const {
    const int lo = u + 1, hi = v - 1;
    if (lo > hi) {
      g1 = g2 = 0.0;
      return;
    }
    const double denom = std::max(delta_, 1.0 - acc_.pair_prob(p, q));
    g1 = std::min(1.0, fwd_.range_max(q, lo, hi) / denom);
    g2 = std::min(1.0, mir_.range_max(n_ - 1 - p, n_ - 1 - hi, n_ - 1 - lo) / denom);
  }

  double g(int p, int u, int v, int q) const {
    double g1, g2;
    roles(p, u, v, q, g1, g2);
    return both_ ? std::max(g1, g2) : g1;
  }

 private:
  int n_;
  const AccessRec& acc_;
  double delta_;
  bool both_;
  gate::PersistentGate fwd_, mir_;
};

// Per gated table, summed over pruning events: candidates seen, eligible (g >= tau),
// ranked (scored by the pruning score), kept, and kept through the reserve.
struct GateStats {
  long long cand[2] = {}, elig[2] = {}, ranked[2] = {}, kept[2] = {}, reserved[2] = {};
};

// --beam-stats=<file>: opt-in TSV logging only (no change to ranking keys or keep sets).
// pruned_key_share = sum_{pruned} exp(key) / sum_{all finite} exp(key) (0 if denom is 0).
// One P/Q row per (pass, table, position); pass is stem|pass1|pass2|scan (set by the driver).
struct BeamStats {
  FILE* f = nullptr;
  std::string pass = "-";
  struct Agg {
    long long n_pos = 0, n_full = 0;
    double sum_share = 0.0, max_share = 0.0;
  };
  // Key = pass + "\t" + table so stem/pass1/scan summaries stay separate.
  std::map<std::string, Agg> by_table;
  void open(const std::string& path) {
    f = std::fopen(path.c_str(), "w");
    if (!f) return;
    std::fprintf(f, "# beam-stats TSV: Q pass t logQ_1_t | P pass table t candidates kept "
                    "pruned_key_share full | S pass table n_pos n_full max_share mean_share\n");
    std::fprintf(f, "kind\tpass\tname\tt\tcandidates\tkept\tpruned_key_share\tfull\n");
    std::fflush(f);
  }
  void set_pass(const char* p) { pass = p ? p : "-"; }
  static double mass_share(const std::vector<double>& keys, size_t keep_n) {
    // keys sorted descending (best first); first keep_n kept, rest pruned.
    double mx = NEG_INF;
    for (double k : keys)
      if (k > mx && k > HALF_NEG) mx = k;
    if (!(mx > HALF_NEG)) return 0.0;
    double all = 0.0, pr = 0.0;
    for (size_t i = 0; i < keys.size(); ++i) {
      if (!(keys[i] > HALF_NEG)) continue;
      const double w = std::exp(keys[i] - mx);
      all += w;
      if (i >= keep_n) pr += w;
    }
    return all > 0.0 ? pr / all : 0.0;
  }
  void log_Q(int t, double logq) {
    if (!f) return;
    std::fprintf(f, "Q\t%s\tQ1t\t%d\t-\t-\t%.17g\t-\n", pass.c_str(), t, logq);
  }
  void log_prune(const char* table, int t, size_t cand, size_t kept, double share, int full) {
    if (!f) return;
    std::fprintf(f, "P\t%s\t%s\t%d\t%zu\t%zu\t%.17g\t%d\n", pass.c_str(), table, t, cand, kept,
                 share, full);
    Agg& a = by_table[pass + "\t" + table];
    ++a.n_pos;
    if (full) ++a.n_full;
    a.sum_share += share;
    if (share > a.max_share) a.max_share = share;
  }
  // Unsorted keys + keep mask: compute share of pruned mass.
  void log_prune_mask(const char* table, int t, const std::vector<double>& keys,
                      const std::vector<char>& keep, int beam) {
    if (!f) return;
    const size_t cand = keys.size();
    size_t kept = 0;
    for (char k : keep)
      if (k) ++kept;
    double mx = NEG_INF;
    for (double k : keys)
      if (k > mx && k > HALF_NEG) mx = k;
    double all = 0.0, pr = 0.0;
    if (mx > HALF_NEG) {
      for (size_t i = 0; i < keys.size(); ++i) {
        if (!(keys[i] > HALF_NEG)) continue;
        const double w = std::exp(keys[i] - mx);
        all += w;
        if (!keep[i]) pr += w;
      }
    }
    const double share = all > 0.0 ? pr / all : 0.0;
    const int full = (beam > 0 && static_cast<int>(cand) > beam) ? 1 : 0;
    log_prune(table, t, cand, kept, share, full);
  }
  void flush_summary() {
    if (!f) return;
    for (const auto& kv : by_table) {
      const Agg& a = kv.second;
      const double mean = a.n_pos ? a.sum_share / static_cast<double>(a.n_pos) : 0.0;
      // kv.first is "pass\ttable"
      const size_t tab = kv.first.find('\t');
      const std::string p = tab == std::string::npos ? "-" : kv.first.substr(0, tab);
      const std::string name = tab == std::string::npos ? kv.first : kv.first.substr(tab + 1);
      std::fprintf(f, "S\t%s\t%s\t%lld\t%lld\t%.17g\t%.17g\t-\n", p.c_str(), name.c_str(),
                   static_cast<long long>(a.n_pos), static_cast<long long>(a.n_full), a.max_share,
                   mean);
    }
    std::fflush(f);
  }
  void close() {
    if (!f) return;
    flush_summary();
    std::fclose(f);
    f = nullptr;
  }
};
static BeamStats g_beam_stats;

// Beamed tables, for --beam-table (per-table beams) and the state / drop dumps.
enum BeamTable { BT_QB, BT_QM, BT_QM2, BT_QMP, BT_PKM, BT_QP, BT_QG, BT_QGLS, BT_QGL, BT_QGR,
                 BT_QGRS, BT_N };
const char* const BT_NAMES[BT_N] = {"Qb", "Qm",  "Qm2", "Qmp", "PKM", "Qp",
                                    "Qg", "Qgls", "Qgl", "Qgr", "Qgrs"};

int beam_table_id(const std::string& s) {
  for (int t = 0; t < BT_N; ++t)
    if (s == BT_NAMES[t]) return t;
  return -1;
}

// State key (table, a, b, c, d) of the dumps: 4-index tables (p, u, v, q); Qb, Qm, Qm2,
// Qmp, Qp (i, j, -1, -1); PKM (d, j, -1, -1).  Indices < 16383.
inline uint64_t state_key(int bt, int a, int b, int c, int d) {
  return (static_cast<uint64_t>(bt) << 56) | (static_cast<uint64_t>(a + 1) << 42) |
         (static_cast<uint64_t>(b + 1) << 28) | (static_cast<uint64_t>(c + 1) << 14) |
         static_cast<uint64_t>(d + 1);
}

// --dump-pruned: every candidate dropped by a pruning event (beamed run), optionally
// restricted to the states of an exact --dump-states file (w = their exact weight).
struct DropLog {
  FILE* f = nullptr;
  const std::unordered_map<uint64_t, double>* exact = nullptr;
  std::string id = "-";
  long long events[BT_N] = {}, dropped[BT_N] = {}, matched[BT_N] = {};
  double w[BT_N] = {};
  // with exact: matched w by relative scan position, bin = floor(10 j / n); D lines only for
  // w >= min_w
  int n = 1;
  double min_w = 0.0;
  double wpos[BT_N][10] = {};
  // --dump-kept: K lines of the kept states of every stem / helper ranking event
  FILE* kf = nullptr;
};

// --stem-est: stem-2 completion estimate from a nested-only run (dp03/stem_est.py).
// C0 = nested prefix C, S[y] = nested exterior suffix [y, n-1] (S[n] = 0), and the index of
// w(x', y') = U(x', y') C0[y'-1]; H(u, v, q) = log sum of w over x' in [u+1, v-1], y' > q.
constexpr double STEM_LOG_SPARSE = 690.77552789821368;  // entries below 1e-300 of the largest

// Mirror (--stem-est-role all): w2(x', y') = U2(x', y') S[x'+1], U2 with the nested prefix
// C0[x0-1] before the outermost pair, at (n-1-y', n-1-x'); H2(p, u, v) = log sum of w2 over
// y' in [u+1, v-1], x' < p.
struct StemEst {
  int n = 0, loop = 0;
  std::vector<double> C0, S;
  stemidx::PersistentSum index, index2;
  bool mirror = false;
  size_t entries = 0, dropped = 0, entries2 = 0, dropped2 = 0;
  // --stem-est-helpers nested: log outside weight of pair (x, y) in the nested-only run (its
  // Qb outside value; pairs without one get the exterior-loop value C0[x-1] + at + S[y+1]).
  std::vector<double> beta;
  // --stem-est-chain: kept pairs (nested pair probability >= floor), the kept right partners of
  // each x ascending, and mbar(i, j) = the nested multiloop interior Qm (Qmp with stem_ml).
  bool chain = false;
  std::vector<char> kept;
  std::vector<std::vector<int>> kept_right;
  std::vector<double> mbar;
  size_t n_kept = 0;
  // --stem-est-hr av: log phi_gamma(x, y) of every pair (availability from the nested pair
  // probabilities) and the index of Uav(x', y') = E_in Tav (no prefix / suffix factors).
  // --stem-est-hr phi: the same lphi, index_av holds the original entries with phi_gamma
  // multiplied in (hr_phi).
  bool hr_av = false, hr_phi = false;
  std::vector<double> lphi;
  stemidx::PersistentSum index_av;
  size_t entries_av = 0, dropped_av = 0;
  // --stem-est-principled: from the nested pre-pass (exact for n <= 1000, else the nested
  // beam, every pair span): exterior suffix prS[y] = Sbar[y, n-1] (prS[n] = 0), gap filler
  // prQz = Qzbar, multiloop interior prMbar (Qmp with stem_ml), prBml[d][e] = every way the
  // closed segment [d, e] sits as one pseudoknot branch of a multiloop (DP03's pseudoknot
  // branch edges, beta1m excluded), kept pairs (pair probability >= --pr-floor);
  // --stem-est-principled-nested: the outsides of the closed Qb / Qm / Qm2 values.
  bool pr = false, pr_nested = false;
  std::vector<double> prS, prQz, prMbar, prBml, prQb_b, prQm_b, prQm2_b;
  // --sweep: finite PrMbar cells. prMbarAt[j] = left ends at right end j;
  // prMbarRight[i] = right ends at left end i, sorted.
  std::vector<std::vector<int>> prMbarAt, prMbarRight;
  long long pr_mbar_n = 0;
  int pr_mbar_max = 0;
  // Hyperedge sweep. prGamma / prMu are pair maps (Qm and Qm2 maps are freed after they push).
  // prMu keys are key2(b, 0). prMuFloor is the smallest kept mu entry (--v-floor).
  std::unordered_map<uint64_t, std::unordered_map<uint64_t, double>> prGamma, prMu;
  std::unordered_map<uint64_t, double> prMuFloor;
  long long pr_sw_states_g = 0, pr_sw_states_v = 0, pr_sw_pushes = 0;
  double pr_sw_sec = 0;
  // --pr-v-diag / --pr-v-ratio: the nested inside log Qb of every pair.
  std::vector<double> prQbIn;
  double QbIn(int i, int j) const {
    if (skip_qbin) lean_abort("QbIn");
    return prQbIn[static_cast<size_t>(i) * n + j];
  }
  std::vector<char> prKept;
  // --l2r-lean: the structures not built (index, S, C0, beta, chain; prKept / prQbIn); any
  // accessor of one aborts
  bool lean = false, skip_kept = false, skip_qbin = false;
  [[noreturn]] static void lean_abort(const char* what) {
    std::fprintf(stderr, "--l2r-lean: %s read but not built\n", what);
    std::fflush(nullptr);
    std::abort();
  }
  std::vector<std::vector<int>> prKeptRight;
  size_t pr_n_kept = 0;
  // --pr-partner-floor: the pairs of the future-stem / partner profiles (pair probability
  // >= that floor); empty: the kept pairs above.
  bool pr_partner = false;
  std::vector<std::vector<int>> prKeptRightP;
  size_t pr_n_partner = 0;
  // --pr-oracle: nested log pair probability of every pair (NEG_INF if none).
  bool want_lp = false;
  std::vector<double> prLP;
  long long bml_check_cells = 0, bml_check_mis = 0;
  double bml_check_rel = 0.0, pr_prep_sec = 0.0, prof_qz_pre = 0.0;
  // Q^z_pre kept per right end l (starts k, sorted by C(k-1)+Q^z). qz_on == 0 is the dense copy.
  struct QzEnt {
    int k;
    double q, p, u1, u2;
  };
  int qz_on = 0;
  std::vector<std::vector<QzEnt>> qz_col;
  std::vector<double> prC, prRight;
  std::vector<std::vector<int>> prKeptRightV, prKeptRightG;
  long long p2_qz_n = 0;
  double p2_qz_sec = 0.0;
  double Qzbar(int a, int b) const {
    if (a > b) return 0.0;
    if (a < 0 || b >= n) return NEG_INF;
    if (qz_on) {
      if (static_cast<size_t>(b) >= qz_col.size()) return NEG_INF;
      const std::vector<QzEnt>& col = qz_col[static_cast<size_t>(b)];
      auto it = std::lower_bound(col.begin(), col.end(), a,
                                 [](const QzEnt& e, int k) { return e.k < k; });
      if (it != col.end() && it->k == a) return it->q;
      return NEG_INF;
    }
    return prQz[static_cast<size_t>(a) * n + b];
  }
  double Bml(int d, int e) const {
    return prBml.empty() ? NEG_INF : prBml[static_cast<size_t>(d) * n + e];
  }
  double PrMbar(int i, int j) const { return prMbar[static_cast<size_t>(i) * n + j]; }
  double H(int u, int v, int q) const {
    if (lean) lean_abort("H");
    return index.range_logsum(q, u + 1, v - 1);
  }
  double Hav(int u, int v, int q) const { return index_av.range_logsum(q, u + 1, v - 1); }
  double H2(int p, int u, int v) const {
    if (lean) lean_abort("H2");
    return index2.range_logsum(n - 1 - p, n - v, n - 2 - u);
  }
  double B(int x, int y) const {
    if (lean) lean_abort("B");
    return beta[static_cast<size_t>(x) * n + y];
  }
  double Mbar(int i, int j) const {
    if (lean) lean_abort("Mbar");
    return mbar[static_cast<size_t>(i) * n + j];
  }
  size_t bytes() const {
    size_t b = 12 * entries + 16 * index.nodes() + 8 * (2 * n + 1);
    if (mirror) b += 12 * entries2 + 16 * index2.nodes();
    b += 8 * beta.size();
    if (chain) b += kept.size() + 8 * mbar.size() + 4 * n_kept;
    if (hr_av) b += 8 * lphi.size() + 12 * entries_av + 16 * index_av.nodes();
    return b;
  }
};

// Entries with log w >= max - STEM_LOG_SPARSE; returns how many were dropped.
inline size_t stem_sparse(std::vector<stemidx::Entry>& ents) {
  double mx = NEG_INF;
  for (const auto& e : ents) mx = std::max(mx, e.lw);
  const double floor = mx - STEM_LOG_SPARSE;
  std::vector<stemidx::Entry> kept;
  for (const auto& e : ents)
    if (e.lw >= floor) kept.push_back(e);
  size_t dropped = ents.size() - kept.size();
  ents.swap(kept);
  return dropped;
}

// --mem-log on: stage checkpoints on stderr (RSS, live heap = malloc in use, CPU) with per
// structure entries and bytes.  Bytes: container capacities plus 16 per heap block and the
// hash nodes / buckets of the unordered maps (libstdc++ layout, no cached hash).
struct MemItem {
  std::string name;
  long long ents;
  double bytes;
};
namespace msz {
inline double node(size_t payload) {
  return static_cast<double>(std::max<size_t>(32, (payload + 16 + 15) & ~static_cast<size_t>(15)));
}
template <class T>
double vec(const std::vector<T>& v) {
  return v.capacity() ? static_cast<double>(v.capacity() * sizeof(T) + 16) : 0.0;
}
template <class T>
double vec2(const std::vector<std::vector<T>>& v) {
  double b = vec(v);
  for (const auto& x : v) b += vec(x);
  return b;
}
template <class M>
double umap(const M& m) {
  return static_cast<double>(m.bucket_count() * sizeof(void*)) +
         static_cast<double>(m.size()) * node(sizeof(typename M::value_type));
}
inline long long fin(const std::vector<double>& v) {
  long long c = 0;
  for (double x : v) c += x > HALF_NEG;
  return c;
}
}  // namespace msz
static double mem_rss_mb() {
  long pages = 0, res = 0;
  if (FILE* f = std::fopen("/proc/self/statm", "r")) {
    if (std::fscanf(f, "%ld %ld", &pages, &res) != 2) res = 0;
    std::fclose(f);
  }
  return static_cast<double>(res) * static_cast<double>(sysconf(_SC_PAGESIZE)) / 1048576.0;
}
static double mem_heap_mb() {
  const struct mallinfo2 mi = mallinfo2();
  return static_cast<double>(mi.uordblks + mi.hblkhd) / 1048576.0;
}
inline void stem_mem_items(const StemEst& s, std::vector<MemItem>& out) {
  const size_t se = s.entries + s.entries2 + s.entries_av;
  out.push_back({"stem_idx", static_cast<long long>(se),
                 12.0 * static_cast<double>(se) +
                     16.0 * static_cast<double>(s.index.nodes() + s.index2.nodes() + s.index_av.nodes())});
  out.push_back({"stem_S_C0", static_cast<long long>(s.S.size() + s.C0.size()), msz::vec(s.S) + msz::vec(s.C0)});
  out.push_back({"nested_beta", msz::fin(s.beta), msz::vec(s.beta)});
  out.push_back({"chain", static_cast<long long>(s.n_kept),
                 msz::vec(s.kept) + msz::vec2(s.kept_right) + msz::vec(s.mbar)});
  out.push_back({"hr_lphi", msz::fin(s.lphi), msz::vec(s.lphi)});
  out.push_back({"pr_Qz", msz::fin(s.prQz), msz::vec(s.prQz)});
  out.push_back({"pr_Mbar", msz::fin(s.prMbar), msz::vec(s.prMbar)});
  out.push_back({"pr_QbIn", msz::fin(s.prQbIn), msz::vec(s.prQbIn)});
  out.push_back({"Bml", msz::fin(s.prBml), msz::vec(s.prBml)});
  out.push_back({"pr_nested_b", msz::fin(s.prQb_b) + msz::fin(s.prQm_b) + msz::fin(s.prQm2_b),
                 msz::vec(s.prQb_b) + msz::vec(s.prQm_b) + msz::vec(s.prQm2_b)});
  out.push_back({"pr_kept", static_cast<long long>(s.pr_n_kept + s.pr_n_partner),
                 msz::vec(s.prKept) + msz::vec2(s.prKeptRight) + msz::vec2(s.prKeptRightP)});
  out.push_back({"pr_LP", msz::fin(s.prLP), msz::vec(s.prLP)});
  out.push_back({"pr_S", static_cast<long long>(s.prS.size()), msz::vec(s.prS)});
}
// One [mem] line: stage, CPU (total, since the previous checkpoint), RSS, peak RSS, live heap,
// then name:entries:MB:bytes_per_entry per item.
static bool g_mem_log = false;
static std::string g_mem_sid = "-";
static double g_mem_cpu = 0.0;
static void mem_log(const char* stage, const std::vector<MemItem>& items = {}) {
  if (!g_mem_log) return;
  const double c = static_cast<double>(std::clock()) / CLOCKS_PER_SEC;
  double hwm = 0.0;
  if (FILE* f = std::fopen("/proc/self/status", "r")) {
    char line[256];
    while (std::fgets(line, sizeof line, f))
      if (std::strncmp(line, "VmHWM:", 6) == 0) hwm = std::atof(line + 6) / 1024.0;
    std::fclose(f);
  }
  std::string s;
  char b[256];
  for (const MemItem& it : items) {
    if (it.ents == 0 && it.bytes == 0.0) continue;
    std::snprintf(b, sizeof b, " %s:%lld:%.1f:%.1f", it.name.c_str(), it.ents, it.bytes / 1048576.0,
                  it.ents > 0 ? it.bytes / static_cast<double>(it.ents) : 0.0);
    s += b;
  }
  std::fprintf(stderr, "[mem] %s stage %s cpu %.2f dcpu %.2f rss_mb %.0f hwm_mb %.0f heap_mb %.0f |%s\n",
               g_mem_sid.c_str(), stage, c, c - g_mem_cpu, mem_rss_mb(), hwm, mem_heap_mb(), s.c_str());
  std::fflush(stderr);
  g_mem_cpu = c;
}

// --sample-prof PATH: SIGPROF samples of the scan (dp.pfunc), 100 Hz. perf and
// gperftools are not installed here. Frames need a binary built with
// -fno-omit-frame-pointer. Off unless the flag is passed.
static int g_sample_fd = -1;
static void sample_on_sig(int) {
  void* stack[48];
  int n = backtrace(stack, 48);
  if (n < 1) n = 1;
  if (n > 48) n = 48;
  unsigned char nb = static_cast<unsigned char>(n);
  if (g_sample_fd >= 0) {
    ssize_t w = write(g_sample_fd, &nb, 1);
    if (w == 1) {
      const size_t bytes = static_cast<size_t>(n) * sizeof(void*);
      ssize_t left = static_cast<ssize_t>(bytes);
      const char* p = reinterpret_cast<const char*>(stack);
      while (left > 0) {
        ssize_t k = write(g_sample_fd, p, static_cast<size_t>(left));
        if (k <= 0) break;
        p += k;
        left -= k;
      }
    }
  }
}
static void sample_start(const char* path) {
  g_sample_fd = open(path, O_CREAT | O_TRUNC | O_WRONLY, 0644);
  if (g_sample_fd < 0) {
    std::perror(path);
    return;
  }
  const char magic[4] = {'S', 'A', 'M', 'P'};
  const uint64_t anchor = static_cast<uint64_t>(reinterpret_cast<uintptr_t>(&sample_on_sig));
  if (write(g_sample_fd, magic, 4) != 4 || write(g_sample_fd, &anchor, sizeof anchor) != static_cast<ssize_t>(sizeof anchor)) {
    std::perror(path);
    return;
  }
  struct sigaction sa;
  std::memset(&sa, 0, sizeof sa);
  sa.sa_handler = sample_on_sig;
  sigemptyset(&sa.sa_mask);
  sa.sa_flags = SA_RESTART;
  sigaction(SIGPROF, &sa, nullptr);
  struct itimerval it;
  std::memset(&it, 0, sizeof it);
  it.it_interval.tv_usec = 10000;
  it.it_value.tv_usec = 10000;
  setitimer(ITIMER_PROF, &it, nullptr);
  std::fprintf(stderr, "[sample] on %s 100 Hz\n", path);
}
static void sample_stop() {
  struct itimerval it;
  std::memset(&it, 0, sizeof it);
  setitimer(ITIMER_PROF, &it, nullptr);
  if (g_sample_fd >= 0) {
    close(g_sample_fd);
    g_sample_fd = -1;
  }
}

// ------------------------------------------------------------------- DP03
class DP03 {
 public:
  DP03(const std::string& seq, int beam, bool allow_pk, int loop_cap, bool score_qz,
       bool arm_score = false, const AccessRec* access = nullptr, int max_gap_scan = 200,
       int max_arm_len = 0, bool gap_step = false, const GateCfg* gate_cfg = nullptr,
       const Gate* gate = nullptr, const int* table_beam = nullptr,
       bool keep_runs = false)
      : E(seq), n(E.n), beam(beam > 0 ? beam : 0), allow_pk(allow_pk),
        cap(loop_cap), score_qz(score_qz), arm_score(arm_score), access(access),
        max_gap_scan(max_gap_scan),
        max_arm_len(max_arm_len <= 0 ? E.n : std::min(max_arm_len, E.n)), gap_step(gap_step),
        gate_(gate), keep_runs(keep_runs) {
    if (gate_cfg) gcfg = *gate_cfg;
    // table_beam[t] < 0: the global beam; 0: infinite
    for (int t = 0; t < BT_N; ++t)
      tbeam[t] = table_beam && table_beam[t] >= 0 ? table_beam[t] : this->beam;
    alloc();
  }

  DropLog* drops = nullptr;
  bool mem_scan = false;  // --mem-log: checkpoints at the quarters of this run's scan
  // --remove-state: states (state_key) erased at their table's pruning point, any beam
  const std::unordered_set<uint64_t>* remove_keys = nullptr;

  ArmStats arm_stats;
  GateStats gate_stats;

  double pfunc() {
    // The nested pass (no pseudoknots) walks mates and does not build the pair list.
    if (l2r_mode != 0 || allow_pk) ensure_pairs();
    for (int j = 0; j < n; ++j) {
      cur_j = j;
      cpu_stop_check(0);
      if (allow_pk && l2r_mode != 1) {
        fill_Qg(j);
        prune4(Qg[j], BT_QG, j, true, GATE_QG);
        fill_Qgls(j);
        prune4(Qgls[j], BT_QGLS, j);
        fill_Qgl(j);
        prune4(Qgl[j], BT_QGL, j, true, GATE_QGL);
        fill_Qgr(j);
        prune4(Qgr[j], BT_QGR, j);
        fill_Qp(j);
        prune2(Qp, Qp_at[j], j, BT_QP, nullptr, nullptr);
      }
      if (allow_pk && l2r_mode != 0) {
        l2r_step_pk(j);
        if (l2r_mode == 1) {
          const int pre = Qp_at[j].live;
          prune2(Qp, Qp_at[j], j, BT_QP, nullptr, nullptr);
          const int kept = Qp_at[j].live;
          qp_pre_n += pre;
          qp_kept_n += kept;
          qp_pre_max = std::max(qp_pre_max, pre);
          qp_kept_max = std::max(qp_kept_max, kept);
          qp_beam = tbeam[BT_QP];
          if (tbeam[BT_QP] > 0 && pre > tbeam[BT_QP]) ++qp_bind_n;
        }
      }
      if (l2r_mode != 0) l2r_done[LT_QP] = j;
      {
        const auto tq = std::chrono::steady_clock::now();
        fill_Qb(j);
        sec_qb += std::chrono::duration<double>(std::chrono::steady_clock::now() - tq).count();
      }
      prune2(Qb, Qb_at[j], j, BT_QB, nullptr, nullptr);
      if (l2r_mode != 0) l2r_done[LT_QB] = j;
      if (l2r_mode == 0 && !allow_pk) {
        const auto ti = std::chrono::steady_clock::now();
        lp_push_interior(j);
        sec_qb += std::chrono::duration<double>(std::chrono::steady_clock::now() - ti).count();
      } else if (l2r_mode != 0) {
        build_qlead(j);
      }
      fill_Qm(j);
      // Qm2 is dropped with Qm (its keys are a subset of Qm's), unless Qm2's own beam is
      // infinite; then it gets its own beam.  With equal beams the second call is a no-op.
      {
        const int pre = Qm_at[j].live;
        if (tbeam[BT_QM2] > 0) prune2(Qm, Qm_at[j], j, BT_QM, &Qm2, &Qm2_at[j]);
        else prune2(Qm, Qm_at[j], j, BT_QM, nullptr, nullptr);
        const int kept = Qm_at[j].live;
        qm_pre_n += pre;
        qm_kept_n += kept;
        qm_pre_max = std::max(qm_pre_max, pre);
        qm_kept_max = std::max(qm_kept_max, kept);
        qm_beam = tbeam[BT_QM];
        if (tbeam[BT_QM] > 0 && pre > tbeam[BT_QM]) ++qm_bind_n;
      }
      prune2(Qm2, Qm2_at[j], j, BT_QM2, nullptr, nullptr);
      if (l2r_mode == 0 && !allow_pk) lp_emit_multi(j);
      if ((allow_pk || chain_qmp) && E.stem_ml) {
        fill_Qmp(j);
        prune2(Qmp, Qmp_at[j], j, BT_QMP, nullptr, nullptr);
      }
      qm_reindex(j);
      if (l2r_mode != 0) l2r_done[LT_QM] = j;
      {
        const auto tq = std::chrono::steady_clock::now();
        if (qz_beam_fill) fill_Qz_kept(j);
        else fill_Qz(j);
        const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - tq).count();
        sec_qz += dt;
        if (prof_on) prof_qz_sec += dt;
      }
      if (l2r_mode != 0) l2r_done[LT_QZ] = j;
      {
        const auto tq = std::chrono::steady_clock::now();
        fill_C(j);
        sec_c += std::chrono::duration<double>(std::chrono::steady_clock::now() - tq).count();
      }
      advance_pk_inside(j);
      if (keep_runs && dz) save_runs(j);
      if (allow_pk && l2r_mode != 1) {
        fill_Qgrs(j);
        prune4(Qgrs[j], BT_QGRS, j);
      }
      if (allow_pk && l2r_mode != 0) l2r_step_end(j);
      if (l2r_mode != 0) note_col();
      if (g_beam_stats.f) g_beam_stats.log_Q(j + 1, C[j]);  // log Q_{1,t} = C[t-1] (1-based t=j+1)
      if (mem_scan && n >= 4 && (j + 1) % (n / 4) == 0 && j + 1 < n) {
        std::vector<MemItem> it;
        mem_items(it, "", false, true);
        mem_log(("scan_j" + std::to_string(j + 1)).c_str(), it);
      }
    }
    inside_done = true;
    logZ = n > 0 ? C[n - 1] : 0.0;
    if (n > 0 && l2r_mode == 0) {
      long long qb_n = 0, qm_n = 0;
      for (int j = 0; j < n; ++j) {
        qb_n += Qb_at[j].live;
        qm_n += Qm_at[j].live;
      }
      const double dn = static_cast<double>(n);
      std::fprintf(stderr,
                   "[p1rule] n %d lead %d hair %.3f int %.3f m2xP %.3f m2multi %.3f multi_ext %.3f "
                   "m_m_u %.3f m_p %.3f extC %.3f qz_cand %.3f qb_live %.3f qm_live %.3f\n",
                   n, ml_lead_cap, r_hair / dn, r_int / dn, r_m2xp / dn, r_m2multi / dn,
                   r_muext / dn, r_mmu / dn, r_mp / dn, r_extc / dn,
                   static_cast<double>(qz_cand) / dn, qb_n / dn, qm_n / dn);
      std::fprintf(stderr, "[p1sec] qb %.6f qm %.6f branch %.6f qz %.6f c %.6f\n", sec_qb, sec_qm,
                   sec_branch, sec_qz, sec_c);
      const double ms = 1000.0 / dn;
      const double rule = sec_hair + sec_int + sec_m2xp + sec_emit + sec_mclose + sec_mp + sec_c + sec_qz;
      const double pfunc_s = sec_qb + sec_qm + sec_branch + sec_qz + sec_c;
      std::fprintf(stderr,
                   "[p1ms] hair %.4f int %.4f m2xp %.4f emit %.4f mclose %.4f m %.4f c %.4f "
                   "qz %.4f sum %.4f pfunc %.4f\n",
                   sec_hair * ms, sec_int * ms, sec_m2xp * ms, sec_emit * ms, sec_mclose * ms,
                   sec_mp * ms, sec_c * ms, sec_qz * ms, rule * ms, pfunc_s * ms);
      std::fprintf(stderr, "[p1z] %.17g\n", logZ);
    }
    if (grammar_simple && n > 0) {
      const double dn = static_cast<double>(n);
      std::fprintf(stderr,
                   "[grammar] simple n %d fill %lld int %lld M1 %lld M2 %lld M3 %lld M4 %lld "
                   "M5 %lld M6 %lld M7 %lld M8 %lld phase %lld lead %lld overcap %lld\n",
                   n, gs_fill, it_int[0] + it_int[1], gs_m[0], gs_m[1], gs_m[2], gs_m[3], gs_m[4],
                   gs_m[5], gs_m[6], gs_m[7], gs_phase, gs_lead, gs_over);
      std::fprintf(stderr,
                   "[grammar] per_nt fill %.3f int %.3f M1 %.3f M2 %.3f M3 %.3f M4 %.3f M5 %.3f "
                   "M6 %.3f M7 %.3f M8 %.3f phase %.3f lead %.3f\n",
                   gs_fill / dn, (it_int[0] + it_int[1]) / dn, gs_m[0] / dn, gs_m[1] / dn,
                   gs_m[2] / dn, gs_m[3] / dn, gs_m[4] / dn, gs_m[5] / dn, gs_m[6] / dn,
                   gs_m[7] / dn, gs_phase / dn, gs_lead / dn);
      static const char* tn[6] = {"R^M", "R^M2", "R^Multi", "U^M", "U^M2", "U^Multi"};
      std::fprintf(stderr, "[grammar-kept] n %d", n);
      for (int t = 0; t < 6; ++t)
        std::fprintf(stderr, " %s %.3f", tn[t], gs_tab_kept[t] / dn);
      std::fprintf(stderr, "\n");
    }
    if (sw_check && l2r_mode != 0) sw_growth_check();
    if (!g_eq_path.empty() && l2r_mode != 0) g_eq_write();
    return logZ;
  }

  std::vector<std::tuple<int, int, double>> bpp(double cutoff) {
    outside();
    return collect_bpp(cutoff);
  }
  // Surviving inside Q^b, 1-based. Used by --dump-qb. Does not change Z.
  void dump_qb(FILE* f) const {
    for (int j = 0; j < n; ++j) {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        const int i = at.keys[r];
        if (i < 0 || at.vals[r] < HALF_NEG) continue;
        std::fprintf(f, "qb %d %d %.17g\n", i + 1, j + 1, at.vals[r]);
      }
    }
  }

  // One JSON line of per-state gate diagnostics (see write_gate_diag).
  std::string gate_diag(const Gate& gd, const std::vector<double>& taus, const std::string& id,
                        int span, double delta, int long_span);

  // Exact posterior weight of every stored state of the beamed tables (format at the
  // definition).
  void dump_states(FILE* f, const std::string& id, double min_w);

  // Rank of every candidate under the default and the --stem-est scores (format at the
  // definition); exact run.
  void dump_scores(FILE* f, const std::string& id, const StemEst& st, double min_w);
  void dump_ext(FILE* f, const std::string& id, const StemEst& st,
                const std::unordered_map<uint64_t, double>& w_of, int kmax);

  // --pr-oracle: the principled Out_hat of every 4-index state of this exact run with each
  // nested piece (H, Qzbar, Mbar, Ctx) nested or exact, against the exact outside (format
  // at the definition).
  struct PrVar {
    std::string name;
    bool hx, zx, mx, cx;  // H, Qzbar, Mbar, Ctx exact
    int fi;               // kept-pair floor: 0 = 1e-7, 1 = 0
    bool ml;              // multiloop stem steps (f, g, h)
    int lc;               // interior steps: unpaired per side <= min(lc, cap)
  };
  void pr_oracle_diag(StemEst& st, FILE* f, const std::string& id, const std::string& sets);
  void pr_hybrid_outside(const PrVar& v, const StemEst& st, const std::vector<char>& km,
                         const std::vector<std::vector<int>>& kr, std::vector<Tab4>* hb,
                         std::vector<double> own[2], const std::vector<size_t> off[5]);
  void pr_build_hbar(const std::vector<char>& km, const std::vector<std::vector<int>>& kl,
                     const std::vector<double>& mb, bool ml, int lc, std::vector<Tab4>& hb);

  // --stem-est tables from this finished nested-only run (StemEst above).
  void build_stem_est(int L, StemEst& st, bool mirror = false) const {
    if (!mates_ready) const_cast<DP03*>(this)->lp_build_mates();
    st.n = n;
    st.loop = L;
    st.C0 = C;
    std::vector<double>& S = st.S;
    S.assign(n + 1, NEG_INF);
    S[n] = 0.0;
    if (dz) {
      std::vector<double> P(n + 1, NEG_INF);
      P[n] = 0.0;
      for (int y = n - 1; y >= 0; --y) {
        double acc = NEG_INF;
        for (int j = y + 1; j < n; ++j) {
          double qb = Qb[ix(y, j)];
          if (qb > HALF_NEG && E.end_ok(y, j)) acc = log_add(acc, qb + E.at(y, j) + S[j + 1]);
        }
        P[y] = acc;
        double s = P[y];
        for (int t = y; t < n; ++t) s = log_add(s, E.run(y, t) + P[t + 1]);
        S[y] = s;
      }
    } else {
      for (int y = n - 1; y >= 0; --y) {
        double acc = S[y + 1];
        for (int j = y + 1; j < n; ++j) {
          double qb = Qb[ix(y, j)];
          if (qb > HALF_NEG && E.end_ok(y, j)) acc = log_add(acc, qb + E.at(y, j) + S[j + 1]);
        }
        S[y] = acc;
      }
    }
    const int lc = std::min(L, cap);
    std::vector<double> T(static_cast<size_t>(n) * n, NEG_INF);
    std::vector<stemidx::Entry> ents;
    for (int x = 0; x < n; ++x)
      for (int y = next_mate(x, x); y >= 0; y = next_mate(x, y)) {
        double acc = E.end_ok(x, y) ? E.at(x, y) + S[y + 1] : NEG_INF;
        for (int x0 = std::max(0, x - 1 - lc); x0 < x; ++x0)
          for (int y0 = y + 1; y0 <= std::min(n - 1, y + 1 + lc); ++y0) {
            double t0 = T[ix(x0, y0)];
            if (t0 < HALF_NEG) continue;
            double w = E.stem_interior(x0, x, y, y0);
            if (w < HALF_NEG) continue;
            acc = log_add(acc, w + t0);
          }
        T[ix(x, y)] = acc;
        if (acc > HALF_NEG && E.end_ok(x, y)) ents.push_back({x, y, (E.at(x, y) + acc) + C[y - 1]});
      }
    st.dropped = stem_sparse(ents);
    st.entries = ents.size();
    st.index = stemidx::PersistentSum(n, ents);
    st.mirror = mirror;
    if (!mirror) return;
    std::fill(T.begin(), T.end(), NEG_INF);
    std::vector<stemidx::Entry> ents2;
    for (int x = 0; x < n; ++x) {
      const double pre = x > 0 ? C[x - 1] : 0.0;
      for (int y = next_mate(x, x); y >= 0; y = next_mate(x, y)) {
        double acc = E.end_ok(x, y) ? E.at(x, y) + pre : NEG_INF;
        for (int x0 = std::max(0, x - 1 - lc); x0 < x; ++x0)
          for (int y0 = y + 1; y0 <= std::min(n - 1, y + 1 + lc); ++y0) {
            double t0 = T[ix(x0, y0)];
            if (t0 < HALF_NEG) continue;
            double w = E.stem_interior(x0, x, y, y0);
            if (w < HALF_NEG) continue;
            acc = log_add(acc, w + t0);
          }
        T[ix(x, y)] = acc;
        if (acc > HALF_NEG && E.end_ok(x, y))
          ents2.push_back({n - 1 - y, n - 1 - x, (E.at(x, y) + acc) + S[x + 1]});
      }
    }
    st.dropped2 = stem_sparse(ents2);
    st.entries2 = ents2.size();
    st.index2 = stemidx::PersistentSum(n, ents2);
  }

  // st.beta from this nested-only run (after build_stem_est, which sets st.S): runs the outside
  // pass and keeps the Qb outside value of every pair.
  void build_nested_beta(StemEst& st) {
    if (!mates_ready) lp_build_mates();
    outside();
    st.beta.assign(static_cast<size_t>(n) * n, NEG_INF);
    for (int x = 0; x < n; ++x)
      for (int y = next_mate(x, x); y >= 0; y = next_mate(x, y)) {
        double b = Qb_b[ix(x, y)];
        if (b < HALF_NEG && E.end_ok(x, y)) b = prefix_log(x) + E.at(x, y) + st.S[y + 1];
        st.beta[ix(x, y)] = b;
      }
  }

  // Nested-only run: also fill Qmp (stem_ml), the chain's mbar.  Set before pfunc().
  bool chain_qmp = false;

  // st's --stem-est-chain data from this nested-only run: runs the outside pass; kept pairs
  // have (Qb + Qb outside) - logZ >= log(floor).
  void build_chain(StemEst& st, double floor) {
    outside();
    const size_t nn = static_cast<size_t>(n) * n;
    st.chain = true;
    st.kept.assign(nn, 0);
    st.kept_right.assign(n, {});
    st.mbar.assign(nn, NEG_INF);
    st.n_kept = 0;
    const double lf = floor > 0.0 ? std::log(floor) : NEG_INF;
    for (int j = 0; j < n; ++j) {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        const int i = at.keys[r];
        if (i < 0) continue;
        const double b = Qb_b[ix(i, j)];
        if (at.vals[r] < HALF_NEG || b < HALF_NEG) continue;
        if (at.vals[r] + b - logZ < lf) continue;
        st.kept[ix(i, j)] = 1;
        st.kept_right[i].push_back(j);
        ++st.n_kept;
      }
    }
    for (std::vector<int>& kr : st.kept_right) std::sort(kr.begin(), kr.end());
    for (int j = 0; j < n; ++j) {
      const OMap& at = qmp_at(j);
      for (size_t r = 0; r < at.slots(); ++r)
        if (at.keys[r] >= 0) st.mbar[ix(at.keys[r], j)] = at.vals[r];
    }
  }

  // Closed-value outsides of Qm / Qm2 with dangles (before run_back_at); set before outside().
  bool keep_closed_b = false;
  SpMat QmC_b, Qm2C_b;

  // st's --stem-est-principled data from this nested-only run (chain_qmp set; runs the outside
  // pass): prS, prQz, prMbar, kept pairs (P >= floor), prBml; nested: the outsides of Qb / Qm /
  // Qm2 (closed values).  Bml[d][e] = (sum of
  //   (a) Qm_b[i][e] w3^(d-i) run(i, d-1)          over cells (i, e), d - i <= ml_left,
  //   (b) Qm_b[i][e] Qm[i][d-1]                     over cells (i, d-1),
  // times alpha2^ml_pk_k: rev_branch's pseudoknot branch) + (c) the PKM closing of every pair
  // i.j' (i < d, j' > e): Qb_b[i][j'] alpha1 alpha2^(1+ml_pk_k) at_close(i, j') w3^(j'-1-e)
  // run(e+1, j'-1) [w3^(d-i-1) run(i+1, d-1) + Qm[i+1][d-1]] (fill_Qb's PKM terms); beta1m
  // is left to Ctx.  check: every cell recomputed directly (n <= 200).
  void build_principled(StemEst& st, double floor, bool nested, bool check,
                        double partner_floor = -1.0) {
    outside();
    const size_t nn = static_cast<size_t>(n) * n;
    st.n = n;
    st.pr = true;
    std::vector<double>& S = st.prS;
    S.assign(n + 1, NEG_INF);
    S[n] = 0.0;
    if (dz) {
      std::vector<double> P(n + 1, NEG_INF);
      P[n] = 0.0;
      for (int y = n - 1; y >= 0; --y) {
        double acc = NEG_INF;
        for (int j = y + 1; j < n; ++j) {
          double qb = Qb[ix(y, j)];
          if (qb > HALF_NEG && E.end_ok(y, j)) acc = log_add(acc, qb + E.at(y, j) + S[j + 1]);
        }
        P[y] = acc;
        double s = P[y];
        for (int t = y; t < n; ++t) s = log_add(s, E.run(y, t) + P[t + 1]);
        S[y] = s;
      }
    } else {
      for (int y = n - 1; y >= 0; --y) {
        double acc = S[y + 1];
        for (int j = y + 1; j < n; ++j) {
          double qb = Qb[ix(y, j)];
          if (qb > HALF_NEG && E.end_ok(y, j)) acc = log_add(acc, qb + E.at(y, j) + S[j + 1]);
        }
        S[y] = acc;
      }
    }
    // --qzbar-old: always dense-copy nested Qz into prQz (snap/incperf), even if qz_beam_fill.
    if (!qz_beam_fill || qzbar_old) {
      st.prQz.assign(nn, NEG_INF);
      for (int i = 0; i < n; ++i)
        for (int j = i; j < n; ++j) st.prQz[ix(i, j)] = Qz[ix(i, j)];
      if (qzbar_old) {
        st.qz_on = 0;
        st.qz_col.clear();
        st.qz_col.shrink_to_fit();
      }
    }
    if (!st.skip_qbin) st.prQbIn = Qb.dense();
    st.prMbar.assign(nn, NEG_INF);
    if (p2_sweep) {
      st.prMbarAt.assign(n, {});
      st.prMbarRight.assign(n, {});
    }
    for (int j = 0; j < n; ++j) {
      const OMap& at = qmp_at(j);
      for (size_t r = 0; r < at.slots(); ++r) {
        const int i = at.keys[r];
        if (i < 0) continue;
        st.prMbar[ix(i, j)] = at.vals[r];
        if (p2_sweep) {
          st.prMbarAt[j].push_back(i);
          st.prMbarRight[i].push_back(j);
        }
      }
    }
    if (p2_sweep) {
      st.pr_mbar_n = 0;
      st.pr_mbar_max = 0;
      for (int j = 0; j < n; ++j) {
        st.pr_mbar_n += static_cast<long long>(st.prMbarAt[j].size());
        st.pr_mbar_max = std::max(st.pr_mbar_max, static_cast<int>(st.prMbarAt[j].size()));
        std::sort(st.prMbarRight[j].begin(), st.prMbarRight[j].end());
      }
      if (beam > 0) {
        auto stat = [&](const std::vector<OMap>& tabs) {
          long long c = 0;
          int m = 0;
          for (int j = 0; j < n; ++j) {
            int k = 0;
            for (size_t r = 0; r < tabs[j].slots(); ++r)
              if (tabs[j].keys[r] >= 0) ++k;
            c += k;
            m = std::max(m, k);
          }
          return std::pair<long long, int>(c, m);
        };
        const auto qb = stat(Qb_at), qm = stat(Qm_at), q2 = stat(Qm2_at);
        std::fprintf(stderr,
                     "[pre-beam] b %d key %s Qb %lld max %d Qm %lld max %d Qm2 %lld max %d\n",
                     beam, pre_key_inside ? "inside" : "prefix", qb.first, qb.second, qm.first,
                     qm.second, q2.first, q2.second);
      }
    }
    if (prof_on && !(p2_sweep && beam > 0)) {
      int qm_max = 0;
      for (int j = 0; j < n; ++j) {
        int k = 0;
        const OMap& at = qmp_at(j);
        for (size_t r = 0; r < at.slots(); ++r)
          if (at.keys[r] >= 0) ++k;
        qm_max = std::max(qm_max, k);
      }
      std::fprintf(stderr, "[qm-col] max %d\n", qm_max);
    }
    if (!st.skip_kept) st.prKept.assign(nn, 0);
    st.prKeptRight.assign(n, {});
    st.pr_n_kept = 0;
    const double lf = floor > 0.0 ? std::log(floor) : NEG_INF;
    st.pr_partner = partner_floor >= 0.0;
    const double lfp = partner_floor > 0.0 ? std::log(partner_floor) : NEG_INF;
    st.prKeptRightP.assign(st.pr_partner ? n : 0, {});
    st.pr_n_partner = 0;
    if (st.want_lp) st.prLP.assign(nn, NEG_INF);
    for (int j = 0; j < n; ++j) {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        const int i = at.keys[r];
        if (i < 0) continue;
        const double b = Qb_b[ix(i, j)];
        if (at.vals[r] < HALF_NEG || b < HALF_NEG) continue;
        if (st.want_lp) st.prLP[ix(i, j)] = at.vals[r] + b - logZ;
        if (st.pr_partner && at.vals[r] + b - logZ >= lfp) {
          st.prKeptRightP[i].push_back(j);
          ++st.pr_n_partner;
        }
        if (at.vals[r] + b - logZ < lf) continue;
        if (!st.skip_kept) st.prKept[ix(i, j)] = 1;
        st.prKeptRight[i].push_back(j);
        ++st.pr_n_kept;
      }
    }
    for (std::vector<int>& kr : st.prKeptRight) std::sort(kr.begin(), kr.end());
    for (std::vector<int>& kr : st.prKeptRightP) std::sort(kr.begin(), kr.end());
    st.prBml.clear();
    if (g_skip_bml) {
      if (nested) {
        st.pr_nested = true;
        st.prQb_b = Qb_b.dense();
        st.prQm_b = dz ? QmC_b.dense() : Qm_b.dense();
        st.prQm2_b = dz ? Qm2C_b.dense() : Qm2_b.dense();
      }
      pass2_apply(st);
      return;
    }

    const double w2 = E.w_alpha2, w3 = E.w_alpha3, kw = E.ml_pk_k * w2;
    std::vector<double> Bab(nn, NEG_INF);
    for (int e = 0; e < n; ++e) {
      const OMap& amj = Qm_at[e];
      for (size_t r = 0; r < amj.slots(); ++r) {
        const int i = amj.keys[r];
        if (i < 0) continue;
        const double mb = Qm_b[ix(i, e)];
        if (mb < HALF_NEG) continue;
        for (int d = i; d <= e && d - i <= E.ml_left && lead_gap_ok(d - i); ++d) {
          double w = w3 * (d - i);
          if (dz) w += E.run(i, d - 1);
          Bab[ix(d, e)] = log_add(Bab[ix(d, e)], mb + w);
        }
      }
    }
    for (int d = 1; d < n; ++d) {
      const OMap& amd = Qm_at[d - 1];
      for (size_t r = 0; r < amd.slots(); ++r) {
        const int i = amd.keys[r];
        if (i < 0) continue;
        const double qm = amd.vals[r];
        for (int e = d; e < n; ++e) {
          const double mb = Qm_b[ix(i, e)];
          if (mb > HALF_NEG) Bab[ix(d, e)] = log_add(Bab[ix(d, e)], mb + qm);
        }
      }
    }
    // (c): R[i][e] = sum over closings i.j' > e of their weight up to the PKM run end e
    std::vector<double> R(nn, NEG_INF);
    const double wc = E.w_alpha1 + (1.0 + E.ml_pk_k) * w2;
    for (int j = 0; j < n; ++j) {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        const int i = at.keys[r];
        if (i < 0 || !E.end_ok(i, j)) continue;
        const double bb = Qb_b[ix(i, j)];
        if (bb < HALF_NEG) continue;
        const double base = bb + wc + E.at_close(i, j);
        for (int e = i + 1; e <= j - 1; ++e) {
          double w = base + w3 * (j - 1 - e);
          if (dz) w += E.run(e + 1, j - 1);
          R[ix(i, e)] = log_add(R[ix(i, e)], w);
        }
      }
    }
    st.prBml.assign(nn, NEG_INF);
    std::vector<double> Lc(n, NEG_INF);
    for (int d = 1; d < n; ++d) {
      for (int i = 0; i < d; ++i) {
        double l = w3 * (d - i - 1);
        if (dz) l += E.run(i + 1, d - 1);
        if (i + 1 <= d - 1 && Qm[ix(i + 1, d - 1)] > HALF_NEG) l = log_add(l, Qm[ix(i + 1, d - 1)]);
        Lc[i] = l;
      }
      for (int e = d; e < n; ++e) {
        double c = NEG_INF;
        for (int i = 0; i < d; ++i) {
          const double rv = R[ix(i, e)];
          if (rv > HALF_NEG) c = log_add(c, rv + Lc[i]);
        }
        const double ab = Bab[ix(d, e)];
        st.prBml[ix(d, e)] = log_add(ab > HALF_NEG ? ab + kw : NEG_INF, c);
      }
    }
    if (nested) {
      st.pr_nested = true;
      st.prQb_b = Qb_b.dense();
      st.prQm_b = dz ? QmC_b.dense() : Qm_b.dense();
      st.prQm2_b = dz ? Qm2C_b.dense() : Qm2_b.dense();
    }
    if (check && n <= 200) bml_check(st);
    pass2_apply(st);
  }

  // Pass 2 limits. Every flag defaults off, and then this function returns without touching tables.
  // Q^z_pre(k,l): DP03 Q^z split on the last base, Q^p dropped. Kept starts ranked by
  // C_pre(k-1) * Q^z_pre (key1) or that times Pr(k-1 is a pair right end) (key2).
  // Pair sets: P1 = Q^b states that survived pass 1; own:b_P = top b_P left ends per right
  // end by the inside Q^b.
  // Outside push of gamma (growth) and mu (V) along the nested Qb / Qm / Qm2 hyperedges.
  // Pseudoknots are off on this pre-pass. Right ends run n-1..0. At j the order is the
  // reverse of the inside: extension from j+1, interior from kept outer pairs, the exterior
  // edge, the Qm/Qm2 branch split over kept states, then the pair pushes into Qm2.
  // Each map is cut before it pushes. Qm and Qm2 maps are freed at the end of j.
  struct LeadEnt {
    int x0;
    double w;
  };
  void sw_hyper_push(StemEst& st) {
    using Ga = std::unordered_map<uint64_t, double>;
    struct M {
      Ga g, v, b;
      bool empty() const { return g.empty() && v.empty() && b.empty(); }
    };
    // Nested pre-pass has l2r_mode==0 (scan's --pk-l2r does not apply here). Multi edges E3/E4
    // multi fire on that path. --g-sweep states turns G on beside V.
    const bool do_g = (g_sweep == G_SWEEP_STATES);
    const bool do_v = p2_vsweep;
    const bool do_b = sw_check && p2_vsweep && p2_v_cap <= 0 && g_ntr_K <= 0;
    const int Kg = g_ntr_K, Kv = p2_v_cap;
    const double a1 = E.w_alpha1, a2 = E.w_alpha2, a3 = E.w_alpha3;
    long long pushes = 0, st_g = 0, st_v = 0;
    int qm_max = 0;
    struct Hist {
      std::vector<uint32_t> c;
      long long n = 0, sum = 0;
      void add(size_t s) {
        if (c.size() <= s) c.resize(s + 1);
        ++c[s];
        ++n;
        sum += static_cast<long long>(s);
      }
      double mean() const { return n ? static_cast<double>(sum) / static_cast<double>(n) : 0.0; }
      long long max() const { return n && !c.empty() ? static_cast<long long>(c.size() - 1) : 0; }
      long long p99() const {
        if (!n) return 0;
        long long target = (n * 99 + 99) / 100;
        if (target < 1) target = 1;
        long long acc = 0;
        for (size_t i = 0; i < c.size(); ++i) {
          acc += c[i];
          if (acc >= target) return static_cast<long long>(i);
        }
        return c.empty() ? 0 : static_cast<long long>(c.size() - 1);
      }
    } g_at, v_at, g_pre, v_pre, g_kept, v_kept;
    auto scale_into = [&](M& d, const M& s, double a) {
      if (a < HALF_NEG || s.empty()) return;
      auto one = [&](Ga& dd, const Ga& ss, double& sec, long long* cnt, Hist* hist) {
        if (ss.empty()) return;
        if (hist) hist->add(ss.size());
        const auto t0 = std::chrono::steady_clock::now();
        for (const auto& [k, v] : ss) {
          if (v < HALF_NEG) continue;
          const double w = a + v;
          auto it = dd.find(k);
          if (it == dd.end()) dd.emplace(k, w);
          else it->second = log_add(it->second, w);
        }
        sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (cnt) ++*cnt;
      };
      if (do_g) one(d.g, s.g, g_merge_sec, &g_pushes, &g_at);
      if (do_v) one(d.v, s.v, v_merge_sec, &v_pushes, &v_at);
      if (do_b) {
        double ignore = 0;
        one(d.b, s.b, ignore, nullptr, nullptr);
      }
    };
    auto push_to = [&](M& d, const M& s, double a) {
      if (a < HALF_NEG || s.empty()) return;
      scale_into(d, s, a);
      ++pushes;
    };
    auto cap_ga = [&](Ga& m) {
      const long long pre = static_cast<long long>(m.size());
      g_pre.add(m.size());
      ++g_pre_n;
      g_pre_sum += pre;
      if (pre > g_pre_max) g_pre_max = pre;
      const auto tc0 = std::chrono::steady_clock::now();
      // --g-sweep states: always entries cut by C_pre * G * S_pre (bG = Kg).
      if (Kg > 0 && static_cast<int>(m.size()) > Kg) {
        struct E {
          uint64_t k;
          double w, s;
        };
        std::vector<E> v;
        v.reserve(m.size());
        for (const auto& [k, w] : m) {
          if (w < HALF_NEG) continue;
          const int p = k_hi(k), q = k_lo(k);
          double s = w;
          if (p > 0 && !st.prC.empty()) s += st.prC[static_cast<size_t>(p - 1)];
          if (q + 1 < 0 || q + 1 >= static_cast<int>(st.prS.size())) continue;
          const double sb = st.prS[static_cast<size_t>(q + 1)];
          if (sb < HALF_NEG) continue;
          s += sb;
          v.push_back({k, w, s});
        }
        if (static_cast<int>(v.size()) > Kg) {
          std::nth_element(v.begin(), v.begin() + Kg, v.end(), [](const E& a, const E& b) {
            if (a.s != b.s) return a.s > b.s;
            return a.k < b.k;
          });
          v.resize(static_cast<size_t>(Kg));
        }
        m.clear();
        for (const E& e : v) m.emplace(e.k, e.w);
      }
      g_cut_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - tc0).count();
      const long long post = static_cast<long long>(m.size());
      ++g_post_n;
      g_post_sum += post;
      if (post > g_post_max) g_post_max = post;
    };
    auto cap_mu = [&](Ga& m) {
      const long long pre = static_cast<long long>(m.size());
      v_pre.add(m.size());
      ++v_pre_n;
      v_pre_sum += pre;
      if (pre > v_pre_max) v_pre_max = pre;
      const auto tc0 = std::chrono::steady_clock::now();
      if (Kv > 0 && static_cast<int>(m.size()) > Kv) {
        struct E {
          uint64_t k;
          double w, s;
        };
        std::vector<E> v;
        v.reserve(m.size());
        for (const auto& [k, w] : m) {
          if (w < HALF_NEG) continue;
          const int b = k_hi(k);
          double s = w;
          if (b > 0 && !st.prC.empty()) s += st.prC[static_cast<size_t>(b - 1)];
          v.push_back({k, w, s});
        }
        if (static_cast<int>(v.size()) > Kv) {
          std::nth_element(v.begin(), v.begin() + Kv, v.end(), [](const E& a, const E& b) {
            if (a.s != b.s) return a.s > b.s;
            return a.k < b.k;
          });
          v.resize(static_cast<size_t>(Kv));
        }
        m.clear();
        for (const E& e : v) m.emplace(e.k, e.w);
      }
      v_cut_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - tc0).count();
      const long long post = static_cast<long long>(m.size());
      ++v_post_n;
      v_post_sum += post;
      if (post > v_post_max) v_post_max = post;
    };
    auto cap_m = [&](M& s) {
      if (do_g && !s.g.empty()) {
        cap_ga(s.g);
        if (!s.g.empty()) {
          ++st_g;
          g_kept.add(s.g.size());
        }
      }
      if (do_v && !s.v.empty()) {
        cap_mu(s.v);
        if (!s.v.empty()) {
          ++st_v;
          v_kept.add(s.v.size());
        }
      }
    };
    auto add1 = [&](Ga& d, uint64_t k, double w) {
      if (w < HALF_NEG) return;
      auto it = d.find(k);
      if (it == d.end()) d.emplace(k, w);
      else it->second = log_add(it->second, w);
    };
    // A pruned target is not zeroed by the edge coefficient: dQ_s/dQ_c does not contain Q_c.
    long long hk_map[3] = {}, hk_bad[3] = {}, hk_e[12] = {}, hk_eb[12] = {};
    auto kept_at = [&](const OMap& at, int i) { return i >= 0 && at.has(i); };
    auto note_push = [&](int e, bool ok) {
      ++hk_e[e];
      if (!ok) ++hk_eb[e];
    };
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<double> Cb(static_cast<size_t>(n), NEG_INF);
    if (n > 0) Cb[static_cast<size_t>(n - 1)] = 0.0;
    double bc1 = NEG_INF, bc2 = NEG_INF;
    std::vector<std::unordered_map<int, M>> qm(static_cast<size_t>(n)), qm2(static_cast<size_t>(n));
    std::vector<std::unordered_map<int, M>> u1(static_cast<size_t>(n)), u2(static_cast<size_t>(n));
    std::vector<std::unordered_map<int, M>> u12(static_cast<size_t>(n)), u22(static_cast<size_t>(n));
    std::vector<std::unordered_map<int, M>> lpm(static_cast<size_t>(n));
    std::unordered_map<uint64_t, M> pairs;
    auto run_mix = [&](M& closed, const M* u1n, const M* u2n, int j, double w, M& ou1, M& ou2) {
      M p = closed;
      if (u1n) scale_into(p, *u1n, w);
      M a, b;
      scale_into(a, closed, E.s1w[j]);
      if (u2n) scale_into(a, *u2n, E.d5w[j] + w);
      scale_into(b, closed, E.d3w[j]);
      if (u2n) scale_into(b, *u2n, w);
      if ((u1n && !u1n->empty()) || (u2n && !u2n->empty())) ++pushes;
      closed = std::move(p);
      ou1 = std::move(a);
      ou2 = std::move(b);
    };
    for (int j = n - 1; j >= 0; --j) {
      qm_max = std::max(qm_max, static_cast<int>(Qm_at[j].live));
      if (dz) {
        double bP, b1, b2;
        run_back(Cb[static_cast<size_t>(j)], bc1, bc2, j, 0.0, bP, b1, b2);
        Cb[static_cast<size_t>(j)] = bP;
        bc1 = b1;
        bc2 = b2;
      } else if (j + 1 < n) {
        Cb[static_cast<size_t>(j)] = log_add(Cb[static_cast<size_t>(j)], Cb[static_cast<size_t>(j + 1)]);
      }
      auto extend_col = [&](int edge, std::unordered_map<int, M>& col, const OMap& at,
                            std::unordered_map<int, M>& U1, std::unordered_map<int, M>& U2,
                            std::unordered_map<int, M>* U1n, std::unordered_map<int, M>* U2n) {
        std::unordered_map<int, M> nu1, nu2;
        for (size_t r = 0; r < at.slots(); ++r) {
          const int i = at.keys[r];
          if (i < 0) continue;
          note_push(edge, kept_at(at, i));
          const M* a = nullptr;
          const M* b = nullptr;
          if (U1n) {
            auto it = U1n->find(i);
            if (it != U1n->end()) a = &it->second;
          }
          if (U2n) {
            auto it = U2n->find(i);
            if (it != U2n->end()) b = &it->second;
          }
          M& closed = col[i];
          if (closed.empty() && !a && !b) continue;
          run_mix(closed, a, b, j, a3, nu1[i], nu2[i]);
        }
        U1.swap(nu1);
        U2.swap(nu2);
      };
      if (dz) {
        extend_col(0, qm[static_cast<size_t>(j)], Qm_at[j], u1[static_cast<size_t>(j)],
                   u2[static_cast<size_t>(j)], j + 1 < n ? &u1[static_cast<size_t>(j + 1)] : nullptr,
                   j + 1 < n ? &u2[static_cast<size_t>(j + 1)] : nullptr);
        extend_col(1, qm2[static_cast<size_t>(j)], Qm2_at[j], u12[static_cast<size_t>(j)],
                   u22[static_cast<size_t>(j)], j + 1 < n ? &u12[static_cast<size_t>(j + 1)] : nullptr,
                   j + 1 < n ? &u22[static_cast<size_t>(j + 1)] : nullptr);
      } else if (j + 1 < n) {
        for (size_t r = 0; r < Qm_at[j].slots(); ++r) {
          const int i = Qm_at[j].keys[r];
          if (i < 0) continue;
          auto it = qm[static_cast<size_t>(j + 1)].find(i);
          if (it != qm[static_cast<size_t>(j + 1)].end()) {
            note_push(0, kept_at(Qm_at[j], i));
            push_to(qm[static_cast<size_t>(j)][i], it->second, a3);
          }
        }
        for (size_t r = 0; r < Qm2_at[j].slots(); ++r) {
          const int i = Qm2_at[j].keys[r];
          if (i < 0) continue;
          auto it = qm2[static_cast<size_t>(j + 1)].find(i);
          if (it != qm2[static_cast<size_t>(j + 1)].end()) {
            note_push(1, kept_at(Qm2_at[j], i));
            push_to(qm2[static_cast<size_t>(j)][i], it->second, a3);
          }
        }
      }
      if (l2r_mode == 0 && !lpP.empty() && !lpM.empty()) {
        const int L = lead_span();
        const OMap& src = lpP[static_cast<size_t>(j)];
        for (size_t s = 0; s < src.slots(); ++s) {
          const int ib = src.keys[s];
          if (ib < 0 || !Qm2_at[j].has(ib)) continue;
          const int p_lo = std::max(0, ib - L);
          for (int p = ib - 1; p >= p_lo; --p) {
            const int q = next_mate(p, j);
            if (q < 0 || q >= n) continue;
            const int k = p * 3 + lp_kind(q - j - 1);
            if (!lpM[static_cast<size_t>(q)].has(k)) continue;
            auto it = lpm[static_cast<size_t>(q)].find(k);
            if (it == lpm[static_cast<size_t>(q)].end() || it->second.empty()) continue;
            const double factor = unpaired_gap(p + 1, ib - 1) + unpaired_gap(j + 1, q - 1);
            note_push(10, true);
            push_to(qm2[static_cast<size_t>(j)][ib], it->second, factor);
          }
        }
      }
      const int qhi = std::min(n, j + 2 + cap);
      for (size_t r = 0; r < Qb_at[j].slots(); ++r) {
        const int d = Qb_at[j].keys[r];
        if (d < 0) continue;
        M& inner = pairs[key2(d, j)];
        if (l2r_mode == 0) {
          const int p_lo = std::max(0, d - 1 - cap);
          for (int p = d - 1; p >= p_lo; --p) {
            const int budget = cap - (d - p - 1);
            if (budget < 0) break;
            for (int q = next_mate(p, j); q != -1 && q - j - 1 <= budget; q = next_mate(p, q)) {
              if (!Qb_at[q].has(p)) continue;
              const double w = E.interior(p, d, j, q);
              if (w < HALF_NEG) continue;
              auto it = pairs.find(key2(p, q));
              note_push(2, it != pairs.end());
              if (it != pairs.end()) push_to(inner, it->second, w);
            }
          }
          continue;
        }
        for (int q = j + 1; q < qhi; ++q) {
          for (size_t s = 0; s < Qb_at[q].slots(); ++s) {
            const int i = Qb_at[q].keys[s];
            if (i < 0 || d <= i) continue;
            if ((d - i - 1) + (q - j - 1) > cap) continue;
            const double w = E.interior(i, d, j, q);
            if (w < HALF_NEG) continue;
            auto it = pairs.find(key2(i, q));
            if (it != pairs.end()) {
              note_push(2, kept_at(Qb_at[j], d) && kept_at(Qb_at[q], i));
              push_to(inner, it->second, w);
            }
          }
        }
      }
      for (size_t r = 0; r < Qb_at[j].slots(); ++r) {
        const int i = Qb_at[j].keys[r];
        if (i < 0 || !E.end_ok(i, j)) continue;
        note_push(3, kept_at(Qb_at[j], i));
        const double at = E.at(i, j);
        const double cb = Cb[static_cast<size_t>(j)];
        if (do_b && cb > HALF_NEG && at > HALF_NEG) add1(pairs[key2(i, j)].b, key2(i, 0), cb + at);
        const double left = i == 0 ? 0.0 : C[static_cast<size_t>(i - 1)];
        const double qb = Qb[ix(i, j)];
        if (i > 0 && cb > HALF_NEG && qb > HALF_NEG)
          Cb[static_cast<size_t>(i - 1)] = log_add(Cb[static_cast<size_t>(i - 1)], (cb + qb) + at);
        (void)left;
      }
      for (auto& kv : qm[static_cast<size_t>(j)]) cap_m(kv.second);
      for (auto& kv : qm2[static_cast<size_t>(j)]) cap_m(kv.second);
      for (size_t r = 0; r < Qb_at[j].slots(); ++r) {
        const int d = Qb_at[j].keys[r];
        if (d < 0 || !E.end_ok(d, j)) continue;
        const double w_at = E.at(d, j);
        const double w2 = a2 + w_at;
        const double obj = Qb[ix(d, j)];
        M& dst = pairs[key2(d, j)];
        for (size_t s = 0; s < Qm_at[j].slots(); ++s) {
          const int i = Qm_at[j].keys[s];
          if (i < 0 || i > d || d - i > E.ml_left || !lead_gap_ok(d - i)) continue;
          const LeadEnt* lead = nullptr;
          if (l2r_mode != 0 && d < static_cast<int>(qlead.size())) {
            const auto& col = qlead[static_cast<size_t>(d)];
            auto it = std::lower_bound(col.begin(), col.end(), i,
                                       [](const LeadEnt& e, int x) { return e.x0 < x; });
            if (it != col.end() && it->x0 == i) lead = &*it;
          }
          // Pass 1's M = P has no left unpaired run; the fragment starts on the pair.
          // --pk-l2r off is the DP03 reference and keeps the leading-run window.
          if (l2r_mode == 0 && !allow_pk) {
            if (i != d) continue;
          } else if (!lead) {
            continue;
          }
          const double w = (l2r_mode == 0 && !allow_pk) ? w2 : w2 + lead->w;
          auto it = qm[static_cast<size_t>(j)].find(i);
          if (it != qm[static_cast<size_t>(j)].end()) {
            note_push(4, kept_at(Qb_at[j], d) && kept_at(Qm_at[j], i));
            push_to(dst, it->second, w);
          }
        }
        if (d == 0) continue;
        for (size_t s = 0; s < Qm_at[d - 1].slots(); ++s) {
          const int i = Qm_at[d - 1].keys[s];
          if (i < 0) continue;
          const double qm_in = Qm[ix(i, d - 1)];
          if (qm_in < HALF_NEG) continue;
          auto it = qm[static_cast<size_t>(j)].find(i);
          auto it2 = qm2[static_cast<size_t>(j)].find(i);
          if (it != qm[static_cast<size_t>(j)].end()) {
            note_push(5, kept_at(Qb_at[j], d) && kept_at(Qm_at[j], i) && kept_at(Qm_at[d - 1], i));
            push_to(dst, it->second, qm_in + w2);
            if (obj > HALF_NEG) {
              note_push(6, kept_at(Qm_at[d - 1], i) && kept_at(Qm_at[j], i));
              push_to(qm[static_cast<size_t>(d - 1)][i], it->second, obj + w2);
            }
          }
          if (it2 != qm2[static_cast<size_t>(j)].end()) {
            note_push(7, kept_at(Qb_at[j], d) && kept_at(Qm2_at[j], i) && kept_at(Qm_at[d - 1], i));
            push_to(dst, it2->second, qm_in + w2);
            if (obj > HALF_NEG) {
              note_push(8, kept_at(Qm_at[d - 1], i) && kept_at(Qm2_at[j], i));
              push_to(qm[static_cast<size_t>(d - 1)][i], it2->second, obj + w2);
            }
          }
        }
      }
      const double w_ml = a1 + a2;
      for (size_t r = 0; r < Qb_at[j].slots(); ++r) {
        const int i = Qb_at[j].keys[r];
        if (i < 0) continue;
        M& m = pairs[key2(i, j)];
        if (E.end_ok(i, j)) {
          if (do_g) add1(m.g, key2(i, j), E.at(i, j));
          if (do_v) {
            const double sb = st.prS[static_cast<size_t>(j + 1)];
            if (sb > HALF_NEG) add1(m.v, key2(i, 0), E.at(i, j) + sb);
          }
        }
        cap_m(m);
        if (l2r_mode != 0 && E.end_ok(i, j) && i + 1 <= j - 1 && Qm2_at[j - 1].has(i + 1)) {
          note_push(9, kept_at(Qb_at[j], i) && kept_at(Qm2_at[j - 1], i + 1));
          push_to(qm2[static_cast<size_t>(j - 1)][i + 1], m, w_ml + E.at_close(i, j));
        }
      }
      if (l2r_mode == 0 && !lpM.empty()) {
        const OMap& col = lpM[static_cast<size_t>(j)];
        for (size_t r = 0; r < col.slots(); ++r) {
          const int k = col.keys[r];
          if (k < 0) continue;
          const int p = k / 3;
          const int kind = k % 3;
          M& mu = lpm[static_cast<size_t>(j)][k];
          if (Qb_at[j].has(p) && E.end_ok(p, j)) {
            auto it = pairs.find(key2(p, j));
            if (it != pairs.end() && !it->second.empty()) {
              note_push(9, true);
              push_to(mu, it->second, E.w_alpha1 + E.w_alpha2 + E.at_close(p, j));
            }
          }
          int qn, nk;
          double delta;
          if (!lp_step(p, j, kind, qn, nk, delta)) continue;
          const int kn = p * 3 + nk;
          if (qn < 0 || qn >= n || !lpM[static_cast<size_t>(qn)].has(kn)) continue;
          auto it = lpm[static_cast<size_t>(qn)].find(kn);
          if (it == lpm[static_cast<size_t>(qn)].end() || it->second.empty()) continue;
          note_push(11, true);
          push_to(mu, it->second, delta);
        }
      }
      auto tally = [&](const std::unordered_map<int, M>& col, const OMap& at, int slot) {
        for (const auto& kv : col) {
          ++hk_map[slot];
          if (!at.has(kv.first)) ++hk_bad[slot];
        }
      };
      tally(qm[static_cast<size_t>(j)], Qm_at[j], 0);
      tally(qm2[static_cast<size_t>(j)], Qm2_at[j], 1);
      std::unordered_map<int, M>().swap(qm[static_cast<size_t>(j)]);
      std::unordered_map<int, M>().swap(qm2[static_cast<size_t>(j)]);
      if (j + 1 < n) {
        std::unordered_map<int, M>().swap(u1[static_cast<size_t>(j + 1)]);
        std::unordered_map<int, M>().swap(u2[static_cast<size_t>(j + 1)]);
        std::unordered_map<int, M>().swap(u12[static_cast<size_t>(j + 1)]);
        std::unordered_map<int, M>().swap(u22[static_cast<size_t>(j + 1)]);
      }
    }
    st.prGamma.clear();
    st.prMu.clear();
    st.prMuFloor.clear();
    for (auto& [k, m] : pairs) {
      const int px = k_hi(k), py = k_lo(k);
      ++hk_map[2];
      if (py < 0 || py >= n || !Qb_at[py].has(px)) ++hk_bad[2];
      if (do_g && !m.g.empty()) st.prGamma.emplace(k, std::move(m.g));
      if (do_v && !m.v.empty()) {
        double mn = NEG_INF;
        bool any = false;
        for (const auto& [bk, w] : m.v) {
          if (w < HALF_NEG) continue;
          if (!any || w < mn) {
            mn = w;
            any = true;
          }
        }
        if (any) st.prMuFloor.emplace(k, mn);
        st.prMu.emplace(k, std::move(m.v));
      }
    }
    if (!g_mu_dump.empty()) {
      std::ofstream md(g_mu_dump);
      for (const auto& [k, mv] : st.prMu) {
        const int x = k_hi(k), y = k_lo(k);
        for (const auto& [bk, w] : mv) {
          if (!(w > HALF_NEG)) continue;
          md << "V " << x << " " << y << " " << k_hi(bk) << " " << w << "\n";
        }
      }
      for (const auto& [k, mg] : st.prGamma) {
        const int x = k_hi(k), y = k_lo(k);
        std::unordered_map<int, double> byq;
        for (const auto& [bk, w] : mg) {
          if (!(w > HALF_NEG)) continue;
          const int p = k_hi(bk), q = k_lo(bk);
          md << "G " << x << " " << y << " " << p << " " << q << " " << w << "\n";
          const double g = w + prefix_log(p);
          auto it = byq.find(q);
          if (it == byq.end()) byq.emplace(q, g);
          else it->second = log_add(it->second, g);
        }
        for (const auto& [q, g] : byq) md << "Q " << x << " " << y << " " << q << " " << g << "\n";
      }
    }
    st.pr_sw_states_g = st_g;
    st.pr_sw_states_v = st_v;
    st.pr_sw_pushes = pushes;
    st.pr_mbar_max = std::max(st.pr_mbar_max, qm_max);
    st.pr_sw_sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::fprintf(stderr,
                 "[sweep-kept] l2r_mode %d allow_pk %d maps qm %lld bad %lld qm2 %lld bad %lld "
                 "pair %lld bad %lld\n",
                 l2r_mode, allow_pk ? 1 : 0, hk_map[0], hk_bad[0], hk_map[1], hk_bad[1], hk_map[2],
                 hk_bad[2]);
    static const char* ename[12] = {"E1 qm-ext", "E2 qm2-ext", "E3 interior", "E4 exterior",
                                    "E5 branch1", "E6 branch-qm-pair", "E7 branch-qm-tail",
                                    "E8 branch-qm2-pair", "E8 branch-qm2-tail", "E9 close",
                                    "E3 multi", "E4 multi-ext"};
    for (int e = 0; e < 12; ++e)
      std::fprintf(stderr, "[sweep-kept] %s push %lld bad %lld\n", ename[e], hk_e[e], hk_eb[e]);
    // E3 multi / E4 multi-ext are the nested pre-pass Multi tables (l2r_mode==0). The scan's
    // --pk-l2r on does not change np.l2r_mode; those edges are expected here.
    if (do_v && p2_v_cap <= 0 && g_ntr_K <= 0 && !Qb_b.empty()) {
      auto rel = [](double a, double b) {
        if (a < HALF_NEG && b < HALF_NEG) return 0.0;
        if (a < HALF_NEG || b < HALF_NEG) return 1e300;
        return std::fabs(std::expm1(a - b));
      };
      double worst = 0, dmin = 1e300, dmax = -1e300, wgot = 0, wref = 0;
      int wx = -1, wy = -1, wn = 0;
      long long nent = 0, nbad = 0;
      for (int j = 0; j < n; ++j) {
        for (size_t r = 0; r < Qb_at[j].slots(); ++r) {
          const int i = Qb_at[j].keys[r];
          if (i < 0 || !E.end_ok(i, j)) continue;
          const double ref = Qb_b[ix(i, j)];
          double got = NEG_INF;
          int nk = 0;
          auto it = st.prMu.find(key2(i, j));
          if (it != st.prMu.end()) {
            for (const auto& [bk, w] : it->second) {
              if (w < HALF_NEG) continue;
              ++nk;
              const int b = k_hi(bk);
              const double c = b > 0 ? C[static_cast<size_t>(b - 1)] : 0.0;
              got = log_add(got, c + w);
            }
          }
          if (got < HALF_NEG && ref < HALF_NEG) continue;
          ++nent;
          const double d = rel(got, ref);
          const double dl = (got > HALF_NEG && ref > HALF_NEG) ? got - ref : 1e300;
          if (dl < dmin) dmin = dl;
          if (dl > dmax) dmax = dl;
          if (d > 1e-9) ++nbad;
          if (d > worst) {
            worst = d;
            wx = i;
            wy = j;
            wgot = got;
            wref = ref;
            wn = nk;
          }
        }
      }
      std::fprintf(stderr,
                   "[v-out] pairs %lld bad %lld max_rel %.3e at %d.%d got %.6f ref %.6f dlog %.6f "
                   "keys %d dlog_min %.6f dlog_max %.6f\n",
                   nent, nbad, worst, wx, wy, wgot, wref, wgot - wref, wn, dmin, dmax);
    }
    std::fprintf(stderr, "[sweep-h] g_states %lld v_states %lld pushes %lld qm_max %d sec %.3f\n",
                 st_g, st_v, pushes, qm_max, st.pr_sw_sec);
    std::fprintf(stderr,
                 "[map] v_maps %lld v_mean %.4f v_max %lld g_maps %lld g_mean %.4f g_max %lld "
                 "pushes %lld\n",
                 v_kept.n, v_kept.mean(), v_kept.max(), g_kept.n, g_kept.mean(), g_kept.max(),
                 pushes);
    const double usec = pushes > 0 ? st.pr_sw_sec * 1e6 / static_cast<double>(pushes) : 0.0;
    std::fprintf(stderr,
                 "[sweep-push] g_n %lld g_mean %.3f g_p99 %lld g_pre_mean %.3f "
                 "v_n %lld v_mean %.3f v_p99 %lld v_pre_mean %.3f usec %.4f\n",
                 g_at.n, g_at.mean(), g_at.p99(), g_pre.mean(), v_at.n, v_at.mean(), v_at.p99(),
                 v_pre.mean(), usec);
    if (do_b && !Qb_b.empty()) {
      auto rel = [](double a, double b) {
        if (a < HALF_NEG && b < HALF_NEG) return 0.0;
        if (a < HALF_NEG || b < HALF_NEG) return 1e300;
        return std::fabs(std::expm1(a - b));
      };
      double worst = 0;
      int wx = -1, wy = -1;
      long long nent = 0;
      for (int j = 0; j < n; ++j) {
        for (size_t r = 0; r < Qb_at[j].slots(); ++r) {
          const int i = Qb_at[j].keys[r];
          if (i < 0 || !E.end_ok(i, j)) continue;
          const double ref = Qb_b[ix(i, j)];
          double got = NEG_INF;
          auto it = pairs.find(key2(i, j));
          if (it != pairs.end()) {
            for (const auto& [bk, w] : it->second.b) {
              if (w < HALF_NEG) continue;
              const int b = k_hi(bk);
              const double c = b > 0 ? C[static_cast<size_t>(b - 1)] : 0.0;
              got = log_add(got, c + w);
            }
          }
          if (got < HALF_NEG && ref < HALF_NEG) continue;
          ++nent;
          const double d = rel(got, ref);
          if (d > worst) {
            worst = d;
            wx = i;
            wy = j;
          }
        }
      }
      std::fprintf(stderr, "[sweep-4c] pairs %lld max_rel %.3e at %d.%d\n", nent, worst, wx, wy);
    }
  }

  void pass2_apply(StemEst& st) {
    if (p2_v_cap > 0 || p2_g_eager || p2_qz_beam > 0 || p2_sweep || p2_vsweep) st.prC = C;
    if (sw_check && !Qb_b.empty()) st.prQb_b = Qb_b.dense();
    if (p2_sweep || p2_vsweep) sw_hyper_push(st);
    if (p2_qz_beam <= 0 && p2_v_mode == 0 && p2_g_mode == 0) return;
    if (st.prC.empty()) st.prC = C;
    auto take = [&](int mode, int bp, std::vector<std::vector<int>>& byL) {
      byL.assign(n, {});
      for (int j = 0; j < n; ++j) {
        const OMap& at = Qb_at[j];
        std::vector<std::pair<double, int>> v;
        v.reserve(at.slots());
        for (size_t r = 0; r < at.slots(); ++r) {
          const int i = at.keys[r];
          if (i < 0 || at.vals[r] < HALF_NEG) continue;
          // mode 2: inside Q^b. mode 3: C_pre(i-1) * Q^b (prefix * inside).
          double key = at.vals[r];
          if (mode == 3 && i > 0) key += C[static_cast<size_t>(i - 1)];
          v.emplace_back(key, i);
        }
        if ((mode == 2 || mode == 3) && bp > 0 && static_cast<int>(v.size()) > bp) {
          std::nth_element(v.begin(), v.begin() + bp, v.end(),
                           [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                             return a.first > b.first;
                           });
          v.resize(static_cast<size_t>(bp));
        }
        for (const auto& e : v) byL[e.second].push_back(j);
      }
      for (std::vector<int>& kr : byL) {
        std::sort(kr.begin(), kr.end());
        kr.erase(std::unique(kr.begin(), kr.end()), kr.end());
      }
    };
    if (p2_v_mode > 0) take(p2_v_mode, p2_v_bp, st.prKeptRightV);
    if (p2_g_mode > 0) {
      if (p2_v_mode == p2_g_mode && p2_v_bp == p2_g_bp)
        st.prKeptRightG = st.prKeptRightV;
      else
        take(p2_g_mode, p2_g_bp, st.prKeptRightG);
    }
    if (p2_qz_beam <= 0) return;
    // --qzbar-old: keep dense prQz from build_principled; do not replace with beamed qz_col.
    if (qzbar_old) {
      st.qz_on = 0;
      st.qz_col.clear();
      st.qz_col.shrink_to_fit();
      return;
    }
    if (qz_beam_fill && p2_qz_key != 2 && !qz_keep.empty()) {
      st.qz_col.assign(static_cast<size_t>(n), {});
      st.p2_qz_n = 0;
      for (int l = 0; l < n; ++l) {
        std::vector<StemEst::QzEnt>& col = st.qz_col[static_cast<size_t>(l)];
        col.reserve(qz_keep[static_cast<size_t>(l)].size());
        for (const QkEnt& e : qz_keep[static_cast<size_t>(l)])
          col.push_back({e.k, e.q, e.p, e.u1, e.u2});
        st.p2_qz_n += static_cast<long long>(col.size());
      }
      st.qz_on = 1;
      st.prQz.clear();
      st.prQz.shrink_to_fit();
      st.p2_qz_sec = 0.0;
      qz_keep.clear();
      qz_keep.shrink_to_fit();
      return;
    }
    if (p2_qz_key == 2) {
      st.prRight.assign(n, NEG_INF);
      for (int j = 0; j < n; ++j) {
        const OMap& at = Qb_at[j];
        double acc = NEG_INF;
        for (size_t r = 0; r < at.slots(); ++r) {
          const int i = at.keys[r];
          if (i < 0 || at.vals[r] < HALF_NEG) continue;
          const double b = Qb_b[ix(i, j)];
          if (b < HALF_NEG) continue;
          acc = log_add(acc, at.vals[r] + b - logZ);
        }
        st.prRight[j] = acc;
      }
    }
    const auto t0 = std::chrono::steady_clock::now();
    struct Acc {
      double acc, u1, u2;
    };
    std::vector<int> seen(static_cast<size_t>(n), 0), touch;
    std::vector<Acc> cell(static_cast<size_t>(n));
    int gen = 0;
    st.qz_col.assign(static_cast<size_t>(n), {});
    st.p2_qz_n = 0;
    auto cell_at = [&](int k) -> Acc& {
      if (seen[static_cast<size_t>(k)] != gen) {
        seen[static_cast<size_t>(k)] = gen;
        cell[static_cast<size_t>(k)] = {NEG_INF, NEG_INF, NEG_INF};
        touch.push_back(k);
      }
      return cell[static_cast<size_t>(k)];
    };
    for (int l = 0; l < n; ++l) {
      if (++gen == 0) {
        std::fill(seen.begin(), seen.end(), 0);
        gen = 1;
      }
      touch.clear();
      const double d5 = dz && l > 0 ? E.d5w[l - 1] : 0.0;
      const double s1 = dz ? E.s1w[l] : 0.0, d3 = dz ? E.d3w[l] : 0.0;
      if (dz) cell_at(l).u1 = E.wu;
      else cell_at(l).acc = E.wu;
      if (l > 0) {
        for (const StemEst::QzEnt& e : st.qz_col[static_cast<size_t>(l - 1)]) {
          Acc& a = cell_at(e.k);
          if (dz) {
            a.u1 = e.p > HALF_NEG ? e.p + E.wu : NEG_INF;
            const double t = log_add(e.u2, e.u1 > HALF_NEG ? e.u1 + d5 : NEG_INF);
            a.u2 = t > HALF_NEG ? t + E.wu : NEG_INF;
          } else {
            a.acc = e.q + E.wu;
          }
        }
      }
      const OMap& atb = Qb_at[l];
      for (size_t r = 0; r < atb.slots(); ++r) {
        const int d = atb.keys[r];
        if (d < 0 || !E.end_ok(d, l)) continue;
        const double w = atb.vals[r] + (E.w_nest + E.at(d, l));
        if (w < HALF_NEG) continue;
        Acc& self = cell_at(d);
        self.acc = log_add(self.acc, w);
        if (d > 0) {
          for (const StemEst::QzEnt& e : st.qz_col[static_cast<size_t>(d - 1)]) {
            if (e.q < HALF_NEG) continue;
            Acc& a = cell_at(e.k);
            a.acc = log_add(a.acc, e.q + w);
          }
        }
      }
      struct Cand {
        int k;
        double q, p, u1, u2, key;
      };
      std::vector<Cand> cs;
      cs.reserve(touch.size());
      for (int k : touch) {
        const Acc& a = cell[static_cast<size_t>(k)];
        const double q = dz ? run_close(a.acc, a.u1, a.u2, s1, d3) : a.acc;
        if (q < HALF_NEG) continue;
        double key = q + (k > 0 ? st.prC[static_cast<size_t>(k - 1)] : 0.0);
        // key2 = key1 * Pr(k-1 is a pair right end). Pr = 0 ranks last and is kept
        // only when the beam still has room (a missing Q^z_pre read stays 0).
        if (p2_qz_key == 2 && k > 0) {
          const double pr = st.prRight[static_cast<size_t>(k - 1)];
          key = pr > HALF_NEG ? key + pr : NEG_INF;
        }
        cs.push_back({k, q, dz ? a.acc : 0.0, a.u1, a.u2, key});
      }
      if (g_beam_stats.f) {
        std::vector<double> keys;
        keys.reserve(cs.size());
        for (const Cand& c : cs) keys.push_back(c.key);
        std::sort(keys.begin(), keys.end(), [](double a, double b) { return a > b; });
        const size_t kept_n =
            static_cast<int>(cs.size()) > p2_qz_beam ? static_cast<size_t>(p2_qz_beam) : keys.size();
        const int full = static_cast<int>(cs.size()) > p2_qz_beam ? 1 : 0;
        g_beam_stats.log_prune("Qz_pre", l, keys.size(), kept_n, BeamStats::mass_share(keys, kept_n),
                               full);
      }
      if (static_cast<int>(cs.size()) > p2_qz_beam) {
        std::nth_element(cs.begin(), cs.begin() + p2_qz_beam, cs.end(),
                         [](const Cand& a, const Cand& b) { return a.key > b.key; });
        cs.resize(static_cast<size_t>(p2_qz_beam));
      }
      std::sort(cs.begin(), cs.end(), [](const Cand& a, const Cand& b) { return a.k < b.k; });
      std::vector<StemEst::QzEnt>& col = st.qz_col[static_cast<size_t>(l)];
      col.reserve(cs.size());
      for (const Cand& c : cs) col.push_back({c.k, c.q, c.p, c.u1, c.u2});
      st.p2_qz_n += static_cast<long long>(col.size());
    }
    st.qz_on = 1;
    st.prQz.clear();
    st.prQz.shrink_to_fit();
    st.p2_qz_sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }

  // build_principled's Bml cell by cell from the three edge lists (long double).
  void bml_check(StemEst& st) const {
    const long double w2 = E.w_alpha2, w3 = E.w_alpha3;
    for (int d = 0; d < n; ++d)
      for (int e = d; e < n; ++e) {
        std::vector<long double> t;
        for (int i = std::max(0, d - std::min(E.ml_left, n)); i <= d; ++i) {
          if (!lead_gap_ok(d - i)) continue;
          if (!Qm_at[e].has(i) || Qm_b[ix(i, e)] < HALF_NEG) continue;
          long double w = Qm_b[ix(i, e)] + E.ml_pk_k * w2 + w3 * (d - i);
          if (dz) w += E.run(i, d - 1);
          t.push_back(w);
        }
        if (d >= 1)
          for (int i = 0; i <= d - 1; ++i) {
            if (!Qm_at[d - 1].has(i) || Qm_b[ix(i, e)] < HALF_NEG) continue;
            t.push_back(static_cast<long double>(Qm_b[ix(i, e)]) + Qm[ix(i, d - 1)] +
                        E.ml_pk_k * w2);
          }
        for (int i = 0; i < d; ++i)
          for (int j = e + 1; j < n; ++j) {
            if (!Qb_at[j].has(i) || !E.end_ok(i, j) || Qb_b[ix(i, j)] < HALF_NEG) continue;
            long double base = static_cast<long double>(Qb_b[ix(i, j)]) + E.w_alpha1 +
                               (1.0L + E.ml_pk_k) * w2 + E.at_close(i, j) + w3 * (j - 1 - e);
            if (dz) base += E.run(e + 1, j - 1);
            long double l1 = w3 * (d - i - 1);
            if (dz) l1 += E.run(i + 1, d - 1);
            t.push_back(base + l1);
            if (i + 1 <= d - 1 && Qm[ix(i + 1, d - 1)] > HALF_NEG)
              t.push_back(base + Qm[ix(i + 1, d - 1)]);
          }
        const double got = st.Bml(d, e);
        ++st.bml_check_cells;
        if (t.empty()) {
          if (got > HALF_NEG) ++st.bml_check_mis;
          continue;
        }
        long double m = *std::max_element(t.begin(), t.end()), s = 0.0L;
        for (long double x : t) s += std::exp(x - m);
        const double ref = static_cast<double>(m + std::log(s));
        if (got < HALF_NEG) ++st.bml_check_mis;
        else st.bml_check_rel = std::max(st.bml_check_rel, std::fabs(std::expm1(got - ref)));
      }
  }

  // st's --stem-est-hr av data from this nested-only run (runs the outside pass):
  // phi(i, j) = (pu(i) + P(i, j)) (pu(j) + P(i, j)), pu(k) = max(0, 1 - sum_j P(k, j)),
  // raised to gamma; Tav(x, y) = phi [E_out(x, y) + w_stack(x-1, y+1; x, y) Tav(x-1, y+1)]
  // over the pairs that can form, outermost first; entries Uav(x', y') = E_in Tav for the
  // end_ok pairs (E_out = E_in = at if end_ok, else 0).
  void build_hr_av(StemEst& st, double gamma) {
    nested_lphi(st, gamma);
    const size_t nn = static_cast<size_t>(n) * n;
    std::vector<double> T(nn, NEG_INF);
    std::vector<stemidx::Entry> ents;
    for (int x = 0; x < n; ++x)
      for (int y : partners_right[x]) {
        const double lp = st.lphi[ix(x, y)];
        if (lp < HALF_NEG) continue;
        double acc = E.end_ok(x, y) ? E.at(x, y) : NEG_INF;
        if (x >= 1 && y + 1 <= n - 1) {
          const double t0 = T[ix(x - 1, y + 1)];
          if (t0 > HALF_NEG) {
            const double w = E.stem_interior(x - 1, x, y, y + 1);
            if (w > HALF_NEG) acc = log_add(acc, w + t0);
          }
        }
        if (acc < HALF_NEG) continue;
        T[ix(x, y)] = lp + acc;
        if (E.end_ok(x, y)) ents.push_back({x, y, E.at(x, y) + T[ix(x, y)]});
      }
    st.dropped_av = stem_sparse(ents);
    st.entries_av = ents.size();
    st.index_av = stemidx::PersistentSum(n, ents);
    st.hr_av = true;
  }

  // st's --stem-est-hr phi data (after build_stem_est on the stem-table run, whose S, C0 and
  // loop st keeps): the original stem table with phi_gamma(x, y) multiplied in once per pair,
  //   T'(x, y) = phi_gamma(x, y) [E_out(x, y) S[y+1] + sum w(x0, y0 | x, y) T'(x0, y0)],
  // entries E_in(x', y') T'(x', y') C0[y'-1], same loop order and index as build_stem_est.
  // gamma = 0 adds nothing, so the entries equal the original ones bit for bit.
  void build_hr_phi(StemEst& st, double gamma) {
    nested_lphi(st, gamma);
    const bool one = gamma == 0.0;
    if (one)
      for (int x = 0; x < n; ++x)
        for (int y : partners_right[x]) st.lphi[ix(x, y)] = 0.0;
    const std::vector<double>& S = st.S;
    const std::vector<double>& C0 = st.C0;
    const int lc = std::min(st.loop, cap);
    std::vector<double> T(static_cast<size_t>(n) * n, NEG_INF);
    std::vector<stemidx::Entry> ents;
    for (int x = 0; x < n; ++x)
      for (int y : partners_right[x]) {
        double acc = E.end_ok(x, y) ? E.at(x, y) + S[y + 1] : NEG_INF;
        for (int x0 = std::max(0, x - 1 - lc); x0 < x; ++x0)
          for (int y0 = y + 1; y0 <= std::min(n - 1, y + 1 + lc); ++y0) {
            double t0 = T[ix(x0, y0)];
            if (t0 < HALF_NEG) continue;
            double w = E.stem_interior(x0, x, y, y0);
            if (w < HALF_NEG) continue;
            acc = log_add(acc, w + t0);
          }
        const double t = one ? acc : st.lphi[ix(x, y)] + acc;
        T[ix(x, y)] = t;
        if (t > HALF_NEG && E.end_ok(x, y)) ents.push_back({x, y, (E.at(x, y) + t) + C0[y - 1]});
      }
    st.dropped_av = stem_sparse(ents);
    st.entries_av = ents.size();
    st.index_av = stemidx::PersistentSum(n, ents);
    st.hr_av = true;
    st.hr_phi = true;
  }

  // st.lphi = gamma (log(pu(x) + P(x, y)) + log(pu(y) + P(x, y))) over the pairs that can form
  // (NEG_INF where a factor is 0), from this nested-only run's outside pass.
  void nested_lphi(StemEst& st, double gamma) {
    outside();
    const size_t nn = static_cast<size_t>(n) * n;
    std::vector<double> P(nn, 0.0), sp(n, 0.0);
    for (int j = 0; j < n; ++j) {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        const int i = at.keys[r];
        if (i < 0) continue;
        const double b = Qb_b[ix(i, j)];
        if (at.vals[r] < HALF_NEG || b < HALF_NEG) continue;
        const double p = std::exp(at.vals[r] + b - logZ);
        P[ix(i, j)] = p;
        sp[i] += p;
        sp[j] += p;
      }
    }
    st.lphi.assign(nn, NEG_INF);
    for (int x = 0; x < n; ++x)
      for (int y : partners_right[x]) {
        const double p = P[ix(x, y)];
        const double a = std::max(0.0, 1.0 - sp[x]) + p, b = std::max(0.0, 1.0 - sp[y]) + p;
        if (a > 0.0 && b > 0.0) st.lphi[ix(x, y)] = gamma * (std::log(a) + std::log(b));
      }
  }

  // --stem-est rank (role both: Qgl / Qgr also scored as stem 2).
  // qg2 (role all): Qg also scored as stem 2 (needs st->mirror).  qg_partner
  // (--stem-est-qg-r partner): Qg's stem-2 score from the partner tables (stem_Bg) instead of
  // the mirrored index.  consistent (--stem-est-consts consistent): every role adds beta1, the
  // closing at(p, q) and only the stem terms not yet in alpha (StemEval).
  void set_loop_cap(int lc) { cap = lc; }
  void set_stem(const StemEst* st, bool both, bool qg2 = false, bool qg_partner = false,
                bool consistent = false) {
    stem_ = st;
    stem_both = both;
    stem_qg2 = qg2;
    stem_qg_partner = qg_partner;
    stem_consistent = consistent;
    first_f = stemidx::MinTree(n + 1);
  }
  // --stem-est-qg-ext (mode 1 stack, 2 bulge1; needs the partner Qg score), K_max, split.
  void set_stem_ext(int mode, int kmax, bool split, bool check) {
    stem_ext = mode;
    stem_ext_max = kmax;
    stem_split = split;
    ext_check = check;
  }
  // --stem-est-qg-ext loop parameters (Lg, cut c, node cap N) and --stem-est-helpers grow
  // helpers: 0 off, 1 grow, 2 nested
  void set_stem_grow(int lg, double cut, int nodes, int helpers, int ml_span) {
    grow_lg = lg;
    grow_cut = cut;
    grow_nodes = nodes;
    stem_helpers = helpers;
    stem_ml_span = ml_span;
  }
  // --stem-est-profile on, --stem-est-profile-eps, --stem-est-ml-tail nested,
  // --stem-profile-check
  void set_stem_profile(bool on, double eps, bool tail, bool check) {
    stem_profile = on;
    prof_eps = eps;
    ml_tail = tail;
    prof_check = check;
  }
  // --stem-est-chain on (needs set_stem_profile on and st->chain), -span, -cut, -ml, -loop
  void set_stem_chain(bool on, int span, double cut, bool ml, int loop) {
    stem_chain = on;
    chain_span = span;
    chain_cut = cut;
    chain_ml = ml;
    chain_loop = loop;
  }
  // --stem-est-hr av / phi (needs st->hr_av): every L-stem score uses Hav(u, v, q) S[q+1]
  // (av) or Hav(u, v, q) / C0[q] (phi) in place of H(u, v, q) / C0[q]; feasibility
  // (StemEval::ok) keeps H.  check: --stem-hr-check.
  void set_stem_hr(bool av, bool sq, bool check) {
    hr_av = av;
    hr_sq = sq;
    hr_check = check;
  }
  // --stem-est-principled on (needs st->pr): every pseudoknot table ranked by
  // alpha + log Out_hat; growth / future stems use the interior loop cap loop and multiloop
  // steps (ml) of span <= ml_span; Qgls closings run <= ml_close; profiles cut at cut x max;
  // arms_k arms per gap (0 = all); check: --stem-pr-check.
  void set_stem_principled(bool on, int loop, bool ml, int ml_span, int ml_close, double cut,
                           int arms_k, bool check) {
    pr_on = on;
    pr_loop = loop;
    pr_ml = ml;
    pr_ml_span = ml_span;
    pr_ml_close = ml_close;
    pr_cut = cut;
    pr_arms_k = arms_k;
    pr_check = check;
  }
  // --pr-arms-rank ratio: arms by phi / Z_arm (Z_arm = wu per arm base x'-xs+1), added in
  // decreasing order until the rest of the gap's total < arms_eps x the selected sum, at most
  // arms_max (replaces --pr-arms-k).  --pr-cut-by gapfree: profile entries cut by their
  // gap-free contribution (growth: x C[p'-1] Sbar[q'+1]; future stems: x Sbar[ys+1] /
  // e^(wu (x - xs))); cut_roots: the same cut on the root and key profiles.
  void set_stem_pr_select(bool ratio, double arms_eps, int arms_max, bool gapfree,
                          bool cut_roots) {
    pr_ratio = ratio;
    pr_arms_eps = arms_eps;
    pr_arms_max = arms_max;
    pr_gapfree = gapfree;
    pr_cut_roots = cut_roots;
  }
  void set_stem_pr_sep(int mode, int arm_z, bool qz_unp, bool check) {
    pr_arms_mode = mode;
    pr_arm_z = arm_z;
    pr_qz_unp = qz_unp;
    pr_sep_check = check;
  }
  // --pr-grow-select bound: growth node, root and key profiles keep their endpoints by
  // iota = w U(p', q') (largest first, until 1 - eps of the total, at most B; B = 0: no cap;
  // eps = 0 and B = 0: everything); check: --pr-grow-check; dump: --pr-score-dump.
  // --pr-armq-sep on: phi(arm) at q' with Qzbar[q'+1][y'-1] replaced by the upper bound
  // Qzbar[0][y'-1] / Qzbar[0][q'], from per-arm suffix sums built once.
  void set_stem_pr_armq_sep(bool on, size_t aq_max) {
    pr_aq_sep = on;
    pr_aq_max = aq_max;
  }
  // --pr-arm-persist on: the top-k arm rankings of every q' from one persistent argmax tree per
  // order (DP03::pidx_build).  --pr-arm-gaponly on: Step A over the arms whose x' lies in the
  // union of the gaps of the states being scored (DP03::pr_gap_begin).  --pr-partner-span W:
  // only partner entries with y' - q' <= W (0: all).  --pr-v-eps e: V maps and arm maps keep
  // their entries by V(xs, ys) Sbar[ys+1] / Qzbar[xs][x] until 1 - e of the total (default
  // 1e-8; 0: all).
  void set_stem_pr_arms(bool persist, bool gaponly, int pspan, double veps) {
    pr_arm_persist = persist;
    pr_arm_gaponly = gaponly;
    pr_pspan = pspan;
    pr_veps = veps;
  }
  // --pr-v-ratio on: arm entries from Vr; --pr-v-diag PREFIX: PREFIX.ventries.tsv.gz and, with
  // --pr-v-sample P, PREFIX.sel.tsv (each L-lookup sampled with probability P, at most 10000).
  void set_stem_pr_v(bool ratio, const std::string& diag, double sample_p) {
    pr_vratio = ratio;
    pr_vdiag = diag;
    pr_v2dump = g_v2dump;
    pr_vsample_p = diag.empty() ? 0.0 : sample_p;
    if (pr_vsample_p > 0.0) {
      pr_vsel = std::fopen((diag + ".sel.tsv").c_str(), "w");
      std::fprintf(pr_vsel, "qp\tu\tv\ttwo\tnarmsV\tnarmsVr\toverlap\tvr_in_v\tshare_vr\tshare_vtop\tshare_eng\t"
                            "overlap_r2\tr2_in_v\tshare_r2\toverlap_sep\tsep_in_v\tshare_sep\n");
      pr_varms = popen(("gzip -1 > '" + diag + ".arms.tsv.gz'").c_str(), "w");
      std::fprintf(pr_varms, "sample\txs\txp\tlogV\tlogVr\tlogVr2\tlogVsep\n");
    }
  }
  void pr_v_finish(const std::string& sid) {
    pr_vdump(sid);
    if (pr_vsel) std::fclose(pr_vsel);
    pr_vsel = nullptr;
    if (pr_varms) pclose(pr_varms);
    pr_varms = nullptr;
  }
  void set_stem_pr_grow(bool bound, double eps, int B, bool check, FILE* dump) {
    pr_gsel = bound;
    pr_geps = eps;
    pr_gB = B;
    pr_gcheck = check;
    pr_dump = dump;
  }
  // sep: distinct queries, columns of the range-sum tables built, (step, q') groups queried,
  // their finite arms, steps with a query, table time; check: queries / arm terms compared,
  // max relative error of S against brute force, max |log sep - log exact|, mismatches
  // (finite against infinite), histograms of the log ratio separable / exact per arm term
  // (weighted by its share of the query's exact sum) and per query, on [-10, 10] in 0.25 bins
  // plus two tails
  void pr_sep_report(FILE* f, const std::string& sid) const {
    std::string ht, hq;
    char b[64];
    double wt = 0.0;
    for (int k = 0; k < PSC_NB; ++k) {
      std::snprintf(b, sizeof b, "%s%.6g", k ? "," : "", psc_ht[k]);
      ht += b;
      std::snprintf(b, sizeof b, "%s%.6g", k ? "," : "", psc_hq[k]);
      hq += b;
      wt += psc_ht[k];
    }
    std::fprintf(f,
                 "[pr-sep-check] %s queries %lld terms %lld tab_rel %.3e tab_mis %lld "
                 "sep_max %.6e sep_mis %lld q_mean %.6e q_abs %.6e term_wmean %.6e "
                 "term_wabs %.6e\n[pr-sep-hist] %s term %s\n[pr-sep-hist] %s query %s\n",
                 sid.c_str(), psc_q, psc_terms, psc_tab_rel, psc_tab_mis, psc_sep_max, psc_sep_mis,
                 psc_q ? psc_d_sum / psc_q : 0.0, psc_q ? psc_d_abs / psc_q : 0.0,
                 wt ? psc_lr_w / wt : 0.0, wt ? psc_lr_wabs / wt : 0.0, sid.c_str(), ht.c_str(),
                 sid.c_str(), hq.c_str());
  }
  long long pr_sep_steps_all() const { return pr_sep_steps + (prSepTouched.empty() ? 0 : 1); }
  static constexpr int PSC_NB = 82;
  long long pr_sep_q = 0, pr_sep_cols = 0, pr_sep_groups = 0, pr_sep_garms = 0,
            pr_sep_steps = 0, psc_q = 0, psc_terms = 0, psc_tab_mis = 0, psc_sep_mis = 0;
  double pr_sep_sec = 0.0, psc_tab_rel = 0.0, psc_sep_max = 0.0, psc_lr_w = 0.0,
         psc_lr_wabs = 0.0, psc_d_sum = 0.0, psc_d_abs = 0.0;
  double psc_ht[PSC_NB] = {}, psc_hq[PSC_NB] = {};
  static int psc_bin(double x) {
    if (x < -10.0) return 0;
    if (x >= 10.0) return PSC_NB - 1;
    return 1 + std::min(PSC_NB - 3, static_cast<int>((x + 10.0) / 0.25));
  }
  // --stem-est-principled-nested on (needs st->pr_nested): Qb / Qm / Qm2 ranked by
  // alpha + their nested outside.
  void set_stem_pr_nested(bool on) { pr_nested = on; }
  // principled: per table (Qg, Qgl, Qgr, Qgls, Qgrs, Qp) ranked states and scoring time;
  // key-profile entries evaluated; T builds (cache hits), arms selected / skipped (outside
  // the gap) over them, ys entries of the T's; K builds (hits) and their entries; SL+SR
  // misses (hits); per-q' arm structures built; static part: V entries, arms, arm entries,
  // G profiles and entries, entries cut; check: states compared per table, finite / infinite
  // mismatches, max relative error
  long long pr_states[6] = {}, pr_evals = 0, pr_T_built = 0, pr_T_hits = 0, pr_T_arms = 0,
            pr_T_skip = 0, pr_T_ys = 0, pr_K_built = 0, pr_K_hits = 0, pr_K_is = 0,
            pr_slr_miss = 0, pr_slr_hits = 0, pr_aq_built = 0, pr_V_entries = 0, pr_arms = 0,
            pr_arm_entries = 0, pr_G_nodes = 0, pr_G_entries = 0, pr_cut_drop = 0,
            pr_roots = 0, pr_keys = 0, pr_key_entries = 0;
  double pr_sec[6] = {}, pr_static_sec = 0.0;
  // grow select: profiles selected (0 growth nodes, 1 evaluated roots / keys), their entries
  // before / after, summed kept iota shares, profiles with zero total; U_L columns added,
  // partner entries added, U evaluations; growth and selection time (all modes: growth).
  // check: endpoints, bound violations (L, R), max log(L / (Qz U_L)), max log(R / (Qz U_R)),
  // endpoints whose selected-arm L exceeds the all-arm L
  long long gs_prof[2] = {}, gs_before[2] = {}, gs_after[2] = {}, gs_zero[2] = {}, gs_cols = 0,
            gs_rent = 0, gs_U = 0, gc_n = 0, gc_vL = 0, gc_vR = 0, gc_sel_gt = 0;
  double gs_share[2] = {}, pr_grow_sec = 0.0, pr_sel_sec = 0.0, gc_maxL = NEG_INF,
         gc_maxR = NEG_INF;
  // g-read profile: one-time fill vs per-root sw_as_rn aggregation.
  double sw_ntr_init_sec = 0.0, sw_fill_sec = 0.0, sw_as_rn_sec = 0.0;
  long long sw_as_rn_calls = 0, sw_as_rn_ent_in = 0, sw_as_rn_ent_out = 0;
  long long sw_ga_nodes = 0, sw_ga_ent_sum = 0, sw_ga_ent_max = 0;
  // V/G fill split (pass-2 V sweep vs scan sw_fill).
  double v_merge_sec = 0.0, v_cut_sec = 0.0, g_merge_sec = 0.0, g_cut_sec = 0.0;
  long long v_pushes = 0, g_pushes = 0;
  long long v_pre_n = 0, v_pre_sum = 0, v_pre_max = 0, v_post_n = 0, v_post_sum = 0, v_post_max = 0;
  long long g_pre_n = 0, g_pre_sum = 0, g_pre_max = 0, g_post_n = 0, g_post_sum = 0, g_post_max = 0;
  // sw_fill multiloop linearity: right-gap yo visits attributed to each pair at column y.
  long long ml_yo_qm_n = 0, ml_yo_qm_sum = 0, ml_yo_qm_max = 0;
  long long ml_yo_sc_n = 0, ml_yo_sc_sum = 0, ml_yo_sc_max = 0;
  int ml_yo_S = 0;
  bool ml_yo_qm_ix = false;
  // K_L inside growth reads (lsc_L) and per-arm (b,c) V/K stats.
  double kl_in_read_sec = 0.0;
  long long arm_read_n = 0;
  long long arm_Vf_sum = 0, arm_Vf_max = 0;
  std::vector<int> arm_Vf_samples;  // for median
  long long arm_Kq_sum = 0, arm_Kq_max = 0;
  std::vector<char> arm_seen;
  // steps of the built nodes by kind: [0] growth nodes + roots, [1] future stems;
  // [.][0] interior, [.][1] multiloop.  Partner terms: wall time of the outermost T / K /
  // Stop build.  stop_j: the step at which a CPU-limit stop was taken (-1: none).
  long long st_kind[2][2] = {}, st_nodes[2] = {};
  double pr_partner_sec = 0.0;
  bool pr_aq_sep = false;
  double pr_aq_sec = 0.0;
  std::vector<double> prSF;
  // inclusive seconds of the pr_T / pr_K / pr_SLR / pr_W misses; W misses and entries read
  double pr_sec_T = 0.0, pr_sec_K = 0.0, pr_sec_SLR = 0.0, pr_sec_W = 0.0;
  long long pr_W_built = 0, pr_W_ent = 0;
  size_t pr_aq_max = 128;  // --pr-aq-max
  bool pr_arm_persist = false, pr_arm_gaponly = false;
  int pr_pspan = 0;
  double pr_veps = 0.0;
  // --pr-v-trim suffix: drop xs.ys only against the entries with xs' >= xs; --pr-v-cap: entries
  // per key kept (largest share of their suffix first), 0 = no cap
  bool pr_vsuffix = false;
  int pr_vcap = 0;
  // Pass 2 limits. All default off: 0 leaves the current tables and ranking keys.
  int p2_qz_beam = 0, p2_qz_key = 1, p2_v_mode = 0, p2_v_bp = 0, p2_v_cap = 0, p2_v_trim = 0;
  int p2_g_mode = 0, p2_g_bp = 0, p2_g_cap = 0;
  bool sw_check = false, p2_sweep = false, p2_vsweep = false, geq_pre = false;
  int p2_v_floor = 0;  // 0 off, 1 smallest kept mu, 2 log-weight 0
  bool p2_v_cut_fb = false;  // --v-cut-fb: kind-(i) miss = kappa_min / C_pre(b-1)
  bool lazy_ml_cap = false;  // --lazy-ml-cap: lazy leading run uses ml_lead_cap
  bool p2_g_suffix = false, p2_g_eager = false, p2_g_fb = false;
  // --pre-key inside: beam rank is the inside weight. Default is C(i-1) + inside.
  bool pre_key_inside = false;
  // --qz-beam on the nested pre-pass: fill Q^z_pre beamed instead of the dense table.
  // The scan's own Q^z (part of Z) is never beamed.
  bool qz_beam_fill = false, qzbar_old = false, prof_on = false;
  long long qz_cand = 0, qz_kept = 0;
  // Pass-1 rule iterations (per nested pfunc). Printed as [p1rule] per nucleotide.
  long long r_hair = 0, r_int = 0, r_m2xp = 0, r_m2multi = 0, r_muext = 0, r_mmu = 0, r_mp = 0,
            r_extc = 0;
  double sec_qb = 0, sec_qm = 0, sec_branch = 0, sec_qz = 0, sec_c = 0;
  double sec_hair = 0, sec_int = 0, sec_m2xp = 0, sec_emit = 0, sec_mclose = 0, sec_mp = 0;
  double prof_qz_sec = 0.0;
  int prof_region = 0;             // 0 other, 1 R prune, 2 spawn rank, 3 H / P3
  double prof_grow_reg[4] = {};
  double sx_k_sec = 0, sx_rkey_sec = 0;
  long long sx_k_n = 0, sx_m_sum = 0;
  // 0 other, 1 R root (rec_root), 2 R key (rec_key), 3 spawn rank
  int sx_ksrc = 0;
  long long sx_k_src[4] = {};
  long long sx_rkey_n = 0;
  long long sx_m_h[8] = {};
  // query-weighted length of the node-trunc list that the K call walks (no layer m on that path)
  long long sx_list_h[65] = {};
  // --k-memo: lK(arm, q') = arm_phi, a function of the pass-2 arm tables only
  bool k_memo = false;
  std::unordered_map<uint64_t, double> lKm;
  long long lKm_calls = 0, lKm_hits = 0;
  // --k-bound: compare lK against the exact gap sum once, after the arms exist.
  // --k-exact-ub: R beams rank by the exact gap sum, in upper-bound order, and stop early.
  bool k_bound = false, k_exact_ub = false, k_exact_check = false, sx_k_exact = false;
  std::unordered_map<uint64_t, double> lKxM;
  long long lKx_calls = 0, lKx_hits = 0, sx_kx_sigma = 0, k_exact_cols = 0, k_exact_mis = 0;
  double sx_kx_sec = 0.0;
  long long prof_help_steps = 0;
  struct QkEnt {
    int k;
    double q, p, u1, u2;
  };
  std::vector<std::vector<QkEnt>> qz_keep;
  long long v_capped = 0;
  // cap: pairs capped; per capped pair the dropped share of the trim weights (sum, max) and of
  // V itself (max); log sums of V over the capped pairs' entries and over the dropped ones
  long long v_cap_pairs = 0;
  double v_cap_fsum = 0.0, v_cap_fmax = 0.0, v_cap_fwmax = 0.0, v_cap_all = NEG_INF,
         v_cap_drop = NEG_INF;
  // --pr-v-mlk K (0 = all): a V node (scan = false, pr_fast 2) keeps its K multiloop helper
  // terms of largest mass w + log sum(helper) and scales them by (kept + dropped) / kept, so
  // the node's mass is unchanged
  int pr_vmlk = 0;
  // --pr-v-mlsteps on: the V nodes enumerate their multiloop steps under pr_fast 2 as well (only
  // the step counts use them; the helpers carry those terms)
  bool pr_vmlsteps = false;
  long long vk_nodes = 0, vk_cand = 0, vk_drop = 0;
  double vk_rmax = 0.0, vk_rsum = 0.0, vk_mass_s = 0.0;
  // --pr-v-key xs: V keyed by xs alone (key2(xs, 0)), the value the sum over ys of
  // V e^beta1 Sbar[ys+1]; arm entries one per (xs, x', y') (prYw), no per-ys entries
  // (prEy / prEys / prEw empty, prArmE0 = prArmY0).
  bool pr_vxs = false;
  // --pr-v-stream on: pass 0 and pass 1 in one y-descending sweep; V(x', y') freed once
  // y' - y - 1 exceeds every step's reach; arm entries bucketed by xs, each bucket freed as
  // it is appended to the arm tables.  v_live_max: most V entries alive at once.
  bool pr_vstream = false;
  long long v_live_max = 0, v_pairs = 0;
  // persistent index: build seconds, nodes, leaf updates; gaponly: Step A arms indexed, union
  // epochs; V trimming: entries before / after
  double pidx_sec = 0.0, pidx_mb = 0.0;
  long long pidx_nodes = 0, pidx_upd = 0, aq_arms = 0, gap_epochs = 0, gap_miss = 0, v_before = 0,
            v_after = 0;
  struct SecTimer {
    double& acc;
    std::chrono::steady_clock::time_point t0;
    explicit SecTimer(double& a) : acc(a), t0(std::chrono::steady_clock::now()) {}
    ~SecTimer() { acc += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); }
  };
  // stop_phase: where the stop was taken (0 scan, 1 the one-time V build); cur_j: current step.
  int pr_pdepth = 0, stop_j = -1, stop_phase = 0, cur_j = -1;
  void cpu_stop_check(int phase) {
    if (!g_cpu_stop) return;
    stop_j = cur_j;
    stop_phase = phase;
    throw CpuStop{};
  }
  struct PrPartnerTimer {
    DP03& d;
    std::chrono::steady_clock::time_point t0;
    explicit PrPartnerTimer(DP03& x) : d(x), t0(std::chrono::steady_clock::now()) { ++d.pr_pdepth; }
    ~PrPartnerTimer() {
      if (--d.pr_pdepth == 0)
        d.pr_partner_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }
  };
  long long prc_states[6] = {}, prc_mis[6] = {};
  double prc_rel[6] = {};
  // av / phi: ranked Qg states with both L-stem factors finite, the sum / min / max of
  // (log Hav + hr_c(q)) - (log H - log C0[q]) over them, states finite under the new mode only
  // / orig only; check: Qg / Qgl states compared with hr_naive, finite / infinite mismatches,
  // max relative error
  long long hr_n = 0, hr_av_only = 0, hr_orig_only = 0, hrc_states = 0, hrc_mis = 0;
  double hr_dsum = 0.0, hr_dmin = 0.0, hr_dmax = 0.0, hrc_rel = 0.0;
  // chain: chain nodes built, edges by step (b interior, f, g, h multiloop), root profiles,
  // their L / R entries, entries dropped by the cut; chain time (inside pair_prof)
  long long ch_nodes = 0, ch_edges[4] = {}, ch_roots = 0, ch_ly = 0, ch_rx = 0, ch_cut_drop = 0;
  double ch_sec = 0.0;
  // profiles: steps with a profile, keys built (Qg pairs, Qgls (i, seed), Qgrs c), pair
  // profiles built, L / R entries of the key profiles, entries dropped by eps, entries
  // evaluated by the ranked states, tail closings; Bg lookups, raw computations and the
  // partner entries those visit; build / total scoring time
  long long pf_steps = 0, pf_keys[3] = {}, pf_pairs = 0, pf_ly = 0, pf_rx = 0, pf_trunc = 0,
            pf_evals = 0, pf_tail_cl = 0, pf_states = 0;
  long long bg_lookups = 0, bg_misses = 0, bg_visits = 0;
  double pf_build_sec = 0.0, pf_total_sec = 0.0;
  // --stem-est-bg-reorder on: M rows built (key, v), the (x', e) entries they visit, the e
  // entries visited by the scored states; check: states compared, finite / infinite
  // mismatches, max relative error, ranking events / states whose keep decision differs
  bool bg_reorder = false, ro_check = false;
  long long ro_rows = 0, ro_build = 0, ro_eval = 0, roc_states = 0, roc_mis = 0, roc_events = 0,
            roc_flip_events = 0, roc_flip_states = 0;
  double roc_rel = 0.0, ro_off = NEG_INF;
  void set_stem_bg_reorder(bool on, bool check) {
    bg_reorder = on;
    ro_check = check;
  }
  // --stem-profile-check: Qg / helper states compared, finite / infinite mismatches, max
  // relative error against the Result-19 score
  long long pc_states[2] = {}, pc_mis[2] = {};
  double pc_rel[2] = {};
  long long stem_removed[3] = {};  // Qg, Qgl, Qgr
  double stem_qg_sec = 0.0;  // time in the Qg stem scoring (feasibility + ranking)
  // sigma_ext: ranked Qg states, extension nodes over them, states with a finite sigma_ext,
  // time; --stem-ext-check: states compared, finite / infinite mismatches, max relative error
  long long ext_states = 0, ext_nodes = 0, ext_finite = 0;
  double ext_sec = 0.0;
  long long xc_states = 0, xc_mis = 0;
  double xc_rel = 0.0;
  // growth node lists built (one per outer pair and step), their nodes, layer nodes reached
  // and dropped by the cut / the node cap; helpers: ranked Qgls / Qgrs states, closings and
  // growth nodes over them, finite scores, time; --stem-ext-check on the helper scores
  long long grow_lists = 0, grow_list_nodes = 0, grow_layer = 0, grow_cut_drop = 0,
            grow_cap_drop = 0;
  long long hlp_states[2] = {}, hlp_closings = 0, hlp_nodes = 0, hlp_finite = 0;
  double hlp_sec = 0.0;
  long long hc_states = 0, hc_mis = 0;
  double hc_rel = 0.0;

  // --stem-bg-check (after the outside pass, every table at infinite beam): stem_Bg of every
  // stored Qg state against a brute-force sum over the stored Qgl states, and the posterior
  // weight of the Qg states a d-dependent removal (H_R = 0 and Bg(d, u, v) = 0) would drop.
  std::string stem_bg_check(const std::string& id, double logZ) {
    struct L { int i, e, f; double val; };
    std::unordered_map<int, std::vector<L>> by_d;
    for (int f = 0; f < n; ++f)
      for (const Grp& g : Qgl[f].grps)
        for (const Ent& s : g.ents)
          if (E.end_ok(s.p, f))
            by_d[s.u + 1].push_back({s.p, g.v, f, prefix_log(s.p) + s.a + E.at(s.p, f)});
    stBc.clear();
    stBgc.clear();
    stc_q = -1;
    long long states = 0, zero_mis = 0, rule_n = 0, ff_only = 0, ff_viol = 0;
    double max_rel = 0.0, rule_w = 0.0;
    std::vector<long double> terms;
    for (int q = 0; q < n; ++q)
      for (const Grp& g : Qg[q].grps)
        for (const Ent& s : g.ents) {
          const int d = s.p, u = s.u, v = g.v;
          ++states;
          terms.clear();
          auto it = by_d.find(d);
          if (E.end_ok(u, v) && it != by_d.end())
            for (const L& x : it->second) {
              if (x.e < u + 1 || x.e > v - 1 || x.f > v - 1) continue;
              const double t = x.val + qz_log(x.f + 1, v - 1) + qz_log(u + 1, x.e - 1);
              if (t > HALF_NEG) terms.push_back(t);
            }
          double brute = NEG_INF;
          if (!terms.empty()) {
            long double m = *std::max_element(terms.begin(), terms.end()), sum = 0.0L;
            for (long double t : terms) sum += std::exp(t - m);
            brute = static_cast<double>(m + std::log(sum) +
                                        (static_cast<long double>(E.w_stem) + E.at(u, v)));
          }
          const double bg = stem_Bg_raw(d, u, v);
          if ((bg > HALF_NEG) != (brute > HALF_NEG)) ++zero_mis;
          else if (bg > HALF_NEG) max_rel = std::max(max_rel, std::fabs(std::expm1(bg - brute)));
          const bool h0 = stem_->H(u, v, q) < HALF_NEG;
          const bool ff = first_f.query(u + 1, v) <= v - 1;
          if (bg > HALF_NEG && !ff) ++ff_viol;
          if (h0 && bg < HALF_NEG) {
            ++rule_n;
            if (s.b > HALF_NEG) rule_w += std::exp(s.a + s.b - logZ);
            if (ff) ++ff_only;
          }
        }
    char buf[512];
    std::snprintf(buf, sizeof buf,
                  "[bg-check] %s n %d states %lld max_rel %.3e zero_mismatch %lld ff_violations "
                  "%lld rule_drop %lld rule_drop_ff_ok %lld rule_drop_w %.6e",
                  id.c_str(), n, states, max_rel, zero_mis, ff_viol, rule_n, ff_only, rule_w);
    return buf;
  }

  // --stem-diag: one JSON line of stem-estimate diagnostics (definition at stem_diag).
  std::string stem_diag(const std::vector<const StemEst*>& ests, const std::string& id);
  void tagged_outside(int tag);

 private:
  Energy E;
  int n;
  int beam;
  bool allow_pk;
  int cap;
  bool score_qz;
  bool arm_score;
  const AccessRec* access;
  int max_gap_scan;
  int max_arm_len;
  bool gap_step;
  GateCfg gcfg;
  const Gate* gate_;
  // --stem-est: partner table A(d, e, f) of the surviving stem-1 Qgl states, per key2(d, e)
  // with f ascending; first_f[e] = first f of any entry (., e, f)
  struct StemA {
    std::vector<int> fs;
    std::vector<double> vals;
    std::vector<std::vector<std::pair<int, double>>> il;  // pr_on: (p, alpha + at(p, f)) per f
  };
  const StemEst* stem_ = nullptr;
  bool stem_both = false, stem_qg2 = false, stem_qg_partner = false, stem_consistent = false;
  std::unordered_map<uint64_t, StemA> stA;
  stemidx::MinTree first_f;
  // e values with an entry (d, e, .) per d, ascending; B(d, e, v) and Bg(d, u, v) caches of
  // the step stc_q (all their terms end before the step, so they are final)
  std::unordered_map<int, std::vector<int>> stE;
  std::unordered_map<uint64_t, double> stBc, stBgc;
  int stc_q = -1;
  // sigma_ext: growth nodes per outer pair key2(p, q) and H(u, v, q') per key3, of step
  // stc_q; helper outside O(x, y; u, v) per key4 (mode 3 = loop)
  struct ExtNode {
    int d, q;
    double w;
  };
  int stem_ext = 0, stem_ext_max = 30;
  bool stem_split = false, ext_check = false;
  int grow_lg = 3, grow_nodes = 2000, stem_ml_span = 30;
  double grow_cut = 1e-8;
  int stem_helpers = 0;
  // nested helper sums for step hn_q: Qgls by 2 * i + (v == q), Qgrs by c; count < 0 = unset
  int hn_q = -1;
  std::vector<std::pair<double, int>> hnL, hnR;
  std::unordered_map<uint64_t, std::vector<ExtNode>> stX;
  std::unordered_map<uint64_t, double> stHc, stOc;
  // --stem-est-profile: sparse log profiles of a key over the final right ends y' (L part,
  // coefficient of Qz[u+1][v-1] H_R(u, v, y')) and left ends x' (R part, coefficient of
  // Bg(x', u, v)), index ascending; tail: the Qgls ml tail (coefficient of Qz[u+1][v-1]);
  // ncl: closings of a helper key.  Per step: untruncated pair profiles by key2(x, y), and
  // the scored key profiles of Qg (by p), Qgls (by 2 i + seed) and Qgrs (by c).
  struct Prof {
    std::vector<int> ly, rx;
    std::vector<double> lw, rw;
    double tail = NEG_INF;
    int ncl = 0;
  };
  bool stem_profile = false, ml_tail = false, prof_check = false;
  double prof_eps = 0.0;
  int pf_q = -1;
  std::unordered_map<uint64_t, Prof> stP;
  std::unordered_map<int, Prof> stPg, stPls, stPrs;
  // --stem-est-bg-reorder: per key profile and v, M(e, v) = sum_x' P_R(x') B(x', e, v) over
  // the e >= lo with a finite term, ascending (log); dense accumulator of the row being built
  struct MRow {
    int lo = 0;
    std::vector<int> es;
    std::vector<double> ms;
  };
  std::unordered_map<const Prof*, std::unordered_map<int, MRow>> stM;
  std::vector<double> mAcc;
  std::vector<int> mStamp, mTouched;
  int m_stamp = 0;
  // prof_eval's dense H / Bg rows of the current gap (u, v) (stamped).
  std::vector<double> rowH, rowB;
  std::vector<int> rowHs, rowBs;
  int row_u = -1, row_v = -1, row_q = -1, row_stamp = 0;
  int tbeam[BT_N];
  // --dump-states with dangles: the run states (P, U1, U2) of Qm, Qm2, Qmp, PKM per
  // cell (i, j) (PKM: (d, j)), and their posterior log weight from the outside pass
  bool keep_runs;
  enum { RS_QM, RS_QM2, RS_QMP, RS_PKM, RS_N };
  std::vector<double> runP[RS_N], runU1[RS_N], runU2[RS_N], runW[RS_N];
  std::unordered_map<uint64_t, std::vector<ArmCand>> arm_cache;
  int arm_cache_q = -1;
  bool inside_done = false;
  bool outside_done = false;
  double logZ = NEG_INF;

  AtView Qb, Qm, Qm2, Qp;
  SpMat Qz;
  std::vector<double> C;
  std::vector<OMap> Qb_at, Qm_at, Qm2_at, Qp_at;
  // Qmp: Qm with the stem-multiloop penalties b_p, c_p; qmp_at() / qmp_b() fall back to
  // Qm unless E.stem_ml
  AtView Qmp;
  SpMat Qmp_b;
  std::vector<OMap> Qmp_at;
  const OMap& qmp_at(int j) const { return E.stem_ml ? Qmp_at[j] : Qm_at[j]; }
  SpMat& qmp_b() { return E.stem_ml ? Qmp_b : Qm_b; }
  // Per column, Qm slot indices bucketed by the base at key-1 (the left end of a case-(h) pair).
  mutable std::vector<std::array<std::vector<int>, 4>> qm_bkt;
  mutable std::vector<char> qm_bkt_ok;
  void qm_reindex(int j) {
    if (static_cast<int>(qm_bkt.size()) != n) {
      qm_bkt.assign(n, {});
      qm_bkt_ok.assign(n, 0);
    }
    auto& b = qm_bkt[j];
    for (auto& v : b) v.clear();
    const OMap& at = qmp_at(j);
    for (size_t r = 0; r < at.slots(); ++r) {
      const int i = at.keys[r];
      if (i < 1) continue;
      const int base = E.code[i - 1];
      if (base < 0 || base > 3) continue;
      b[base].push_back(static_cast<int>(r));
    }
    qm_bkt_ok[j] = 1;
  }
  // Slots of Qm at `col` whose base at key-1 can pair with the base at `right`.
  // Falls back to every slot when the column has not been bucketed.
  template <class F>
  void qm_pair_slots(int col, int right, F&& fn) const {
    const OMap& at = qmp_at(col);
    if (static_cast<int>(qm_bkt_ok.size()) == n && qm_bkt_ok[col]) {
      const int bj = E.code[right];
      for (int bi = 0; bi < 4; ++bi) {
        const int t = pair_type(bi, bj);
        if (t < 0 || (!allow_gu && (t == 3 || t == 5))) continue;
        for (int r : qm_bkt[col][bi]) fn(r);
      }
      return;
    }
    for (size_t r = 0; r < at.slots(); ++r) fn(static_cast<int>(r));
  }
  std::vector<OMap> PKM_at;  // PKM_at[j][d] = sum_{e<=j} Qp[d][e] e^{-α3 (j-e)}, beamed
  std::vector<Tab4> Qg, Qgls, Qgrs, Qgl, Qgr;
  std::vector<std::vector<int>> partners_left, partners_right;
  std::vector<double> qb_acc;  // fill_Qb scratch, one row, reset on the partners just written
  std::vector<std::vector<LeadEnt>> qlead;
  // Pass 1 (l2r_mode == 0) enumerates like LinearPartition, via mate_next only.
  // Interior loops use the energy cap. The multiloop leading run uses ml_lead_cap.
  // Multi key is (p, q, kind): 0 empty right run, 1 one base, 2 length >= 2.
  // A length >= 2 extension does not depend on where the run started, so those
  // arrivals merge. Empty and one-base runs do: the dangle is s1w of that base.
  std::vector<int> mate_after[4];
  bool mates_ready = false;
  std::vector<OMap> lpH;   // hairpin, column = right end, key = left end
  std::vector<OMap> lpM;   // Multi inside, key = p * 3 + kind
  std::vector<OMap> lpMb;  // Multi outside
  std::vector<OMap> lpP;   // M2 weight that ends on a branch, no trailing run

  SpMat Qb_b, Qm_b, Qm2_b, Qp_b, PKM_b;
  SpMat Qz_b;
  std::vector<double> C_b;

  // Dangles: trailing-run states (P, U1, U2) of the tables that grow by unpaired
  // bases, indexed by column parity; the stored table value closes the run for a
  // paired j+1: X = P + U1 s1w[j] + U2 d3w[j].  Outside keeps the U1, U2 outsides.
  bool dz = false;
  RunSt rsQm[2], rsQm2[2], rsQmp[2], rsPKM[2];
  std::vector<double> zP[2], zU1[2], zU2[2];
  double cP = 0.0, cU1 = NEG_INF, cU2 = NEG_INF;
  RunSt rbQm[2], rbQm2[2], rbQmp[2], rbPKM[2];
  std::vector<double> bz1[2], bz2[2];
  double bc1 = NEG_INF, bc2 = NEG_INF;

  static std::vector<int> live_keys(const OMap& m) {
    std::vector<int> out;
    for (size_t r = 0; r < m.slots(); ++r)
      if (m.keys[r] >= 0) out.push_back(m.keys[r]);
    return out;
  }

  // U1, U2 at j from the run states pv at j-1 of the survivors surv.
  void run_ext(const RunSt& pv, RunSt& cu, const std::vector<int>& surv, int j, double w) {
    const double d5 = E.d5w[j - 1];
    for (int i : surv) {
      double p = pv.P.get(i);
      if (p > HALF_NEG) cu.U1.put(i, p + w);
      double t = log_add(pv.U2.get(i), pv.U1.get(i) + d5);
      if (t > HALF_NEG) cu.U2.put(i, t + w);
    }
  }

  static double run_close(double p, double u1, double u2, double s1, double d3) {
    double v = p;
    if (u1 > HALF_NEG) v = log_add(v, u1 + s1);
    if (u2 > HALF_NEG) v = log_add(v, u2 + d3);
    return v;
  }

  // Closed values at j into at (and tbl unless null); keys first, then new keys of P.
  void run_fin(AtView* tbl, OMap& at, const RunSt& c, const std::vector<int>& first,
               int j) {
    const double s1 = E.s1w[j], d3 = E.d3w[j];
    auto one = [&](int i) {
      if (at.has(i)) return;
      double v = run_close(c.P.get(i), c.U1.get(i), c.U2.get(i), s1, d3);
      if (v > HALF_NEG) {
        if (tbl) (*tbl)[ix(i, j)] = v;
        at.set(i, v);
      }
    };
    for (int i : first) one(i);
    for (int i : c.P.keys) one(i);
  }

  void run_back(double bF, double bU1n, double bU2n, int j, double w, double& bP, double& bU1,
                double& bU2) const {
    bP = log_add(bF, bU1n + w);
    bU1 = log_add(bF + E.s1w[j], bU2n + (E.d5w[j] + w));
    bU2 = log_add(bF + E.d3w[j], bU2n + w);
  }

  // Replaces the outside of the closed values at j by that of P; keeps U1's, U2's.
  // rs >= 0 (keep_runs): also the cell's posterior log weight into runW[rs].
  template <class Beta>
  void run_back_at(Beta& beta, const OMap& at, int j, double w, RunSt* rb, int rs = -1) {
    const RunSt& nx = rb[(j + 1) & 1];
    RunSt& cu = rb[j & 1];
    cu.clear();
    for (size_t r = 0; r < at.slots(); ++r) {
      int i = at.keys[r];
      if (i < 0) continue;
      double bP, b1, b2;
      run_back(static_cast<double>(beta[ix(i, j)]), nx.U1.get(i), nx.U2.get(i), j, w, bP, b1, b2);
      beta[ix(i, j)] = bP;
      cu.U1.put(i, b1);
      cu.U2.put(i, b2);
      if (rs >= 0 && keep_runs) run_weight(rs, i, j, bP, b1, b2);
    }
  }

  // Every derivation through a run cell is in exactly one of its states P, U1, U2 at j.
  void run_weight(int rs, int i, int j, double bP, double b1, double b2) {
    const size_t k = ix(i, j);
    double v = NEG_INF;
    if (runP[rs][k] > HALF_NEG && bP > HALF_NEG) v = log_add(v, runP[rs][k] + bP);
    if (runU1[rs][k] > HALF_NEG && b1 > HALF_NEG) v = log_add(v, runU1[rs][k] + b1);
    if (runU2[rs][k] > HALF_NEG && b2 > HALF_NEG) v = log_add(v, runU2[rs][k] + b2);
    runW[rs][k] = v;
  }

  void save_runs(int j) {
    auto save = [&](int rs, const RunSt& c, const OMap& at) {
      for (size_t r = 0; r < at.slots(); ++r) {
        int i = at.keys[r];
        if (i < 0) continue;
        runP[rs][ix(i, j)] = c.P.get(i);
        runU1[rs][ix(i, j)] = c.U1.get(i);
        runU2[rs][ix(i, j)] = c.U2.get(i);
      }
    };
    save(RS_QM, rsQm[j & 1], Qm_at[j]);
    save(RS_QM2, rsQm2[j & 1], Qm2_at[j]);
    if (allow_pk && E.stem_ml) save(RS_QMP, rsQmp[j & 1], Qmp_at[j]);
    if (allow_pk) save(RS_PKM, rsPKM[j & 1], PKM_at[j]);
  }

  size_t ix(int i, int j) const { return static_cast<size_t>(i) * n + j; }

  void alloc() {
    size_t nn = static_cast<size_t>(n) * n;
    Qz.init(n);
    C.assign(n, NEG_INF);
    auto init_maps = [&](std::vector<OMap>& v) {
      v.resize(n);
      for (auto& m : v) m.init(n);
    };
    init_maps(Qb_at);
    init_maps(Qm_at);
    init_maps(Qm2_at);
    init_maps(Qp_at);
    init_maps(PKM_at);
    Qb.bind(&Qb_at, n);
    Qm.bind(&Qm_at, n);
    Qm2.bind(&Qm2_at, n);
    Qp.bind(&Qp_at, n);
    if (E.stem_ml) {
      init_maps(Qmp_at);
      Qmp.bind(&Qmp_at, n);
    }
    dz = E.dangle;
    if (dz) {
      for (int k = 0; k < 2; ++k) {
        rsQm[k].init(n);
        rsQm2[k].init(n);
        rsQmp[k].init(n);
        rsPKM[k].init(n);
        zP[k].assign(n, NEG_INF);
        zU1[k].assign(n, NEG_INF);
        zU2[k].assign(n, NEG_INF);
      }
      if (keep_runs)
        for (int r = 0; r < RS_N; ++r) {
          runP[r].assign(nn, NEG_INF);
          runU1[r].assign(nn, NEG_INF);
          runU2[r].assign(nn, NEG_INF);
          runW[r].assign(nn, NEG_INF);
        }
    }
    Qg.resize(n);
    Qgls.resize(n);
    Qgrs.resize(n);
    Qgl.resize(n);
    Qgr.resize(n);
  }

  void build_pairs() {
    partners_left.assign(n, {});
    partners_right.assign(n, {});
    for (int j = 0; j < n; ++j)
      for (int i = 0; i <= j - MIN_PAIR; ++i)
        if (E.pairs(i, j)) {
          partners_left[j].push_back(i);
          partners_right[i].push_back(j);
        }
  }
  // Pass 1 does not build this list. The scan and the stem index do, on first use.
  void ensure_pairs() const {
    if (n <= 0 || !partners_left.empty()) return;
    const_cast<DP03*>(this)->build_pairs();
  }

  // ---------------------------------------------------------------- scores
  double qz_log(int left, int right) const {
    if (left > right) return 0.0;
    if (left < 0 || right >= n) return NEG_INF;
    ++qz_reads;
    if (qz_audit && right > l2r_done[LT_QZ]) ++qz_unf;
    return Qz[ix(left, right)];
  }

  double prefix_log(int p) const {
    if (p <= 0) return 0.0;
    return C[p - 1];
  }

  double four_score(int p, int u, int v, double log_x) const {
    double s = prefix_log(p) + log_x;
    if (score_qz) s += qz_log(u + 1, v - 1);
    return s;
  }

  // Arm windows of gap [u+1, v-1] at step q with D_q > 0 as (log F-part,
  // log D_q - log acc, acc from floor), Pareto-pruned (see DP03._arm_front in engine.py).
  const std::vector<ArmCand>& arm_front(int u, int v, int q) {
    if (q != arm_cache_q) {
      arm_cache.clear();
      arm_cache_q = q;
    }
    uint64_t key = key2(u, v);
    auto it = arm_cache.find(key);
    if (it != arm_cache.end()) return it->second;
    const int lo = u + 1, hi = v - 1, M = max_gap_scan, L = max_arm_len;
    std::vector<ArmCand> cand;
    auto scan = [&](int s) {
      const bool near_s = s - lo < M;
      for (const AccWin& w : access->by_start[s]) {
        const int t = w.t;
        if (t > hi || t - s + 1 > L) break;
        if (!near_s && hi - t >= M) continue;
        if (w.ys.back() <= q) continue;
        size_t k = std::upper_bound(w.ys.begin(), w.ys.end(), q) - w.ys.begin();
        const bool floor = !access->has_acc(s, t);
        ++arm_stats.windows;
        arm_stats.floor_windows += floor;
        double z = w.lsuf[k] - access->log_acc(s, t);
        double a = qz_log(lo, s - 1) + (t - s + 1) * E.wu + qz_log(t + 1, hi);
        cand.push_back({a, z, floor});
      }
    };
    if (lo <= hi) {
      if (hi - lo + 1 > M) {
        for (int s = lo; s < lo + M; ++s) scan(s);
        for (int s = std::max(lo + M, hi - M + 2 - L); s <= hi; ++s) scan(s);
      } else {
        for (int s = lo; s <= hi; ++s) scan(s);
      }
    }
    std::stable_sort(cand.begin(), cand.end(), [](const ArmCand& x, const ArmCand& y) {
      if (x.a != y.a) return x.a > y.a;
      return x.z > y.z;
    });
    std::vector<ArmCand> front;
    double best_z = NEG_INF;
    for (const auto& c : cand)
      if (c.z > best_z) {
        front.push_back(c);
        best_z = c.z;
      }
    return arm_cache.emplace(key, std::move(front)).first->second;
  }

  static double softplus(double z) {
    if (z > 0.0) return z + std::log1p(std::exp(-z));
    return std::log1p(std::exp(z));
  }

  double arm_score_of(int p, int u, int v, int q, double log_x) {
    const double base = qz_log(u + 1, v - 1);
    double best = base, chosen = -1.0;
    bool chosen_floor = false;
    const auto& front = arm_front(u, v, q);
    if (!front.empty()) {
      const double log_omp = std::log(std::max(1e-3, 1.0 - access->pair_prob(p, q)));
      for (const ArmCand& c : front) {
        double lk = softplus(c.z - log_omp);
#ifdef DP03_DEBUG
        if (!(lk >= 0.0)) { std::fprintf(stderr, "kappa < 1\n"); std::abort(); }
#endif
        double val = c.a + lk;
        if (val > best) {
          best = val;
          chosen = lk;
          chosen_floor = c.floor;
        }
      }
    }
#ifdef DP03_DEBUG
    if (!(best >= base)) { std::fprintf(stderr, "G < F(gap)\n"); std::abort(); }
#endif
    ++arm_stats.scored;
    if (chosen >= 0.0) {
      ++arm_stats.arm_wins;
      arm_stats.sum_log_kappa += chosen;
      if (chosen > arm_stats.max_log_kappa) arm_stats.max_log_kappa = chosen;
      arm_stats.chosen_floor += chosen_floor;
      int k = 0;
      while (k < LK_HIST_BINS - 1 && chosen >= LK_HIST_EDGES[k]) ++k;
      ++arm_stats.lk_hist[k];
    }
    return prefix_log(p) + log_x + best;
  }

  // gtab: GATE_QG / GATE_QGL for the tables the gate may act on (stats are kept for them
  // either way), -1 for the others.
  void prune4(Tab4& t, int bt, int q, bool arm = false, int gtab = -1) {
    if (stem_ && stem_->lean) StemEst::lean_abort("prune4");
    const bool use_arm = arm && arm_score;
    const int beam = tbeam[bt];
    if (remove_keys) {
      std::vector<char> keep;
      bool any = false;
      for (const Grp& g : t.grps)
        for (const Ent& e : g.ents) {
          keep.push_back(!remove_keys->count(state_key(bt, e.p, e.u, g.v, q)));
          any |= !keep.back();
        }
      if (any) filter_step(t, keep);
    }
    if (stem_ && (bt == BT_QG || bt == BT_QGL || bt == BT_QGR)) {
      prune4_stem(t, bt, q, gtab, beam);
      return;
    }
    if (stem_ && (stem_helpers || pr_on) && (bt == BT_QGLS || bt == BT_QGRS)) {
      prune4_helper(t, bt, q, beam);
      return;
    }
    if (gtab >= 0 && gcfg.on && gcfg.tables[gtab] && gate_) {
      prune4_gated(t, q, use_arm, gtab, beam);
      return;
    }
    if (beam <= 0) {
      if (gtab >= 0)
        for (const Grp& g : t.grps) {
          gate_stats.cand[gtab] += g.ents.size();
          gate_stats.elig[gtab] += g.ents.size();
          gate_stats.kept[gtab] += g.ents.size();
        }
      return;
    }
    if (gap_step) {
      prune4_step(t, bt, q, use_arm, gtab, beam);
      return;
    }
    std::vector<double> scores;
    std::vector<char> keep;
    for (auto& g : t.grps) {
      if (gtab >= 0) {
        gate_stats.cand[gtab] += g.ents.size();
        gate_stats.elig[gtab] += g.ents.size();
      }
      if (static_cast<int>(g.ents.size()) <= beam) {
        if (gtab >= 0) gate_stats.kept[gtab] += g.ents.size();
        continue;
      }
      scores.clear();
      for (const Ent& e : g.ents)
        scores.push_back(use_arm ? arm_score_of(e.p, e.u, g.v, q, e.a)
                                 : four_score(e.p, e.u, g.v, e.a));
      keep_mask(scores, keep, beam);
      if (drops) {
        std::vector<size_t> rank = score_ranks(scores);
        ++drops->events[bt];
        for (size_t k = 0; k < g.ents.size(); ++k)
          if (!keep[k])
            note_drop(bt, bt, q, g.v, g.ents[k].p, g.ents[k].u, g.v, q, rank[k], scores.size(),
                      scores[k]);
      }
      filter_group(g, keep, 0);
      if (gtab >= 0) {
        gate_stats.ranked[gtab] += scores.size();
        gate_stats.kept[gtab] += g.ents.size();
      }
    }
  }

  // ------------------------------------------------------- stem estimate
  // log B(d, e, e2) = log sum over f <= e2 - 1 of A(d, e, f) Qz[f+1][e2-1].
  double stem_B(int d, int e, int e2) const {
    auto it = stA.find(key2(d, e));
    double acc = NEG_INF;
    if (it == stA.end()) return acc;
    const StemA& a = it->second;
    for (size_t k = 0; k < a.fs.size(); ++k) {
      if (a.fs[k] > e2 - 1) break;
      acc = log_add(acc, a.vals[k] + qz_log(a.fs[k] + 1, e2 - 1));
    }
    return acc;
  }

  static uint64_t key3(int a, int b, int c) {
    return (static_cast<uint64_t>(a) << 42) | (static_cast<uint64_t>(b) << 21) |
           static_cast<uint64_t>(c);
  }

  void stem_cache_step(int q) {
    if (q == stc_q) return;
    stc_q = q;
    std::unique_ptr<SecTimer> ptm(lpf ? new SecTimer(lpf_stc_sec) : nullptr);
    if (lpf) {
      ++lpf_stc_n;
      lpf_stc_buckets += static_cast<long long>(stBc.bucket_count() + stBgc.bucket_count() + stX.bucket_count() +
                                                stHc.bucket_count() + stOc.bucket_count());
      if (pr_on)
        lpf_stc_scan += static_cast<long long>(prTc.size() + prAQ.size() + prKls.bucket_count() +
                                               prKrs.bucket_count() + prRKls.bucket_count() +
                                               prRKrs.bucket_count() + prRT.bucket_count() + prAc.bucket_count());
    }
    stBc.clear();
    stBgc.clear();
    stX.clear();
    stHc.clear();
    stOc.clear();
    if (stem_profile) {
      stP.clear();
      stPg.clear();
      stPls.clear();
      stPrs.clear();
      stM.clear();
    }
    if (pr_on) pr_step(q);
  }

  double stem_B_cached(int d, int e, int v) {
    const uint64_t key = key3(d, e, v);
    auto it = stBc.find(key);
    if (it != stBc.end()) return it->second;
    double b = stem_B(d, e, v);
    stBc.emplace(key, b);
    return b;
  }

  // Partner weight of a Qg state (d, u, v) as the right stem (DP03._stem_Bg): the engine's
  // Qg -> Qgl edge (beta2 + at(u, v), end_ok(u, v), filler Qz[u+1][e-1]) times
  // B(d, e, v) (Qgl -> Qgr filler, Qp with the surviving left stems), over the e of stE[d]
  // in [u+1, v-1], ascending.
  double stem_Bg_raw(int d, int u, int v) {
    if (!E.end_ok(u, v)) return NEG_INF;
    auto it = stE.find(d);
    if (it == stE.end()) return NEG_INF;
    const std::vector<int>& es = it->second;
    double acc = NEG_INF;
    for (auto k = std::lower_bound(es.begin(), es.end(), u + 1); k != es.end() && *k <= v - 1; ++k) {
      const int e = *k;
      ++bg_visits;
      double b = stem_B_cached(d, e, v);
      if (b < HALF_NEG) continue;
      double qz = qz_log(u + 1, e - 1);
      if (qz < HALF_NEG) continue;
      acc = log_add(acc, qz + b);
    }
    if (acc < HALF_NEG) return NEG_INF;
    return acc + (E.w_stem + E.at(u, v));
  }

  double stem_Bg(int d, int u, int v) {
    ++bg_lookups;
    const uint64_t key = key3(d, u, v);
    auto it = stBgc.find(key);
    if (it != stBgc.end()) return it->second;
    ++bg_misses;
    double bg = stem_Bg_raw(d, u, v);
    stBgc.emplace(key, bg);
    return bg;
  }

  // L-stem factor of a final right end q: log H(u, v, q) - log C0[q]; with hr_av the index
  // is st->index_av (av or phi entries), and with hr_sq (av) the constant is log S[q+1]
  // (h + hr_c(q) equals h - C0[q] bit for bit otherwise).
  bool hr_av = false, hr_sq = false, hr_check = false;
  double hr_h(int u, int v, int q) const { return hr_av ? stem_->Hav(u, v, q) : stem_->H(u, v, q); }
  double hr_c(int q) const { return hr_sq ? stem_->S[q + 1] : -stem_->C0[q]; }

  // --stem-hr-check: Hav(u, v, q) by listing every partner stem (innermost x'.y', x' in
  // [u+1, v-1], y' > q, stacked outward while the pairs can form) and summing
  // E_out E_in prod w_stack prod phi_gamma in long double.
  double hr_naive(int u, int v, int q) const {
    if (stem_->hr_phi) return hr_naive_phi(u, v, q);
    std::vector<long double> terms;
    for (int x = u + 1; x <= v - 1; ++x)
      for (int y = q + 1; y < n; ++y) {
        if (y - x < MIN_PAIR || !E.pairs(x, y) || !E.end_ok(x, y)) continue;
        const double lp0 = stem_->lphi[ix(x, y)];
        if (lp0 < HALF_NEG) continue;
        long double w = static_cast<long double>(E.at(x, y)) + lp0;
        for (int i = x, j = y;;) {
          if (E.end_ok(i, j)) terms.push_back(w + E.at(i, j));
          const int i0 = i - 1, j0 = j + 1;
          if (i0 < 0 || j0 > n - 1 || j0 - i0 < MIN_PAIR || !E.pairs(i0, j0)) break;
          const double s = E.stem_interior(i0, i, j, j0), lp = stem_->lphi[ix(i0, j0)];
          if (s < HALF_NEG || lp < HALF_NEG) break;
          w += static_cast<long double>(s) + lp;
          i = i0;
          j = j0;
        }
      }
    if (terms.empty()) return NEG_INF;
    long double m = *std::max_element(terms.begin(), terms.end()), sum = 0.0L;
    for (long double t : terms) sum += std::exp(t - m);
    return static_cast<double>(m + std::log(sum));
  }

  // --stem-hr-check with phi: H_R'(u, v, q) by listing every partner stem of the original
  // table (innermost x'.y' an end pair, x' in [u+1, v-1], y' > q; outward steps with at most
  // the table's loop unpaired bases per side through pairs that can form; outermost pair an
  // end pair) and summing E_in C0[y'-1] prod phi_gamma prod w E_out S[j_1+1] in long double.
  double hr_naive_phi(int u, int v, int q) const {
    const int lc = std::min(stem_->loop, cap);
    auto partner = [&](int i, int j) {
      if (i < 0 || j > n - 1) return false;
      const auto& pr = partners_right[i];
      return std::find(pr.begin(), pr.end(), j) != pr.end();
    };
    std::vector<long double> terms;
    std::vector<std::tuple<int, int, long double>> stack;
    for (int x = u + 1; x <= v - 1; ++x)
      for (int y : partners_right[x]) {
        if (y <= q || !E.end_ok(x, y)) continue;
        const double lp = stem_->lphi[ix(x, y)];
        if (lp < HALF_NEG) continue;
        stack.emplace_back(x, y, static_cast<long double>(E.at(x, y)) + stem_->C0[y - 1] + lp);
        while (!stack.empty()) {
          auto [i, j, w] = stack.back();
          stack.pop_back();
          if (E.end_ok(i, j)) terms.push_back(w + E.at(i, j) + stem_->S[j + 1]);
          for (int i0 = std::max(0, i - 1 - lc); i0 < i; ++i0)
            for (int j0 = j + 1; j0 <= std::min(n - 1, j + 1 + lc); ++j0) {
              if (!partner(i0, j0)) continue;
              const double s = E.stem_interior(i0, i, j, j0), lp0 = stem_->lphi[ix(i0, j0)];
              if (s < HALF_NEG || lp0 < HALF_NEG) continue;
              stack.emplace_back(i0, j0, w + s + lp0);
            }
        }
      }
    if (terms.empty()) return NEG_INF;
    long double m = *std::max_element(terms.begin(), terms.end()), sum = 0.0L;
    for (long double t : terms) sum += std::exp(t - m);
    return static_cast<double>(m + std::log(sum));
  }

  double stem_H_cached(int u, int v, int q2) {
    const uint64_t key = key3(u, v, q2);
    auto it = stHc.find(key);
    if (it != stHc.end()) return it->second;
    double h = hr_h(u, v, q2);
    stHc.emplace(key, h);
    return h;
  }

  // Extension nodes of a Qg state with outer pair p.stc_q (DP03._stem_ext_nodes).
  const std::vector<ExtNode>& stem_ext_nodes(int p) { return stem_grow_nodes(p, stc_q); }

  // Growth nodes of outer pair p.q (DP03._stem_grow_nodes): the outer pairs (d', q') reached
  // by 1..stem_ext_max case-(b) wraps (weight stem_interior + pair_bonus each; stack: d'-1,
  // q'+1; bulge1 also d'-2, q'+1 and d'-1, q'+2; loop: d'-1-a, q'+1+b for a, b <=
  // min(grow_lg, cap), a outer), with log W summed over the paths (layer by layer; nodes in
  // first-reached order).  loop: after each layer, nodes with log W < log(grow_cut) + the
  // layer max are dropped, then new nodes beyond grow_nodes in total are dropped (lowest W
  // first, ties by order); dropped nodes are neither kept nor expanded.  Only end_ok nodes
  // are kept: every exit of a Qg state needs end_ok of its outer pair.
  const std::vector<ExtNode>& stem_grow_nodes(int p, int q) {
    const uint64_t key = key2(p, q);
    auto it = stX.find(key);
    if (it != stX.end()) return it->second;
    std::vector<ExtNode> out;
    if (stem_ext == 3) {
      const int lg = std::min(grow_lg, cap);
      const double lcut = grow_cut > 0.0 ? std::log(grow_cut) : NEG_INF;
      std::vector<ExtNode> cur{{p, q, 0.0}}, nxt, tot;
      std::unordered_map<uint64_t, size_t> nidx, tidx;
      std::vector<size_t> fresh, order;
      std::vector<char> keep;
      for (int t = 1; t <= stem_ext_max && !cur.empty(); ++t) {
        nxt.clear();
        nidx.clear();
        for (const ExtNode& c : cur)
          for (int a = 0; a <= lg; ++a)
            for (int b = 0; b <= lg; ++b) {
              const int i = c.d - 1 - a, j = c.q + 1 + b;
              if (i < 0 || j > n - 1) continue;
              const double s = E.stem_interior(i, c.d, c.q, j);
              if (s < HALF_NEG) continue;
              const double val = c.w + (s + pair_bonus);
              auto r = nidx.emplace(key2(i, j), nxt.size());
              if (r.second) nxt.push_back({i, j, val});
              else nxt[r.first->second].w = log_add(nxt[r.first->second].w, val);
            }
        grow_layer += static_cast<long long>(nxt.size());
        keep.assign(nxt.size(), 1);
        if (grow_cut > 0.0 && !nxt.empty()) {
          double mx = NEG_INF;
          for (const ExtNode& x : nxt) mx = std::max(mx, x.w);
          const double thr = lcut + mx;
          for (size_t k = 0; k < nxt.size(); ++k)
            if (nxt[k].w < thr) {
              keep[k] = 0;
              ++grow_cut_drop;
            }
        }
        if (grow_nodes > 0) {
          fresh.clear();
          for (size_t k = 0; k < nxt.size(); ++k)
            if (keep[k] && !tidx.count(key2(nxt[k].d, nxt[k].q))) fresh.push_back(k);
          const size_t room = static_cast<size_t>(grow_nodes) > tot.size()
                                  ? static_cast<size_t>(grow_nodes) - tot.size() : 0;
          if (fresh.size() > room) {
            order.resize(fresh.size());
            for (size_t m = 0; m < fresh.size(); ++m) order[m] = m;
            std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
              return nxt[fresh[a]].w > nxt[fresh[b]].w;
            });
            for (size_t m = room; m < order.size(); ++m) keep[fresh[order[m]]] = 0;
            grow_cap_drop += static_cast<long long>(fresh.size() - room);
          }
        }
        cur.clear();
        for (size_t k = 0; k < nxt.size(); ++k) {
          if (!keep[k]) continue;
          const ExtNode& x = nxt[k];
          cur.push_back(x);
          auto r = tidx.emplace(key2(x.d, x.q), tot.size());
          if (r.second) tot.push_back(x);
          else tot[r.first->second].w = log_add(tot[r.first->second].w, x.w);
        }
      }
      grow_list_nodes += static_cast<long long>(tot.size());
      for (const ExtNode& x : tot)
        if (E.end_ok(x.d, x.q)) out.push_back(x);
    } else if (stem_ext == 1) {
      double w = 0.0;
      for (int t = 1; t <= stem_ext_max; ++t) {
        const int i = p - t, j = q + t;
        if (i < 0 || j > n - 1) break;
        const double s = E.stem_interior(i, i + 1, j - 1, j);
        if (s < HALF_NEG) break;
        w = w + (s + pair_bonus);
        ++grow_list_nodes;
        if (E.end_ok(i, j)) out.push_back({i, j, w});
      }
    } else {
      static const int MV[3][2] = {{1, 1}, {2, 1}, {1, 2}};
      std::vector<ExtNode> cur{{p, q, 0.0}}, nxt, tot;
      std::unordered_map<uint64_t, size_t> nidx, tidx;
      for (int t = 1; t <= stem_ext_max && !cur.empty(); ++t) {
        nxt.clear();
        nidx.clear();
        for (const ExtNode& c : cur)
          for (const auto& m : MV) {
            const int i = c.d - m[0], j = c.q + m[1];
            if (i < 0 || j > n - 1) continue;
            const double s = E.stem_interior(i, c.d, c.q, j);
            if (s < HALF_NEG) continue;
            const double val = c.w + (s + pair_bonus);
            auto r = nidx.emplace(key2(i, j), nxt.size());
            if (r.second) nxt.push_back({i, j, val});
            else nxt[r.first->second].w = log_add(nxt[r.first->second].w, val);
          }
        for (const ExtNode& x : nxt) {
          auto r = tidx.emplace(key2(x.d, x.q), tot.size());
          if (r.second) tot.push_back(x);
          else tot[r.first->second].w = log_add(tot[r.first->second].w, x.w);
        }
        cur.swap(nxt);
      }
      grow_list_nodes += static_cast<long long>(tot.size());
      for (const ExtNode& x : tot)
        if (E.end_ok(x.d, x.q)) out.push_back(x);
    }
    ++grow_lists;
    return stX.emplace(key, std::move(out)).first->second;
  }

  // sigma_ext of the Qg state (p, u, v, stc_q) with log weight a (DP03._stem_ext_score): log
  // sum over its extension nodes (d', q', W) of the stem-1 score four_score(d', u, v, a) + W +
  // H(u, v, q') - C0[q'] + w_pk and the stem-2 score Bg(d', u, v) + a + W + S[q'+1].
  double stem_ext_score(int p, int u, int v, double a) {
    const std::vector<ExtNode>& xs = stem_ext_nodes(p);
    const double w_pk = E.w_beta1 + 2.0 * E.w_stem;
    double acc = NEG_INF;
    for (const ExtNode& x : xs) {
      double sl = NEG_INF, sr = NEG_INF;
      const double h = stem_H_cached(u, v, x.q);
      if (h > HALF_NEG) {
        sl = four_score(x.d, u, v, a) + x.w + h + hr_c(x.q) + w_pk;
        if (stem_consistent) sl += (E.at(u, v) + E.at(x.d, x.q));
      }
      const double bg = stem_Bg(x.d, u, v);
      if (bg > HALF_NEG) {
        sr = bg + a + x.w + stem_->S[x.q + 1];
        if (stem_consistent) sr += (E.at(x.d, x.q) + E.w_beta1);
      }
      acc = log_add(acc, log_add(sl, sr));
    }
    ext_nodes += static_cast<long long>(xs.size());
    return acc;
  }

  // --stem-ext-check: sigma_ext by walking every wrap path separately (no shared path sums;
  // H and Bg recomputed), every term summed in long double.
  void ext_paths(int d, int q, int t, long double w, std::vector<std::tuple<int, int, long double>>& out) {
    if (t > 0 && E.end_ok(d, q)) out.emplace_back(d, q, w);
    if (t == stem_ext_max) return;
    std::vector<std::pair<int, int>> mv;
    if (stem_ext == 3) {
      const int lg = std::min(grow_lg, cap);
      for (int a = 0; a <= lg; ++a)
        for (int b = 0; b <= lg; ++b) mv.emplace_back(a + 1, b + 1);
    } else {
      mv = {{1, 1}, {2, 1}, {1, 2}};
      if (stem_ext == 1) mv.resize(1);
    }
    for (const auto& m : mv) {
      const int i = d - m.first, j = q + m.second;
      if (i < 0 || j > n - 1) continue;
      const double s = E.stem_interior(i, d, q, j);
      if (s < HALF_NEG) continue;
      ext_paths(i, j, t + 1, w + static_cast<long double>(s) + pair_bonus, out);
    }
  }

  double stem_ext_naive(int p, int u, int v, double a) {
    std::vector<std::tuple<int, int, long double>> paths;
    ext_paths(p, stc_q, 0, 0.0L, paths);
    const long double w_pk = static_cast<long double>(E.w_beta1) + 2.0L * E.w_stem;
    std::vector<long double> terms;
    for (const auto& [d, q, w] : paths) {
      const double h = hr_h(u, v, q);
      if (h > HALF_NEG) {
        long double s = static_cast<long double>(prefix_log(d)) + a + w + h + hr_c(q) + w_pk;
        if (score_qz) s += qz_log(u + 1, v - 1);
        if (stem_consistent) s += static_cast<long double>(E.at(u, v)) + E.at(d, q);
        terms.push_back(s);
      }
      const double bg = stem_Bg_raw(d, u, v);
      if (bg > HALF_NEG) {
        long double s = static_cast<long double>(bg) + a + w + stem_->S[q + 1];
        if (stem_consistent) s += static_cast<long double>(E.at(d, q)) + E.w_beta1;
        terms.push_back(s);
      }
    }
    if (terms.empty()) return NEG_INF;
    long double m = *std::max_element(terms.begin(), terms.end()), sum = 0.0L;
    for (long double t : terms) sum += std::exp(t - m);
    return static_cast<double>(m + std::log(sum));
  }

  // ------------------------------------------------- helper scores (Qgls, Qgrs)
  // Stop term of a stem with final outer pair d.q2 and inner pair (u, v), without alpha
  // (DP03._helper_stop): logsumexp of C[d-1] Qz[u+1][v-1] H(u, v, q2) / C0[q2] w_pk and
  // Bg(d, u, v) S[q2+1] (consistent constants as stem_ext_score).
  double helper_stop(int d, int q2, int u, int v) {
    double sl = NEG_INF, sr = NEG_INF;
    const double h = stem_H_cached(u, v, q2);
    if (h > HALF_NEG) {
      double s = prefix_log(d);
      if (score_qz) s += qz_log(u + 1, v - 1);
      sl = s + h + hr_c(q2) + (E.w_beta1 + 2.0 * E.w_stem);
      if (stem_consistent) sl += (E.at(u, v) + E.at(d, q2));
    }
    const double bg = stem_Bg(d, u, v);
    if (bg > HALF_NEG) {
      sr = bg + stem_->S[q2 + 1];
      if (stem_consistent) sr += (E.at(d, q2) + E.w_beta1);
    }
    return log_add(sl, sr);
  }

  // O(x, y; u, v) (DP03._helper_O): the stem closed at x.y stops there or grows on (growth
  // nodes of x.y with their path weights), cached per step.
  double helper_O(int x, int y, int u, int v) {
    const uint64_t key = (static_cast<uint64_t>(x) << 48) | (static_cast<uint64_t>(y) << 32) |
                         (static_cast<uint64_t>(u) << 16) | static_cast<uint64_t>(v);
    auto it = stOc.find(key);
    if (it != stOc.end()) return it->second;
    double acc = helper_stop(x, y, u, v);
    const std::vector<ExtNode>& xs = stem_grow_nodes(x, y);
    for (const ExtNode& nd : xs) acc = log_add(acc, nd.w + helper_stop(nd.d, nd.q, u, v));
    hlp_nodes += static_cast<long long>(xs.size());
    stOc.emplace(key, acc);
    return acc;
  }

  // A closing of a helper state into a Qg with outer pair x.y and closing weight w (alpha
  // not included); filled by helper_closings.
  struct Closing {
    int x, y;
    double w;
  };

  // Closings of Qgls[q] state (i, u, v): case (f) into outer pair (i-1, j), j = q+1 ..
  // min(n-1, q+1+stem_ml_span) (seed group v == q: run capped as in fill_Qg).
  void closings_ls(int i, int v, int q, std::vector<Closing>& out) const {
    out.clear();
    const int x = i - 1;
    if (x < 0) return;
    const int j_hi = std::min(n - 1, q + 1 + stem_ml_span);
    for (int j = q + 1; j <= j_hi; ++j) {
      const int unpaired = j - q - 1;
      if (v == q && unpaired > cap) continue;
      if (!E.pairs(x, j) || !E.end_ok(x, j)) continue;
      double w_ext = E.w_ap + E.w_bp + E.w_cp * unpaired + pair_bonus;
      if (dz) w_ext += E.run(q + 1, j - 1);
      out.push_back({x, j, w_ext + E.at(x, j)});
    }
  }

  // Closings of Qgrs[q] state (c, u, v) into outer pair (i, q+1): case (g) (unpaired run on
  // the left, i ascending), then case (h) (Qm[i+1][c-1] on the left, i ascending).
  void closings_rs(int c, int q, std::vector<Closing>& out) const {
    out.clear();
    const int j = q + 1;
    if (j > n - 1) return;
    for (int i = lead_i_lo(c); i < c; ++i) {
      if (!E.pairs(i, j) || !E.end_ok(i, j)) continue;
      double w = E.w_ap + E.w_bp + E.w_cp * (c - i - 1) + E.at(i, j) + pair_bonus;
      if (dz) w += E.run(i + 1, c - 1);
      out.push_back({i, j, w});
    }
    if (c - 1 < 0) return;
    const OMap& at = qmp_at(c - 1);
    std::vector<std::pair<int, double>> qm;
    qm_pair_slots(c - 1, j, [&](int r) {
      if (at.keys[r] >= 0) qm.emplace_back(at.keys[r], at.vals[r]);
    });
    std::sort(qm.begin(), qm.end());
    for (const auto& [i, val] : qm) {
      const int outer = i - 1;
      if (outer < 0 || outer >= c || !E.pairs(outer, j) || !E.end_ok(outer, j)) continue;
      out.push_back({outer, j, val + E.w_ap + E.w_bp + E.at(outer, j) + pair_bonus});
    }
  }

  // sigma of a helper state with inner pair (u, v) and log weight a: logsumexp over its
  // closings of w + a + O(x, y; u, v).
  double helper_sigma(const std::vector<Closing>& cl, int u, int v, double a) {
    double acc = NEG_INF;
    for (const Closing& c : cl) acc = log_add(acc, c.w + a + helper_O(c.x, c.y, u, v));
    return acc;
  }

  // --stem-ext-check: the helper score by walking every closing and growth path separately
  // (H and Bg recomputed), every term summed in long double.
  double helper_naive(const std::vector<Closing>& cl, int u, int v, double a) {
    std::vector<long double> terms;
    std::vector<std::tuple<int, int, long double>> paths;
    const long double w_pk = static_cast<long double>(E.w_beta1) + 2.0L * E.w_stem;
    for (const Closing& c : cl) {
      paths.clear();
      paths.emplace_back(c.x, c.y, 0.0L);
      ext_paths(c.x, c.y, 0, 0.0L, paths);
      for (const auto& [d, q, w] : paths) {
        const long double base = static_cast<long double>(c.w) + a + w;
        const double h = hr_h(u, v, q);
        if (h > HALF_NEG) {
          long double s = base + prefix_log(d) + h + hr_c(q) + w_pk;
          if (score_qz) s += qz_log(u + 1, v - 1);
          if (stem_consistent) s += static_cast<long double>(E.at(u, v)) + E.at(d, q);
          terms.push_back(s);
        }
        const double bg = stem_Bg_raw(d, u, v);
        if (bg > HALF_NEG) {
          long double s = base + bg + stem_->S[q + 1];
          if (stem_consistent) s += static_cast<long double>(E.at(d, q)) + E.w_beta1;
          terms.push_back(s);
        }
      }
    }
    if (terms.empty()) return NEG_INF;
    long double m = *std::max_element(terms.begin(), terms.end()), sum = 0.0L;
    for (long double t : terms) sum += std::exp(t - m);
    return static_cast<double>(m + std::log(sum));
  }

  // --stem-est-helpers nested: log sum over the closings (x, y) of w + beta_nested(x, y), in
  // closing order.  The closings depend only on (i, v == q) for Qgls and on c for Qgrs, so the
  // sum is cached per step; returns it and sets ncl to the number of closings.
  double helper_nested_sum(bool ls, int p, int v, int q, int& ncl) {
    if (hn_q != q) {
      hn_q = q;
      hnL.assign(2 * static_cast<size_t>(n), {NEG_INF, -1});
      hnR.assign(static_cast<size_t>(n), {NEG_INF, -1});
    }
    std::pair<double, int>& slot = ls ? hnL[2 * static_cast<size_t>(p) + (v == q)] : hnR[p];
    if (slot.second < 0) {
      std::vector<Closing> cl;
      if (ls) closings_ls(p, v, q, cl);
      else closings_rs(p, q, cl);
      double acc = NEG_INF;
      for (const Closing& c : cl) acc = log_add(acc, c.w + stem_->B(c.x, c.y));
      slot = {acc, static_cast<int>(cl.size())};
    }
    ncl = slot.second;
    return slot.first;
  }

  // sigma of a helper state (nested): alpha + its closing sum + the gap filler, as four_score.
  double helper_sigma_nested(bool ls, const Ent& e, int v, int q, int& ncl) {
    const double s = helper_nested_sum(ls, e.p, v, q, ncl);
    if (s < HALF_NEG) return NEG_INF;
    double r = e.a + s;
    if (score_qz) r += qz_log(e.u + 1, v - 1);
    return r;
  }

  // --stem-ext-check (nested): the same score without the per-step cache, in long double.
  double helper_nested_naive(bool ls, const Ent& e, int v, int q) {
    std::vector<Closing> cl;
    if (ls) closings_ls(e.p, v, q, cl);
    else closings_rs(e.p, q, cl);
    std::vector<long double> terms;
    for (const Closing& c : cl) {
      const double b = stem_->B(c.x, c.y);
      if (b < HALF_NEG) continue;
      long double t = static_cast<long double>(e.a) + c.w + b;
      if (score_qz) t += qz_log(e.u + 1, v - 1);
      terms.push_back(t);
    }
    if (terms.empty()) return NEG_INF;
    long double m = *std::max_element(terms.begin(), terms.end()), sum = 0.0L;
    for (long double t : terms) sum += std::exp(t - m);
    return static_cast<double>(m + std::log(sum));
  }

  // ------------------------------------------------ profiles (--stem-est-profile)
  // Sparse profile accumulator: entries in first-reached order, each summed in log space in the
  // order added; finish() sorts them by index (DP03._ProfAcc).
  struct ProfAcc {
    std::unordered_map<int, size_t> li, ri;
    std::vector<std::pair<int, double>> l, r;
    static void add(std::unordered_map<int, size_t>& ix, std::vector<std::pair<int, double>>& v,
                    int key, double w) {
      auto res = ix.emplace(key, v.size());
      if (res.second) v.emplace_back(key, w);
      else v[res.first->second].second = log_add(v[res.first->second].second, w);
    }
    void addL(int y, double w) { add(li, l, y, w); }
    void addR(int x, double w) { add(ri, r, x, w); }
    void finish(Prof& p) {
      auto by_ix = [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
        return a.first < b.first;
      };
      std::sort(l.begin(), l.end(), by_ix);
      std::sort(r.begin(), r.end(), by_ix);
      for (const auto& e : l) {
        p.ly.push_back(e.first);
        p.lw.push_back(e.second);
      }
      for (const auto& e : r) {
        p.rx.push_back(e.first);
        p.rw.push_back(e.second);
      }
    }
  };

  void pf_note_step() {
    if (pf_q != stc_q) {
      pf_q = stc_q;
      ++pf_steps;
    }
  }

  // Profile of the stem closed at outer pair x.y (DP03._pair_prof): the stop-here node (x, y)
  // with G = 1, then its growth nodes in list order.  Node (d', q', G) adds
  // G C[d'-1] w_pk / C0[q'] to L at q' and G S[q'+1] to R at d' (consistent: at(d', q') to
  // both and beta1 to R).  Untruncated; cached per step.
  const Prof& pair_prof(int x, int y) {
    if (stem_chain) return chain_root(x, y);
    const uint64_t key = key2(x, y);
    auto it = stP.find(key);
    if (it != stP.end()) return it->second;
    const double w_pk = E.w_beta1 + 2.0 * E.w_stem;
    ProfAcc acc;
    auto add = [&](int d, int q2, double g) {
      double cl = ((g + prefix_log(d)) + hr_c(q2)) + w_pk;
      double cr = g + stem_->S[q2 + 1];
      if (stem_consistent) {
        cl += E.at(d, q2);
        cr += (E.at(d, q2) + E.w_beta1);
      }
      acc.addL(q2, cl);
      acc.addR(d, cr);
    };
    add(x, y, 0.0);
    for (const ExtNode& nd : stem_grow_nodes(x, y)) add(nd.d, nd.q, nd.w);
    Prof p;
    acc.finish(p);
    ++pf_pairs;
    return stP.emplace(key, std::move(p)).first->second;
  }

  // ------------------------------------------------ arm chains (--stem-est-chain)
  // A step from pair x.y out to a kept pair x'.y' enclosing it, with log weight w; t = 0 (b)
  // interior, 1 (f), 2 (g), 3 (h) multiloop.
  struct ChainNbr {
    int x, y, t;
    double w;
  };
  void count_steps(int which, const std::vector<ChainNbr>& nb) {
    ++st_nodes[which];
    for (const ChainNbr& s : nb) ++st_kind[which][s.t == 0 ? 0 : 1];
  }
  bool stem_chain = false, chain_ml = true;
  int chain_span = 150, chain_loop = 30;
  double chain_cut = 1e-8;
  // chain node profiles Prof*(x, y) and root profiles, by key2; valid for the whole run
  std::unordered_map<uint64_t, Prof> chM, chR;

  // Steps out of x.y (DP03._chain_nbrs): (b) x' = x-1 down to x-1-Lc, y' ascending in
  // (y, y+1+Lc], Lc = min(chain_loop, cap); with chain_ml and end_ok(x, y): (f) x' = i-1 over
  // the keys i of qmp_at(x-1) ascending (x - i <= S), y' in (y, y+1+ml_span]; (g) x' ascending
  // in [x-1-cap, x), y' with mbar(y+1, y'-1) and y'-y-1 <= S; (h) x' as (f), y' as (g).
  // x'.y' kept (and end_ok for f-h).
  void chain_nbrs(int x, int y, std::vector<ChainNbr>& out) {
    out.clear();
    const StemEst& st = *stem_;
    const int lc = std::min(chain_loop, cap);
    for (int xo = x - 1; xo >= std::max(0, x - 1 - lc); --xo) {
      const std::vector<int>& kr = st.kept_right[xo];
      for (auto it = std::upper_bound(kr.begin(), kr.end(), y);
           it != kr.end() && *it <= y + 1 + lc; ++it) {
        const double s = E.stem_interior(xo, x, y, *it);
        if (s < HALF_NEG) continue;
        out.push_back({xo, *it, 0, s + pair_bonus});
      }
    }
    if (!chain_ml || !E.end_ok(x, y)) return;
    const double br = E.w_bp + E.at(x, y);
    std::vector<std::pair<int, double>> qm;
    if (x >= 1) {
      const OMap& at = qmp_at(x - 1);
      for (size_t r = 0; r < at.slots(); ++r) {
        const int i = at.keys[r];
        if (i >= 1 && x - i <= chain_span) qm.emplace_back(i - 1, at.vals[r]);
      }
      std::sort(qm.begin(), qm.end());
    }
    auto close = [&](int xo, int yo) { return (E.w_ap + E.w_bp + E.at(xo, yo)) + pair_bonus; };
    for (const auto& [xo, qv] : qm) {
      const std::vector<int>& kr = st.kept_right[xo];
      for (auto it = std::upper_bound(kr.begin(), kr.end(), y);
           it != kr.end() && *it <= y + 1 + stem_ml_span; ++it) {
        const int yo = *it;
        if (!E.end_ok(xo, yo)) continue;
        double w = qv + br + E.w_cp * (yo - y - 1);
        if (dz) w += E.run(y + 1, yo - 1);
        out.push_back({xo, yo, 1, w + close(xo, yo)});
      }
    }
    for (int xo = ml_run_lo(x); xo < x; ++xo) {
      double rl = E.w_cp * (x - xo - 1);
      if (dz) rl += E.run(xo + 1, x - 1);
      const std::vector<int>& kr = st.kept_right[xo];
      for (auto it = std::upper_bound(kr.begin(), kr.end(), y + 1);
           it != kr.end() && *it - y - 1 <= chain_span; ++it) {
        const int yo = *it;
        const double mb = st.Mbar(y + 1, yo - 1);
        if (mb < HALF_NEG || !E.end_ok(xo, yo)) continue;
        out.push_back({xo, yo, 2, ((mb + br) + rl) + close(xo, yo)});
      }
    }
    for (const auto& [xo, qv] : qm) {
      const std::vector<int>& kr = st.kept_right[xo];
      for (auto it = std::upper_bound(kr.begin(), kr.end(), y + 1);
           it != kr.end() && *it - y - 1 <= chain_span; ++it) {
        const int yo = *it;
        const double mb = st.Mbar(y + 1, yo - 1);
        if (mb < HALF_NEG || !E.end_ok(xo, yo)) continue;
        out.push_back({xo, yo, 3, ((qv + mb) + br) + close(xo, yo)});
      }
    }
  }

  // Stop-here entries of node x.y with G = 1 (as pair_prof's add).
  void chain_stop(ProfAcc& acc, int x, int y) const {
    double cl = (prefix_log(x) + hr_c(y)) + (E.w_beta1 + 2.0 * E.w_stem);
    double cr = stem_->S[y + 1];
    if (stem_consistent) {
      cl += E.at(x, y);
      cr += (E.at(x, y) + E.w_beta1);
    }
    acc.addL(y, cl);
    acc.addR(x, cr);
  }

  // Adds w + Prof*(x', y') of every step in nb (step order; all in chM); counts the edges.
  void chain_merge(ProfAcc& acc, const std::vector<ChainNbr>& nb) {
    for (const ChainNbr& t : nb) {
      const Prof& q = chM.find(key2(t.x, t.y))->second;
      for (size_t k = 0; k < q.ly.size(); ++k) acc.addL(q.ly[k], t.w + q.lw[k]);
      for (size_t k = 0; k < q.rx.size(); ++k) acc.addR(q.rx[k], t.w + q.rw[k]);
      ++ch_edges[t.t];
    }
  }

  // Drop entries below log(chain_cut) + the part's max (order kept).
  long long chain_trunc(std::vector<int>& ix, std::vector<double>& w) const {
    if (chain_cut <= 0.0 || w.empty()) return 0;
    const double thr = *std::max_element(w.begin(), w.end()) + std::log(chain_cut);
    size_t o = 0;
    for (size_t k = 0; k < w.size(); ++k)
      if (w[k] >= thr) {
        ix[o] = ix[k];
        w[o] = w[k];
        ++o;
      }
    const long long d = static_cast<long long>(w.size() - o);
    ix.resize(o);
    w.resize(o);
    return d;
  }

  // Prof*(x, y) = [end_ok(x, y)] stop(x, y) + sum over the steps of w Prof*(x', y'), cut;
  // built bottom-up with an explicit stack (steps strictly enclose, so no cycles).
  void chain_build(int x, int y) {
    if (chM.count(key2(x, y))) return;
    std::vector<std::pair<int, int>> stk{{x, y}};
    std::vector<ChainNbr> nb;
    while (!stk.empty()) {
      const auto [a, b] = stk.back();
      if (chM.count(key2(a, b))) {
        stk.pop_back();
        continue;
      }
      chain_nbrs(a, b, nb);
      bool ready = true;
      for (const ChainNbr& t : nb)
        if (!chM.count(key2(t.x, t.y))) {
          stk.emplace_back(t.x, t.y);
          ready = false;
        }
      if (!ready) continue;
      ProfAcc acc;
      if (E.end_ok(a, b)) chain_stop(acc, a, b);
      chain_merge(acc, nb);
      Prof p;
      acc.finish(p);
      ch_cut_drop += chain_trunc(p.ly, p.lw) + chain_trunc(p.rx, p.rw);
      ++ch_nodes;
      chM.emplace(key2(a, b), std::move(p));
      stk.pop_back();
    }
  }

  // Root profile of the stem closed at x.y (DP03._chain_root): stop(x, y) unconditionally plus
  // the steps out of x.y; untruncated, cached for the run.
  const Prof& chain_root(int x, int y) {
    const uint64_t key = key2(x, y);
    auto it = chR.find(key);
    if (it != chR.end()) return it->second;
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<ChainNbr> nb;
    chain_nbrs(x, y, nb);
    for (const ChainNbr& t : nb) chain_build(t.x, t.y);
    ProfAcc acc;
    chain_stop(acc, x, y);
    chain_merge(acc, nb);
    Prof p;
    acc.finish(p);
    ++ch_roots;
    ch_ly += static_cast<long long>(p.ly.size());
    ch_rx += static_cast<long long>(p.rx.size());
    ch_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return chR.emplace(key, std::move(p)).first->second;
  }

  // Every chain path from x.y (root: stop unconditional; later nodes need end_ok) with its
  // weight summed in long double.
  void chain_paths(int x, int y, bool root, long double w,
                   std::vector<std::tuple<int, int, long double>>& out) {
    if (root || E.end_ok(x, y)) out.emplace_back(x, y, w);
    std::vector<ChainNbr> nb;
    chain_nbrs(x, y, nb);
    for (const ChainNbr& t : nb) chain_paths(t.x, t.y, false, w + t.w, out);
  }

  // --stem-profile-check with chain: the score of a state with inner pair (u, v) and log
  // weight a closed by cl, every chain path's stop terms summed separately (as helper_naive);
  // the path list of each root is enumerated once per run.
  std::unordered_map<uint64_t, std::vector<std::tuple<int, int, long double>>> chPaths;
  double chain_naive(const std::vector<Closing>& cl, int u, int v, double a) {
    std::vector<long double> terms;
    const long double w_pk = static_cast<long double>(E.w_beta1) + 2.0L * E.w_stem;
    for (const Closing& c : cl) {
      auto it = chPaths.find(key2(c.x, c.y));
      if (it == chPaths.end()) {
        it = chPaths.emplace(key2(c.x, c.y), std::vector<std::tuple<int, int, long double>>())
                 .first;
        chain_paths(c.x, c.y, true, 0.0L, it->second);
      }
      for (const auto& [d, q, w] : it->second) {
        const long double base = static_cast<long double>(c.w) + a + w;
        const double h = hr_h(u, v, q);
        if (h > HALF_NEG) {
          long double s = base + prefix_log(d) + h + hr_c(q) + w_pk;
          if (score_qz) s += qz_log(u + 1, v - 1);
          if (stem_consistent) s += static_cast<long double>(E.at(u, v)) + E.at(d, q);
          terms.push_back(s);
        }
        const double bg = stem_Bg_raw(d, u, v);
        if (bg > HALF_NEG) {
          long double s = base + bg + stem_->S[q + 1];
          if (stem_consistent) s += static_cast<long double>(E.at(d, q)) + E.w_beta1;
          terms.push_back(s);
        }
      }
    }
    if (terms.empty()) return NEG_INF;
    long double m = *std::max_element(terms.begin(), terms.end()), sum = 0.0L;
    for (long double t : terms) sum += std::exp(t - m);
    return static_cast<double>(m + std::log(sum));
  }

  // ------------------------------------------ principled score (--stem-est-principled)
  // sigma(X) = alpha_X + log Out_hat(X); Out_hat follows DP03's edges outward from X: exact
  // scan tables at positions <= q, the nested pre-pass (StemEst pr*) beyond q, and outside
  // the whole pseudoknot Ctx(i, j) = e^beta1 C[i-1] Sbar[j+1] + e^beta1m Bml[i][j].
  // Stem steps from a pair outward (pr_steps) run over the pre-pass kept pairs.  The growth G
  // of the scored stem takes Qm from the scan on its left; future stems H (profiles over
  // their outer ends) take Mbar on both sides.
  bool pr_on = false, pr_check = false, pr_nested = false, pr_ml = true, pr_built = false;
  int pr_loop = 30, pr_ml_span = 150, pr_ml_close = 30, pr_arms_k = 64;
  bool pr_ratio = false, pr_gapfree = false, pr_cut_roots = false;
  double pr_arms_eps = 1e-3;
  int pr_arms_max = 256;
  double pr_cut = 1e-8;
  // --pr-arms sep: the exterior term of the L-stem arm part over all arms, the gap filler made
  // separable: Qz[u+1][v-1] S_q'(u, v), S = sum over arms xs >= u+1 (two: = u+1), x' <= v-1
  // of e^phi / Zarm (arm_z 0: Qz[xs][x'], 1: e^(wu (x' - xs + 1))).  The multiloop-context
  // term (Bml(p', ys), not one factor per arm) keeps the exact-placement arms of pr_select.
  // qz_unp: the score's gap filler unpaired-only (test); sep_check: every distinct query
  // against brute force and exact placement.
  int pr_arms_mode = 0, pr_arm_z = 0;
  bool pr_qz_unp = false, pr_sep_check = false;
  bool pr_gsel = false, pr_gcheck = false;
  double pr_geps = 1e-6;
  int pr_gB = 512;
  FILE* pr_dump = nullptr;

 public:
  // Growth instrumentation, per profile kind (0 node, 1 root, 2 key): entries before / after
  // selection, distinct p' and q' after (and before) selection; under --grow-rec the nonzero
  // q' of g_L and p' of g_R (all built layers).  Memory in bytes of the stored growth data:
  // nodes and roots (kept for the run), helpers (cleared at PR_HG_MAX), the per-step root /
  // key caches; peaks over the run.
  bool rec_noml = false;
  std::vector<int> gz_before[3], gz_after[3], gz_dp[3], gz_dq[3], gz_dpb[3], gz_dqb[3];
  int gs_who = 0;
  double gm_nodes = 0, gm_roots = 0, gm_help = 0, gm_help_peak = 0, gm_step = 0, gm_step_peak = 0,
         gm_total_peak = 0;
  long long gm_help_entries_peak = 0;
  void gm_update() {
    gm_help_peak = std::max(gm_help_peak, gm_help);
    gm_step_peak = std::max(gm_step_peak, gm_step);
    gm_total_peak = std::max(gm_total_peak, gm_nodes + gm_roots + gm_help + gm_step);
  }
  std::vector<int> gz_mark;
  int gz_stamp = 0;
  std::pair<int, int> gz_distinct(const std::vector<uint64_t>& k) {
    if (gz_mark.size() < static_cast<size_t>(n)) gz_mark.assign(n, 0);
    ++gz_stamp;
    int dp = 0, dq = 0, last = -1;
    for (uint64_t key : k) {
      const int a = k_hi(key), b = k_lo(key);
      if (a != last) {
        ++dp;
        last = a;
      }
      if (gz_mark[b] != gz_stamp) {
        gz_mark[b] = gz_stamp;
        ++dq;
      }
    }
    return {dp, dq};
  }

 private:
  double pr_qz(int a, int b) const {
    if (qz_audit && a <= b && b < n && b > l2r_done[LT_QZ]) ++qz_unf_pr;
    if (!pr_qz_unp) return qz_log(a, b);
    if (a > b) return 0.0;
    return E.wu * (b - a + 1);
  }
  // Sparse log profile over key2(a, b), ascending.  --pr-fast: lin[k] = e^(w[k] - mx) as well
  // (filled by pr_linify for the profiles that are merged into others).
  struct PProf {
    std::vector<uint64_t> k;
    std::vector<double> w;
    double mx = NEG_INF;
    std::vector<double> lin;
    double tot = NEG_INF;  // log sum of w, set by pr_linify
  };
  // --sweep-check: double V profiles (pass 0) and the pair list, saved before the float readback.
  std::unordered_map<uint64_t, PProf> sw_ref;
  std::vector<std::pair<int, int>> sw_pairs;
  // --sweep: gamma of each pair, without the prefix. 0 = full right gap, 1 = gap <= loop cap.
  std::unordered_map<uint64_t, std::unordered_map<uint64_t, double>> sw_ga[2];
  bool sw_filled = false;
  static double pbytes(const PProf& p) {
    return 8.0 * static_cast<double>(p.k.capacity() + p.w.capacity() + p.lin.capacity()) +
           static_cast<double>(sizeof(PProf) + 24);
  }
  static void pr_linify(PProf& p) {
    p.mx = NEG_INF;
    for (double x : p.w) p.mx = std::max(p.mx, x);
    p.lin.resize(p.w.size());
    double s = 0.0;
    for (size_t k = 0; k < p.w.size(); ++k) s += p.lin[k] = std::exp(p.w[k] - p.mx);
    p.tot = s > 0.0 ? p.mx + std::log(s) : NEG_INF;
  }

  // --pr-fast: dense accumulator over pair keys a * n + b (key2 order), sums in linear space
  // relative to the log scale dRef (the largest term added, set by dacc_begin).
 public:
  long long pr_adds = 0;  // --pr-fast: profile entries added into the dense accumulator
  int pr_fast = 0;  // 0 off, 1 dense accumulator, 2 also the multiloop helper profiles
  long long pr_hbuilt[2] = {}, pr_hentries[2] = {}, pr_hterms[2] = {}, pr_nterms[2] = {};
  long long pr_hg_ent = 0, pr_hg_clears = 0;
  // --pr-abl: bits of PA_* (couplings of sigma_pr replaced by separable forms); the scored
  // state's p and q (set by pr_score4).
  static constexpr int PA_BML = 1, PA_Q = 2, PA_P = 4, PA_XS = 8, PA_KS = 16, PA_SEP = 32,
                       PA_XB = 64, PA_QR = 128;
  // log Qz[q+1][q'] as --pr-abl q uses it (qr: the prefix ratio Qz[0][q'] / Qz[0][q])
  double pr_abl_qz(int q, int qp) const {
    if (qp <= q) return 0.0;
    if (!(pr_abl & PA_QR)) return pr_qz(q + 1, qp);
    const double a = qz_log(0, qp), b = qz_log(0, q);
    return a > HALF_NEG && b > HALF_NEG ? a - b : NEG_INF;
  }
  int pr_abl = 0, pr_abl_p = -1, pr_abl_q = -1, pr_abl_pk = 0;
  long long pr_sep_evals = 0, pr_sep_pp = 0, pr_sep_ppk = 0;

 private:
  std::vector<double> dA;
  std::vector<char> dMark;
  std::vector<uint32_t> dT;
  double dRef = 0.0;
  void dacc_begin(double ref) {
    if (dA.empty()) {
      dA.assign(static_cast<size_t>(n) * n, 0.0);
      dMark.assign(static_cast<size_t>(n) * n, 0);
    }
    dRef = ref;
  }
  uint32_t dslot(uint64_t key) const {
    return static_cast<uint32_t>(k_hi(key)) * static_cast<uint32_t>(n) + static_cast<uint32_t>(k_lo(key));
  }
  void dacc_add1(uint64_t key, double w) {
    if (w < HALF_NEG) return;
    const uint32_t s = dslot(key);
    if (!dMark[s]) {
      dMark[s] = 1;
      dT.push_back(s);
    }
    dA[s] += std::exp(w - dRef);
  }
  void dacc_add(const PProf& q, double w) {
    if (q.k.empty() || q.mx < HALF_NEG || w < HALF_NEG) return;
    const double sc = std::exp((w + q.mx) - dRef);
    if (sc == 0.0) return;
    const size_t m = q.k.size();
    const uint64_t* kk = q.k.data();
    const double* ll = q.lin.data();
    pr_adds += static_cast<long long>(m);
    for (size_t k = 0; k < m; ++k) {
      const uint32_t s = dslot(kk[k]);
      if (!dMark[s]) {
        dMark[s] = 1;
        dT.push_back(s);
      }
      dA[s] += sc * ll[k];
    }
  }
  // As PAcc::finish (key order, cut on w + adj(key) against the max), then resets.
  template <class F>
  long long dacc_finish(PProf& p, double cut, F adj) {
    std::sort(dT.begin(), dT.end());
    double thr = NEG_INF;
    if (cut > 0.0 && !dT.empty()) {
      double mx = NEG_INF;
      for (uint32_t s : dT)
        if (dA[s] > 0.0) mx = std::max(mx, dRef + std::log(dA[s]) + adj(dkey(s)));
      thr = mx + std::log(cut);
    }
    long long d = 0;
    p.k.reserve(dT.size());
    p.w.reserve(dT.size());
    for (uint32_t s : dT) {
      if (dA[s] > 0.0) {
        const uint64_t key = dkey(s);
        const double w = dRef + std::log(dA[s]);
        if (w + adj(key) >= thr) {
          p.k.push_back(key);
          p.w.push_back(w);
        } else {
          ++d;
        }
      }
      dA[s] = 0.0;
      dMark[s] = 0;
    }
    dT.clear();
    return d;
  }
  long long dacc_finish(PProf& p, double cut) {
    return dacc_finish(p, cut, [](uint64_t) { return 0.0; });
  }
  uint64_t dkey(uint32_t s) const {
    return key2(static_cast<int>(s / static_cast<uint32_t>(n)), static_cast<int>(s % static_cast<uint32_t>(n)));
  }
  long long pr_finish_out_d(PProf& p, double cut) {
    if (!pr_gapfree) return dacc_finish(p, cut);
    const StemEst& st = *stem_;
    return dacc_finish(p, cut, [&](uint64_t k) { return prefix_log(k_hi(k)) + st.prS[k_lo(k) + 1]; });
  }
  long long pr_finish_fut_d(PProf& p, int x) {
    if (!pr_gapfree) return dacc_finish(p, pr_cut);
    const StemEst& st = *stem_;
    return dacc_finish(p, pr_cut, [&](uint64_t k) {
      return (pr_vxs ? 0.0 : st.prS[k_lo(k) + 1]) - E.wu * (x - k_hi(k));
    });
  }
  // Entries summed in log space in the order added; finish() sorts them by key and drops
  // those below max + log(cut) (cut > 0); returns the drops.
  struct PAcc {
    std::unordered_map<uint64_t, size_t> ix;
    std::vector<std::pair<uint64_t, double>> v;
    void add(uint64_t key, double w) {
      if (w < HALF_NEG) return;
      auto res = ix.emplace(key, v.size());
      if (res.second) v.emplace_back(key, w);
      else v[res.first->second].second = log_add(v[res.first->second].second, w);
    }
    long long finish(PProf& p, double cut) {
      return finish(p, cut, [](uint64_t) { return 0.0; });
    }
    // cut on w + adj(key) (the entry's contribution), keeping w
    template <class F>
    long long finish(PProf& p, double cut, F adj) {
      std::sort(v.begin(), v.end(), [](const std::pair<uint64_t, double>& a,
                                       const std::pair<uint64_t, double>& b) {
        return a.first < b.first;
      });
      double thr = NEG_INF;
      if (cut > 0.0 && !v.empty()) {
        double mx = NEG_INF;
        for (const auto& e : v) mx = std::max(mx, e.second + adj(e.first));
        thr = mx + std::log(cut);
      }
      long long d = 0;
      for (const auto& e : v)
        if (e.second + adj(e.first) >= thr) {
          p.k.push_back(e.first);
          p.w.push_back(e.second);
        } else {
          ++d;
        }
      return d;
    }
  };
  static int k_hi(uint64_t k) { return static_cast<int>(k >> 32); }
  static int k_lo(uint64_t k) { return static_cast<int>(k & 0xffffffffu); }

  // Profile cuts.  Outer-pair profiles (growth, roots, keys; entries p'.q'): gap-free factor
  // C[p'-1] Sbar[q'+1].  Future-stem profiles of a pair with left end x (entries xs.ys, the
  // left part xs..x-1 inside a gap): Sbar[ys+1] / e^(wu (x - xs)).
  long long pr_finish_out(PAcc& acc, PProf& p, double cut) const {
    if (!pr_gapfree) return acc.finish(p, cut);
    const StemEst& st = *stem_;
    return acc.finish(p, cut, [&](uint64_t k) { return prefix_log(k_hi(k)) + st.prS[k_lo(k) + 1]; });
  }
  long long pr_finish_fut(PAcc& acc, PProf& p, int x) const {
    if (!pr_gapfree) return acc.finish(p, pr_cut);
    const StemEst& st = *stem_;
    return acc.finish(p, pr_cut, [&](uint64_t k) {
      return st.prS[k_lo(k) + 1] - E.wu * (x - k_hi(k));
    });
  }

  // DP03 stem steps out of x.y onto kept pairs x'.y' (the engine's weights): (b) interior,
  // x - x' - 1 and y' - y - 1 <= min(pr_loop, cap); with pr_ml and end_ok(x, y): (f) Qm on the
  // left (scan: qmp_at(x-1); else Mbar), unpaired run on the right (<= cap from a single pair,
  // the engine's seed cap; else <= S); (g) unpaired run (<= cap) left, Mbar right; (h) Qm left,
  // Mbar right; multiloop parts of span <= S = pr_ml_span, x'.y' end_ok for f-h.  Not scan
  // (future stems / partners): onto the partner pairs if --pr-partner-floor is set.  ml = false:
  // (b) only.
  // Pair list for a pass-2 step. Flags off: the current kept pairs (partner list when set).
  const std::vector<std::vector<int>>& pr_pair_list(bool scan) const {
    const StemEst& st = *stem_;
    if (!scan && p2_v_mode > 0) return st.prKeptRightV;
    if (scan && p2_g_mode > 0) return st.prKeptRightG;
    if (!scan && st.pr_partner) return st.prKeptRightP;
    return st.prKeptRight;
  }
  void pr_steps(int x, int y, bool inner, bool scan, std::vector<ChainNbr>& out, bool ml = true) const {
    out.clear();
    const StemEst& st = *stem_;
    const std::vector<std::vector<int>>& KR = pr_pair_list(scan);
    const int lc = std::min(pr_loop, cap);
    for (int xo = x - 1; xo >= std::max(0, x - 1 - lc); --xo) {
      const std::vector<int>& kr = KR[xo];
      for (auto it = std::upper_bound(kr.begin(), kr.end(), y);
           it != kr.end() && *it <= y + 1 + lc; ++it) {
        const double s = E.stem_interior(xo, x, y, *it);
        if (s < HALF_NEG) continue;
        out.push_back({xo, *it, 0, s + pair_bonus});
      }
    }
    if (!ml || !pr_ml || (!scan && !pr_partner_ml) || (scan && rec_noml) || !E.end_ok(x, y)) return;
    const int S = pr_ml_span;
    const double br = E.w_bp + E.at(x, y);
    std::vector<std::pair<int, double>> qm;
    if (x >= 1) {
      if (scan) {
        const OMap& at = qmp_at(x - 1);
        for (size_t r = 0; r < at.slots(); ++r) {
          const int i = at.keys[r];
          if (i >= 1 && x - i <= S) qm.emplace_back(i - 1, at.vals[r]);
        }
        std::sort(qm.begin(), qm.end());
      } else {
        for (int i = std::max(1, x - S); i <= x - 1; ++i) {
          const double mb = st.PrMbar(i, x - 1);
          if (mb > HALF_NEG) qm.emplace_back(i - 1, mb);
        }
      }
    }
    auto close = [&](int xo, int yo) { return (E.w_ap + E.w_bp + E.at(xo, yo)) + pair_bonus; };
    const int rf = inner ? cap : S;
    for (const auto& [xo, qv] : qm) {
      const std::vector<int>& kr = KR[xo];
      for (auto it = std::upper_bound(kr.begin(), kr.end(), y);
           it != kr.end() && *it - y - 1 <= rf; ++it) {
        const int yo = *it;
        if (!E.end_ok(xo, yo)) continue;
        double w = qv + br + E.w_cp * (yo - y - 1);
        if (dz) w += E.run(y + 1, yo - 1);
        out.push_back({xo, yo, 1, w + close(xo, yo)});
      }
    }
    for (int xo = ml_run_lo(x); xo < x; ++xo) {
      double rl = E.w_cp * (x - xo - 1);
      if (dz) rl += E.run(xo + 1, x - 1);
      const std::vector<int>& kr = KR[xo];
      for (auto it = std::upper_bound(kr.begin(), kr.end(), y + 1);
           it != kr.end() && *it - y - 1 <= S; ++it) {
        const int yo = *it;
        const double mb = st.PrMbar(y + 1, yo - 1);
        if (mb < HALF_NEG || !E.end_ok(xo, yo)) continue;
        out.push_back({xo, yo, 2, ((mb + br) + rl) + close(xo, yo)});
      }
    }
    for (const auto& [xo, qv] : qm) {
      const std::vector<int>& kr = KR[xo];
      for (auto it = std::upper_bound(kr.begin(), kr.end(), y + 1);
           it != kr.end() && *it - y - 1 <= S; ++it) {
        const int yo = *it;
        const double mb = st.PrMbar(y + 1, yo - 1);
        if (mb < HALF_NEG || !E.end_ok(xo, yo)) continue;
        out.push_back({xo, yo, 3, ((qv + mb) + br) + close(xo, yo)});
      }
    }
  }

  // --pr-fast 2: pr_steps' multiloop steps grouped by the new outer left end xo.  The steps
  // of x.y are br(x, y) times
  //   sum over the Qm keys (xo, qv) of pr_steps (f, h):  qv HF(xo, y, rf)
  //   sum over xo in [x-1-cap, x-1] (g):                 e^(wcp (x-xo-1)) run(xo+1, x-1) HM(xo, y)
  // with the right-side sums over the kept partners yo of xo (end_ok), C the child profile:
  //   HM(xo, y)     = sum_{yo > y+1, yo-y-1 <= S} Mbar(y+1, yo-1) cl(xo, yo) C(xo, yo)
  //   HF(xo, y, rf) = HM(xo, y) + sum_{yo > y, yo-y-1 <= rf} e^(wcp (yo-y-1)) run(y+1, yo-1) cl C
  // HF and HM depend on x only through xo, so every node with the same y shares them.
  struct MlSplit {
    bool on = false;
    double br = 0.0;
    std::vector<std::pair<int, double>> qm, lr;
  };
  void pr_split(int x, int y, bool scan, MlSplit& s) const {
    s.on = false;
    s.qm.clear();
    s.lr.clear();
    if (!pr_ml || (!scan && !pr_partner_ml) || (scan && rec_noml) || !E.end_ok(x, y)) return;
    const StemEst& st = *stem_;
    const int S = pr_ml_span;
    s.on = true;
    s.br = E.w_bp + E.at(x, y);
    if (x >= 1) {
      if (scan) {
        const OMap& at = qmp_at(x - 1);
        for (size_t r = 0; r < at.slots(); ++r) {
          const int i = at.keys[r];
          if (i >= 1 && x - i <= S) s.qm.emplace_back(i - 1, at.vals[r]);
        }
        std::sort(s.qm.begin(), s.qm.end());
      } else {
        for (int i = std::max(1, x - S); i <= x - 1; ++i) {
          const double mb = st.PrMbar(i, x - 1);
          if (mb > HALF_NEG) s.qm.emplace_back(i - 1, mb);
        }
      }
    }
    for (int xo = ml_run_lo(x); xo < x; ++xo) {
      double rl = E.w_cp * (x - xo - 1);
      if (dz) rl += E.run(xo + 1, x - 1);
      s.lr.emplace_back(xo, rl);
    }
  }

  // Helper profile HM (kind 0), HF with rf = S (1) or rf = cap (2) of (xo, y); cached by key3.
  std::vector<std::pair<const PProf*, double>> prHT;
  template <class Get>
  const PProf& pr_helper(std::unordered_map<uint64_t, PProf>& cache, int xo, int y, int kind,
                         bool scan, Get child) {
    const uint64_t key = key3(xo, y, kind);
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;
    std::unique_ptr<SecTimer> ptm(lpf && !scan ? new SecTimer(vi_help_s) : nullptr);
    const long long a0 = pr_adds;
    const StemEst& st = *stem_;
    const std::vector<std::vector<int>>& KR = pr_pair_list(scan);
    const std::vector<int>& kr = KR[xo];
    const int S = pr_ml_span;
    auto close = [&](int a, int b) { return (E.w_ap + E.w_bp + E.at(a, b)) + pair_bonus; };
    prHT.clear();
    if (kind != 0) {
      const int rf = kind == 1 ? S : cap;
      for (auto jt = std::upper_bound(kr.begin(), kr.end(), y); jt != kr.end() && *jt - y - 1 <= rf; ++jt) {
        const int yo = *jt;
        if (!E.end_ok(xo, yo)) continue;
        double w = E.w_cp * (yo - y - 1);
        if (dz) w += E.run(y + 1, yo - 1);
        prHT.emplace_back(child(xo, yo), w + close(xo, yo));
      }
    }
    for (auto jt = std::upper_bound(kr.begin(), kr.end(), y + 1); jt != kr.end() && *jt - y - 1 <= S; ++jt) {
      const int yo = *jt;
      const double mb = st.PrMbar(y + 1, yo - 1);
      if (mb < HALF_NEG || !E.end_ok(xo, yo)) continue;
      prHT.emplace_back(child(xo, yo), mb + close(xo, yo));
    }
    double ref = NEG_INF;
    for (const auto& [q, w] : prHT)
      if (q->mx > HALF_NEG) ref = std::max(ref, w + q->mx);
    PProf p;
    if (ref > HALF_NEG) {
      dacc_begin(ref);
      for (const auto& [q, w] : prHT) dacc_add(*q, w);
      dacc_finish(p, 0.0);
      pr_linify(p);
    }
    const int sc = scan ? 0 : 1;
    if (!scan) vi_help_adds += pr_adds - a0;
    ++pr_hbuilt[sc];
    pr_hentries[sc] += static_cast<long long>(p.k.size());
    pr_hterms[sc] += static_cast<long long>(prHT.size());
    if (scan) {
      gm_help += pbytes(p);
      gm_help_entries_peak = std::max(gm_help_entries_peak, pr_hentries[0] - pr_hg_ent);
      gm_update();
    }
    return cache.emplace(key, std::move(p)).first->second;
  }

  // log of the total of pr_helper(xo, y, kind) from its terms and the children's totals (the
  // helper is uncut, so this is its exact mass); cached for the current y.
  std::unordered_map<uint64_t, double> prHM;
  int prHMy = -1;
  template <class Get>
  double pr_hmass(int xo, int y, int kind, Get child) {
    if (y != prHMy) {
      prHM.clear();
      prHMy = y;
    }
    const uint64_t key = key3(xo, y, kind);
    auto it = prHM.find(key);
    if (it != prHM.end()) return it->second;
    const StemEst& st = *stem_;
    const std::vector<int>& kr = (st.pr_partner ? st.prKeptRightP : st.prKeptRight)[xo];
    const int S = pr_ml_span;
    auto close = [&](int a, int b) { return (E.w_ap + E.w_bp + E.at(a, b)) + pair_bonus; };
    double m = NEG_INF;
    if (kind != 0) {
      const int rf = kind == 1 ? S : cap;
      for (auto jt = std::upper_bound(kr.begin(), kr.end(), y); jt != kr.end() && *jt - y - 1 <= rf; ++jt) {
        const int yo = *jt;
        if (!E.end_ok(xo, yo)) continue;
        double w = E.w_cp * (yo - y - 1);
        if (dz) w += E.run(y + 1, yo - 1);
        const double t = child(xo, yo)->tot;
        if (t > HALF_NEG) m = log_add(m, w + close(xo, yo) + t);
      }
    }
    for (auto jt = std::upper_bound(kr.begin(), kr.end(), y + 1); jt != kr.end() && *jt - y - 1 <= S; ++jt) {
      const int yo = *jt;
      const double mb = st.PrMbar(y + 1, yo - 1);
      if (mb < HALF_NEG || !E.end_ok(xo, yo)) continue;
      const double t = child(xo, yo)->tot;
      if (t > HALF_NEG) m = log_add(m, mb + close(xo, yo) + t);
    }
    return prHM.emplace(key, m).first->second;
  }

  // --pr-fast: starts the dense accumulator on the node x.y with steps nb (pr_steps(x, y, inner,
  // scan)) and adds the steps' w C(x', y') plus, if self_w is finite, self_w at key2(x, y).
  // pr_fast 2 adds the multiloop steps through the helpers in hcache.  The caller finishes.
  template <class Get>
  void pr_node_acc(int x, int y, bool inner, bool scan, const std::vector<ChainNbr>& nb,
                   double self_w, Get child, std::unordered_map<uint64_t, PProf>& hcache) {
    std::vector<std::pair<const PProf*, double>> terms;
    terms.reserve(nb.size());
    const bool vp = lpf && !scan;
    std::chrono::steady_clock::time_point t0;
    if (vp) t0 = std::chrono::steady_clock::now();
    size_t n_int = 0;
    if (pr_fast >= 2) {
      for (const ChainNbr& t : nb)
        if (t.t == 0) terms.emplace_back(child(t.x, t.y), t.w);
      n_int = terms.size();
      MlSplit s;
      pr_split(x, y, scan, s);
      if (s.on) {
        if (vp) {
          const long long c = static_cast<long long>(s.qm.size() + s.lr.size());
          vi_mlc += c;
          vi_mlc_max = std::max(vi_mlc_max, c);
          ++vi_ml_nodes;
        }
        const size_t nc = s.qm.size() + s.lr.size();
        if (!scan && pr_vmlk > 0 && nc > static_cast<size_t>(pr_vmlk)) {
          std::chrono::steady_clock::time_point tk = std::chrono::steady_clock::now();
          std::vector<double> w(nc), ms(nc);
          for (size_t c = 0; c < nc; ++c) {
            const bool q = c < s.qm.size();
            const int xo = q ? s.qm[c].first : s.lr[c - s.qm.size()].first;
            w[c] = (q ? s.qm[c].second : s.lr[c - s.qm.size()].second) + s.br;
            ms[c] = w[c] + pr_hmass(xo, y, q ? (inner ? 2 : 1) : 0, child);
          }
          std::vector<std::pair<double, int>> rk(nc);
          for (size_t c = 0; c < nc; ++c) rk[c] = {ms[c], static_cast<int>(c)};
          std::nth_element(rk.begin(), rk.begin() + pr_vmlk, rk.end(),
                           [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                             return a.first > b.first || (a.first == b.first && a.second < b.second);
                           });
          std::vector<char> keep(nc, 0);
          double mk = NEG_INF, md = NEG_INF;
          for (size_t i = 0; i < nc; ++i) {
            if (i < static_cast<size_t>(pr_vmlk)) {
              keep[rk[i].second] = 1;
              mk = log_add(mk, rk[i].first);
            } else {
              md = log_add(md, rk[i].first);
            }
          }
          const double sc = mk > HALF_NEG && md > HALF_NEG ? log_add(mk, md) - mk : 0.0;
          ++vk_nodes;
          vk_cand += static_cast<long long>(nc);
          vk_drop += static_cast<long long>(nc) - pr_vmlk;
          if (mk > HALF_NEG && md > HALF_NEG) {
            const double r = std::exp(md - log_add(mk, md));
            vk_rmax = std::max(vk_rmax, r);
            vk_rsum += r;
          }
          vk_mass_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - tk).count();
          for (size_t c = 0; c < nc; ++c) {
            if (!keep[c]) continue;
            const bool q = c < s.qm.size();
            const int xo = q ? s.qm[c].first : s.lr[c - s.qm.size()].first;
            terms.emplace_back(&pr_helper(hcache, xo, y, q ? (inner ? 2 : 1) : 0, scan, child), w[c] + sc);
          }
        } else {
          for (const auto& [xo, qv] : s.qm)
            terms.emplace_back(&pr_helper(hcache, xo, y, inner ? 2 : 1, scan, child), qv + s.br);
          for (const auto& [xo, rl] : s.lr)
            terms.emplace_back(&pr_helper(hcache, xo, y, 0, scan, child), rl + s.br);
        }
      }
    } else {
      for (const ChainNbr& t : nb) terms.emplace_back(child(t.x, t.y), t.w);
      n_int = terms.size();
    }
    if (vp) vi_lap(t0, vi_terms_s);
    double ref = self_w;
    for (const auto& [q, w] : terms)
      if (q->mx > HALF_NEG) ref = std::max(ref, w + q->mx);
    dacc_begin(ref > HALF_NEG ? ref : 0.0);
    if (self_w > HALF_NEG) dacc_add1(!scan && pr_vxs ? key2(x, 0) : key2(x, y), self_w);
    const long long a0 = pr_adds;
    for (size_t i = 0; i < n_int; ++i) dacc_add(*terms[i].first, terms[i].second);
    if (vp) {
      vi_add_int += pr_adds - a0;
      vi_lap(t0, vi_add_int_s);
    }
    const long long a1 = pr_adds;
    for (size_t i = n_int; i < terms.size(); ++i) dacc_add(*terms[i].first, terms[i].second);
    if (vp) {
      vi_add_ml += pr_adds - a1;
      vi_lap(t0, vi_add_ml_s);
    }
    pr_nterms[scan ? 0 : 1] += static_cast<long long>(terms.size());
  }

  // Future stems (nested left), built once: V(x, y) = [end_ok] at(x, y) at key (x, y) plus the
  // steps' w V(x', y'), cut, outermost pairs first.  Inner pair x'.y' (end_ok): its steps
  // with the seed cap plus (single_stems) the pair alone, times its own beta2 + at +
  // pair_bonus, gives the arm entries (xs, x'; y', ys, w).  Arms (xs, x') by id in (xs, x')
  // order with their entries by (y', ys) and, per y', log sum_ys w e^beta1 Sbar[ys+1];
  // prByXp: arm ids in (x', xs) order.
  int prA = 0;
  std::vector<int> prArmXs, prArmXp, prArmY0, prYy, prArmE0, prEy, prEys, prXsFirst, prByXp,
      prXpFirst;
  std::vector<double> prYw, prEw;

  struct AT {
    int xs, xp, yp, ys;
    double w;
  };
  static bool at_less(const AT& a, const AT& b) {
    if (a.xs != b.xs) return a.xs < b.xs;
    if (a.xp != b.xp) return a.xp < b.xp;
    if (a.yp != b.yp) return a.yp < b.yp;
    return a.ys < b.ys;
  }

  // --pr-partner-maxpairs K (0 = off): future partner stems of at most K chain pairs, counted
  // innermost to outermost (the inner pair included; branch helices of their multiloops not
  // counted).  V_c(x, y): the profile of exactly c pairs from x.y out (V_1 = the stop at x.y,
  // V_c = sum of the steps' w V_{c-1}(x', y')), keys (xs, ys, c), c <= K - 1; the arm entries
  // of an inner pair are its own stop plus its steps' w sum_c V_c.  The cut acts on sum_c.
  // check: the uncapped V built as well, arm entries compared.
 public:
  int pr_pmax = 0;
  bool pr_pcheck = false;
  long long pr_V_pairs = 0, pr_V_layer_entries = 0, pc_n = 0, pc_key_mis = 0, pe_n = 0,
            pe_key_mis = 0;
  double pc_maxd = 0.0, pe_maxd = 0.0;
  void set_stem_pr_partner_cap(int K, bool check) {
    pr_pmax = K;
    pr_pcheck = check;
  }
  // --pr-partner-ml off: the future partner stems (V and the arm entries, pr_steps with
  // scan = false) take interior steps only; the growth chains (scan = true) keep their
  // multiloop steps.  check (--pr-partner-ml-check, n <= 40, cut 0, v-eps 0): the arm entries
  // against every multiloop-free partner stem enumerated over all enclosing kept pairs.
  bool pr_partner_ml = true, pr_pmlcheck = false;
  long long pm_n = 0, pm_key_mis = 0;
  double pm_maxd = 0.0;
  void set_stem_pr_partner_ml(bool ml, bool check) {
    pr_partner_ml = ml;
    pr_pmlcheck = check;
  }

 private:
  static uint64_t kcl(int xs, int ys, int c) {
    return (static_cast<uint64_t>(static_cast<uint32_t>(xs)) << 40) |
           (static_cast<uint64_t>(static_cast<uint32_t>(ys)) << 8) | static_cast<uint64_t>(c);
  }
  static int kcl_xs(uint64_t k) { return static_cast<int>(k >> 40); }
  static int kcl_ys(uint64_t k) { return static_cast<int>((k >> 8) & 0xffffffffu); }

  // Layered profile of a future stem with left end x: entries sorted, cut by the per-(xs, ys)
  // sum over c (pr_finish_fut's rule); returns the dropped (xs, ys); groups: kept (xs, ys).
  long long pr_finish_layers(PAcc& acc, PProf& p, int x, long long& groups) const {
    std::vector<std::pair<uint64_t, double>>& v = acc.v;
    std::sort(v.begin(), v.end(), [](const std::pair<uint64_t, double>& a,
                                     const std::pair<uint64_t, double>& b) {
      return a.first < b.first;
    });
    const StemEst& st = *stem_;
    auto adj = [&](uint64_t k) {
      return pr_gapfree ? st.prS[kcl_ys(k) + 1] - E.wu * (x - kcl_xs(k)) : 0.0;
    };
    std::vector<size_t> g0;
    std::vector<double> col;
    for (size_t i = 0; i < v.size(); ++i) {
      if (i == 0 || (v[i].first >> 8) != (v[i - 1].first >> 8)) {
        g0.push_back(i);
        col.push_back(v[i].second);
      } else {
        col.back() = log_add(col.back(), v[i].second);
      }
    }
    g0.push_back(v.size());
    double thr = NEG_INF;
    if (pr_cut > 0.0 && !col.empty()) {
      double mx = NEG_INF;
      for (size_t g = 0; g < col.size(); ++g) mx = std::max(mx, col[g] + adj(v[g0[g]].first));
      thr = mx + std::log(pr_cut);
    }
    long long d = 0;
    for (size_t g = 0; g < col.size(); ++g) {
      if (col[g] + adj(v[g0[g]].first) >= thr) {
        ++groups;
        for (size_t i = g0[g]; i < g0[g + 1]; ++i) {
          p.k.push_back(v[i].first);
          p.w.push_back(v[i].second);
        }
      } else {
        ++d;
      }
    }
    return d;
  }

  // --pr-v-eps: keep the entries xs.ys (layered: xs.ys.c) of the map of a pair with left end x
  // by V Sbar[ys+1] / Qzbar[xs][x], largest first, until 1 - eps of their total (key order kept).
  // --pr-v-trim suffix: the same weights, but xs.ys is dropped only if below eps times the sum
  // over the entries with xs' >= xs (the lookups take xs > a2, so only those always come with
  // it); then at most pr_vcap entries, by weight over that suffix sum.
  void pr_vtrim_suffix(PProf& p, int x, bool layered) {
    const StemEst& st = *stem_;
    const size_t m = p.k.size();
    std::vector<double> io(m);
    std::vector<int> xsv(m);
    for (size_t i = 0; i < m; ++i) {
      const int xs = layered ? kcl_xs(p.k[i]) : k_hi(p.k[i]);
      const int ys = layered ? kcl_ys(p.k[i]) : k_lo(p.k[i]);
      const double z = st.Qzbar(xs, x);
      io[i] = p.w[i] + (pr_vxs ? 0.0 : st.prS[ys + 1]) - (z > HALF_NEG ? z : 0.0);
      xsv[i] = xs;
    }
    std::vector<size_t> ord(m);
    for (size_t i = 0; i < m; ++i) ord[i] = i;
    std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return xsv[a] > xsv[b]; });
    std::vector<double> rel(m, NEG_INF);
    const double le = pr_veps > 0.0 ? std::log(pr_veps) : NEG_INF;
    double suf = NEG_INF;
    for (size_t a = 0; a < m;) {
      size_t b = a;
      while (b < m && xsv[ord[b]] == xsv[ord[a]]) suf = log_add(suf, io[ord[b++]]);
      for (size_t c = a; c < b; ++c) rel[ord[c]] = io[ord[c]] - suf;
      a = b;
    }
    std::vector<char> keep(m, 0);
    std::vector<size_t> kept;
    for (size_t i = 0; i < m; ++i)
      if (io[i] < HALF_NEG ? false : rel[i] >= le) kept.push_back(i);
    if (pr_vcap > 0 && kept.size() > static_cast<size_t>(pr_vcap)) {
      v_capped += static_cast<long long>(kept.size()) - pr_vcap;
      std::stable_sort(kept.begin(), kept.end(), [&](size_t a, size_t b) { return rel[a] > rel[b]; });
      double ai = NEG_INF, di = NEG_INF, aw = NEG_INF, dw = NEG_INF;
      for (size_t c = 0; c < kept.size(); ++c) {
        ai = log_add(ai, io[kept[c]]);
        aw = log_add(aw, p.w[kept[c]]);
        if (c >= static_cast<size_t>(pr_vcap)) {
          di = log_add(di, io[kept[c]]);
          dw = log_add(dw, p.w[kept[c]]);
        }
      }
      ++v_cap_pairs;
      const double fi = std::exp(di - ai), fw = std::exp(dw - aw);
      v_cap_fsum += fi;
      v_cap_fmax = std::max(v_cap_fmax, fi);
      v_cap_fwmax = std::max(v_cap_fwmax, fw);
      v_cap_all = log_add(v_cap_all, aw);
      v_cap_drop = log_add(v_cap_drop, dw);
      kept.resize(static_cast<size_t>(pr_vcap));
    }
    for (size_t i : kept) keep[i] = 1;
    size_t o = 0;
    for (size_t i = 0; i < m; ++i)
      if (keep[i]) {
        p.k[o] = p.k[i];
        p.w[o] = p.w[i];
        ++o;
      }
    p.k.resize(o);
    p.w.resize(o);
  }
  // --v-cap n_V: keep the top n_V left ends b of one V map, ranked by C_pre(b-1) * V.
  void pr_vtrim_cpre(PProf& p, bool layered) {
    if (p2_v_cap <= 0 || p.k.empty() || !stem_ || stem_->prC.empty()) return;
    const size_t m = p.k.size();
    std::vector<int> xs(m);
    for (size_t i = 0; i < m; ++i) xs[i] = layered ? kcl_xs(p.k[i]) : k_hi(p.k[i]);
    std::vector<size_t> ord(m);
    for (size_t i = 0; i < m; ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return xs[a] < xs[b]; });
    std::vector<std::pair<int, double>> grp;
    for (size_t a = 0; a < m;) {
      size_t b = a;
      double s = NEG_INF;
      const int x0 = xs[ord[a]];
      while (b < m && xs[ord[b]] == x0) s = log_add(s, p.w[ord[b++]]);
      const double c = x0 > 0 ? stem_->prC[static_cast<size_t>(x0 - 1)] : 0.0;
      grp.emplace_back(x0, s + c);
      a = b;
    }
    if (static_cast<int>(grp.size()) <= p2_v_cap) return;
    std::nth_element(grp.begin(), grp.begin() + p2_v_cap, grp.end(),
                     [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                       return a.second > b.second;
                     });
    std::vector<char> keepx(static_cast<size_t>(n), 0);
    for (int i = 0; i < p2_v_cap; ++i)
      if (grp[static_cast<size_t>(i)].first >= 0 && grp[static_cast<size_t>(i)].first < n)
        keepx[static_cast<size_t>(grp[static_cast<size_t>(i)].first)] = 1;
    size_t o = 0;
    for (size_t i = 0; i < m; ++i)
      if (xs[i] >= 0 && xs[i] < n && keepx[static_cast<size_t>(xs[i])]) {
        p.k[o] = p.k[i];
        p.w[o] = p.w[i];
        ++o;
      }
    p.k.resize(o);
    p.w.resize(o);
  }
  void pr_vtrim(PProf& p, int x, bool layered = false) {
    v_before += static_cast<long long>(p.k.size());
    if (!(p2_v_cap > 0 && p2_v_trim == 2)) {
    if (pr_vsuffix && p.k.size() > 1) {
      pr_vtrim_suffix(p, x, layered);
    } else if (pr_veps > 0.0 && p.k.size() > 1) {
      const StemEst& st = *stem_;
      const size_t m = p.k.size();
      std::vector<double> io(m);
      double tot = NEG_INF;
      for (size_t i = 0; i < m; ++i) {
        const int xs = layered ? kcl_xs(p.k[i]) : k_hi(p.k[i]);
        const int ys = layered ? kcl_ys(p.k[i]) : k_lo(p.k[i]);
        const double z = st.Qzbar(xs, x);
        io[i] = p.w[i] + (pr_vxs ? 0.0 : st.prS[ys + 1]) - (z > HALF_NEG ? z : 0.0);
        tot = log_add(tot, io[i]);
      }
      if (tot > HALF_NEG) {
        std::vector<size_t> ord(m);
        for (size_t i = 0; i < m; ++i) ord[i] = i;
        std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return io[a] > io[b]; });
        const double goal = tot + std::log1p(-pr_veps);
        std::vector<char> keep(m, 0);
        double s = NEG_INF;
        for (size_t i : ord) {
          keep[i] = 1;
          s = log_add(s, io[i]);
          if (s >= goal) break;
        }
        size_t o = 0;
        for (size_t i = 0; i < m; ++i)
          if (keep[i]) {
            p.k[o] = p.k[i];
            p.w[o] = p.w[i];
            ++o;
          }
        p.k.resize(o);
        p.w.resize(o);
      }
    }
    }
    if (p2_v_cap > 0 && p2_v_trim != 1) pr_vtrim_cpre(p, layered);
    v_after += static_cast<long long>(p.k.size());
  }

  // --pr-v-diag / --pr-v-ratio: Vr(xs; x', y') = sum over partner pairs xs.ys, ys > y', of
  // Qb(xs, ys) / Qb(x', y') with the entry's end factors (base(x', y') at(xs, ys) e^beta1
  // Sbar[ys+1]); xs = x': the pair alone.  prVsuf[xs][k]: log suffix sum over the k-th and
  // later partners ys of xs of Qb(xs, ys) at(xs, ys) e^beta1 Sbar[ys+1] (end_ok).
  bool pr_vratio = false;
  std::string pr_vdiag;
  // --pr-v2-dump PATH: per full-support group (xs; x', y') log Vr and log Vr2, where Vr2 takes
  // per partner ys min(Qb(xs,ys)/Qb(x',y'), Qout(x',y')/Qout(xs,ys)) = the inside ratio times
  // min(1, P(x'.y') / P(xs.ys)); then exits.
  std::string pr_v2dump;
  void pr_v2_dump(const std::vector<std::pair<int, int>>& pairs) {
    const StemEst& st = *stem_;
    const auto& KR = pr_KRp();
    FILE* f = std::fopen(pr_v2dump.c_str(), "w");
    if (!f) throw std::runtime_error("cannot open --pr-v2-dump file");
    std::fprintf(f, "xs\txp\typ\tlogVr\tlogVr2\n");
    long long nin = 0, nout = 0, groups = 0;
    for (const auto& [x, y] : pairs) {
      if (!E.end_ok(x, y)) continue;
      const double qi = st.QbIn(x, y);
      const double lpi = st.prLP[static_cast<size_t>(x) * n + y];
      if (qi < HALF_NEG || lpi < HALF_NEG) continue;
      const double base = (E.at(x, y) + E.w_stem) + pair_bonus;
      for (int xs = 0; xs < x; ++xs) {
        double a = NEG_INF, b = NEG_INF;
        for (auto it = std::upper_bound(KR[xs].begin(), KR[xs].end(), y); it != KR[xs].end(); ++it) {
          const int ys = *it;
          const double qb = st.QbIn(xs, ys);
          if (qb < HALF_NEG || !E.end_ok(xs, ys)) continue;
          const double t = qb + E.at(xs, ys) + (E.w_beta1 + st.prS[ys + 1]);
          const double g = std::min(0.0, lpi - st.prLP[static_cast<size_t>(xs) * n + ys]);
          ++(g < 0.0 ? nout : nin);
          a = log_add(a, t);
          b = log_add(b, t + g);
        }
        if (a < HALF_NEG) continue;
        ++groups;
        std::fprintf(f, "%d\t%d\t%d\t%.9g\t%.9g\n", xs, x, y, base - qi + a, base - qi + b);
      }
    }
    std::fclose(f);
    std::fprintf(stderr, "[v2] groups %lld terms_inside %lld terms_outside %lld\n", groups, nin, nout);
  }
  double pr_vsample_p = 0.0;
  long long pr_vsample_max = 10000;
  std::vector<std::vector<double>> prVsuf, prVsufO;
  std::vector<double> prVrg, prVint, prVlen, prVcon, prVshare;
  // --pr-v-diag, per group: Vr2 (per partner min of the inside and outside ratios), Vsep
  // (min of the two sums), Vo (the outside-ratio sum), LB (sum of the per-partner lower bounds
  // Qb(xs,ys)/Qb(x',y') max(0, P(x'.y') + P(xs.ys) - 1) / P(xs.ys)).
  std::vector<double> prVr2g, prVsepg, prVog, prVlbg;
  std::vector<long long> prVuse;
  std::unordered_map<uint64_t, int> prArmOf;
  FILE* pr_vsel = nullptr;
  // --pr-v-sample: PREFIX.arms.tsv.gz, per sampled lookup and arm (xs, x') of the gap, the
  // partner value with the gap fillers under exact V, Vr, Vr2, Vsep (log; -inf: none)
  FILE* pr_varms = nullptr;
  long long pr_vsel_n = 0, pr_vT_n = 0, pr_vr_groups = 0, pr_vr_extra = 0;
  double pr_vr_sec = 0.0, pr_vint_sec = 0.0;
  uint64_t pr_vrng = 0x9E3779B97F4A7C15ull;
  double pr_vrand() {
    pr_vrng ^= pr_vrng << 13;
    pr_vrng ^= pr_vrng >> 7;
    pr_vrng ^= pr_vrng << 17;
    return static_cast<double>(pr_vrng >> 11) * (1.0 / 9007199254740992.0);
  }
  const std::vector<std::vector<int>>& pr_KRp() const {
    const StemEst& st = *stem_;
    return st.pr_partner ? st.prKeptRightP : st.prKeptRight;
  }
  void pr_vsuf_build() {
    if (!prVsuf.empty()) return;
    const StemEst& st = *stem_;
    const auto& KR = pr_KRp();
    prVsuf.assign(n, {});
    prVsufO.assign(n, {});
    for (int xs = 0; xs < n; ++xs) {
      const std::vector<int>& kr = KR[xs];
      std::vector<double>& s = prVsuf[xs];
      std::vector<double>& so = prVsufO[xs];
      s.assign(kr.size() + 1, NEG_INF);
      so.assign(kr.size() + 1, NEG_INF);
      for (int k = static_cast<int>(kr.size()) - 1; k >= 0; --k) {
        const int ys = kr[k];
        const double qb = st.QbIn(xs, ys);
        double t = NEG_INF;
        if (qb > HALF_NEG && E.end_ok(xs, ys)) t = qb + E.at(xs, ys) + (E.w_beta1 + st.prS[ys + 1]);
        s[k] = log_add(s[k + 1], t);
        const double lps = st.prLP.empty() ? NEG_INF : st.prLP[static_cast<size_t>(xs) * n + ys];
        so[k] = log_add(so[k + 1], t > HALF_NEG && lps > HALF_NEG ? t - lps : NEG_INF);
      }
    }
  }
  // kind 0: Vsep = min(Vr, Vo); 1: Vr2; 2: LB; 3: Vo = P(x'.y') sum_ys t / P(xs.ys) / Qb(x',y').
  // xs = x': the pair alone (exact).  Same units as pr_vr_group.
  double pr_vb_group(int xs, int xp, int yp, int kind) const {
    const StemEst& st = *stem_;
    if (xs == xp || st.prLP.empty()) return pr_vr_group(xs, xp, yp);
    const double base = (E.at(xp, yp) + E.w_stem) + pair_bonus;
    const double qi = st.QbIn(xp, yp);
    const double lpi = st.prLP[static_cast<size_t>(xp) * n + yp];
    if (qi < HALF_NEG || lpi < HALF_NEG) return NEG_INF;
    const std::vector<int>& kr = pr_KRp()[xs];
    const size_t k0 = std::upper_bound(kr.begin(), kr.end(), yp) - kr.begin();
    if (kind == 0 || kind == 3) {
      const double so = prVsufO[xs][k0];
      const double vo = so > HALF_NEG ? base - qi + lpi + so : NEG_INF;
      if (kind == 3) return vo;
      const double s = prVsuf[xs][k0];
      const double vr = s > HALF_NEG ? base - qi + s : NEG_INF;
      return std::min(vr, vo);
    }
    const uint64_t key = key3(xs, xp, yp);
    auto& memo = prVbc[kind - 1];
    if (auto it = memo.find(key); it != memo.end()) return it->second;
    const double pi = std::exp(lpi);
    double acc = NEG_INF;
    for (size_t k = k0; k < kr.size(); ++k) {
      const int ys = kr[k];
      const double qb = st.QbIn(xs, ys);
      const double lps = st.prLP[static_cast<size_t>(xs) * n + ys];
      if (qb < HALF_NEG || lps < HALF_NEG || !E.end_ok(xs, ys)) continue;
      const double t = qb + E.at(xs, ys) + (E.w_beta1 + st.prS[ys + 1]);
      if (kind == 1) {
        acc = log_add(acc, t + std::min(0.0, lpi - lps));
      } else {
        const double j = pi + std::expm1(lps);  // P(x'.y') + P(xs.ys) - 1
        if (j > 0.0) acc = log_add(acc, t + std::log(j) - lps);
      }
    }
    const double out = acc > HALF_NEG ? base - qi + acc : NEG_INF;
    memo.emplace(key, out);
    return out;
  }
  mutable std::unordered_map<uint64_t, double> prVbc[2];
  // log Vr(xs; x', y') (the V group value, prYw's units); x'.y' end_ok.
  double pr_vr_group(int xs, int xp, int yp) const {
    const StemEst& st = *stem_;
    const double base = (E.at(xp, yp) + E.w_stem) + pair_bonus;
    if (xs == xp)
      return single_stems ? base + E.at(xp, yp) + (E.w_beta1 + st.prS[yp + 1]) : NEG_INF;
    const double qi = st.QbIn(xp, yp);
    if (qi < HALF_NEG) return NEG_INF;
    const std::vector<int>& kr = pr_KRp()[xs];
    const size_t k = std::upper_bound(kr.begin(), kr.end(), yp) - kr.begin();
    const double s = prVsuf[xs][k];
    return s > HALF_NEG ? base - qi + s : NEG_INF;
  }

  void pr_static_ats(int K, bool count, const std::vector<std::pair<int, int>>& pairs,
                     std::vector<AT>& ats, std::vector<std::vector<AT>>* bk = nullptr) {
    const StemEst& st = *stem_;
    if (pr_vratio && K == 0) {
      const auto& KR = pr_KRp();
      for (const auto& [x, y] : pairs) {
        if (!E.end_ok(x, y)) continue;
        const double base = (E.at(x, y) + E.w_stem) + pair_bonus;
        if (single_stems) ats.push_back({x, x, y, y, base + E.at(x, y)});
        const double qi = st.QbIn(x, y);
        if (qi < HALF_NEG) continue;
        for (int xs = 0; xs < x; ++xs)
          for (auto it = std::upper_bound(KR[xs].begin(), KR[xs].end(), y); it != KR[xs].end(); ++it) {
            const int ys = *it;
            const double qb = st.QbIn(xs, ys);
            if (qb < HALF_NEG || !E.end_ok(xs, ys)) continue;
            const double g =
                g_v2ratio ? std::min(0.0, st.prLP[static_cast<size_t>(x) * n + y] -
                                              st.prLP[static_cast<size_t>(xs) * n + ys])
                          : 0.0;
            ats.push_back({xs, x, y, ys, base + E.at(xs, ys) + qb - qi + g});
            if (count) ++pr_V_entries;
          }
      }
      return;
    }
    std::unordered_map<uint64_t, PProf> V;
    V.reserve(pairs.size() * 2 + 1);
    std::vector<ChainNbr> nb;
    if (K == 0 && pr_fast > 0) {
      // y descending: every step's x'.y' has y' > y, so it is built first, and the helpers
      // of one y are dropped when y changes.
      std::vector<std::pair<int, int>> ord(pairs);
      std::stable_sort(ord.begin(), ord.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
        return a.second > b.second;
      });
      std::unordered_map<uint64_t, PProf> H;
      auto child = [&](int a, int b) -> const PProf* { return &V.find(key2(a, b))->second; };
      // pass 0: the stop at x.y (end_ok); pass 1: the pair alone (single_stems)
      auto node = [&](int x, int y, bool inner, std::unordered_map<uint64_t, PProf>& hc, PProf& p) {
        std::chrono::steady_clock::time_point t0;
        if (lpf) t0 = std::chrono::steady_clock::now();
        pr_steps(x, y, inner, false, nb, pr_fast < 2 || pr_vmlsteps);
        if (count) count_steps(1, nb);
        if (lpf) {
          ++vi_nodes[inner ? 1 : 0];
          for (const ChainNbr& t : nb) ++(t.t == 0 ? vi_nb_int : vi_nb_ml);
          vi_lap(t0, vi_steps_s);
        }
        double self = !inner ? (E.end_ok(x, y) ? E.at(x, y) : NEG_INF)
                             : (single_stems ? E.at(x, y) : NEG_INF);
        if (pr_vxs && self > HALF_NEG) self += E.w_beta1 + st.prS[y + 1];
        pr_node_acc(x, y, inner, false, nb, self, child, hc);
        if (lpf) t0 = std::chrono::steady_clock::now();
        const long long d = pr_finish_fut_d(p, x);
        if (lpf) vi_ent_fin += static_cast<long long>(p.k.size());
        pr_vtrim(p, x);
        if (lpf) {
          vi_ent_kept += static_cast<long long>(p.k.size());
          vi_lap(t0, vi_fin_s);
        }
        if (count) pr_cut_drop += d;
      };
      auto put_v = [&](int x, int y, PProf& p) {
        if (sw_check) sw_ref.emplace(key2(x, y), p);
        if (count) {
          pr_V_entries += static_cast<long long>(p.k.size());
          ++v_pairs;
        }
        pr_linify(p);
        V.emplace(key2(x, y), std::move(p));
      };
      auto put_ats = [&](int x, int y, const PProf& p) {
        std::unique_ptr<SecTimer> ptm(lpf ? new SecTimer(vi_ats_s) : nullptr);
        if (lpf) vi_ats += static_cast<long long>(p.k.size());
        const double base = (E.at(x, y) + E.w_stem) + pair_bonus;
        for (size_t k = 0; k < p.k.size(); ++k) {
          const AT a{k_hi(p.k[k]), x, y, k_lo(p.k[k]), base + p.w[k]};
          if (bk) (*bk)[a.xs].push_back(a);
          else ats.push_back(a);
        }
      };
      auto ats_sz = [&](long long& e) {
        if (!bk) {
          e = static_cast<long long>(ats.size());
          return msz::vec(ats);
        }
        double b = 0.0;
        e = 0;
        for (const std::vector<AT>& v : *bk) {
          e += static_cast<long long>(v.size());
          b += msz::vec(v);
        }
        return b;
      };
      if (pr_vstream) {
        // V(x', y') is read only by pairs with y' - y - 1 <= R (pr_steps / pr_helper reach)
        const int lc = std::min(pr_loop, cap);
        const int R = pr_ml && pr_partner_ml ? std::max({lc, cap, pr_ml_span}) : lc;
        std::vector<std::vector<int>> byY(n);
        std::unordered_map<uint64_t, PProf> H1;
        long long live = 0;
        int yfree = n - 1, ylast = -1;
        for (const auto& [x, y] : ord) {
          cpu_stop_check(1);
          if (y != ylast) {
            H.clear();
            H1.clear();
            ylast = y;
            for (; yfree > y + 1 + R; --yfree) {
              for (int a : byY[yfree]) {
                auto it = V.find(key2(a, yfree));
                if (g_rfv) rfv_put(it->first, it->second);
                live -= static_cast<long long>(it->second.k.size());
                V.erase(it);
              }
              std::vector<int>().swap(byY[yfree]);
            }
          }
          PProf p;
          node(x, y, false, H, p);
          live += static_cast<long long>(p.k.size());
          v_live_max = std::max(v_live_max, live);
          put_v(x, y, p);
          byY[y].push_back(x);
          if (E.end_ok(x, y)) {
            PProf q;
            node(x, y, true, H1, q);
            put_ats(x, y, q);
          }
        }
        if (g_mem_log) {
          long long ve = 0, he = 0, ae = 0;
          const double vb = pmap_sz(V, ve), hb = pmap_sz(H, he) + pmap_sz(H1, he);
          const double ab = ats_sz(ae);
          mem_log("V_stream", {{"V_prof", ve, vb},
                               {"V_help", he, hb},
                               {"V_ats", ae, ab},
                               {"V_live_max", v_live_max, 0.0},
                               {"V_untrimmed", v_before, 0.0},
                               {"V_trimmed", v_after, 0.0}});
        }
        vquery_dump(V);
        mu_dump_lazy(V);
        if (g_rfv) rfv_keep(V);
        return;
      }
      for (int pass = 0; pass < 2; ++pass) {
        H.clear();
        int ylast = -1;
        for (const auto& [x, y] : ord) {
          cpu_stop_check(1);
          if (pass == 1 && !E.end_ok(x, y)) continue;
          if (y != ylast) {
            H.clear();
            ylast = y;
          }
          PProf p;
          node(x, y, pass == 1, H, p);
          if (pass == 0) put_v(x, y, p);
          else put_ats(x, y, p);
        }
        if (pass == 0) v_live_max = pr_V_entries;
        if (g_mem_log) {
          long long ve = 0, he = 0, ae = 0;
          const double vb = pmap_sz(V, ve), hb = pmap_sz(H, he);
          const double ab = ats_sz(ae);
          mem_log(pass == 0 ? "V_pass0" : "V_pass1",
                  {{"V_prof", ve, vb},
                   {"V_help", he, hb},
                   {"V_ats", ae, ab},
                   {"V_untrimmed", v_before, 0.0},
                   {"V_trimmed", v_after, 0.0}});
        }
      }
      vquery_dump(V);
      if (g_rfv) rfv_keep(V);
      return;
    }
    if (K == 0) {
      auto merge = [&](PAcc& acc) {
        for (const ChainNbr& t : nb) {
          const PProf& q = V.find(key2(t.x, t.y))->second;
          for (size_t k = 0; k < q.k.size(); ++k) acc.add(q.k[k], t.w + q.w[k]);
        }
      };
      for (const auto& [x, y] : pairs) {
        cpu_stop_check(1);
        pr_steps(x, y, false, false, nb);
        if (count) count_steps(1, nb);
        PAcc acc;
        if (E.end_ok(x, y)) acc.add(key2(x, y), E.at(x, y));
        merge(acc);
        PProf p;
        const long long d = pr_finish_fut(acc, p, x);
        pr_vtrim(p, x);
        if (count) {
          pr_cut_drop += d;
          pr_V_entries += static_cast<long long>(p.k.size());
        }
        V.emplace(key2(x, y), std::move(p));
      }
      for (const auto& [x, y] : pairs) {
        cpu_stop_check(1);
        if (!E.end_ok(x, y)) continue;
        pr_steps(x, y, true, false, nb);
        if (count) count_steps(1, nb);
        PAcc acc;
        if (single_stems) acc.add(key2(x, y), E.at(x, y));
        merge(acc);
        PProf p;
        const long long d = pr_finish_fut(acc, p, x);
        pr_vtrim(p, x);
        if (count) pr_cut_drop += d;
        const double base = (E.at(x, y) + E.w_stem) + pair_bonus;
        for (size_t k = 0; k < p.k.size(); ++k)
          ats.push_back({k_hi(p.k[k]), x, y, k_lo(p.k[k]), base + p.w[k]});
      }
      return;
    }
    const uint64_t cmax = static_cast<uint64_t>(K - 1);
    for (const auto& [x, y] : pairs) {
      cpu_stop_check(1);
      pr_steps(x, y, false, false, nb);
      if (count) count_steps(1, nb);
      PAcc acc;
      if (cmax >= 1) {
        if (E.end_ok(x, y)) acc.add(kcl(x, y, 1), E.at(x, y));
        for (const ChainNbr& t : nb) {
          const PProf& q = V.find(key2(t.x, t.y))->second;
          for (size_t k = 0; k < q.k.size(); ++k)
            if ((q.k[k] & 0xffu) < cmax) acc.add(q.k[k] + 1, t.w + q.w[k]);
        }
      }
      PProf p;
      long long groups = 0;
      const long long d = pr_finish_layers(acc, p, x, groups);
      pr_vtrim(p, x, true);
      if (count) {
        pr_cut_drop += d;
        pr_V_entries += groups;
        pr_V_layer_entries += static_cast<long long>(p.k.size());
        ++pr_V_pairs;
      }
      V.emplace(key2(x, y), std::move(p));
    }
    for (const auto& [x, y] : pairs) {
      cpu_stop_check(1);
      if (!E.end_ok(x, y)) continue;
      pr_steps(x, y, true, false, nb);
      if (count) count_steps(1, nb);
      PAcc acc;
      if (single_stems) acc.add(key2(x, y), E.at(x, y));
      for (const ChainNbr& t : nb) {
        const PProf& q = V.find(key2(t.x, t.y))->second;
        for (size_t k = 0; k < q.k.size(); ++k)
          acc.add(key2(kcl_xs(q.k[k]), kcl_ys(q.k[k])), t.w + q.w[k]);
      }
      PProf p;
      const long long d = pr_finish_fut(acc, p, x);
      pr_vtrim(p, x);
      if (count) pr_cut_drop += d;
      const double base = (E.at(x, y) + E.w_stem) + pair_bonus;
      for (size_t k = 0; k < p.k.size(); ++k)
        ats.push_back({k_hi(p.k[k]), x, y, k_lo(p.k[k]), base + p.w[k]});
    }
  }

  // --grow-r-from-v: V kept after the build as flat arrays (xs ascending, float log V) per pair
  // (prVkI: key2(x, y) -> offset, length), filled as the build frees each pair's profile;
  // prVod: V of pairs outside V's pair set (key3(x, y, inner)), built from the kept V of their
  // steps with the growth vectors' stop at x.y; prVodH: their helpers; prVch: the kept V of a
  // node's steps as profiles while it is built; prVkey: per step, the U keys' sums over their
  // closings (index order).  Counters: lookups, lookups with no entry, seed lookups through the
  // arm entries, on-demand nodes and entries, steps onto pairs not kept, keys and entries.
 public:
  struct VSpan {
    const int32_t* xs = nullptr;
    const float* lw = nullptr;
    uint32_t len = 0;
  };
  struct VOwn {
    std::vector<int32_t> xs;
    std::vector<float> lw;
    VSpan span() const { return {xs.data(), lw.data(), static_cast<uint32_t>(xs.size())}; }
  };
  std::unordered_map<uint64_t, std::pair<uint32_t, uint32_t>> prVkI;
  std::vector<int32_t> prVkX;
  std::vector<float> prVkW;
  std::unordered_map<uint64_t, VOwn> prVod;
  std::unordered_map<uint64_t, PProf> prVodH, prVch;
  std::unordered_map<int, std::vector<std::pair<int, double>>> prVkey[2];
  std::vector<double> rfvAcc;
  std::vector<int> rfvT;
  long long rfv_look = 0, rfv_miss = 0, rfv_seed = 0, rfv_od = 0, rfv_od_ent = 0,
            rfv_child_miss = 0, rfv_keys = 0, rfv_key_ent = 0;
  double rfv_sec = 0;
  // v-sweep lookups split by U keys (rfv_key) and R/K terms (rfv_root). key0 is a state
  // scored in a column over its beam whose key is -inf because the V read missed.
  long long rfv_look_u = 0, rfv_miss_u = 0, rfv_look_r = 0, rfv_miss_r = 0;
  long long rfv_key0_u = 0, rfv_part_u = 0, rfv_key0_r = 0;
  long long rfv_look_s = 0, rfv_miss_s = 0, vfb_n = 0;
  // Miss kind: (i) map full (b_V), entry cut; (ii) map exists, not full; (iii) no map.
  enum { ST_P1 = 0, ST_P2, ST_P3, ST_R, ST_S, ST_U, ST_N };
  long long vm_i[ST_N] = {}, vm_ii[ST_N] = {}, vm_iii[ST_N] = {};
  long long km_i[ST_N] = {}, km_ii[ST_N] = {}, km_iii[ST_N] = {};
  long long hm_i[ST_N] = {}, hm_ii[ST_N] = {}, hm_iii[ST_N] = {};
  int v_map_kind(int x, int y, int xs) const {
    if (!stem_) return 3;
    auto it = stem_->prMu.find(key2(x, y));
    if (it == stem_->prMu.end()) return 3;
    auto jt = it->second.find(key2(xs, 0));
    if (jt != it->second.end() && jt->second > HALF_NEG) return 0;
    if (p2_v_cap > 0 && static_cast<int>(it->second.size()) >= p2_v_cap) return 1;
    return 2;
  }
  double v_cut_fb_mu(int x, int y, int xs) const {
    if (!stem_) return NEG_INF;
    auto ft = stem_->prMuFloor.find(key2(x, y));
    if (ft == stem_->prMuFloor.end()) return NEG_INF;
    const double c = (xs > 0 && !stem_->prC.empty()) ? stem_->prC[static_cast<size_t>(xs - 1)] : 0.0;
    if (c < HALF_NEG) return NEG_INF;
    return ft->second - c;
  }
  void note_kind(long long* a, long long* b, long long* c, int st, int k) {
    if (st < 0 || st >= ST_N) return;
    if (k == 1) ++a[st];
    else if (k == 2) ++b[st];
    else if (k == 3) ++c[st];
  }
  int k_map_kind(int r, int qp) const {
    if (r < 0 || qp + 1 > n - 1 || prArmY0.empty()) return 3;
    const int xs = prArmXs[r], xp = prArmXp[r], hi = prArmY0[r + 1];
    const int* yy = prYy.data();
    const int k = static_cast<int>(std::upper_bound(yy + prArmY0[r], yy + hi, qp) - yy);
    if (k >= hi) return 3;
    const int mk = v_map_kind(xp, yy[k], xs);
    return mk == 0 ? 3 : mk;
  }
  int h_map_kind(int a2, int j) const {
    if (a2 < 0 || a2 >= n || !stem_) return 3;
    bool full = false, any = false;
    for (int c1 : partners_right[a2]) {
      if (c1 <= j) continue;
      auto it = stem_->prMu.find(key2(a2, c1));
      if (it == stem_->prMu.end()) continue;
      any = true;
      if (p2_v_cap > 0 && static_cast<int>(it->second.size()) >= p2_v_cap) full = true;
    }
    if (!any) return 3;
    return full ? 1 : 2;
  }
  void rfv_put(uint64_t key, const PProf& p) {
    prVkI[key] = {static_cast<uint32_t>(prVkX.size()), static_cast<uint32_t>(p.k.size())};
    for (size_t i = 0; i < p.k.size(); ++i) {
      prVkX.push_back(k_hi(p.k[i]));
      prVkW.push_back(static_cast<float>(p.w[i]));
    }
  }
  void rfv_keep(std::unordered_map<uint64_t, PProf>& V) {
    for (const auto& kv : V) rfv_put(kv.first, kv.second);
    prVkX.shrink_to_fit();
    prVkW.shrink_to_fit();
  }
  bool rfv_get(int x, int y, VSpan& s) const {
    auto it = prVkI.find(key2(x, y));
    if (it == prVkI.end()) return false;
    s = {prVkX.data() + it->second.first, prVkW.data() + it->second.first, it->second.second};
    return true;
  }
  static double rfv_at(const VSpan& s, int xs) {
    const int32_t* e = s.xs + s.len;
    const int32_t* it = std::lower_bound(s.xs, e, xs);
    return it != e && *it == xs ? static_cast<double>(s.lw[it - s.xs]) : NEG_INF;
  }
  static void rfv_prof(const VSpan& s, PProf& p) {
    p = PProf();
    for (uint32_t i = 0; i < s.len; ++i) p.mx = std::max(p.mx, static_cast<double>(s.lw[i]));
    p.k.resize(s.len);
    p.lin.resize(s.len);
    for (uint32_t i = 0; i < s.len; ++i) {
      p.k[i] = key2(s.xs[i], 0);
      p.lin[i] = std::exp(static_cast<double>(s.lw[i]) - p.mx);
    }
  }
  VSpan rfv_node(int x, int y, bool inner) {
    const uint64_t key = key3(x, y, inner ? 1 : 0);
    auto it = prVod.find(key);
    if (it != prVod.end()) return it->second.span();
    const auto rfv_t0 = std::chrono::steady_clock::now();
    const StemEst& st = *stem_;
    std::vector<ChainNbr> nb;
    pr_steps(x, y, inner, false, nb, pr_fast < 2 || pr_vmlsteps);
    const double self = E.end_ok(x, y) ? E.at(x, y) + E.w_beta1 + st.prS[y + 1] : NEG_INF;
    static const PProf empty;
    auto child = [&](int a, int b) -> const PProf* {
      const uint64_t k = key2(a, b);
      auto jt = prVch.find(k);
      if (jt != prVch.end()) return &jt->second;
      VSpan s;
      if (!rfv_get(a, b, s)) {
        ++rfv_child_miss;
        return &empty;
      }
      PProf& p = prVch[k];
      rfv_prof(s, p);
      return &p;
    };
    if (prVodH.size() > 100000) prVodH.clear();
    pr_node_acc(x, y, inner, false, nb, self, child, prVodH);
    prVch.clear();
    PProf p;
    pr_finish_fut_d(p, x);
    VOwn o;
    o.xs.reserve(p.k.size());
    o.lw.reserve(p.k.size());
    for (size_t i = 0; i < p.k.size(); ++i) {
      o.xs.push_back(k_hi(p.k[i]));
      o.lw.push_back(static_cast<float>(p.w[i]));
    }
    ++rfv_od;
    rfv_od_ent += static_cast<long long>(o.xs.size());
    rfv_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - rfv_t0).count();
    return prVod.emplace(key, std::move(o)).first->second.span();
  }
  // --v-sweep: mu[b] of pair x.y. Missing b: --v-cut-fb uses kappa_min/C_pre(b-1) on a
  // full map; --v-floor uses the smallest kept mu or 0; else -inf.
  void mu_dump_lazy(const std::unordered_map<uint64_t, PProf>& V) const {
    if (g_mu_dump.empty()) return;
    std::ofstream md(g_mu_dump);
    for (const auto& kv : V) {
      const int x = k_hi(kv.first), y = k_lo(kv.first);
      std::unordered_map<int, double> grp;
      for (size_t i = 0; i < kv.second.k.size(); ++i) {
        if (!(kv.second.w[i] > HALF_NEG)) continue;
        const int b = k_hi(kv.second.k[i]);
        auto it = grp.find(b);
        if (it == grp.end()) grp.emplace(b, kv.second.w[i]);
        else it->second = log_add(it->second, kv.second.w[i]);
      }
      for (const auto& g : grp) md << "V " << x << " " << y << " " << g.first << " " << g.second << "\n";
    }
  }
  static const char* trace_edge(int t, int trail) {
    if (t == 0) return "E3 interior";
    if (t == 1) return trail == 0 ? "E6 branch-qm-pair" : "E7 branch-qm-tail";
    if (t == 2) return "E5 branch1";
    return trail == 0 ? "E8 branch-qm2-pair" : "E8 branch-qm2-tail";
  }
  void trace_v(int x, int y, int b, const std::unordered_map<uint64_t, PProf>& V,
               std::ostream& o) const {
    auto at = [&](int px, int py, int pb, double& val) {
      val = NEG_INF;
      auto it = V.find(key2(px, py));
      if (it == V.end()) return;
      for (size_t i = 0; i < it->second.k.size(); ++i)
        if (k_hi(it->second.k[i]) == pb && it->second.w[i] > HALF_NEG)
          val = log_add(val, it->second.w[i]);
    };
    o << "TRACE " << x << " " << y << " " << b << "\n";
    int cx = x, cy = y;
    std::vector<ChainNbr> nb;
    for (int hop = 0; hop < n; ++hop) {
      if (cx == b) {
        o << "  step " << hop << " born " << cx << " " << cy
          << " edge V-seed/E4 exterior lead 0 trail 0 pk none\n";
        return;
      }
      pr_steps(cx, cy, false, false, nb, true);
      int bx = -1, by = -1, bt = -1, bl = 0, br = 0;
      double bw = NEG_INF;
      for (const ChainNbr& t : nb) {
        double cv = NEG_INF;
        if (t.x != b) at(t.x, t.y, b, cv);
        if (t.x != b && !(cv > HALF_NEG)) continue;
        const double s = t.w + (cv > HALF_NEG ? cv : 0.0);
        if (!(s > bw)) continue;
        bw = s;
        bx = t.x;
        by = t.y;
        bt = t.t;
        bl = cx - t.x - 1;
        br = t.y - cy - 1;
      }
      if (bx < 0) {
        o << "  step " << hop << " no-parent " << cx << " " << cy << "\n";
        return;
      }
      const char* side = bt == 0 ? "interior" : bt == 1 ? "branch-left" : bt == 2 ? "branch-right" : "branch-both";
      o << "  step " << hop << " " << cx << " " << cy << " -> " << bx << " " << by << " " << side
        << " l1 " << (bt == 0 ? bl : 0) << " l2 " << (bt == 0 ? br : 0) << " lead "
        << (bt == 2 || bt == 3 ? bl : 0) << " trail " << (bt == 1 || bt == 3 ? br : 0)
        << " edge " << trace_edge(bt, bt == 1 || bt == 3 ? br : 0) << " pk none\n";
      cx = bx;
      cy = by;
    }
  }
  void vquery_dump(const std::unordered_map<uint64_t, PProf>& V) const {
    if (g_vquery_out.empty()) return;
    std::ofstream out(g_vquery_out);
    auto at_b = [&](const PProf& p, int b, double& val, double& mx) {
      val = NEG_INF;
      mx = NEG_INF;
      std::unordered_map<int, double> grp;
      grp.reserve(p.k.size() * 2 + 1);
      for (size_t i = 0; i < p.k.size(); ++i) {
        if (!(p.w[i] > HALF_NEG)) continue;
        const int xs = k_hi(p.k[i]);
        auto it = grp.find(xs);
        if (it == grp.end()) grp.emplace(xs, p.w[i]);
        else it->second = log_add(it->second, p.w[i]);
      }
      for (const auto& kv : grp) {
        if (kv.second > mx) mx = kv.second;
        if (kv.first == b) val = kv.second;
      }
    };
    std::ofstream tr;
    if (!g_vtrace_path.empty()) tr.open(g_vtrace_path);
    for (const VQ& q : g_vquery) {
      if (q.b >= 0) {
        double val = NEG_INF, mx = NEG_INF;
        auto it = V.find(key2(q.x, q.y));
        if (it != V.end()) at_b(it->second, q.b, val, mx);
        out << (q.tag == 1 ? "K " : "V ") << q.st << " " << q.x << " " << q.y << " " << q.b
            << " " << val << " " << mx << "\n";
        if (tr.is_open() && val > HALF_NEG) trace_v(q.x, q.y, q.b, V, tr);
        continue;
      }
      double best = NEG_INF, best_mx = NEG_INF;
      int nz = 0, maps = 0;
      for (const auto& kv : V) {
        if (k_hi(kv.first) != q.x || k_lo(kv.first) <= q.y) continue;
        ++maps;
        double val = NEG_INF, mx = NEG_INF;
        at_b(kv.second, -1, val, mx);
        if (mx > HALF_NEG) {
          ++nz;
          if (mx > best) {
            best = mx;
            best_mx = mx;
          }
        }
      }
      out << "H " << q.st << " " << q.x << " " << q.y << " " << maps << " " << nz << " " << best
          << " " << best_mx << "\n";
    }
  }
  void miss_k(int r, int qp, int st) {
    if (g_miss_cap <= 0 || (st != ST_P3 && st != ST_S)) return;
    if (r < 0 || qp + 1 > n - 1 || prArmY0.empty()) return;
    const int xs = prArmXs[r], xp = prArmXp[r], hi = prArmY0[r + 1];
    const int* yy = prYy.data();
    const int k = static_cast<int>(std::upper_bound(yy + prArmY0[r], yy + hi, qp) - yy);
    if (k >= hi || v_map_kind(xp, yy[k], xs) != 2) return;
    miss_take(1, st, xp, yy[k], xs);
  }
  double mu_at(int x, int y, int xs, int st = -1) {
    const int kind = v_map_kind(x, y, xs);
    if (kind == 2) miss_take(0, st, x, y, xs);
    if (kind == 0) {
      auto it = stem_->prMu.find(key2(x, y));
      return it->second.find(key2(xs, 0))->second;
    }
    note_kind(vm_i, vm_ii, vm_iii, st, kind);
    if (p2_v_cut_fb && kind == 1) {
      const double fb = v_cut_fb_mu(x, y, xs);
      if (fb > HALF_NEG) {
        ++vfb_n;
        return fb;
      }
    }
    if (p2_v_floor == 2) return 0.0;
    if (p2_v_floor == 1) {
      auto ft = stem_->prMuFloor.find(key2(x, y));
      return ft == stem_->prMuFloor.end() ? NEG_INF : ft->second;
    }
    return NEG_INF;
  }
  // V of the root x.y (seed: x.y is the stem's inner pair) at the outer left end xs
  double rfv_root(int x, int y, bool seed, int xs) {
    if (p2_vsweep) {
      ++rfv_look;
      ++rfv_look_r;
      const double v = mu_at(x, y, xs, ST_R);
      if (v < HALF_NEG) {
        ++rfv_miss;
        ++rfv_miss_r;
        ++rfv_key0_r;
      }
      return v;
    }
    ++rfv_look;
    VSpan p;
    const bool kept = rfv_get(x, y, p);
    double v = NEG_INF;
    if (!kept || (seed && !E.end_ok(x, y))) {
      v = rfv_at(rfv_node(x, y, seed), xs);
    } else if (!seed) {
      v = rfv_at(p, xs);
    } else if (xs == x) {
      v = E.at(x, y) + E.w_beta1 + stem_->prS[y + 1];
    } else {
      ++rfv_seed;
      const int r = larm(xs, x);
      if (r >= 0) {
        const auto b = prYy.begin();
        auto it = std::lower_bound(b + prArmY0[r], b + prArmY0[r + 1], y);
        if (it != b + prArmY0[r + 1] && *it == y)
          v = prYw[static_cast<size_t>(it - b)] - ((E.at(x, y) + E.w_stem) + pair_bonus);
      }
    }
    if (v < HALF_NEG) ++rfv_miss;
    return v;
  }
  // U key (ls / rs of p at step q): sum over its closings c of w_c V(c.x, c.y) at xs
  double rfv_key(bool ls, int p, bool seed, int q, int xs) {
    if (p2_vsweep) {
      ++rfv_look;
      ++rfv_look_u;
      std::vector<Closing> cl;
      if (ls) pr_ls_closings(p, seed, q, cl);
      else closings_rs(p, q, cl);
      double s = NEG_INF;
      int missed = 0;
      for (const Closing& c : cl) {
        const double v = mu_at(c.x, c.y, xs, ST_U);
        if (v > HALF_NEG) s = log_add(s, c.w + v);
        else ++missed;
      }
      if (s < HALF_NEG) {
        ++rfv_miss;
        ++rfv_miss_u;
        ++rfv_key0_u;
      } else if (missed > 0) {
        ++rfv_part_u;
      }
      return s;
    }
    auto& m = prVkey[ls ? 0 : 1];
    const int key = ls ? 2 * p + (seed ? 1 : 0) : p;
    auto it = m.find(key);
    if (it == m.end()) {
      std::vector<Closing> cl;
      if (ls) pr_ls_closings(p, seed, q, cl);
      else closings_rs(p, q, cl);
      if (rfvAcc.size() < static_cast<size_t>(n)) rfvAcc.assign(n, NEG_INF);
      rfvT.clear();
      for (const Closing& c : cl) {
        VSpan s;
        if (!rfv_get(c.x, c.y, s)) s = rfv_node(c.x, c.y, false);
        for (uint32_t i = 0; i < s.len; ++i) {
          const int a = s.xs[i];
          if (rfvAcc[a] < HALF_NEG) rfvT.push_back(a);
          rfvAcc[a] = log_add(rfvAcc[a], c.w + static_cast<double>(s.lw[i]));
        }
      }
      std::sort(rfvT.begin(), rfvT.end());
      std::vector<std::pair<int, double>> o;
      o.reserve(rfvT.size());
      for (int a : rfvT) {
        o.emplace_back(a, rfvAcc[a]);
        rfvAcc[a] = NEG_INF;
      }
      ++rfv_keys;
      rfv_key_ent += static_cast<long long>(o.size());
      it = m.emplace(key, std::move(o)).first;
    }
    ++rfv_look;
    const auto& v = it->second;
    auto jt = std::lower_bound(v.begin(), v.end(), xs,
                               [](const std::pair<int, double>& e, int x) { return e.first < x; });
    if (jt != v.end() && jt->first == xs) return jt->second;
    ++rfv_miss;
    return NEG_INF;
  }
  double rfv_bytes(long long& e) const {
    e = static_cast<long long>(prVkX.size());
    double b = msz::umap(prVkI) + msz::vec(prVkX) + msz::vec(prVkW) + msz::umap(prVod);
    for (const auto& kv : prVod) {
      e += static_cast<long long>(kv.second.xs.size());
      b += msz::vec(kv.second.xs) + msz::vec(kv.second.lw);
    }
    return b;
  }

 private:
  // Outside sweep of the pass-0 future-stem profile. Multiloop steps of pr_steps are not
  // enumerated: the right side of a closing pair is a per-left-end chain (empty gap, one
  // base, then a long gap whose dangle factors), and the left side is a kept Qm. Weights
  // are the ones pr_steps uses (PrMbar, w_ap/w_bp/w_cp, at, run, pair_bonus).
  void sw_verify() {
    if (!pr_vxs) {
      std::fprintf(stderr, "[sweep-check] needs --pr-v-key xs\n");
      return;
    }
    if (sw_ref.empty()) {
      std::fprintf(stderr, "[sweep-check] no pass-0 profiles saved\n");
      return;
    }
    using Mu = std::unordered_map<int, double>;
    auto add_mu = [](Mu& d, const Mu& s, double a) {
      if (a < HALF_NEG) return;
      for (const auto& [k, v] : s) {
        if (v < HALF_NEG) continue;
        const double w = a + v;
        auto it = d.find(k);
        if (it == d.end()) d.emplace(k, w);
        else it->second = log_add(it->second, w);
      }
    };
    auto add1 = [](Mu& d, int k, double w) {
      if (w < HALF_NEG) return;
      auto it = d.find(k);
      if (it == d.end()) d.emplace(k, w);
      else it->second = log_add(it->second, w);
    };
    const StemEst& st = *stem_;
    const int S = pr_ml_span;
    const int lc = std::min(pr_loop, cap);
    std::unordered_set<uint64_t> is_pair;
    std::vector<std::vector<int>> xs_at(n);
    std::unordered_map<uint64_t, Mu> pmu;
    is_pair.reserve(sw_pairs.size() * 2 + 1);
    pmu.reserve(sw_pairs.size() * 2 + 1);
    for (const auto& [x, y] : sw_pairs) {
      is_pair.insert(key2(x, y));
      xs_at[y].push_back(x);
      pmu.emplace(key2(x, y), Mu{});
    }
    auto close_of = [&](int xo, int yo) {
      if (!E.end_ok(xo, yo)) return NEG_INF;
      return (E.w_ap + E.w_bp + E.at(xo, yo)) + pair_bonus;
    };
    std::unordered_map<int, Mu> glong;
    for (int y = n - 1; y >= 0; --y) {
      std::unordered_map<int, Mu> gap0, gap1, rq;
      auto fold = [&](int yo, std::unordered_map<int, Mu>& dst, double extra_base) {
        if (yo < 0 || yo >= n) return;
        for (int xo : xs_at[yo]) {
          const double cl = close_of(xo, yo);
          if (cl < HALF_NEG) continue;
          add_mu(dst[xo], pmu[key2(xo, yo)], cl + extra_base);
        }
      };
      fold(y + 1, gap0, 0.0);
      if (y + 2 < n)
        fold(y + 2, gap1, (dz ? E.s1w[y + 1] : 0.0) + E.w_cp);
      if (y + 3 < n) {
        const int yo = y + 3;
        const double extra = (dz ? E.d3w[yo - 1] : 0.0) + E.w_cp * (yo - 1);
        fold(yo, glong, extra);
      }
      for (int yo = y + 2; yo < n; ++yo) {
        if (yo - y - 1 > S) continue;
        const double mb = st.PrMbar(y + 1, yo - 1);
        if (mb < HALF_NEG) continue;
        for (int xo : xs_at[yo]) {
          const double cl = close_of(xo, yo);
          if (cl < HALF_NEG) continue;
          add_mu(rq[xo], pmu[key2(xo, yo)], mb + cl);
        }
      }
      const double ls = (dz && y + 1 < n ? E.d5w[y + 1] : 0.0) - E.w_cp * y;
      for (int x : xs_at[y]) {
        if (!E.end_ok(x, y)) continue;
        const double br = E.w_bp + E.at(x, y);
        if (br < HALF_NEG) continue;
        Mu& dst = pmu[key2(x, y)];
        for (int xo = ml_run_lo(x); xo < x; ++xo) {
          auto it = rq.find(xo);
          if (it == rq.end()) continue;
          double lw = E.w_cp * (x - xo - 1);
          if (dz && xo + 1 <= x - 1) lw += E.run(xo + 1, x - 1);
          add_mu(dst, it->second, br + lw);
        }
        for (int i = std::max(1, x - S); i <= x - 1; ++i) {
          const double qv = st.PrMbar(i, x - 1);
          if (qv < HALF_NEG) continue;
          const int xo = i - 1;
          const double sc = br + qv;
          if (auto it = gap0.find(xo); it != gap0.end()) add_mu(dst, it->second, sc);
          if (auto it = gap1.find(xo); it != gap1.end()) add_mu(dst, it->second, sc);
          if (auto it = glong.find(xo); it != glong.end()) add_mu(dst, it->second, sc + ls);
          if (auto it = rq.find(xo); it != rq.end()) add_mu(dst, it->second, sc);
        }
      }
      for (int x : xs_at[y]) {
        Mu& m = pmu[key2(x, y)];
        if (E.end_ok(x, y)) {
          double self = E.at(x, y);
          if (self > HALF_NEG) add1(m, x, self + E.w_beta1 + st.prS[y + 1]);
        }
        for (int xi = x + 1; xi <= std::min(n - 1, x + 1 + lc); ++xi) {
          for (int yi = y - 1; yi >= std::max(0, y - 1 - lc); --yi) {
            const uint64_t ik = key2(xi, yi);
            if (!is_pair.count(ik)) continue;
            const double s = E.stem_interior(x, xi, yi, y);
            if (s < HALF_NEG) continue;
            add_mu(pmu[ik], m, s + pair_bonus);
          }
        }
      }
    }
    auto rel = [](double a, double b) {
      if (a < HALF_NEG && b < HALF_NEG) return 0.0;
      if (a < HALF_NEG || b < HALF_NEG) return 1e300;
      return std::expm1(std::fabs(a - b));
    };
    double worst = 0.0;
    int wx = -1, wy = -1, wb = -1;
    double ww = 0, wr = 0;
    long long nent = 0, nmiss = 0;
    for (const auto& [x, y] : sw_pairs) {
      const Mu& got = pmu[key2(x, y)];
      auto rit = sw_ref.find(key2(x, y));
      if (rit == sw_ref.end()) {
        ++nmiss;
        continue;
      }
      std::unordered_set<int> seen;
      for (size_t t = 0; t < rit->second.k.size(); ++t) {
        const int b = k_hi(rit->second.k[t]);
        seen.insert(b);
        ++nent;
        auto it = got.find(b);
        const double g = it == got.end() ? NEG_INF : it->second;
        const double r = rel(g, rit->second.w[t]);
        if (r > worst) {
          worst = r;
          wx = x;
          wy = y;
          wb = b;
          ww = g;
          wr = rit->second.w[t];
        }
      }
      for (const auto& [b, g] : got) {
        if (seen.count(b)) continue;
        ++nent;
        const double r = rel(g, NEG_INF);
        if (r > worst) {
          worst = r;
          wx = x;
          wy = y;
          wb = b;
          ww = g;
          wr = NEG_INF;
        }
      }
    }
    std::fprintf(stderr,
                 "[sweep-check] pairs %zu entries %lld missing_pairs %lld max_rel %.3e "
                 "at %d.%d b %d sweep %.9g ref %.9g\n",
                 sw_pairs.size(), nent, nmiss, worst, wx, wy, wb, ww, wr);
  }

  // Gamma half of the same chain. Left Qm is PrMbar or the scan's qmp_at. Self is at(x, y)
  // at key (x, y), with no beta1; the scan prefix is applied when summing onto q'.
  void sw_growth_check() {
    if (g_ntr_K > 0 || g_ml_c > 0) {
      std::fprintf(stderr, "[sweep-g] skipped: needs --grow-node-trunc 0 and --grow-ml-c 0\n");
      return;
    }
    using Ga = std::unordered_map<uint64_t, double>;
    auto add_ga = [](Ga& d, const Ga& s, double a) {
      if (a < HALF_NEG) return;
      for (const auto& [k, v] : s) {
        if (v < HALF_NEG) continue;
        const double w = a + v;
        auto it = d.find(k);
        if (it == d.end()) d.emplace(k, w);
        else it->second = log_add(it->second, w);
      }
    };
    auto add1 = [](Ga& d, uint64_t k, double w) {
      if (w < HALF_NEG) return;
      auto it = d.find(k);
      if (it == d.end()) d.emplace(k, w);
      else it->second = log_add(it->second, w);
    };
    auto add1_q = [](std::unordered_map<int, double>& d, int k, double w) {
      if (w < HALF_NEG) return;
      auto it = d.find(k);
      if (it == d.end()) d.emplace(k, w);
      else it->second = log_add(it->second, w);
    };
    const StemEst& st = *stem_;
    const int S = pr_ml_span;
    const int lc = std::min(pr_loop, cap);
    std::vector<std::pair<int, int>> pairs;
    const std::vector<std::vector<int>>& src = pr_pair_list(true);
    for (int x = 0; x < n; ++x)
      for (int y : src[x]) pairs.emplace_back(x, y);
    auto once = [&](bool scan_qm, const char* tag) {
      std::unordered_set<uint64_t> is_pair;
      std::vector<std::vector<int>> xs_at(n);
      std::unordered_map<uint64_t, Ga> pg;
      is_pair.reserve(pairs.size() * 2 + 1);
      pg.reserve(pairs.size() * 2 + 1);
      for (const auto& [x, y] : pairs) {
        is_pair.insert(key2(x, y));
        xs_at[y].push_back(x);
        pg.emplace(key2(x, y), Ga{});
      }
      auto lq = [&](int i, int j) {
        if (i < 0 || j < i || j >= n) return NEG_INF;
        return scan_qm ? qmp_at(j).get(i) : st.PrMbar(i, j);
      };
      auto close_of = [&](int xo, int yo) {
        if (!E.end_ok(xo, yo)) return NEG_INF;
        return (E.w_ap + E.w_bp + E.at(xo, yo)) + pair_bonus;
      };
      std::unordered_map<int, Ga> glong;
      for (int y = n - 1; y >= 0; --y) {
        std::unordered_map<int, Ga> gap0, gap1, rq;
        auto fold = [&](int yo, std::unordered_map<int, Ga>& dst, double extra) {
          if (yo < 0 || yo >= n) return;
          for (int xo : xs_at[yo]) {
            const double cl = close_of(xo, yo);
            if (cl < HALF_NEG) continue;
            add_ga(dst[xo], pg[key2(xo, yo)], cl + extra);
          }
        };
        fold(y + 1, gap0, 0.0);
        if (y + 2 < n) fold(y + 2, gap1, (dz ? E.s1w[y + 1] : 0.0) + E.w_cp);
        if (y + 3 < n) {
          const int yo = y + 3;
          fold(yo, glong, (dz ? E.d3w[yo - 1] : 0.0) + E.w_cp * (yo - 1));
        }
        for (int yo = y + 2; yo < n; ++yo) {
          if (yo - y - 1 > S) continue;
          const double mb = st.PrMbar(y + 1, yo - 1);
          if (mb < HALF_NEG) continue;
          for (int xo : xs_at[yo]) {
            const double cl = close_of(xo, yo);
            if (cl < HALF_NEG) continue;
            add_ga(rq[xo], pg[key2(xo, yo)], mb + cl);
          }
        }
        const double ls = (dz && y + 1 < n ? E.d5w[y + 1] : 0.0) - E.w_cp * y;
        for (int x : xs_at[y]) {
          if (!E.end_ok(x, y)) continue;
          const double br = E.w_bp + E.at(x, y);
          if (br < HALF_NEG) continue;
          Ga& dst = pg[key2(x, y)];
          for (int xo = ml_run_lo(x); xo < x; ++xo) {
            auto it = rq.find(xo);
            if (it == rq.end()) continue;
            double lw = E.w_cp * (x - xo - 1);
            if (dz && xo + 1 <= x - 1) lw += E.run(xo + 1, x - 1);
            add_ga(dst, it->second, br + lw);
          }
          for (int i = std::max(1, x - S); i <= x - 1; ++i) {
            const double qv = lq(i, x - 1);
            if (qv < HALF_NEG) continue;
            const int xo = i - 1;
            const double sc = br + qv;
            if (auto it = gap0.find(xo); it != gap0.end()) add_ga(dst, it->second, sc);
            if (auto it = gap1.find(xo); it != gap1.end()) add_ga(dst, it->second, sc);
            if (auto it = glong.find(xo); it != glong.end()) add_ga(dst, it->second, sc + ls);
            if (auto it = rq.find(xo); it != rq.end()) add_ga(dst, it->second, sc);
          }
        }
        for (int x : xs_at[y]) {
          Ga& m = pg[key2(x, y)];
          if (E.end_ok(x, y)) add1(m, key2(x, y), E.at(x, y));
          for (int xi = x + 1; xi <= std::min(n - 1, x + 1 + lc); ++xi) {
            for (int yi = y - 1; yi >= std::max(0, y - 1 - lc); --yi) {
              const uint64_t ik = key2(xi, yi);
              if (!is_pair.count(ik)) continue;
              const double s = E.stem_interior(x, xi, yi, y);
              if (s < HALF_NEG) continue;
              add_ga(pg[ik], m, s + pair_bonus);
            }
          }
        }
      }
      auto rel = [](double a, double b) {
        if (a < HALF_NEG && b < HALF_NEG) return 0.0;
        if (a < HALF_NEG || b < HALF_NEG) return 1e300;
        return std::expm1(std::fabs(a - b));
      };
      double worst = 0.0;
      int wx = -1, wy = -1, wq = -1;
      double ww = 0, wr = 0;
      int fx = -1, fy = -1, fq = -1;
      double fw = 0, fr = 0;
      long long nroot = 0, nent = 0;
      for (const auto& [key, rn] : prRR) {
        const int inner = static_cast<int>(key & ((1ull << 21) - 1));
        if (inner != 0) continue;
        const int y = static_cast<int>((key >> 21) & ((1ull << 21) - 1));
        const int x = static_cast<int>(key >> 42);
        ++nroot;
        std::unordered_map<int, double> byq;
        auto git = pg.find(key2(x, y));
        if (git != pg.end()) {
          for (const auto& [k, w] : git->second) {
            const int p = k_hi(k), q = k_lo(k);
            add1_q(byq, q, w + prefix_log(p));
          }
        }
        for (int q = 0; q < n; ++q) {
          double ref = NEG_INF;
          for (int m = 0; m < rnL; ++m) {
            const RV& a = rn.L[m];
            if (a.v.empty() || q < a.lo || q >= a.lo + static_cast<int>(a.v.size())) continue;
            const double v = a.v[static_cast<size_t>(q - a.lo)];
            if (v > 0.0) ref = log_add(ref, a.mx + std::log(v));
          }
          auto it = byq.find(q);
          const double got = it == byq.end() ? NEG_INF : it->second;
          if (got < HALF_NEG && ref < HALF_NEG) continue;
          ++nent;
          const double r = rel(got, ref);
          if (r > worst) {
            worst = r;
            wx = x;
            wy = y;
            wq = q;
            ww = got;
            wr = ref;
          }
          if (r > 1e-9 && (fx < 0 || y < fy || (y == fy && x < fx) || (y == fy && x == fx && q < fq))) {
            fy = y;
            fx = x;
            fq = q;
            fw = got;
            fr = ref;
          }
        }
      }
      std::fprintf(stderr,
                   "[sweep-g] %s roots %lld entries %lld max_rel %.3e at %d.%d q %d sweep %.9g ref %.9g "
                   "first %d.%d q %d sweep %.9g ref %.9g\n",
                   tag, nroot, nent, worst, wx, wy, wq, ww, wr, fx, fy, fq, fw, fr);
      if (!scan_qm && fx >= 0) {
        for (int i = 1; i <= fx - 1; ++i) {
          const double a = st.PrMbar(i, fx - 1), b = qmp_at(fx - 1).get(i);
          const double r = rel(a, b);
          if (r > 1e-9)
            std::fprintf(stderr, "[sweep-g] Qm %d..%d pre %.9g scan %.9g rel %.3e\n", i, fx - 1, a, b, r);
        }
      }
    };
    once(false, "pre-Qm");
    once(true, "scan-Qm");
  }

  // Gamma for --sweep. Right-gap limit rf (full span, or the loop cap for an inner root).
  // Left Qm is the nested pre-pass. Self is at(x, y), prefix not stored.
  void sw_fill(int slot, int rf) {
    using Ga = std::unordered_map<uint64_t, double>;
    auto add_ga = [&](Ga& d, const Ga& s, double a) {
      if (a < HALF_NEG || s.empty()) return;
      const auto tm0 = std::chrono::steady_clock::now();
      for (const auto& [k, v] : s) {
        if (v < HALF_NEG) continue;
        const double w = a + v;
        auto it = d.find(k);
        if (it == d.end()) d.emplace(k, w);
        else it->second = log_add(it->second, w);
      }
      g_merge_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - tm0).count();
      ++g_pushes;
    };
    auto add1 = [](Ga& d, uint64_t k, double w) {
      if (w < HALF_NEG) return;
      auto it = d.find(k);
      if (it == d.end()) d.emplace(k, w);
      else it->second = log_add(it->second, w);
    };
    const StemEst& st = *stem_;
    // Cap a pair's map before it is pushed (--g-cut). See g_cut enum at file scope.
    auto cap_ga = [&](Ga& m) {
      const long long pre = static_cast<long long>(m.size());
      ++g_pre_n;
      g_pre_sum += pre;
      if (pre > g_pre_max) g_pre_max = pre;
      const auto tc0 = std::chrono::steady_clock::now();
      auto finish = [&]() {
        g_cut_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - tc0).count();
        const long long post = static_cast<long long>(m.size());
        ++g_post_n;
        g_post_sum += post;
        if (post > g_post_max) g_post_max = post;
      };
      const int K = g_ntr_K;
      if (K <= 0 || m.empty()) {
        finish();
        return;
      }
      if (g_cut == G_CUT_ENTRIES) {
        // Top-K raw (p',q') by C_pre(p'-1) * G * S_pre(q'+1). Scan still applies Chat at read.
        struct E {
          uint64_t k;
          double w, s;
        };
        std::vector<E> v;
        v.reserve(m.size());
        for (const auto& [k, w] : m) {
          if (w < HALF_NEG) continue;
          const int p = k_hi(k), q = k_lo(k);
          double s = w;
          if (p > 0 && !st.prC.empty()) s += st.prC[static_cast<size_t>(p - 1)];
          if (q + 1 < 0 || q + 1 >= static_cast<int>(st.prS.size())) continue;
          const double sb = st.prS[static_cast<size_t>(q + 1)];
          if (sb < HALF_NEG) continue;
          s += sb;
          v.push_back({k, w, s});
        }
        if (static_cast<int>(v.size()) > K) {
          std::nth_element(v.begin(), v.begin() + K, v.end(), [](const E& a, const E& b) {
            if (a.s != b.s) return a.s > b.s;
            return a.k < b.k;
          });
          v.resize(static_cast<size_t>(K));
        }
        Ga out;
        out.reserve(v.size());
        for (const E& e : v) out.emplace(e.k, e.w);
        m.swap(out);
        finish();
        return;
      }
      // topq and collapse: rank by q' with Ghat[q'] = sum_p' G*C_pre (+ ntrUL / Sbar).
      std::unordered_map<int, double> byq;
      byq.reserve(m.size());
      for (const auto& [k, w] : m) {
        if (w < HALF_NEG) continue;
        const int p = k_hi(k), q = k_lo(k);
        double g = w;
        // Collapse stores key2(0,q') with Chat already applied.
        if (p > 0 && !st.prC.empty()) g += st.prC[static_cast<size_t>(p - 1)];
        auto it = byq.find(q);
        if (it == byq.end()) byq.emplace(q, g);
        else it->second = log_add(it->second, g);
      }
      struct E {
        int q;
        double g, s;
      };
      std::vector<E> v;
      v.reserve(byq.size());
      for (const auto& [q, g] : byq) {
        double s = NEG_INF;
        if (p2_g_suffix) {
          const double sb = st.prS[static_cast<size_t>(q + 1)];
          if (sb > HALF_NEG) s = g + sb;
        } else if (q >= 0 && q < static_cast<int>(ntrUL.size()) && ntrUL[q] > HALF_NEG) {
          s = g + ntrUL[q];
        }
        if (s > HALF_NEG) v.push_back({q, g, s});
      }
      if (static_cast<int>(v.size()) > K) {
        std::nth_element(v.begin(), v.begin() + K, v.end(), [](const E& a, const E& b) {
          if (a.s != b.s) return a.s > b.s;
          return a.q < b.q;
        });
        v.resize(static_cast<size_t>(K));
      }
      if (g_cut == G_CUT_COLLAPSE) {
        Ga out;
        out.reserve(v.size());
        for (const E& e : v) out.emplace(key2(0, e.q), e.g);
        m.swap(out);
        finish();
        return;
      }
      // topq: retain all raw (p',q') mass for the kept q'.
      std::unordered_set<int> keep;
      keep.reserve(v.size());
      for (const E& e : v) keep.insert(e.q);
      Ga out;
      out.reserve(m.size());
      for (const auto& [k, w] : m) {
        if (w < HALF_NEG) continue;
        if (keep.count(k_lo(k))) out.emplace(k, w);
      }
      m.swap(out);
      finish();
    };
    std::vector<char> fin(static_cast<size_t>(n), 0);
    const int S = pr_ml_span;
    const int lc = std::min(pr_loop, cap);
    std::unordered_set<uint64_t> is_pair;
    std::vector<std::vector<int>> xs_at(n);
    std::unordered_map<uint64_t, Ga> pg;
    const std::vector<std::vector<int>>& src = pr_pair_list(true);
    for (int x = 0; x < n; ++x)
      for (int y : src[x]) {
        is_pair.insert(key2(x, y));
        xs_at[y].push_back(x);
        pg.emplace(key2(x, y), Ga{});
      }
    auto close_of = [&](int xo, int yo) {
      if (!E.end_ok(xo, yo)) return NEG_INF;
      return (E.w_ap + E.w_bp + E.at(xo, yo)) + pair_bonus;
    };
    std::unordered_map<int, Ga> long_acc;
    const bool full_rf = rf >= n - 1;
    const bool qm_ix = static_cast<int>(st.prMbarRight.size()) == n;
    ml_yo_S = S;
    ml_yo_qm_ix = qm_ix;
    for (int y = n - 1; y >= 0; --y) {
      std::unordered_map<int, Ga> gap0, gap1, glong_win, rq;
      auto fold = [&](int yo, std::unordered_map<int, Ga>& dst, double extra) {
        if (yo < 0 || yo >= n) return;
        if (!fin[static_cast<size_t>(yo)])
          throw std::runtime_error("sweep: read of a pair column before it was capped");
        for (int xo : xs_at[yo]) {
          const double cl = close_of(xo, yo);
          if (cl < HALF_NEG) continue;
          add_ga(dst[xo], pg[key2(xo, yo)], cl + extra);
        }
      };
      auto add_rq = [&](int yo, double mb) {
        if (yo < 0 || yo >= n || mb < HALF_NEG) return;
        if (!fin[static_cast<size_t>(yo)])
          throw std::runtime_error("sweep: read of a pair column before it was capped");
        for (int xo : xs_at[yo]) {
          const double cl = close_of(xo, yo);
          if (cl < HALF_NEG) continue;
          add_ga(rq[xo], pg[key2(xo, yo)], mb + cl);
        }
      };
      if (rf >= 0) fold(y + 1, gap0, 0.0);
      if (rf >= 1 && y + 2 < n) fold(y + 2, gap1, (dz ? E.s1w[y + 1] : 0.0) + E.w_cp);
      if (full_rf) {
        if (y + 3 < n)
          fold(y + 3, long_acc, (dz ? E.d3w[y + 2] : 0.0) + E.w_cp * (y + 2));
      } else {
        for (int yo = y + 3; yo < n && yo - y - 1 <= rf; ++yo)
          fold(yo, glong_win, (dz ? E.d3w[yo - 1] : 0.0) + E.w_cp * (yo - 1));
      }
      auto& glong = full_rf ? long_acc : glong_win;
      // Right-gap yo count for this column (each pair ending at y receives this many).
      long long n_yo_qm = 0, n_yo_sc = 0;
      if (qm_ix) {
        if (y + 1 < n)
          for (int j : st.prMbarRight[static_cast<size_t>(y + 1)]) {
            if (j < y + 1) continue;
            if (j > y + S) break;
            ++n_yo_qm;
            add_rq(j + 1, st.PrMbar(y + 1, j));
          }
      } else {
        for (int yo = y + 2; yo < n && yo - y - 1 <= S; ++yo) {
          ++n_yo_sc;
          add_rq(yo, st.PrMbar(y + 1, yo - 1));
        }
      }
      // The right-gap sums are not pair states, but each inner pair copies them.
      // Cut each one to n_G first, with the same key as a pair map. long_acc stays cut,
      // so the next column adds into an already capped bucket.
      auto cap_bucket = [&](std::unordered_map<int, Ga>& b) {
        for (auto& kv : b) cap_ga(kv.second);
      };
      cap_bucket(gap0);
      cap_bucket(gap1);
      cap_bucket(glong);
      cap_bucket(rq);
      const double ls = (dz && y + 1 < n ? E.d5w[y + 1] : 0.0) - E.w_cp * y;
      for (int x : xs_at[y]) {
        if (!E.end_ok(x, y)) continue;
        // Full-span slot only (the second sw_fill is the short inner window).
        if (full_rf && E.end_ok(x, y)) {
          if (qm_ix) {
            ++ml_yo_qm_n;
            ml_yo_qm_sum += n_yo_qm;
            if (n_yo_qm > ml_yo_qm_max) ml_yo_qm_max = n_yo_qm;
          } else {
            ++ml_yo_sc_n;
            ml_yo_sc_sum += n_yo_sc;
            if (n_yo_sc > ml_yo_sc_max) ml_yo_sc_max = n_yo_sc;
          }
        }
        const double br = E.w_bp + E.at(x, y);
        if (br < HALF_NEG) continue;
        Ga& dst = pg[key2(x, y)];
        for (int xo = ml_run_lo(x); xo < x; ++xo) {
          auto it = rq.find(xo);
          if (it == rq.end()) continue;
          double lw = E.w_cp * (x - xo - 1);
          if (dz && xo + 1 <= x - 1) lw += E.run(xo + 1, x - 1);
          add_ga(dst, it->second, br + lw);
        }
        auto use_left = [&](int i, double qv) {
          if (i < 1 || x - i > S || qv < HALF_NEG) return;
          const int xo = i - 1;
          const double sc = br + qv;
          if (auto it = gap0.find(xo); it != gap0.end()) add_ga(dst, it->second, sc);
          if (auto it = gap1.find(xo); it != gap1.end()) add_ga(dst, it->second, sc);
          if (auto it = glong.find(xo); it != glong.end()) add_ga(dst, it->second, sc + ls);
          if (auto it = rq.find(xo); it != rq.end()) add_ga(dst, it->second, sc);
        };
        if (qm_ix && x >= 1) {
          for (int i : st.prMbarAt[static_cast<size_t>(x - 1)]) use_left(i, st.PrMbar(i, x - 1));
        } else {
          for (int i = std::max(1, x - S); i <= x - 1; ++i) use_left(i, st.PrMbar(i, x - 1));
        }
      }
      for (int x : xs_at[y]) {
        Ga& m = pg[key2(x, y)];
        if (E.end_ok(x, y)) add1(m, key2(x, y), E.at(x, y));
        cap_ga(m);
        for (int xi = x + 1; xi <= std::min(n - 1, x + 1 + lc); ++xi)
          for (int yi = y - 1; yi >= std::max(0, y - 1 - lc); --yi) {
            const uint64_t ik = key2(xi, yi);
            if (!is_pair.count(ik)) continue;
            if (fin[static_cast<size_t>(yi)])
              throw std::runtime_error("sweep: push into a pair column after it was capped");
            const double s = E.stem_interior(x, xi, yi, y);
            if (s < HALF_NEG) continue;
            add_ga(pg[ik], m, s + pair_bonus);
          }
      }
      fin[static_cast<size_t>(y)] = 1;
    }
    sw_ga[slot] = std::move(pg);
  }

  // Sweep prGamma (and the unused sw_fill) against a fresh lazy growth vector.
  // Both sides use the final scan prefix. Writes --g-eq.
  void g_eq_write() {
    std::ofstream out(g_eq_path);
    if (!out) {
      std::fprintf(stderr, "[g-eq] cannot write %s\n", g_eq_path.c_str());
      return;
    }
    auto rel = [](double a, double b) {
      if (a < HALF_NEG && b < HALF_NEG) return 0.0;
      if (a < HALF_NEG || b < HALF_NEG) return 1e300;
      return std::fabs(std::expm1(a - b));
    };
    auto agg = [&](const std::unordered_map<uint64_t, double>* src, std::vector<double>& v) {
      v.assign(static_cast<size_t>(n), NEG_INF);
      if (!src) return;
      for (const auto& [k, w] : *src) {
        if (w < HALF_NEG) continue;
        const int q = k_lo(k);
        if (q < 0 || q >= n) continue;
        v[static_cast<size_t>(q)] = log_add(v[static_cast<size_t>(q)], w + prefix_log(k_hi(k)));
      }
    };
    auto from_rn = [&](const RN& r, std::vector<double>& v) {
      v.assign(static_cast<size_t>(n), NEG_INF);
      if (g_ntr_K > 0) {
        for (const auto& [q, g] : r.T[0])
          if (q >= 0 && q < n && g > HALF_NEG) v[static_cast<size_t>(q)] = log_add(v[static_cast<size_t>(q)], g);
        return;
      }
      for (int m = 0; m < rnL; ++m) {
        const RV& a = r.L[m];
        for (size_t i = 0; i < a.v.size(); ++i) {
          if (!(a.v[i] > 0.0)) continue;
          const int q = a.lo + static_cast<int>(i);
          if (q < 0 || q >= n) continue;
          v[static_cast<size_t>(q)] = log_add(v[static_cast<size_t>(q)], a.mx + std::log(a.v[i]));
        }
      }
    };
    try {
      sw_fill(0, n);
    } catch (const std::exception& e) {
      out << "fill_error " << e.what() << "\n";
      std::fprintf(stderr, "[g-eq] sw_fill %s\n", e.what());
    }
    const bool saved = p2_sweep;
    p2_sweep = false;
    prRT.clear();
    prRR.clear();
    prRG.clear();
    prRH.clear();
    std::vector<std::pair<int, int>> pairs;
    {
      std::unordered_set<uint64_t> seen;
      const std::vector<std::vector<int>>& src = pr_pair_list(true);
      for (int x = 0; x < n; ++x)
        for (int y : src[x]) {
          const uint64_t k = key2(x, y);
          if (seen.insert(k).second) pairs.emplace_back(x, y);
        }
      if (stem_)
        for (const auto& kv : stem_->prGamma) {
          const uint64_t k = kv.first;
          if (seen.insert(k).second) pairs.emplace_back(k_hi(k), k_lo(k));
        }
    }
    std::fprintf(stderr, "[g-eq] pairs %zu\n", pairs.size());
    struct Bad {
      int x, y, q;
      double hy, gr, ad, r;
    };
    std::vector<Bad> top;
    auto consider = [&](int x, int y, int q, double hy, double gr) {
      const double r = rel(hy, gr);
      if (!(r > 1e-9)) return;
      const double ad = (hy < HALF_NEG || gr < HALF_NEG) ? 1e300 : std::fabs(hy - gr);
      top.push_back({x, y, q, hy, gr, ad, r});
    };
    long long ncmp = 0, nbad = 0, ncmp_f = 0, nbad_f = 0;
    std::unordered_set<uint64_t> bad_pair;
    std::vector<double> hy, gr, fl, spec_h, spec_g;
    bool have_spec = false;
    int built = 0;
    for (const auto& [x, y] : pairs) {
      const RN& rn = rec_root(x, y, false);
      from_rn(rn, gr);
      const std::unordered_map<uint64_t, double>* hs = nullptr;
      {
        auto it = sw_ga[0].find(key2(x, y));
        if (it != sw_ga[0].end()) hs = &it->second;
      }
      agg(hs, hy);
      const std::unordered_map<uint64_t, double>* fs = nullptr;
      auto fit = sw_ga[0].find(key2(x, y));
      if (fit != sw_ga[0].end()) fs = &fit->second;
      agg(fs, fl);
      bool any = false;
      for (int q = 0; q < n; ++q) {
        const double a = hy[static_cast<size_t>(q)], b = gr[static_cast<size_t>(q)];
        if (a < HALF_NEG && b < HALF_NEG) continue;
        ++ncmp;
        if (rel(a, b) > 1e-9) {
          ++nbad;
          any = true;
          consider(x, y, q, a, b);
        }
        const double c = fl[static_cast<size_t>(q)];
        if (c < HALF_NEG && b < HALF_NEG) continue;
        ++ncmp_f;
        if (rel(c, b) > 1e-9) ++nbad_f;
      }
      if (any) bad_pair.insert(key2(x, y));
      if (x == 4 && y == 45) {
        spec_h = hy;
        spec_g = gr;
        have_spec = true;
      }
      if (++built % 400 == 0) std::fprintf(stderr, "[g-eq] built %d\n", built);
    }
    p2_sweep = saved;
    long long ncmp_p = 0, nbad_p = 0, nbad_pf = 0;
    {
      geq_pre = true;
      p2_sweep = false;
      prRT.clear();
      prRR.clear();
      prRG.clear();
      prRH.clear();
      int built_p = 0;
      for (const auto& [x, y] : pairs) {
        const RN& rn = rec_root(x, y, false);
        from_rn(rn, gr);
        const std::unordered_map<uint64_t, double>* hs = nullptr;
        {
          auto it = sw_ga[0].find(key2(x, y));
          if (it != sw_ga[0].end()) hs = &it->second;
        }
        agg(hs, hy);
        const std::unordered_map<uint64_t, double>* fs = nullptr;
        auto fit = sw_ga[0].find(key2(x, y));
        if (fit != sw_ga[0].end()) fs = &fit->second;
        agg(fs, fl);
        for (int q = 0; q < n; ++q) {
          const double a = hy[static_cast<size_t>(q)], b = gr[static_cast<size_t>(q)],
                       c = fl[static_cast<size_t>(q)];
          if (a < HALF_NEG && b < HALF_NEG) continue;
          ++ncmp_p;
          if (rel(a, b) > 1e-9) ++nbad_p;
          if (!(c < HALF_NEG && b < HALF_NEG) && rel(c, b) > 1e-9) ++nbad_pf;
        }
        if (x == 1 && y == 5) {
          out << "pre15\n";
          for (int q = 0; q < n; ++q) {
            const double a = hy[static_cast<size_t>(q)], b = gr[static_cast<size_t>(q)],
                         c = fl[static_cast<size_t>(q)];
            if (a < HALF_NEG && b < HALF_NEG && c < HALF_NEG) continue;
            out << "pq " << q << " " << std::setprecision(17) << a << " " << c << " " << b << "\n";
          }
        }
        if (++built_p % 400 == 0) std::fprintf(stderr, "[g-eq] pre %d\n", built_p);
      }
      geq_pre = false;
      prRT.clear();
      prRR.clear();
      prRG.clear();
      prRH.clear();
      p2_sweep = saved;
    }
    out << "pre " << ncmp_p << " " << nbad_p << " " << nbad_pf << "\n";
    std::sort(top.begin(), top.end(), [](const Bad& a, const Bad& b) {
      if (a.ad != b.ad) return a.ad > b.ad;
      if (a.y != b.y) return a.y < b.y;
      if (a.x != b.x) return a.x < b.x;
      return a.q < b.q;
    });
    if (top.size() > 10) top.resize(10);
    out << "cmp " << ncmp << " " << nbad << "\n";
    out << "fill " << ncmp_f << " " << nbad_f << "\n";
    out << "pairs " << pairs.size() << " bad_pairs " << bad_pair.size() << "\n";
    for (size_t i = 0; i < top.size(); ++i)
      out << "top " << i << " " << top[i].x << " " << top[i].y << " " << top[i].q << " " << std::setprecision(17)
          << top[i].hy << " " << top[i].gr << " " << top[i].ad << " " << top[i].r << "\n";
    out << "vec " << (have_spec ? 1 : 0) << "\n";
    if (have_spec) {
      for (int q = 0; q < n; ++q) {
        const double a = spec_h[static_cast<size_t>(q)], b = spec_g[static_cast<size_t>(q)];
        if (a < HALF_NEG && b < HALF_NEG) continue;
        out << "vq " << q << " " << std::setprecision(17) << a << " " << b << "\n";
      }
    }
    std::unordered_map<uint64_t, int> depth;
    std::function<int(int, int)> dep = [&](int x, int y) -> int {
      const uint64_t k = key2(x, y);
      auto it = depth.find(k);
      if (it != depth.end()) return it->second;
      depth.emplace(k, 0);
      std::vector<ChainNbr> nb;
      pr_steps(x, y, false, true, nb);
      int d = 0;
      if (!nb.empty()) {
        d = n + 1;
        for (const ChainNbr& t : nb) d = std::min(d, 1 + dep(t.x, t.y));
      }
      depth[k] = d;
      return d;
    };
    int bx = -1, by = -1, bd = n + 2;
    for (uint64_t k : bad_pair) {
      const int x = k_hi(k), y = k_lo(k);
      const int d = dep(x, y);
      if (d < bd || (d == bd && (y < by || (y == by && x < bx)))) {
        bd = d;
        bx = x;
        by = y;
      }
    }
    out << "short " << bx << " " << by << " " << (bx < 0 ? -1 : bd) << "\n";
    p2_sweep = false;
    long long ibad = 0;
    int ix = -1, iy = -1, idp = n + 2;
    for (uint64_t k : bad_pair) {
      const int x = k_hi(k), y = k_lo(k);
      std::vector<ChainNbr> nb;
      pr_steps(x, y, false, true, nb);
      bool ml = false;
      for (const ChainNbr& t : nb)
        if (t.t != 0) ml = true;
      if (ml) continue;
      ++ibad;
      const int d = dep(x, y);
      if (d < idp || (d == idp && (y < iy || (y == iy && x < ix)))) {
        idp = d;
        ix = x;
        iy = y;
      }
    }
    out << "ibad " << ibad << " ishort " << ix << " " << iy << " " << (ix < 0 ? -1 : idp) << "\n";
    if (bx >= 0) {
      int x = bx, y = by;
      for (int hop = 0; hop < n && x >= 0; ++hop) {
        const RN& rn = rec_root(x, y, false);
        from_rn(rn, gr);
        const std::unordered_map<uint64_t, double>* hs = nullptr;
        {
          auto it = sw_ga[0].find(key2(x, y));
          if (it != sw_ga[0].end()) hs = &it->second;
        }
        agg(hs, hy);
        int nd = 0;
        for (int q = 0; q < n; ++q)
          if (rel(hy[static_cast<size_t>(q)], gr[static_cast<size_t>(q)]) > 1e-9) ++nd;
        std::vector<ChainNbr> nb;
        pr_steps(x, y, false, true, nb);
        out << "node " << hop << " " << x << " " << y << " depth " << dep(x, y) << " ndiff " << nd << " nstep "
            << nb.size() << " at " << std::setprecision(17) << (E.end_ok(x, y) ? E.at(x, y) : NEG_INF) << " pref "
            << prefix_log(x) << "\n";
        if (hop == 0) {
          for (int q = 0; q < n; ++q) {
            const double a = hy[static_cast<size_t>(q)], b = gr[static_cast<size_t>(q)];
            if (a < HALF_NEG && b < HALF_NEG) continue;
            out << "sq " << q << " " << std::setprecision(17) << a << " " << b << "\n";
            if (!(rel(a, b) > 1e-9)) continue;
            out << "diffq " << q << " " << std::setprecision(17) << a << " " << b << "\n";
            if (hs) {
              for (const auto& [k, w] : *hs)
                if (k_lo(k) == q)
                  out << "raw " << k_hi(k) << " " << q << " " << std::setprecision(17) << w << "\n";
            }
            for (size_t s = 0; s < nb.size(); ++s) {
              const RN& ch = rec_root(nb[s].x, nb[s].y, false);
              std::vector<double> cg;
              from_rn(ch, cg);
              const double cv = cg[static_cast<size_t>(q)];
              if (cv < HALF_NEG) continue;
              out << "via " << s << " " << q << " " << std::setprecision(17) << (cv + nb[s].w) << "\n";
            }
          }
        }
        int nxt = -1, ny = -1, bestd = n + 3;
        for (size_t s = 0; s < nb.size(); ++s) {
          const int d = dep(nb[s].x, nb[s].y);
          out << "step " << s << " t " << nb[s].t << " " << nb[s].x << " " << nb[s].y << " " << std::setprecision(17)
              << nb[s].w << " depth " << d << " bad " << (bad_pair.count(key2(nb[s].x, nb[s].y)) ? 1 : 0) << "\n";
          if (d < bestd) {
            bestd = d;
            nxt = nb[s].x;
            ny = nb[s].y;
          }
        }
        if (nxt < 0 || dep(x, y) == 0) break;
        x = nxt;
        y = ny;
      }
    }
    out.flush();
    std::fprintf(stderr, "[g-eq] cmp %lld bad %lld fill %lld bad %lld\n", ncmp, nbad, ncmp_f, nbad_f);
  }

  void pr_build_static() {
    if (pr_built) return;
    pr_built = true;
    SecTimer ts(pr_static_sec);
    const StemEst& st = *stem_;
    std::vector<std::pair<int, int>> pairs;
    const std::vector<std::vector<int>>& vsrc = pr_pair_list(false);
    for (int x = 0; x < n; ++x)
      for (int y : vsrc[x]) pairs.emplace_back(x, y);
    if (!pr_v2dump.empty()) {
      pr_v2_dump(pairs);
      std::fflush(nullptr);
      std::_Exit(0);
    }
    std::stable_sort(pairs.begin(), pairs.end(),
                     [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
                       return a.second - a.first > b.second - b.first;
                     });
    std::vector<AT> ats;
    std::vector<std::vector<AT>> bk;
    if (pr_vstream) bk.resize(n);
    if (lpf) {
      vi_pairs = static_cast<long long>(pairs.size());
      vi_pairs_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - ts.t0).count();
    }
    {
      std::unique_ptr<SecTimer> vtm(lpf ? new SecTimer(vi_vbuild_s) : nullptr);
      pr_static_ats(pr_pmax, true, pairs, ats, pr_vstream ? &bk : nullptr);
    }
    if (sw_check) {
      sw_pairs = pairs;
      sw_verify();
      if (p2_vsweep && stem_ && !stem_->prMu.empty() && !sw_ref.empty()) {
        auto rel = [](double a, double b) {
          if (a < HALF_NEG && b < HALF_NEG) return 0.0;
          if (a < HALF_NEG || b < HALF_NEG) return 1e300;
          return std::fabs(std::expm1(a - b));
        };
        double worst = 0, wgot = 0, wref = 0;
        int wx = -1, wy = -1, wb = -1;
        int fx = -1, fy = -1, fb = -1;
        double fgot = 0, fref = 0;
        long long nent = 0, miss = 0;
        for (const auto& [pk, p] : sw_ref) {
          const int x = k_hi(pk), y = k_lo(pk);
          auto it = stem_->prMu.find(pk);
          for (size_t i = 0; i < p.k.size(); ++i) {
            const int b = k_hi(p.k[i]);
            const double ref = p.w[i];
            double got = NEG_INF;
            if (it != stem_->prMu.end()) {
              auto jt = it->second.find(key2(b, 0));
              if (jt != it->second.end()) got = jt->second;
            }
            if (got < HALF_NEG && ref < HALF_NEG) continue;
            ++nent;
            if (got < HALF_NEG || ref < HALF_NEG) ++miss;
            const double d = rel(got, ref);
            if (d > worst) {
              worst = d;
              wx = x;
              wy = y;
              wb = b;
              wgot = got;
              wref = ref;
            }
            if (d > 1e-12 && (fx < 0 || y < fy || (y == fy && x < fx) || (y == fy && x == fx && b < fb))) {
              fx = x;
              fy = y;
              fb = b;
              fgot = got;
              fref = ref;
            }
          }
        }
        std::fprintf(stderr,
                     "[v-sweep-check] entries %lld one_sided %lld max_rel %.3e at %d.%d b %d "
                     "sweep %.9g ref %.9g first %d.%d b %d sweep %.9g ref %.9g\n",
                     nent, miss, worst, wx, wy, wb, wgot, wref, fx, fy, fb, fgot, fref);
      }
    }
    std::unique_ptr<SecTimer> atm(lpf ? new SecTimer(vi_arm_s) : nullptr);
    mtrim();
    long long bk_n = 0;
    double bk_b = 0.0;
    for (const std::vector<AT>& v : bk) {
      bk_n += static_cast<long long>(v.size());
      bk_b += msz::vec(v);
    }
    mem_log("V_ats_done", {{"V_ats", static_cast<long long>(ats.size()) + bk_n, msz::vec(ats) + bk_b}});
    std::sort(ats.begin(), ats.end(), at_less);
    if (pr_pmax > 0 && pr_pcheck && n <= 40 && pr_cut <= 0.0) {
      std::map<std::tuple<int, int, int, int>, double> en;
      std::vector<ChainNbr> nb0;
      std::function<void(int, int, int, int, int, double)> dfs = [&](int x0, int y0, int x,
                                                                    int y, int c, double w) {
        const bool inner = c == 1;
        if (inner ? single_stems : E.end_ok(x, y)) {
          auto key = std::make_tuple(x, x0, y0, y);
          auto it = en.find(key);
          const double v = w + E.at(x, y);
          if (it == en.end()) en.emplace(key, v);
          else it->second = log_add(it->second, v);
        }
        if (c >= pr_pmax) return;
        std::vector<ChainNbr> nb;
        pr_steps(x, y, inner, false, nb);
        for (const ChainNbr& t : nb) dfs(x0, y0, t.x, t.y, c + 1, w + t.w);
      };
      for (const auto& [x, y] : pairs)
        if (E.end_ok(x, y)) dfs(x, y, x, y, 1, (E.at(x, y) + E.w_stem) + pair_bonus);
      std::vector<AT> ref;
      for (const auto& [k, w] : en)
        ref.push_back({std::get<0>(k), std::get<1>(k), std::get<2>(k), std::get<3>(k), w});
      std::sort(ref.begin(), ref.end(), at_less);
      size_t a = 0, b = 0;
      while (a < ats.size() || b < ref.size()) {
        if (b == ref.size() || (a < ats.size() && at_less(ats[a], ref[b]))) {
          ++pe_key_mis;
          ++a;
        } else if (a == ats.size() || at_less(ref[b], ats[a])) {
          ++pe_key_mis;
          ++b;
        } else {
          ++pe_n;
          pe_maxd = std::max(pe_maxd, std::fabs(ats[a].w - ref[b].w));
          ++a;
          ++b;
        }
      }
    }
    if (pr_pmax > 0 && pr_pcheck) {
      std::vector<AT> ref;
      pr_static_ats(0, false, pairs, ref);
      std::sort(ref.begin(), ref.end(), at_less);
      size_t a = 0, b = 0;
      while (a < ats.size() || b < ref.size()) {
        if (b == ref.size() || (a < ats.size() && at_less(ats[a], ref[b]))) {
          ++pc_key_mis;
          ++a;
        } else if (a == ats.size() || at_less(ref[b], ats[a])) {
          ++pc_key_mis;
          ++b;
        } else {
          ++pc_n;
          pc_maxd = std::max(pc_maxd, std::fabs(ats[a].w - ref[b].w));
          ++a;
          ++b;
        }
      }
    }
    if (!pr_partner_ml && pr_pmlcheck && pr_pmax == 0 && n <= 40 && pr_cut <= 0.0 &&
        pr_veps <= 0.0) {
      std::map<std::tuple<int, int, int, int>, double> en;
      std::function<void(int, int, int, int, bool, double)> dfs = [&](int x0, int y0, int x,
                                                                     int y, bool inner, double w) {
        if (inner ? single_stems : E.end_ok(x, y)) {
          auto key = std::make_tuple(x, x0, y0, y);
          auto it = en.find(key);
          const double v = w + E.at(x, y);
          if (it == en.end()) en.emplace(key, v);
          else it->second = log_add(it->second, v);
        }
        for (const auto& [xo, yo] : pairs) {
          if (!(xo < x && yo > y)) continue;
          const double s = E.stem_interior(xo, x, y, yo);
          if (s > HALF_NEG) dfs(x0, y0, xo, yo, false, w + (s + pair_bonus));
        }
      };
      for (const auto& [x, y] : pairs)
        if (E.end_ok(x, y)) dfs(x, y, x, y, true, (E.at(x, y) + E.w_stem) + pair_bonus);
      std::vector<AT> ref;
      for (const auto& [k, w] : en)
        ref.push_back({std::get<0>(k), std::get<1>(k), std::get<2>(k), std::get<3>(k), w});
      std::sort(ref.begin(), ref.end(), at_less);
      size_t a = 0, b = 0;
      while (a < ats.size() || b < ref.size()) {
        if (b == ref.size() || (a < ats.size() && at_less(ats[a], ref[b]))) {
          ++pm_key_mis;
          ++a;
        } else if (a == ats.size() || at_less(ref[b], ats[a])) {
          ++pm_key_mis;
          ++b;
        } else {
          ++pm_n;
          pm_maxd = std::max(pm_maxd, std::fabs(ats[a].w - ref[b].w));
          ++a;
          ++b;
        }
      }
    }
    // pv: the previous entry in at_less order (null at the first)
    auto add = [&](const AT& a, const AT* pv) {
      const bool new_arm = !pv || a.xs != pv->xs || a.xp != pv->xp;
      if (new_arm) {
        prArmXs.push_back(a.xs);
        prArmXp.push_back(a.xp);
        prArmE0.push_back(static_cast<int>(pr_vxs ? prYy.size() : prEy.size()));
        prArmY0.push_back(static_cast<int>(prYy.size()));
      }
      if (!pr_vxs) {
        prEy.push_back(a.yp);
        prEys.push_back(a.ys);
        prEw.push_back(a.w);
      }
      const double yw = pr_vxs ? a.w : a.w + (E.w_beta1 + st.prS[a.ys + 1]);
      if (new_arm || a.yp != pv->yp) {
        prYy.push_back(a.yp);
        prYw.push_back(yw);
      } else {
        prYw.back() = log_add(prYw.back(), yw);
      }
    };
    for (size_t k = 0; k < ats.size(); ++k) add(ats[k], k == 0 ? nullptr : &ats[k - 1]);
    if (pr_vstream) {
      if (pr_vxs) {
        prYy.reserve(static_cast<size_t>(bk_n));
        prYw.reserve(static_cast<size_t>(bk_n));
      } else {
        prEy.reserve(static_cast<size_t>(bk_n));
        prEys.reserve(static_cast<size_t>(bk_n));
        prEw.reserve(static_cast<size_t>(bk_n));
      }
      for (int xs = 0; xs < n; ++xs) {
        std::vector<AT>& v = bk[xs];
        std::sort(v.begin(), v.end(), at_less);
        for (size_t k = 0; k < v.size(); ++k) add(v[k], k == 0 ? nullptr : &v[k - 1]);
        std::vector<AT>().swap(v);
        if ((xs & 31) == 31) mtrim();
      }
      std::vector<std::vector<AT>>().swap(bk);
      mtrim();
    }
    prA = static_cast<int>(prArmXs.size());
    prArmE0.push_back(static_cast<int>(pr_vxs ? prYy.size() : prEy.size()));
    prArmY0.push_back(static_cast<int>(prYy.size()));
    prXsFirst.assign(n + 1, prA);
    for (int r = prA - 1; r >= 0; --r) prXsFirst[prArmXs[r]] = r;
    for (int x = n - 1; x >= 0; --x) prXsFirst[x] = std::min(prXsFirst[x], prXsFirst[x + 1]);
    prByXp.resize(prA);
    for (int r = 0; r < prA; ++r) prByXp[r] = r;
    std::stable_sort(prByXp.begin(), prByXp.end(),
                     [&](int a, int b) { return prArmXp[a] < prArmXp[b]; });
    prXpFirst.assign(n + 1, prA);
    for (int r = prA - 1; r >= 0; --r) prXpFirst[prArmXp[prByXp[r]]] = r;
    for (int x = n - 1; x >= 0; --x) prXpFirst[x] = std::min(prXpFirst[x], prXpFirst[x + 1]);
    pr_arms = prA;
    pr_arm_entries = static_cast<long long>(pr_vxs ? prYy.size() : prEy.size());
    pr_yw_from_sweep();
    if (k_bound) k_bound_check();
    if (g_mem_log) {
      std::vector<MemItem> it{{"V_ats", static_cast<long long>(ats.size()), msz::vec(ats)}};
      mem_items(it, "", false, true);
      it.resize(2);
      mem_log("V_arms_built", it);
    }
    if (pr_vstream || g_mtrim) {
      std::vector<AT>().swap(ats);
      mtrim();
    }
    if (!pr_vdiag.empty() || pr_vsample_p > 0.0) pr_vdiag_static(pairs);
    if (g_static_only) {
      const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - ts.t0).count();
      std::fprintf(stderr,
                   "[vstatic] n %d pairs %zu V_entries %lld arm_entries %lld arms %d adds %lld "
                   "node_terms %lld hv_built %lld hv_entries %lld hv_terms %lld fut_built %lld "
                   "fut_int %lld fut_ml %lld v_before %lld sec %.6f\n",
                   n, pairs.size(), pr_V_entries, pr_arm_entries, prA, pr_adds, pr_nterms[1],
                   pr_hbuilt[1], pr_hentries[1], pr_hterms[1], st_nodes[1], st_kind[1][0],
                   st_kind[1][1], v_before, sec);
      std::fflush(nullptr);
      std::_Exit(0);
    }
    if (p2_sweep || p2_vsweep) {
      long long ent = 0;
      if (stem_)
        for (const auto& kv : stem_->prGamma) ent += static_cast<long long>(kv.second.size());
      std::fprintf(stderr,
                   "[sweep] pairs %zu gamma %lld g_states %lld v_states %lld pushes %lld "
                   "cap %d qm_max %d sec %.3f\n",
                   pairs.size(), ent, stem_ ? stem_->pr_sw_states_g : 0,
                   stem_ ? stem_->pr_sw_states_v : 0, stem_ ? stem_->pr_sw_pushes : 0, g_ntr_K,
                   stem_ ? stem_->pr_mbar_max : 0, stem_ ? stem_->pr_sw_sec : 0.0);
    }
  }

  // Per V group g (prYy index): log Vr, the interior-only V and log sum w len (len = loops
  // from xs.ys down to x'.y'); arm ids by (xs, x'); the full-support Vr groups not in V.
  void pr_vdiag_static(const std::vector<std::pair<int, int>>& pairs) {
    const StemEst& st = *stem_;
    auto t0 = std::chrono::steady_clock::now();
    pr_vsuf_build();
    const size_t G = prYy.size();
    prVrg.assign(G, NEG_INF);
    prVr2g.assign(G, NEG_INF);
    prVsepg.assign(G, NEG_INF);
    prVog.assign(G, NEG_INF);
    prVlbg.assign(G, NEG_INF);
    std::unordered_map<uint64_t, int> gix;
    gix.reserve(G * 2 + 1);
    for (int r = 0; r < prA; ++r) {
      prArmOf[key2(prArmXs[r], prArmXp[r])] = r;
      for (int g = prArmY0[r]; g < prArmY0[r + 1]; ++g) {
        prVrg[g] = pr_vr_group(prArmXs[r], prArmXp[r], prYy[g]);
        prVsepg[g] = pr_vb_group(prArmXs[r], prArmXp[r], prYy[g], 0);
        prVr2g[g] = pr_vb_group(prArmXs[r], prArmXp[r], prYy[g], 1);
        prVlbg[g] = pr_vb_group(prArmXs[r], prArmXp[r], prYy[g], 2);
        prVog[g] = pr_vb_group(prArmXs[r], prArmXp[r], prYy[g], 3);
        gix[key3(prArmXs[r], prArmXp[r], prYy[g])] = g;
      }
    }
    for (const auto& [x, y] : pairs) {
      if (!E.end_ok(x, y)) continue;
      for (int xs = 0; xs <= x; ++xs) {
        if (pr_vr_group(xs, x, y) < HALF_NEG) continue;
        ++pr_vr_groups;
        if (!gix.count(key3(xs, x, y))) ++pr_vr_extra;
      }
    }
    pr_vr_sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    t0 = std::chrono::steady_clock::now();
    // interior-only V (pr_ml off): profiles of w and of log sum w len
    const bool ml0 = pr_ml;
    pr_ml = false;
    std::unordered_map<uint64_t, std::pair<PProf, PProf>> VI;
    VI.reserve(pairs.size() * 2 + 1);
    std::vector<ChainNbr> nb;
    auto merge2 = [&](PAcc& aw, PAcc& al) {
      for (const ChainNbr& t : nb) {
        const auto& q = VI.find(key2(t.x, t.y))->second;
        for (size_t k = 0; k < q.first.k.size(); ++k) aw.add(q.first.k[k], t.w + q.first.w[k]);
        for (size_t k = 0; k < q.first.k.size(); ++k) al.add(q.first.k[k], t.w + q.first.w[k]);
        for (size_t k = 0; k < q.second.k.size(); ++k) al.add(q.second.k[k], t.w + q.second.w[k]);
      }
    };
    for (const auto& [x, y] : pairs) {
      pr_steps(x, y, false, false, nb);
      PAcc aw, al;
      if (E.end_ok(x, y)) aw.add(key2(x, y), E.at(x, y));
      merge2(aw, al);
      std::pair<PProf, PProf> p;
      aw.finish(p.first, 0.0);
      al.finish(p.second, 0.0);
      VI.emplace(key2(x, y), std::move(p));
    }
    prVint.assign(G, NEG_INF);
    prVlen.assign(G, NEG_INF);
    for (const auto& [x, y] : pairs) {
      if (!E.end_ok(x, y)) continue;
      pr_steps(x, y, true, false, nb);
      PAcc aw, al;
      if (single_stems) aw.add(key2(x, y), E.at(x, y));
      merge2(aw, al);
      PProf pw, pl;
      aw.finish(pw, 0.0);
      al.finish(pl, 0.0);
      const double base = (E.at(x, y) + E.w_stem) + pair_bonus;
      for (size_t k = 0; k < pw.k.size(); ++k) {
        auto it = gix.find(key3(k_hi(pw.k[k]), x, y));
        if (it == gix.end()) continue;
        const double f = base + (E.w_beta1 + st.prS[k_lo(pw.k[k]) + 1]);
        prVint[it->second] = log_add(prVint[it->second], pw.w[k] + f);
      }
      for (size_t k = 0; k < pl.k.size(); ++k) {
        auto it = gix.find(key3(k_hi(pl.k[k]), x, y));
        if (it == gix.end()) continue;
        const double f = base + (E.w_beta1 + st.prS[k_lo(pl.k[k]) + 1]);
        prVlen[it->second] = log_add(prVlen[it->second], pl.w[k] + f);
      }
    }
    pr_ml = ml0;
    pr_vint_sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    prVuse.assign(G, 0);
    prVcon.assign(G, NEG_INF);
    prVshare.assign(G, 0.0);
  }

  // A chosen arm set of an L-lookup (q'; u, v): per group of the selected arms with y' > q',
  // the use, log sum of Qzbar[q'+1][y'-1] V (no gap fillers), and its share of the lookup's
  // exterior total (with the fillers).  Sampled lookups: every arm of the gap ranked by its
  // partner value with V and with Vr (full support).
  std::vector<std::pair<int, double>> prVtmp;
  void pr_vdiag_T(int qp, int u, int v, bool two, const std::vector<int>& sel) {
    const StemEst& st = *stem_;
    ++pr_vT_n;
    if (!prVuse.empty()) {
      prVtmp.clear();
      double tot = NEG_INF;
      for (int r : sel) {
        const double gl = two ? 0.0 : pr_qz(u + 1, prArmXs[r] - 1);
        const double gr = pr_qz(prArmXp[r] + 1, v - 1);
        if (gl < HALF_NEG || gr < HALF_NEG) continue;
        const int hi = prArmY0[r + 1];
        for (int g = static_cast<int>(std::upper_bound(prYy.begin() + prArmY0[r], prYy.begin() + hi, qp) -
                                      prYy.begin());
             g < hi; ++g) {
          const double z = st.Qzbar(qp + 1, prYy[g] - 1);
          if (z < HALF_NEG) continue;
          const double c = z + prYw[g];
          prVtmp.emplace_back(g, c + gl + gr);
          tot = log_add(tot, c + gl + gr);
          ++prVuse[g];
          prVcon[g] = log_add(prVcon[g], c);
        }
      }
      if (tot > HALF_NEG)
        for (const auto& [g, t] : prVtmp) prVshare[g] += std::exp(t - tot);
    }
    if (!pr_vsel || pr_vsel_n >= pr_vsample_max || pr_vrand() >= pr_vsample_p) return;
    ++pr_vsel_n;
    // exact partner value per arm of the gap
    std::vector<std::pair<double, int>> pv;
    std::unordered_map<int, double> pvOf;
    int lo, hi;
    if (two) {
      lo = prXsFirst[u + 1];
      hi = prXsFirst[u + 2];
    } else {
      lo = prXpFirst[u + 1];
      hi = prXpFirst[v];
    }
    double totV = NEG_INF;
    for (int pos = lo; pos < hi; ++pos) {
      const int r = two ? pos : prByXp[pos];
      if (prArmXp[r] > v - 1 || (!two && prArmXs[r] < u + 1)) continue;
      const double gl = two ? 0.0 : pr_qz(u + 1, prArmXs[r] - 1);
      const double gr = pr_qz(prArmXp[r] + 1, v - 1);
      if (gl < HALF_NEG || gr < HALF_NEG) continue;
      double s = NEG_INF;
      const int gh = prArmY0[r + 1];
      for (int g = static_cast<int>(std::upper_bound(prYy.begin() + prArmY0[r], prYy.begin() + gh, qp) -
                                    prYy.begin());
           g < gh; ++g) {
        const double z = st.Qzbar(qp + 1, prYy[g] - 1);
        if (z > HALF_NEG) s = log_add(s, z + prYw[g]);
      }
      if (s < HALF_NEG) continue;
      s += gl + gr;
      pv.emplace_back(s, r);
      pvOf[r] = s;
      totV = log_add(totV, s);
    }
    // partner value per arm (xs, x') of the gap, full support: kind -1 Vr, else pr_vb_group's
    const auto& KR = pr_KRp();
    auto rank = [&](int kind, std::vector<std::tuple<double, int, int>>& pr) {
      pr.clear();
      for (int xp = u + 1; xp <= v - 1; ++xp) {
        const int xs_lo = u + 1, xs_hi = two ? std::min(u + 1, xp) : xp;
        if (xs_lo > xs_hi) continue;
        const double gr = pr_qz(xp + 1, v - 1);
        if (gr < HALF_NEG) continue;
        std::vector<double> acc(xs_hi - xs_lo + 1, NEG_INF);
        for (auto it = std::upper_bound(KR[xp].begin(), KR[xp].end(), qp); it != KR[xp].end(); ++it) {
          const int yp = *it;
          if (!E.end_ok(xp, yp)) continue;
          const double z = st.Qzbar(qp + 1, yp - 1);
          if (z < HALF_NEG) continue;
          for (int xs = xs_lo; xs <= xs_hi; ++xs) {
            const double w = kind < 0 ? pr_vr_group(xs, xp, yp) : pr_vb_group(xs, xp, yp, kind);
            if (w > HALF_NEG) acc[xs - xs_lo] = log_add(acc[xs - xs_lo], z + w);
          }
        }
        for (int xs = xs_lo; xs <= xs_hi; ++xs) {
          if (acc[xs - xs_lo] < HALF_NEG) continue;
          const double gl = two ? 0.0 : pr_qz(u + 1, xs - 1);
          if (gl < HALF_NEG) continue;
          pr.emplace_back(acc[xs - xs_lo] + gl + gr, xs, xp);
        }
      }
      std::sort(pr.begin(), pr.end(), [](const auto& a, const auto& b) { return std::get<0>(a) > std::get<0>(b); });
    };
    const size_t K = static_cast<size_t>(pr_arms_k > 0 ? pr_arms_k : 64);
    std::sort(pv.begin(), pv.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
    std::unordered_set<int> topV;
    double keepV = NEG_INF, keepE = NEG_INF;
    for (size_t k = 0; k < std::min(K, pv.size()); ++k) {
      topV.insert(pv[k].second);
      keepV = log_add(keepV, pv[k].first);
    }
    // per ranking: arms ranked, overlap with V's top K, of those in V, log kept exact V
    std::vector<std::tuple<double, int, int>> pr;
    size_t nr[3], ov[3], inv[3];
    double keep[3];
    const int KINDS[3] = {-1, 1, 0};
    std::unordered_map<uint64_t, std::array<double, 4>> armv;  // exact, Vr, Vr2, Vsep
    if (pr_varms)
      for (const auto& [s, r] : pv) armv[key2(prArmXs[r], prArmXp[r])] = {s, NEG_INF, NEG_INF, NEG_INF};
    for (int c = 0; c < 3; ++c) {
      rank(KINDS[c], pr);
      if (pr_varms)
        for (const auto& [val, a, b] : pr) {
          auto it = armv.try_emplace(key2(a, b), std::array<double, 4>{NEG_INF, NEG_INF, NEG_INF, NEG_INF}).first;
          it->second[1 + c] = val;
        }
      nr[c] = pr.size();
      ov[c] = inv[c] = 0;
      keep[c] = NEG_INF;
      for (size_t k = 0; k < std::min(K, pr.size()); ++k) {
        auto ia = prArmOf.find(key2(std::get<1>(pr[k]), std::get<2>(pr[k])));
        if (ia == prArmOf.end()) continue;
        auto ip = pvOf.find(ia->second);
        if (ip == pvOf.end()) continue;
        ++inv[c];
        keep[c] = log_add(keep[c], ip->second);
        if (topV.count(ia->second)) ++ov[c];
      }
    }
    for (int r : sel) {
      auto ip = pvOf.find(r);
      if (ip != pvOf.end()) keepE = log_add(keepE, ip->second);
    }
    auto sh = [&](double x) { return totV > HALF_NEG && x > HALF_NEG ? std::exp(x - totV) : 0.0; };
    std::fprintf(pr_vsel, "%d\t%d\t%d\t%d\t%zu\t%zu\t%zu\t%zu\t%.9g\t%.9g\t%.9g\t%zu\t%zu\t%.9g\t%zu\t%zu\t%.9g\n",
                 qp, u, v, two ? 1 : 0, pv.size(), nr[0], ov[0], inv[0], sh(keep[0]), sh(keepV), sh(keepE),
                 ov[1], inv[1], sh(keep[1]), ov[2], inv[2], sh(keep[2]));
    if (pr_varms)
      for (const auto& [k, a] : armv)
        std::fprintf(pr_varms, "%lld\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\n", pr_vsel_n, k_hi(k), k_lo(k), a[0],
                     a[1], a[2], a[3]);
  }

  void pr_vdump(const std::string& sid) {
    if (prVrg.empty()) return;
    const StemEst& st = *stem_;
    if (!pr_vdiag.empty()) {
      const std::string cmd = "gzip -1 > '" + pr_vdiag + ".ventries.tsv.gz'";
      FILE* f = popen(cmd.c_str(), "w");
      std::fprintf(f, "xs\txp\typ\tlogV\tlogVr\tlogVint\tlogVlen\tlpn\tuse\tlogcon\tshare\t"
                      "logVr2\tlogVsep\tlogVo\tlogLB\n");
      for (int r = 0; r < prA; ++r)
        for (int g = prArmY0[r]; g < prArmY0[r + 1]; ++g) {
          const double lp = st.prLP.empty() ? NAN : st.prLP[static_cast<size_t>(prArmXp[r]) * n + prYy[g]];
          std::fprintf(f, "%d\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.6g\t%lld\t%.9g\t%.6g\t%.9g\t%.9g\t%.9g\t%.9g\n",
                       prArmXs[r], prArmXp[r], prYy[g], prYw[g], prVrg[g], prVint[g], prVlen[g], lp,
                       prVuse[g], prVcon[g], prVshare[g], prVr2g[g], prVsepg[g], prVog[g], prVlbg[g]);
        }
      pclose(f);
    }
    std::fprintf(stderr,
                 "[vratio] %s groups %zu vr_groups %lld vr_extra %lld used %lld T_seen %lld "
                 "sel_samples %lld vr_sec %.3f vint_sec %.3f\n",
                 sid.c_str(), prYy.size(), pr_vr_groups, pr_vr_extra,
                 static_cast<long long>(std::count_if(prVuse.begin(), prVuse.end(),
                                                      [](long long u) { return u > 0; })),
                 pr_vT_n, pr_vsel_n, pr_vr_sec, pr_vint_sec);
  }

  // Per future end q': phi(arm) = log sum over its entries with y' > q' of
  // Qzbar[q'+1][y'-1] w e^beta1 Sbar[ys+1] (the gap-free future factor), argmax trees over the
  // (xs, x') and (x', xs) orders (k > 0), and W(arm)(ys) = sum_{y' > q'} Qzbar[q'+1][y'-1] w
  // (lazy).  At most PR_AQ_MAX kept, least recently used rebuilt on demand.
  // Argmax tree over a sorted list of positions (ties: lower position, as aq_better).
  struct CTree {
    std::vector<int> pos, t;
    std::vector<double> val;
    int sz = 0;
    int better(int i, int j) const {
      if (i < 0) return j;
      if (j < 0) return i;
      if (val[i] != val[j]) return val[i] > val[j] ? i : j;
      return std::min(i, j);
    }
    void build() {
      sz = 1;
      while (sz < std::max<int>(1, static_cast<int>(pos.size()))) sz <<= 1;
      t.assign(2 * static_cast<size_t>(sz), -1);
      for (int i = 0; i < static_cast<int>(pos.size()); ++i) t[sz + i] = i;
      for (int k = sz - 1; k >= 1; --k) t[k] = better(t[2 * k], t[2 * k + 1]);
    }
    // index (into pos) of the best position in [l, r], -1 if none
    int query(int l, int r) const {
      const int a = static_cast<int>(std::lower_bound(pos.begin(), pos.end(), l) - pos.begin());
      const int b = static_cast<int>(std::upper_bound(pos.begin(), pos.end(), r) - pos.begin()) - 1;
      int best = -1;
      for (int lo = a + sz, hi = b + sz + 1; lo < hi; lo >>= 1, hi >>= 1) {
        if (lo & 1) best = better(best, t[lo++]);
        if (hi & 1) best = better(best, t[--hi]);
      }
      return best;
    }
  };
  struct ArmQ {
    std::vector<double> phi, rs;  // rs: --pr-arms-rank ratio, phi - wu (x' - xs + 1)
    std::vector<int> t1, t2;
    int sz = 0;
    // sep: ra = phi - log Zarm; D[b][x] = log sum ra over arms xs >= x, x' <= b; P1[r] = log
    // sum ra over the arms of r's xs with x' <= x'_r; columns 0..known built (Qz final there)
    std::vector<double> ra, P1;
    std::vector<std::vector<double>> D;
    int known = -1, nfin = 0;
    std::unordered_map<int, std::vector<std::pair<int, double>>> W;
    long long use = 0;
    // gaponly: argmax trees over the indexed arms only (positions of order 1 / 2), built for
    // the gap union `epoch`
    CTree g1, g2;
    long long epoch = -1;
  };
  static constexpr size_t PR_AQ_MAX = 128;
  std::unordered_map<int, ArmQ> prAQ;
  long long pr_use = 0;

  // gaponly: merged x' intervals [l, r] of the gaps (u, v) scored since the last
  // pr_gap_begin (u+1..v-1), and its epoch.
  std::vector<std::pair<int, int>> prGapIv;
  long long gap_epoch = 0;
  void pr_gap_begin(std::vector<std::pair<int, int>>& iv) {
    if (!pr_arm_gaponly) return;
    std::sort(iv.begin(), iv.end());
    prGapIv.clear();
    for (const auto& [l, r] : iv) {
      if (l > r) continue;
      if (!prGapIv.empty() && l <= prGapIv.back().second + 1)
        prGapIv.back().second = std::max(prGapIv.back().second, r);
      else
        prGapIv.emplace_back(l, r);
    }
    ++gap_epoch;
    ++gap_epochs;
  }
  bool pr_gap_covered(int l, int r) const {
    auto it = std::upper_bound(prGapIv.begin(), prGapIv.end(), std::make_pair(l, INT32_MAX));
    return it != prGapIv.begin() && std::prev(it)->second >= r;
  }

  // phi(arm r) at q' as Step A computes it (separable: from the suffix sums prSF; with
  // --pr-partner-span W only the entries q' < y' <= q' + W, summed in y' order).
  // --v-sweep: lK reads the swept b_V entries. mu is the profile weight; the arm
  // base (pair energy, stem weight, pair_bonus) stays the one prYw was built with.
  // A partner the sweep did not keep is -inf. The floor is not applied here.
  void pr_yw_from_sweep() {
    if (!p2_vsweep || !stem_ || stem_->prMu.empty()) return;
    long long rep = 0, miss = 0, ent = 0;
    for (int r = 0; r < prA; ++r) {
      const int xs = prArmXs[r], xp = prArmXp[r];
      for (int k = prArmY0[r]; k < prArmY0[r + 1]; ++k) {
        ++ent;
        const int yp = prYy[k];
        double mu = NEG_INF;
        auto it = stem_->prMu.find(key2(xp, yp));
        if (it != stem_->prMu.end()) {
          auto jt = it->second.find(key2(xs, 0));
          if (jt != it->second.end()) mu = jt->second;
        }
        if (mu < HALF_NEG && p2_v_cut_fb &&
            p2_v_cap > 0 && it != stem_->prMu.end() &&
            static_cast<int>(it->second.size()) >= p2_v_cap) {
          mu = v_cut_fb_mu(xp, yp, xs);
          if (mu > HALF_NEG) ++vfb_n;
        }
        if (mu < HALF_NEG) {
          prYw[k] = NEG_INF;
          ++miss;
          continue;
        }
        prYw[k] = (mu + (E.at(xp, yp) + E.w_stem)) + pair_bonus;
        ++rep;
      }
    }
    prSF.clear();
    std::fprintf(stderr, "[v-sweep-lk] arms %d entries %lld replaced %lld missing %lld\n", prA, ent,
                 rep, miss);
  }
  // lK(arm, q') >= lKx(arm, q') on the weights lK is about to use. One pass, not on the hot path.
  void k_bound_check() {
    pr_sf_build();
    long long checks = 0, bad = 0, bad_inf = 0, bad_fin = 0;
    double worst = 0.0, wub = 0.0, wex = 0.0;
    int wr = -1, wq = -1;
    for (int r = 0; r < prA; ++r) {
      for (int q = 0; q + 1 <= n - 1; ++q) {
        const double ub = lK(r, q);
        const double ex = lKx_sum(r, q);
        ++checks;
        const bool ub_fin = ub > HALF_NEG, ex_fin = ex > HALF_NEG;
        double gap = 0.0;
        if (ex_fin && !ub_fin) gap = 1e300;
        else if (ex_fin && ub_fin) gap = ex - ub;
        if (gap > 1e-8) {
          ++bad;
          if (!ub_fin) ++bad_inf;
          else ++bad_fin;
          if (gap > worst) {
            worst = gap;
            wub = ub;
            wex = ex;
            wr = r;
            wq = q;
          }
        }
      }
    }
    std::fprintf(stderr,
                 "[k-bound] checks %lld violations %lld inf %lld finite %lld worst %.6e ub %.6e ex %.6e arm %d q %d\n",
                 checks, bad, bad_inf, bad_fin, worst, wub, wex, wr, wq);
  }
  void pr_sf_build() {
    if (!pr_aq_sep || !prSF.empty()) return;
    const StemEst& st = *stem_;
    prSF.assign(prYy.size(), NEG_INF);
    for (int r = 0; r < prA; ++r) {
      double acc = NEG_INF;
      for (int k = prArmY0[r + 1] - 1; k >= prArmY0[r]; --k) {
        const double z = st.Qzbar(0, prYy[k] - 1);
        if (z > HALF_NEG && prYw[k] > HALF_NEG) acc = log_add(acc, z + prYw[k]);
        prSF[k] = acc;
      }
    }
  }
  // raw: phi + log Qzbar[0][q'] (separable) or phi (not); ranking only depends on raw.
  double arm_raw(int r, int qp) const {
    const StemEst& st = *stem_;
    const int hi = prArmY0[r + 1];
    int k = static_cast<int>(std::upper_bound(prYy.begin() + prArmY0[r], prYy.begin() + hi, qp) -
                             prYy.begin());
    if (pr_pspan > 0 || !pr_aq_sep) {
      double acc = NEG_INF;
      for (; k < hi && (pr_pspan <= 0 || prYy[k] - qp <= pr_pspan); ++k) {
        const double z = pr_aq_sep ? st.Qzbar(0, prYy[k] - 1) : st.Qzbar(qp + 1, prYy[k] - 1);
        if (z > HALF_NEG && prYw[k] > HALF_NEG) acc = log_add(acc, z + prYw[k]);
      }
      return acc;
    }
    return k < hi ? prSF[k] : NEG_INF;
  }
  double arm_phi(int r, int qp, double z0) const {
    const double raw = arm_raw(r, qp);
    if (!pr_aq_sep) return raw;
    return raw > HALF_NEG && z0 > HALF_NEG ? raw - z0 : NEG_INF;
  }

  // Persistent argmax trees over the arm positions of order 1 (arm id) and order 2 (prByXp),
  // one version per q' (pRoot[t][q']): sweeping q' down from n-1, an arm's leaf is set to its
  // raw value whenever an entry y' = q' + 1 enters (or, with W, y' = q' + W + 1 leaves) its
  // window.  Ties: lower position, as aq_better.  Node 0: the empty tree.
  struct PNode {
    int l, r, best, ver;
    double val;
  };
  std::vector<PNode> pnd;
  std::vector<int> pRoot[2], prXpPos;
  int pver = 0;
  bool pidx_built = false;
  void pidx_pull(int nd) {
    const PNode& L = pnd[pnd[nd].l];
    const PNode& R = pnd[pnd[nd].r];
    int b;
    double v;
    if (L.best < 0 || (R.best >= 0 && R.val > L.val)) {
      b = R.best;
      v = R.val;
    } else {
      b = L.best;
      v = L.val;
    }
    pnd[nd].best = b;
    pnd[nd].val = v;
  }
  int pidx_set(int nd, int lo, int hi, int pos, double val) {
    if (pnd[nd].ver != pver) {
      PNode c = pnd[nd];
      c.ver = pver;
      pnd.push_back(c);
      nd = static_cast<int>(pnd.size()) - 1;
    }
    if (lo == hi) {
      pnd[nd].best = val > HALF_NEG ? pos : -1;
      pnd[nd].val = val > HALF_NEG ? val : NEG_INF;
      return nd;
    }
    const int mid = (lo + hi) / 2;
    if (pos <= mid) {
      const int c = pidx_set(pnd[nd].l, lo, mid, pos, val);
      pnd[nd].l = c;
    } else {
      const int c = pidx_set(pnd[nd].r, mid + 1, hi, pos, val);
      pnd[nd].r = c;
    }
    pidx_pull(nd);
    return nd;
  }
  void pidx_q(int nd, int lo, int hi, int l, int r, int& bp, double& bv) const {
    if (nd == 0 || r < lo || hi < l) return;
    const PNode& x = pnd[nd];
    if (x.best < 0) return;
    if (l <= lo && hi <= r) {
      if (bp < 0 || x.val > bv || (x.val == bv && x.best < bp)) {
        bp = x.best;
        bv = x.val;
      }
      return;
    }
    const int mid = (lo + hi) / 2;
    pidx_q(x.l, lo, mid, l, r, bp, bv);
    pidx_q(x.r, mid + 1, hi, l, r, bp, bv);
  }
  void pidx_build() {
    if (pidx_built) return;
    pidx_built = true;
    SecTimer ts(pidx_sec);
    pr_sf_build();
    pnd.assign(1, PNode{0, 0, -1, -1, NEG_INF});
    for (auto& v : pRoot) v.assign(n, 0);
    prXpPos.assign(prA, 0);
    for (int pos = 0; pos < prA; ++pos) prXpPos[prByXp[pos]] = pos;
    std::vector<std::vector<int>> ev(n);
    for (int r = 0; r < prA; ++r)
      for (int k = prArmY0[r]; k < prArmY0[r + 1]; ++k) {
        if (prYy[k] - 1 >= 0) ev[prYy[k] - 1].push_back(r);
        if (pr_pspan > 0 && prYy[k] - pr_pspan - 1 >= 0) ev[prYy[k] - pr_pspan - 1].push_back(r);
      }
    std::vector<int> stamp(prA, -1);
    int root[2] = {0, 0};
    for (int qp = n - 1; qp >= 0; --qp) {
      ++pver;
      for (int r : ev[qp]) {
        if (stamp[r] == qp) continue;
        stamp[r] = qp;
        const double raw = arm_raw(r, qp);
        root[0] = pidx_set(root[0], 0, prA - 1, r, raw);
        root[1] = pidx_set(root[1], 0, prA - 1, prXpPos[r], raw);
        ++pidx_upd;
      }
      std::vector<int>().swap(ev[qp]);
      pRoot[0][qp] = root[0];
      pRoot[1][qp] = root[1];
    }
    pidx_nodes = static_cast<long long>(pnd.size());
    pidx_mb = (static_cast<double>(pnd.capacity()) * sizeof(PNode) +
               static_cast<double>(2 * n + prA) * sizeof(int)) / 1048576.0;
  }

  double aq_val(const ArmQ& a, int ord, int pos) const {
    return (pr_ratio ? a.rs : a.phi)[ord == 1 ? pos : prByXp[pos]];
  }
  int aq_better(const ArmQ& a, int ord, int p1, int p2) const {
    if (p1 < 0) return p2;
    if (p2 < 0) return p1;
    const double v1 = aq_val(a, ord, p1), v2 = aq_val(a, ord, p2);
    if (v1 != v2) return v1 > v2 ? p1 : p2;
    return std::min(p1, p2);
  }
  void aq_tree(ArmQ& a, int ord, std::vector<int>& t) const {
    t.assign(2 * static_cast<size_t>(a.sz), -1);
    for (int p = 0; p < prA; ++p) t[a.sz + p] = p;
    for (int k = a.sz - 1; k >= 1; --k) t[k] = aq_better(a, ord, t[2 * k], t[2 * k + 1]);
  }
  int aq_query(const ArmQ& a, int ord, int l, int r) const {
    const std::vector<int>& t = ord == 1 ? a.t1 : a.t2;
    int best = -1;
    for (int lo = l + a.sz, hi = r + a.sz + 1; lo < hi; lo >>= 1, hi >>= 1) {
      if (lo & 1) best = aq_better(a, ord, best, t[lo++]);
      if (hi & 1) best = aq_better(a, ord, best, t[--hi]);
    }
    return best;
  }

  ArmQ& pr_armq(int qp) {
    auto it = prAQ.find(qp);
    if (it != prAQ.end()) {
      if (!pr_arm_gaponly || it->second.epoch == gap_epoch) {
        it->second.use = ++pr_use;
        return it->second;
      }
      prAQ.erase(it);
    }
    if (prAQ.size() >= pr_aq_max) {
      auto old = prAQ.begin();
      for (auto jt = prAQ.begin(); jt != prAQ.end(); ++jt)
        if (jt->second.use < old->second.use) old = jt;
      prAQ.erase(old);
    }
    const auto t0 = std::chrono::steady_clock::now();
    const StemEst& st = *stem_;
    ArmQ a;
    if (pr_arm_persist) {
      pidx_build();
      a.use = ++pr_use;
      ++pr_aq_built;
      pr_aq_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      return prAQ.emplace(qp, std::move(a)).first->second;
    }
    pr_sf_build();
    const double z0 = pr_aq_sep ? st.Qzbar(0, qp) : 0.0;
    if (pr_arm_gaponly) {
      a.epoch = gap_epoch;
      for (const auto& [l, r] : prGapIv)
        for (int pos = prXpFirst[l]; pos < prXpFirst[std::min(r, n - 1) + 1]; ++pos) {
          a.g2.pos.push_back(pos);
          a.g1.pos.push_back(prByXp[pos]);
        }
      std::sort(a.g1.pos.begin(), a.g1.pos.end());
      for (int p : a.g1.pos) a.g1.val.push_back(arm_phi(p, qp, z0));
      for (int p : a.g2.pos) a.g2.val.push_back(arm_phi(prByXp[p], qp, z0));
      a.g1.build();
      a.g2.build();
      aq_arms += static_cast<long long>(a.g2.pos.size());
      a.use = ++pr_use;
      ++pr_aq_built;
      pr_aq_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      return prAQ.emplace(qp, std::move(a)).first->second;
    }
    a.phi.assign(prA, NEG_INF);
    aq_arms += prA;
    for (int r = 0; r < prA; ++r) {
      if (pr_pspan > 0) {
        a.phi[r] = arm_phi(r, qp, z0);
        continue;
      }
      const int hi = prArmY0[r + 1];
      int k = static_cast<int>(std::upper_bound(prYy.begin() + prArmY0[r], prYy.begin() + hi, qp) -
                               prYy.begin());
      if (pr_aq_sep) {
        if (k < hi && prSF[k] > HALF_NEG && z0 > HALF_NEG) a.phi[r] = prSF[k] - z0;
        continue;
      }
      double acc = NEG_INF;
      for (; k < hi; ++k) {
        const double z = st.Qzbar(qp + 1, prYy[k] - 1);
        if (z > HALF_NEG) acc = log_add(acc, z + prYw[k]);
      }
      a.phi[r] = acc;
    }
    if (pr_arms_mode == 1)
      for (int r = 0; r < prA; ++r) a.nfin += a.phi[r] > HALF_NEG;
    if (pr_ratio) {
      a.rs.resize(prA);
      for (int r = 0; r < prA; ++r)
        a.rs[r] = a.phi[r] > HALF_NEG ? a.phi[r] - E.wu * (prArmXp[r] - prArmXs[r] + 1) : NEG_INF;
    }
    if (pr_arms_k > 0 || pr_ratio) {
      a.sz = 1;
      while (a.sz < std::max(1, prA)) a.sz <<= 1;
      aq_tree(a, 1, a.t1);
      aq_tree(a, 2, a.t2);
    }
    a.use = ++pr_use;
    ++pr_aq_built;
    pr_aq_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return prAQ.emplace(qp, std::move(a)).first->second;
  }

  // Arms of gap (u, v) at q' in selection order: xs >= u+1 (two: xs = u+1), x' <= v-1, finite
  // phi; the top pr_arms_k by phi (ties: lower position), all if 0; ratio: by rs until the
  // adaptive stop (set_stem_pr_select).
  void pr_select(const ArmQ& a, int qp, int u, int v, bool two, std::vector<int>& out) {
    out.clear();
    if (u + 1 > v - 1) return;
    if (pr_arm_gaponly && !pr_gap_covered(u + 1, v - 1)) ++gap_miss;
    int lo, hi, ord;
    if (two) {
      lo = prXsFirst[u + 1];
      const int end = prXsFirst[u + 2];
      hi = static_cast<int>(std::upper_bound(prArmXp.begin() + lo, prArmXp.begin() + end, v - 1) -
                            prArmXp.begin());
      ord = 1;
    } else {
      lo = prXpFirst[u + 1];
      hi = prXpFirst[v];
      ord = 2;
    }
    auto arm_of = [&](int pos) { return ord == 1 ? pos : prByXp[pos]; };
    auto ok = [&](int r) { return two || prArmXs[r] >= u + 1; };
    if (pr_arms_k == 0 && !pr_ratio) {
      for (int pos = lo; pos < hi; ++pos) {
        const int r = arm_of(pos);
        if (ok(r) && a.phi[r] > HALF_NEG) out.push_back(r);
      }
      return;
    }
    // ratio: the gap's total r over its arms (scaled by the largest), for the adaptive stop
    double smax = NEG_INF, tot = 0.0, sel = 0.0;
    if (pr_ratio) {
      for (int pos = lo; pos < hi; ++pos) {
        const int r = arm_of(pos);
        if (ok(r)) smax = std::max(smax, a.rs[r]);
      }
      if (smax < HALF_NEG) return;
      for (int pos = lo; pos < hi; ++pos) {
        const int r = arm_of(pos);
        if (ok(r) && a.rs[r] > HALF_NEG) tot += std::exp(a.rs[r] - smax);
      }
    }
    const int n_max = pr_ratio ? pr_arms_max : pr_arms_k;
    struct HE {
      double v;
      int pos, l, r;
      bool operator<(const HE& o) const { return v != o.v ? v < o.v : pos > o.pos; }
    };
    std::priority_queue<HE> pq;
    double z0 = 0.0;
    int root = 0;
    if (pr_arm_persist) {
      z0 = stem_->Qzbar(0, qp);
      if (z0 < HALF_NEG) return;
      root = pRoot[ord - 1][qp];
    }
    auto push = [&](int l, int r) {
      if (l > r) return;
      int p;
      double val;
      if (pr_arm_persist) {
        p = -1;
        double raw = NEG_INF;
        pidx_q(root, 0, prA - 1, l, r, p, raw);
        if (p < 0) return;
        val = raw - z0;
      } else if (pr_arm_gaponly) {
        const CTree& c = ord == 1 ? a.g1 : a.g2;
        const int i = c.query(l, r);
        if (i < 0) return;
        p = c.pos[i];
        val = c.val[i];
      } else {
        p = aq_query(a, ord, l, r);
        if (p < 0) return;
        val = aq_val(a, ord, p);
      }
      if (val < HALF_NEG) return;
      pq.push({val, p, l, r});
    };
    push(lo, hi - 1);
    while (!pq.empty() && static_cast<int>(out.size()) < n_max) {
      const HE t = pq.top();
      pq.pop();
      const int r = arm_of(t.pos);
      if (ok(r)) {
        out.push_back(r);
        if (pr_ratio) {
          sel += std::exp(t.v - smax);
          if (tot - sel < pr_arms_eps * sel) break;
        }
      } else {
        ++pr_T_skip;
      }
      push(t.l, t.pos - 1);
      push(t.pos + 1, t.r);
    }
  }

  // sep: columns known+1..b of a's tables (arms with x' = c added to every D[c][x], x <= xs).
  void pr_sep_extend(ArmQ& a, int b) {
    const auto t0 = std::chrono::steady_clock::now();
    if (a.ra.empty()) {
      a.ra.assign(prA, NEG_INF);
      a.P1.assign(prA, NEG_INF);
    }
    std::vector<double> col;
    for (int c = a.known + 1; c <= b; ++c) {
      col.assign(c + 1, NEG_INF);
      for (int pos = prXpFirst[c]; pos < prXpFirst[c + 1]; ++pos) {
        const int r = prByXp[pos], xs = prArmXs[r];
        if (a.phi[r] > HALF_NEG) {
          const double z = pr_arm_z == 1 ? E.wu * (c - xs + 1) : pr_qz(xs, c);
          if (z > HALF_NEG) a.ra[r] = a.phi[r] - z;
        }
        col[xs] = log_add(col[xs], a.ra[r]);
        const bool cont = r > 0 && prArmXs[r - 1] == xs;
        a.P1[r] = log_add(cont ? a.P1[r - 1] : NEG_INF, a.ra[r]);
      }
      std::vector<double> d(c + 1, NEG_INF);
      double suf = NEG_INF;
      for (int x = c; x >= 0; --x) {
        suf = log_add(suf, col[x]);
        d[x] = log_add(x <= c - 1 ? a.D[c - 1][x] : NEG_INF, suf);
      }
      a.D.push_back(std::move(d));
      ++pr_sep_cols;
    }
    a.known = std::max(a.known, b);
    pr_sep_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }

  std::unordered_set<int> prSepTouched;
  std::unordered_set<uint64_t> pscSeen;

  // sep: log of the exterior arm term without Ctx's prefix, Qz[u+1][v-1] S_q'(u, v).
  double pr_sep_ext(int qp, int u, int v, bool two) {
    if (u + 1 > v - 1 || qp + 1 > n - 1) return NEG_INF;
    ArmQ& a = pr_armq(qp);
    if (prSepTouched.insert(qp).second) {
      ++pr_sep_groups;
      pr_sep_garms += a.nfin;
    }
    ++pr_sep_q;
    if (a.known < v - 1) pr_sep_extend(a, v - 1);
    double S = NEG_INF;
    if (two) {
      const int lo = prXsFirst[u + 1], end = prXsFirst[u + 2];
      const int hi = static_cast<int>(
          std::upper_bound(prArmXp.begin() + lo, prArmXp.begin() + end, v - 1) - prArmXp.begin());
      if (hi > lo) S = a.P1[hi - 1];
    } else {
      S = a.D[v - 1][u + 1];
    }
    const double g = pr_qz(u + 1, v - 1);
    const double se = S > HALF_NEG && g > HALF_NEG ? g + S : NEG_INF;
    if (pr_sep_check) pr_sep_diag(a, qp, u, v, two, S, se);
    return se;
  }

  // sep_check: S against the brute-force sum of ra, se against the exact placement over all
  // arms; the per-term log ratio (Qz[u+1][v-1] / Zarm) / (Qz[u+1][xs-1] Qz[x'+1][v-1]).
  void pr_sep_diag(const ArmQ& a, int qp, int u, int v, bool two, double S, double se) {
    const uint64_t key = (static_cast<uint64_t>(qp) << 43) | (static_cast<uint64_t>(u) << 22) |
                         (static_cast<uint64_t>(v) << 1) | (two ? 1u : 0u);
    if (!pscSeen.insert(key).second) return;
    ++psc_q;
    int lo, hi;
    if (two) {
      lo = prXsFirst[u + 1];
      hi = prXsFirst[u + 2];
    } else {
      lo = prXpFirst[u + 1];
      hi = prXpFirst[v];
    }
    const double g = pr_qz(u + 1, v - 1);
    std::vector<std::pair<double, double>> tm;  // (exact term, log ratio)
    double bs = NEG_INF, ex = NEG_INF;
    for (int pos = lo; pos < hi; ++pos) {
      const int r = two ? pos : prByXp[pos];
      if (prArmXs[r] < u + 1 || prArmXp[r] > v - 1 || a.phi[r] < HALF_NEG) continue;
      bs = log_add(bs, a.ra[r]);
      const double gl = two ? 0.0 : pr_qz(u + 1, prArmXs[r] - 1);
      const double gr = pr_qz(prArmXp[r] + 1, v - 1);
      if (gl < HALF_NEG || gr < HALF_NEG) continue;
      const double t = (gl + gr) + a.phi[r];
      ex = log_add(ex, t);
      if (a.ra[r] > HALF_NEG && g > HALF_NEG) tm.emplace_back(t, (g + (a.ra[r] - a.phi[r])) - (gl + gr));
    }
    if ((S > HALF_NEG) != (bs > HALF_NEG)) ++psc_tab_mis;
    else if (S > HALF_NEG) psc_tab_rel = std::max(psc_tab_rel, std::fabs(std::expm1(S - bs)));
    if ((se > HALF_NEG) != (ex > HALF_NEG)) {
      ++psc_sep_mis;
    } else if (se > HALF_NEG) {
      const double d = se - ex;
      psc_sep_max = std::max(psc_sep_max, std::fabs(d));
      psc_d_sum += d;
      psc_d_abs += std::fabs(d);
      psc_hq[psc_bin(d)] += 1.0;
      for (const auto& [t, lr] : tm) {
        const double w = std::exp(t - ex);
        psc_ht[psc_bin(lr)] += w;
        psc_lr_w += w * lr;
        psc_lr_wabs += w * std::fabs(lr);
        ++psc_terms;
      }
    }
  }

  const std::vector<std::pair<int, double>>& pr_W(ArmQ& a, int r, int qp) {
    auto it = a.W.find(r);
    if (it != a.W.end()) return it->second;
    if (pr_vxs) throw std::runtime_error("--pr-v-key xs: pr_W needs per-ys arm entries");
    SecTimer tw(pr_sec_W);
    ++pr_W_built;
    const StemEst& st = *stem_;
    std::vector<std::pair<int, double>> w;
    const int hi = prArmE0[r + 1];
    for (int k = static_cast<int>(std::upper_bound(prEy.begin() + prArmE0[r], prEy.begin() + hi, qp) -
                                  prEy.begin());
         k < hi; ++k) {
      if (pr_pspan > 0 && prEy[k] - qp > pr_pspan) break;
      const double z = st.Qzbar(qp + 1, prEy[k] - 1);
      if (z > HALF_NEG) w.emplace_back(prEys[k], z + prEw[k]);
      ++pr_W_ent;
    }
    std::stable_sort(w.begin(), w.end(), [](const std::pair<int, double>& x,
                                            const std::pair<int, double>& y) {
      return x.first < y.first;
    });
    std::vector<std::pair<int, double>> m;
    for (const auto& e : w) {
      if (!m.empty() && m.back().first == e.first) m.back().second = log_add(m.back().second, e.second);
      else m.push_back(e);
    }
    return a.W.emplace(r, std::move(m)).first->second;
  }

  // Dense accumulator over [0, n) (stamped).
  std::vector<double> prAcc;
  std::vector<int> prStamp, prTouched;
  int pr_stamp = 0;
  void acc_begin() {
    if (prAcc.size() != static_cast<size_t>(n)) {
      prAcc.assign(n, NEG_INF);
      prStamp.assign(n, 0);
    }
    ++pr_stamp;
    prTouched.clear();
  }
  void acc_add(int i, double w) {
    if (prStamp[i] != pr_stamp) {
      prStamp[i] = pr_stamp;
      prAcc[i] = NEG_INF;
      prTouched.push_back(i);
    }
    prAcc[i] = log_add(prAcc[i], w);
  }

  // T(q'; u, v)(ys) = sum over the selected arms (xs, x') of Qz[u+1][xs-1] Qz[x'+1][v-1]
  // W(arm)(ys) (two: xs = u+1, no left filler); ext = log sum_ys T e^beta1 Sbar[ys+1].
  // Cached by q' until the scan passes q' (no term depends on the step).
  struct PrTV {
    std::vector<int> ys;
    std::vector<double> w;
    double ext = NEG_INF;
  };
  std::unordered_map<int, std::unordered_map<uint64_t, PrTV>> prTc;
  std::vector<int> prSel;

  const PrTV& pr_T(int qp, int u, int v, bool two) {
    std::unordered_map<uint64_t, PrTV>& bucket = prTc[qp];
    const uint64_t key = key3(u, v, two ? 1 : 0);
    auto it = bucket.find(key);
    if (it != bucket.end()) {
      ++pr_T_hits;
      return it->second;
    }
    ++pr_T_built;
    PrPartnerTimer pt(*this);
    SecTimer tt(pr_sec_T);
    const StemEst& st = *stem_;
    PrTV t;
    if (qp + 1 <= n - 1) {
      ArmQ& a = pr_armq(qp);
      const bool xs = pr_abl & (PA_XS | PA_XB);
      if (pr_abl & PA_XS) {
        const int k0 = pr_arms_k;
        const bool r0 = pr_ratio;
        pr_arms_k = 0;
        pr_ratio = false;
        pr_select(a, qp, u, v, two, prSel);
        pr_arms_k = k0;
        pr_ratio = r0;
      } else {
        pr_select(a, qp, u, v, two, prSel);
      }
      pr_T_arms += static_cast<long long>(prSel.size());
      if (!prVrg.empty()) pr_vdiag_T(qp, u, v, two, prSel);
      const double gz = xs ? pr_qz(u + 1, v - 1) : 0.0;
      acc_begin();
      for (int r : prSel) {
        double g;
        if (xs) {
          const double za = pr_qz(prArmXs[r], prArmXp[r]);
          if (gz < HALF_NEG || za < HALF_NEG) continue;
          g = gz - za;
        } else {
          const double gl = two ? 0.0 : pr_qz(u + 1, prArmXs[r] - 1);
          const double gr = pr_qz(prArmXp[r] + 1, v - 1);
          if (gl < HALF_NEG || gr < HALF_NEG) continue;
          g = gl + gr;
        }
        for (const auto& [ys, w] : pr_W(a, r, qp)) acc_add(ys, g + w);
      }
      std::sort(prTouched.begin(), prTouched.end());
      for (int ys : prTouched) {
        t.ys.push_back(ys);
        t.w.push_back(prAcc[ys]);
        t.ext = log_add(t.ext, prAcc[ys] + (E.w_beta1 + st.prS[ys + 1]));
      }
      pr_T_ys += static_cast<long long>(t.ys.size());
    }
    return bucket.emplace(key, std::move(t)).first->second;
  }

  // K(p'; u, v)(i) = sum over the surviving L-stems Qgl[i][p'-1][e][f], u+1 <= e <= f <= v-1,
  // of Qz[u+1][e-1] alpha at(i, f) Qz[f+1][v-1]; ext = log sum_i K C[i-1].  Final once the
  // scan passes v - 1, so kept across steps (cleared when it outgrows PR_KC_MAX).
  struct PrKV {
    std::vector<int> is;
    std::vector<double> w;
    double ext = NEG_INF;
  };
  static constexpr size_t PR_KC_MAX = size_t(1) << 21;
  std::unordered_map<uint64_t, PrKV> prKc;

  const PrKV& pr_K(int pp, int u, int v) {
    const uint64_t key = key3(pp, u, v);
    auto it = prKc.find(key);
    if (it != prKc.end()) {
      ++pr_K_hits;
      return it->second;
    }
    if (prKc.size() >= PR_KC_MAX) prKc.clear();
    ++pr_K_built;
    PrPartnerTimer pt(*this);
    SecTimer tk(pr_sec_K);
    PrKV kv;
    auto ie = stE.find(pp);
    if (ie != stE.end()) {
      acc_begin();
      const std::vector<int>& es = ie->second;
      for (auto ke = std::lower_bound(es.begin(), es.end(), u + 1); ke != es.end() && *ke <= v - 1;
           ++ke) {
        const int e = *ke;
        const double gl = pr_qz(u + 1, e - 1);
        if (gl < HALF_NEG) continue;
        auto ia = stA.find(key2(pp, e));
        if (ia == stA.end()) continue;
        const StemA& sa = ia->second;
        for (size_t kf = 0; kf < sa.fs.size(); ++kf) {
          const int f = sa.fs[kf];
          if (f > v - 1) break;
          const double gr = pr_qz(f + 1, v - 1);
          if (gr < HALF_NEG) continue;
          const double g = (pr_abl & PA_KS) ? pr_qz(u + 1, v - 1) - pr_qz(e, f) : gl + gr;
        if (!(g > HALF_NEG)) continue;
        for (const auto& [i, lw] : sa.il[kf]) acc_add(i, g + lw);
        }
      }
      std::sort(prTouched.begin(), prTouched.end());
      for (int i : prTouched) {
        kv.is.push_back(i);
        kv.w.push_back(prAcc[i]);
        kv.ext = log_add(kv.ext, prAcc[i] + prefix_log(i));
      }
      pr_K_is += static_cast<long long>(kv.is.size());
    }
    return prKc.emplace(key, std::move(kv)).first->second;
  }

  double pr_ctx(int i, int j) const {
    const StemEst& st = *stem_;
    return log_add(E.w_beta1 + prefix_log(i) + st.prS[j + 1], E.w_beta1m + st.Bml(i, j));
  }

  // Stop_L + Stop_R of an outer pair p'.q' for a stem with inner pair u.v, without that
  // stem's own beta2 + at(u, v) and the growth factor:
  //   Stop_L = sum_ys T(q'; u, v)(ys) Ctx(p', ys),  Stop_R = sum_i K(p'; u, v)(i) Ctx(i, q').
  // Final once built: cached by (u, v), then (p', q'); cleared when it outgrows PR_SC_MAX.
  static constexpr size_t PR_SC_MAX = size_t(1) << 22;
  std::unordered_map<uint64_t, std::unordered_map<uint64_t, double>> prSc;
  size_t pr_sc_n = 0;

  // --pr-abl q: T(q; u, v) / Qz[q+1][q'] for T(q'; u, v); p: K(p; u, v) for K(p'; u, v)
  // (p, q: the scored state's); bml: Ctx without its multiloop term.
  double pr_SLR_abl(int pp, int qp, int u, int v) {
    const StemEst& st = *stem_;
    const bool bml = !(pr_abl & PA_BML);
    const int tq = (pr_abl & PA_Q) && pr_abl_q >= 0 ? pr_abl_q : qp;
    const double tz = -pr_abl_qz(tq, qp);
    double sl = NEG_INF;
    if (tz > HALF_NEG && -tz > HALF_NEG) {
      const PrTV& t = pr_T(tq, u, v, false);
      if (!t.ys.empty()) {
        double bm = NEG_INF;
        if (bml)
          for (size_t k = 0; k < t.ys.size(); ++k) bm = log_add(bm, t.w[k] + st.Bml(pp, t.ys[k]));
        sl = tz + log_add(prefix_log(pp) + t.ext, E.w_beta1m + bm);
      }
    }
    const int kp = (pr_abl & PA_P) && pr_abl_p >= 0 ? pr_abl_p : pp;
    double sr = NEG_INF;
    const PrKV& kv = pr_K(kp, u, v);
    if (!kv.is.empty()) {
      double bm = NEG_INF;
      if (bml)
        for (size_t k = 0; k < kv.is.size(); ++k) bm = log_add(bm, kv.w[k] + st.Bml(kv.is[k], qp));
      sr = log_add((E.w_beta1 + st.prS[qp + 1]) + kv.ext, E.w_beta1m + bm);
    }
    return log_add(sl, sr);
  }

  double pr_SLR(std::unordered_map<uint64_t, double>& b, int pp, int qp, int u, int v) {
    if (pr_abl & (PA_Q | PA_P)) {
      ++pr_slr_miss;
      PrPartnerTimer pt(*this);
      SecTimer ts(pr_sec_SLR);
      return pr_SLR_abl(pp, qp, u, v);
    }
    const uint64_t key = key2(pp, qp);
    auto it = b.find(key);
    if (it != b.end()) {
      ++pr_slr_hits;
      return it->second;
    }
    ++pr_slr_miss;
    PrPartnerTimer pt(*this);
    SecTimer ts(pr_sec_SLR);
    if (pr_abl & PA_BML) {
      const double s = pr_SLR_abl(pp, qp, u, v);
      b.emplace(key, s);
      ++pr_sc_n;
      return s;
    }
    const StemEst& st = *stem_;
    double sl = NEG_INF;
    if (pr_arms_mode == 1) {
      const double se = pr_sep_ext(qp, u, v, false);
      const PrTV& t = pr_T(qp, u, v, false);
      double bm = NEG_INF;
      for (size_t k = 0; k < t.ys.size(); ++k) bm = log_add(bm, t.w[k] + st.Bml(pp, t.ys[k]));
      sl = log_add(prefix_log(pp) + se, E.w_beta1m + bm);
    } else {
      const PrTV& t = pr_T(qp, u, v, false);
      if (!t.ys.empty()) {
        double bm = NEG_INF;
        for (size_t k = 0; k < t.ys.size(); ++k) bm = log_add(bm, t.w[k] + st.Bml(pp, t.ys[k]));
        sl = log_add(prefix_log(pp) + t.ext, E.w_beta1m + bm);
      }
    }
    double sr = NEG_INF;
    const PrKV& kv = pr_K(pp, u, v);
    if (!kv.is.empty()) {
      double bm = NEG_INF;
      for (size_t k = 0; k < kv.is.size(); ++k) bm = log_add(bm, kv.w[k] + st.Bml(kv.is[k], qp));
      sr = log_add((E.w_beta1 + st.prS[qp + 1]) + kv.ext, E.w_beta1m + bm);
    }
    const double s = log_add(sl, sr);
    b.emplace(key, s);
    ++pr_sc_n;
    return s;
  }

  // --pr-grow-select bound.  At step q every gap lies in [0, q-1], every arm (xs, x') with
  // x' <= q-1 has its final Qz[xs][x'], and every partner entry with f <= q-1 is in stA.  An
  // arm inside the gap has Qz[u+1][xs-1] Qz[x'+1][v-1] <= Qz[u+1][v-1] / Qz[xs][x'] (a filler
  // of [u+1, v-1] contains every concatenation), so for every gap of the step
  //   Stop_L(p', q'; u, v) <= Qz[u+1][v-1] U_L(p', q'),
  //     U_L = sum over those arms and their ys of W(arm)(ys) Ctx(p', ys) / Qz[xs][x'],
  //   Stop_R(p', q'; u, v) <= Qz[u+1][v-1] U_R(p', q'),
  //     U_R = sum over those entries (p', e, f; i) of alpha at(i, f) Ctx(i, q') / Qz[e][f]
  // (the selected arms are a subset of all arms).  With dangles a cut b | b+1 can split an
  // unpaired run, adding d3(b) + d5(b+1), so each cut of the arm (xs - 1 | xs, x' | x' + 1)
  // or partner entry multiplies its 1 / Qz by e^kappa(b), kappa(b) = max(0, d3w[b] + d5w[b+1]).
  // Per q': the ys-vector of sum W / Qz, columns added as the scan passes them; per p': the
  // i-vector of sum alpha at / Qz.  The Bml parts of Ctx are dot products with exp(Bml) scaled
  // per row (U_L) or per column (U_R), in long double (the product of the two scalings stays
  // in range).
  struct PrPend {
    int d, e, f, i;
    double w;
  };
  std::vector<PrPend> prPend;
  size_t prPendDone = 0;
  int pr_q = -1, pr_sync_q = -1;
  struct PrUL {
    int known = -1;
    bool dirty = false;
    std::vector<double> w;
    std::vector<long double> lin;
    double phi = NEG_INF, lmax = NEG_INF;
    double acc = NEG_INF;  // --pr-v-key xs: phi accumulated directly (w / lin unused)
  };
  std::unordered_map<int, PrUL> prUL;
  struct PrUR {
    bool dirty = false;
    std::vector<double> w;
    std::vector<long double> lin;
    double r1 = NEG_INF, lmax = NEG_INF;
  };
  std::vector<PrUR> prUR;
  std::vector<long double> prEBr, prEBc;
  std::vector<double> prMBr, prMBc, prKap;
  std::unordered_map<uint64_t, std::pair<double, double>> prUc;
  // log 1 / Qz[a][b] with the cut factors of both ends (a >= 1)
  double pr_iqz(int a, int b) const {
    const double z = pr_qz(a, b);
    if (z < HALF_NEG) return NEG_INF;
    return (prKap[a - 1] + prKap[b]) - z;
  }

  void pr_ub_init() {
    if (!prMBr.empty()) return;
    const StemEst& st = *stem_;
    const size_t nn = static_cast<size_t>(n) * n;
    prKap.assign(n + 1, 0.0);
    if (dz)
      for (int b = 0; b + 1 < n; ++b) prKap[b] = std::max(0.0, E.d3w[b] + E.d5w[b + 1]);
    if (g_U_nobml || g_grow_rec) {
      prMBr.assign(n, NEG_INF);
      prMBc.assign(n, NEG_INF);
      prUR.assign(n, PrUR());
      return;
    }
    prEBr.assign(nn, 0.0L);
    prEBc.assign(nn, 0.0L);
    prMBr.assign(n, NEG_INF);
    prMBc.assign(n, NEG_INF);
    for (int d = 0; d < n; ++d)
      for (int e = 0; e < n; ++e) {
        const double b = st.Bml(d, e);
        prMBr[d] = std::max(prMBr[d], b);
        prMBc[e] = std::max(prMBc[e], b);
      }
    for (int d = 0; d < n; ++d)
      for (int e = 0; e < n; ++e) {
        const double b = st.Bml(d, e);
        if (b < HALF_NEG) continue;
        prEBr[ix(d, e)] = std::exp(static_cast<long double>(b) - prMBr[d]);
        prEBc[ix(e, d)] = std::exp(static_cast<long double>(b) - prMBc[e]);
      }
    prUR.assign(n, PrUR());
  }

  // Partner entries with f <= q-1 into their p' vectors (entries arrive in f order).
  void pr_ub_sync() {
    if (pr_sync_q == pr_q) return;
    pr_sync_q = pr_q;
    while (prPendDone < prPend.size() && prPend[prPendDone].f <= pr_q - 1) {
      const PrPend& e = prPend[prPendDone++];
      const double iz = e.e >= 1 ? pr_iqz(e.e, e.f) : NEG_INF;
      if (iz < HALF_NEG) continue;
      PrUR& r = prUR[e.d];
      if (r.w.empty()) r.w.assign(n, NEG_INF);
      r.w[e.i] = log_add(r.w[e.i], e.w + iz);
      r.dirty = true;
      ++gs_rent;
    }
  }

  PrUL& pr_ul(int qp) {
    const StemEst& st = *stem_;
    PrUL& u = prUL[qp];
    if (u.w.empty() && !pr_vxs) u.w.assign(n, NEG_INF);
    for (int c = u.known + 1; c <= pr_q - 1; ++c) {
      for (int pos = prXpFirst[c]; pos < prXpFirst[c + 1]; ++pos) {
        const int r = prByXp[pos];
        const double iz = prArmXs[r] >= 1 ? pr_iqz(prArmXs[r], c) : NEG_INF;
        if (iz < HALF_NEG) continue;
        if (pr_vxs) {
          const int hy = prArmY0[r + 1];
          for (int k = static_cast<int>(std::upper_bound(prYy.begin() + prArmY0[r], prYy.begin() + hy, qp) -
                                        prYy.begin());
               k < hy; ++k) {
            if (pr_pspan > 0 && prYy[k] - qp > pr_pspan) break;
            const double zb = st.Qzbar(qp + 1, prYy[k] - 1);
            if (zb < HALF_NEG) continue;
            u.acc = log_add(u.acc, (zb + prYw[k]) + iz);
            u.dirty = true;
          }
          continue;
        }
        const int hi = prArmE0[r + 1];
        for (int k = static_cast<int>(std::upper_bound(prEy.begin() + prArmE0[r], prEy.begin() + hi, qp) -
                                      prEy.begin());
             k < hi; ++k) {
          if (pr_pspan > 0 && prEy[k] - qp > pr_pspan) break;
          const double zb = st.Qzbar(qp + 1, prEy[k] - 1);
          if (zb < HALF_NEG) continue;
          double& t = u.w[prEys[k]];
          t = log_add(t, (zb + prEw[k]) + iz);
          u.dirty = true;
        }
      }
      ++gs_cols;
    }
    u.known = std::max(u.known, pr_q - 1);
    if (u.dirty && pr_vxs) {
      u.dirty = false;
      u.phi = u.lmax = u.acc;
    } else if (u.dirty) {
      u.dirty = false;
      u.phi = NEG_INF;
      u.lmax = NEG_INF;
      for (int y = qp + 1; y < n; ++y)
        if (u.w[y] > HALF_NEG) {
          u.phi = log_add(u.phi, u.w[y] + (E.w_beta1 + st.prS[y + 1]));
          u.lmax = std::max(u.lmax, u.w[y]);
        }
      u.lin.assign(n, 0.0L);
      for (int y = qp + 1; y < n; ++y)
        if (u.w[y] > HALF_NEG) u.lin[y] = std::exp(static_cast<long double>(u.w[y]) - u.lmax);
    }
    return u;
  }

  PrUR& pr_ur(int d) {
    PrUR& r = prUR[d];
    if (r.dirty) {
      r.dirty = false;
      r.r1 = NEG_INF;
      r.lmax = NEG_INF;
      for (int i = 0; i < d; ++i)
        if (r.w[i] > HALF_NEG) {
          r.r1 = log_add(r.r1, r.w[i] + prefix_log(i));
          r.lmax = std::max(r.lmax, r.w[i]);
        }
      r.lin.assign(n, 0.0L);
      for (int i = 0; i < d; ++i)
        if (r.w[i] > HALF_NEG) r.lin[i] = std::exp(static_cast<long double>(r.w[i]) - r.lmax);
    }
    return r;
  }

  // (log U_L, log U_R) of the endpoint p'.q' at the current step (cached per step).
  std::pair<double, double> pr_U(int pp, int qp) {
    const uint64_t key = key2(pp, qp);
    auto it = prUc.find(key);
    if (it != prUc.end()) return it->second;
    if (pr_vxs) throw std::runtime_error("--pr-v-key xs: pr_U needs per-ys arm entries");
    pr_ub_init();
    pr_ub_sync();
    ++gs_U;
    const StemEst& st = *stem_;
    double ul = NEG_INF, ur = NEG_INF;
    if (qp + 1 <= n - 1) {
      const PrUL& u = pr_ul(qp);
      if (u.lmax > HALF_NEG) {
        ul = prefix_log(pp) + u.phi;
        if (!g_U_nobml) {
          long double dot = 0.0L;
          const long double* eb = &prEBr[ix(pp, 0)];
          for (int y = qp + 1; y < n; ++y) dot += u.lin[y] * eb[y];
          if (dot > 0.0L)
            ul = log_add(ul, (E.w_beta1m + u.lmax) + (prMBr[pp] + static_cast<double>(std::log(dot))));
        }
      }
    }
    if (pp >= 1) {
      const PrUR& r = pr_ur(pp);
      if (r.lmax > HALF_NEG) {
        ur = (E.w_beta1 + st.prS[qp + 1]) + r.r1;
        if (!g_U_nobml) {
          long double dot = 0.0L;
          const long double* eb = &prEBc[ix(qp, 0)];
          for (int i = 0; i < pp; ++i) dot += r.lin[i] * eb[i];
          if (dot > 0.0L)
            ur = log_add(ur, (E.w_beta1m + r.lmax) + (prMBc[qp] + static_cast<double>(std::log(dot))));
        }
      }
    }
    return prUc.emplace(key, std::make_pair(ul, ur)).first->second;
  }

  // Keeps the entries of `all` (key order) with the largest iota = w (U_L + U_R) until their
  // sum reaches (1 - eps) of the total, at most B (0: no cap); a zero total keeps the B
  // largest by w.  kind 0: growth node, 1: evaluated root / key.
  void pr_gselect(PProf& all, PProf& out, int kind) {
    const auto t0 = std::chrono::steady_clock::now();
    const size_t m = all.k.size();
    ++gs_prof[kind];
    gs_before[kind] += static_cast<long long>(m);
    out.k.clear();
    out.w.clear();
    const int who = gs_who;
    {
      const auto [dp, dq] = gz_distinct(all.k);
      gz_before[who].push_back(static_cast<int>(m));
      gz_dpb[who].push_back(dp);
      gz_dqb[who].push_back(dq);
    }
    auto rec_after = [&]() {
      const auto [dp, dq] = gz_distinct(out.k);
      gz_after[who].push_back(static_cast<int>(out.k.size()));
      gz_dp[who].push_back(dp);
      gz_dq[who].push_back(dq);
    };
    if (pr_geps <= 0.0 && pr_gB <= 0) {
      out = std::move(all);
      gs_after[kind] += static_cast<long long>(m);
      gs_share[kind] += 1.0;
      rec_after();
      pr_sel_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      return;
    }
    std::vector<double> io(m);
    double tot = NEG_INF;
    for (size_t k = 0; k < m; ++k) {
      const auto [ul, ur] = pr_U(k_hi(all.k[k]), k_lo(all.k[k]));
      const double u = log_add(ul, ur);
      io[k] = u > HALF_NEG ? all.w[k] + u : NEG_INF;
      tot = log_add(tot, io[k]);
    }
    std::vector<size_t> ord(m);
    for (size_t k = 0; k < m; ++k) ord[k] = k;
    std::vector<char> keep(m, 0);
    const size_t cap = pr_gB > 0 ? std::min(m, static_cast<size_t>(pr_gB)) : m;
    if (tot < HALF_NEG) {
      ++gs_zero[kind];
      std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return all.w[a] > all.w[b]; });
      for (size_t c = 0; c < cap; ++c) keep[ord[c]] = 1;
      gs_share[kind] += 1.0;
    } else {
      std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return io[a] > io[b]; });
      const double target = pr_geps > 0.0 ? tot + std::log1p(-pr_geps) : tot;
      double kept = NEG_INF;
      size_t c = 0;
      for (size_t k : ord) {
        if (c >= cap || io[k] < HALF_NEG) break;
        keep[k] = 1;
        ++c;
        kept = log_add(kept, io[k]);
        if (pr_geps > 0.0 && kept >= target) break;
      }
      gs_share[kind] += std::exp(kept - tot);
    }
    for (size_t k = 0; k < m; ++k)
      if (keep[k]) {
        out.k.push_back(all.k[k]);
        out.w.push_back(all.w[k]);
      }
    gs_after[kind] += static_cast<long long>(out.k.size());
    rec_after();
    pr_sel_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }

  // Wall time of the outermost growth / root / key build, less the selection inside it.
  int pr_gdepth = 0;
  struct PrGrowTimer {
    DP03& d;
    std::chrono::steady_clock::time_point t0;
    double s0;
    explicit PrGrowTimer(DP03& x) : d(x), t0(std::chrono::steady_clock::now()), s0(x.pr_sel_sec) {
      ++d.pr_gdepth;
    }
    ~PrGrowTimer() {
      if (--d.pr_gdepth == 0) {
        const double dt =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() -
            (d.pr_sel_sec - s0);
        d.pr_grow_sec += dt;
        if (d.prof_on) {
          const int r = d.prof_region;
          if (r >= 0 && r < 4) d.prof_grow_reg[r] += dt;
        }
      }
    }
  };

  // Growth profiles (scan left, final once built): G(x, y) = [end_ok] at(x, y) at key (x, y)
  // plus the steps' w G(x', y'), cut (bound: selected); roots (x, y, inner) the same,
  // untruncated (bound: selected per step when a Qg state evaluates them, pr_root_eval).
  std::unordered_map<uint64_t, PProf> prG, prGR;
  // --pr-fast 2: growth helpers HF / HM over the G profiles (final, so the helpers are too);
  // dropped whenever they hold more than PR_HG_MAX entries (only between node builds).
  std::unordered_map<uint64_t, PProf> prHg;
  static constexpr long long PR_HG_MAX = 8000000;
  void pr_hg_trim() {
    if (pr_hentries[0] - pr_hg_ent <= PR_HG_MAX) return;
    prHg.clear();
    pr_hg_ent = pr_hentries[0];
    ++pr_hg_clears;
    gm_help = 0;
  }
  auto pr_gchild() {
    return [this](int a, int b) -> const PProf* { return &prG.find(key2(a, b))->second; };
  }

  void pr_gbuild(int x, int y) {
    if (prG.count(key2(x, y))) return;
    std::vector<std::pair<int, int>> stk{{x, y}};
    std::vector<ChainNbr> nb;
    while (!stk.empty()) {
      const auto [a, b] = stk.back();
      if (prG.count(key2(a, b))) {
        stk.pop_back();
        continue;
      }
      pr_steps(a, b, false, true, nb);
      if (nb.empty() && p2_g_fb && p2_g_mode > 0 && !g_in_set(a, b)) g_one_step(a, b, nb);
      bool ready = true;
      for (const ChainNbr& t : nb)
        if (!prG.count(key2(t.x, t.y))) {
          stk.emplace_back(t.x, t.y);
          ready = false;
        }
      if (!ready) continue;
      count_steps(0, nb);
      PProf p;
      if (pr_fast > 0) {
        pr_hg_trim();
        pr_node_acc(a, b, false, true, nb, E.end_ok(a, b) ? E.at(a, b) : NEG_INF, pr_gchild(), prHg);
        if (pr_gsel) {
          PProf all;
          dacc_finish(all, 0.0);
          gs_who = 0;
          pr_gselect(all, p, 0);
        } else {
          pr_cut_drop += pr_finish_out_d(p, pr_cut);
        }
        pr_linify(p);
        ++pr_G_nodes;
        pr_G_entries += static_cast<long long>(p.k.size());
        gm_nodes += pbytes(p);
        gm_update();
        prG.emplace(key2(a, b), std::move(p));
        stk.pop_back();
        continue;
      }
      PAcc acc;
      if (E.end_ok(a, b)) acc.add(key2(a, b), E.at(a, b));
      for (const ChainNbr& t : nb) {
        const PProf& q = prG.find(key2(t.x, t.y))->second;
        for (size_t k = 0; k < q.k.size(); ++k) acc.add(q.k[k], t.w + q.w[k]);
      }
      if (pr_gsel) {
        PProf all;
        acc.finish(all, 0.0);
        gs_who = 0;
        pr_gselect(all, p, 0);
      } else {
        pr_cut_drop += pr_finish_out(acc, p, pr_cut);
      }
      ++pr_G_nodes;
      pr_G_entries += static_cast<long long>(p.k.size());
      gm_nodes += pbytes(p);
      gm_update();
      prG.emplace(key2(a, b), std::move(p));
      stk.pop_back();
    }
  }

  const PProf& pr_root(int x, int y, bool inner) {
    const uint64_t key = key3(x, y, inner ? 1 : 0);
    auto it = prGR.find(key);
    if (it != prGR.end()) return it->second;
    PrGrowTimer gt(*this);
    std::vector<ChainNbr> nb;
    pr_steps(x, y, inner, true, nb);
    if (nb.empty() && p2_g_fb && p2_g_mode > 0 && !g_in_set(x, y)) g_one_step(x, y, nb);
    count_steps(0, nb);
    for (const ChainNbr& t : nb) pr_gbuild(t.x, t.y);
    PProf p;
    if (pr_fast > 0) {
      pr_hg_trim();
      pr_node_acc(x, y, inner, true, nb, E.end_ok(x, y) ? E.at(x, y) : NEG_INF, pr_gchild(), prHg);
      if (pr_gapfree && pr_cut_roots) pr_cut_drop += pr_finish_out_d(p, pr_cut);
      else dacc_finish(p, 0.0);
      pr_linify(p);
      ++pr_roots;
      gm_roots += pbytes(p);
      gm_update();
      return prGR.emplace(key, std::move(p)).first->second;
    }
    PAcc acc;
    if (E.end_ok(x, y)) acc.add(key2(x, y), E.at(x, y));
    for (const ChainNbr& t : nb) {
      const PProf& q = prG.find(key2(t.x, t.y))->second;
      for (size_t k = 0; k < q.k.size(); ++k) acc.add(q.k[k], t.w + q.w[k]);
    }
    if (pr_gapfree && pr_cut_roots) pr_cut_drop += pr_finish_out(acc, p, pr_cut);
    else acc.finish(p, 0.0);
    ++pr_roots;
    gm_roots += pbytes(p);
    gm_update();
    return prGR.emplace(key, std::move(p)).first->second;
  }

  // The root profile a Qg state at this step evaluates (bound: selected, cached per step).
  std::unordered_map<uint64_t, PProf> prGRs;
  const PProf& pr_root_eval(int x, int y, bool inner) {
    if (!pr_gsel) return pr_root(x, y, inner);
    const uint64_t key = key3(x, y, inner ? 1 : 0);
    auto it = prGRs.find(key);
    if (it != prGRs.end()) return it->second;
    PrGrowTimer gt(*this);
    PProf all = pr_root(x, y, inner), out;
    gs_who = 1;
    pr_gselect(all, out, 1);
    gm_step += pbytes(out);
    gm_update();
    return prGRs.emplace(key, std::move(out)).first->second;
  }

  // Key profiles of the helper states at step q: Qgls (i, seed) closes by (f) into
  // (i-1, j), j - q - 1 <= cap (seed) or pr_ml_close; Qgrs c by closings_rs.  Per step.
  std::unordered_map<int, PProf> prKls, prKrs;

  void pr_ls_closings(int i, bool seed, int q, std::vector<Closing>& out) const {
    out.clear();
    const int x = i - 1;
    if (x < 0) return;
    const int j_hi = std::min(n - 1, q + 1 + pr_ml_close);
    for (int j = q + 1; j <= j_hi; ++j) {
      const int unpaired = j - q - 1;
      if (seed && unpaired > cap) continue;
      if (!E.pairs(x, j) || !E.end_ok(x, j)) continue;
      double w_ext = E.w_ap + E.w_bp + E.w_cp * unpaired + pair_bonus;
      if (dz) w_ext += E.run(q + 1, j - 1);
      out.push_back({x, j, w_ext + E.at(x, j)});
    }
  }

  const PProf& pr_key(bool ls, int p, bool seed, int q) {
    std::unordered_map<int, PProf>& m = ls ? prKls : prKrs;
    const int key = ls ? 2 * p + (seed ? 1 : 0) : p;
    auto it = m.find(key);
    if (it != m.end()) return it->second;
    PrGrowTimer gt(*this);
    std::vector<Closing> cl;
    if (ls) pr_ls_closings(p, seed, q, cl);
    else closings_rs(p, q, cl);
    PProf out;
    if (pr_fast > 0) {
      std::vector<std::pair<const PProf*, double>> rs;
      rs.reserve(cl.size());
      for (const Closing& c : cl) rs.emplace_back(&pr_root(c.x, c.y, false), c.w);
      double ref = NEG_INF;
      for (const auto& [r, w] : rs)
        if (r->mx > HALF_NEG) ref = std::max(ref, w + r->mx);
      dacc_begin(ref > HALF_NEG ? ref : 0.0);
      for (const auto& [r, w] : rs) dacc_add(*r, w);
      if (pr_gsel) {
        PProf all;
        dacc_finish(all, 0.0);
        gs_who = 2;
        pr_gselect(all, out, 1);
      } else if (pr_gapfree && pr_cut_roots) {
        pr_cut_drop += pr_finish_out_d(out, pr_cut);
      } else {
        dacc_finish(out, 0.0);
      }
      ++pr_keys;
      pr_key_entries += static_cast<long long>(out.k.size());
      gm_step += pbytes(out);
      gm_update();
      return m.emplace(key, std::move(out)).first->second;
    }
    PAcc acc;
    for (const Closing& c : cl) {
      const PProf& r = pr_root(c.x, c.y, false);
      for (size_t k = 0; k < r.k.size(); ++k) acc.add(r.k[k], c.w + r.w[k]);
    }
    if (pr_gsel) {
      PProf all;
      acc.finish(all, 0.0);
      gs_who = 2;
      pr_gselect(all, out, 1);
    } else if (pr_gapfree && pr_cut_roots) {
      pr_cut_drop += pr_finish_out(acc, out, pr_cut);
    } else {
      acc.finish(out, 0.0);
    }
    ++pr_keys;
    pr_key_entries += static_cast<long long>(out.k.size());
    gm_step += pbytes(out);
    gm_update();
    return m.emplace(key, std::move(out)).first->second;
  }

 public:
  // ------------------------------------------------ --grow-rec: layered marginal vectors
  // RN of a growth node / root / key / helper: L[m][q'] = sum_p' G_m[p', q'] C[p'-1] and
  // R[m][p'] = sum_q' G_m[p', q'] e^beta1 Sbar[q'+1], G_m the growth weight with m multiloop
  // steps (m = 2: two or more).  Same steps, tables and caps as the profiles (pr_steps,
  // pr_split, pr_helper's ranges); interior steps keep m, multiloop steps add one.  Layers
  // m >= rnL (L) or >= rnR (R) are not built.  Values: e^mx v[i - lo].
  struct RV {
    int lo = 0;
    double mx = NEG_INF;
    std::vector<double> v;
  };
  struct RN {
    RV L[3], R[3];
    double tot = NEG_INF;  // log of the summed built layers (L and R), the --grow-ml-c rank
    // --grow-node-trunc: kept (index, log value) of L (0) and R (1), layers summed, index order;
    // rest: log of the R mass dropped in this node's subtree (bounds any missing R entry)
    std::vector<std::pair<int, double>> T[2];
    double rest = NEG_INF;
    double rsum = NEG_INF;  // --grow-r-from-v: log of the R side's sum (R itself not built)
  };
  static double rbytes(const RN& r) {
    size_t c = 0;
    for (int m = 0; m < 3; ++m) c += r.L[m].v.capacity() + r.R[m].v.capacity();
    c += 2 * (r.T[0].capacity() + r.T[1].capacity());
    return 8.0 * static_cast<double>(c) + static_cast<double>(sizeof(RN) + 24);
  }
  int rnL = 3, rnR = 3;
  std::unordered_map<uint64_t, RN> prRG, prRR, prRH;
  std::unordered_map<int, RN> prRKls, prRKrs;
  // R side (layer vectors and kept R lists) of the stored growth nodes, roots, helpers, keys
  double rside_bytes(long long& e) const {
    double b = 0.0;
    e = 0;
    auto add = [&](const RN& r) {
      for (int m = 0; m < 3; ++m) {
        e += static_cast<long long>(r.R[m].v.size());
        b += 8.0 * static_cast<double>(r.R[m].v.capacity());
      }
      e += static_cast<long long>(r.T[1].size());
      b += 16.0 * static_cast<double>(r.T[1].capacity());
    };
    for (const auto* M : {&prRG, &prRR, &prRH})
      for (const auto& kv : *M) add(kv.second);
    for (const auto* M : {&prRKls, &prRKrs})
      for (const auto& kv : *M) add(kv.second);
    return b;
  }
  long long rh_elems = 0, rh_elems0 = 0, rh_built = 0, rh_terms = 0, rh_clears = 0, rn_nodes = 0,
            rn_roots = 0, rn_keys = 0, rn_terms = 0, rn_mlc_drop_l = 0, rn_mlc_drop_r = 0;
  static long long relems(const RN& r) {
    long long c = 0;
    for (int m = 0; m < 3; ++m) c += static_cast<long long>(r.L[m].v.size() + r.R[m].v.size());
    return c + static_cast<long long>(r.T[0].size() + r.T[1].size());
  }
  // dense accumulators of the six output vectors (rel. to rRef), touched ranges
  std::vector<double> rA[6];
  double rRef[6];
  int rLo[6], rHi[6];
  std::vector<int> rTouch[6];
  bool rSparse = false;
  struct RTerm {
    const RN* c;
    double w;
    int sh;
  };
  std::vector<RTerm> rTerms;
  double rSelfR = NEG_INF;  // --grow-r-from-v: the node's own R stop (log)
  void racc_begin() {
    rSparse = g_ntr_K > 0 && (p2_g_cap > 0 || p2_g_mode > 0 || p2_g_eager);
    for (int s = 0; s < 6; ++s) {
      if (rA[s].size() < static_cast<size_t>(n)) rA[s].assign(n, 0.0);
      rRef[s] = NEG_INF;
      rLo[s] = n;
      rHi[s] = -1;
      if (rSparse) rTouch[s].clear();
    }
    rSelfR = NEG_INF;
  }
  // slot of layer m on side (0 L, 1 R), -1 if not built
  int rslot(int side, int m) const { return m < (side ? rnR : rnL) ? side * 3 + m : -1; }
  static int rdst(int m, int sh) { return sh ? std::min(m + 1, 2) : m; }
  // Eager growth freezes the prefix at C_pre. Lazy (the default) uses the scan's C.
  double g_pref(int x) const {
    if (p2_g_eager && stem_ && !stem_->prC.empty())
      return x > 0 ? stem_->prC[static_cast<size_t>(x - 1)] : 0.0;
    return prefix_log(x);
  }
  void racc_self(int x, int y, double at) {
    const StemEst& st = *stem_;
    const double px = g_pref(x);
    rRef[0] = std::max(rRef[0], at + (lLam.empty() ? px : px + lLam[x]));
    if (g_rfv) return;
    rRef[3] = std::max(rRef[3], at + E.w_beta1 + st.prS[y + 1]);
  }
  void racc_refs() {
    if (g_ntr_K > 0) {
      for (const RTerm& t : rTerms) {
        for (const auto& [i, g] : t.c->T[0]) rRef[0] = std::max(rRef[0], t.w + g);
        if (g_rfv) continue;
        if (g_ntr_KR >= 0) {
          for (const auto& [i, g] : t.c->T[1]) rRef[3] = std::max(rRef[3], t.w + g);
        } else if (!t.c->R[0].v.empty()) {
          rRef[3] = std::max(rRef[3], t.w + t.c->R[0].mx);
        }
      }
      return;
    }
    for (const RTerm& t : rTerms)
      for (int side = 0; side < 2; ++side)
        for (int m = 0; m < 3; ++m) {
          const RV& a = side ? t.c->R[m] : t.c->L[m];
          if (a.v.empty()) continue;
          const int s = rslot(side, rdst(m, t.sh));
          if (s >= 0) rRef[s] = std::max(rRef[s], t.w + a.mx);
        }
  }
  void racc_add_self(int x, int y, double at) {
    const StemEst& st = *stem_;
    auto put = [&](int s, int i, double w) {
      if (w < HALF_NEG || rRef[s] < HALF_NEG) return;
      if (rSparse && rA[s][i] == 0.0) rTouch[s].push_back(i);
      rA[s][i] += std::exp(w - rRef[s]);
      rLo[s] = std::min(rLo[s], i);
      rHi[s] = std::max(rHi[s], i);
    };
    const double px = g_pref(x);
    put(0, y, at + (lLam.empty() ? px : px + lLam[x]));
    if (g_rfv) {
      rSelfR = at + E.w_beta1 + st.prS[y + 1];
      return;
    }
    put(3, x, at + E.w_beta1 + st.prS[y + 1]);
  }
  void racc_add_terms() {
    if (g_ntr_K > 0) {
      for (const RTerm& t : rTerms)
        for (int side = 0; side < 2; ++side) {
          const int s = side * 3;
          if (side == 1 && g_rfv) continue;
          if (side == 1 && g_ntr_KR < 0) {
            const RV& a = t.c->R[0];
            if (a.v.empty()) continue;
            const double sc = std::exp((t.w + a.mx) - rRef[s]);
            if (sc == 0.0) continue;
            double* d = rA[s].data() + a.lo;
            for (size_t i = 0; i < a.v.size(); ++i) d[i] += sc * a.v[i];
            if (rSparse)
              for (size_t i = 0; i < a.v.size(); ++i)
                if (a.v[i] != 0.0) rTouch[s].push_back(a.lo + static_cast<int>(i));
            if (lpf) lpf_racc_dense += static_cast<long long>(a.v.size());
            rLo[s] = std::min(rLo[s], a.lo);
            rHi[s] = std::max(rHi[s], a.lo + static_cast<int>(a.v.size()) - 1);
            continue;
          }
          for (const auto& [i, g] : t.c->T[side]) {
            const double v = std::exp((t.w + g) - rRef[s]);
            if (v == 0.0) continue;
            if (rSparse && rA[s][i] == 0.0) rTouch[s].push_back(i);
            rA[s][i] += v;
            rLo[s] = std::min(rLo[s], i);
            rHi[s] = std::max(rHi[s], i);
          }
        }
      return;
    }
    for (const RTerm& t : rTerms)
      for (int side = 0; side < 2; ++side)
        for (int m = 0; m < 3; ++m) {
          const RV& a = side ? t.c->R[m] : t.c->L[m];
          if (a.v.empty()) continue;
          const int s = rslot(side, rdst(m, t.sh));
          if (s < 0) continue;
          const double sc = std::exp((t.w + a.mx) - rRef[s]);
          if (sc == 0.0) continue;
          double* d = rA[s].data() + a.lo;
          const double* src = a.v.data();
          const size_t len = a.v.size();
          for (size_t i = 0; i < len; ++i) d[i] += sc * src[i];
          rLo[s] = std::min(rLo[s], a.lo);
          rHi[s] = std::max(rHi[s], a.lo + static_cast<int>(len) - 1);
        }
  }
  // --grow-node-trunc ranking (static): L entry q' by g uL(q'), uL(q') = sum over the V groups
  // (arm r, y' > q') of Qzbar(0, y'-1) prYw / Qz[xs][x'] / Qzbar(0, q') (tr_uL with the
  // separable fill); entries with uL = 0 dropped (no arm's lK is nonzero there).  R entry p'
  // by g, the p' that are some arm's xs first (lsc_R reads R at an arm's xs).
  std::vector<double> ntrUL;
  std::vector<char> ntrXs;
  long long ntr_in[2] = {}, ntr_kept[2] = {}, ntr_nodes = 0;
  mutable long long ntr_rmiss = 0;  // lsc_R lookups answered by the rest bound
  struct NE {
    int i;
    double g, s;
  };
  std::vector<NE> ntrE;
  int ntr_first_j = -2;
  long long ntr_arm_skip = 0;
  // --l2r-qz-audit on: reads of the scan's Qz past the filled column, total and via pr_qz
  bool qz_audit = false, qz_audit_req = false;
  mutable long long qz_unf = 0, qz_unf_pr = 0, qz_reads = 0;
  void ntr_init() {
    if (!ntrUL.empty()) return;
    const StemEst& st = *stem_;
    pr_ub_init();
    ntr_first_j = cur_j;
    std::vector<double> acc(n + 1, NEG_INF);
    for (int r = 0; r < prA; ++r) {
      const int xs = prArmXs[r], xp = prArmXp[r];
      double iz = NEG_INF;
      if (xs >= 1) {
        if (g_ntr_zbar && !pr_qz_unp) {
          const double z = st.Qzbar(xs, xp);
          if (z > HALF_NEG) iz = (prKap[xs - 1] + prKap[xp]) - z;
        } else {
          iz = pr_iqz(xs, xp);
        }
      }
      if (iz < HALF_NEG) {
        ++ntr_arm_skip;
        continue;
      }
      for (int k = prArmY0[r]; k < prArmY0[r + 1]; ++k) {
        const double z = st.Qzbar(0, prYy[k] - 1);
        if (z > HALF_NEG) acc[prYy[k]] = log_add(acc[prYy[k]], (z + prYw[k]) + iz);
      }
    }
    ntrUL.assign(n, NEG_INF);
    double suf = NEG_INF;
    for (int q = n - 1; q >= 0; --q) {
      suf = log_add(suf, acc[q + 1]);
      const double z = st.Qzbar(0, q);
      ntrUL[q] = suf > HALF_NEG ? suf - (z > HALF_NEG ? z : 0.0) : NEG_INF;
    }
    ntrXs.assign(n, 0);
    for (int r = 0; r < prA; ++r) ntrXs[prArmXs[r]] = 1;
  }
  void racc_finish_trunc(RN& r) {
    ntr_init();
    double tot = NEG_INF;
    r.rest = NEG_INF;
    for (const RTerm& t : rTerms)
      if (t.c->rest > HALF_NEG) r.rest = log_add(r.rest, t.w + t.c->rest);
    for (int side = 0; side < 2; ++side) {
      const size_t K = static_cast<size_t>(side && g_ntr_KR > 0 ? g_ntr_KR : g_ntr_K);
      const int s = side * 3;
      std::vector<std::pair<int, double>>& o = r.T[side];
      o.clear();
      if (side == 1 && g_rfv) {
        double rs = rSelfR;
        for (const RTerm& t : rTerms)
          if (t.c->rsum > HALF_NEG) rs = log_add(rs, t.w + t.c->rsum);
        r.R[0] = RV();
        r.rsum = rs;
        tot = log_add(tot, rs);
        continue;
      }
      if (side == 1 && g_ntr_KR < 0) {
        tot = log_add(tot, rv_finish(s, r.R[0]));
        continue;
      }
      if (rHi[s] < rLo[s]) continue;
      double* a = rA[s].data();
      ntrE.clear();
      auto one = [&](int i) {
        if (a[i] > 0.0) {
          const double g = rRef[s] + std::log(a[i]);
          tot = log_add(tot, g);
          if (side == 0) {
            if (p2_g_suffix) {
              const double sb = stem_->prS[static_cast<size_t>(i + 1)];
              if (sb > HALF_NEG) ntrE.push_back({i, g, g + sb});
            } else if (ntrUL[i] > HALF_NEG) {
              ntrE.push_back({i, g, g + ntrUL[i]});
            }
          } else {
            ntrE.push_back({i, g, ntrXs[i] ? 0.0 : NEG_INF});
          }
        }
        a[i] = 0.0;
      };
      if (rSparse) {
        std::vector<int>& tt = rTouch[s];
        std::sort(tt.begin(), tt.end());
        tt.erase(std::unique(tt.begin(), tt.end()), tt.end());
        for (int i : tt) one(i);
      } else {
        for (int i = rLo[s]; i <= rHi[s]; ++i) one(i);
      }
      ntr_in[side] += static_cast<long long>(ntrE.size());
      auto better = [side](const NE& x, const NE& y) {
        if (side == 1 && x.s != y.s) return x.s > y.s;
        const double vx = side ? x.g : x.s, vy = side ? y.g : y.s;
        if (vx != vy) return vx > vy;
        return x.i < y.i;
      };
      if (ntrE.size() > K) {
        std::nth_element(ntrE.begin(), ntrE.begin() + static_cast<long>(K), ntrE.end(), better);
        if (side == 1)
          for (size_t k = K; k < ntrE.size(); ++k) r.rest = log_add(r.rest, ntrE[k].g);
        ntrE.resize(K);
      }
      std::sort(ntrE.begin(), ntrE.end(), [](const NE& x, const NE& y) { return x.i < y.i; });
      o.reserve(ntrE.size());
      for (const NE& e : ntrE) o.emplace_back(e.i, e.g);
      ntr_kept[side] += static_cast<long long>(o.size());
    }
    ++ntr_nodes;
    r.tot = tot;
    rTerms.clear();
  }
  void racc_finish(RN& r) {
    if (g_ntr_K > 0) {
      racc_finish_trunc(r);
      return;
    }
    double tot = NEG_INF;
    for (int s = 0; s < 6; ++s) {
      RV& o = s < 3 ? r.L[s] : r.R[s - 3];
      const double t = rv_finish(s, o);
      if (t > HALF_NEG) tot = log_add(tot, t);
    }
    r.tot = tot;
    rTerms.clear();
  }
  // slot s of the accumulator into o (dense, relative to its max), reset; log of its sum
  double rv_finish(int s, RV& o) {
    o = RV();
    if (rHi[s] < rLo[s]) return NEG_INF;
    if (lpf) {
      ++lpf_racc_n;
      lpf_racc_fin += rHi[s] - rLo[s] + 1;
    }
    int lo = rLo[s], hi = rHi[s];
    double* a = rA[s].data();
    while (lo <= hi && a[lo] == 0.0) ++lo;
    while (hi >= lo && a[hi] == 0.0) --hi;
    double t = NEG_INF;
    if (lo <= hi) {
      double mxv = 0.0, sum = 0.0;
      for (int i = lo; i <= hi; ++i) {
        mxv = std::max(mxv, a[i]);
        sum += a[i];
      }
      o.lo = lo;
      o.mx = rRef[s] + std::log(mxv);
      o.v.resize(static_cast<size_t>(hi - lo + 1));
      for (int i = lo; i <= hi; ++i) o.v[i - lo] = a[i] / mxv;
      t = rRef[s] + std::log(sum);
    }
    for (int i = rLo[s]; i <= rHi[s]; ++i) a[i] = 0.0;
    return t;
  }
  void rec_hg_trim() {
    if (rh_elems - rh_elems0 <= PR_HG_MAX) return;
    std::unique_ptr<SecTimer> ptm(lpf ? new SecTimer(lpf_hg_sec) : nullptr);
    if (lpf) ++lpf_hg_clears;
    prRH.clear();
    {
      std::unique_ptr<SecTimer> mtm(lpf ? new SecTimer(lpf_mtrim_sec) : nullptr);
      mtrim();
    }
    rh_elems0 = rh_elems;
    ++rh_clears;
    gm_help = 0;
  }
  // --grow-ml-c: keep the C terms with the largest w + tot (ties: first listed)
  void rec_topc(std::vector<RTerm>& t, size_t from, long long& drop) {
    const size_t C = static_cast<size_t>(g_ml_c);
    if (C == 0 || t.size() - from <= C) return;
    std::stable_sort(t.begin() + static_cast<long>(from), t.end(), [](const RTerm& a, const RTerm& b) {
      return a.w + a.c->tot > b.w + b.c->tot;
    });
    drop += static_cast<long long>(t.size() - from - C);
    t.resize(from + C);
  }
  // terms of the helper (xo, y, kind) after its --grow-ml-c cut
  void rec_helper_terms(int xo, int y, int kind, std::vector<RTerm>& t, long long& drop) {
    const StemEst& st = *stem_;
    // Pair set P when --g-pairs is on. Flags off keeps prKeptRight, the current list.
    const std::vector<int>& kr =
        (p2_g_mode > 0 ? st.prKeptRightG[xo] : st.prKeptRight[xo]);
    const int S = pr_ml_span;
    auto close = [&](int a, int b) { return (E.w_ap + E.w_bp + E.at(a, b)) + pair_bonus; };
    t.clear();
    if (kind != 0) {
      const int rf = kind == 1 ? S : cap;
      for (auto jt = std::upper_bound(kr.begin(), kr.end(), y); jt != kr.end() && *jt - y - 1 <= rf; ++jt) {
        const int yo = *jt;
        if (!E.end_ok(xo, yo)) continue;
        double w = E.w_cp * (yo - y - 1);
        if (dz) w += E.run(y + 1, yo - 1);
        t.push_back({&prRG.find(key2(xo, yo))->second, w + close(xo, yo), 0});
      }
    }
    for (auto jt = std::upper_bound(kr.begin(), kr.end(), y + 1); jt != kr.end() && *jt - y - 1 <= S; ++jt) {
      const int yo = *jt;
      const double mb = st.PrMbar(y + 1, yo - 1);
      if (mb < HALF_NEG || !E.end_ok(xo, yo)) continue;
      t.push_back({&prRG.find(key2(xo, yo))->second, mb + close(xo, yo), 0});
    }
    if (prof_on) prof_help_steps += static_cast<long long>(t.size());
    rec_topc(t, 0, drop);
  }
  bool g_in_set(int x, int y) const {
    if (p2_g_mode <= 0 || !stem_ || x < 0 || x >= n) return true;
    const std::vector<int>& kr = stem_->prKeptRightG[static_cast<size_t>(x)];
    return std::binary_search(kr.begin(), kr.end(), y);
  }
  // One enclosing pair in P, the nearest left end that has a right end past y.
  void g_one_step(int x, int y, std::vector<ChainNbr>& nb) const {
    if (!stem_ || p2_g_mode <= 0) return;
    const StemEst& st = *stem_;
    const std::vector<std::vector<int>>& KR = st.prKeptRightG;
    for (int xo = x - 1; xo >= 0; --xo) {
      const std::vector<int>& kr = KR[static_cast<size_t>(xo)];
      auto it = std::upper_bound(kr.begin(), kr.end(), y);
      if (it == kr.end()) continue;
      const int yo = *it;
      const int gapL = x - xo - 1, gapR = yo - y - 1;
      if (gapL + gapR <= cap) {
        const double s = E.stem_interior(xo, x, y, yo);
        if (s > HALF_NEG) {
          nb.push_back({xo, yo, 0, s + pair_bonus});
          return;
        }
      }
      if (E.end_ok(x, y) && E.end_ok(xo, yo) && gapR <= pr_ml_span && y + 1 <= yo - 1) {
        const double mb = st.PrMbar(y + 1, yo - 1);
        if (mb > HALF_NEG) {
          double rl = E.w_cp * gapL;
          if (dz) rl += E.run(xo + 1, x - 1);
          const double br = E.w_bp + E.at(x, y);
          const double close = (E.w_ap + E.w_bp + E.at(xo, yo)) + pair_bonus;
          nb.push_back({xo, yo, 2, ((mb + br) + rl) + close});
          return;
        }
      }
    }
  }
  // --grow-ml-ub: log sum over the helper's terms of e^(w + child tot), cached for the run.
  // Bounds the helper's tot from above: the helper adds, per term, the child's T[0] (a subset of
  // its L entries) and its R side, while the child's tot sums all its L entries and all of R.
  static constexpr double kMlUbMargin = 1e-7;
  std::unordered_map<uint64_t, double> prRHub;
  long long ml_ub_cand = 0, ml_ub_built = 0, ml_ub_nodes = 0;
  double rec_helper_ub(int xo, int y, int kind) {
    const uint64_t key = key3(xo, y, kind);
    auto it = prRHub.find(key);
    if (it != prRHub.end()) return it->second;
    std::vector<RTerm> t;
    long long drop = 0;
    rec_helper_terms(xo, y, kind, t, drop);
    double u = NEG_INF;
    for (const RTerm& e : t)
      if (e.c->tot > HALF_NEG) u = log_add(u, e.w + e.c->tot);
    if (u > HALF_NEG) u += kMlUbMargin * std::max(1.0, std::fabs(u));
    return prRHub.emplace(key, u).first->second;
  }
  // helper of pr_helper's (xo, y, kind) over the RN children (layers unshifted)
  const RN& rec_helper(int xo, int y, int kind) {
    const uint64_t key = key3(xo, y, kind);
    auto it = prRH.find(key);
    if (it != prRH.end()) return it->second;
    std::unique_ptr<SecTimer> ptm(lpf ? new SecTimer(lpf_help_sec) : nullptr);
    if (lpf) {
      ++lpf_help_built;
      if (!lpf_help_keys.insert(key).second) ++lpf_help_rebuilt;
    }
    std::vector<RTerm> t;
    rec_helper_terms(xo, y, kind, t, rn_mlc_drop_r);
    std::vector<RTerm> keep;
    keep.swap(rTerms);
    rTerms = std::move(t);
    racc_begin();
    racc_refs();
    racc_add_terms();
    ++rh_built;
    rh_terms += static_cast<long long>(rTerms.size());
    RN r;
    racc_finish(r);
    rTerms.swap(keep);
    // --grow-r-from-v: helpers hold no R side, so their fixed size counts toward the clear limit
    rh_elems += relems(r) + (g_rfv ? static_cast<long long>(sizeof(RN) + 64) / 8 : 0);
    gm_help += rbytes(r);
    gm_help_entries_peak = std::max(gm_help_entries_peak, rh_elems - rh_elems0);
    gm_update();
    return prRH.emplace(key, std::move(r)).first->second;
  }
  // The multiloop terms rec_topc keeps from s (same terms, same order), building a helper only
  // while its bound (cached tot, else rec_helper_ub) is not below the C-th best tot so far.
  void rec_ml_ub(int y, bool inner, const MlSplit& s, std::vector<RTerm>& ml) {
    struct Cand {
      int xo, kind, idx;
      double w, ub;
    };
    std::vector<Cand> cd;
    cd.reserve(s.qm.size() + s.lr.size());
    for (const auto& [xo, qv] : s.qm) cd.push_back({xo, inner ? 2 : 1, static_cast<int>(cd.size()), qv + s.br, 0.0});
    for (const auto& [xo, rl] : s.lr) cd.push_back({xo, 0, static_cast<int>(cd.size()), rl + s.br, 0.0});
    for (Cand& c : cd) {
      auto it = prRH.find(key3(c.xo, y, c.kind));
      c.ub = c.w + (it != prRH.end() ? it->second.tot : rec_helper_ub(c.xo, y, c.kind));
    }
    std::vector<int> ord(cd.size());
    for (size_t k = 0; k < cd.size(); ++k) ord[k] = static_cast<int>(k);
    std::sort(ord.begin(), ord.end(), [&](int a, int b) { return cd[a].ub > cd[b].ub || (cd[a].ub == cd[b].ub && a < b); });
    const size_t C = static_cast<size_t>(g_ml_c);
    struct Got {
      RTerm t;
      int idx;
      double v;
    };
    std::vector<Got> got;
    std::vector<double> top;  // the C largest v so far, ascending
    for (int k : ord) {
      const Cand& c = cd[k];
      if (top.size() == C && c.ub < top.front()) break;
      const RN& h = rec_helper(c.xo, y, c.kind);
      const RTerm t{&h, c.w, 1};
      const double v = t.w + t.c->tot;
      got.push_back({t, c.idx, v});
      if (top.size() < C) top.insert(std::upper_bound(top.begin(), top.end(), v), v);
      else if (v > top.front()) {
        top.erase(top.begin());
        top.insert(std::upper_bound(top.begin(), top.end(), v), v);
      }
    }
    std::sort(got.begin(), got.end(), [](const Got& a, const Got& b) { return a.v > b.v || (a.v == b.v && a.idx < b.idx); });
    ++ml_ub_nodes;
    ml_ub_cand += static_cast<long long>(cd.size());
    ml_ub_built += static_cast<long long>(got.size());
    rn_mlc_drop_l += static_cast<long long>(cd.size() - C);
    for (size_t k = 0; k < C && k < got.size(); ++k) ml.push_back(got[k].t);
  }
  // RN of x.y with steps nb (pr_steps(x, y, inner, true)): the stop at x.y, interior steps
  // into the children, multiloop steps through the helpers as pr_node_acc with --pr-fast 2.
  void rec_node(int x, int y, bool inner, const std::vector<ChainNbr>& nb, RN& out) {
    rec_hg_trim();
    rTerms.clear();
    for (const ChainNbr& t : nb)
      if (t.t == 0) rTerms.push_back({&prRG.find(key2(t.x, t.y))->second, t.w, 0});
    MlSplit s;
    pr_split(x, y, !geq_pre, s);
    if (s.on) {
      if (lpf) {
        lpf_ml_cand += static_cast<long long>(s.qm.size() + s.lr.size());
        lpf_ml_max = std::max(lpf_ml_max, static_cast<long long>(s.qm.size() + s.lr.size()));
        ++lpf_ml_nodes;
      }
      std::vector<RTerm> ml;
      if (g_ml_ub && g_ml_c > 0 && s.qm.size() + s.lr.size() > static_cast<size_t>(g_ml_c)) rec_ml_ub(y, inner, s, ml);
      else {
        for (const auto& [xo, qv] : s.qm) ml.push_back({&rec_helper(xo, y, inner ? 2 : 1), qv + s.br, 1});
        for (const auto& [xo, rl] : s.lr) ml.push_back({&rec_helper(xo, y, 0), rl + s.br, 1});
        rec_topc(ml, 0, rn_mlc_drop_l);
      }
      rTerms.insert(rTerms.end(), ml.begin(), ml.end());
    }
    rn_terms += static_cast<long long>(rTerms.size());
    const double at = E.end_ok(x, y) ? E.at(x, y) : NEG_INF;
    racc_begin();
    if (at > HALF_NEG) racc_self(x, y, at);
    racc_refs();
    if (at > HALF_NEG) racc_add_self(x, y, at);
    racc_add_terms();
    racc_finish(out);
  }
  // R side entries per built node / root / key: count, sum, max
  long long rz_n[3] = {}, rz_sum[3] = {}, rz_max[3] = {};
  void rec_stat(int who, const RN& r) {
    gz_after[who].push_back(static_cast<int>(relems(r)));
    {
      long long c = static_cast<long long>(r.T[1].size());
      for (int m = 0; m < 3; ++m) c += static_cast<long long>(r.R[m].v.size());
      ++rz_n[who];
      rz_sum[who] += c;
      rz_max[who] = std::max(rz_max[who], c);
    }
    if (g_ntr_K > 0) {
      gz_dq[who].push_back(static_cast<int>(r.T[0].size()));
      gz_dp[who].push_back(static_cast<int>(r.T[1].size()));
      return;
    }
    auto nz = [&](const RV* a, int nl) {
      int lo = n, hi = -1;
      for (int m = 0; m < nl; ++m)
        if (!a[m].v.empty()) {
          lo = std::min(lo, a[m].lo);
          hi = std::max(hi, a[m].lo + static_cast<int>(a[m].v.size()) - 1);
        }
      int c = 0;
      for (int i = lo; i <= hi; ++i) {
        bool z = true;
        for (int m = 0; m < nl && z; ++m)
          if (!a[m].v.empty() && i >= a[m].lo && i < a[m].lo + static_cast<int>(a[m].v.size()) &&
              a[m].v[i - a[m].lo] > 0.0)
            z = false;
        if (!z) ++c;
      }
      return c;
    };
    gz_dq[who].push_back(nz(r.L, rnL));
    gz_dp[who].push_back(nz(r.R, rnR));
  }
  void rec_gbuild(int x, int y) {
    if (prRG.count(key2(x, y))) return;
    std::vector<std::pair<int, int>> stk{{x, y}};
    std::vector<ChainNbr> nb;
    while (!stk.empty()) {
      const auto [a, b] = stk.back();
      if (prRG.count(key2(a, b))) {
        stk.pop_back();
        continue;
      }
      pr_steps(a, b, false, !geq_pre, nb);
      if (nb.empty() && p2_g_fb && p2_g_mode > 0 && !g_in_set(a, b)) g_one_step(a, b, nb);
      bool ready = true;
      for (const ChainNbr& t : nb)
        if (!prRG.count(key2(t.x, t.y))) {
          stk.emplace_back(t.x, t.y);
          ready = false;
        }
      if (!ready) continue;
      count_steps(0, nb);
      RN r;
      rec_node(a, b, false, nb, r);
      ++rn_nodes;
      rec_stat(0, r);
      gm_nodes += rbytes(r);
      gm_update();
      prRG.emplace(key2(a, b), std::move(r));
      stk.pop_back();
    }
  }
  // One root from the swept gamma. g at q' already includes the scan prefix, as rec_root stores it.
  bool sw_g_ready = false;  // one-shot gate for G availability (direct fill or states prGamma)
  RN sw_as_rn(int x, int y, bool inner) {
    if (p2_sweep && l2r_mode != 0 && stem_ && !sw_g_ready) {
      sw_g_ready = true;
      if (g_sweep == G_SWEEP_DIRECT) {
        sw_filled = true;
        // ntrUL must exist before sw_fill's cap_ga (same L trunc key as racc_finish_trunc).
        if (g_ntr_K > 0 && !p2_g_suffix) {
          const auto t0 = std::chrono::steady_clock::now();
          ntr_init();
          sw_ntr_init_sec +=
              std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        }
        {
          const auto t0 = std::chrono::steady_clock::now();
          sw_fill(0, n);
          sw_fill(1, std::min(pr_loop, cap));
          sw_fill_sec +=
              std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        }
        // Entries per G node after fill (slot 0 = full right gap).
        for (const auto& kv : sw_ga[0]) {
          const long long e = static_cast<long long>(kv.second.size());
          ++sw_ga_nodes;
          sw_ga_ent_sum += e;
          if (e > sw_ga_ent_max) sw_ga_ent_max = e;
        }
      } else {
        // --g-sweep states: G already in stem_->prGamma from sw_hyper_push.
        // Fill time is the G merge+cut copied onto this DP from the nested pass.
        sw_filled = false;
        sw_fill_sec = g_merge_sec + g_cut_sec;
        if (stem_) {
          for (const auto& kv : stem_->prGamma) {
            const long long e = static_cast<long long>(kv.second.size());
            ++sw_ga_nodes;
            sw_ga_ent_sum += e;
            if (e > sw_ga_ent_max) sw_ga_ent_max = e;
          }
        }
      }
    }
    const auto t_read0 = std::chrono::steady_clock::now();
    RN r;
    std::unordered_map<int, double> byq;
    const std::unordered_map<uint64_t, double>* src = nullptr;
    std::unordered_map<uint64_t, std::unordered_map<uint64_t, double>>::const_iterator git;
    if (sw_filled) {
      auto it = sw_ga[inner ? 1 : 0].find(key2(x, y));
      if (it != sw_ga[inner ? 1 : 0].end()) src = &it->second;
    } else if (stem_ && !stem_->prGamma.empty()) {
      git = stem_->prGamma.find(key2(x, y));
      if (git != stem_->prGamma.end()) src = &git->second;
    } else {
      auto it = sw_ga[inner ? 1 : 0].find(key2(x, y));
      if (it != sw_ga[inner ? 1 : 0].end()) src = &it->second;
    }
    if (src) {
      sw_as_rn_ent_in += static_cast<long long>(src->size());
      for (const auto& [k, w] : *src) {
        const int q = k_lo(k);
        const double g = w + prefix_log(k_hi(k));
        auto jt = byq.find(q);
        if (jt == byq.end()) byq.emplace(q, g);
        else jt->second = log_add(jt->second, g);
      }
    }
    double tot = NEG_INF;
    for (const auto& [q, g] : byq) tot = log_add(tot, g);
    r.tot = tot;
    if (g_ntr_K > 0) {
      ntr_init();
      struct E {
        int i;
        double g, s;
      };
      std::vector<E> v;
      v.reserve(byq.size());
      for (const auto& [q, g] : byq) {
        double s = NEG_INF;
        if (p2_g_suffix) {
          const double sb = stem_->prS[static_cast<size_t>(q + 1)];
          if (sb > HALF_NEG) s = g + sb;
        } else if (q < static_cast<int>(ntrUL.size()) && ntrUL[q] > HALF_NEG) {
          s = g + ntrUL[q];
        }
        if (s > HALF_NEG) v.push_back({q, g, s});
      }
      const size_t K = static_cast<size_t>(g_ntr_K);
      if (v.size() > K) {
        std::nth_element(v.begin(), v.begin() + static_cast<long>(K), v.end(),
                         [](const E& a, const E& b) {
                           if (a.s != b.s) return a.s > b.s;
                           return a.i < b.i;
                         });
        v.resize(K);
      }
      std::sort(v.begin(), v.end(), [](const E& a, const E& b) { return a.i < b.i; });
      r.T[0].reserve(v.size());
      for (const E& e : v) r.T[0].emplace_back(e.i, e.g);
    } else if (!byq.empty()) {
      int lo = n, hi = -1;
      double mx = NEG_INF;
      for (const auto& [q, g] : byq) {
        lo = std::min(lo, q);
        hi = std::max(hi, q);
        mx = std::max(mx, g);
      }
      RV& a = r.L[0];
      a.lo = lo;
      a.mx = mx;
      a.v.assign(static_cast<size_t>(hi - lo + 1), 0.0);
      for (const auto& [q, g] : byq)
        if (g > HALF_NEG) a.v[static_cast<size_t>(q - lo)] = std::exp(g - mx);
    }
    ++sw_as_rn_calls;
    sw_as_rn_ent_out += static_cast<long long>(r.T[0].size());
    sw_as_rn_sec +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t_read0).count();
    return r;
  }
  const RN& rec_root(int x, int y, bool inner) {
    const uint64_t key = key3(x, y, inner ? 1 : 0);
    auto it = prRR.find(key);
    const bool dump_focus = std::getenv("DUMP_ROOT") &&
                            ((x == 4 && y == 35) || (x == 9 && y == 37) || (x == 9 && y == 41));
    auto dump_t0 = [&](const RN& placed) {
      if (!dump_focus) return;
      std::fprintf(stderr,
                   "DUMP_ROOT_CREATE\tx\t%d\ty\t%d\tinner\t%d\tp2_sweep\t%d\tK\t%d\tntr\t%d\n", x, y,
                   inner ? 1 : 0, p2_sweep ? 1 : 0, g_ntr_K, ntrUL.empty() ? 0 : 1);
      if (!placed.T[0].empty()) {
        for (const auto& e : placed.T[0])
          std::fprintf(stderr, "DUMP_ROOT_T0\tx\t%d\ty\t%d\tq\t%d\tg\t%.17g\n", x, y, e.first,
                       e.second);
        std::fprintf(stderr, "DUMP_ROOT_T0META\tx\t%d\ty\t%d\tn\t%zu\tsrc\tT0\n", x, y,
                     placed.T[0].size());
      } else {
        auto lam = lamL(placed);
        for (const auto& e : lam)
          std::fprintf(stderr, "DUMP_ROOT_T0\tx\t%d\ty\t%d\tq\t%d\tg\t%.17g\n", x, y, e.first,
                       e.second);
        std::fprintf(stderr, "DUMP_ROOT_T0META\tx\t%d\ty\t%d\tn\t%zu\tsrc\tlamL\n", x, y, lam.size());
      }
    };
    if (it != prRR.end()) {
      if (dump_focus)
        std::fprintf(stderr, "DUMP_ROOT_HIT\tx\t%d\ty\t%d\tinner\t%d\tp2_sweep\t%d\n", x, y,
                     inner ? 1 : 0, p2_sweep ? 1 : 0);
      return it->second;
    }
    if (p2_sweep) {
      PrGrowTimer gt(*this);
      RN r = sw_as_rn(x, y, inner);
      ++rn_roots;
      auto& placed = prRR.emplace(key, std::move(r)).first->second;
      dump_t0(placed);
      return placed;
    }
    PrGrowTimer gt(*this);
    std::vector<ChainNbr> nb;
    pr_steps(x, y, inner, !geq_pre, nb);
    if (nb.empty() && p2_g_fb && p2_g_mode > 0 && !g_in_set(x, y)) g_one_step(x, y, nb);
    count_steps(0, nb);
    for (const ChainNbr& t : nb) rec_gbuild(t.x, t.y);
    RN r;
    rec_node(x, y, inner, nb, r);
    ++rn_roots;
    rec_stat(1, r);
    gm_roots += rbytes(r);
    gm_update();
    auto& placed = prRR.emplace(key, std::move(r)).first->second;
    dump_t0(placed);
    return placed;
  }
  const RN& rec_key(bool ls, int p, bool seed, int q) {
    const auto tick = []() { return std::chrono::steady_clock::now(); };
    auto t0 = tick();
    auto charge = [&](std::chrono::steady_clock::time_point t1) {
      sx_rkey_sec += std::chrono::duration<double>(t1 - t0).count();
    };
    std::unordered_map<int, RN>& m = ls ? prRKls : prRKrs;
    const int key = ls ? 2 * p + (seed ? 1 : 0) : p;
    auto it = m.find(key);
    if (it != m.end()) {
      charge(tick());
      return it->second;
    }
    PrGrowTimer gt(*this);
    std::vector<Closing> cl;
    if (ls) pr_ls_closings(p, seed, q, cl);
    else closings_rs(p, q, cl);
    charge(tick());
    std::vector<RTerm> t;
    for (const Closing& c : cl) t.push_back({&rec_root(c.x, c.y, false), c.w, 0});
    t0 = tick();
    rTerms = std::move(t);
    racc_begin();
    racc_refs();
    racc_add_terms();
    RN r;
    racc_finish(r);
    ++rn_keys;
    rec_stat(2, r);
    gm_step += rbytes(r);
    gm_update();
    auto placed = m.emplace(key, std::move(r)).first;
    charge(tick());
    return placed->second;
  }

  // DUMP_P2_ALL=prefix: write Lambda_L/R (roots+nodes+keys), V-key lists, Qz_pre for offline
  // sweep-vs-lazy comparison (Ther |d| gap diagnosis).
  void dump_p2_all(const char* prefix) const {
    if (!prefix || !*prefix) return;
    auto dec3 = [](uint64_t k, int& a, int& b, int& c) {
      a = static_cast<int>(k >> 42);
      b = static_cast<int>((k >> 21) & 0x1fffff);
      c = static_cast<int>(k & 0x1fffff);
    };
    auto wmap = [&](const char* tag, const std::unordered_map<uint64_t, RN>& M, int side) {
      std::string path = std::string(prefix) + "." + tag + ".tsv";
      std::ofstream out(path);
      out << "key\tx\ty\tinner\ti\tg\n";
      for (const auto& [k, r] : M) {
        int x, y, inn;
        dec3(k, x, y, inn);
        if (g_ntr_K > 0) {
          for (const auto& e : r.T[side])
            out << k << '\t' << x << '\t' << y << '\t' << inn << '\t' << e.first << '\t'
                << std::setprecision(17) << e.second << '\n';
        } else {
          const RV& a = side ? r.R[0] : r.L[0];
          for (size_t i = 0; i < a.v.size(); ++i) {
            if (!(a.v[i] > 0.0)) continue;
            out << k << '\t' << x << '\t' << y << '\t' << inn << '\t' << (a.lo + (int)i) << '\t'
                << std::setprecision(17) << (a.mx + std::log(a.v[i])) << '\n';
          }
        }
      }
      std::fprintf(stderr, "[dump-p2] %s n=%zu\n", path.c_str(), M.size());
    };
    auto wkey = [&](const char* tag, const std::unordered_map<int, RN>& M, int side) {
      std::string path = std::string(prefix) + "." + tag + ".tsv";
      std::ofstream out(path);
      out << "key\ti\tg\n";
      for (const auto& [k, r] : M) {
        if (g_ntr_K > 0) {
          for (const auto& e : r.T[side])
            out << k << '\t' << e.first << '\t' << std::setprecision(17) << e.second << '\n';
        }
      }
      std::fprintf(stderr, "[dump-p2] %s n=%zu\n", path.c_str(), M.size());
    };
    wmap("Lambda_L_root", prRR, 0);
    wmap("Lambda_R_root", prRR, 1);
    wmap("Lambda_L_node", prRG, 0);
    wmap("Lambda_R_node", prRG, 1);
    wkey("K_L_ls", prRKls, 0);
    wkey("K_R_ls", prRKls, 1);
    wkey("K_L_rs", prRKrs, 0);
    wkey("K_R_rs", prRKrs, 1);
    for (int s = 0; s < 2; ++s) {
      std::string path = std::string(prefix) + ".V_s" + std::to_string(s) + ".tsv";
      std::ofstream out(path);
      out << "xs\ty\tw\n";
      for (const auto& [xs, vec] : prVkey[s])
        for (const auto& e : vec)
          out << xs << '\t' << e.first << '\t' << std::setprecision(17) << e.second << '\n';
      std::fprintf(stderr, "[dump-p2] %s keys=%zu\n", path.c_str(), prVkey[s].size());
    }
    if (stem_) {
      std::string path = std::string(prefix) + ".Qz_pre.tsv";
      std::ofstream out(path);
      out << "i\tj\tz\n";
      long long nq = 0;
      if (!stem_->qz_col.empty()) {
        for (int j = 0; j < n; ++j)
          for (const auto& e : stem_->qz_col[static_cast<size_t>(j)]) {
            out << e.k << '\t' << j << '\t' << std::setprecision(17) << e.q << '\n';
            ++nq;
          }
      } else if (!stem_->prQz.empty()) {
        for (int i = 0; i < n; ++i)
          for (int j = i; j < n; ++j) {
            const double z = stem_->prQz[static_cast<size_t>(i) * n + j];
            if (z < HALF_NEG) continue;
            out << i << '\t' << j << '\t' << std::setprecision(17) << z << '\n';
            ++nq;
          }
      }
      std::fprintf(stderr, "[dump-p2] %s n=%lld\n", path.c_str(), nq);
    }
  }

  // --grow-trunc: the kept (index, log g) of g_L (q') and g_R (p') of an evaluated root / key,
  // summed over the used layers; per step (u_L, u_R move with q).
  struct RT {
    std::vector<std::pair<int, double>> L, R;
  };
  std::unordered_map<const RN*, RT> prRT;
  std::vector<int> tr_before[2], tr_after[2], tr_kL[2], tr_kR[2];
  double tr_share[2] = {}, tr_zero[2] = {};
  double tr_uL(int qp) {
    if (qp + 1 > n - 1) return NEG_INF;
    const PrUL& u = pr_ul(qp);
    return u.lmax > HALF_NEG ? u.phi : NEG_INF;
  }
  double tr_uR(int pp) {
    if (pp < 1) return NEG_INF;
    const PrUR& r = pr_ur(pp);
    return r.lmax > HALF_NEG ? r.r1 : NEG_INF;
  }
  const RT& rec_trunc(const RN& r, int who) {
    auto it = prRT.find(&r);
    if (it != prRT.end()) return it->second;
    if (g_ntr_K > 0) throw std::runtime_error("--grow-node-trunc: no layer vectors for rec_trunc");
    const auto t0 = std::chrono::steady_clock::now();
    pr_ub_init();
    pr_ub_sync();
    // entries: side, index, log g, log score
    struct E4 {
      int side, i;
      double g, s;
    };
    std::vector<E4> e[2];
    for (int side = 0; side < 2; ++side) {
      const RV* a = side ? r.R : r.L;
      const int nl = std::min(side ? g_ml_R + 1 : g_ml_L + 1, side ? rnR : rnL);
      int lo = n, hi = -1;
      for (int m = 0; m < nl; ++m)
        if (!a[m].v.empty()) {
          lo = std::min(lo, a[m].lo);
          hi = std::max(hi, a[m].lo + static_cast<int>(a[m].v.size()) - 1);
        }
      for (int i = lo; i <= hi; ++i) {
        double g = NEG_INF;
        for (int m = 0; m < nl; ++m) {
          const int k = i - a[m].lo;
          if (k >= 0 && k < static_cast<int>(a[m].v.size()) && a[m].v[k] > 0.0)
            g = log_add(g, a[m].mx + std::log(a[m].v[k]));
        }
        if (g < HALF_NEG) continue;
        // --g-key suffix: rank g_L(q') by G_L * S̄_pre(q'+1). Beam-DP03 uses this
        // path (--grow-trunc); LR-DP03 ranks in racc_finish_trunc.
        double u;
        if (p2_g_suffix && side == 0 && stem_)
          u = stem_->prS[static_cast<size_t>(i + 1)];
        else
          u = side ? tr_uR(i) : tr_uL(i);
        e[side].push_back({side, i, g, u > HALF_NEG ? g + u : NEG_INF});
      }
    }
    RT t;
    const size_t before = e[0].size() + e[1].size();
    double share_num = NEG_INF, share_den = NEG_INF;
    auto pick = [&](std::vector<E4>& v) {
      double tot = NEG_INF;
      for (const E4& x : v) tot = log_add(tot, x.s);
      const size_t cap = std::min(v.size(), static_cast<size_t>(g_tr_K));
      std::vector<size_t> ord(v.size());
      for (size_t k = 0; k < v.size(); ++k) ord[k] = k;
      double kept = NEG_INF;
      if (tot < HALF_NEG) {
        tr_zero[who] += 1.0;
        std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return v[a].g > v[b].g; });
        for (size_t c = 0; c < cap; ++c) (v[ord[c]].side ? t.R : t.L).emplace_back(v[ord[c]].i, v[ord[c]].g);
      } else {
        std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return v[a].s > v[b].s; });
        const double target = g_tr_eps > 0.0 ? tot + std::log1p(-g_tr_eps) : tot;
        size_t c = 0;
        for (size_t k : ord) {
          if (c >= cap || v[k].s < HALF_NEG) break;
          (v[k].side ? t.R : t.L).emplace_back(v[k].i, v[k].g);
          ++c;
          kept = log_add(kept, v[k].s);
          if (g_tr_eps > 0.0 && kept >= target) break;
        }
        share_num = log_add(share_num, kept);
        share_den = log_add(share_den, tot);
      }
    };
    if (g_tr_pool) {
      std::vector<E4> all = e[0];
      all.insert(all.end(), e[1].begin(), e[1].end());
      pick(all);
    } else {
      pick(e[0]);
      pick(e[1]);
    }
    tr_before[who].push_back(static_cast<int>(before));
    tr_after[who].push_back(static_cast<int>(t.L.size() + t.R.size()));
    tr_kL[who].push_back(static_cast<int>(t.L.size()));
    tr_kR[who].push_back(static_cast<int>(t.R.size()));
    tr_share[who] += share_den > HALF_NEG ? std::exp(share_num - share_den) : 1.0;
    if (p2_g_suffix || p2_g_mode > 0) ntr_kept[0] += static_cast<long long>(t.L.size());
    gm_step += 16.0 * static_cast<double>(t.L.capacity() + t.R.capacity()) + sizeof(RT) + 32;
    gm_update();
    pr_sel_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return prRT.emplace(&r, std::move(t)).first->second;
  }

  // Score of a stem with inner pair u.v from the RN r of its outer pairs: per layer m,
  //   SL[m] = sum_q' L[m][q'] tau_q'(u, v)     (sep: sum_q' L[m][q'] / Qz[q+1][q'] times tau_q)
  //   SR[m] = sum_p' R[m][p'] kappa_p'(u, v)   (--pr-abl-pk K: the top K p' by the used g_R)
  // and sigma = sum_{m <= ml_L} SL[m] + sum_{m <= ml_R} SR[m] (pr_eval_sep with the layer sums
  // for g_L, g_R).  --grow-ml-diag: every layer evaluated, mass kept of each configuration.
  std::vector<double> rTau;
  std::vector<int> rTauStamp;
  int rTauSt = 0;
  static constexpr int RD_CFG = 5, RD_BINS = 2000;
  static constexpr int RD_L[RD_CFG] = {2, 2, 1, 1, 0}, RD_R[RD_CFG] = {2, 0, 0, 1, 0};
  long long rd_n = 0, rd_hist[RD_CFG][RD_BINS + 1] = {};
  double rd_sum[RD_CFG] = {}, rd_layer_sum[2][3] = {};
  double rd_min[RD_CFG] = {1, 1, 1, 1, 1};
  double rec_tau(int qp, int u, int v) {
    if (rTauStamp.size() < static_cast<size_t>(n)) {
      rTauStamp.assign(n, 0);
      rTau.assign(n, NEG_INF);
    }
    if (rTauStamp[qp] == rTauSt) return rTau[qp];
    rTauStamp[qp] = rTauSt;
    PrPartnerTimer pt(*this);
    return rTau[qp] = pr_T(qp, u, v, false).ext;
  }
  double pr_eval_rec(const RN& r, int u, int v, int who) {
    if (g_ntr_K > 0) throw std::runtime_error("--grow-node-trunc: pr_eval_rec needs the layer vectors");
    if (!E.end_ok(u, v)) return NEG_INF;
    ++rTauSt;
    ++pr_sep_evals;
    const int q = pr_abl_q;
    const bool vecL = !(pr_abl & PA_Q);
    if (g_tr_K > 0) {
      const RT& t = rec_trunc(r, who);
      double sl = NEG_INF, sr = NEG_INF;
      for (const auto& [qp, g] : t.L) {
        if (vecL) {
          const double ta = rec_tau(qp, u, v);
          if (ta > HALF_NEG) sl = log_add(sl, g + ta);
        } else {
          const double z = pr_abl_qz(q, qp);
          if (z > HALF_NEG) sl = log_add(sl, g - z);
        }
      }
      if (!vecL && sl > HALF_NEG) {
        const double ta = rec_tau(q, u, v);
        sl = ta > HALF_NEG ? sl + ta : NEG_INF;
      }
      for (const auto& [pp, g] : t.R) {
        PrPartnerTimer pt(*this);
        const double kx = pr_K(pp, u, v).ext;
        if (kx > HALF_NEG) sr = log_add(sr, g + kx);
      }
      const double s = log_add(sl, sr);
      if (s < HALF_NEG) return NEG_INF;
      return s + (E.w_stem + E.at(u, v));
    }
    const int nL = g_ml_diag ? rnL : g_ml_L + 1, nR = g_ml_diag ? rnR : g_ml_R + 1;
    double SL[3] = {NEG_INF, NEG_INF, NEG_INF}, SR[3] = {NEG_INF, NEG_INF, NEG_INF};
    for (int m = 0; m < nL; ++m) {
      const RV& a = r.L[m];
      double acc = NEG_INF;
      for (size_t i = 0; i < a.v.size(); ++i) {
        if (a.v[i] <= 0.0) continue;
        const int qp = a.lo + static_cast<int>(i);
        const double g = a.mx + std::log(a.v[i]);
        if (vecL) {
          const double t = rec_tau(qp, u, v);
          if (t > HALF_NEG) acc = log_add(acc, g + t);
        } else {
          const double z = pr_abl_qz(q, qp);
          if (z > HALF_NEG) acc = log_add(acc, g - z);
        }
      }
      if (!vecL && acc > HALF_NEG) {
        const double t = rec_tau(q, u, v);
        acc = t > HALF_NEG ? acc + t : NEG_INF;
      }
      SL[m] = acc;
    }
    if (pr_abl_pk > 0) {
      // top K p' by the configured g_R, then every layer over those p'
      std::vector<std::pair<double, int>> gr;
      std::vector<double> g(n, NEG_INF);
      for (int m = 0; m <= g_ml_R && m < rnR; ++m) {
        const RV& a = r.R[m];
        for (size_t i = 0; i < a.v.size(); ++i)
          if (a.v[i] > 0.0) g[a.lo + i] = log_add(g[a.lo + i], a.mx + std::log(a.v[i]));
      }
      for (int p = 0; p < n; ++p)
        if (g[p] > HALF_NEG) gr.emplace_back(g[p], p);
      std::stable_sort(gr.begin(), gr.end(), [](const auto& x, const auto& y) { return x.first > y.first; });
      if (gr.size() > static_cast<size_t>(pr_abl_pk)) gr.resize(static_cast<size_t>(pr_abl_pk));
      for (const auto& [gw, pp] : gr) {
        PrPartnerTimer pt(*this);
        const double kx = pr_K(pp, u, v).ext;
        if (kx < HALF_NEG) continue;
        for (int m = 0; m < nR; ++m) {
          const RV& a = r.R[m];
          const int i = pp - a.lo;
          if (i < 0 || i >= static_cast<int>(a.v.size()) || a.v[i] <= 0.0) continue;
          SR[m] = log_add(SR[m], a.mx + std::log(a.v[i]) + kx);
        }
      }
    } else {
      for (int m = 0; m < nR; ++m) {
        const RV& a = r.R[m];
        double acc = NEG_INF;
        for (size_t i = 0; i < a.v.size(); ++i) {
          if (a.v[i] <= 0.0) continue;
          PrPartnerTimer pt(*this);
          const double kx = pr_K(a.lo + static_cast<int>(i), u, v).ext;
          if (kx > HALF_NEG) acc = log_add(acc, a.mx + std::log(a.v[i]) + kx);
        }
        SR[m] = acc;
      }
    }
    double s = NEG_INF;
    for (int m = 0; m <= g_ml_L && m < nL; ++m) s = log_add(s, SL[m]);
    for (int m = 0; m <= g_ml_R && m < nR; ++m) s = log_add(s, SR[m]);
    if (g_ml_diag) {
      double all = NEG_INF;
      for (int m = 0; m < 3; ++m) all = log_add(all, log_add(SL[m], SR[m]));
      if (all > HALF_NEG) {
        ++rd_n;
        for (int m = 0; m < 3; ++m) {
          if (SL[m] > HALF_NEG) rd_layer_sum[0][m] += std::exp(SL[m] - all);
          if (SR[m] > HALF_NEG) rd_layer_sum[1][m] += std::exp(SR[m] - all);
        }
        for (int c = 0; c < RD_CFG; ++c) {
          double sc = NEG_INF;
          for (int m = 0; m <= RD_L[c]; ++m) sc = log_add(sc, SL[m]);
          for (int m = 0; m <= RD_R[c]; ++m) sc = log_add(sc, SR[m]);
          const double f = sc > HALF_NEG ? std::min(1.0, std::exp(sc - all)) : 0.0;
          rd_sum[c] += f;
          rd_min[c] = std::min(rd_min[c], f);
          // bins of -log10 f, width 0.01; the last bin holds f < 1e-20 and f = 0
          const double d = f > 0.0 ? -std::log10(f) : 1e9;
          const int b = std::min(RD_BINS, static_cast<int>(d * 100.0));
          ++rd_hist[c][b];
        }
      }
    }
    if (s < HALF_NEG) return NEG_INF;
    return s + (E.w_stem + E.at(u, v));
  }
  // --grow-rec-check: recursion score and summed vectors against the profile of the state.
  long long rc_n = 0, rc_vec_n = 0;
  double rc_max_score = 0.0, rc_max_vec = 0.0;
  void rec_check(const RN& r, const PProf& kp, double s_rec, double s_prof) {
    ++rc_n;
    if (s_rec > HALF_NEG || s_prof > HALF_NEG) {
      const double d = (s_rec > HALF_NEG && s_prof > HALF_NEG) ? std::fabs(s_rec - s_prof) : INFINITY;
      rc_max_score = std::max(rc_max_score, d);
    }
    const StemEst& st = *stem_;
    std::vector<double> L(n, NEG_INF), R(n, NEG_INF);
    for (size_t k = 0; k < kp.k.size(); ++k) {
      const int pp = k_hi(kp.k[k]), qp = k_lo(kp.k[k]);
      L[qp] = log_add(L[qp], kp.w[k] + prefix_log(pp));
      R[pp] = log_add(R[pp], kp.w[k] + E.w_beta1 + st.prS[qp + 1]);
    }
    auto cmp = [&](const RV* a, const std::vector<double>& ref) {
      std::vector<double> x(n, NEG_INF);
      double mx = NEG_INF;
      for (int m = 0; m < 3; ++m)
        for (size_t i = 0; i < a[m].v.size(); ++i)
          if (a[m].v[i] > 0.0) x[a[m].lo + i] = log_add(x[a[m].lo + i], a[m].mx + std::log(a[m].v[i]));
      for (int i = 0; i < n; ++i) mx = std::max(mx, std::max(x[i], ref[i]));
      if (mx < HALF_NEG) return;
      for (int i = 0; i < n; ++i) {
        if (std::max(x[i], ref[i]) < mx - 27.6) continue;  // entries below 1e-12 of the max
        ++rc_vec_n;
        const double d = (x[i] > HALF_NEG && ref[i] > HALF_NEG) ? std::fabs(x[i] - ref[i]) : INFINITY;
        rc_max_vec = std::max(rc_max_vec, d);
      }
    };
    cmp(r.L, L);
    cmp(r.R, R);
  }

  // log Out of a stem with inner pair u.v whose outer pairs p'.q' carry the profile kp:
  // [end_ok(u, v)] e^beta2 at(u, v) sum kp(p', q') (Stop_L + Stop_R)(p', q').
  // --pr-abl sep (with bml, q): the separable form of pr_eval,
  //   Gamma_L tau_q(u, v) + sum over the top --pr-abl-pk p' (by g_R; all if 0) g_R(p') kappa_p'(u, v),
  //   Gamma_L = sum kp(p', q') C[p'-1] / Qz[q+1][q'],  g_R(p') = sum_q' kp(p', q') e^beta1 Sbar[q'+1],
  //   tau_q = T(q; u, v).ext,  kappa_p' = K(p'; u, v).ext.
  // Without q (sepv: bml only, exact): Stop_L as sum_q' g_L(q') tau_q'(u, v) with
  // g_L(q') = sum_p' kp(p', q') C[p'-1].
  std::vector<std::pair<double, int>> prSepR, prSepL;
  double pr_eval_sep(const PProf& kp, int u, int v) {
    const StemEst& st = *stem_;
    const int q = pr_abl_q;
    const bool vecL = !(pr_abl & PA_Q);
    double gl = NEG_INF;
    prSepR.clear();
    prSepL.clear();
    for (size_t k = 0; k < kp.k.size(); ++k) {
      const int pp = k_hi(kp.k[k]), qp = k_lo(kp.k[k]);
      if (vecL) {
        prSepL.emplace_back(kp.w[k] + prefix_log(pp), qp);
      } else {
        const double z = pr_abl_qz(q, qp);
        if (z > HALF_NEG) gl = log_add(gl, kp.w[k] + prefix_log(pp) - z);
      }
      prSepR.emplace_back(kp.w[k] + (E.w_beta1 + st.prS[qp + 1]), pp);
    }
    double sv = NEG_INF;
    if (vecL) {
      std::sort(prSepL.begin(), prSepL.end(), [](const auto& a, const auto& b) { return a.second < b.second; });
      for (size_t k = 0; k < prSepL.size();) {
        double g = NEG_INF;
        const int qp = prSepL[k].second;
        for (; k < prSepL.size() && prSepL[k].second == qp; ++k) g = log_add(g, prSepL[k].first);
        if (g < HALF_NEG) continue;
        PrPartnerTimer pt(*this);
        const PrTV& t = pr_T(qp, u, v, false);
        if (t.ext > HALF_NEG) sv = log_add(sv, g + t.ext);
      }
    }
    std::sort(prSepR.begin(), prSepR.end(), [](const auto& a, const auto& b) { return a.second < b.second; });
    size_t m = 0;
    for (size_t k = 0; k < prSepR.size(); ++k) {
      if (m > 0 && prSepR[m - 1].second == prSepR[k].second)
        prSepR[m - 1].first = log_add(prSepR[m - 1].first, prSepR[k].first);
      else prSepR[m++] = prSepR[k];
    }
    prSepR.resize(m);
    std::sort(prSepR.begin(), prSepR.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
    const size_t K = pr_abl_pk > 0 ? std::min(m, static_cast<size_t>(pr_abl_pk)) : m;
    ++pr_sep_evals;
    pr_sep_pp += static_cast<long long>(m);
    pr_sep_ppk += static_cast<long long>(K);
    double s = sv;
    if (gl > HALF_NEG) {
      PrPartnerTimer pt(*this);
      const PrTV& t = pr_T(q, u, v, false);
      if (t.ext > HALF_NEG) s = gl + t.ext;
    }
    for (size_t k = 0; k < K; ++k) {
      if (prSepR[k].first < HALF_NEG) continue;
      PrPartnerTimer pt(*this);
      const PrKV& kv = pr_K(prSepR[k].second, u, v);
      if (kv.ext > HALF_NEG) s = log_add(s, prSepR[k].first + kv.ext);
    }
    if (s < HALF_NEG) return NEG_INF;
    return s + (E.w_stem + E.at(u, v));
  }

  double pr_eval(const PProf& kp, int u, int v) {
    if (!E.end_ok(u, v)) return NEG_INF;
    if (pr_abl & PA_SEP) return pr_eval_sep(kp, u, v);
    if (pr_sc_n >= PR_SC_MAX) {
      prSc.clear();
      pr_sc_n = 0;
    }
    std::unordered_map<uint64_t, double>& b = prSc[key2(u, v)];
    double acc = NEG_INF;
    for (size_t k = 0; k < kp.k.size(); ++k) {
      const double s = pr_SLR(b, k_hi(kp.k[k]), k_lo(kp.k[k]), u, v);
      if (s > HALF_NEG) acc = log_add(acc, kp.w[k] + s);
    }
    pr_evals += static_cast<long long>(kp.k.size());
    if (pr_gcheck)
      for (size_t k = 0; k < kp.k.size(); ++k) pr_gcheck_ep(k_hi(kp.k[k]), k_lo(kp.k[k]), u, v);
    if (acc < HALF_NEG) return NEG_INF;
    return acc + (E.w_stem + E.at(u, v));
  }

  // --pr-grow-check: Stop_L over every arm in the gap (exact placement) and as scored (the
  // selected arms), Stop_R as scored, against Qz[u+1][v-1] U_L and Qz[u+1][v-1] U_R.
  void pr_gcheck_ep(int pp, int qp, int u, int v) {
    const double g = pr_qz(u + 1, v - 1);
    if (g < HALF_NEG) return;
    const auto [ul, ur] = pr_U(pp, qp);
    double sl_all = NEG_INF, sl_sel = NEG_INF, sr = NEG_INF;
    if (qp + 1 <= n - 1 && u + 1 <= v - 1) {
      ArmQ& a = pr_armq(qp);
      acc_begin();
      for (int pos = prXpFirst[u + 1]; pos < prXpFirst[v]; ++pos) {
        const int r = prByXp[pos];
        if (prArmXs[r] < u + 1 || a.phi[r] < HALF_NEG) continue;
        const double gl = pr_qz(u + 1, prArmXs[r] - 1), gr = pr_qz(prArmXp[r] + 1, v - 1);
        if (gl < HALF_NEG || gr < HALF_NEG) continue;
        for (const auto& [ys, w] : pr_W(a, r, qp)) acc_add(ys, (gl + gr) + w);
      }
      for (int ys : prTouched) sl_all = log_add(sl_all, prAcc[ys] + pr_ctx(pp, ys));
      const PrTV& t = pr_T(qp, u, v, false);
      for (size_t k = 0; k < t.ys.size(); ++k) sl_sel = log_add(sl_sel, t.w[k] + pr_ctx(pp, t.ys[k]));
    }
    const PrKV& kv = pr_K(pp, u, v);
    for (size_t k = 0; k < kv.is.size(); ++k) sr = log_add(sr, kv.w[k] + pr_ctx(kv.is[k], qp));
    const double slack = std::log1p(1e-12);
    ++gc_n;
    if (sl_all > HALF_NEG) {
      const double d = ul > HALF_NEG ? sl_all - (g + ul) : INFINITY;
      gc_maxL = std::max(gc_maxL, d);
      if (d > slack) ++gc_vL;
    }
    if (sr > HALF_NEG) {
      const double d = ur > HALF_NEG ? sr - (g + ur) : INFINITY;
      gc_maxR = std::max(gc_maxR, d);
      if (d > slack) ++gc_vR;
    }
    if (sl_sel > HALF_NEG && !(sl_all > HALF_NEG && sl_sel <= sl_all + slack)) ++gc_sel_gt;
  }

  // log sum_i alpha(Qgl[i][d-1][e][f]) at(i, f) Ctx(i, q) over the surviving Qgl; per step.
  std::unordered_map<uint64_t, double> prAc;
  double pr_actx(int d, int e, int f, int q) {
    const uint64_t key = key3(d, e, f);
    auto it = prAc.find(key);
    if (it != prAc.end()) return it->second;
    const StemEst& st = *stem_;
    double s = NEG_INF;
    auto ia = stA.find(key2(d, e));
    if (ia != stA.end()) {
      const StemA& sa = ia->second;
      const auto kf = std::upper_bound(sa.fs.begin(), sa.fs.end(), f) - sa.fs.begin() - 1;
      if (kf >= 0 && sa.fs[kf] == f && sa.vals[kf] > HALF_NEG) {
        double bm = NEG_INF;
        if (!(pr_abl & PA_BML))
          for (const auto& [i, lw] : sa.il[kf]) bm = log_add(bm, lw + st.Bml(i, q));
        s = log_add((E.w_beta1 + st.prS[q + 1]) + sa.vals[kf], E.w_beta1m + bm);
      }
    }
    prAc.emplace(key, s);
    return s;
  }

  // sigma of a state of Qg / Qgl / Qgr / Qgls / Qgrs (outer p, inner u, group v) at step q.
  double pr_score4(int bt, const Ent& e, int v, int q) {
    pr_build_static();
    const int p = e.p, u = e.u;
    pr_abl_p = p;
    pr_abl_q = q;
    if (g_grow_rec && (bt == BT_QG || bt == BT_QGLS || bt == BT_QGRS)) {
      const RN& r = bt == BT_QG ? rec_root(p, q, p == u)
                                : rec_key(bt == BT_QGLS, p, bt == BT_QGLS && v == q, q);
      const double s = pr_eval_rec(r, u, v, bt == BT_QG ? 0 : 1);
      if (g_rec_check) {
        const PProf& kp = bt == BT_QG ? pr_root_eval(p, q, p == u)
                                      : pr_key(bt == BT_QGLS, p, bt == BT_QGLS && v == q, q);
        rec_check(r, kp, s, pr_eval(kp, u, v));
      }
      return e.a + s;
    }
    if (bt == BT_QG) return e.a + pr_eval(pr_root_eval(p, q, p == u), u, v);
    if (bt == BT_QGLS) return e.a + pr_eval(pr_key(true, p, v == q, q), u, v);
    if (bt == BT_QGRS) return e.a + pr_eval(pr_key(false, p, false, q), u, v);
    if (!E.end_ok(p, q)) return NEG_INF;
    const StemEst& st = *stem_;
    double s = NEG_INF;
    if (bt == BT_QGL) {
      if (pr_arms_mode == 1) {
        const double se = pr_sep_ext(q, u, v, true);
        const PrTV& t = pr_T(q, u, v, true);
        double bm = NEG_INF;
        for (size_t k = 0; k < t.ys.size(); ++k) bm = log_add(bm, t.w[k] + st.Bml(p, t.ys[k]));
        s = log_add(prefix_log(p) + se, E.w_beta1m + bm);
      } else {
        const PrTV& t = pr_T(q, u, v, true);
        if (!t.ys.empty()) {
          double bm = NEG_INF;
          if (!(pr_abl & PA_BML))
            for (size_t k = 0; k < t.ys.size(); ++k) bm = log_add(bm, t.w[k] + st.Bml(p, t.ys[k]));
          s = log_add(prefix_log(p) + t.ext, E.w_beta1m + bm);
        }
      }
      auto ia = stA.find(key2(p, u + 1));
      if (ia != stA.end())
        for (int f : ia->second.fs) {
          if (f > v - 1) break;
          const double gr = pr_qz(f + 1, v - 1);
          if (gr < HALF_NEG) continue;
          const double c = pr_actx(p, u + 1, f, q);
          if (c > HALF_NEG) s = log_add(s, gr + c);
        }
    } else {
      s = pr_actx(p, u + 1, v - 1, q);
    }
    if (s < HALF_NEG) return NEG_INF;
    return e.a + E.at(p, q) + s;
  }

  static int pr_tab(int bt) {
    return bt == BT_QG ? 0 : bt == BT_QGL ? 1 : bt == BT_QGR ? 2 : bt == BT_QGLS ? 3 : bt == BT_QGRS ? 4 : 5;
  }

  // At a new step: per-step caches cleared, q'-keyed ones below q dropped.
  void pr_step(int q) {
    if (!prSepTouched.empty()) {
      ++pr_sep_steps;
      prSepTouched.clear();
    }
    prKls.clear();
    prKrs.clear();
    prRKls.clear();
    prRKrs.clear();
    prVkey[0].clear();
    prVkey[1].clear();
    prRT.clear();
    gm_step = 0;
    prAc.clear();
    pr_q = q;
    if (pr_gsel || pr_gcheck) {
      prGRs.clear();
      prUc.clear();
      for (auto it = prUL.begin(); it != prUL.end();) it = it->first < q ? prUL.erase(it) : std::next(it);
    }
    for (auto it = prTc.begin(); it != prTc.end();) it = it->first < q ? prTc.erase(it) : std::next(it);
    for (auto it = prAQ.begin(); it != prAQ.end();) it = it->first < q ? prAQ.erase(it) : std::next(it);
  }

  // --stem-pr-check: Out_hat of a ranked state with every DP03 path (growth, future stem,
  // partner, closing, context) as its own term, long double.
  struct PrPath {
    int x, y;
    long double w;
  };
  std::unordered_map<uint64_t, std::vector<PrPath>> prnG, prnV;

  static double lse(const std::vector<long double>& t) {
    if (t.empty()) return NEG_INF;
    long double m = *std::max_element(t.begin(), t.end()), s = 0.0L;
    for (long double x : t) s += std::exp(x - m);
    return static_cast<double>(m + std::log(s));
  }

  // Weight of the step x.y -> xo.yo, every step kind tested directly.
  double prn_step(int x, int y, int xo, int yo, bool inner, bool scan) const {
    const StemEst& st = *stem_;
    if (st.skip_kept) StemEst::lean_abort("prKept");
    if (!(xo < x && yo > y) || !st.prKept[ix(xo, yo)]) return NEG_INF;
    std::vector<long double> t;
    const int lc = std::min(pr_loop, cap);
    if (x - xo - 1 <= lc && yo - y - 1 <= lc) {
      const double s = E.stem_interior(xo, x, y, yo);
      if (s > HALF_NEG) t.push_back(static_cast<long double>(s) + pair_bonus);
    }
    if (pr_ml && (scan || pr_partner_ml) && E.end_ok(x, y) && E.end_ok(xo, yo)) {
      const long double br = static_cast<long double>(E.w_bp) + E.at(x, y);
      const long double cl = static_cast<long double>(E.w_ap) + E.w_bp + E.at(xo, yo) + pair_bonus;
      double lq = NEG_INF, rq = NEG_INF;
      if (xo + 1 <= x - 1 && x - (xo + 1) <= pr_ml_span)
        lq = scan ? qmp_at(x - 1).get(xo + 1) : st.PrMbar(xo + 1, x - 1);
      if (y + 1 <= yo - 1 && yo - y - 1 <= pr_ml_span) rq = st.PrMbar(y + 1, yo - 1);
      if (lq > HALF_NEG && yo - y - 1 <= (inner ? cap : pr_ml_span)) {
        long double w = lq + br + static_cast<long double>(E.w_cp) * (yo - y - 1) + cl;
        if (dz) w += E.run(y + 1, yo - 1);
        t.push_back(w);
      }
      if (ml_run_ok(x - xo - 1) && rq > HALF_NEG) {
        long double w = rq + br + static_cast<long double>(E.w_cp) * (x - xo - 1) + cl;
        if (dz) w += E.run(xo + 1, x - 1);
        t.push_back(w);
      }
      if (lq > HALF_NEG && rq > HALF_NEG) t.push_back(lq + rq + br + cl);
    }
    return lse(t);
  }

  void prn_paths(int x, int y, bool inner, bool scan, bool start_ok, long double w,
                 std::vector<PrPath>& out) const {
    if (start_ok && E.end_ok(x, y)) out.push_back({x, y, w});
    for (int xo = 0; xo < x; ++xo)
      for (int yo : stem_->prKeptRight[xo]) {
        if (yo <= y) continue;
        const double s = prn_step(x, y, xo, yo, inner, scan);
        if (s > HALF_NEG) prn_paths(xo, yo, false, scan, true, w + s, out);
      }
  }

  const std::vector<PrPath>& prn_G(int x, int y, bool inner) {
    const uint64_t key = key3(x, y, inner ? 1 : 0);
    auto it = prnG.find(key);
    if (it != prnG.end()) return it->second;
    std::vector<PrPath> v;
    prn_paths(x, y, inner, true, true, 0.0L, v);
    return prnG.emplace(key, std::move(v)).first->second;
  }

  const std::vector<PrPath>& prn_V(int x, int y) {
    const uint64_t key = key2(x, y);
    auto it = prnV.find(key);
    if (it != prnV.end()) return it->second;
    std::vector<PrPath> v;
    prn_paths(x, y, true, false, single_stems, 0.0L, v);
    return prnV.emplace(key, std::move(v)).first->second;
  }

  long double prn_ctx(int i, int j) const {
    const StemEst& st = *stem_;
    std::vector<long double> t;
    if (st.prS[j + 1] > HALF_NEG)
      t.push_back(static_cast<long double>(E.w_beta1) + prefix_log(i) + st.prS[j + 1]);
    if (st.Bml(i, j) > HALF_NEG) t.push_back(static_cast<long double>(E.w_beta1m) + st.Bml(i, j));
    return lse(t);
  }

  // Future R-stem terms: for the L-stem's outer pair pp.qp and inner gap (u, v), every kept
  // inner pair x'.y' (x' <= v-1, y' > qp) and stem path to xs.ys (xs >= u+1; xs == u+1 only
  // if exact_xs), with base added.
  void prn_future(int pp, int qp, int u, int v, bool exact_xs, long double base,
                  std::vector<long double>& t) {
    const StemEst& st = *stem_;
    for (int xp = u + 1; xp <= v - 1; ++xp)
      for (int yp : st.prKeptRight[xp]) {
        if (yp <= qp || !E.end_ok(xp, yp)) continue;
        const double gr = qz_log(xp + 1, v - 1), zb = st.Qzbar(qp + 1, yp - 1);
        if (gr < HALF_NEG || zb < HALF_NEG) continue;
        for (const PrPath& h : prn_V(xp, yp)) {
          if (h.x < u + 1 || (exact_xs && h.x != u + 1)) continue;
          const double gl = exact_xs ? 0.0 : qz_log(u + 1, h.x - 1);
          if (gl < HALF_NEG) continue;
          const long double c = prn_ctx(pp, h.y);
          if (c < HALF_NEG) continue;
          t.push_back(base + gl + gr + zb + E.w_stem + pair_bonus + E.at(xp, yp) + h.w +
                      E.at(h.x, h.y) + c);
        }
      }
  }

  // Surviving L-stems Qgl[i][pp-1][e][f] (e in [e_lo, e_hi], f <= v-1) closing a pseudoknot
  // ending at qp, read from the Qgl tables, with base added.
  void prn_partner(int pp, int qp, int u, int v, bool exact_e, long double base,
                   std::vector<long double>& t) {
    for (int f = u + 1; f <= v - 1; ++f)
      for (const Grp& g : Qgl[f].grps) {
        const int e = g.v;
        if (e < u + 1 || e > v - 1 || (exact_e && e != u + 1)) continue;
        const double gl = exact_e ? 0.0 : qz_log(u + 1, e - 1), gr = qz_log(f + 1, v - 1);
        if (gl < HALF_NEG || gr < HALF_NEG) continue;
        for (const Ent& s : g.ents) {
          if (s.u != pp - 1 || !E.end_ok(s.p, f)) continue;
          const long double c = prn_ctx(s.p, qp);
          if (c < HALF_NEG) continue;
          t.push_back(base + gl + s.a + E.at(s.p, f) + gr + c);
        }
      }
  }

  double prn_outqg(int x, int y, bool inner, int u, int v, long double w0) {
    std::vector<long double> t;
    if (!E.end_ok(u, v)) return NEG_INF;
    const long double own = static_cast<long double>(E.w_stem) + E.at(u, v);
    for (const PrPath& g : prn_G(x, y, inner)) {
      const long double base = w0 + g.w + E.at(g.x, g.y) + own;
      prn_future(g.x, g.y, u, v, false, base, t);
      prn_partner(g.x, g.y, u, v, false, base, t);
    }
    return lse(t);
  }

  double pr_naive4(int bt, const Ent& e, int v, int q) {
    const int p = e.p, u = e.u;
    if (bt == BT_QG) return prn_outqg(p, q, p == u, u, v, e.a);
    if (bt == BT_QGLS || bt == BT_QGRS) {
      std::vector<Closing> cl;
      if (bt == BT_QGLS) pr_ls_closings(p, v == q, q, cl);
      else closings_rs(p, q, cl);
      std::vector<long double> t;
      for (const Closing& c : cl) {
        const double o = prn_outqg(c.x, c.y, false, u, v, static_cast<long double>(e.a) + c.w);
        if (o > HALF_NEG) t.push_back(o);
      }
      return lse(t);
    }
    if (!E.end_ok(p, q)) return NEG_INF;
    std::vector<long double> t;
    const long double base = static_cast<long double>(e.a) + E.at(p, q);
    if (bt == BT_QGL) {
      prn_future(p, q, u, v, true, base, t);
      prn_partner(p, q, u, v, true, base, t);
    } else {
      for (const Grp& g : Qgl[v - 1].grps) {
        if (g.v != u + 1) continue;
        for (const Ent& s : g.ents) {
          if (s.u != p - 1 || !E.end_ok(s.p, v - 1)) continue;
          const long double c = prn_ctx(s.p, q);
          if (c > HALF_NEG) t.push_back(base + s.a + E.at(s.p, v - 1) + c);
        }
      }
    }
    return lse(t);
  }

  void prc_cmp(int t, double s, double ref) {
    ++prc_states[t];
    if ((ref > HALF_NEG) != (s > HALF_NEG)) ++prc_mis[t];
    else if (ref > HALF_NEG) prc_rel[t] = std::max(prc_rel[t], std::fabs(std::expm1(s - ref)));
  }

  // --stem-est-profile-eps: keep the smallest set of entries (largest first, ties by index)
  // holding 1 - eps of the part's weight, in index order; returns the entries dropped.
  static long long prof_trunc(std::vector<int>& ix, std::vector<double>& w, double eps) {
    const size_t m = w.size();
    if (m <= 1) return 0;
    double tot = NEG_INF;
    for (double x : w) tot = log_add(tot, x);
    const double target = tot + std::log1p(-eps);
    std::vector<size_t> order(m);
    for (size_t k = 0; k < m; ++k) order[k] = k;
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return w[a] > w[b]; });
    std::vector<char> keep(m, 0);
    double acc = NEG_INF;
    for (size_t k : order) {
      keep[k] = 1;
      acc = log_add(acc, w[k]);
      if (acc >= target) break;
    }
    size_t o = 0;
    for (size_t k = 0; k < m; ++k)
      if (keep[k]) {
        ix[o] = ix[k];
        w[o] = w[k];
        ++o;
      }
    ix.resize(o);
    w.resize(o);
    return static_cast<long long>(m - o);
  }

  void prof_finish_key(Prof& p, int kind) {
    if (prof_eps > 0.0) pf_trunc += prof_trunc(p.ly, p.lw, prof_eps) + prof_trunc(p.rx, p.rw, prof_eps);
    pf_ly += static_cast<long long>(p.ly.size());
    pf_rx += static_cast<long long>(p.rx.size());
    ++pf_keys[kind];
  }

  // Key profile of the Qg states with outer pair p.stc_q: the pair profile, truncated.
  const Prof& qg_prof(int p) {
    pf_note_step();
    auto it = stPg.find(p);
    if (it != stPg.end()) return it->second;
    const auto t0 = std::chrono::steady_clock::now();
    Prof k = pair_prof(p, stc_q);
    prof_finish_key(k, 0);
    pf_build_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return stPg.emplace(p, std::move(k)).first->second;
  }

  // Sum over the closings (closing order) of w + the pair profile of the closed pair.
  void prof_closings(const std::vector<Closing>& cl, Prof& out) {
    ProfAcc acc;
    for (const Closing& c : cl) {
      const Prof& pp = pair_prof(c.x, c.y);
      for (size_t k = 0; k < pp.ly.size(); ++k) acc.addL(pp.ly[k], c.w + pp.lw[k]);
      for (size_t k = 0; k < pp.rx.size(); ++k) acc.addR(pp.rx[k], c.w + pp.rw[k]);
    }
    acc.finish(out);
    out.ncl = static_cast<int>(cl.size());
  }

  // --stem-est-ml-tail nested: log sum over the case-(f) closings j > q+1+L_ml of the Qgls key
  // (i, seed) of w_f(i-1, j) + beta_N(i-1, j) (Result 20's term), j ascending.
  double ml_tail_sum(int i, bool seed, int q) {
    const int x = i - 1;
    if (x < 0) return NEG_INF;
    double acc = NEG_INF;
    for (int j = q + 2 + stem_ml_span; j <= n - 1; ++j) {
      const int unpaired = j - q - 1;
      if (seed && unpaired > cap) continue;
      if (!E.pairs(x, j) || !E.end_ok(x, j)) continue;
      const double b = stem_->B(x, j);
      if (b < HALF_NEG) continue;
      double w_ext = E.w_ap + E.w_bp + E.w_cp * unpaired + pair_bonus;
      if (dz) w_ext += E.run(q + 1, j - 1);
      acc = log_add(acc, (w_ext + E.at(x, j)) + b);
      ++pf_tail_cl;
    }
    return acc;
  }

  // Key profile of the Qgls[q] states (i, ., v) with seed = (v == q).
  const Prof& ls_prof(int i, bool seed, int q) {
    pf_note_step();
    const int key = 2 * i + (seed ? 1 : 0);
    auto it = stPls.find(key);
    if (it != stPls.end()) return it->second;
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<Closing> cl;
    closings_ls(i, seed ? q : -1, q, cl);
    Prof k;
    prof_closings(cl, k);
    if (ml_tail) k.tail = ml_tail_sum(i, seed, q);
    prof_finish_key(k, 1);
    pf_build_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return stPls.emplace(key, std::move(k)).first->second;
  }

  // Key profile of the Qgrs[q] states (c, ., .).
  const Prof& rs_prof(int c, int q) {
    pf_note_step();
    auto it = stPrs.find(c);
    if (it != stPrs.end()) return it->second;
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<Closing> cl;
    closings_rs(c, q, cl);
    Prof k;
    prof_closings(cl, k);
    prof_finish_key(k, 2);
    pf_build_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return stPrs.emplace(c, std::move(k)).first->second;
  }

  // Gap part of a profile score, alpha excluded (DP03._prof_eval): logsumexp of
  // Qz[u+1][v-1] sum_y' H_R(u, v, y') P_L(y') (consistent: times at(u, v)),
  // sum_x' Bg(x', u, v) P_R(x') and Qz[u+1][v-1] tail.
  double prof_eval(const Prof& p, int u, int v) {
    if (u != row_u || v != row_v || stc_q != row_q) {
      if (rowH.size() != static_cast<size_t>(n + 1)) {
        rowH.assign(n + 1, NEG_INF);
        rowB.assign(n + 1, NEG_INF);
        rowHs.assign(n + 1, 0);
        rowBs.assign(n + 1, 0);
      }
      row_u = u;
      row_v = v;
      row_q = stc_q;
      ++row_stamp;
    }
    double L = NEG_INF;
    for (size_t k = 0; k < p.ly.size(); ++k) {
      const int y = p.ly[k];
      if (rowHs[y] != row_stamp) {
        rowH[y] = stem_H_cached(u, v, y);
        rowHs[y] = row_stamp;
      }
      const double h = rowH[y];
      if (h > HALF_NEG) L = log_add(L, h + p.lw[k]);
    }
    const double qz = score_qz ? qz_log(u + 1, v - 1) : 0.0;
    if (L > HALF_NEG) {
      L += qz;
      if (stem_consistent) L += E.at(u, v);
    }
    double R = NEG_INF, Ro = NEG_INF;
    if (E.end_ok(u, v) && bg_reorder) R = prof_R_reorder(p, u, v);
    if (E.end_ok(u, v) && (!bg_reorder || ro_check))
      for (size_t k = 0; k < p.rx.size(); ++k) {
        const int x = p.rx[k];
        if (rowBs[x] != row_stamp) {
          rowB[x] = stem_Bg(x, u, v);
          rowBs[x] = row_stamp;
        } else {
          ++bg_lookups;
        }
        const double bg = rowB[x];
        if (bg > HALF_NEG) Ro = log_add(Ro, bg + p.rw[k]);
      }
    if (!bg_reorder) R = Ro;
    double s = log_add(L, R);
    if (p.tail > HALF_NEG) s = log_add(s, p.tail + qz);
    if (ro_check) {
      ro_off = log_add(L, Ro);
      if (p.tail > HALF_NEG) ro_off = log_add(ro_off, p.tail + qz);
    }
    pf_evals += static_cast<long long>(p.ly.size() + p.rx.size());
    return s;
  }

  // --stem-est-bg-reorder: the M row of key profile p at v (cached per step) from e = lo up,
  // built over the R entries x' in index order and, per x', the e of stE[x'] in [lo, v - 1]
  // ascending.  States are scored by gap (v, u), so the first lo = u + 1 serves the rest; a
  // smaller u rebuilds the row (M(e, v) of each e does not depend on lo).
  const MRow& prof_mrow(const Prof& p, int v, int lo) {
    std::unordered_map<int, MRow>& byv = stM[&p];
    auto it = byv.find(v);
    if (it != byv.end() && it->second.lo <= lo) return it->second;
    ++ro_rows;
    if (mAcc.size() != static_cast<size_t>(n + 1)) {
      mAcc.assign(n + 1, NEG_INF);
      mStamp.assign(n + 1, 0);
    }
    ++m_stamp;
    mTouched.clear();
    for (size_t k = 0; k < p.rx.size(); ++k) {
      const int x = p.rx[k];
      auto ie = stE.find(x);
      if (ie == stE.end()) continue;
      const std::vector<int>& es = ie->second;
      for (auto ke = std::lower_bound(es.begin(), es.end(), lo); ke != es.end(); ++ke) {
        const int e = *ke;
        if (e > v - 1) break;
        ++ro_build;
        const double b = stem_B_cached(x, e, v);
        if (b < HALF_NEG) continue;
        if (mStamp[e] != m_stamp) {
          mStamp[e] = m_stamp;
          mAcc[e] = NEG_INF;
          mTouched.push_back(e);
        }
        mAcc[e] = log_add(mAcc[e], p.rw[k] + b);
      }
    }
    std::sort(mTouched.begin(), mTouched.end());
    MRow r;
    r.lo = lo;
    r.es = mTouched;
    r.ms.reserve(mTouched.size());
    for (int e : mTouched) r.ms.push_back(mAcc[e]);
    MRow& slot = byv[v];
    slot = std::move(r);
    return slot;
  }

  // log sum_x' P_R(x') Bg(x', u, v) = beta2 + at(u, v) + log sum_{e in [u+1, v-1]}
  // Qz[u+1][e-1] M(e, v) (end_ok(u, v) checked by the caller).
  double prof_R_reorder(const Prof& p, int u, int v) {
    if (p.rx.empty()) return NEG_INF;
    const MRow& r = prof_mrow(p, v, u + 1);
    double acc = NEG_INF;
    for (size_t k = std::lower_bound(r.es.begin(), r.es.end(), u + 1) - r.es.begin();
         k < r.es.size(); ++k) {
      ++ro_eval;
      const double qz = qz_log(u + 1, r.es[k] - 1);
      if (qz < HALF_NEG) continue;
      acc = log_add(acc, qz + r.ms[k]);
    }
    if (acc < HALF_NEG) return NEG_INF;
    return acc + (E.w_stem + E.at(u, v));
  }

  // --stem-bg-reorder-check: one scored state (reordered s, original ro_off, alpha a), stored
  // at off[m]
  void ro_cmp(double s, double a, std::vector<double>& off, size_t m) {
    const double so = a + ro_off;
    off[m] = so;
    ++roc_states;
    if ((s > HALF_NEG) != (so > HALF_NEG)) ++roc_mis;
    else if (s > HALF_NEG) roc_rel = std::max(roc_rel, std::fabs(std::expm1(s - so)));
  }

  // --stem-bg-reorder-check: keep decisions of one ranking event under the original scores
  // (order of the event's score vector) against km.
  void ro_flips(std::vector<double>& off, const std::vector<char>& km, int beam) {
    std::vector<char> ko;
    keep_mask(off, ko, beam);
    ++roc_events;
    long long f = 0;
    for (size_t r = 0; r < km.size(); ++r) f += ko[r] != km[r];
    roc_flip_states += f;
    roc_flip_events += f > 0;
  }

  // Scoring order of m profile states: by gap (v, u), stable, so prof_eval's rows are reused.
  // Scores do not depend on the order.
  template <class F>
  static std::vector<size_t> gap_order(size_t m, F key) {
    std::vector<std::pair<std::pair<int, int>, size_t>> ks(m);
    for (size_t r = 0; r < m; ++r) ks[r] = {key(r), r};
    std::sort(ks.begin(), ks.end());
    std::vector<size_t> out(m);
    for (size_t r = 0; r < m; ++r) out[r] = ks[r].second;
    return out;
  }

  // --stem-profile-check (t = 0 Qg, 1 helpers)
  void prof_cmp(int t, double s, double ref) {
    ++pc_states[t];
    if ((ref > HALF_NEG) != (s > HALF_NEG)) ++pc_mis[t];
    else if (ref > HALF_NEG) pc_rel[t] = std::max(pc_rel[t], std::fabs(std::expm1(s - ref)));
  }

  // --stem-est-helpers grow / nested: rank Qgls / Qgrs by helper_sigma / helper_sigma_nested
  // (step: whole table, group: per group v; threshold rule of keep_mask).  No state is removed.
  void prune4_helper(Tab4& t, int bt, int q, int beam) {
    if (beam <= 0) return;
    const auto t0 = std::chrono::steady_clock::now();
    stem_cache_step(q);
    std::vector<const Ent*> ent;
    std::vector<int> gv;
    std::vector<std::pair<size_t, size_t>> spans;
    for (const Grp& g : t.grps) {
      const size_t lo = ent.size();
      for (const Ent& e : g.ents) {
        ent.push_back(&e);
        gv.push_back(g.v);
      }
      if (!gap_step) spans.emplace_back(lo, ent.size());
    }
    if (gap_step) spans.emplace_back(0, ent.size());
    if (pr_on && pr_arm_gaponly) {
      std::vector<std::pair<int, int>> iv;
      for (const auto& sp : spans)
        if (sp.second - sp.first > static_cast<size_t>(beam))
          for (size_t k = sp.first; k < sp.second; ++k) iv.emplace_back(ent[k]->u + 1, gv[k] - 1);
      pr_gap_begin(iv);
    }
    const bool ls = bt == BT_QGLS;
    std::vector<char> keep(ent.size(), 1), km;
    std::vector<double> sc, so;
    std::vector<Closing> cl;
    bool any = false;
    for (const auto& sp : spans) {
      const size_t m = sp.second - sp.first;
      if (m <= static_cast<size_t>(beam)) continue;
      sc.clear();
      if (pr_on) {
        const auto tp = std::chrono::steady_clock::now();
        const int ti = pr_tab(bt);
        sc.assign(m, NEG_INF);
        for (size_t r : gap_order(m, [&](size_t r) {
               return std::make_pair(gv[sp.first + r], ent[sp.first + r]->u);
             })) {
          const size_t k = sp.first + r;
          sc[r] = pr_score4(bt, *ent[k], gv[k], q);
          if (pr_dump)
            std::fprintf(pr_dump, "%d %d %d %d %d %.17g\n", bt, q, ent[k]->p, ent[k]->u, gv[k], sc[r]);
          hlp_finite += sc[r] > HALF_NEG;
          if (pr_check) prc_cmp(ti, sc[r], pr_naive4(bt, *ent[k], gv[k], q));
        }
        pr_states[ti] += static_cast<long long>(m);
        pr_sec[ti] += std::chrono::duration<double>(std::chrono::steady_clock::now() - tp).count();
      } else if (stem_profile) {
        sc.assign(m, NEG_INF);
        if (ro_check) so.assign(m, NEG_INF);
        for (size_t r : gap_order(m, [&](size_t r) {
               return std::make_pair(gv[sp.first + r], ent[sp.first + r]->u);
             })) {
          const size_t k = sp.first + r;
          const Ent& e = *ent[k];
          const Prof& kp = ls ? ls_prof(e.p, gv[k] == q, q) : rs_prof(e.p, q);
          const double s = e.a + prof_eval(kp, e.u, gv[k]);
          sc[r] = s;
          if (ro_check) ro_cmp(s, e.a, so, r);
          hlp_closings += kp.ncl;
          hlp_finite += s > HALF_NEG;
          ++pf_states;
          if (prof_check) {
            if (ls) closings_ls(e.p, gv[k], q, cl);
            else closings_rs(e.p, q, cl);
            prof_cmp(1, s, stem_chain ? chain_naive(cl, e.u, gv[k], e.a)
                                      : helper_sigma(cl, e.u, gv[k], e.a));
          }
        }
      }
      for (size_t k = sp.first; k < sp.second && !stem_profile && !pr_on; ++k) {
        const Ent& e = *ent[k];
        if (stem_helpers == 2) {
          int ncl = 0;
          const double s = helper_sigma_nested(ls, e, gv[k], q, ncl);
          sc.push_back(s);
          hlp_closings += ncl;
          hlp_finite += s > HALF_NEG;
          if (ext_check) {
            const double nv = helper_nested_naive(ls, e, gv[k], q);
            ++hc_states;
            if ((nv > HALF_NEG) != (s > HALF_NEG)) ++hc_mis;
            else if (nv > HALF_NEG) hc_rel = std::max(hc_rel, std::fabs(std::expm1(s - nv)));
          }
          continue;
        }
        if (ls) closings_ls(e.p, gv[k], q, cl);
        else closings_rs(e.p, q, cl);
        const double s = helper_sigma(cl, e.u, gv[k], e.a);
        sc.push_back(s);
        hlp_closings += static_cast<long long>(cl.size());
        hlp_finite += s > HALF_NEG;
        if (ext_check) {
          const double nv = helper_naive(cl, e.u, gv[k], e.a);
          ++hc_states;
          if ((nv > HALF_NEG) != (s > HALF_NEG)) ++hc_mis;
          else if (nv > HALF_NEG) hc_rel = std::max(hc_rel, std::fabs(std::expm1(s - nv)));
        }
      }
      hlp_states[ls ? 0 : 1] += static_cast<long long>(m);
      keep_mask(sc, km, beam);
      if (ro_check && stem_profile) ro_flips(so, km, beam);
      for (size_t r = 0; r < m; ++r) {
        keep[sp.first + r] = km[r];
        any |= !km[r];
      }
      if (drops) {
        std::vector<size_t> rank = score_ranks(sc);
        ++drops->events[bt];
        for (size_t r = 0; r < m; ++r)
          if (!km[r]) {
            const Ent& e = *ent[sp.first + r];
            const int v = gv[sp.first + r];
            note_drop(bt, bt, q, gap_step ? -1 : v, e.p, e.u, v, q, rank[r], m, sc[r]);
          } else if (drops->kf) {
            const Ent& e = *ent[sp.first + r];
            const int v = gv[sp.first + r];
            note_kept(bt, q, gap_step ? -1 : v, e.p, e.u, v, q, rank[r], m, sc[r]);
          }
      }
    }
    if (any) filter_step(t, keep);
    const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    hlp_sec += dt;
    if (stem_profile) pf_total_sec += dt;
  }

  // Ranking order of keep_mask: descending, earliest first on ties.
  static std::vector<size_t> order_desc(const std::vector<double>& s) {
    std::vector<size_t> order(s.size());
    for (size_t k = 0; k < s.size(); ++k) order[k] = k;
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return s[a] > s[b]; });
    return order;
  }

  // --stem-est-qg-combine split (DP03._split_mask): for each of parts (sigma_L, sigma_R,
  // sigma_ext) in turn, up to ceil(beam/3) more states by that score in order_desc (finite
  // scores only, chosen states skipped), then up to beam states by the combined score sc.
  static void split_mask(const std::vector<double>* const parts[3], const std::vector<double>& sc,
                         std::vector<char>& keep, int beam) {
    keep.assign(sc.size(), 0);
    const int m = (beam + 2) / 3;
    int taken = 0;
    for (int r = 0; r < 3; ++r) {
      const std::vector<double>& s = *parts[r];
      int got = 0;
      for (size_t k : order_desc(s)) {
        if (got == m || taken == beam || s[k] < HALF_NEG) break;
        if (keep[k]) continue;
        keep[k] = 1;
        ++got;
        ++taken;
      }
    }
    for (size_t k : order_desc(sc)) {
      if (taken == beam) break;
      if (keep[k]) continue;
      keep[k] = 1;
      ++taken;
    }
  }

  double stem_A(int d, int e, int f) const {
    auto it = stA.find(key2(d, e));
    if (it == stA.end()) return NEG_INF;
    const StemA& a = it->second;
    auto k = std::upper_bound(a.fs.begin(), a.fs.end(), f) - a.fs.begin() - 1;
    return k >= 0 && a.fs[k] == f ? a.vals[k] : NEG_INF;
  }

  // A(u+1, v, q) += C[p-1] Qgl at(p, q) for every surviving Qgl state at step q (with ok:
  // only the flattened states k with ok[k]).
  void stem_add_A(int q, const std::vector<char>* ok = nullptr) {
    size_t k = 0;
    for (const Grp& g : Qgl[q].grps)
      for (const Ent& s : g.ents) {
        if (ok && !(*ok)[k++]) continue;
        if (!E.end_ok(s.p, q)) continue;
        double val = prefix_log(s.p) + s.a + E.at(s.p, q);
        uint64_t key = key2(s.u + 1, g.v);
        auto it = stA.find(key);
        if (it == stA.end()) {
          it = stA.emplace(key, StemA{{q}, {val}, {}}).first;
          if (pr_on) it->second.il.emplace_back();
          first_f.set_min(g.v, q);
          std::vector<int>& es = stE[s.u + 1];
          es.insert(std::lower_bound(es.begin(), es.end(), g.v), g.v);
        } else if (it->second.fs.back() == q) {
          it->second.vals.back() = log_add(it->second.vals.back(), val);
        } else {
          it->second.fs.push_back(q);
          it->second.vals.push_back(val);
          if (pr_on) it->second.il.emplace_back();
        }
        if (pr_on) it->second.il.back().emplace_back(s.p, s.a + E.at(s.p, q));
        if (pr_on && (pr_gsel || pr_gcheck))
          prPend.push_back({s.u + 1, g.v, q, s.p, s.a + E.at(s.p, q)});
      }
  }

  // Feasibility and stem score of the flattened states of Qg / Qgl / Qgr at step q
  // (DP03._stem_eval in engine.py).
  struct StemEval {
    DP03& dp;
    int bt, q;
    std::vector<const Grp*> grp;
    std::vector<const Ent*> ent;
    std::vector<char> ok;
    std::vector<double> bval, aval;
    std::unordered_map<uint64_t, double> hcache, hacache;
    double w_pk, c0q, sq, cq;

    StemEval(DP03& d, const Tab4& t, int bt_, int q_) : dp(d), bt(bt_), q(q_) {
      const Energy& E = dp.E;
      dp.stem_cache_step(q);
      w_pk = E.w_beta1 + 2.0 * E.w_stem;
      c0q = dp.stem_->C0[q];
      sq = dp.stem_->S[q + 1];
      cq = dp.hr_c(q);
      for (const Grp& g : t.grps)
        for (const Ent& e : g.ents) {
          grp.push_back(&g);
          ent.push_back(&e);
        }
      const size_t N = ent.size();
      ok.assign(N, 0);
      bval.assign(N, std::nan(""));
      if (bt == BT_QG) {
        for (size_t k = 0; k < N; ++k) {
          int u = ent[k]->u, v = grp[k]->v;
          ok[k] = H(u, v) > HALF_NEG || dp.first_f.query(u + 1, v) <= v - 1;
        }
      } else if (bt == BT_QGL) {
        for (size_t k = 0; k < N; ++k) {
          if (!E.end_ok(ent[k]->p, q)) continue;
          ok[k] = H(ent[k]->u, grp[k]->v) > HALF_NEG || B(k) > HALF_NEG;
        }
      } else {
        aval.assign(N, NEG_INF);
        for (size_t k = 0; k < N; ++k) {
          int p = ent[k]->p, u = ent[k]->u, v = grp[k]->v;
          if (E.end_ok(p, q)) aval[k] = dp.stem_A(p, u + 1, v - 1);
          ok[k] = aval[k] > HALF_NEG;
        }
      }
    }

    double H(int u, int v) {
      uint64_t key = key2(u, v);
      auto it = hcache.find(key);
      if (it != hcache.end()) return it->second;
      double h = dp.stem_->H(u, v, q);
      hcache.emplace(key, h);
      return h;
    }

    // The L-stem score's H: H itself, or Hav with --stem-est-hr av (ok[] keeps H).
    double Hs(int u, int v) {
      if (!dp.hr_av) return H(u, v);
      uint64_t key = key2(u, v);
      auto it = hacache.find(key);
      if (it != hacache.end()) return it->second;
      double h = dp.stem_->Hav(u, v, q);
      hacache.emplace(key, h);
      return h;
    }

    double B(size_t k) {
      if (std::isnan(bval[k])) {
        const int d = ent[k]->p, e = ent[k]->u + 1, v = grp[k]->v;
        bval[k] = dp.stem_qg_partner ? dp.stem_B_cached(d, e, v) : dp.stem_B(d, e, v);
      }
      return bval[k];
    }

    double s1h(size_t k) {
      const Ent& e = *ent[k];
      double h = Hs(e.u, grp[k]->v);
      if (h < HALF_NEG) return NEG_INF;
      return dp.four_score(e.p, e.u, grp[k]->v, e.a) + h + cq + w_pk;
    }

    // Stem-1 score (state as the left stem, H_R); NEG_INF if none.
    double sigL(size_t k) {
      const Energy& E = dp.E;
      const Ent& e = *ent[k];
      if (bt == BT_QGR) return NEG_INF;
      double s = s1h(k);
      if (!dp.stem_consistent || s < HALF_NEG) return s;
      if (bt == BT_QG) return s + (E.at(e.u, grp[k]->v) + E.at(e.p, q));
      return s + (E.at(e.p, q) - E.w_stem);
    }

    // Stem-2 score (state as the right stem); NEG_INF if none or the role is off.
    double sigR(size_t k) {
      const Energy& E = dp.E;
      const Ent& e = *ent[k];
      const int v = grp[k]->v;
      if (bt == BT_QG) {
        if (!dp.stem_qg2) return NEG_INF;
        if (dp.stem_qg_partner) {
          double bg = dp.stem_Bg(e.p, e.u, v);
          if (bg < HALF_NEG) return NEG_INF;
          double s2 = bg + e.a + sq;
          if (dp.stem_consistent) s2 += (E.at(e.p, q) + E.w_beta1);
          return s2;
        }
        double h2 = dp.stem_->H2(e.p, e.u, v);
        if (h2 < HALF_NEG) return NEG_INF;
        double s2 = sq + e.a;
        if (dp.score_qz) s2 += dp.qz_log(e.u + 1, v - 1);
        s2 = s2 + h2 - dp.stem_->S[e.p] + w_pk;
        if (dp.stem_consistent) s2 += (E.at(e.u, v) + E.at(e.p, q));
        return s2;
      }
      if (!dp.stem_both) return NEG_INF;
      if (bt == BT_QGL) {
        double b = B(k);
        return b < HALF_NEG ? NEG_INF : b + e.a + sq + E.at(e.p, q) + E.w_beta1;
      }
      return aval[k] + e.a + sq + E.at(e.p, q) + E.w_beta1;
    }

    double score(size_t k) {
      const Ent& e = *ent[k];
      if (bt == BT_QG) {
        double s = sigL(k);
        if (!dp.stem_qg2) return s;
        double s2 = sigR(k);
        if (s2 < HALF_NEG) return s;
        return log_add(s, s2);
      }
      if (bt == BT_QGL) {
        double s = sigL(k);
        if (!dp.stem_both) return s;
        return log_add(s, sigR(k));
      }
      if (!dp.stem_both) return dp.four_score(e.p, e.u, grp[k]->v, e.a);
      return sigR(k);
    }
  };

  // --stem-est rank: drop infeasible states, rank the rest by the stem score (step: whole
  // table, group: per group v; threshold rule of keep_mask) (DP03._prune_stem).
  void prune4_stem(Tab4& t, int bt, int q, int gtab, int beam) {
    const auto t_qg = std::chrono::steady_clock::now();
    StemEval ev(*this, t, bt, q);
    const size_t N = ev.ent.size();
    std::vector<char> keep(ev.ok);
    size_t ranked = 0;
    if (beam > 0) {
      std::vector<std::pair<size_t, size_t>> spans;
      if (gap_step) {
        spans.emplace_back(0, N);
      } else {
        size_t lo = 0;
        for (const Grp& g : t.grps) {
          spans.emplace_back(lo, lo + g.ents.size());
          lo += g.ents.size();
        }
      }
      if (pr_on && pr_arm_gaponly) {
        std::vector<std::pair<int, int>> iv;
        for (const auto& sp : spans) {
          size_t m = 0;
          for (size_t k = sp.first; k < sp.second; ++k) m += ev.ok[k] ? 1 : 0;
          if (m <= static_cast<size_t>(beam)) continue;
          for (size_t k = sp.first; k < sp.second; ++k)
            if (ev.ok[k]) iv.emplace_back(ev.ent[k]->u + 1, ev.grp[k]->v - 1);
        }
        pr_gap_begin(iv);
      }
      std::vector<size_t> idx;
      std::vector<double> sc, xs, sl, sr, so;
      std::vector<char> km;
      const bool ext = bt == BT_QG && stem_ext;
      for (const auto& sp : spans) {
        idx.clear();
        for (size_t k = sp.first; k < sp.second; ++k)
          if (ev.ok[k]) idx.push_back(k);
        if (idx.size() <= static_cast<size_t>(beam)) continue;
        ranked += idx.size();
        sc.clear();
        if (pr_on) {
          const auto tp = std::chrono::steady_clock::now();
          const int ti = pr_tab(bt);
          sc.assign(idx.size(), NEG_INF);
          for (size_t m : gap_order(idx.size(), [&](size_t m) {
                 return std::make_pair(ev.grp[idx[m]]->v, ev.ent[idx[m]]->u);
               })) {
            const size_t k = idx[m];
            sc[m] = pr_score4(bt, *ev.ent[k], ev.grp[k]->v, q);
            if (pr_dump)
              std::fprintf(pr_dump, "%d %d %d %d %d %.17g\n", bt, q, ev.ent[k]->p, ev.ent[k]->u,
                           ev.grp[k]->v, sc[m]);
            if (pr_check) prc_cmp(ti, sc[m], pr_naive4(bt, *ev.ent[k], ev.grp[k]->v, q));
          }
          pr_states[ti] += static_cast<long long>(idx.size());
          pr_sec[ti] += std::chrono::duration<double>(std::chrono::steady_clock::now() - tp).count();
        } else if (ext && stem_profile) {
          const auto t_x = std::chrono::steady_clock::now();
          sc.assign(idx.size(), NEG_INF);
          if (ro_check) so.assign(idx.size(), NEG_INF);
          for (size_t m : gap_order(idx.size(), [&](size_t m) {
                 return std::make_pair(ev.grp[idx[m]]->v, ev.ent[idx[m]]->u);
               })) {
            const size_t k = idx[m];
            const Ent& e = *ev.ent[k];
            const int v = ev.grp[k]->v;
            const double s = e.a + prof_eval(qg_prof(e.p), e.u, v);
            sc[m] = s;
            if (ro_check) ro_cmp(s, e.a, so, m);
            ext_finite += s > HALF_NEG;
            if (prof_check && stem_chain)
              prof_cmp(0, s, chain_naive({{e.p, stc_q, 0.0}}, e.u, v, e.a));
            else if (prof_check)
              prof_cmp(0, s, log_add(ev.score(k), stem_ext_score(e.p, e.u, v, e.a)));
          }
          const double dt =
              std::chrono::duration<double>(std::chrono::steady_clock::now() - t_x).count();
          ext_sec += dt;
          pf_total_sec += dt;
          ext_states += static_cast<long long>(idx.size());
          pf_states += static_cast<long long>(idx.size());
        } else {
          for (size_t k : idx) sc.push_back(ev.score(k));
        }
        if ((hr_av || hr_check) && bt != BT_QGR)
          for (size_t k : idx) {
            const int u = ev.ent[k]->u, v = ev.grp[k]->v;
            if (hr_check) {
              const double hs = ev.Hs(u, v), nv = hr_naive(u, v, q);
              ++hrc_states;
              if ((hs > HALF_NEG) != (nv > HALF_NEG)) ++hrc_mis;
              else if (nv > HALF_NEG) hrc_rel = std::max(hrc_rel, std::fabs(std::expm1(hs - nv)));
            }
            if (hr_av && bt == BT_QG) {
              const double lo = ev.H(u, v), la = ev.Hs(u, v);
              if (lo > HALF_NEG && la > HALF_NEG) {
                const double d = (la + ev.cq) - (lo - ev.c0q);
                hr_dmin = hr_n == 0 ? d : std::min(hr_dmin, d);
                hr_dmax = hr_n == 0 ? d : std::max(hr_dmax, d);
                hr_dsum += d;
                ++hr_n;
              } else if (la > HALF_NEG) {
                ++hr_av_only;
              } else if (lo > HALF_NEG) {
                ++hr_orig_only;
              }
            }
          }
        if (ext && !stem_profile && !pr_on) {
          const auto t_x = std::chrono::steady_clock::now();
          xs.clear();
          for (size_t k : idx) {
            const Ent& e = *ev.ent[k];
            xs.push_back(stem_ext_score(e.p, e.u, ev.grp[k]->v, e.a));
          }
          ext_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_x).count();
          ext_states += static_cast<long long>(idx.size());
          for (size_t m = 0; m < idx.size(); ++m) {
            ext_finite += xs[m] > HALF_NEG;
            sc[m] = log_add(sc[m], xs[m]);
          }
          if (ext_check)
            for (size_t m = 0; m < idx.size(); ++m) {
              const Ent& e = *ev.ent[idx[m]];
              const double nv = stem_ext_naive(e.p, e.u, ev.grp[idx[m]]->v, e.a);
              ++xc_states;
              if ((nv > HALF_NEG) != (xs[m] > HALF_NEG)) ++xc_mis;
              else if (nv > HALF_NEG) xc_rel = std::max(xc_rel, std::fabs(std::expm1(xs[m] - nv)));
            }
        }
        if (ext && stem_split && !pr_on) {
          sl.clear();
          sr.clear();
          for (size_t k : idx) {
            sl.push_back(ev.sigL(k));
            sr.push_back(ev.sigR(k));
          }
          const std::vector<double>* parts[3] = {&sl, &sr, &xs};
          split_mask(parts, sc, km, beam);
        } else {
          keep_mask(sc, km, beam);
          if (ro_check && ext && stem_profile) ro_flips(so, km, beam);
        }
        for (size_t m = 0; m < idx.size(); ++m) keep[idx[m]] = km[m];
        if (drops) {
          std::vector<size_t> rank = score_ranks(sc);
          ++drops->events[bt];
          for (size_t m = 0; m < idx.size(); ++m)
            if (!km[m]) {
              const Ent& e = *ev.ent[idx[m]];
              const int v = ev.grp[idx[m]]->v;
              note_drop(bt, bt, q, gap_step ? -1 : v, e.p, e.u, v, q, rank[m], idx.size(), sc[m]);
            } else if (drops->kf) {
              const Ent& e = *ev.ent[idx[m]];
              const int v = ev.grp[idx[m]]->v;
              note_kept(bt, q, gap_step ? -1 : v, e.p, e.u, v, q, rank[m], idx.size(), sc[m]);
            }
        }
      }
    }
    if (bt == BT_QG)
      stem_qg_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_qg).count();
    size_t n_ok = 0, n_keep = 0;
    for (size_t k = 0; k < N; ++k) {
      n_ok += ev.ok[k];
      n_keep += keep[k];
    }
    stem_removed[bt == BT_QG ? 0 : bt == BT_QGL ? 1 : 2] += N - n_ok;
    if (gtab >= 0) {
      gate_stats.cand[gtab] += N;
      gate_stats.elig[gtab] += n_ok;
      gate_stats.ranked[gtab] += ranked;
      gate_stats.kept[gtab] += n_keep;
    }
    if (n_keep < N) filter_step(t, keep);
    if (bt == BT_QGL) stem_add_A(q);
  }

  // Gate selection over `n_all` candidates (one group, or the whole table in step mode):
  // (a) drop g < tau; (b) if more than `beam` remain, keep the top `beam` of them by the
  // pruning score (threshold rule of keep_mask); (c) add up to `reserve` more eligible
  // states by g, descending, ties by insertion order.  Returns false if every candidate
  // is kept (the caller then leaves the states untouched).
  template <class ScoreAt>
  bool gate_select(const std::vector<double>& gv, ScoreAt score_at, std::vector<char>& keep,
                   int gtab, int beam) {
    const size_t n_all = gv.size();
    std::vector<char> elig(n_all, 0);
    size_t ne = 0;
    for (size_t k = 0; k < n_all; ++k) {
      elig[k] = gv[k] >= gcfg.tau;
      ne += elig[k];
    }
    gate_stats.cand[gtab] += n_all;
    gate_stats.elig[gtab] += ne;
    if (beam <= 0 || ne <= static_cast<size_t>(beam)) {
      gate_stats.kept[gtab] += ne;
      if (ne == n_all) return false;
      keep.swap(elig);
      return true;
    }
    std::vector<double> scores;
    scores.reserve(ne);
    for (size_t k = 0; k < n_all; ++k)
      if (elig[k]) scores.push_back(score_at(k));
    gate_stats.ranked[gtab] += ne;
    std::vector<char> km;
    keep_mask(scores, km, beam);
    keep.assign(n_all, 0);
    size_t m = 0, kept = 0;
    for (size_t k = 0; k < n_all; ++k)
      if (elig[k]) {
        keep[k] = km[m++];
        kept += keep[k];
      }
    if (gcfg.reserve > 0) {
      std::vector<size_t> rest;
      for (size_t k = 0; k < n_all; ++k)
        if (elig[k] && !keep[k]) rest.push_back(k);
      std::stable_sort(rest.begin(), rest.end(),
                       [&](size_t a, size_t b) { return gv[a] > gv[b]; });
      size_t add = std::min(rest.size(), static_cast<size_t>(gcfg.reserve));
      for (size_t r = 0; r < add; ++r) keep[rest[r]] = 1;
      kept += add;
      gate_stats.reserved[gtab] += add;
    }
    gate_stats.kept[gtab] += kept;
    return true;
  }

  void prune4_gated(Tab4& t, int q, bool use_arm, int gtab, int beam) {
    const Gate& G = *gate_;
    std::vector<double> gv;
    std::vector<char> keep;
    if (!gap_step) {
      for (auto& g : t.grps) {
        gv.clear();
        for (const Ent& e : g.ents) gv.push_back(G.g(e.p, e.u, g.v, q));
        auto score_at = [&](size_t k) {
          const Ent& e = g.ents[k];
          return use_arm ? arm_score_of(e.p, e.u, g.v, q, e.a) : four_score(e.p, e.u, g.v, e.a);
        };
        if (gate_select(gv, score_at, keep, gtab, beam)) filter_group(g, keep, 0);
      }
      return;
    }
    std::vector<std::pair<int, int>> where;  // (group, entry) of each flattened candidate
    for (int gi = 0; gi < static_cast<int>(t.grps.size()); ++gi) {
      const Grp& g = t.grps[gi];
      for (int r = 0; r < static_cast<int>(g.ents.size()); ++r) {
        const Ent& e = g.ents[r];
        gv.push_back(G.g(e.p, e.u, g.v, q));
        where.emplace_back(gi, r);
      }
    }
    auto score_at = [&](size_t k) {
      const Grp& g = t.grps[where[k].first];
      const Ent& e = g.ents[where[k].second];
      return use_arm ? arm_score_of(e.p, e.u, g.v, q, e.a) : four_score(e.p, e.u, g.v, e.a);
    };
    if (gate_select(gv, score_at, keep, gtab, beam)) filter_step(t, keep);
  }

  // Top `beam` of `scores`; ties at the threshold go to the earliest (DP03._keep_mask).
  static void keep_mask(const std::vector<double>& scores, std::vector<char>& keep, int beam) {
    std::vector<double> tmp(scores);
    std::nth_element(tmp.begin(), tmp.begin() + (beam - 1), tmp.end(), std::greater<double>());
    const double thr = tmp[beam - 1];
    int room = beam;
    for (double s : scores) room -= s > thr;
    keep.assign(scores.size(), 0);
    for (size_t k = 0; k < scores.size(); ++k) {
      if (scores[k] > thr) {
        keep[k] = 1;
      } else if (scores[k] == thr && room > 0) {
        keep[k] = 1;
        --room;
      }
    }
  }

  // Keeps g.ents[r] with keep[off + r], in insertion order; rebuilds the index.
  // Releases the pre-prune capacity: a group can hold many times the beam
  // before pruning and would otherwise keep that memory for the whole run.
  static void filter_group(Grp& g, const std::vector<char>& keep, size_t off) {
    std::vector<Ent> ents;
    for (size_t r = 0; r < g.ents.size(); ++r)
      if (keep[off + r]) ents.push_back(g.ents[r]);
    ents.shrink_to_fit();
    g.ents.swap(ents);
    std::unordered_map<uint64_t, int> idx;
    idx.reserve(g.ents.size());
    for (int r = 0; r < static_cast<int>(g.ents.size()); ++r)
      idx.emplace(key2(g.ents[r].p, g.ents[r].u), r);
    g.idx.swap(idx);
  }

  // Top `beam` states of the whole table at this step (DP03._prune, gap_beam="step").
  // Kept states stay in insertion order; emptied groups are dropped.
  void prune4_step(Tab4& t, int bt, int q, bool use_arm, int gtab, int beam) {
    size_t total = 0;
    for (const Grp& g : t.grps) total += g.ents.size();
    if (gtab >= 0) {
      gate_stats.cand[gtab] += total;
      gate_stats.elig[gtab] += total;
    }
    if (total <= static_cast<size_t>(beam)) {
      if (gtab >= 0) gate_stats.kept[gtab] += total;
      return;
    }
    std::vector<double> scores;
    scores.reserve(total);
    for (const Grp& g : t.grps)
      for (const Ent& e : g.ents)
        scores.push_back(use_arm ? arm_score_of(e.p, e.u, g.v, q, e.a)
                                 : four_score(e.p, e.u, g.v, e.a));
    std::vector<char> keep;
    keep_mask(scores, keep, beam);
    if (gtab >= 0) {
      gate_stats.ranked[gtab] += total;
      gate_stats.kept[gtab] += std::min(total, static_cast<size_t>(beam));
    }
    if (drops) {
      std::vector<size_t> rank = score_ranks(scores);
      ++drops->events[bt];
      size_t k = 0;
      for (const Grp& g : t.grps)
        for (const Ent& e : g.ents) {
          if (!keep[k])
            note_drop(bt, bt, q, -1, e.p, e.u, g.v, q, rank[k], total, scores[k]);
          ++k;
        }
    }
    filter_step(t, keep);
  }

  // Rank of each score (0 = best) in keep_mask's order: descending, earliest first on ties.
  static std::vector<size_t> score_ranks(const std::vector<double>& s) {
    std::vector<size_t> order(s.size()), rank(s.size());
    for (size_t k = 0; k < s.size(); ++k) order[k] = k;
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return s[a] > s[b]; });
    for (size_t k = 0; k < order.size(); ++k) rank[order[k]] = k;
    return rank;
  }

  // One dropped candidate (table bt, dropped by the ranking of table src at step j; v = its
  // group in group mode, else -1).  With an exact state file only its states are written.
  void note_drop(int bt, int src, int j, int v, int a, int b, int c, int d, size_t rank,
                 size_t cands, double score) {
    DropLog& L = *drops;
    ++L.dropped[bt];
    double w = std::nan("");
    if (L.exact) {
      auto it = L.exact->find(state_key(bt, a, b, c, d));
      if (it == L.exact->end()) return;
      w = it->second;
      ++L.matched[bt];
      L.w[bt] += w;
      L.wpos[bt][std::min(9, 10 * j / L.n)] += w;
      if (w < L.min_w) return;
    }
    std::fprintf(L.f, "D %s %s %s %d %d %d %d %d %d %zu %zu %.9g %.6e\n", L.id.c_str(),
                 BT_NAMES[bt], BT_NAMES[src], j, v, a, b, c, d, rank, cands, score, w);
  }

  // Pass-2 / l2r table drop (R,U,S,P1-P3,…): same D line shape; table/src are l2r names.
  void note_l2r_drop(const char* tname, int j, int v, int a, int b, int c, int d, size_t rank,
                     size_t cands, double score) {
    if (!drops) return;
    std::fprintf(drops->f, "D %s %s %s %d %d %d %d %d %d %zu %zu %.9g nan\n", drops->id.c_str(),
                 tname, tname, j, v, a, b, c, d, rank, cands, score);
  }

  // --dump-kept: one kept candidate of a stem / helper ranking event (fields of note_drop;
  // w = exact w, nan if not in the exact state file or without one).
  void note_kept(int bt, int j, int v, int a, int b, int c, int d, size_t rank, size_t cands,
                 double score) {
    DropLog& L = *drops;
    double w = std::nan("");
    if (L.exact) {
      auto it = L.exact->find(state_key(bt, a, b, c, d));
      if (it != L.exact->end()) w = it->second;
    }
    std::fprintf(L.kf, "K %s %s %s %d %d %d %d %d %d %zu %zu %.9g %.6e\n", L.id.c_str(),
                 BT_NAMES[bt], BT_NAMES[bt], j, v, a, b, c, d, rank, cands, score, w);
  }

  // Keeps the flattened states with keep[k] (groups in order, entries in insertion order);
  // emptied groups are dropped.
  void filter_step(Tab4& t, const std::vector<char>& keep) {
    std::vector<Grp> grps;
    t.vidx.clear();
    size_t off = 0;
    for (Grp& g : t.grps) {
      size_t sz = g.ents.size();
      filter_group(g, keep, off);
      off += sz;
      if (g.ents.empty()) continue;
      t.vidx.emplace(g.v, static_cast<int>(grps.size()));
      grps.push_back(std::move(g));
    }
    t.grps.swap(grps);
  }

  // Beam on a right-end map {left: value} by prefix score (DP03._prune_1d).
  void prune1(OMap& at, int j) {
    const int beam = tbeam[BT_PKM];
    if (remove_keys) {
      for (size_t r = 0; r < at.slots(); ++r) {
        int d = at.keys[r];
        if (d >= 0 && remove_keys->count(state_key(BT_PKM, d, j, -1, -1))) at.erase(d);
      }
      at.compact();
    }
    if (beam <= 0 || at.live <= beam) return;
    std::vector<std::pair<double, int>> ranked;
    ranked.reserve(at.live);
    for (size_t r = 0; r < at.slots(); ++r) {
      int i = at.keys[r];
      if (i < 0) continue;
      ranked.emplace_back(prefix_log(i) + at.vals[r], i);
    }
    std::sort(ranked.begin(), ranked.end(),
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                if (a.first != b.first) return a.first > b.first;
                return a.second > b.second;
              });
    std::vector<char> keep(n, 0);
    for (int k = 0; k < beam; ++k) keep[ranked[k].second] = 1;
    if (drops) {
      ++drops->events[BT_PKM];
      for (size_t k = beam; k < ranked.size(); ++k)
        note_drop(BT_PKM, BT_PKM, j, -1, ranked[k].second, j, -1, -1, k, ranked.size(),
                  ranked[k].first);
    }
    for (size_t r = 0; r < at.slots(); ++r) {
      int i = at.keys[r];
      if (i >= 0 && !keep[i]) at.erase(i);
    }
    at.compact();
  }

  // extra_row / extra_at: Qm2, dropped with Qm (bt = BT_QM).
  void prune2(AtView& table, OMap& at, int j, int bt, AtView* extra_row, OMap* extra_at) {
    const int beam = tbeam[bt];
    if (remove_keys) {
      for (size_t r = 0; r < at.slots(); ++r) {
        int i = at.keys[r];
        if (i < 0 || !remove_keys->count(state_key(bt, i, j, -1, -1))) continue;
        at.erase(i);
        table[ix(i, j)] = NEG_INF;
      }
      at.compact();
    }
    if (beam <= 0 || at.live <= beam) {
      if (g_beam_stats.f && beam > 0)
        g_beam_stats.log_prune(BT_NAMES[bt], j, static_cast<size_t>(at.live),
                               static_cast<size_t>(at.live), 0.0, 0);
      return;
    }
    std::vector<std::pair<double, int>> ranked;
    ranked.reserve(at.live);
    // --stem-est-principled: Qp by alpha + log Ctx(i, j); -nested: Qb / Qm / Qm2 by alpha +
    // their nested outside
    const std::vector<double>* nb = nullptr;
    if (pr_nested && stem_)
      nb = bt == BT_QB ? &stem_->prQb_b : bt == BT_QM ? &stem_->prQm_b
                                         : bt == BT_QM2 ? &stem_->prQm2_b : nullptr;
    const bool pq = pr_on && bt == BT_QP;
    const auto tp = std::chrono::steady_clock::now();
    for (size_t r = 0; r < at.slots(); ++r) {
      int i = at.keys[r];
      if (i < 0) continue;
      if (pq) {
        const double s = at.vals[r] + pr_ctx(i, j);
        ranked.emplace_back(s < HALF_NEG ? NEG_INF : s, i);
        if (pr_check) {
          const long double c = prn_ctx(i, j);
          prc_cmp(5, ranked.back().first,
                  c < HALF_NEG ? NEG_INF : static_cast<double>(at.vals[r] + c));
        }
      } else if (nb) {
        const double s = at.vals[r] + (*nb)[ix(i, j)];
        ranked.emplace_back(s < HALF_NEG ? NEG_INF : s, i);
      } else if (pre_key_inside) {
        ranked.emplace_back(at.vals[r], i);
      } else {
        ranked.emplace_back(prefix_log(i) + at.vals[r], i);
      }
    }
    if (pq) {
      pr_states[5] += at.live;
      pr_sec[5] += std::chrono::duration<double>(std::chrono::steady_clock::now() - tp).count();
    }
    std::sort(ranked.begin(), ranked.end(),
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                if (a.first != b.first) return a.first > b.first;
                return a.second > b.second;
              });
    if (g_beam_stats.f) {
      std::vector<double> keys;
      keys.reserve(ranked.size());
      for (const auto& e : ranked) keys.push_back(e.first);
      const size_t kept_n = static_cast<size_t>(beam);
      g_beam_stats.log_prune(BT_NAMES[bt], j, ranked.size(), kept_n,
                             BeamStats::mass_share(keys, kept_n), 1);
    }
    std::vector<char> keep(n, 0);
    for (int k = 0; k < beam; ++k) keep[ranked[k].second] = 1;
    if (drops) {
      ++drops->events[bt];
      for (size_t k = beam; k < ranked.size(); ++k) {
        const int i = ranked[k].second;
        note_drop(bt, bt, j, -1, i, j, -1, -1, k, ranked.size(), ranked[k].first);
        if (extra_at && extra_at->has(i))
          note_drop(BT_QM2, bt, j, -1, i, j, -1, -1, k, ranked.size(), ranked[k].first);
      }
    }
    for (size_t r = 0; r < at.slots(); ++r) {
      int i = at.keys[r];
      if (i < 0 || keep[i]) continue;
      at.erase(i);
      table[ix(i, j)] = NEG_INF;
      if (extra_row) (*extra_row)[ix(i, j)] = NEG_INF;
      if (extra_at) extra_at->erase(i);
    }
    at.compact();
    if (extra_at) extra_at->compact();
  }

  void set_pair(AtView& table, std::vector<OMap>& at, int i, int j, double val) {
    if (val < HALF_NEG) return;
    double t = log_add(static_cast<double>(table[ix(i, j)]), val);
    table[ix(i, j)] = t;
    at[j].set(i, t);
  }

  void set_Qm(int i, int j, double val, bool is_m2) {
    if (val < HALF_NEG) return;
    if (dz) {
      rsQm[j & 1].P.add(i, val);
      if (is_m2) rsQm2[j & 1].P.add(i, val);
      return;
    }
    double t = log_add(static_cast<double>(Qm[ix(i, j)]), val);
    Qm[ix(i, j)] = t;
    Qm_at[j].set(i, t);
    if (is_m2) {
      double t2 = log_add(static_cast<double>(Qm2[ix(i, j)]), val);
      Qm2[ix(i, j)] = t2;
      Qm2_at[j].set(i, t2);
    }
  }

  // PKM[j] = PKM[j-1] + one unpaired, plus Qp ending at j; then beam.
  void advance_pk_inside(int j) {
    if (!allow_pk) return;
    const double w3 = E.w_alpha3;
    OMap& at = PKM_at[j];
    if (dz) {
      RunSt& cu = rsPKM[j & 1];
      cu.clear();
      std::vector<int> surv;
      if (j > 0) {
        surv = live_keys(PKM_at[j - 1]);
        run_ext(rsPKM[(j - 1) & 1], cu, surv, j, w3);
      }
      const OMap& qp = Qp_at[j];
      for (size_t r = 0; r < qp.slots(); ++r)
        if (qp.keys[r] >= 0) cu.P.put(qp.keys[r], qp.vals[r]);
      run_fin(nullptr, at, cu, surv, j);
      prune1(at, j);
      return;
    }
    if (j > 0) {
      const OMap& prev = PKM_at[j - 1];
      for (size_t r = 0; r < prev.slots(); ++r) {
        int d = prev.keys[r];
        if (d >= 0) at.set(d, prev.vals[r] + w3);
      }
    }
    const OMap& qp = Qp_at[j];
    for (size_t r = 0; r < qp.slots(); ++r) {
      int d = qp.keys[r];
      if (d >= 0) at.set(d, log_add(at.get(d), qp.vals[r]));
    }
    prune1(at, j);
  }

  // ============================================ --pk-l2r: H-type pseudoknots in scan order
  // Qp(a', ys) = sum A F1 F2 F3 B (stemest/l2r/r_stage0.md) with the partner's pieces as states:
  //   A (LA): seeds Gs(a2, xs, x', c1) = pb + w_stem + at(a2, c1) + F1 + F2, grown G(a, xs, x', j);
  //     Tab4 group x', entry (p = a, u = xs).
  //   S(a', xs, x', j) = sum_q' [Afull(q') + at(a', q')] Qz(q'+1, j) through the Qz run automaton
  //     (P, U1, U2; c = the close-A term at j; F = closed + c).
  //   B (LB): seeds Gs(a', xs, x', y') = S_F(y'-1) + w_stem + at(x', y') + pb, grown G(a', xs, x, y),
  //     x >= xs; Tab4 group xs, entry (p = x, u = a').  Qp(a', y) += Bfull(xs, y) + at(xs, y).
  //   Gls / Glss (seed-derived) / Grs: Qgls / Qgrs with the carried keys.
  // l2r_mode 1: Qp from the recursion; 2 (check): both, the engine's Qp feeds the nested tables.
  int l2r_mode = 0;
  bool l2r_check_R = false;
  struct LFam {
    bool isU = false;
    std::vector<Tab4> G, Gs, Gls, Glss, Grs, Glt;
    // --grammar simple: M, M2, Multi. M2 groups with v >= n are M5 (stem then nested);
    // v < n are M4 (nested then stem). The logical group key is v % n.
    std::vector<Tab4> SM, SM2, SMulti;
  };
  LFam LA, LB;
  bool grammar_simple = false;
  bool c_trace = false;
  long long gs_fill = 0, gs_phase = 0, gs_lead = 0, gs_over = 0, gs_m[8] = {};
  // Kept after each type's own beam: R^M, R^M2, R^Multi, U^M, U^M2, U^Multi.
  long long gs_tab_pre[6] = {}, gs_tab_kept[6] = {}, gs_tab_max[6] = {};
  struct LeadCol {
    std::vector<double> w;
  };
  std::vector<LeadCol> gsLead;
  std::unordered_set<uint64_t> gs_seedA, gs_seedB;
  struct SSt {
    double P, U1, U2, c, F;
    double bP, bU1, bU2, bc, bF;
  };
  struct SMap {
    std::vector<uint64_t> keys;
    std::vector<SSt> v;
    std::unordered_map<uint64_t, int> idx;
    int at(uint64_t k) {
      auto it = idx.find(k);
      if (it != idx.end()) return it->second;
      idx.emplace(k, static_cast<int>(keys.size()));
      keys.push_back(k);
      v.push_back(SSt{NEG_INF, NEG_INF, NEG_INF, NEG_INF, NEG_INF, NEG_INF, NEG_INF, NEG_INF,
                      NEG_INF, NEG_INF});
      return static_cast<int>(keys.size()) - 1;
    }
    int find(uint64_t k) const {
      auto it = idx.find(k);
      return it == idx.end() ? -1 : it->second;
    }
  };
  std::vector<SMap> LS;
  std::vector<double> LQp;
  std::vector<int> lmaxR;  // largest y with pairs(x, y) and end_ok(x, y), -1 if none
  enum { LT_QZ, LT_QM, LT_QB, LT_QP, LT_AG, LT_AGS, LT_ALS, LT_ARS, LT_BG, LT_BGS, LT_BLS,
         LT_BRS, LT_S, LT_P1, LT_P2, LT_P3, LT_N };
  // --l2r-spawn inc: the A seeds built in scan order from three filler tables. The column
  // is the right end, so the index counts are P1(a, t) = 2, P2(a, b, t) = 3, P3(a, b, c, t) = 4.
  // key3 packs the indices other than t; unused slots stay 0.
  //   P1(a2, j) = Qz(a2+1, j), opened with c = 1 at j = a2, extended as S (fill_Qz automaton);
  //   P2(a2, xs, j) = P1(a2, xs-1) carried (value in F): B's arm xs.. is weighed later by U;
  //   P3(a2, xs, x', j) = P2(a2, xs, x') Qz(x'+1, j), opened with c = P2 at j = x', extended as S;
  //   Rs(a2, xs, x', c1) = P3(a2, xs, x', c1-1) + pb + w_stem + at(a2, c1) at a pair a2.c1.
  bool l2r_inc = false;
  std::vector<SMap> LP1, LP2, LP3;
  int l2r_done[LT_N];
  long long l2r_tr[LT_N] = {};
  // l2r_grow edges that land in R or U: 0..3 R interior/(f)/(g)/(h), 4..7 the same for U.
  int ext_slot = -1;
  long long ext_n[8] = {};
  // Loop iterations, including ones that fail the pair check. [0] R, [1] U.
  long long it_int[2] = {}, it_h[2] = {}, it_ls[2] = {}, it_rs[2] = {};
  // Cumulative counters at the end of each scan column, for the last-20% means.
  struct ColCnt {
    long long iti, ith, itls, itrs, tail, edge, qz, k, qm;
  };
  std::vector<ColCnt> col_cnt;
  bool next_pair = false;
  void note_col() {
    long long edge = 0;
    for (int s = 0; s < 8; ++s) edge += ext_n[s];
    col_cnt.push_back({it_int[0] + it_int[1], it_h[0] + it_h[1], it_ls[0] + it_ls[1],
                       it_rs[0] + it_rs[1], rlt_ext, edge, qz_reads, sx_k_n, qm_kept_n});
  }
  // Qp after the per-position beam (preset beam). pre is the column before prune2.
  long long qp_pre_n = 0, qp_kept_n = 0, qp_bind_n = 0;
  int qp_pre_max = 0, qp_kept_max = 0, qp_beam = 0;
  // Qm keys after the per-column beam.
  long long qm_pre_n = 0, qm_kept_n = 0, qm_bind_n = 0;
  int qm_pre_max = 0, qm_kept_max = 0, qm_beam = 0;
  // --rlt: case (f) right flank kept as a one-base tail. rlt_beam 0 does not prune.
  bool rlt = false;
  int rlt_beam = 0;
  long long rlt_kept = 0;
  int rlt_max = 0;
  long long rlt_ext = 0;
  long long rlt_kept_f[2] = {};
  int rlt_max_f[2] = {};
  std::vector<int> rlt_pos[2];
  // seconds: A grow+helpers, S, B grow+helpers, closes, spawn, checks
  double l2r_sec[6] = {};
  long long lc_R_n = 0, lc_R_mis = 0, lc_Qp_n = 0, lc_Qp_mis = 0;
  double lc_R_rel = 0.0, lc_Qp_rel = 0.0;

  static int k3a(uint64_t k) { return static_cast<int>(k >> 42); }
  static int k3b(uint64_t k) { return static_cast<int>((k >> 21) & 0x1FFFFF); }
  static int k3c(uint64_t k) { return static_cast<int>(k & 0x1FFFFF); }

  void l2r_init() {
    LA.isU = false;
    LB.isU = true;
    for (LFam* F : {&LA, &LB}) {
      F->G.assign(n, Tab4{});
      F->Gs.assign(n, Tab4{});
      F->Gls.assign(n, Tab4{});
      F->Glss.assign(n, Tab4{});
      F->Grs.assign(n, Tab4{});
      if (rlt) F->Glt.assign(n, Tab4{});
      if (grammar_simple) {
        F->SM.assign(n, Tab4{});
        F->SM2.assign(n, Tab4{});
        F->SMulti.assign(n, Tab4{});
      }
    }
    if (grammar_simple) gsLead.assign(n, LeadCol{});
    LS.assign(n, SMap{});
    if (l2r_inc) {
      LP1.assign(n, SMap{});
      LP2.assign(n, SMap{});
      LP3.assign(n, SMap{});
    }
    if (l2r_mode == 2) LQp.assign(static_cast<size_t>(n) * n, NEG_INF);
    ensure_pairs();
    lmaxR.assign(n, -1);
    for (int x = 0; x < n; ++x)
      for (int y : partners_right[x])
        if (E.end_ok(x, y)) lmaxR[x] = std::max(lmaxR[x], y);
    std::fill(l2r_done, l2r_done + LT_N, -1);
    qz_audit = qz_audit_req;
  }

  void lchk(int t, int pos) const {
    if (pos <= l2r_done[t]) return;
    static const char* nm[LT_N] = {"Qz", "Qm", "Qb", "Qp", "R", "Rs", "Rls", "Rrs",
                                   "U",  "Us", "Uls", "Urs", "S", "P1", "P2", "P3"};
    std::fprintf(stderr, "[l2r] abort: read of %s at %d, finished through %d (scan %d)\n", nm[t],
                 pos, l2r_done[t], cur_j);
    std::abort();
  }
  // B with left arm xs..x' can still close right of j
  bool l2r_arm(int xs, int xp, int j) const { return lmaxR[xs] > j && lmaxR[xp] > j; }
  void lacc(Tab4& t, int lt, int v, int p, int u, double val) {
    if (val < HALF_NEG) return;
    ++l2r_tr[lt];
    if (ext_slot >= 0) ++ext_n[ext_slot];
    t.acc(v, p, u, val);
  }
  // Every transition of the grow / helper tables is target (v, p, u) += s.a + w (+ Qm at m
  // when qm_i >= 0).  l2r_adj: the same enumeration pushes the target's outside value back.
  bool l2r_adj = false;
  void ledge(Tab4& T, int lt, int v, int p, int u, Ent& s, double w, int qm_i = -1,
             int qm_j = -1, double qm = 0.0) {
    if (!l2r_adj) {
      lacc(T, lt, v, p, u, s.a + (w + qm));
      return;
    }
    Ent* t = T.find(v, p, u);
    if (!t || t->b < HALF_NEG) return;
    add4(s, t->b + (w + qm));
    if (qm_i >= 0) add2b(qmp_b(), qm_i, qm_j, t->b + s.a + w);
  }

  // Tail key packs the run's left end s with the carried u (xs for R, a' for U).
  int tail_pack(int s, int u) const { return s * n + u; }
  void tail_unpack(int packed, int& s, int& u) const {
    s = packed / n;
    u = packed - s * n;
  }
  // Log weight added when the unpaired flank's right end becomes `right` (s == right: a new
  // one-base run). Dangles some: one base is s1, a longer run is d5(left) + d3(right).
  double tail_extend(int s, int right) const {
    if (!dz) return E.w_cp;
    if (s == right) return E.w_cp + E.s1w[right];
    if (right - s == 1) return E.w_cp + (E.d5w[s] + E.d3w[right]) - E.s1w[s];
    return E.w_cp + E.d3w[right] - E.d3w[right - 1];
  }
  void l2r_prune_by_inside(Tab4& t, int beam) {
    size_t tot = 0;
    for (const Grp& g : t.grps) tot += g.ents.size();
    if (beam <= 0 || tot <= static_cast<size_t>(beam)) return;
    std::vector<double> sc;
    sc.reserve(tot);
    for (const Grp& g : t.grps)
      for (const Ent& e : g.ents) sc.push_back(e.a);
    std::vector<char> keep;
    keep_mask(sc, keep, beam);
    filter_step(t, keep);
  }
  // Rank by C(p-1) + inside. Keys only; the stored inside is unchanged.
  void l2r_prune_prefix(Tab4& t, int b, int slot = -1) {
    size_t tot = 0;
    for (const Grp& g : t.grps) tot += g.ents.size();
    size_t kept = tot;
    if (b > 0 && tot > static_cast<size_t>(b)) {
      std::vector<double> sc;
      sc.reserve(tot);
      for (const Grp& g : t.grps)
        for (const Ent& e : g.ents) sc.push_back(prefix_log(e.p) + e.a);
      std::vector<char> keep;
      keep_mask(sc, keep, b);
      filter_step(t, keep);
      kept = 0;
      for (const Grp& g : t.grps) kept += g.ents.size();
    }
    if (slot >= 0) {
      gs_tab_pre[slot] += static_cast<long long>(tot);
      gs_tab_kept[slot] += static_cast<long long>(kept);
      if (static_cast<long long>(kept) > gs_tab_max[slot])
        gs_tab_max[slot] = static_cast<long long>(kept);
    }
  }
  // Same score as Rls / Uls (kind 1, not a seed). The packed run start is not part of the key.
  void l2r_prune_tail(LFam& F, Tab4& t, int j) {
    size_t tot = 0;
    for (const Grp& g : t.grps) tot += g.ents.size();
    if (rlt_beam <= 0 || tot <= static_cast<size_t>(rlt_beam)) return;
    std::vector<double> sc;
    sc.reserve(tot);
    for (const Grp& g : t.grps)
      for (const Ent& e : g.ents) {
        int run_s, uu;
        tail_unpack(e.u, run_s, uu);
        (void)run_s;
        if (!F.isU) {
          const RN& r = rec_key(true, e.p, false, j);
          const int arm = larm(uu, g.v);
          const double s = lsc_L(r, arm, 1);
          sc.push_back(s > HALF_NEG ? e.a + s + lkb_arm(e.p, arm) : NEG_INF);
        } else {
          double s;
          if (g_rfv) s = rfv_key(true, e.p, false, j, g.v);
          else s = lsc_R(rec_key(true, e.p, false, j), g.v);
          const double add = prefix_log(uu) + (lsub_c ? lrho(uu) : lkb_xs(uu, g.v));
          sc.push_back(s > HALF_NEG ? e.a + s + add + (lLam.empty() ? 0.0 : lLam[uu]) : NEG_INF);
        }
      }
    std::vector<char> keep;
    keep_mask(sc, keep, rlt_beam);
    filter_step(t, keep);
  }
  // Rlt/Ult at j: Rls/Uls at j-1 plus one unpaired base, and every previous tail plus one base.
  void l2r_tail(LFam& F, int j) {
    if (!rlt || j < 1) return;
    Tab4& dst = F.Glt[j];
    const double add_new = tail_extend(j, j);
    for (const Grp& g : F.Gls[j - 1].grps)
      for (const Ent& e : g.ents) {
        ++rlt_ext;
        dst.acc(g.v, e.p, tail_pack(j, e.u), e.a + add_new);
      }
    for (const Grp& g : F.Glt[j - 1].grps)
      for (const Ent& e : g.ents) {
        ++rlt_ext;
        int s, uu;
        tail_unpack(e.u, s, uu);
        dst.acc(g.v, e.p, e.u, e.a + tail_extend(s, j));
      }
    if (rlt_beam > 0) l2r_prune_tail(F, dst, j);
    size_t tot = 0;
    for (const Grp& g : dst.grps) tot += g.ents.size();
    const int fam = F.isU ? 1 : 0;
    if (rlt_pos[fam].size() != static_cast<size_t>(n)) rlt_pos[fam].assign(n, 0);
    rlt_pos[fam][j] = static_cast<int>(tot);
    rlt_kept_f[fam] += static_cast<long long>(tot);
    rlt_max_f[fam] = std::max(rlt_max_f[fam], static_cast<int>(tot));
    rlt_kept += static_cast<long long>(tot);
    rlt_max = std::max(rlt_max, static_cast<int>(tot));
  }
  // Adjoint of l2r_tail(j): Glt[j] back into Gls[j-1] (new one-base runs) and Glt[j-1].
  void l2r_tail_back(LFam& F, int j) {
    if (!rlt || j < 1) return;
    const double add_new = tail_extend(j, j);
    for (Grp& g : F.Glt[j].grps)
      for (Ent& e : g.ents) {
        if (e.b < HALF_NEG) continue;
        int s, uu;
        tail_unpack(e.u, s, uu);
        Tab4& src_tab = s == j ? F.Gls[j - 1] : F.Glt[j - 1];
        Grp* src = src_tab.get(g.v);
        if (!src) continue;
        const uint64_t key = s == j ? key2(e.p, uu) : key2(e.p, e.u);
        auto it = src->idx.find(key);
        if (it == src->idx.end()) continue;
        add4(src->ents[it->second], e.b + (s == j ? add_new : tail_extend(s, j)));
      }
  }
  // Case (f) from the tail at j-1, the empty run on Rls/Rlss at j-1, and Rlss for runs 1..cap.
  // Rlss is not copied into the tail, so a seed derivation is closed once.
  void l2r_close_tails(LFam& F, int j, int tg) {
    const int o = F.isU ? 4 : 0;
    ext_slot = (F.isU ? 4 : 0) + 1;
    auto close = [&](Grp& g, Ent& s, int uu, double w_pre) {
      const int i = s.p - 1;
      const int lo = F.isU ? g.v : 0;
      if (i < lo || !E.pairs(i, j) || !E.end_ok(i, j)) return;
      if (!F.isU && !l2r_arm(uu, g.v, j)) return;
      ledge(F.G[j], tg, g.v, i, uu, s, w_pre + E.at(i, j));
    };
    if (j < 1) return;
    lchk(LT_ALS + o, j - 1);
    const double w_empty = (E.w_ap + E.w_bp) + pair_bonus;
    for (Grp& g : F.Glt[j - 1].grps)
      for (Ent& s : g.ents) {
        int run_s, uu;
        tail_unpack(s.u, run_s, uu);
        (void)run_s;
        close(g, s, uu, w_empty);
      }
    for (int src = 0; src < 2; ++src)
      for (Grp& g : (src ? F.Glss[j - 1] : F.Gls[j - 1]).grps)
        for (Ent& s : g.ents) close(g, s, s.u, w_empty);
    const int flo = std::max(0, j - 1 - cap);
    for (int f = flo; f < j - 1; ++f) {
      lchk(LT_ALS + o, f);
      const int unpaired = j - f - 1;
      double w_ext = E.w_ap + E.w_bp + E.w_cp * unpaired + pair_bonus;
      if (dz) w_ext += E.run(f + 1, j - 1);
      for (Grp& g : F.Glss[f].grps)
        for (Ent& s : g.ents) close(g, s, s.u, w_ext);
    }
  }
  // New outer pairs i.j of family F (fill_Qg (b), (f), (g), (h); (a) = (b) on the seeds Gs).
  // B: every new left end >= the reserved xs (group key).
  void l2r_grow(LFam& F, int j) {
    const int o = F.isU ? 4 : 0;
    Tab4& G = F.G[j];
    const int tg = LT_AG + o;
    auto alive = [&](const Grp& g, const Ent& s) { return F.isU || l2r_arm(s.u, g.v, j); };
    const int base = F.isU ? 4 : 0;
    ext_slot = base;
    for (int f = std::max(0, j - 1 - cap); f < j; ++f) {
      lchk(LT_AG + o, f);
      lchk(LT_AGS + o, f);
      for (int src = 0; src < 2; ++src)
        for (Grp& g : (src ? F.Gs[f] : F.G[f]).grps) {
          const int lo = F.isU ? g.v : 0;
          for (Ent& s : g.ents) {
            if (!alive(g, s)) continue;
            const int c = s.p;
            const int iu = F.isU ? 1 : 0;
            const int i0 = std::max(lo, c - 1 - cap);
            auto take = [&](int i) {
              const double w = E.stem_interior(i, c, f, j);
              if (w < HALF_NEG) return;
              ledge(G, tg, g.v, i, s.u, s, w + pair_bonus);
            };
            if (next_pair) {
              const std::vector<int>& pl = partners_left[static_cast<size_t>(j)];
              for (auto it = std::lower_bound(pl.begin(), pl.end(), i0); it != pl.end() && *it < c;
                   ++it) {
                if (!l2r_adj) ++it_int[iu];
                take(*it);
              }
            } else {
              for (int i = i0; i < c; ++i) {
                if (!l2r_adj) ++it_int[iu];
                if (!E.pairs(i, j)) continue;
                take(i);
              }
            }
        }
      }
    }
    ext_slot = base + 1;
    // --grammar simple: the multiloop is M1..M8, not cases (f)(g)(h) and not the tail.
    if (grammar_simple) {
      ext_slot = -1;
      return;
    }
    if (rlt) {
      l2r_close_tails(F, j, tg);
    } else {
      for (int f = 0; f < j; ++f) {
        lchk(LT_ALS + o, f);
        const int unpaired = j - f - 1;
        double w_ext = E.w_ap + E.w_bp + E.w_cp * unpaired + pair_bonus;
        if (dz) w_ext += E.run(f + 1, j - 1);
        for (int src = 0; src < 2; ++src) {
          if (src == 1 && unpaired > cap) continue;
          for (Grp& g : (src ? F.Glss[f] : F.Gls[f]).grps) {
            const int lo = F.isU ? g.v : 0;
            for (Ent& s : g.ents) {
              const int i = s.p - 1;
              if (i < lo || !alive(g, s) || !E.pairs(i, j) || !E.end_ok(i, j)) continue;
              ledge(G, tg, g.v, i, s.u, s, w_ext + E.at(i, j));
            }
          }
        }
      }
    }
    if (j < 1) {
      ext_slot = -1;
      return;
    }
    lchk(LT_ARS + o, j - 1);
    for (Grp& g : F.Grs[j - 1].grps) {
      const int lo = F.isU ? g.v : 0;
      for (Ent& s : g.ents) {
        if (!alive(g, s)) continue;
        const int c = s.p;
        ext_slot = base + 2;
        for (int i = std::max(lo, lead_i_lo(c)); i < c; ++i) {
          if (!E.pairs(i, j) || !E.end_ok(i, j)) continue;
          double w = E.w_ap + E.w_bp + E.w_cp * (c - i - 1) + E.at(i, j) + pair_bonus;
          if (dz) w += E.run(i + 1, c - 1);
          ledge(G, tg, g.v, i, s.u, s, w);
        }
        if (c - 1 < 0) continue;
        lchk(LT_QM, c - 1);
        ext_slot = base + 3;
        const OMap& at = qmp_at(c - 1);
        const int iu = F.isU ? 1 : 0;
        qm_pair_slots(c - 1, j, [&](int r) {
          if (!l2r_adj) ++it_h[iu];
          const int ii = at.keys[r];
          if (ii < 0) return;
          const int outer = ii - 1;
          if (outer < lo || outer >= c || !E.pairs(outer, j) || !E.end_ok(outer, j)) return;
          ledge(G, tg, g.v, outer, s.u, s, E.w_ap + E.w_bp + E.at(outer, j) + pair_bonus, ii,
                c - 1, at.vals[r]);
        });
      }
    }
    ext_slot = -1;
  }

  // Gls / Glss at j: Qm(i, c-1) left of a branch c.j (fill_Qgls; seeds from Gs).
  void l2r_ls(LFam& F, int j) {
    const int o = F.isU ? 4 : 0;
    for (int src = 0; src < 2; ++src) {
      Tab4& out = src ? F.Glss[j] : F.Gls[j];
      for (Grp& g : (src ? F.Gs[j] : F.G[j]).grps) {
        const int lo = F.isU ? g.v + 1 : 1;
        for (Ent& s : g.ents) {
          const int c = s.p;
          if (c - 1 < 0 || !E.end_ok(c, j)) continue;
          lchk(LT_QM, c - 1);
          const double w = E.w_bp + E.at(c, j);
          const OMap& at = qmp_at(c - 1);
          const int iu = F.isU ? 1 : 0;
          for (size_t r = 0; r < at.slots(); ++r) {
            if (!l2r_adj) ++it_ls[iu];
            const int i = at.keys[r];
            if (i < lo) continue;
            ledge(out, LT_ALS + o, g.v, i, s.u, s, w, i, c - 1, at.vals[r]);
          }
        }
      }
    }
  }

  // Grs at the end of step j: a branch c.f, f < j, then Qm(f+1, j) (fill_Qgrs).
  void l2r_rs(LFam& F, int j) {
    const int o = F.isU ? 4 : 0;
    lchk(LT_QM, j);
    Tab4& out = F.Grs[j];
    const OMap& at = qmp_at(j);
    const int iu = F.isU ? 1 : 0;
    for (size_t r = 0; r < at.slots(); ++r) {
      if (!l2r_adj) ++it_rs[iu];
      const int p = at.keys[r];
      if (p < 1) continue;
      const int f = p - 1;
      lchk(LT_AG + o, f);
      lchk(LT_AGS + o, f);
      for (int src = 0; src < 2; ++src)
        for (Grp& g : (src ? F.Gs[f] : F.G[f]).grps)
          for (Ent& s : g.ents) {
            if (!E.end_ok(s.p, f) || (!F.isU && !l2r_arm(s.u, g.v, j))) continue;
            ledge(out, LT_ARS + o, g.v, s.p, s.u, s, E.w_bp + E.at(s.p, f), p, j, at.vals[r]);
          }
    }
  }

  // Outside of the A seeds: Rs(a2, xs, x', j) -> Qz(a2+1, xs-1), Qz(x'+1, j-1).
  void l2r_spawn_back(int j) {
    for (Grp& g : LA.Gs[j].grps) {
      const int xp = g.v;
      for (Ent& s : g.ents) {
        if (s.b < HALF_NEG) continue;
        const int a2 = s.p, xs = s.u;
        const double k0 = pair_bonus + E.w_stem + E.at(a2, j);
        const double f1 = qz_log(a2 + 1, xs - 1), f2 = qz_log(xp + 1, j - 1);
        if (xs - 1 >= a2 + 1) add2b(Qz_b, a2 + 1, xs - 1, s.b + k0 + f2);
        if (j - 1 >= xp + 1) add2b(Qz_b, xp + 1, j - 1, s.b + k0 + f1);
      }
    }
  }

  // --l2r-spawn cube (Huang & Chiang 2007 cube pruning) at the inner pair a2.c1 (c1 = j).
  // Rows: xs > a2 by Qzbar(a2+1, xs-1) + Wmax(xs) (pre-pass Qzbar: the scan's Qz right of j is
  // not filled yet), built on a2's first spawn and cached.  Columns: xs's V arms by the static
  // weight log sum_{y'} E_in V (lYwS at the arm's first y').  Priority of cell (i, c): exact
  // F1(a2, xs_i) + F2(x'_ic, c1) + W_c1(arm), W_c1 = lS(arm, c1) = log sum_{y' > c1} E_in V
  // (--l2r-spawn-w suf, default) or lK(arm, c1), the same sum with the gap Qz(c1+1, y'-1) (k):
  // R's K at q' = c1.  Pop the best; push (i, c+1),
  // and the next row's (i', 0) when c == 0; invalid cells (x' >= c1, no y' > c1, a B that
  // cannot close right of c1, F = 0) are counted as pops and expanded at once.  Stop after
  // sp_pool valid pops or sp_maxpop pops; the valid ones re-ranked by the Rs score (alpha +
  // lsc_L + context, as l2r_pruneA) and the top l2r_spawn_k spawned with value k0 F1 F2.
  // --l2r-spawn cube (default form): per step c1 = j every arm valid at c1 gets G = F2(x', c1) +
  // W_c1 (a2-free), columns are xs's arms sorted by G, and per pair every row head (exact
  // F1(a2, xs) + best G of xs) enters one heap: the axes are exactly monotone, so the pops are
  // the exact top sp_pool by F1 + F2 + W_c1.  cube-static: the cached rows / static columns
  // above (sp_static).
  bool sp_cube = false, sp_wk = false, sp_static = false;
  int sp_pool = 128, sp_maxpop = 0;
  std::vector<std::vector<std::pair<double, int>>> spG;
  std::vector<int> spGx;
  int spG_j = -1;
  size_t sp_col_ents = 0, sp_col_max = 0;
  void l2r_spawn_step_cols(int j) {
    SecTimer t(sp_col_sec);
    if (spG.empty()) spG.assign(n, {});
    for (int xs : spGx) spG[xs].clear();
    spGx.clear();
    size_t tot = 0;
    for (int r = 0; r < prXsFirst[j]; ++r) {
      const int xs = prArmXs[r], xp = prArmXp[r];
      if (xp >= j || lmaxR[xs] <= j || lmaxR[xp] <= j) continue;
      const double f2 = qz_log(xp + 1, j - 1);
      if (f2 < HALF_NEG) continue;
      const double w = sp_wk ? lK(r, j) : lS(r, j);
      if (w < HALF_NEG) continue;
      if (spG[xs].empty()) spGx.push_back(xs);
      spG[xs].emplace_back(f2 + w, r);
      ++tot;
    }
    for (int xs : spGx)
      std::sort(spG[xs].begin(), spG[xs].end(),
                [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                  return a.first > b.first || (a.first == b.first && a.second < b.second);
                });
    std::sort(spGx.begin(), spGx.end());
    sp_col_ents += tot;
    sp_col_max = std::max(sp_col_max, tot);
    spG_j = j;
  }
  void l2r_spawn_cube_step(int a2, int j, double k0) {
    if (spG_j != j) l2r_spawn_step_cols(j);
    struct Cell {
      double pr;
      int xs, c;
      bool operator<(const Cell& o) const { return pr < o.pr; }
    };
    std::vector<Cell> heap;
    std::vector<std::pair<double, int>> pool;  // (F1 + F2, arm)
    long long pops = 0;
    {
      SecTimer t(sp_heap_sec);
      for (auto it = std::upper_bound(spGx.begin(), spGx.end(), a2); it != spGx.end(); ++it) {
        const int xs = *it;
        const double f1 = qz_log(a2 + 1, xs - 1);
        if (f1 > HALF_NEG) heap.push_back({f1 + spG[xs][0].first, xs, 0});
      }
      std::make_heap(heap.begin(), heap.end());
      const int maxpop = sp_maxpop > 0 ? sp_maxpop : 4 * sp_pool;
      while (!heap.empty() && static_cast<int>(pool.size()) < sp_pool && pops < maxpop) {
        std::pop_heap(heap.begin(), heap.end());
        const Cell x = heap.back();
        heap.pop_back();
        ++pops;
        const auto& col = spG[x.xs];
        const double f1 = x.pr - col[x.c].first;
        const int r = col[x.c].second;
        pool.emplace_back(f1 + qz_log(prArmXp[r] + 1, j - 1), r);
        if (x.c + 1 < static_cast<int>(col.size())) {
          heap.push_back({f1 + col[x.c + 1].first, x.xs, x.c + 1});
          std::push_heap(heap.begin(), heap.end());
        }
      }
    }
    sp_pops += pops;
    sp_cand += static_cast<long long>(pool.size());
    if (ltr_f && a2 == ltr_a2 && j == ltr_c1) ltr_pool(a2, j, pool, pops, heap.size());
    l2r_spawn_rank(a2, j, k0, pool);
  }
  int ltr_a2 = -1, ltr_c1 = -1;
  void ltr_pool(int a2, int j, const std::vector<std::pair<double, int>>& pool, long long pops, size_t left) {
    const int rt = larm(ltr_xs, ltr_xp);
    const double f1 = qz_log(a2 + 1, ltr_xs - 1), f2 = qz_log(ltr_xp + 1, j - 1);
    const double w = rt >= 0 ? (sp_wk ? lK(rt, j) : lS(rt, j)) : NEG_INF;
    int col = -1;
    for (size_t c = 0; c < spG[ltr_xs].size(); ++c)
      if (spG[ltr_xs][c].second == rt) col = static_cast<int>(c);
    const double pt = f1 + f2 + w;
    int in = -1, rank = 1;
    double lo = INFINITY;
    for (size_t q = 0; q < pool.size(); ++q) {
      const int r = pool[q].second;
      const double p = pool[q].first + (sp_wk ? lK(r, j) : lS(r, j));
      lo = std::min(lo, p);
      rank += p > pt;
      if (r == rt) in = static_cast<int>(q);
    }
    std::fprintf(ltr_f, "spawn-pool a2 %d c1 %d arm %d f1 %.6f f2 %.6f w %.6f pr %.6f col %d/%zu in_pool %d "
                 "rank %d pool %zu pool_min %.6f pops %lld heap_left %zu sp_pool %d\n",
                 a2, j, rt, f1, f2, w, pt, col, spG[ltr_xs].size(), in, rank, pool.size(), lo, pops, left,
                 sp_pool);
  }
  void ltr_rank(int a2, int j, double k0, const RN& rn, const std::vector<std::pair<double, int>>& pool,
                const std::vector<std::pair<double, int>>& rk) {
    const int rt = larm(ltr_xs, ltr_xp);
    const double f1 = qz_log(a2 + 1, ltr_xs - 1), f2 = qz_log(ltr_xp + 1, j - 1);
    const double a = k0 + f1 + f2, s = lsc_L(rn, rt, 0), kb = lkb_arm(a2, rt);
    const double tot = a + s + kb;
    int rank = 1;
    std::vector<double> v;
    for (const auto& e : rk) {
      rank += e.first + k0 > tot;
      v.push_back(e.first + k0);
    }
    const size_t k = std::min(rk.size(), static_cast<size_t>(l2r_spawn_k));
    double cut = NEG_INF;
    if (k > 0) {
      std::nth_element(v.begin(), v.begin() + (k - 1), v.end(), std::greater<double>());
      cut = v[k - 1];
    }
    double full = NEG_INF;
    for (int m = 0; m <= g_ml_L && m < rnL; ++m) {
      const RV& x = rn.L[m];
      for (size_t i = 0; i < x.v.size(); ++i)
        if (x.v[i] > 0.0) {
          const double kk = lK(rt, x.lo + static_cast<int>(i));
          if (kk > HALF_NEG) full = log_add(full, x.mx + std::log(x.v[i]) + kk);
        }
    }
    std::fprintf(ltr_f, "spawn-rank a2 %d c1 %d arm %d a %.6f lsc_L %.6f ctx %.6f score %.6f rank %d of %zu k %zu "
                 "cut %.6f lsc_L_untrunc %.6f pool %zu\n",
                 a2, j, rt, a, s, kb, tot, rank, rk.size(), k, cut, full, pool.size());
    for (const auto& [qp, g] : rn.T[0])
      std::fprintf(ltr_f, "spawn-T0 qp %d lam %.6f K %.6f Kx %.6f prefix_a? zbar0q %.6f\n", qp, g, lK(rt, qp),
                   lKx(rt, qp), stem_->Qzbar(0, qp));
    if (rt >= 0)
      for (int k = prArmY0[rt]; k < prArmY0[rt + 1]; ++k)
        std::fprintf(ltr_f, "spawn-V y %d Yw %.6f zbar0y %.6f\n", prYy[k], prYw[k], stem_->Qzbar(0, prYy[k] - 1));
    for (int m = 0; m <= g_ml_L && m < rnL; ++m) {
      const RV& x = rn.L[m];
      for (size_t i = 0; i < x.v.size(); ++i)
        if (x.v[i] > 0.0)
          std::fprintf(ltr_f, "spawn-L m %d qp %d lam %.6f K %.6f\n", m, x.lo + static_cast<int>(i),
                       x.mx + std::log(x.v[i]), lK(rt, x.lo + static_cast<int>(i)));
    }
  }
  void l2r_spawn_rank(int a2, int j, double k0, const std::vector<std::pair<double, int>>& pool) {
    if (pool.empty()) return;
    ProfReg reg(*this, 2);
    std::vector<std::pair<double, int>> rk;
    {
      SecTimer t(sp_rank_sec);
      const RN& rn = rec_root(a2, j, true);
      for (size_t q = 0; q < pool.size(); ++q) {
        const int r = pool[q].second;
        if (lsub_e) {
          const double s = lors(prArmXs[r], prArmXp[r], a2, j);
          if (s > HALF_NEG) rk.emplace_back(pool[q].first + s, static_cast<int>(q));
          continue;
        }
        sx_ksrc = 3;
        const double s = lsc_L(rn, r, 0);
        sx_ksrc = 0;
        if (s > HALF_NEG) rk.emplace_back(pool[q].first + s + lkb_arm(a2, r), static_cast<int>(q));
      }
      if (ltr_f && a2 == ltr_a2 && j == ltr_c1) ltr_rank(a2, j, k0, rn, pool, rk);
    }
    const size_t k = std::min(rk.size(), static_cast<size_t>(l2r_spawn_k));
    std::partial_sort(rk.begin(), rk.begin() + static_cast<long>(k), rk.end(),
                      [](const std::pair<double, int>& x, const std::pair<double, int>& y) {
                        return x.first > y.first || (x.first == y.first && x.second < y.second);
                      });
    sp_kept += static_cast<long long>(k);
    for (size_t c = 0; c < k; ++c) {
      const auto& p = pool[rk[c].second];
      lacc(LA.Gs[j], LT_AGS, prArmXp[p.second], a2, prArmXs[p.second], k0 + p.first);
    }
  }
  long long sp_pops = 0;
  double sp_row_sec = 0.0, sp_col_sec = 0.0, sp_heap_sec = 0.0, sp_rank_sec = 0.0;
  std::vector<std::vector<int>> spCol, spRow;
  std::vector<char> spRowDone;
  std::vector<double> spWmax;
  size_t sp_row_ents = 0;
  void l2r_spawn_cols() {
    SecTimer t(sp_col_sec);
    spCol.assign(n, {});
    spWmax.assign(n, NEG_INF);
    spRow.assign(n, {});
    spRowDone.assign(n, 0);
    std::vector<double> w(prA, NEG_INF);
    for (int r = 0; r < prA; ++r) {
      if (prArmY0[r] < prArmY0[r + 1]) w[r] = lYwS[prArmY0[r]];
      if (w[r] < HALF_NEG) continue;
      spCol[prArmXs[r]].push_back(r);
      spWmax[prArmXs[r]] = std::max(spWmax[prArmXs[r]], w[r]);
    }
    for (auto& c : spCol)
      std::stable_sort(c.begin(), c.end(), [&](int a, int b) { return w[a] > w[b]; });
  }
  const std::vector<int>& l2r_spawn_row(int a2) {
    if (spRowDone[a2]) return spRow[a2];
    SecTimer t(sp_row_sec);
    spRowDone[a2] = 1;
    std::vector<std::pair<double, int>> v;
    for (int xs = a2 + 1; xs < n; ++xs) {
      if (spCol[xs].empty()) continue;
      const double z = stem_->Qzbar(a2 + 1, xs - 1);
      if (z > HALF_NEG) v.emplace_back(z + spWmax[xs], xs);
    }
    std::stable_sort(v.begin(), v.end(), [](const std::pair<double, int>& x,
                                           const std::pair<double, int>& y) { return x.first > y.first; });
    std::vector<int>& row = spRow[a2];
    row.reserve(v.size());
    for (const auto& e : v) row.push_back(e.second);
    sp_row_ents += row.size();
    return row;
  }
  void l2r_spawn_cube(int a2, int j, double k0) {
    const std::vector<int>& row = l2r_spawn_row(a2);
    const int R = static_cast<int>(row.size());
    const int maxpop = sp_maxpop > 0 ? sp_maxpop : 4 * sp_pool;
    std::vector<double> f1(R, 1.0);  // 1.0: not yet computed
    auto next_row = [&](int i) {
      for (; i < R; ++i) {
        const int xs = row[i];
        if (xs < j && lmaxR[xs] > j) return i;
      }
      return -1;
    };
    struct Cell {
      double pr;
      int i, c;
      bool operator<(const Cell& o) const { return pr < o.pr; }
    };
    std::priority_queue<Cell> heap;
    std::vector<std::pair<int, int>> todo;
    std::vector<std::pair<double, int>> pool;  // (F1 + F2, arm)
    long long pops = 0;
    auto offer = [&](int i0, int c0) {
      todo.emplace_back(i0, c0);
      while (!todo.empty()) {
        const auto [i, c] = todo.back();
        todo.pop_back();
        if (i < 0 || c >= static_cast<int>(spCol[row[i]].size())) continue;
        const int xs = row[i], r = spCol[xs][c], xp = prArmXp[r];
        if (f1[i] > 0.5) f1[i] = qz_log(a2 + 1, xs - 1);
        double pr = NEG_INF;
        if (xp < j && lmaxR[xp] > j && f1[i] > HALF_NEG) {
          const double f2 = qz_log(xp + 1, j - 1), w = sp_wk ? lK(r, j) : lS(r, j);
          if (f2 > HALF_NEG && w > HALF_NEG) pr = f1[i] + f2 + w;
        }
        if (pr > HALF_NEG) {
          heap.push({pr, i, c});
          continue;
        }
        if (++pops >= maxpop) return;
        todo.emplace_back(i, c + 1);
        if (c == 0) todo.emplace_back(next_row(i + 1), 0);
      }
    };
    {
      SecTimer t(sp_heap_sec);
      offer(next_row(0), 0);
      while (!heap.empty() && static_cast<int>(pool.size()) < sp_pool && pops < maxpop) {
        const Cell x = heap.top();
        heap.pop();
        ++pops;
        const int r = spCol[row[x.i]][x.c];
        pool.emplace_back(f1[x.i] + qz_log(prArmXp[r] + 1, j - 1), r);
        offer(x.i, x.c + 1);
        if (x.c == 0) offer(next_row(x.i + 1), 0);
      }
    }
    sp_pops += pops;
    sp_cand += static_cast<long long>(pool.size());
    l2r_spawn_rank(a2, j, k0, pool);
  }

  // A seeds at the inner pair a2.j: every arm a2 < xs <= x' < j whose B can close right of j.
  void l2r_spawn(int j) {
    std::vector<std::pair<double, int>> cand;
    for (int a2 : partners_left[j]) {
      if (!E.end_ok(a2, j)) continue;
      const double k0 = pair_bonus + E.w_stem + E.at(a2, j);
      ++sp_events;
      if (sp_cube) {
        lchk(LT_QZ, j - 1);
        if (sp_static) l2r_spawn_cube(a2, j, k0);
        else l2r_spawn_cube_step(a2, j, k0);
        continue;
      }
      if (l2r_spawn_k > 0) {
        cand.clear();
        lchk(LT_QZ, j - 1);
        for (int r = prXsFirst[a2 + 1]; r < prXsFirst[j]; ++r) {
          const int xs = prArmXs[r], xp = prArmXp[r];
          if (xp >= j || lmaxR[xs] <= j || lmaxR[xp] <= j) continue;
          const double f1 = qz_log(a2 + 1, xs - 1), f2 = qz_log(xp + 1, j - 1);
          if (f1 < HALF_NEG || f2 < HALF_NEG) continue;
          const double s = lS(r, j);
          if (s < HALF_NEG) continue;
          cand.emplace_back(f1 + s + f2, r);
        }
        sp_cand += static_cast<long long>(cand.size());
        const size_t k = std::min(cand.size(), static_cast<size_t>(l2r_spawn_k));
        std::partial_sort(cand.begin(), cand.begin() + static_cast<long>(k), cand.end(),
                          [](const std::pair<double, int>& x, const std::pair<double, int>& y) {
                            return x.first > y.first || (x.first == y.first && x.second < y.second);
                          });
        sp_kept += static_cast<long long>(k);
        for (size_t c = 0; c < k; ++c) {
          const int r = cand[c].second, xs = prArmXs[r], xp = prArmXp[r];
          lacc(LA.Gs[j], LT_AGS, xp, a2, xs, k0 + qz_log(a2 + 1, xs - 1) + qz_log(xp + 1, j - 1));
        }
        continue;
      }
      const long long before = l2r_tr[LT_AGS];
      for (int xs = a2 + 1; xs < j; ++xs) {
        if (lmaxR[xs] <= j) continue;
        lchk(LT_QZ, xs - 1);
        const double f1 = qz_log(a2 + 1, xs - 1);
        if (f1 < HALF_NEG) continue;
        for (int xp = xs; xp < j; ++xp) {
          if (lmaxR[xp] <= j) continue;
          lchk(LT_QZ, j - 1);
          const double f2 = qz_log(xp + 1, j - 1);
          if (f2 < HALF_NEG) continue;
          lacc(LA.Gs[j], LT_AGS, xp, a2, xs, k0 + f1 + f2);
        }
      }
      sp_cand += l2r_tr[LT_AGS] - before;
      sp_kept += l2r_tr[LT_AGS] - before;
    }
  }

  void l2r_closeA(int j) {
    SMap& S = LS[j];
    for (int src = 0; src < 2; ++src) {
      if (src == 1 && !single_stems) continue;
      for (const Grp& g : (src ? LA.Gs[j] : LA.G[j]).grps)
        for (const Ent& s : g.ents) {
          if (!E.end_ok(s.p, j)) continue;
          if (grammar_simple && !single_stems && gs_seedA.count(key3(s.p, s.u, g.v))) continue;
          SSt& st = S.v[S.at(key3(s.p, s.u, g.v))];
          st.c = log_add(st.c, s.a + E.at(s.p, j));
          ++l2r_tr[LT_S];
          if (grammar_simple) ++gs_phase;
        }
    }
  }

  // B seeds x'.j from S at j-1.
  void l2r_bseed(int j) {
    if (j < 1) return;
    lchk(LT_S, j - 1);
    const SMap& pv = LS[j - 1];
    for (size_t k = 0; k < pv.keys.size(); ++k) {
      const double f = pv.v[k].F;
      if (f < HALF_NEG) continue;
      const uint64_t key = pv.keys[k];
      const int xp = k3c(key);
      if (!E.pairs(xp, j) || !E.end_ok(xp, j)) continue;
      if (grammar_simple) {
        lacc(LB.Gs[j], LT_BGS, k3b(key), xp, k3a(key), f + E.w_stem + E.at(xp, j) + pair_bonus);
        gs_seedB.insert(key3(xp, k3a(key), k3b(key)));
        ++gs_phase;
      } else {
        lacc(LB.Gs[j], LT_BGS, k3b(key), xp, k3a(key), f + E.w_stem + E.at(xp, j) + pair_bonus);
      }
    }
  }

  void lc_cmp(double x, double y, long long& cnt, long long& mis, double& rel) {
    ++cnt;
    const bool fx = x > HALF_NEG, fy = y > HALF_NEG;
    if (fx != fy) {
      ++mis;
      rel = INFINITY;
      return;
    }
    if (!fx) return;
    const double r = std::expm1(std::fabs(x - y));
    rel = std::max(rel, r);
    if (r > 1e-10) ++mis;
  }

  void l2r_closeB(int j) {
    for (int src = 0; src < 2; ++src) {
      if (src == 1 && !single_stems) continue;
      for (const Grp& g : (src ? LB.Gs[j] : LB.G[j]).grps) {
        const int xs = g.v;
        if (!E.end_ok(xs, j)) continue;
        for (const Ent& s : g.ents) {
          if (s.p != xs) continue;
          if (grammar_simple && !single_stems && gs_seedB.count(key3(s.p, s.u, xs))) continue;
          const double val = s.a + E.at(xs, j);
          if (grammar_simple) ++gs_phase;
          if (l2r_mode == 1) set_pair(Qp, Qp_at, s.u, j, val);
          else LQp[ix(s.u, j)] = log_add(LQp[ix(s.u, j)], val);
        }
      }
    }
    if (l2r_mode == 2)
      for (int i = 0; i <= j; ++i) {
        const double x = Qp[ix(i, j)], y = LQp[ix(i, j)];
        if (x > HALF_NEG || y > HALF_NEG) lc_cmp(x, y, lc_Qp_n, lc_Qp_mis, lc_Qp_rel);
      }
  }

  // Check 1.1: Afull(a, xs, x', j) against sum_{a2, c1} Qg(a, a2, c1, j) w_stem at(a2, c1) F1 F2
  // from the engine's Qg[j] (same live-arm filter).
  void l2r_checkR(int j) {
    std::unordered_map<uint64_t, double> def, rec;
    for (const Grp& g : Qg[j].grps) {
      const int c1 = g.v;
      for (const Ent& s : g.ents) {
        const int a = s.p, a2 = s.u;
        if (!E.end_ok(a2, c1)) continue;
        const double base = s.a + E.w_stem + E.at(a2, c1);
        for (int xs = a2 + 1; xs < c1; ++xs) {
          if (lmaxR[xs] <= j) continue;
          const double f1 = qz_log(a2 + 1, xs - 1);
          if (f1 < HALF_NEG) continue;
          for (int xp = xs; xp < c1; ++xp) {
            if (lmaxR[xp] <= j) continue;
            const double f2 = qz_log(xp + 1, c1 - 1);
            if (f2 < HALF_NEG) continue;
            double& d = def.try_emplace(key3(a, xs, xp), NEG_INF).first->second;
            d = log_add(d, base + f1 + f2);
          }
        }
      }
    }
    for (int src = 0; src < 2; ++src) {
      if (src == 1 && !single_stems) continue;
      for (const Grp& g : (src ? LA.Gs[j] : LA.G[j]).grps)
        for (const Ent& s : g.ents) {
          double& d = rec.try_emplace(key3(s.p, s.u, g.v), NEG_INF).first->second;
          d = log_add(d, s.a);
        }
    }
    for (const auto& [k, x] : def) {
      auto it = rec.find(k);
      lc_cmp(x, it == rec.end() ? NEG_INF : it->second, lc_R_n, lc_R_mis, lc_R_rel);
    }
    for (const auto& [k, y] : rec)
      if (!def.count(k)) lc_cmp(NEG_INF, y, lc_R_n, lc_R_mis, lc_R_rel);
  }

  // S at the end of step j (after fill_Qz(j)): the Qz recursion of fill_Qz summed over q'.
  void l2r_sext(int j) {
    l2r_fext(LS, LT_S, j, [&](uint64_t key) { return l2r_arm(k3b(key), k3c(key), j); });
  }
  // One fill_Qz step of the filler table T at j (S, P1, P3): every state of T[j-1] / T[d-1]
  // with alive(key) at j takes one unpaired base / the item d.j; F = closed run + c (the
  // filler that starts after j, entered by the caller before this call).
  template <class Alive>
  void l2r_fext(std::vector<SMap>& T, int lt, int j, Alive alive) {
    const int fi = lt == LT_S ? 0 : lt == LT_P1 ? 1 : 2;
    std::unique_ptr<SecTimer> ptm(lpf ? new SecTimer(lpf_fx_sec[fi]) : nullptr);
    long long fpair = 0;
    if (lpf) ++lpf_fx_pos[fi];
    lchk(LT_QB, j);
    if (allow_pk) lchk(LT_QP, j);
    SMap& cur = T[j];
    auto branch = [&](int d, double bw) {
      if (d < 1) return;
      lchk(lt, d - 1);
      const SMap& pv = T[d - 1];
      if (lpf) {
        ++lpf_fx_items[fi];
        fpair += static_cast<long long>(pv.keys.size());
      }
      for (size_t k = 0; k < pv.keys.size(); ++k) {
        const double f = pv.v[k].F;
        const uint64_t key = pv.keys[k];
        if (f < HALF_NEG || !alive(key)) continue;
        SSt& s = cur.v[cur.at(key)];
        s.P = log_add(s.P, f + bw);
        ++l2r_tr[lt];
        if (grammar_simple) ++gs_fill;
      }
    };
    const OMap& atb = Qb_at[j];
    for (size_t r = 0; r < atb.slots(); ++r) {
      const int d = atb.keys[r];
      if (d < 0 || !E.end_ok(d, j)) continue;
      branch(d, atb.vals[r] + (E.w_nest + E.at(d, j)));
    }
    if (allow_pk) {
      const OMap& atp = Qp_at[j];
      for (size_t r = 0; r < atp.slots(); ++r) {
        const int d = atp.keys[r];
        if (d < 0) continue;
        branch(d, atp.vals[r] + E.w_beta1p + E.nest_pk_k * E.w_nest);
      }
    }
    if (lpf) {
      lpf_fx_pair[fi] += fpair;
      lpf_fx_pairmax[fi] = std::max(lpf_fx_pairmax[fi], fpair);
    }
    if (j >= 1) {
      lchk(lt, j - 1);
      const SMap& pv = T[j - 1];
      if (lpf) lpf_fx_unp[fi] += static_cast<long long>(pv.keys.size());
      for (size_t k = 0; k < pv.keys.size(); ++k) {
        const uint64_t key = pv.keys[k];
        if (!alive(key)) continue;
        const SSt& o = pv.v[k];
        if (dz) {
          const double u1 = log_add(o.P, o.c), u2 = log_add(o.U2, o.U1 + E.d5w[j - 1]);
          if (u1 < HALF_NEG && u2 < HALF_NEG) continue;
          SSt& s = cur.v[cur.at(key)];
          s.U1 = u1 > HALF_NEG ? u1 + E.wu : NEG_INF;
          s.U2 = u2 > HALF_NEG ? u2 + E.wu : NEG_INF;
        } else {
          if (o.F < HALF_NEG) continue;
          SSt& s = cur.v[cur.at(key)];
          s.P = log_add(s.P, o.F + E.wu);
        }
        ++l2r_tr[lt];
        if (grammar_simple) ++gs_fill;
      }
    }
    for (SSt& s : cur.v)
      s.F = log_add(dz ? run_close(s.P, s.U1, s.U2, E.s1w[j], E.d3w[j]) : s.P, s.c);
  }

  // ---- --l2r-spawn inc: P1 / P2 / P3 (forward).  alive: A (a2) and B (xs, x') can still
  // pair right of every c1 the state can reach (c1 >= j+1; the arm needs lmaxR > c1, as spawn).
  bool lpa1(int a2, int j) const { return lmaxR[a2] > j + 1; }
  bool lpa2(int a2, int xs, int j) const { return lmaxR[a2] > j && lmaxR[xs] > j + 1; }
  bool lpa3(int a2, int xs, int xp, int j) const {
    return lmaxR[a2] > j && lmaxR[xs] > j + 1 && lmaxR[xp] > j + 1;
  }
  bool l2r_pairable(int a2, int c1) const {
    const std::vector<int>& pl = partners_left[c1];
    return std::binary_search(pl.begin(), pl.end(), a2) && E.end_ok(a2, c1);
  }
  long long lp_seed_checks = 0, lp_seeds = 0;
  void l2r_p1_step(int j) {
    if (lpa1(j, j)) {
      SSt& s = LP1[j].v[LP1[j].at(key3(j, 0, 0))];
      s.c = 0.0;
      ++l2r_tr[LT_P1];
    }
    l2r_fext(LP1, LT_P1, j, [&](uint64_t key) { return lpa1(k3a(key), j); });
  }
  void l2r_p2_step(int j) {
    if (j < 1) return;
    SMap& cur = LP2[j];
    lchk(LT_P2, j - 1);
    const SMap& pv = LP2[j - 1];
    for (size_t k = 0; k < pv.keys.size(); ++k) {
      const uint64_t key = pv.keys[k];
      const double f = pv.v[k].F;
      if (f < HALF_NEG || !lpa2(k3a(key), k3b(key), j)) continue;
      cur.v[cur.at(key)].F = f;
      ++l2r_tr[LT_P2];
      if (grammar_simple) ++gs_fill;
    }
    lchk(LT_P1, j - 1);
    const SMap& p1 = LP1[j - 1];
    for (size_t k = 0; k < p1.keys.size(); ++k) {
      const int a2 = k3a(p1.keys[k]);
      const double f = p1.v[k].F;
      if (f < HALF_NEG || !lpa2(a2, j, j)) continue;
      cur.v[cur.at(key3(a2, j, 0))].F = f;
      ++l2r_tr[LT_P2];
      if (grammar_simple) ++gs_phase;
    }
  }
  void l2r_p3_step(int j) {
    SMap& cur = LP3[j];
    lchk(LT_P2, j);
    const SMap& p2 = LP2[j];
    for (size_t k = 0; k < p2.keys.size(); ++k) {
      const int a2 = k3a(p2.keys[k]), xs = k3b(p2.keys[k]);
      const double f = p2.v[k].F;
      if (f < HALF_NEG || !lpa3(a2, xs, j, j)) continue;
      SSt& s = cur.v[cur.at(key3(a2, xs, j))];
      s.c = log_add(s.c, f);
      ++l2r_tr[LT_P3];
      if (grammar_simple) ++gs_phase;
    }
    l2r_fext(LP3, LT_P3, j, [&](uint64_t key) { return lpa3(k3a(key), k3b(key), k3c(key), j); });
  }
  // A seeds at the inner pair a2.j from P3 at j-1 (spawn's arm rule: lmaxR[xs], lmaxR[x'] > j).
  void l2r_inc_seed(int j) {
    if (j < 1) return;
    lchk(LT_P3, j - 1);
    const SMap& pv = LP3[j - 1];
    for (size_t k = 0; k < pv.keys.size(); ++k) {
      const uint64_t key = pv.keys[k];
      const double f = pv.v[k].F;
      if (f < HALF_NEG) continue;
      const int a2 = k3a(key), xs = k3b(key), xp = k3c(key);
      if (lmaxR[xs] <= j || lmaxR[xp] <= j) continue;
      ++lp_seed_checks;
      if (!l2r_pairable(a2, j)) continue;
      ++lp_seeds;
      if (grammar_simple) {
        // Seeds are their own 4-index type (Rs), beamed at b4, not mixed into R.
        lacc(LA.Gs[j], LT_AGS, xp, a2, xs, (pair_bonus + E.w_stem + E.at(a2, j)) + f);
        gs_seedA.insert(key3(a2, xs, xp));
        ++gs_phase;
      } else {
        lacc(LA.Gs[j], LT_AGS, xp, a2, xs, (pair_bonus + E.w_stem + E.at(a2, j)) + f);
      }
    }
  }
  // --l2r-inc-check on (unbeamed): every A seed at j against spawn's k0 F1 F2 from the Qz table.
  bool l2r_inc_check = false;
  long long lc_I_n = 0, lc_I_mis = 0;
  double lc_I_rel = 0.0;
  void l2r_inc_check_step(int j) {
    std::unordered_map<uint64_t, double> def;
    for (int a2 : partners_left[j]) {
      if (!E.end_ok(a2, j)) continue;
      const double k0 = pair_bonus + E.w_stem + E.at(a2, j);
      for (int xs = a2 + 1; xs < j; ++xs) {
        if (lmaxR[xs] <= j) continue;
        const double f1 = qz_log(a2 + 1, xs - 1);
        if (f1 < HALF_NEG) continue;
        for (int xp = xs; xp < j; ++xp) {
          if (lmaxR[xp] <= j) continue;
          const double f2 = qz_log(xp + 1, j - 1);
          if (f2 < HALF_NEG) continue;
          def.emplace(key3(a2, xs, xp), k0 + f1 + f2);
        }
      }
    }
    std::unordered_map<uint64_t, double> rec;
    for (const Grp& g : LA.Gs[j].grps)
      for (const Ent& s : g.ents) rec.emplace(key3(s.p, s.u, g.v), s.a);
    for (const auto& [k, x] : def) {
      auto it = rec.find(k);
      lc_cmp(x, it == rec.end() ? NEG_INF : it->second, lc_I_n, lc_I_mis, lc_I_rel);
    }
    for (const auto& [k, y] : rec)
      if (!def.count(k)) lc_cmp(NEG_INF, y, lc_I_n, lc_I_mis, lc_I_rel);
  }

  // ---- beam (--pk-l2r on with a finite --l2r-beam; needs the principled pre-pass, --grow-rec)
  // Scores (log): R-type (R, Rs, Rls, Rlss, Rrs): alpha + log sum_q' g_L[q'] K_q'(arm), g_L the
  // (truncated) left growth vector of rec_root / rec_key, which already holds C[p'-1];
  // K_q'(arm) = sum_{y' > q'} Qz(q'+1, y'-1) E_in V (the engine's arm_phi, separable with
  // --pr-armq-sep).  U-type: alpha + g_R[xs] (V of the outer pair, with e^beta1 Sbar) +
  // C[a'-1].  S: F + K_j(arm) + C[a'-1].  Spawn: arms ranked by F1 + s_j(arm) + F2,
  // s_j = log sum_{y' > j} E_in V; top l2r_spawn_k (0: every arm, unranked).
  enum { PT_R, PT_RS, PT_RLS, PT_RLSS, PT_RRS, PT_U, PT_US, PT_ULS, PT_ULSS, PT_URS, PT_S,
         PT_P1, PT_P2, PT_P3, PT_N };
  int l2r_beam = 0, l2r_spawn_k = 0;
  // -1: l2r_beam; 0: none
  int l2r_tb[PT_N] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
  int ltb(int pt) const { return l2r_tb[pt] < 0 ? l2r_beam : l2r_tb[pt]; }
  std::unordered_map<uint64_t, int> lArm;
  std::vector<double> lYwS;
  long long l2r_pre[PT_N] = {}, l2r_kept[PT_N] = {}, l2r_kmax[PT_N] = {};
  long long l2r_pos[PT_N] = {}, l2r_bind[PT_N] = {};
  // scored states (pruning steps only), scores of -inf, R / S scores with no V arm
  long long l2r_nsc[PT_N] = {}, l2r_ninf[PT_N] = {}, l2r_noarm[PT_N] = {};
  // --l2r-ninf-mass on: per table, the states a beam cut scores -inf: log sum of their inside
  // values, kept / dropped counts, and the kept ones' keys per step (their posterior after the
  // outside pass).  --l2r-keep-ninf on: they are all kept, beyond the beam.
  bool l2r_ninf_mass = false, l2r_keep_ninf = false;
  std::vector<double> ninf_in = std::vector<double>(PT_N, NEG_INF);
  long long ninf_kept[PT_N] = {}, ninf_drop[PT_N] = {};
  std::vector<std::unordered_set<uint64_t>> ninf_keys[PT_N];
  static uint64_t ninf_key(int a, int b, int c) {
    return (static_cast<uint64_t>(a + 1) << 42) | (static_cast<uint64_t>(b + 1) << 21) |
           static_cast<uint64_t>(c + 1);
  }
  void ninf_note(int pt, int j, uint64_t key, double in, char& keep) {
    if (l2r_keep_ninf) keep = 1;
    if (!l2r_ninf_mass) return;
    if (in > HALF_NEG) ninf_in[pt] = log_add(ninf_in[pt], in);
    if (!keep) {
      ++ninf_drop[pt];
      return;
    }
    ++ninf_kept[pt];
    if (ninf_keys[pt].empty()) ninf_keys[pt].resize(n);
    ninf_keys[pt][j].insert(key);
  }
  // After the outside pass: per table "name ninf kept dropped log10(sum inside / Z) posterior"
  // of the -inf-scored states (posterior of the kept ones), then the posterior total.
  std::string l2r_ninf_report() {
    static const char* pn[PT_N] = {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us", "Uls", "Ulss",
                                   "Urs", "S", "P1", "P2", "P3"};
    outside();
    char buf[200];
    std::snprintf(buf, sizeof buf, " logZ %.9f | ninf kept dropped log10_in/Z post", logZ);
    std::string out = buf;
    double all = 0.0;
    for (int pt = 0; pt < PT_N; ++pt) {
      double post = 0.0;
      if (!ninf_keys[pt].empty())
        for (int j = 0; j < n; ++j) {
          const std::unordered_set<uint64_t>& ks = ninf_keys[pt][j];
          if (ks.empty()) continue;
          if (pt >= PT_S) {
            const SMap& m = l2r_smap(pt, j);
            for (size_t k = 0; k < m.keys.size(); ++k)
              if (ks.count(m.keys[k]) && m.v[k].F > HALF_NEG && m.v[k].bF > HALF_NEG)
                post += std::exp(m.v[k].F + m.v[k].bF - logZ);
          } else {
            for (const Grp& g : l2r_tab(pt, j).grps)
              for (const Ent& e : g.ents)
                if (ks.count(ninf_key(g.v, e.u, e.p)) && e.a > HALF_NEG && e.b > HALF_NEG)
                  post += std::exp(e.a + e.b - logZ);
          }
        }
      all += post;
      std::snprintf(buf, sizeof buf, " %s %lld %lld %lld %.4f %.6g", pn[pt], l2r_ninf[pt], ninf_kept[pt],
                    ninf_drop[pt], ninf_in[pt] > HALF_NEG ? (ninf_in[pt] - logZ) / std::log(10.0) : -999.0,
                    post);
      out += buf;
    }
    std::snprintf(buf, sizeof buf, " | post_total %.6g", all);
    return out + buf;
  }
  double l2r_bsec[PT_N] = {}, l2r_psec[PT_N] = {}, l2r_csec = 0.0, l2r_isec = 0.0;
  long long sp_events = 0, sp_cand = 0, sp_kept = 0;

  void l2r_init_scores() {
    if (l2r_beam <= 0) return;
    SecTimer t(l2r_isec);
    pr_build_static();
    std::unique_ptr<SecTimer> rtm(lpf ? new SecTimer(vi_rest_s) : nullptr);
    pr_sf_build();
    for (int r = 0; r < prA; ++r) lArm.emplace(key2(prArmXs[r], prArmXp[r]), r);
    lYwS.assign(prYy.size(), NEG_INF);
    for (int r = 0; r < prA; ++r) {
      double acc = NEG_INF;
      for (int k = prArmY0[r + 1] - 1; k >= prArmY0[r]; --k) {
        acc = log_add(acc, prYw[k]);
        lYwS[k] = acc;
      }
    }
    if (sp_cube && sp_static) l2r_spawn_cols();
    if (l2r_inc) lp_init_scores();
  }
  // --g-mode eager: G_L for every pair in P, before the scan, prefix frozen at C_pre.
  // y descending so an enclosing pair (larger right end) is built first.
  void growth_eager() {
    if (!p2_g_eager || !stem_ || l2r_beam <= 0) return;
    const std::vector<std::vector<int>>& KR = pr_pair_list(true);
    std::vector<std::pair<int, int>> pairs;
    for (int x = 0; x < n; ++x)
      for (int y : KR[x]) pairs.emplace_back(x, y);
    std::sort(pairs.begin(), pairs.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
      if (a.second != b.second) return a.second > b.second;
      return a.first > b.first;
    });
    for (const auto& [x, y] : pairs) {
      rec_root(x, y, false);
      rec_root(x, y, true);
    }
  }
  // Posterior mass of the spawned Rs states, summed over the steps (after the outside pass).
  double l2r_rs_mass() {
    outside();
    double m = 0.0;
    for (int j = 0; j < n; ++j)
      for (const Grp& g : LA.Gs[j].grps)
        for (const Ent& e : g.ents)
          if (e.a > HALF_NEG && e.b > HALF_NEG) m += std::exp(e.a + e.b - logZ);
    return m;
  }
  int l2r_arms() const { return prA; }
  size_t l2r_arm_ents() const { return pr_vxs ? prYy.size() : prEy.size(); }
  size_t l2r_arm_groups() const { return prYy.size(); }
  int larm(int xs, int xp) const {
    auto it = lArm.find(key2(xs, xp));
    return it == lArm.end() ? -1 : it->second;
  }
  // --l2r-arm-fb on: an arm xs..x' not in V scores as the nearest V arm xs''..x', xs < xs'' <= x',
  // times per extra outer pair x the best stack x.y onto (x+1).(y-1) over x's partners (with the
  // terminal penalties); with no such arm, or a score of -inf, the floor kFbFloor (finite).
  bool l2r_arm_fb = false;
  static constexpr double kFbFloor = -1e4;
  long long l2r_fb_hit = 0, l2r_fb_floor = 0;
  std::unordered_map<uint64_t, std::pair<int, double>> lFb;
  std::vector<double> lStk;
  double lstack_best(int x) {
    if (lStk.empty()) lStk.assign(n, std::nan(""));
    if (!std::isnan(lStk[x])) return lStk[x];
    double b = NEG_INF;
    for (int y : partners_right[x]) {
      if (x + 1 >= y - 1 || !E.pairs(x + 1, y - 1)) continue;
      const double s = E.stem_interior(x, x + 1, y - 1, y);
      if (s < HALF_NEG) continue;
      b = std::max(b, s + pair_bonus + E.at(x, y) - E.at(x + 1, y - 1));
    }
    return lStk[x] = b;
  }
  std::pair<int, double> larm_fb(int xs, int xp) {
    const int r = larm(xs, xp);
    if (r >= 0 || !l2r_arm_fb) return {r, 0.0};
    const uint64_t key = key2(xs, xp);
    auto it = lFb.find(key);
    if (it != lFb.end()) return it->second;
    std::pair<int, double> res{-1, 0.0};
    double extra = 0.0;
    for (int x = xs; x < xp; ++x) {
      const double s = lstack_best(x);
      if (s < HALF_NEG) break;
      extra += s;
      const int r2 = larm(x + 1, xp);
      if (r2 >= 0) {
        res = {r2, extra};
        break;
      }
    }
    return lFb.emplace(key, res).first->second;
  }
  // ---- D2 oracle substitutions (--l2r-sub a,c,d; (b) untrimmed V = --pr-v-eps 0)
  // a: K_q'(arm) = sum_{y' > q'} Qzbar(q'+1, y'-1) E_in V (not separable), cached per
  //    (q', arm) for the current step.  c: + rho(a'), the exact run's context correction of
  //    C(a'-1) (S and the U tables).  d: the exact outside of U / Us states of the exact run
  //    (U, Us; S through its B seeds), with rho divided out unless c.
  bool lsub_a = false, lsub_c = false, lsub_d = false;
  std::vector<double> lRho;
  std::unordered_map<uint64_t, double> lOU, lOUs;
  std::unordered_map<uint64_t, double> lKc;
  int lKc_q = -1;
  static uint64_t key4(int a, int b, int c, int d) {
    return (static_cast<uint64_t>(a) << 48) | (static_cast<uint64_t>(b) << 32) |
           (static_cast<uint64_t>(c) << 16) | static_cast<uint64_t>(d);
  }
  double lrho(int a) const { return lRho.empty() ? 0.0 : lRho[a]; }
  // --l2r-ctx p1: log lambda(a') from a finished pass 1 (empty: off).  Scores only: g_L's stop
  // term C(p-1) lambda(p), and the S and U-type scores times lambda(a').
  std::vector<double> lLam;
  double llam(int a) const { return lLam.empty() ? 0.0 : lLam[a]; }
  // lambda(a') = [sum_ys I1 O1 + tau W] / [sum_ys I1 E + tau W] over src's Qp items (a', ys),
  // I1 / O1 their inside / outside, E = e^beta1 C1(a'-1) Sbar_pk(ys+1) (C1: src's exterior
  // prefix), W = mean of I1 E at a'; 1 where src has no item at a'.  src.outside() is run here.
  void l2r_ctx_lambda(DP03& src, const std::vector<double>& spk, double tau,
                      std::vector<double>& lam, std::vector<int>& cnt) const {
    src.outside();
    std::vector<double> num(n, NEG_INF), den(n, NEG_INF);
    cnt.assign(n, 0);
    for (int ys = 0; ys < n; ++ys) {
      const OMap& at = src.Qp_at[ys];
      for (size_t r = 0; r < at.slots(); ++r) {
        const int a = at.keys[r];
        if (a < 0) continue;
        const double I = at.vals[r];
        num[a] = log_add(num[a], I + src.Qp_b[ix(a, ys)]);
        den[a] = log_add(den[a], I + (E.w_beta1 + src.prefix_log(a) + spk[ys + 1]));
        ++cnt[a];
      }
    }
    lam.assign(n, 0.0);
    const double lt = std::log(tau);
    for (int a = 0; a < n; ++a) {
      if (cnt[a] == 0 || den[a] < HALF_NEG) continue;
      const double tw = tau > 0.0 ? lt + den[a] - std::log(static_cast<double>(cnt[a])) : NEG_INF;
      lam[a] = log_add(num[a], tw) - log_add(den[a], tw);
    }
  }
  double lKx_sum(int r, int qp) const {
    const StemEst& st = *stem_;
    const int hi = prArmY0[r + 1];
    double acc = NEG_INF;
    for (int k = static_cast<int>(std::upper_bound(prYy.begin() + prArmY0[r], prYy.begin() + hi, qp) -
                                  prYy.begin());
         k < hi; ++k) {
      const double z = st.Qzbar(qp + 1, prYy[k] - 1);
      if (z > HALF_NEG && prYw[k] > HALF_NEG) acc = log_add(acc, z + prYw[k]);
    }
    return acc;
  }
  double lKx(int r, int qp) {
    if (lKc_q != cur_j) {
      lKc.clear();
      lKc_q = cur_j;
    }
    const uint64_t key = key2(qp, r);
    auto it = lKc.find(key);
    if (it != lKc.end()) return it->second;
    const double acc = lKx_sum(r, qp);
    lKc.emplace(key, acc);
    return acc;
  }
  // Persistent exact-sum memo for --k-exact-ub. Same key as --k-memo, value is lKx.
  double lKx_memo(int r, int qp) {
    ++lKx_calls;
    const uint64_t key = key2(qp, r);
    auto it = lKxM.find(key);
    if (it != lKxM.end()) {
      ++lKx_hits;
      return it->second;
    }
    const double v = lKx_sum(r, qp);
    lKxM.emplace(key, v);
    return v;
  }
  double lK(int r, int qp) {
    if (r < 0 || qp + 1 > n - 1) return NEG_INF;
    if (sx_k_exact) return lKx_memo(r, qp);
    if (lsub_a) return lKx(r, qp);
    if (!k_memo) return arm_phi(r, qp, pr_aq_sep ? stem_->Qzbar(0, qp) : 0.0);
    ++lKm_calls;
    const uint64_t key = key2(qp, r);
    auto it = lKm.find(key);
    if (it != lKm.end()) {
      ++lKm_hits;
      return it->second;
    }
    const double v = arm_phi(r, qp, pr_aq_sep ? stem_->Qzbar(0, qp) : 0.0);
    lKm.emplace(key, v);
    return v;
  }
  // --l2r-ctx on: the nested pre-pass's span context pr_ctx(c, ys) = e^beta1 C[c-1] Sbar[ys+1] +
  // e^beta1m Bml(c, ys) in place of the scores' exterior e^beta1 C[c-1] Sbar[ys+1], as one log
  // factor per state: kbar(c) = log sum_ys W(ys) e^(pr_ctx(c, ys) - C[c-1]) -
  // log sum_ys W(ys) e^(beta1 + Sbar[ys+1]), W the pre-pass V weight by ys of the state's arm
  // (S and R types; c = a' for S, the A pair's left end for R) or of every arm with the state's
  // xs (U types, c = a').  Cached per (c, arm) / (c, xs).
  bool l2r_ctx = false;
  std::vector<std::vector<std::pair<int, double>>> lWr, lWx;
  std::unordered_map<uint64_t, double> lKbr, lKbx;
  void l2r_ctx_init() {
    if (!l2r_ctx) return;
    if (pr_vxs) throw std::runtime_error("--pr-v-key xs: --l2r-ctx needs per-ys arm entries");
    auto add = [](std::vector<std::pair<int, double>>& v, int ys, double w) {
      for (auto& [y, s] : v)
        if (y == ys) {
          s = log_add(s, w);
          return;
        }
      v.emplace_back(ys, w);
    };
    lWr.assign(prA, {});
    lWx.assign(n, {});
    for (int r = 0; r < prA; ++r)
      for (int k = prArmE0[r]; k < prArmE0[r + 1]; ++k) {
        add(lWr[r], prEys[k], prEw[k]);
        add(lWx[prArmXs[r]], prEys[k], prEw[k]);
      }
  }
  double lkbar(const std::vector<std::pair<int, double>>& w, int c) const {
    const StemEst& st = *stem_;
    const double pc = prefix_log(c);
    double num = NEG_INF, den = NEG_INF;
    for (const auto& [ys, lw] : w) {
      if (ys <= c) continue;
      const double ex = E.w_beta1 + st.prS[ys + 1];
      den = log_add(den, lw + ex);
      num = log_add(num, lw + log_add(ex, E.w_beta1m + st.Bml(c, ys) - pc));
    }
    return den > HALF_NEG ? num - den : 0.0;
  }
  double lkb_arm(int c, int r) {
    if (!l2r_ctx || r < 0) return 0.0;
    const uint64_t key = key2(c, r);
    auto it = lKbr.find(key);
    return it != lKbr.end() ? it->second : lKbr.emplace(key, lkbar(lWr[r], c)).first->second;
  }
  double lkb_xs(int c, int xs) {
    if (!l2r_ctx) return 0.0;
    const uint64_t key = key2(c, xs);
    auto it = lKbx.find(key);
    return it != lKbx.end() ? it->second : lKbx.emplace(key, lkbar(lWx[xs], c)).first->second;
  }

  double scoreS_k(int r, int j) {
    ++rfv_look_s;
    const double k = lK(r, j);
    if (k < HALF_NEG) {
      ++rfv_miss_s;
      const int mk = k_map_kind(r, j);
      note_kind(vm_i, vm_ii, vm_iii, ST_S, mk);
      note_kind(km_i, km_ii, km_iii, ST_S, mk);
      if (mk == 2) miss_k(r, j, ST_S);
    }
    return k;
  }
  double l2r_scoreS(uint64_t key, double f, int j) {
    if (f < HALF_NEG) return NEG_INF;
    const int a = k3a(key), xs = k3b(key), xp = k3c(key);
    if (!lsub_d && l2r_arm_fb) {
      const auto [r, fx] = larm_fb(xs, xp);
      if (r >= 0 && r != larm(xs, xp)) ++l2r_fb_hit;
      const double k = scoreS_k(r, j);
      if (k < HALF_NEG) {
        ++l2r_fb_floor;
        return f + prefix_log(a) + kFbFloor;
      }
      return f + (k + fx) + prefix_log(a) + (lsub_c ? lrho(a) : lkb_arm(a, r)) + llam(a);
    }
    if (!lsub_d && !lLam.empty()) {
      const int r = larm(xs, xp);
      return f + scoreS_k(r, j) + prefix_log(a) + (lsub_c ? lrho(a) : lkb_arm(a, r)) + lLam[a];
    }
    if (!lsub_d) {
      const int r = larm(xs, xp);
      if (r < 0) ++l2r_noarm[PT_S];
      return f + scoreS_k(r, j) + prefix_log(a) + (lsub_c ? lrho(a) : lkb_arm(a, r));
    }
    const StemEst& st = *stem_;
    double acc = NEG_INF;
    for (int y : partners_right[xp]) {
      if (y <= j || !E.end_ok(xp, y)) continue;
      auto it = lOUs.find(key4(a, xs, xp, y));
      if (it == lOUs.end()) continue;
      const double fill = lsub_a ? st.Qzbar(j + 1, y - 1) : st.Qzbar(0, y - 1) - st.Qzbar(0, j);
      if (!(fill > HALF_NEG) || !(fill < 1e300)) continue;
      acc = log_add(acc, fill + (E.w_stem + E.at(xp, y) + pair_bonus) + it->second);
    }
    return acc > HALF_NEG ? f + acc - (lsub_c ? 0.0 : lrho(a)) : NEG_INF;
  }

  // ---- --l2r-spawn inc scores (log; values exact, scores only rank the states of one step j).
  // Pz(i) = Qzbar(0, i), the pre-pass prefix gap filler: Qz(i, j) ~ Pz(j) / Pz(i-1), an
  // overestimate, used only as a ratio whose 1 / Pz(j) is common to the states at j.
  //   H(a2, j)[q'] = sum_{c1 > j, a2.c1} Pz(c1-1) E_A(a2, c1) Lambda_L(a2, c1)[q'], Lambda_L the
  //   truncated L vector of rec_root(a2, c1, inner) (C[p'-1] folded in); per a2 a suffix list
  //   over c1, each suffix kept as its top lp_hk entries by H uL(q') (ntrUL; q' with uL = 0
  //   dropped).  Built at the scan's step a2 (the roots read C left of a2 only).
  //   P3: F - Pz(j) + log sum_q' H(a2, j)[q'] K_q'(arm) + ctx(a2, arm)
  //   P2: F - Pz(j) + Aopen(a2, j) + ArmOpen(xs, j) + ctx(a2, xs)
  //   P1: F - Pz(j) + Aopen(a2, j) + Bopen(j)
  //   Aopen = log sum_q' H[q'] / Pz(q') (--l2r-p-aopen q; raw: sum_q' H[q']); ArmOpen(xs, j) =
  //   log sum_{V arms xs..x', x' >= j} W(xs, x') / Pz(x'), W = sum_y' Pz(y'-1) E_in V (prYw);
  //   Bopen(j) = log sum_{xs > j} Pz(xs-1) ArmOpen(xs, xs).
  // A score of -inf (no floor) is counted in lp_miss[table][reason]: 0 arm not in V, 1 H(a2, j)
  // empty, 2 no kept q' with K_q'(arm) > 0 (every y' <= q'), 3 ArmOpen / Bopen empty.
  int lp_hk = 16;
  bool lp_aopen_raw = false;
  struct LH {
    std::vector<int> c1;                                 // ascending
    std::vector<std::vector<std::pair<int, double>>> v;  // suffix c1[k..], q' ascending
    std::vector<double> aopen;
    std::vector<const RN*> rt;                           // rec_root(a2, c1[k], true)
    std::vector<std::vector<double>> l;                  // --l2r-p3-lin: e^(v - lpz(q') - lm)
    std::vector<double> lm;
  };
  // --l2r-p3-lin on (separable arms, no span, no sub a): the P3 dot product in linear space,
  // H as e^(v - lpz(q') - m) with one m per list, K as e^(prSF - prSF[first of arm]).
  bool lp_lin = false;
  std::vector<double> lSFl, lSFm;
  bool lp_lin_ok() const { return lp_lin && pr_aq_sep && pr_pspan <= 0 && !lsub_a; }
  void lp_lin_of(const std::vector<std::pair<int, double>>& v, std::vector<double>& l, double& m) const {
    l.resize(v.size());
    m = NEG_INF;
    for (size_t i = 0; i < v.size(); ++i) {
      const int q = v[i].first;
      const double z = q + 1 <= n - 1 ? lpz(q) : NEG_INF;
      l[i] = z > HALF_NEG && v[i].second > HALF_NEG ? v[i].second - z : NEG_INF;
      m = std::max(m, l[i]);
    }
    for (double& x : l) x = x > HALF_NEG ? std::exp(x - m) : 0.0;
  }
  std::vector<LH> lH;
  std::vector<char> lHdone;
  std::vector<std::vector<std::pair<int, double>>> lArmOpen;  // per xs: (x', suffix), x' ascending
  std::vector<double> lBopen;
  long long lp_h_ents = 0, lp_h_roots = 0, lp_h_lists = 0;
  double lp_h_sec = 0.0;
  long long lp_miss[3][4] = {};
  double lpz(int i) const { return stem_->Qzbar(0, i); }
  std::vector<std::pair<int, double>> lamL(const RN& r) {
    if (g_ntr_K > 0) return r.T[0];
    if (g_tr_K > 0) return rec_trunc(r, 0).L;
    std::vector<double> acc(n, NEG_INF);
    for (int m = 0; m <= g_ml_L && m < rnL; ++m) {
      const RV& a = r.L[m];
      for (size_t i = 0; i < a.v.size(); ++i)
        if (a.v[i] > 0.0) acc[a.lo + i] = log_add(acc[a.lo + i], a.mx + std::log(a.v[i]));
    }
    std::vector<std::pair<int, double>> o;
    for (int q = 0; q < n; ++q)
      if (acc[q] > HALF_NEG) o.emplace_back(q, acc[q]);
    return o;
  }
  const LH& lh_get(int a2) {
    if (lHdone[a2]) return lH[a2];
    if (a2 > cur_j + 1) {
      std::fprintf(stderr, "[l2r] abort: H(%d) built at scan %d (needs C left of a2)\n", a2, cur_j);
      std::abort();
    }
    SecTimer t(lp_h_sec);
    ProfReg reg(*this, 3);
    lHdone[a2] = 1;
    ntr_init();
    LH& h = lH[a2];
    for (int c1 : partners_right[a2])
      if (c1 >= a2 + 2 && E.end_ok(a2, c1)) h.c1.push_back(c1);
    const size_t m = h.c1.size();
    h.v.resize(m);
    h.aopen.assign(m, NEG_INF);
    h.rt.assign(m, nullptr);
    std::vector<std::pair<int, double>> cum, mrg;
    std::vector<std::pair<double, int>> rk;
    const long long hb0 = lpf_help_built, hc0 = lpf_hg_clears, rr0 = rn_roots;
    const double hs0 = lpf_help_sec;
    for (size_t k = m; k-- > 0;) {
      const int c1 = h.c1[k];
      const double base = lpz(c1 - 1) + (pair_bonus + E.w_stem + E.at(a2, c1));
      std::chrono::steady_clock::time_point ht0;
      if (lpf) ht0 = std::chrono::steady_clock::now();
      h.rt[k] = &rec_root(a2, c1, true);
      std::vector<std::pair<int, double>> add = lamL(*h.rt[k]);
      if (lpf) {
        const auto ht1 = std::chrono::steady_clock::now();
        lpf_h_root_sec += std::chrono::duration<double>(ht1 - ht0).count();
        ht0 = ht1;
        lpf_h_lam += static_cast<long long>(add.size());
        lpf_h_merge += static_cast<long long>(cum.size() + add.size());
        if (cum.size() + add.size() > static_cast<size_t>(lp_hk)) {
          ++lpf_h_trunc;
          lpf_h_trunc_in += static_cast<long long>(cum.size() + add.size());
        }
      }
      ++lp_h_roots;
      mrg.clear();
      size_t x = 0, y = 0;
      while (x < cum.size() || y < add.size()) {
        if (y < add.size() && ntrUL[add[y].first] < HALF_NEG) {
          ++y;
          continue;
        }
        if (y == add.size() || (x < cum.size() && cum[x].first < add[y].first)) {
          mrg.push_back(cum[x++]);
        } else if (x == cum.size() || add[y].first < cum[x].first) {
          mrg.emplace_back(add[y].first, base + add[y].second);
          ++y;
        } else {
          mrg.emplace_back(cum[x].first, log_add(cum[x].second, base + add[y].second));
          ++x;
          ++y;
        }
      }
      if (mrg.size() > static_cast<size_t>(lp_hk)) {
        rk.clear();
        for (size_t i = 0; i < mrg.size(); ++i) rk.emplace_back(mrg[i].second + ntrUL[mrg[i].first], static_cast<int>(i));
        std::nth_element(rk.begin(), rk.begin() + lp_hk, rk.end(),
                         [](const std::pair<double, int>& p, const std::pair<double, int>& q) {
                           return p.first > q.first || (p.first == q.first && p.second < q.second);
                         });
        std::vector<char> keep(mrg.size(), 0);
        for (int i = 0; i < lp_hk; ++i) keep[rk[i].second] = 1;
        size_t w = 0;
        for (size_t i = 0; i < mrg.size(); ++i)
          if (keep[i]) mrg[w++] = mrg[i];
        mrg.resize(w);
      }
      cum.swap(mrg);
      h.v[k] = cum;
      if (lp_lin_ok()) {
        h.l.resize(m);
        h.lm.resize(m);
        lp_lin_of(cum, h.l[k], h.lm[k]);
      }
      double ao = NEG_INF;
      for (const auto& [qp, hv] : cum) ao = log_add(ao, lp_aopen_raw ? hv : hv - lpz(qp));
      h.aopen[k] = ao;
      lp_h_ents += static_cast<long long>(cum.size());
      if (lpf) lpf_h_merge_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - ht0).count();
    }
    if (lpf) {
      lpf_h_newroots += rn_roots - rr0;
      lpf_h_helpers += lpf_help_built - hb0;
      lpf_h_clears += lpf_hg_clears - hc0;
      lpf_h_help_sec += lpf_help_sec - hs0;
    }
    ++lp_h_lists;
    return h;
  }
  // index of the first c1 > j in H(a2), or -1
  int lh_at(const LH& h, int j) const {
    const auto it = std::upper_bound(h.c1.begin(), h.c1.end(), j);
    return it == h.c1.end() ? -1 : static_cast<int>(it - h.c1.begin());
  }
  // --l2r-p3-hex on: the P3 score's continuation with the pre-pass filler in place of the prefix
  // ratio, H_ex(a2, j)[q'] = sum_{c1 > j} Qzbar(j+1, c1-1) E_A(a2, c1) Lambda_L(a2, c1)[q'] (the
  // backward filler recursion's value), top lp_hk entries as H ranks them; built per (a2, j) at
  // step j, kept for the current j only.
  bool lp_hex = false;
  int lhx_j = -1;
  struct LX {
    std::vector<std::pair<int, double>> v;
    std::vector<double> l;
    double m = NEG_INF;
  };
  std::unordered_map<int, LX> lHx;
  std::vector<double> lhxA;
  std::vector<int> lhxT;
  long long lp_hx_builds = 0, lp_hx_terms = 0, lp_hx_ents = 0;
  double lp_hx_sec = 0.0;
  const LX& lh_ex(const LH& h, int a2, int j) {
    if (j != lhx_j) {
      lHx.clear();
      lhx_j = j;
    }
    auto it = lHx.find(a2);
    if (it != lHx.end()) return it->second;
    SecTimer tm(lp_hx_sec);
    if (lhxA.empty()) lhxA.assign(n, NEG_INF);
    lhxT.clear();
    const int k0 = lh_at(h, j);
    for (int k = k0 < 0 ? static_cast<int>(h.c1.size()) : k0; k < static_cast<int>(h.c1.size()); ++k) {
      const int c1 = h.c1[k];
      const double z = stem_->Qzbar(j + 1, c1 - 1);
      if (z < HALF_NEG) continue;
      const double base = z + (pair_bonus + E.w_stem + E.at(a2, c1));
      std::vector<std::pair<int, double>> own;
      const std::vector<std::pair<int, double>>& lam = g_ntr_K > 0 ? h.rt[k]->T[0] : (own = lamL(*h.rt[k]));
      for (const auto& [q, g] : lam) {
        if (ntrUL[q] < HALF_NEG) continue;
        if (lhxA[q] < HALF_NEG) lhxT.push_back(q);
        lhxA[q] = log_add(lhxA[q], base + g);
        ++lp_hx_terms;
      }
    }
    std::sort(lhxT.begin(), lhxT.end());
    std::vector<std::pair<int, double>> o;
    o.reserve(lhxT.size());
    for (int q : lhxT) {
      o.emplace_back(q, lhxA[q]);
      lhxA[q] = NEG_INF;
    }
    if (o.size() > static_cast<size_t>(lp_hk)) {
      std::vector<std::pair<double, int>> rk;
      for (size_t i = 0; i < o.size(); ++i) rk.emplace_back(o[i].second + ntrUL[o[i].first], static_cast<int>(i));
      std::nth_element(rk.begin(), rk.begin() + lp_hk, rk.end(),
                       [](const std::pair<double, int>& p, const std::pair<double, int>& q) {
                         return p.first > q.first || (p.first == q.first && p.second < q.second);
                       });
      std::vector<char> keep(o.size(), 0);
      for (int i = 0; i < lp_hk; ++i) keep[rk[i].second] = 1;
      size_t w = 0;
      for (size_t i = 0; i < o.size(); ++i)
        if (keep[i]) o[w++] = o[i];
      o.resize(w);
    }
    ++lp_hx_builds;
    lp_hx_ents += static_cast<long long>(o.size());
    LX x;
    if (lp_lin_ok()) lp_lin_of(o, x.l, x.m);
    x.v = std::move(o);
    return lHx.emplace(a2, std::move(x)).first->second;
  }
  double l_armopen(int xs, int j) const {
    if (xs >= n) return NEG_INF;
    const auto& v = lArmOpen[xs];
    const auto it = std::lower_bound(v.begin(), v.end(), j,
                                     [](const std::pair<int, double>& e, int x) { return e.first < x; });
    return it == v.end() ? NEG_INF : it->second;
  }
  void lp_init_scores() {
    if (lp_lin_ok()) {
      lSFl.assign(prSF.size(), 0.0);
      lSFm.assign(prA, NEG_INF);
      for (int r = 0; r < prA; ++r) {
        if (prArmY0[r] >= prArmY0[r + 1]) continue;
        const double m = prSF[prArmY0[r]];
        if (m < HALF_NEG) continue;
        lSFm[r] = m;
        for (int k = prArmY0[r]; k < prArmY0[r + 1]; ++k)
          if (prSF[k] > HALF_NEG) lSFl[k] = std::exp(prSF[k] - m);
      }
    }
    lH.assign(n, LH{});
    lHdone.assign(n, 0);
    lArmOpen.assign(n, {});
    for (int xs = 0; xs < n; ++xs) {
      auto& v = lArmOpen[xs];
      for (int r = prXsFirst[xs]; r < prXsFirst[xs + 1]; ++r) {
        double w = NEG_INF;
        for (int k = prArmY0[r]; k < prArmY0[r + 1]; ++k) {
          const double z = lpz(prYy[k] - 1);
          if (z > HALF_NEG) w = log_add(w, z + prYw[k]);
        }
        if (w > HALF_NEG) v.emplace_back(prArmXp[r], w - lpz(prArmXp[r]));
      }
      for (size_t k = v.size(); k-- > 1;) v[k - 1].second = log_add(v[k - 1].second, v[k].second);
    }
    lBopen.assign(n + 1, NEG_INF);
    for (int j = n - 2; j >= 0; --j) lBopen[j] = log_add(lBopen[j + 1], lpz(j) + l_armopen(j + 1, j + 1));
  }
  std::chrono::steady_clock::time_point lpf_t0;
  void lpf_lap(double& acc) {
    const auto t1 = std::chrono::steady_clock::now();
    acc += std::chrono::duration<double>(t1 - lpf_t0).count();
    lpf_t0 = t1;
  }
  double lp_score(int pt, uint64_t key, double f, int j) {
    if (f < HALF_NEG) return NEG_INF;
    const int t = pt - PT_P1, a2 = k3a(key), xs = k3b(key);
    if (lpf) {
      ++lpf_sc_calls[t];
      lpf_t0 = std::chrono::steady_clock::now();
    }
    int arm = -1;
    if (pt == PT_P3) {
      arm = larm(xs, k3c(key));
      if (lpf) {
        ++lpf_sc_osteps[t];
        lpf_lap(lpf_sc_osec[t]);
      }
      if (arm < 0) {
        ++lp_miss[t][0];
        return NEG_INF;
      }
    }
    const LH& h = lh_get(a2);
    const int k = lh_at(h, j);
    if (lpf) {
      lpf_sc_hsteps[t] += 1 + lpf_log2(h.c1.size());
      lpf_lap(lpf_sc_hsec[t]);
    }
    if (k < 0 || h.v[k].empty()) {
      ++lp_miss[t][1];
      const int hk = h_map_kind(a2, j);
      note_kind(hm_i, hm_ii, hm_iii, t, hk);
      if (hk == 2 && (t == ST_P3 || t == ST_S)) miss_take(2, t, a2, j, -1);
      return NEG_INF;
    }
    if (pt == PT_P3) {
      const LX* hx = lp_hex ? &lh_ex(h, a2, j) : nullptr;
      const std::vector<std::pair<int, double>>& hl = hx ? hx->v : h.v[k];
      double s = NEG_INF;
      if (lp_lin_ok()) {
        const std::vector<double>& ll = hx ? hx->l : h.l[k];
        const int hi = prArmY0[arm + 1];
        const int* yy = prYy.data();
        int kk = prArmY0[arm];
        double a = 0.0;
        for (size_t i = 0; i < hl.size(); ++i) {
          if (ll[i] == 0.0) continue;
          kk = static_cast<int>(std::upper_bound(yy + kk, yy + hi, hl[i].first) - yy);
          if (kk >= hi) break;
          a += ll[i] * lSFl[kk];
        }
        if (a > 0.0) s = (hx ? hx->m : h.lm[k]) + lSFm[arm] + std::log(a);
      } else {
        for (const auto& [qp, hv] : hl) {
          const double kq = lK(arm, qp);
          if (kq > HALF_NEG) s = log_add(s, hv + kq);
        }
      }
      if (lpf) {
        const long long m = static_cast<long long>(hl.size());
        lpf_sc_dot[t] += m;
        lpf_sc_kcalls[t] += m;
        lpf_sc_ksteps[t] += m * lpf_log2(static_cast<size_t>(prArmY0[arm + 1] - prArmY0[arm]));
        lpf_lap(lpf_sc_ksec[t]);
      }
      if (s < HALF_NEG) {
        ++lp_miss[t][2];
        const int qp = hl.empty() ? j : hl.front().first;
        const int mk = k_map_kind(arm, qp);
        note_kind(km_i, km_ii, km_iii, t, mk);
        if (mk == 2) miss_k(arm, qp, t);
        return NEG_INF;
      }
      if (lp_hex) return f + s + lkb_arm(a2, arm);
      return f - lpz(j) + s + lkb_arm(a2, arm);
    }
    if (lpf) {
      lpf_sc_osteps[t] += pt == PT_P2 && xs < n ? lpf_log2(lArmOpen[xs].size()) : 1;
      lpf_sc_dot[t] += 1;
    }
    const double o = pt == PT_P2 ? l_armopen(xs, j) : lBopen[j];
    if (lpf) lpf_lap(lpf_sc_osec[t]);
    if (o < HALF_NEG) {
      ++lp_miss[t][3];
      return NEG_INF;
    }
    return f - lpz(j) + h.aopen[k] + o + (pt == PT_P2 ? lkb_xs(a2, xs) : 0.0);
  }

  // Pre-pass suffix recursion (build_principled's prS) with qp[y] = (j, log Qp(y, j)) items as
  // one more paired element, closed as fill_C does: Qp e^beta1 Sbar[j+1].  qb: (j, log Qb(y, j))
  // items in place of the pre-pass nested Qb.
  using QpItems = std::vector<std::vector<std::pair<int, double>>>;
  std::vector<double> l2r_pk_suffix(const QpItems& qp, const QpItems* qb = nullptr) const {
    const StemEst& st = *stem_;
    std::vector<double> S2(n + 1, NEG_INF);
    S2[n] = 0.0;
    auto paired = [&](int y) {
      double acc = NEG_INF;
      if (qb) {
        for (const auto& it : (*qb)[y])
          if (E.end_ok(y, it.first)) acc = log_add(acc, it.second + E.at(y, it.first) + S2[it.first + 1]);
      } else {
        for (int j = y + 1; j < n; ++j) {
          const double q = st.QbIn(y, j);
          if (q > HALF_NEG && E.end_ok(y, j)) acc = log_add(acc, q + E.at(y, j) + S2[j + 1]);
        }
      }
      for (const auto& it : qp[y]) acc = log_add(acc, it.second + E.w_beta1 + S2[it.first + 1]);
      return acc;
    };
    if (dz) {
      std::vector<double> P(n + 1, NEG_INF);
      P[n] = 0.0;
      for (int y = n - 1; y >= 0; --y) {
        P[y] = paired(y);
        double s = P[y];
        for (int t = y; t < n; ++t) s = log_add(s, E.run(y, t) + P[t + 1]);
        S2[y] = s;
      }
    } else {
      for (int y = n - 1; y >= 0; --y) S2[y] = log_add(S2[y + 1], paired(y));
    }
    return S2;
  }
  // --l2r-p1 R: Sbar_pk from a finished beamed pseudoknot pass src: its stored Qp items with
  // the pre-pass nested Qb (Sn) or with src's own Qb items (Sp; Sp[0] = src's log Z, fill_C).
  static QpItems l2r_items(const std::vector<OMap>& at, int n, long long& cnt) {
    QpItems it(n);
    cnt = 0;
    for (int j = 0; j < n; ++j)
      for (size_t r = 0; r < at[j].slots(); ++r)
        if (at[j].keys[r] >= 0 && at[j].vals[r] > HALF_NEG) {
          it[at[j].keys[r]].emplace_back(j, at[j].vals[r]);
          ++cnt;
        }
    return it;
  }
  void l2r_sfx_pass(const DP03& src, std::vector<double>& Sn, std::vector<double>& Sp,
                    long long& qp_items, long long& qb_items) const {
    const QpItems qp = l2r_items(src.Qp_at, n, qp_items);
    const QpItems qb = l2r_items(src.Qb_at, n, qb_items);
    Sn = l2r_pk_suffix(qp);
    Sp = l2r_pk_suffix(qp, &qb);
  }
  // One pass's time split (s): init (l2r_init_scores: arm entries + prYw + prSF), growth and
  // selection (counted where they run), Rs = spawn (build + score), all tables' score and build
  // (Rs included), close; other = pass CPU - init - build - score - close (nested tables).
  std::string l2r_prof(double pass_cpu) const {
    double bs = 0.0, ps = 0.0;
    for (int t = 0; t < PT_N; ++t) {
      bs += l2r_bsec[t];
      ps += l2r_psec[t];
    }
    char b[512];
    std::snprintf(b, sizeof b,
                  "cpu %.2f init %.2f static %.2f grow %.2f sel %.2f spawn %.2f (cand %lld) "
                  "score %.2f build %.2f close %.2f other %.2f",
                  pass_cpu, l2r_isec, pr_static_sec, pr_grow_sec, pr_sel_sec,
                  l2r_bsec[PT_RS] + l2r_psec[PT_RS], sp_cand, ps, bs, l2r_csec,
                  pass_cpu - l2r_isec - bs - ps - l2r_csec);
    return b;
  }
  static int pct_p99(const long long* h, int nbin) {
    long long tot = 0;
    for (int i = 0; i < nbin; ++i) tot += h[i];
    if (tot <= 0) return -1;
    const long long need = (tot * 99 + 99) / 100;
    long long c = 0;
    for (int i = 0; i < nbin; ++i) {
      c += h[i];
      if (c >= need) return i;
    }
    return nbin - 1;
  }
  void scan_ct_dump() const {
    auto grp = [&](int a, int b, long long& pre, long long& kept, long long& mx, long long& bind,
                   long long& pos, int& beam) {
      pre = kept = mx = bind = pos = 0;
      beam = ltb(a);
      for (int t = a; t <= b; ++t) {
        pre += l2r_pre[t];
        kept += l2r_kept[t];
        mx = std::max(mx, l2r_kmax[t]);
        bind += l2r_bind[t];
        pos += l2r_pos[t];
      }
    };
    auto one = [&](const char* name, int a, int b) {
      long long pre, kept, mx, bind, pos;
      int beam;
      grp(a, b, pre, kept, mx, bind, pos, beam);
      const double per = pos > 0 ? static_cast<double>(pre) / static_cast<double>(pos) : 0.0;
      std::fprintf(stderr,
                   " %s beam %d pre %lld kept %lld max %lld bind %lld pos %lld cand_per_pos %.3f",
                   name, beam, pre, kept, mx, bind, pos, per);
    };
    long long m_n = 0, m_sum = 0, list_n = 0, list_sum = 0;
    for (int i = 0; i < 8; ++i) {
      m_n += sx_m_h[i];
      m_sum += static_cast<long long>(i) * sx_m_h[i];
    }
    for (int i = 0; i < 65; ++i) {
      list_n += sx_list_h[i];
      list_sum += static_cast<long long>(i) * sx_list_h[i];
    }
    const double k_per = sx_rkey_n > 0 ? static_cast<double>(sx_k_src[2]) / static_cast<double>(sx_rkey_n)
                                       : 0.0;
    long long rpre = 0;
    for (int t = PT_R; t <= PT_RRS; ++t) rpre += l2r_pre[t];
    const double rkey_per = rpre > 0 ? static_cast<double>(sx_rkey_n) / static_cast<double>(rpre) : 0.0;
    std::fprintf(stderr, "[scan-ct] n %d", n);
    one("P1", PT_P1, PT_P1);
    one("P2", PT_P2, PT_P2);
    one("P3", PT_P3, PT_P3);
    one("R", PT_R, PT_RRS);
    one("S", PT_S, PT_S);
    one("U", PT_U, PT_URS);
    std::fprintf(stderr,
                 " rkey_calls %lld r_cand %lld rkey_per_cand %.6f k_rkey %lld k_per_rkey %.3f "
                 "k_root %lld k_spawn %lld k_n %lld m_n %lld m_mean %.3f m_p99 %d "
                 "qlist_n %lld qlist_mean %.3f qlist_p99 %d\n",
                 sx_rkey_n, rpre, rkey_per, sx_k_src[2], k_per, sx_k_src[1], sx_k_src[3], sx_k_n,
                 m_n, m_n > 0 ? static_cast<double>(m_sum) / static_cast<double>(m_n) : 0.0,
                 pct_p99(sx_m_h, 8), list_n,
                 list_n > 0 ? static_cast<double>(list_sum) / static_cast<double>(list_n) : 0.0,
                 pct_p99(sx_list_h, 65));
    if (k_memo) {
      const double bytes =
          static_cast<double>(lKm.size()) * static_cast<double>(sizeof(uint64_t) + sizeof(double)) +
          static_cast<double>(lKm.bucket_count()) * static_cast<double>(sizeof(void*));
      std::fprintf(stderr,
                   "[k-memo] calls %lld hits %lld hit_rate %.6f entries %zu bytes %.0f k_s %.4f\n",
                   lKm_calls, lKm_hits,
                   lKm_calls > 0 ? static_cast<double>(lKm_hits) / static_cast<double>(lKm_calls) : 0.0,
                   lKm.size(), bytes, sx_k_sec);
    }
    if (k_exact_ub || sx_k_exact) {
      const double bytes = static_cast<double>(lKxM.size()) *
                               static_cast<double>(sizeof(uint64_t) + sizeof(double)) +
                           static_cast<double>(lKxM.bucket_count()) * static_cast<double>(sizeof(void*));
      const double per = n > 0 ? static_cast<double>(sx_kx_sigma) / static_cast<double>(n) : 0.0;
      std::fprintf(stderr,
                   "[k-exact] sigma %lld per_pos %.3f kx_calls %lld kx_hits %lld hit_rate %.6f "
                   "entries %zu bytes %.0f kx_s %.4f cols %lld mis %lld\n",
                   sx_kx_sigma, per, lKx_calls, lKx_hits,
                   lKx_calls > 0 ? static_cast<double>(lKx_hits) / static_cast<double>(lKx_calls) : 0.0,
                   lKxM.size(), bytes, sx_kx_sec, k_exact_cols, k_exact_mis);
    }
  }
  void profile_dump(const char* sid, double wall, double qz_pre) const {
    if (!prof_on) return;
    auto bucket = [&](int a, int b) {
      double s = 0.0;
      for (int t = a; t <= b; ++t) s += l2r_bsec[t] + l2r_psec[t];
      return s;
    };
    double rscore = 0.0;
    for (int t = PT_R; t <= PT_RRS; ++t) rscore += l2r_psec[t];
    const double rkey = std::max(0.0, rscore - prof_grow_reg[1]);
    long long pairs = 0;
    if (p2_g_mode > 0 && stem_)
      for (const std::vector<int>& kr : stem_->prKeptRightG)
        pairs += static_cast<long long>(kr.size());
    std::fprintf(stderr,
                 "[profile] %s wall %.3f qz_pre %.4f qz_scan %.4f v_s %.4f grow_s %.4f "
                 "grow_R %.4f grow_spawn %.4f grow_H %.4f grow_other %.4f rkey_s %.4f "
                 "P1 %.4f P2 %.4f P3 %.4f R %.4f U %.4f S %.4f close %.4f init %.4f "
                 "vectors %lld roots %lld closure %lld steps %lld steps_ml %lld help_steps %lld "
                 "pairs_P %lld\n",
                 sid, wall, qz_pre, prof_qz_sec, pr_static_sec, pr_grow_sec, prof_grow_reg[1],
                 prof_grow_reg[2], prof_grow_reg[3], prof_grow_reg[0], rkey, bucket(PT_P1, PT_P1),
                 bucket(PT_P2, PT_P2), bucket(PT_P3, PT_P3), bucket(PT_R, PT_RRS),
                 bucket(PT_U, PT_URS), bucket(PT_S, PT_S), l2r_csec, l2r_isec, rn_nodes + rn_roots,
                 rn_roots, static_cast<long long>(prRG.size()), st_kind[0][0] + st_kind[0][1],
                 st_kind[0][1], prof_help_steps, pairs);
    // Exclusive scan split. R and P3 above still include the growth built inside them.
    // Growth charged past a bucket is not removed twice; `over` is that excess and `rest`
    // is whatever of the scan none of these timers owns. v_s and init sit before the scan.
    const double P1 = bucket(PT_P1, PT_P1);
    const double P2 = bucket(PT_P2, PT_P2);
    const double P3b = bucket(PT_P3, PT_P3);
    const double Rb = bucket(PT_R, PT_RRS);
    const double Ub = bucket(PT_U, PT_URS);
    const double Sb = bucket(PT_S, PT_S);
    const double gR = std::min(prof_grow_reg[1], Rb);
    const double gH = std::min(prof_grow_reg[3], P3b);
    const double over = (prof_grow_reg[1] - gR) + (prof_grow_reg[3] - gH);
    const double rest = wall - prof_qz_sec - pr_grow_sec - P1 - P2 - (P3b - gH) - (Rb - gR) - Ub -
                        Sb - l2r_csec;
    std::fprintf(stderr,
                 "[profile-x] %s scan %.3f qz %.4f grow %.4f P1 %.4f P2 %.4f P3 %.4f R %.4f "
                 "U %.4f S %.4f close %.4f rest %.4f over %.4f\n",
                 sid, wall, prof_qz_sec, pr_grow_sec, P1, P2, P3b - gH, Rb - gR, Ub, Sb, l2r_csec,
                 rest, over);
    const double spawn = l2r_bsec[PT_RS];
    const double ext = l2r_bsec[PT_R] + l2r_bsec[PT_U] + l2r_bsec[PT_S];
    const double nest = l2r_bsec[PT_P1] + l2r_bsec[PT_P2] + l2r_bsec[PT_P3];
    const double km = sx_k_n > 0 ? static_cast<double>(sx_m_sum) / static_cast<double>(sx_k_n) : 0.0;
    std::fprintf(stderr,
                 "[profile-s] %s spawn %.4f ext %.4f nest %.4f k %.4f k_n %lld k_m %.3f rkey %.4f "
                 "per_nt %.6f\n",
                 sid, spawn, ext, nest, sx_k_sec, sx_k_n, km, sx_rkey_sec,
                 n > 0 ? wall / static_cast<double>(n) : 0.0);
  }
  // The pass-1 copy: every beam (nested tables and l2r) capped at B.
  void l2r_cap_beam(int B) {
    beam = beam > 0 ? std::min(beam, B) : B;
    for (int t = 0; t < BT_N; ++t)
      if (tbeam[t] <= 0 || tbeam[t] > B) tbeam[t] = B;
    l2r_beam = B;
  }
  // After st.prS is replaced: prYw (per arm and y', log sum_ys w e^beta1 Sbar[ys+1]) again.
  void l2r_yw_rebuild() {
    if (pr_vxs) throw std::runtime_error("--pr-v-key xs: no per-ys arm entries to rebuild prYw from");
    const StemEst& st = *stem_;
    std::fill(prYw.begin(), prYw.end(), NEG_INF);
    for (int r = 0; r < prA; ++r) {
      int g = prArmY0[r];
      for (int k = prArmE0[r]; k < prArmE0[r + 1]; ++k) {
        if (k > prArmE0[r] && prEy[k] != prEy[k - 1]) ++g;
        prYw[g] = log_add(prYw[g], prEw[k] + (E.w_beta1 + st.prS[prEys[k] + 1]));
      }
    }
    prSF.clear();
  }
  bool l2r_stem_is(const StemEst* s) const { return s == stem_; }
  int l2r_arm_entries() const { return static_cast<int>(l2r_arm_ents()); }

  // --l2r-split on (with --l2r-diag): Qp_b by consumer, QpS[0] exterior (beta1), [1] multiloop
  // (beta1m: Qm / Qmp branch, PK inside a pair via PKM), [2] gap (beta1p: Qz gap, S gap).  The
  // l2r states' outsides are linear in the Qp_b read at close B, so an outside pass reading
  // QpS[c] there gives every state's share by consumer (lspB[c]); a last normal pass restores b.
  bool l2r_split = false;
  bool lsp_rec = false;
  int lsp_use = -1;
  std::vector<double> QpS[3];
  std::vector<std::vector<double>> lspB[3][PT_N];
  void qps(int c, int i, int j, double v) {
    if (lsp_rec) add2b(QpS[c], i, j, v);
  }
  SMap& l2r_smap(int pt, int j) {
    return pt == PT_S ? LS[j] : pt == PT_P1 ? LP1[j] : pt == PT_P2 ? LP2[j] : LP3[j];
  }
  Tab4& l2r_tab(int pt, int j) {
    LFam& F = pt < PT_U ? LA : LB;
    const int r = pt - (pt < PT_U ? PT_R : PT_U);
    return r == 0 ? F.G[j] : r == 1 ? F.Gs[j] : r == 2 ? F.Gls[j] : r == 3 ? F.Glss[j] : F.Grs[j];
  }
  void l2r_reset_b() {
    for (int j = 0; j < n; ++j) {
      for (int pt = 0; pt < PT_S; ++pt)
        for (Grp& g : l2r_tab(pt, j).grps)
          for (Ent& e : g.ents) e.b = NEG_INF;
      for (SSt& s : LS[j].v) s.bP = s.bU1 = s.bU2 = s.bc = s.bF = NEG_INF;
      if (l2r_inc)
        for (std::vector<SMap>* T : {&LP1, &LP2, &LP3})
          for (SSt& s : (*T)[j].v) s.bP = s.bU1 = s.bU2 = s.bc = s.bF = NEG_INF;
    }
  }
  void l2r_collect_b(int c) {
    for (int pt = 0; pt < PT_N; ++pt) {
      lspB[c][pt].assign(n, {});
      for (int j = 0; j < n; ++j) {
        std::vector<double>& o = lspB[c][pt][j];
        if (pt >= PT_S) {
          if (pt > PT_S && !l2r_inc) continue;
          for (const SSt& s : l2r_smap(pt, j).v) o.push_back(s.bF);
        } else {
          for (const Grp& g : l2r_tab(pt, j).grps)
            for (const Ent& e : g.ents) o.push_back(e.b);
        }
      }
    }
  }
  void l2r_split_run() {
    lsp_rec = false;
    for (int c = 0; c < 3; ++c) {
      l2r_reset_b();
      outside_done = false;
      lsp_use = c;
      outside();
      l2r_collect_b(c);
    }
    lsp_use = -1;
    l2r_reset_b();
    outside_done = false;
    outside();
  }

  // Exact run (beam inf, after the outside pass): rho(a') = log sum_ys Qp Qp_out -
  // log sum_ys Qp e^(beta1) C(a'-1) Sbar(ys+1); the outside of every U and Us state.
  void l2r_oracle_dump(const std::string& path) {
    const StemEst& st = *stem_;
    std::vector<double> num(n, NEG_INF), den(n, NEG_INF), dcx(n, NEG_INF);
    for (int ys = 0; ys < n; ++ys) {
      const OMap& at = Qp_at[ys];
      for (size_t r = 0; r < at.slots(); ++r) {
        const int a = at.keys[r];
        if (a < 0) continue;
        num[a] = log_add(num[a], at.vals[r] + Qp_b[ix(a, ys)]);
        den[a] = log_add(den[a], at.vals[r] + E.w_beta1 + prefix_log(a) + st.prS[ys + 1]);
        dcx[a] = log_add(dcx[a], at.vals[r] + pr_ctx(a, ys));
      }
    }
    std::string ln = "[l2r-rho] a:rho_ext/rho_ctx/log_mass";
    for (int a = 0; a < n; ++a) {
      if (num[a] < HALF_NEG) continue;
      char b[96];
      std::snprintf(b, sizeof b, " %d:%.2f/%.2f/%.1f", a, num[a] - den[a], num[a] - dcx[a],
                    num[a] - logZ);
      ln += b;
    }
    std::fprintf(stderr, "%s\n", ln.c_str());
    // exterior suffix with pseudoknots (C_b[y]) against the pre-pass's nested Sbar[y+1]
    ln = "[l2r-suffix] y:log(C_b[y]/Sbar[y+1])";
    for (int y = 0; y < n; ++y) {
      char b[48];
      std::snprintf(b, sizeof b, " %d:%.2f", y, C_b[y] - st.prS[y + 1]);
      ln += b;
    }
    std::fprintf(stderr, "%s\n", ln.c_str());
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) {
      std::fprintf(stderr, "--l2r-oracle-dump: cannot write %s\n", path.c_str());
      return;
    }
    const int32_t nn = n;
    std::fwrite(&nn, sizeof nn, 1, f);
    for (int a = 0; a < n; ++a) {
      const double r = num[a] > HALF_NEG && den[a] > HALF_NEG ? num[a] - den[a] : 0.0;
      std::fwrite(&r, sizeof r, 1, f);
    }
    for (int seed = 0; seed < 2; ++seed) {
      std::vector<std::pair<uint64_t, double>> rec;
      for (int j = 0; j < n; ++j)
        for (const Grp& g : (seed ? LB.Gs[j] : LB.G[j]).grps)
          for (const Ent& e : g.ents)
            if (e.b > HALF_NEG) rec.emplace_back(key4(e.u, g.v, e.p, j), e.b);
      const int64_t c = static_cast<int64_t>(rec.size());
      std::fwrite(&c, sizeof c, 1, f);
      for (const auto& [k, b] : rec) {
        std::fwrite(&k, sizeof k, 1, f);
        std::fwrite(&b, sizeof b, 1, f);
      }
    }
    std::vector<std::pair<uint64_t, double>> rec;
    for (int j = 0; j < n; ++j)
      for (const Grp& g : LA.Gs[j].grps)
        for (const Ent& e : g.ents)
          if (e.b > HALF_NEG) rec.emplace_back(key4(e.u, g.v, e.p, j), e.b);
    const int64_t c = static_cast<int64_t>(rec.size());
    std::fwrite(&c, sizeof c, 1, f);
    for (const auto& [k, b] : rec) {
      std::fwrite(&k, sizeof k, 1, f);
      std::fwrite(&b, sizeof b, 1, f);
    }
    std::fclose(f);
  }
  bool l2r_oracle_load(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return false;
    int32_t nn = 0;
    bool ok = std::fread(&nn, sizeof nn, 1, f) == 1 && nn == n;
    lRho.assign(n, 0.0);
    if (ok) ok = std::fread(lRho.data(), sizeof(double), n, f) == static_cast<size_t>(n);
    for (int seed = 0; ok && seed < 2; ++seed) {
      int64_t c = 0;
      ok = std::fread(&c, sizeof c, 1, f) == 1;
      std::unordered_map<uint64_t, double>& m = seed ? lOUs : lOU;
      m.reserve(static_cast<size_t>(c));
      for (int64_t i = 0; ok && i < c; ++i) {
        uint64_t k;
        double b;
        ok = std::fread(&k, sizeof k, 1, f) == 1 && std::fread(&b, sizeof b, 1, f) == 1;
        m.emplace(k, b);
      }
    }
    int64_t c = 0;
    if (ok && std::fread(&c, sizeof c, 1, f) == 1) {
      lORs.reserve(static_cast<size_t>(c));
      for (int64_t i = 0; ok && i < c; ++i) {
        uint64_t k;
        double b;
        ok = std::fread(&k, sizeof k, 1, f) == 1 && std::fread(&b, sizeof b, 1, f) == 1;
        lORs.emplace(k, b);
      }
    } else if (lsub_e) {
      ok = false;
    }
    std::fclose(f);
    return ok;
  }
  // --l2r-sub e (oracle diagnostic): A seeds ranked by alpha + their exact outside (exact run's dump)
  bool lsub_e = false;
  std::unordered_map<uint64_t, double> lORs;
  double lors(int u, int v, int p, int j) const {
    auto it = lORs.find(key4(u, v, p, j));
    return it == lORs.end() ? NEG_INF : it->second;
  }
  // log sum_{y' > j} E_in V of arm r
  double lS(int r, int j) const {
    const int hi = prArmY0[r + 1];
    const int k = static_cast<int>(std::upper_bound(prYy.begin() + prArmY0[r], prYy.begin() + hi, j) -
                                   prYy.begin());
    return k < hi ? lYwS[k] : NEG_INF;
  }
  double lsc_L(const RN& r, int arm, int who) {
    const auto t0 = std::chrono::steady_clock::now();
    auto done = [&](double s) {
      const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      if (sx_k_exact) sx_kx_sec += dt;
      else sx_k_sec += dt;
      kl_in_read_sec += dt;
      return s;
    };
    if (arm < 0) return done(NEG_INF);
    double s = NEG_INF;
    auto note = [&](int m, int list) {
      if (sx_k_exact) return;
      ++sx_k_n;
      ++sx_k_src[sx_ksrc];
      if (m >= 0 && m < 8) ++sx_m_h[m];
      if (list >= 0) ++sx_list_h[std::min(list, 64)];
    };
    auto note_arm = [&]() {
      if (arm < 0) return;
      if (arm_seen.empty()) arm_seen.assign(static_cast<size_t>(std::max(prA, 1)), 0);
      if (arm >= prA || arm_seen[static_cast<size_t>(arm)]) return;
      arm_seen[static_cast<size_t>(arm)] = 1;
      ++arm_read_n;
      // |{f : V(c,f)[b] kept}| for arm (b,c)=(xs,xp)
      int nf = 0;
      if (stem_ && !stem_->prMu.empty()) {
        const int b = prArmXs[arm], c = prArmXp[arm];
        for (int f : pr_pair_list(true)[c]) {
          auto it = stem_->prMu.find(key2(c, f));
          if (it == stem_->prMu.end()) continue;
          if (it->second.find(key2(b, 0)) != it->second.end()) ++nf;
        }
      }
      arm_Vf_sum += nf;
      if (nf > arm_Vf_max) arm_Vf_max = nf;
      arm_Vf_samples.push_back(nf);
      // |{q' : K_L(arm,q') nonzero}|
      int n_kq = 0;
      for (int qp = 0; qp + 1 <= n - 1; ++qp)
        if (lK(arm, qp) > HALF_NEG) ++n_kq;
      arm_Kq_sum += n_kq;
      if (n_kq > arm_Kq_max) arm_Kq_max = n_kq;
    };
    if (g_ntr_K > 0) {
      const int len = static_cast<int>(r.T[0].size());
      for (const auto& [qp, g] : r.T[0]) {
        note(-1, len);
        const double k = lK(arm, qp);
        if (k > HALF_NEG) s = log_add(s, g + k);
      }
      note_arm();
      return done(s);
    }
    if (g_tr_K > 0) {
      const auto& L = rec_trunc(r, who).L;
      const int len = static_cast<int>(L.size());
      for (const auto& [qp, g] : L) {
        note(-1, len);
        const double k = lK(arm, qp);
        if (k > HALF_NEG) s = log_add(s, g + k);
      }
      return done(s);
    }
    for (int m = 0; m <= g_ml_L && m < rnL; ++m) {
      const RV& a = r.L[m];
      for (size_t i = 0; i < a.v.size(); ++i) {
        if (a.v[i] <= 0.0) continue;
        note(m, -1);
        sx_m_sum += m;
        const double k = lK(arm, a.lo + static_cast<int>(i));
        if (k > HALF_NEG) s = log_add(s, a.mx + std::log(a.v[i]) + k);
      }
    }
    return done(s);
  }
  double lsc_R(const RN& r, int xs) const {
    if (g_ntr_K > 0 && g_ntr_KR < 0) {
      const RV& a = r.R[0];
      const int k = xs - a.lo;
      return k >= 0 && k < static_cast<int>(a.v.size()) && a.v[k] > 0.0 ? a.mx + std::log(a.v[k]) : NEG_INF;
    }
    if (g_ntr_K > 0) {
      const auto& t = r.T[1];
      auto it = std::lower_bound(t.begin(), t.end(), xs,
                                 [](const std::pair<int, double>& e, int v) { return e.first < v; });
      if (it != t.end() && it->first == xs) return it->second;
      ++ntr_rmiss;
      return r.rest;
    }
    double s = NEG_INF;
    for (int m = 0; m <= g_ml_R && m < rnR; ++m) {
      const RV& a = r.R[m];
      const int k = xs - a.lo;
      if (k >= 0 && k < static_cast<int>(a.v.size()) && a.v[k] > 0.0)
        s = log_add(s, a.mx + std::log(a.v[k]));
    }
    return s;
  }

  // --l2r-diag: no pruning; every state's score kept per (table, step) in entry order, and
  // after the outside pass the posterior mass a top-b cut by score drops against the cut by
  // the posterior itself (oracle).
  bool l2r_diag = false;
  std::vector<std::vector<double>> ldS[PT_N];
  void ldiag_put(int pt, int j, std::vector<double>&& sc) {
    if (ldS[pt].empty()) ldS[pt].resize(n);
    ldS[pt][j] = std::move(sc);
  }

  // --l2r-trace FILE a',xs,x' (0-based): every state of one pseudoknot path at its step (A
  // tables with arm (xs, x'), B tables with (xs, a'), S with (a', xs, x')); scores only where
  // the beam prunes anyway, so the run is unchanged.  After the outside pass: a, b, posterior.
  FILE* ltr_f = nullptr;
  int ltr_ap = -1, ltr_xs = -1, ltr_xp = -1;
  static const char* ltr_name(int pt) {
    static const char* pn[PT_N] = {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us", "Uls", "Ulss",
                                   "Urs", "S", "P1", "P2", "P3"};
    return pn[pt];
  }
  bool ltr_match4(int pt, const Grp& g, const Ent& e) const {
    return pt < PT_U ? (g.v == ltr_xp && e.u == ltr_xs) : (g.v == ltr_xs && e.u == ltr_ap);
  }
  bool ltr_match3(uint64_t key) const {
    return k3a(key) == ltr_ap && k3b(key) == ltr_xs && k3c(key) == ltr_xp;
  }
  void ltr_line(int pt, int j, int p, int u, int v, double a, const std::vector<double>* sc,
                const std::vector<char>* keep, size_t k, size_t tot, int b) {
    std::fprintf(ltr_f, "scan %s j %d p %d u %d v %d a %.6f", ltr_name(pt), j, p, u, v, a);
    if (sc) {
      const double s = (*sc)[k];
      long long rank = 1;
      for (double x : *sc) rank += x > s;
      std::vector<double> c(*sc);
      std::nth_element(c.begin(), c.begin() + (b - 1), c.end(), std::greater<double>());
      std::fprintf(ltr_f, " s %.6f rank %lld tot %zu b %d cut %.6f kept %d\n", s, rank, tot, b,
                   c[b - 1], static_cast<int>((*keep)[k]));
    } else {
      std::fprintf(ltr_f, " tot %zu b %d noprune\n", tot, b);
    }
  }
  void ltr_step4(const Tab4& t, int pt, int j, const std::vector<double>* sc,
                 const std::vector<char>* keep, size_t tot, int b) {
    size_t k = 0;
    for (const Grp& g : t.grps)
      for (const Ent& e : g.ents) {
        if (ltr_match4(pt, g, e)) ltr_line(pt, j, e.p, e.u, g.v, e.a, sc, keep, k, tot, b);
        ++k;
      }
  }
  void ltr_final() {
    if (!ltr_f) return;
    std::fprintf(ltr_f, "arm %d,%d idx %d logZ %.9f\n", ltr_xs, ltr_xp, larm(ltr_xs, ltr_xp), logZ);
    for (int j = 0; j < n; ++j) {
      for (int pt = 0; pt < PT_S; ++pt)
        for (const Grp& g : l2r_tab(pt, j).grps)
          for (const Ent& e : g.ents)
            if (ltr_match4(pt, g, e))
              std::fprintf(ltr_f, "final %s j %d p %d u %d v %d a %.6f b %.6f post %.6e\n",
                           ltr_name(pt), j, e.p, e.u, g.v, e.a, e.b,
                           e.a > HALF_NEG && e.b > HALF_NEG ? std::exp(e.a + e.b - logZ) : 0.0);
      if (j == ltr_c1)
        for (const Grp& g : LA.Gs[j].grps)
          for (const Ent& e : g.ents)
            if (e.p == ltr_a2)
              std::fprintf(ltr_f, "final-seed j %d p %d u %d v %d arm %d a %.6f b %.6f post %.6e\n", j, e.p, e.u,
                           g.v, larm(e.u, g.v), e.a, e.b,
                           e.a > HALF_NEG && e.b > HALF_NEG ? std::exp(e.a + e.b - logZ) : 0.0);
      const SMap& m = LS[j];
      for (size_t k = 0; k < m.keys.size(); ++k)
        if (ltr_match3(m.keys[k]))
          std::fprintf(ltr_f, "final S j %d p %d u %d v %d a %.6f b %.6f post %.6e\n", j,
                       k3a(m.keys[k]), k3b(m.keys[k]), k3c(m.keys[k]), m.v[k].F, m.v[k].bF,
                       m.v[k].F > HALF_NEG && m.v[k].bF > HALF_NEG
                           ? std::exp(m.v[k].F + m.v[k].bF - logZ) : 0.0);
    }
    std::fclose(ltr_f);
    ltr_f = nullptr;
  }

  // --l2r-perf on: spawn-stage counters and wall-clock timers, reported on stderr ([l2r-perf]);
  // nothing they touch feeds a value or a pruning decision.
  bool lpf = false;
  long long lpf_npos[PT_N] = {}, lpf_full[PT_N] = {}, lpf_nsc[PT_N] = {}, lpf_premax[PT_N] = {};
  double lpf_prune_sec[PT_N] = {};
  // filler (l2r_fext) per table 0 S, 1 P1, 2 P3: items d.j, state x item visits, unpaired visits
  long long lpf_fx_pos[3] = {}, lpf_fx_items[3] = {}, lpf_fx_pair[3] = {}, lpf_fx_unp[3] = {},
            lpf_fx_pairmax[3] = {};
  double lpf_fx_sec[3] = {};
  // lp_score per table 0 P1, 1 P2, 2 P3
  long long lpf_sc_calls[3] = {}, lpf_sc_hsteps[3] = {}, lpf_sc_dot[3] = {}, lpf_sc_kcalls[3] = {},
            lpf_sc_ksteps[3] = {}, lpf_sc_osteps[3] = {};
  double lpf_sc_hsec[3] = {}, lpf_sc_ksec[3] = {}, lpf_sc_osec[3] = {};
  // H construction
  long long lpf_h_merge = 0, lpf_h_trunc = 0, lpf_h_trunc_in = 0, lpf_h_newroots = 0, lpf_h_lam = 0,
            lpf_h_helpers = 0, lpf_h_clears = 0;
  double lpf_h_root_sec = 0.0, lpf_h_merge_sec = 0.0, lpf_h_help_sec = 0.0;
  // growth helpers and per-position overheads
  long long lpf_hg_clears = 0, lpf_help_built = 0, lpf_help_rebuilt = 0, lpf_racc_dense = 0,
            lpf_racc_fin = 0, lpf_racc_n = 0;
  double lpf_hg_sec = 0.0, lpf_mtrim_sec = 0.0, lpf_help_sec = 0.0;
  long long lpf_ml_cand = 0, lpf_ml_max = 0, lpf_ml_nodes = 0, lpf_stc_n = 0, lpf_stc_buckets = 0,
            lpf_stc_scan = 0;
  double lpf_stc_sec = 0.0;
  // init (V build): nodes (outer, inner), pr_steps' interior / multiloop steps, multiloop
  // candidates (helper terms) per node, entries added (interior children, helpers, inside the
  // helpers), entries after the finish / after the trim, arm entries, partner pairs
  long long vi_nodes[2] = {}, vi_nb_int = 0, vi_nb_ml = 0, vi_mlc = 0, vi_mlc_max = 0, vi_ml_nodes = 0,
            vi_add_int = 0, vi_add_ml = 0, vi_help_adds = 0, vi_ent_fin = 0, vi_ent_kept = 0, vi_ats = 0,
            vi_pairs = 0;
  double vi_steps_s = 0.0, vi_terms_s = 0.0, vi_add_int_s = 0.0, vi_add_ml_s = 0.0, vi_help_s = 0.0,
         vi_fin_s = 0.0, vi_ats_s = 0.0, vi_vbuild_s = 0.0, vi_arm_s = 0.0, vi_rest_s = 0.0, vi_pairs_s = 0.0;
  static void vi_lap(std::chrono::steady_clock::time_point& t0, double& acc) {
    const auto t1 = std::chrono::steady_clock::now();
    acc += std::chrono::duration<double>(t1 - t0).count();
    t0 = t1;
  }
  std::unordered_set<uint64_t> lpf_help_keys;
  static long long lpf_log2(size_t m) {
    long long s = 0;
    while (m > 0) {
      ++s;
      m >>= 1;
    }
    return s;
  }
  void lpf_prune_stat(int pt, size_t pre, size_t kept, int b, bool scored) {
    ++lpf_npos[pt];
    lpf_premax[pt] = std::max(lpf_premax[pt], static_cast<long long>(pre));
    if (b > 0 && kept >= static_cast<size_t>(b)) ++lpf_full[pt];
    if (scored) lpf_nsc[pt] += static_cast<long long>(pre);
  }

  // Top l2r_beam entries of t by score(group, entry).
  template <class Score>
  void l2r_prune4(Tab4& t, int pt, Score score, int j) {
    std::unique_ptr<SecTimer> ptm(lpf ? new SecTimer(lpf_prune_sec[pt]) : nullptr);
    size_t tot = 0;
    for (const Grp& g : t.grps) tot += g.ents.size();
    const size_t pre = tot;
    l2r_pre[pt] += static_cast<long long>(tot);
    ++l2r_pos[pt];
    if (l2r_diag) {
      SecTimer tm(l2r_psec[pt]);
      std::vector<double> sc;
      sc.reserve(tot);
      for (const Grp& g : t.grps)
        for (const Ent& e : g.ents) sc.push_back(score(g, e));
      ldiag_put(pt, j, std::move(sc));
      l2r_kept[pt] += static_cast<long long>(tot);
      l2r_kmax[pt] = std::max(l2r_kmax[pt], static_cast<long long>(tot));
      return;
    }
    const int b = ltb(pt);
    if (b > 0 && pre > static_cast<size_t>(b)) ++l2r_bind[pt];
    if (b > 0 && tot > static_cast<size_t>(b)) {
      SecTimer tm(l2r_psec[pt]);
      std::vector<double> sc;
      sc.reserve(tot);
      for (const Grp& g : t.grps)
        for (const Ent& e : g.ents) sc.push_back(score(g, e));
      l2r_nsc[pt] += static_cast<long long>(tot);
      for (double s : sc)
        if (s < HALF_NEG) ++l2r_ninf[pt];
      std::vector<char> keep;
      keep_mask(sc, keep, b);
      if (drops) {
        static const char* pn[PT_N] = {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us", "Uls", "Ulss",
                                       "Urs", "S", "P1", "P2", "P3"};
        std::vector<size_t> rk = score_ranks(sc);
        size_t k = 0;
        for (const Grp& g : t.grps)
          for (const Ent& e : g.ents) {
            if (!keep[k])
              note_l2r_drop(pn[pt], j, g.v, e.p, e.u, g.v, j, rk[k], tot, sc[k]);
            ++k;
          }
      }
      if (l2r_ninf_mass || l2r_keep_ninf) {
        size_t k = 0;
        for (const Grp& g : t.grps)
          for (const Ent& e : g.ents) {
            if (sc[k] < HALF_NEG) ninf_note(pt, j, ninf_key(g.v, e.u, e.p), e.a, keep[k]);
            ++k;
          }
      }
      if (g_beam_stats.f) {
        static const char* pn[PT_N] = {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us", "Uls", "Ulss",
                                       "Urs", "S", "P1", "P2", "P3"};
        g_beam_stats.log_prune_mask(pn[pt], j, sc, keep, b);
      }
      if (ltr_f) ltr_step4(t, pt, j, &sc, &keep, tot, b);
      filter_step(t, keep);
      tot = l2r_keep_ninf ? static_cast<size_t>(std::count(keep.begin(), keep.end(), 1)) : static_cast<size_t>(b);
    } else if (ltr_f) {
      ltr_step4(t, pt, j, nullptr, nullptr, tot, b);
    }
    if (lpf) lpf_prune_stat(pt, pre, tot, b, pre != tot || (b > 0 && pre > static_cast<size_t>(b)));
    l2r_kept[pt] += static_cast<long long>(tot);
    l2r_kmax[pt] = std::max(l2r_kmax[pt], static_cast<long long>(tot));
  }
  void l2r_pruneS(int j) {
    l2r_prune_sm(LS[j], PT_S, j, [&](uint64_t key, double f) { return l2r_scoreS(key, f, j); });
  }
  // --l2r-p3-split Ba,bh: P3 kept as the top Ba open a2 values by log sum of their states'
  // scores, then the top bh holes (xs, x') of each kept a2 by score (ties: earliest entry);
  // applied when P3 holds more than its beam (Ba bh).
  int lp_split_a = 0, lp_split_h = 0;
  long long lp_split_grp = 0, lp_split_grp_kept = 0;
  void keep_split(const std::vector<uint64_t>& keys, const std::vector<double>& sc, std::vector<char>& keep) {
    std::unordered_map<int, int> gi;
    std::vector<std::vector<int>> mem;
    std::vector<double> gs;
    for (size_t k = 0; k < keys.size(); ++k) {
      auto [it, ins] = gi.emplace(k3a(keys[k]), static_cast<int>(mem.size()));
      if (ins) {
        mem.emplace_back();
        gs.push_back(NEG_INF);
      }
      mem[it->second].push_back(static_cast<int>(k));
      if (sc[k] > HALF_NEG) gs[it->second] = log_add(gs[it->second], sc[k]);
    }
    keep.assign(keys.size(), 0);
    const int ng = static_cast<int>(mem.size());
    lp_split_grp += ng;
    std::vector<char> gk;
    if (ng > lp_split_a) keep_mask(gs, gk, lp_split_a);
    else gk.assign(ng, 1);
    std::vector<double> ms;
    std::vector<char> mk;
    for (int g = 0; g < ng; ++g) {
      if (!gk[g]) continue;
      ++lp_split_grp_kept;
      const std::vector<int>& v = mem[g];
      if (static_cast<int>(v.size()) <= lp_split_h) {
        for (int k : v) keep[k] = 1;
        continue;
      }
      ms.clear();
      for (int k : v) ms.push_back(sc[k]);
      keep_mask(ms, mk, lp_split_h);
      for (size_t r = 0; r < v.size(); ++r)
        if (mk[r]) keep[v[r]] = 1;
    }
  }
  void l2r_pruneP(int pt, int j) {
    std::vector<SMap>& T = pt == PT_P1 ? LP1 : pt == PT_P2 ? LP2 : LP3;
    l2r_prune_sm(T[j], pt, j, [&](uint64_t key, double f) { return lp_score(pt, key, f, j); });
  }
  template <class Score>
  void l2r_prune_sm(SMap& m, int pt, int j, Score score) {
    std::unique_ptr<SecTimer> ptm(lpf ? new SecTimer(lpf_prune_sec[pt]) : nullptr);
    size_t tot = m.keys.size();
    l2r_pre[pt] += static_cast<long long>(tot);
    ++l2r_pos[pt];
    const int b = l2r_diag ? 0 : ltb(pt);
    const size_t pre = tot;
    if (b > 0 && pre > static_cast<size_t>(b)) ++l2r_bind[pt];
    if (l2r_diag || (b > 0 && tot > static_cast<size_t>(b))) {
      SecTimer tm(l2r_psec[pt]);
      std::vector<double> sc(tot);
      for (size_t k = 0; k < tot; ++k) {
        const uint64_t key = m.keys[k];
        sc[k] = score(key, m.v[k].F);
      }
      if (l2r_diag) {
        ldiag_put(pt, j, std::move(sc));
        l2r_kept[pt] += static_cast<long long>(tot);
        l2r_kmax[pt] = std::max(l2r_kmax[pt], static_cast<long long>(tot));
        return;
      }
      l2r_nsc[pt] += static_cast<long long>(tot);
      for (double s : sc)
        if (s < HALF_NEG) ++l2r_ninf[pt];
      std::vector<char> keep;
      if (pt == PT_P3 && lp_split_a > 0) keep_split(m.keys, sc, keep);
      else keep_mask(sc, keep, b);
      if (drops) {
        static const char* pn[PT_N] = {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us", "Uls", "Ulss",
                                       "Urs", "S", "P1", "P2", "P3"};
        std::vector<size_t> rk = score_ranks(sc);
        for (size_t k = 0; k < tot; ++k) {
          if (keep[k]) continue;
          const uint64_t key = m.keys[k];
          note_l2r_drop(pn[pt], j, -1, k3a(key), k3b(key), k3c(key), j, rk[k], tot, sc[k]);
        }
      }
      if (l2r_ninf_mass || l2r_keep_ninf)
        for (size_t k = 0; k < tot; ++k)
          if (sc[k] < HALF_NEG) ninf_note(pt, j, m.keys[k], m.v[k].F, keep[k]);
      if (g_beam_stats.f) {
        static const char* pn[PT_N] = {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us", "Uls", "Ulss",
                                       "Urs", "S", "P1", "P2", "P3"};
        g_beam_stats.log_prune_mask(pn[pt], j, sc, keep, b);
      }
      if (ltr_f && pt == PT_S)
        for (size_t k = 0; k < tot; ++k)
          if (ltr_match3(m.keys[k]))
            ltr_line(pt, j, k3a(m.keys[k]), k3b(m.keys[k]), k3c(m.keys[k]), m.v[k].F, &sc, &keep, k,
                     tot, b);
      SMap o;
      for (size_t k = 0; k < tot; ++k)
        if (keep[k]) {
          o.idx.emplace(m.keys[k], static_cast<int>(o.keys.size()));
          o.keys.push_back(m.keys[k]);
          o.v.push_back(m.v[k]);
        }
      m = std::move(o);
      tot = m.keys.size();
    } else if (ltr_f && pt == PT_S) {
      for (size_t k = 0; k < tot; ++k)
        if (ltr_match3(m.keys[k]))
          ltr_line(pt, j, k3a(m.keys[k]), k3b(m.keys[k]), k3c(m.keys[k]), m.v[k].F, nullptr, nullptr,
                   k, tot, b);
    }
    if (lpf) lpf_prune_stat(pt, pre, tot, b, b > 0 && pre > static_cast<size_t>(b));
    l2r_kept[pt] += static_cast<long long>(tot);
    l2r_kmax[pt] = std::max(l2r_kmax[pt], static_cast<long long>(tot));
  }
  // R-type tables at j: kind 0 root (seed: inner), 1 ls key, 2 rs key
  struct ProfReg {
    DP03& d;
    int prev;
    bool on;
    ProfReg(DP03& d, int r) : d(d), prev(d.prof_region), on(d.prof_on) {
      if (on) d.prof_region = r;
    }
    ~ProfReg() {
      if (on) d.prof_region = prev;
    }
  };
  // Rank one R column by the exact gap sum. Every candidate gets sigma_ub from lK.
  // Exact sigma (lKx, memoized) is computed in that order and stops once the next
  // upper bound is below the current b-th exact score. exact <= ub carries the
  // per-state additives, so the comparison is on the full score.
  void l2r_pruneA_exact(Tab4& t, int pt, int kind, bool seed, int j) {
    ProfReg reg(*this, 1);
    std::unique_ptr<SecTimer> ptm(lpf ? new SecTimer(lpf_prune_sec[pt]) : nullptr);
    size_t tot = 0;
    for (const Grp& g : t.grps) tot += g.ents.size();
    const size_t pre = tot;
    l2r_pre[pt] += static_cast<long long>(tot);
    ++l2r_pos[pt];
    const int b = ltb(pt);
    if (b > 0 && pre > static_cast<size_t>(b)) ++l2r_bind[pt];
    if (!(b > 0 && tot > static_cast<size_t>(b))) {
      if (ltr_f) ltr_step4(t, pt, j, nullptr, nullptr, tot, b);
      if (lpf) lpf_prune_stat(pt, pre, tot, b, false);
      l2r_kept[pt] += static_cast<long long>(tot);
      l2r_kmax[pt] = std::max(l2r_kmax[pt], static_cast<long long>(tot));
      return;
    }
    SecTimer tm(l2r_psec[pt]);
    struct Ref {
      const Grp* g;
      const Ent* e;
    };
    std::vector<Ref> ref;
    ref.reserve(tot);
    for (const Grp& g : t.grps)
      for (const Ent& e : g.ents) ref.push_back(Ref{&g, &e});
    auto score = [&](const Ref& c, bool exact) {
      if (kind != 0 && !exact) ++sx_rkey_n;
      const RN& r =
          kind == 0 ? rec_root(c.e->p, j, seed) : rec_key(kind == 1, c.e->p, kind == 1 && seed, j);
      sx_ksrc = kind == 0 ? 1 : 2;
      sx_k_exact = exact;
      const int arm = larm(c.e->u, c.g->v);
      if (!exact && arm < 0) ++l2r_noarm[pt];
      const double s = lsc_L(r, arm, kind == 0 ? 0 : 1);
      sx_k_exact = false;
      sx_ksrc = 0;
      if (!exact && g_mass_f && s > HALF_NEG)
        g_mass_note(kind, seed ? 1 : 0, c.e->p, j, s);
      return s > HALF_NEG ? c.e->a + s + lkb_arm(c.e->p, arm) : NEG_INF;
    };
    std::vector<double> ub(tot), ex(tot, NEG_INF);
    std::vector<char> got(tot, 0);
    for (size_t i = 0; i < tot; ++i) ub[i] = score(ref[i], false);
    l2r_nsc[pt] += static_cast<long long>(tot);
    for (double s : ub)
      if (s < HALF_NEG) ++l2r_ninf[pt];
    std::vector<int> ord(tot);
    for (size_t i = 0; i < tot; ++i) ord[i] = static_cast<int>(i);
    std::sort(ord.begin(), ord.end(), [&](int a, int c) { return ub[a] > ub[c]; });
    std::priority_queue<double, std::vector<double>, std::greater<double>> heap;
    for (int idx : ord) {
      if (static_cast<int>(heap.size()) >= b && ub[idx] < heap.top()) break;
      ex[idx] = score(ref[idx], true);
      got[idx] = 1;
      ++sx_kx_sigma;
      if (static_cast<int>(heap.size()) < b) heap.push(ex[idx]);
      else if (ex[idx] > heap.top()) {
        heap.pop();
        heap.push(ex[idx]);
      }
    }
    std::vector<char> keep;
    keep_mask(ex, keep, b);
    if (k_exact_check) {
      ++k_exact_cols;
      std::vector<double> full(tot);
      for (size_t i = 0; i < tot; ++i) full[i] = got[i] ? ex[i] : score(ref[i], true);
      std::vector<char> keep2;
      keep_mask(full, keep2, b);
      if (keep != keep2) ++k_exact_mis;
    }
    if (l2r_ninf_mass || l2r_keep_ninf) {
      for (size_t k = 0; k < tot; ++k)
        if (ex[k] < HALF_NEG)
          ninf_note(pt, j, ninf_key(ref[k].g->v, ref[k].e->u, ref[k].e->p), ref[k].e->a, keep[k]);
    }
    if (ltr_f) ltr_step4(t, pt, j, &ex, &keep, tot, b);
    filter_step(t, keep);
    tot = l2r_keep_ninf ? static_cast<size_t>(std::count(keep.begin(), keep.end(), 1))
                        : static_cast<size_t>(b);
    if (lpf) lpf_prune_stat(pt, pre, tot, b, true);
    l2r_kept[pt] += static_cast<long long>(tot);
    l2r_kmax[pt] = std::max(l2r_kmax[pt], static_cast<long long>(tot));
  }
  void l2r_pruneA(Tab4& t, int pt, int kind, bool seed, int j) {
    if (k_exact_ub && !l2r_diag && !l2r_arm_fb && !(lsub_e && kind == 0 && seed)) {
      l2r_pruneA_exact(t, pt, kind, seed, j);
      return;
    }
    ProfReg reg(*this, 1);
    l2r_prune4(t, pt, [&](const Grp& g, const Ent& e) {
      if (lsub_e && kind == 0 && seed) {
        const double b = lors(e.u, g.v, e.p, j);
        return b > HALF_NEG ? e.a + b : NEG_INF;
      }
      if (kind != 0) ++sx_rkey_n;
      const RN& r = kind == 0 ? rec_root(e.p, j, seed) : rec_key(kind == 1, e.p, kind == 1 && seed, j);
      sx_ksrc = kind == 0 ? 1 : 2;
      if (l2r_arm_fb) {
        const auto [arm, fx] = larm_fb(e.u, g.v);
        if (arm >= 0 && arm != larm(e.u, g.v)) ++l2r_fb_hit;
        const double s = lsc_L(r, arm, kind == 0 ? 0 : 1);
        sx_ksrc = 0;
        if (s < HALF_NEG) {
          ++l2r_fb_floor;
          return e.a + kFbFloor;
        }
        return e.a + (s + fx) + lkb_arm(e.p, arm);
      }
      const int arm = larm(e.u, g.v);
      if (arm < 0) ++l2r_noarm[pt];
      const double s = lsc_L(r, arm, kind == 0 ? 0 : 1);
      sx_ksrc = 0;
      if (g_mass_f && s > HALF_NEG) g_mass_note(kind, seed ? 1 : 0, e.p, j, s);
      return s > HALF_NEG ? e.a + s + lkb_arm(e.p, arm) : NEG_INF;
    }, j);
  }
  void l2r_pruneB(Tab4& t, int pt, int kind, bool seed, int j) {
    l2r_prune4(t, pt, [&](const Grp& g, const Ent& e) {
      if (lsub_d && kind == 0) {
        const auto& m = seed ? lOUs : lOU;
        auto it = m.find(key4(e.u, g.v, e.p, j));
        return it == m.end() ? NEG_INF : e.a + it->second - (lsub_c ? 0.0 : lrho(e.u));
      }
      double s;
      if (g_rfv) {
        s = kind == 0 ? rfv_root(e.p, j, seed, g.v) : rfv_key(kind == 1, e.p, kind == 1 && seed, j, g.v);
      } else {
        const RN& r = kind == 0 ? rec_root(e.p, j, seed) : rec_key(kind == 1, e.p, kind == 1 && seed, j);
        s = lsc_R(r, g.v);
      }
      if (!lLam.empty())
        return s > HALF_NEG ? e.a + s + prefix_log(e.u) + (lsub_c ? lrho(e.u) : lkb_xs(e.u, g.v)) +
                                  lLam[e.u]
                            : NEG_INF;
      return s > HALF_NEG ? e.a + s + prefix_log(e.u) + (lsub_c ? lrho(e.u) : lkb_xs(e.u, g.v))
                          : NEG_INF;
    }, j);
  }

  // lost[t][0..1]: posterior mass outside the top-b by score / by posterior, b = bs[0], bs[1]
  std::string l2r_diag_report(const std::vector<int>& bs) {
    static const char* pn[PT_N] = {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us", "Uls", "Ulss",
                                   "Urs", "S", "P1", "P2", "P3"};
    std::string out;
    char buf[512];
    for (int pt = 0; pt < PT_N; ++pt) {
      if (ldS[pt].empty()) continue;
      std::vector<double> lostS(bs.size(), 0.0), lostO(bs.size(), 0.0);
      double tot = 0.0, minf = 0.0, dsum = 0.0, dsq = 0.0, dw = 0.0, hdrop[8] = {};
      long long steps_over = 0, ninf = 0, hin[8] = {};
      std::vector<std::pair<double, double>> dq;  // (d, posterior weight) for the quantiles
      std::vector<std::array<double, 5>> ds;      // (d, weight, exterior / ml / gap share of b)
      for (int j = 0; j < n; ++j) {
        const std::vector<double>& sc = ldS[pt][j];
        std::vector<double> po, pbv;
        if (pt >= PT_S) {
          const SMap& m = l2r_smap(pt, j);
          for (size_t k = 0; k < m.keys.size(); ++k) {
            po.push_back(m.v[k].F + m.v[k].bF);
            pbv.push_back(m.v[k].bF);
          }
        } else {
          for (const Grp& g : l2r_tab(pt, j).grps)
            for (const Ent& e : g.ents) {
              po.push_back(e.a + e.b);
              pbv.push_back(e.b);
            }
        }
        if (po.size() != sc.size()) {
          std::snprintf(buf, sizeof buf, " [%s size mismatch at %d]", pn[pt], j);
          out += buf;
          continue;
        }
        std::vector<double> pe(po.size());
        for (size_t k = 0; k < po.size(); ++k) {
          pe[k] = po[k] > HALF_NEG ? std::exp(std::min(0.0, po[k] - logZ)) : 0.0;
          tot += pe[k];
          if (sc[k] < HALF_NEG) {
            ++ninf;
            minf += pe[k];
          } else if (po[k] > HALF_NEG) {
            // score - log posterior weight (the score's outside estimate / the true one)
            const double d = sc[k] - po[k];
            dsum += pe[k] * d;
            dsq += pe[k] * d * d;
            dw += pe[k];
            if (pe[k] > 0.0) {
              dq.emplace_back(d, pe[k]);
              std::array<double, 5> t{d, pe[k], 0.0, 0.0, 0.0};
              const double b = pbv[k];
              for (int c = 0; c < 3 && !lspB[c][pt].empty(); ++c) {
                const double bc = lspB[c][pt][j][k];
                t[2 + c] = bc > HALF_NEG ? std::exp(bc - b) : 0.0;
              }
              ds.push_back(t);
            }
          }
        }
        for (size_t c = 0; c < bs.size(); ++c) {
          const size_t b = static_cast<size_t>(bs[c]);
          if (po.size() <= b) continue;
          if (c == 0) ++steps_over;
          std::vector<char> keep, ko;
          keep_mask(sc, keep, static_cast<int>(b));
          for (size_t k = 0; k < po.size(); ++k)
            if (!keep[k]) lostS[c] += pe[k];
          keep_mask(po, ko, static_cast<int>(b));
          for (size_t k = 0; k < po.size(); ++k)
            if (!ko[k]) lostO[c] += pe[k];
          if (bs[c] != 1000) continue;
          // dropped by score, kept by the oracle: mass by d = score - log post; intruders
          // (kept by score, not by the oracle): count by d
          for (size_t k = 0; k < po.size(); ++k) {
            const double d = sc[k] < HALF_NEG ? -1e9 : po[k] < HALF_NEG ? 1e9 : sc[k] - po[k];
            const int bin = d < -4 ? 0 : d < -2 ? 1 : d < -1 ? 2 : d < 0 ? 3 : d < 1 ? 4 : d < 2 ? 5 : d < 4 ? 6 : 7;
            if (!keep[k] && ko[k]) hdrop[bin] += pe[k];
            if (keep[k] && !ko[k]) ++hin[bin];
          }
        }
      }
      const double mu = dw > 0 ? dsum / dw : 0.0;
      std::snprintf(buf, sizeof buf,
                    " %s mass %.4g steps_over %lld ninf %lld inf_mass %.4g logratio_mean %.3f sd %.3f",
                    pn[pt], tot, steps_over, ninf, minf, mu,
                    dw > 0 ? std::sqrt(std::max(0.0, dsq / dw - mu * mu)) : 0.0);
      out += buf;
      // low mode: peak of the weighted d histogram (bins 0.05) below -3; mean, sd, share of the
      // mass within 0.3 of it; consumer split of b in the mode and at |d| < 0.3
      {
        std::vector<double> h(340, 0.0);
        for (const auto& t : ds)
          if (t[0] >= -20.0 && t[0] < -3.0) h[static_cast<size_t>((t[0] + 20.0) / 0.05)] += t[1];
        const size_t pk = static_cast<size_t>(std::max_element(h.begin(), h.end()) - h.begin());
        const double c0 = -20.0 + 0.05 * (static_cast<double>(pk) + 0.5);
        double w = 0.0, m1 = 0.0, m2 = 0.0, sp[3] = {}, w0 = 0.0, s0[3] = {};
        for (const auto& t : ds) {
          if (std::fabs(t[0] - c0) < 0.3) {
            w += t[1];
            m1 += t[1] * t[0];
            m2 += t[1] * t[0] * t[0];
            for (int c = 0; c < 3; ++c) sp[c] += t[1] * t[2 + c];
          }
          if (std::fabs(t[0]) < 0.3) {
            w0 += t[1];
            for (int c = 0; c < 3; ++c) s0[c] += t[1] * t[2 + c];
          }
        }
        const double mu1 = w > 0 ? m1 / w : 0.0;
        std::snprintf(buf, sizeof buf, " mode %.3f sd %.3f share %.3f", mu1,
                      w > 0 ? std::sqrt(std::max(0.0, m2 / w - mu1 * mu1)) : 0.0, dw > 0 ? w / dw : 0.0);
        out += buf;
        if (!lspB[0][pt].empty()) {
          auto f = [](double a, double b) { return b > 0 ? a / b : 0.0; };
          std::snprintf(buf, sizeof buf, " split_mode(ext/ml/gap) %.4f/%.4f/%.4f split_d0 %.4f/%.4f/%.4f",
                        f(sp[0], w), f(sp[1], w), f(sp[2], w), f(s0[0], w0), f(s0[1], w0), f(s0[2], w0));
          out += buf;
        }
      }
      // posterior-weighted quantiles of d = log(score / (inside x exact outside))
      std::sort(dq.begin(), dq.end());
      out += " q01/05/25/50/75/95/99";
      {
        const double qs[7] = {0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99};
        double c = 0.0;
        size_t k = 0;
        for (double q : qs) {
          while (k < dq.size() && c + dq[k].second < q * dw) c += dq[k++].second;
          std::snprintf(buf, sizeof buf, " %.2f", k < dq.size() ? dq[k].first : 0.0);
          out += buf;
        }
      }
      for (size_t c = 0; c < bs.size(); ++c) {
        std::snprintf(buf, sizeof buf, " b%d score %.4g oracle %.4g", bs[c], lostS[c], lostO[c]);
        out += buf;
      }
      out += " b1000 drop_mass_by_d(<-4,-4,-2,-1,0,1,2,>=4)";
      for (int k = 0; k < 8; ++k) {
        std::snprintf(buf, sizeof buf, " %.2g", hdrop[k]);
        out += buf;
      }
      out += " intruders_by_d";
      for (int k = 0; k < 8; ++k) {
        std::snprintf(buf, sizeof buf, " %lld", hin[k]);
        out += buf;
      }
      out += " |";
    }
    return out;
  }

  // Pseudoknot part of step j (after the engine's Qg .. Qp in check mode).
  // Qhat^Lead(x0, t): empty at x0 == t, else w_cp per base plus the dangle of the whole run.
  // Beam b2 per column, key C(x0-1) + Lead. A dropped start is NEG_INF.
  void l2r_lead(int t) {
    LeadCol col;
    col.w.assign(static_cast<size_t>(t) + 1, NEG_INF);
    struct Cand {
      int x0;
      double w, k;
    };
    std::vector<Cand> cs;
    cs.reserve(static_cast<size_t>(t) + 1);
    for (int x0 = 0; x0 <= t; ++x0) {
      if (!lead_gap_ok(t - x0)) continue;
      double w = 0.0;
      if (x0 < t) {
        w = E.w_cp * (t - x0);
        if (dz) w += E.run(x0 + 1, t);
      }
      cs.push_back(Cand{x0, w, prefix_log(x0) + w});
    }
    int b = tbeam[BT_QM];
    if (b > 0 && static_cast<int>(cs.size()) > b) {
      std::nth_element(cs.begin(), cs.begin() + b, cs.end(),
                       [](const Cand& a, const Cand& c) { return a.k > c.k; });
      cs.resize(static_cast<size_t>(b));
    }
    for (const Cand& c : cs) col.w[static_cast<size_t>(c.x0)] = c.w;
    gs_lead += static_cast<long long>(cs.size());
    gsLead[static_cast<size_t>(t)] = std::move(col);
  }

  // M1 and M4 from the stems kept at j. M3 records each M4 state as an M as well.
  // M4's group key is the carried v. Unpaired bases inside Qhat^M are already in that table.
  void l2r_simple_m4(LFam& F, int j) {
    auto from = [&](const Tab4& tab) {
      for (const Grp& g : tab.grps) {
        const int lo = F.isU ? g.v + 1 : 1;
        for (const Ent& s : g.ents) {
          if (E.end_ok(s.p, j)) {
            F.SM[j].acc(g.v, s.p, s.u, s.a + E.w_bp + E.at(s.p, j));
            ++gs_m[0];
          }
          const int c = s.p;
          if (c - 1 < 0 || !E.end_ok(c, j)) continue;
          lchk(LT_QM, c - 1);
          const double w = E.w_bp + E.at(c, j);
          const OMap& at = qmp_at(c - 1);
          for (size_t r = 0; r < at.slots(); ++r) {
            const int i = at.keys[r];
            if (i < lo) continue;
            F.SM2[j].acc(g.v, i, s.u, s.a + w + at.vals[r]);
            ++gs_m[3];
          }
        }
      }
    };
    from(F.G[j]);
    from(F.Gs[j]);
    for (const Grp& g : F.SM2[j].grps) {
      if (g.v >= n) continue;
      for (const Ent& s : g.ents) {
        F.SM[j].acc(g.v, s.p, s.u, s.a);
        ++gs_m[2];
      }
    }
  }

  // M5: stem ended at f = p-1, nested segment p..j on its right. Group key is v + n.
  void l2r_simple_m5(LFam& F, int j) {
    lchk(LT_QM, j);
    const OMap& at = qmp_at(j);
    const int stem = F.isU ? LT_BG : LT_AG;
    for (size_t r = 0; r < at.slots(); ++r) {
      const int p = at.keys[r];
      if (p < 1) continue;
      const int f = p - 1;
      const double qm = at.vals[r];
      lchk(stem, f);
      auto take = [&](const Tab4& tab) {
        for (const Grp& g : tab.grps)
          for (const Ent& s : g.ents) {
            if (!E.end_ok(s.p, f)) continue;
            if (!F.isU && !l2r_arm(s.u, g.v, j)) continue;
            F.SM2[j].acc(g.v + n, s.p, s.u, s.a + E.w_bp + E.at(s.p, f) + qm);
            ++gs_m[4];
          }
      };
      take(F.G[f]);
      take(F.Gs[f]);
    }
  }

  // Closes into X at j. (f) is M6+M7+M8 from M4. (g) is M6+M8 from M5 through Lead.
  // (h) is M5 plus a nested segment on the left. Survivors of the Multi beam are written.
  void l2r_simple_close(LFam& F, int j) {
    Tab4& Mu = F.SMulti[j];
    const int tg = F.isU ? LT_BG : LT_AG;
    for (int f = 0; f < j; ++f) {
      const int unpaired = j - f - 1;
      double w_ext = E.w_ap + E.w_bp + E.w_cp * unpaired + pair_bonus;
      if (dz) w_ext += E.run(f + 1, j - 1);
      for (const Grp& g : F.SM2[f].grps) {
        if (g.v >= n) continue;
        const int lo = F.isU ? g.v : 0;
        for (const Ent& s : g.ents) {
          const int i = s.p - 1;
          if (i < lo || !E.pairs(i, j) || !E.end_ok(i, j)) continue;
          if (!F.isU && !l2r_arm(s.u, g.v, j)) continue;
          Mu.acc(g.v, i, s.u, s.a + w_ext + E.at(i, j));
          ++gs_m[5];
          gs_m[6] += unpaired;
          ++gs_m[7];
          if (unpaired > cap) ++gs_over;
        }
      }
    }
    if (j >= 1) {
      for (const Grp& g : F.SM2[j - 1].grps) {
        if (g.v < n) continue;
        const int real_v = g.v - n;
        const int lo = F.isU ? real_v : 0;
        for (const Ent& s : g.ents) {
          const int c = s.p;
          if (!F.isU && !l2r_arm(s.u, real_v, j)) continue;
          // Leading run is the same window as pass 1: p in [c - ml_lead_cap, c).
          // ml_lead_cap 0 keeps every p. No Qhat^Lead table and no beam on the run.
          {
            const int i_lo = std::max(lo, lead_i_lo(c));
            for (int i = i_lo; i < c; ++i) {
              if (!E.pairs(i, j) || !E.end_ok(i, j)) continue;
              const int gap = c - i - 1;
              double lw = 0.0;
              if (gap > 0) {
                lw = E.w_cp * gap;
                if (dz) lw += E.run(i + 1, c - 1);
              }
              Mu.acc(real_v, i, s.u, s.a + lw + E.w_ap + E.w_bp + E.at(i, j) + pair_bonus);
              ++gs_m[5];
              ++gs_m[7];
              ++gs_lead;
              if (gap > cap) ++gs_over;
            }
          }
          if (c - 1 < 0) continue;
          lchk(LT_QM, c - 1);
          const OMap& at = qmp_at(c - 1);
          qm_pair_slots(c - 1, j, [&](int r) {
            const int ii = at.keys[r];
            if (ii < 0) return;
            const int outer = ii - 1;
            if (outer < lo || outer >= c || !E.pairs(outer, j) || !E.end_ok(outer, j)) return;
            Mu.acc(real_v, outer, s.u,
                   at.vals[r] + s.a + E.w_ap + E.w_bp + E.at(outer, j) + pair_bonus);
            ++gs_m[7];
          });
        }
      }
    }
    l2r_prune_prefix(Mu, l2r_beam, F.isU ? 5 : 2);
    for (const Grp& g : Mu.grps)
      for (const Ent& s : g.ents) lacc(F.G[j], tg, g.v, s.p, s.u, s.a);
  }

  void l2r_simple_pk(int j) {
    if (!l2r_inc) {
      std::fprintf(stderr, "[grammar] simple needs --l2r-spawn inc\n");
      std::abort();
    }
    if (l2r_beam > 0) stem_cache_step(j);
    gs_seedA.clear();
    gs_seedB.clear();
    {
      SecTimer t(l2r_bsec[PT_R]);
      l2r_grow(LA, j);
    }
    {
      SecTimer t(l2r_bsec[PT_RS]);
      l2r_inc_seed(j);
    }
    l2r_pruneA(LA.Gs[j], PT_RS, 0, true, j);
    l2r_simple_close(LA, j);
    l2r_pruneA(LA.G[j], PT_R, 0, false, j);
    l2r_done[LT_AG] = l2r_done[LT_AGS] = j;
    l2r_simple_m4(LA, j);
    {
      SecTimer t(l2r_csec);
      l2r_closeA(j);
    }
    l2r_done[LT_ALS] = j;
    {
      SecTimer t(l2r_bsec[PT_US]);
      l2r_bseed(j);
    }
    l2r_pruneB(LB.Gs[j], PT_US, 0, true, j);
    {
      SecTimer t(l2r_bsec[PT_U]);
      l2r_grow(LB, j);
    }
    l2r_simple_close(LB, j);
    l2r_pruneB(LB.G[j], PT_U, 0, false, j);
    l2r_done[LT_BG] = l2r_done[LT_BGS] = j;
    l2r_simple_m4(LB, j);
    l2r_done[LT_BLS] = j;
    {
      SecTimer t(l2r_csec);
      l2r_closeB(j);
    }
  }

  void l2r_simple_end(int j) {
    {
      SecTimer t(l2r_bsec[PT_S]);
      l2r_sext(j);
    }
    l2r_pruneS(j);
    l2r_done[LT_S] = j;
    {
      SecTimer t(l2r_bsec[PT_P1]);
      l2r_p1_step(j);
    }
    l2r_pruneP(PT_P1, j);
    l2r_done[LT_P1] = j;
    {
      SecTimer t(l2r_bsec[PT_P2]);
      l2r_p2_step(j);
    }
    l2r_pruneP(PT_P2, j);
    l2r_done[LT_P2] = j;
    {
      SecTimer t(l2r_bsec[PT_P3]);
      l2r_p3_step(j);
    }
    l2r_pruneP(PT_P3, j);
    l2r_done[LT_P3] = j;
    l2r_simple_m5(LA, j);
    l2r_simple_m5(LB, j);
    l2r_prune_prefix(LA.SM[j], l2r_beam, 0);
    l2r_prune_prefix(LB.SM[j], l2r_beam, 3);
    l2r_prune_prefix(LA.SM2[j], l2r_beam, 1);
    l2r_prune_prefix(LB.SM2[j], l2r_beam, 4);
    l2r_done[LT_ARS] = l2r_done[LT_BRS] = j;
  }

  void l2r_step_pk(int j) {
    if (grammar_simple) {
      l2r_simple_pk(j);
      return;
    }
    if (l2r_beam > 0) stem_cache_step(j);
    {
      SecTimer t(l2r_bsec[PT_R]);
      l2r_grow(LA, j);
    }
    {
      SecTimer t(l2r_bsec[PT_RS]);
      if (l2r_inc) l2r_inc_seed(j);
      else l2r_spawn(j);
    }
    if (l2r_inc_check) {
      SecTimer t(l2r_sec[5]);
      l2r_inc_check_step(j);
    }
    l2r_pruneA(LA.G[j], PT_R, 0, false, j);
    l2r_pruneA(LA.Gs[j], PT_RS, 0, true, j);
    l2r_done[LT_AG] = l2r_done[LT_AGS] = j;
    if (l2r_mode == 2 && l2r_check_R) {
      SecTimer t(l2r_sec[5]);
      l2r_checkR(j);
    }
    {
      SecTimer t(l2r_bsec[PT_RLS]);
      l2r_ls(LA, j);
    }
    l2r_pruneA(LA.Gls[j], PT_RLS, 1, false, j);
    l2r_pruneA(LA.Glss[j], PT_RLSS, 1, true, j);
    l2r_done[LT_ALS] = j;
    {
      SecTimer t(l2r_csec);
      l2r_closeA(j);
    }
    {
      SecTimer t(l2r_bsec[PT_US]);
      l2r_bseed(j);
    }
    {
      SecTimer t(l2r_bsec[PT_U]);
      l2r_grow(LB, j);
    }
    l2r_pruneB(LB.G[j], PT_U, 0, false, j);
    l2r_pruneB(LB.Gs[j], PT_US, 0, true, j);
    l2r_done[LT_BG] = l2r_done[LT_BGS] = j;
    {
      SecTimer t(l2r_bsec[PT_ULS]);
      l2r_ls(LB, j);
    }
    l2r_pruneB(LB.Gls[j], PT_ULS, 1, false, j);
    l2r_pruneB(LB.Glss[j], PT_ULSS, 1, true, j);
    l2r_done[LT_BLS] = j;
    {
      SecTimer t(l2r_csec);
      l2r_closeB(j);
    }
    if (rlt) {
      l2r_tail(LA, j);
      l2r_tail(LB, j);
    }
  }

  // End of step j: S (needs Qb, Qp, Qz at j), then Grs (needs Qm(., j)).
  void l2r_step_end(int j) {
    if (grammar_simple) {
      l2r_simple_end(j);
      return;
    }
    {
      SecTimer t(l2r_bsec[PT_S]);
      l2r_sext(j);
    }
    l2r_pruneS(j);
    l2r_done[LT_S] = j;
    if (l2r_inc) {
      {
        SecTimer t(l2r_bsec[PT_P1]);
        l2r_p1_step(j);
      }
      l2r_pruneP(PT_P1, j);
      l2r_done[LT_P1] = j;
      {
        SecTimer t(l2r_bsec[PT_P2]);
        l2r_p2_step(j);
      }
      l2r_pruneP(PT_P2, j);
      l2r_done[LT_P2] = j;
      {
        SecTimer t(l2r_bsec[PT_P3]);
        l2r_p3_step(j);
      }
      l2r_pruneP(PT_P3, j);
      l2r_done[LT_P3] = j;
    }
    {
      SecTimer t(l2r_bsec[PT_RRS]);
      l2r_rs(LA, j);
    }
    l2r_pruneA(LA.Grs[j], PT_RRS, 2, false, j);
    {
      SecTimer t(l2r_bsec[PT_URS]);
      l2r_rs(LB, j);
    }
    l2r_pruneB(LB.Grs[j], PT_URS, 2, false, j);
    l2r_done[LT_ARS] = l2r_done[LT_BRS] = j;
  }

  // ---- outside (--pk-l2r on): each step mirrors the forward one in reverse order.
  void l2r_sext_back(int j) {
    l2r_fext_back(LS, j, [&](uint64_t key) { return l2r_arm(k3b(key), k3c(key), j); });
  }
  // Adjoint of l2r_fext at j: bF of T[j] (complete) into bP / bU1 / bU2 / bc, then into the
  // items d.j and T[j-1], T[d-1].
  template <class Alive>
  void l2r_fext_back(std::vector<SMap>& T, int j, Alive alive) {
    SMap& cur = T[j];
    for (SSt& s : cur.v) {
      if (s.bF < HALF_NEG) continue;
      s.bP = log_add(s.bP, s.bF);
      s.bc = log_add(s.bc, s.bF);
      if (dz) {
        s.bU1 = log_add(s.bU1, s.bF + E.s1w[j]);
        s.bU2 = log_add(s.bU2, s.bF + E.d3w[j]);
      }
    }
    auto branch = [&](int d, double qv, double bw, SpMat& qb) {
      if (d < 1) return;
      SMap& pv = T[d - 1];
      for (size_t k = 0; k < pv.keys.size(); ++k) {
        SSt& o = pv.v[k];
        const uint64_t key = pv.keys[k];
        if (o.F < HALF_NEG || !alive(key)) continue;
        const int c = cur.find(key);
        if (c < 0 || cur.v[c].bP < HALF_NEG) continue;
        const double bp = cur.v[c].bP;
        o.bF = log_add(o.bF, bp + qv + bw);
        add2b(qb, d, j, bp + o.F + bw);
        if (&qb == &Qp_b) qps(2, d, j, bp + o.F + bw);
      }
    };
    const OMap& atb = Qb_at[j];
    for (size_t r = 0; r < atb.slots(); ++r) {
      const int d = atb.keys[r];
      if (d < 0 || !E.end_ok(d, j)) continue;
      branch(d, atb.vals[r], E.w_nest + E.at(d, j), Qb_b);
    }
    if (allow_pk) {
      const OMap& atp = Qp_at[j];
      for (size_t r = 0; r < atp.slots(); ++r) {
        const int d = atp.keys[r];
        if (d < 0) continue;
        branch(d, atp.vals[r], E.w_beta1p + E.nest_pk_k * E.w_nest, Qp_b);
      }
    }
    if (j < 1) return;
    SMap& pv = T[j - 1];
    for (size_t k = 0; k < pv.keys.size(); ++k) {
      const uint64_t key = pv.keys[k];
      if (!alive(key)) continue;
      const int c = cur.find(key);
      if (c < 0) continue;
      const SSt& s = cur.v[c];
      SSt& o = pv.v[k];
      if (dz) {
        if (s.bU1 > HALF_NEG) {
          o.bP = log_add(o.bP, s.bU1 + E.wu);
          o.bc = log_add(o.bc, s.bU1 + E.wu);
        }
        if (s.bU2 > HALF_NEG) {
          o.bU2 = log_add(o.bU2, s.bU2 + E.wu);
          o.bU1 = log_add(o.bU1, s.bU2 + E.wu + E.d5w[j - 1]);
        }
      } else if (s.bP > HALF_NEG && o.F > HALF_NEG) {
        o.bF = log_add(o.bF, s.bP + E.wu);
      }
    }
  }

  // --l2r-spawn inc, outside at the start of step j: P3 (filler, then the hole closed at
  // x' = j into P2), P2 (carried into P2[j-1], or opened at xs = j from P1[j-1]), P1.
  void l2r_inc_back(int j) {
    l2r_fext_back(LP3, j, [&](uint64_t key) { return lpa3(k3a(key), k3b(key), k3c(key), j); });
    SMap& p2 = LP2[j];
    for (size_t k = 0; k < LP3[j].keys.size(); ++k) {
      const uint64_t key = LP3[j].keys[k];
      const SSt& s = LP3[j].v[k];
      if (k3c(key) != j || s.c < HALF_NEG || s.bc < HALF_NEG) continue;
      const int t = p2.find(key3(k3a(key), k3b(key), 0));
      if (t >= 0) p2.v[t].bF = log_add(p2.v[t].bF, s.bc);
    }
    if (j >= 1) {
      SMap& pv = LP2[j - 1];
      SMap& p1 = LP1[j - 1];
      for (size_t k = 0; k < p2.keys.size(); ++k) {
        const double b = p2.v[k].bF;
        if (b < HALF_NEG) continue;
        const uint64_t key = p2.keys[k];
        SMap& src = k3b(key) < j ? pv : p1;
        const int t = src.find(k3b(key) < j ? key : key3(k3a(key), 0, 0));
        if (t >= 0) src.v[t].bF = log_add(src.v[t].bF, b);
      }
    }
    l2r_fext_back(LP1, j, [&](uint64_t key) { return lpa1(k3a(key), j); });
  }

  // Outside of the A seeds from P3: Rs(a2, xs, x', j) -> P3(a2, xs, x', j-1).
  void l2r_inc_seed_back(int j) {
    if (j < 1) return;
    SMap& pv = LP3[j - 1];
    for (Grp& g : LA.Gs[j].grps)
      for (Ent& s : g.ents) {
        if (s.b < HALF_NEG) continue;
        const int k = pv.find(key3(s.p, s.u, g.v));
        if (k < 0) continue;
        SSt& o = pv.v[k];
        o.bF = log_add(o.bF, s.b + (pair_bonus + E.w_stem + E.at(s.p, j)));
      }
  }

  void l2r_closeA_back(int j) {
    const SMap& S = LS[j];
    for (int src = 0; src < 2; ++src) {
      if (src == 1 && !single_stems) continue;
      for (Grp& g : (src ? LA.Gs[j] : LA.G[j]).grps)
        for (Ent& s : g.ents) {
          if (!E.end_ok(s.p, j)) continue;
          const int k = S.find(key3(s.p, s.u, g.v));
          if (k >= 0) add4(s, S.v[k].bc + E.at(s.p, j));
        }
    }
  }

  void l2r_bseed_back(int j) {
    if (j < 1) return;
    SMap& pv = LS[j - 1];
    for (Grp& g : LB.Gs[j].grps)
      for (Ent& s : g.ents) {
        if (s.b < HALF_NEG) continue;
        const int k = pv.find(key3(s.u, g.v, s.p));
        if (k < 0) continue;
        SSt& o = pv.v[k];
        o.bF = log_add(o.bF, s.b + E.w_stem + E.at(s.p, j) + pair_bonus);
      }
  }

  void l2r_closeB_back(int j) {
    for (int src = 0; src < 2; ++src) {
      if (src == 1 && !single_stems) continue;
      for (Grp& g : (src ? LB.Gs[j] : LB.G[j]).grps) {
        const int xs = g.v;
        if (!E.end_ok(xs, j)) continue;
        for (Ent& s : g.ents)
          if (s.p == xs) {
            const double qpb = lsp_use >= 0 ? QpS[lsp_use][ix(s.u, j)]
                                            : static_cast<double>(Qp_b[ix(s.u, j)]);
            add4(s, qpb + E.at(xs, j));
          }
      }
    }
  }

  // Start of the outside step j (before the engine's): Urs, Rrs, S.
  void l2r_back_end(int j) {
    l2r_adj = true;
    l2r_rs(LB, j);
    l2r_rs(LA, j);
    l2r_sext_back(j);
    if (l2r_inc) l2r_inc_back(j);
    l2r_adj = false;
  }

  // At the engine's Qp adjoint (Qp outside at j final): close B, U side, B seeds, close A, R side.
  void l2r_back_pk(int j) {
    l2r_adj = true;
    if (rlt) {
      l2r_tail_back(LB, j);
      l2r_tail_back(LA, j);
    }
    l2r_closeB_back(j);
    l2r_ls(LB, j);
    l2r_grow(LB, j);
    l2r_bseed_back(j);
    l2r_closeA_back(j);
    l2r_ls(LA, j);
    if (l2r_inc) l2r_inc_seed_back(j);
    else l2r_spawn_back(j);
    l2r_grow(LA, j);
    l2r_adj = false;
  }

  // Stem pairs: the outer pair of every R / U state and the inner pair of every seed.
  template <class AddPair>
  void l2r_pairs(AddPair& add_pair) {
    for (LFam* F : {&LA, &LB})
      for (int j = 0; j < n; ++j)
        for (const std::vector<Tab4>* T : {&F->G, &F->Gs})
          for (const Grp& g : (*T)[j].grps)
            for (const Ent& s : g.ents) add_pair(s.p, j, s.a + s.b - logZ);
  }

  // states per table: R, Rs, Rls, Rlss, Rrs, U, Us, Uls, Ulss, Urs, S
  void log_nest(const char* tag) const {
    auto oml = [](const std::vector<OMap>& v) {
      long long e = 0;
      for (const OMap& m : v) e += m.live;
      return e;
    };
    const long long nn = static_cast<long long>(n) * n;
    std::fprintf(stderr,
                 "[nest] %s n %d dense %lld Qb %lld Qm %lld Qm2 %lld Qp %lld Qmp %lld Qz %lld "
                 "qz_mb %.3f\n",
                 tag, n, nn, oml(Qb_at), oml(Qm_at), oml(Qm2_at), oml(Qp_at), oml(Qmp_at),
                 Qz.finite(), static_cast<double>(Qz.bytes()) / 1048576.0);
  }

  // pr_Mbar and pr_QbIn are built later, inside the principled pass. Time and size them
  // here so the pass-1 peak above this call does not include the two dense copies.
  void log_pr_side() const {
    const auto t0 = std::chrono::steady_clock::now();
    const std::vector<double> qbin = Qb.dense();
    const double tq = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const auto t1 = std::chrono::steady_clock::now();
    std::vector<double> mb(static_cast<size_t>(n) * static_cast<size_t>(n), NEG_INF);
    long long nz = 0;
    for (int j = 0; j < n; ++j) {
      const OMap& at = qmp_at(j);
      for (size_t r = 0; r < at.slots(); ++r) {
        const int i = at.keys[r];
        if (i < 0) continue;
        mb[static_cast<size_t>(i) * static_cast<size_t>(n) + static_cast<size_t>(j)] = at.vals[r];
        ++nz;
      }
    }
    const double tm = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
    std::fprintf(stderr,
                 "[pr-side] qbin_s %.6f qbin_mb %.1f mbar_s %.6f mbar_mb %.1f mbar_nz %lld\n",
                 tq, static_cast<double>(qbin.size()) * 8.0 / 1048576.0, tm,
                 static_cast<double>(mb.size()) * 8.0 / 1048576.0, nz);
  }

  std::vector<long long> l2r_states() const {
    std::vector<long long> c(11, 0);
    auto cnt = [](const std::vector<Tab4>& v) {
      long long s = 0;
      for (const Tab4& t : v)
        for (const Grp& g : t.grps) s += static_cast<long long>(g.ents.size());
      return s;
    };
    const LFam* fs[2] = {&LA, &LB};
    for (int k = 0; k < 2; ++k) {
      c[5 * k + 0] = cnt(fs[k]->G);
      c[5 * k + 1] = cnt(fs[k]->Gs);
      c[5 * k + 2] = cnt(fs[k]->Gls);
      c[5 * k + 3] = cnt(fs[k]->Glss);
      c[5 * k + 4] = cnt(fs[k]->Grs);
    }
    for (const SMap& m : LS) c[10] += static_cast<long long>(m.keys.size());
    return c;
  }

  // --mem-log: per structure entries and bytes (MemItem).  tabs: the engine tables (inside
  // dense + kept maps, outside, b16 4-index); pk: V / arms, K, growth, l2r tables, caches.
  static void tab4_sz(const std::vector<Tab4>& v, long long& e, double& b) {
    b += msz::vec(v);
    for (const Tab4& t : v) {
      b += msz::vec(t.grps) + msz::umap(t.vidx);
      for (const Grp& g : t.grps) {
        e += static_cast<long long>(g.ents.size());
        b += msz::vec(g.ents) + msz::umap(g.idx);
      }
    }
  }
  static void omap_sz(const std::vector<OMap>& v, long long& e, double& b) {
    b += msz::vec(v);
    for (const OMap& m : v) {
      e += m.live;
      b += msz::vec(m.keys) + msz::vec(m.vals) + m.pos.bytes();
    }
  }
  template <class M>
  static double pmap_sz(const M& m, long long& e) {
    double b = msz::umap(m);
    for (const auto& kv : m) {
      e += static_cast<long long>(kv.second.k.size());
      b += msz::vec(kv.second.k) + msz::vec(kv.second.w) + msz::vec(kv.second.lin);
    }
    return b;
  }
  template <class M>
  static double rmap_sz(const M& m, long long& e) {
    double b = msz::umap(m);
    for (const auto& kv : m) {
      e += relems(kv.second);
      for (int s = 0; s < 3; ++s) b += msz::vec(kv.second.L[s].v) + msz::vec(kv.second.R[s].v);
      b += msz::vec(kv.second.T[0]) + msz::vec(kv.second.T[1]);
    }
    return b;
  }
  template <class M>
  static double profmap_sz(const M& m, long long& e) {
    double b = msz::umap(m);
    for (const auto& kv : m) {
      const Prof& p = kv.second;
      e += static_cast<long long>(p.ly.size() + p.rx.size());
      b += msz::vec(p.ly) + msz::vec(p.rx) + msz::vec(p.lw) + msz::vec(p.rw);
    }
    return b;
  }
  template <class M>
  static double vmap_sz(const M& m, long long& e) {
    double b = msz::umap(m);
    for (const auto& kv : m) {
      e += static_cast<long long>(kv.second.size());
      b += msz::vec(kv.second);
    }
    return b;
  }
  void mem_items(std::vector<MemItem>& out, const std::string& pre, bool tabs, bool pk) const {
    if (tabs) {
      long long e = Qz.finite();
      double b = Qz.bytes() + msz::vec(C);
      for (int k = 0; k < 2; ++k) b += msz::vec(zP[k]) + msz::vec(zU1[k]) + msz::vec(zU2[k]);
      for (int r = 0; r < RS_N; ++r) b += msz::vec(runP[r]) + msz::vec(runU1[r]) + msz::vec(runU2[r]) + msz::vec(runW[r]);
      for (const std::vector<OMap>* v : {&Qb_at, &Qm_at, &Qm2_at, &Qp_at, &Qmp_at, &PKM_at}) omap_sz(*v, e, b);
      out.push_back({pre + "in", e, b});
      e = 0;
      b = Qz_b.bytes() + msz::vec(C_b);
      e += Qz_b.finite();
      for (const SpMat* v : {&Qb_b, &Qm_b, &Qm2_b, &Qp_b, &PKM_b, &Qmp_b, &QmC_b, &Qm2C_b}) {
        e += v->finite();
        b += v->bytes();
      }
      for (const std::vector<double>* v : {&QpS[0], &QpS[1], &QpS[2]}) {
        e += msz::fin(*v);
        b += msz::vec(*v);
      }
      for (int k = 0; k < 2; ++k) b += msz::vec(bz1[k]) + msz::vec(bz2[k]);
      out.push_back({pre + "out", e, b});
      e = 0;
      b = 0.0;
      for (const std::vector<Tab4>* v : {&Qg, &Qgls, &Qgrs, &Qgl, &Qgr}) tab4_sz(*v, e, b);
      out.push_back({pre + "tab4", e, b});
    }
    if (!pk) return;
    {
      double b = msz::umap(lArm) + msz::vec(lYwS) + msz::vec(prYw) + msz::vec(prEw);
      for (const std::vector<int>* v : {&prArmXs, &prArmXp, &prArmY0, &prYy, &prArmE0, &prEy, &prEys,
                                        &prXsFirst, &prByXp, &prXpFirst})
        b += msz::vec(*v);
      out.push_back({"V_arms", static_cast<long long>(l2r_arm_ents()), b});
    }
    {
      double b = msz::vec(prSF) + msz::umap(prAQ);
      long long e = static_cast<long long>(prSF.size());
      for (const auto& kv : prAQ) {
        const ArmQ& a = kv.second;
        e += static_cast<long long>(a.phi.size());
        b += msz::vec(a.phi) + msz::vec(a.rs) + msz::vec(a.t1) + msz::vec(a.t2) + msz::vec(a.ra) +
             msz::vec(a.P1) + msz::vec2(a.D) + vmap_sz(a.W, e);
      }
      out.push_back({"K_sfx", e, b});
    }
    long long e = 0;
    double b = rmap_sz(prRG, e);
    out.push_back({"grow_nodes", e, b});
    e = 0;
    b = rmap_sz(prRR, e);
    out.push_back({"grow_roots", e, b});
    e = 0;
    b = rmap_sz(prRH, e);
    out.push_back({"grow_help", e, b});
    e = 0;
    b = rmap_sz(prRKls, e) + rmap_sz(prRKrs, e);
    out.push_back({"grow_keys", e, b});
    e = 0;
    b = rside_bytes(e);
    out.push_back({"grow_R_side", e, b});
    e = 0;
    b = rfv_bytes(e);
    out.push_back({"V_keep", e, b});
    e = 0;
    b = msz::umap(prRT);
    for (const auto& kv : prRT) {
      e += static_cast<long long>(kv.second.L.size() + kv.second.R.size());
      b += msz::vec(kv.second.L) + msz::vec(kv.second.R);
    }
    out.push_back({"grow_trunc", e, b});
    e = 0;
    b = 0.0;
    for (const LFam* F : {&LA, &LB})
      for (const std::vector<Tab4>* v : {&F->G, &F->Gs, &F->Gls, &F->Glss, &F->Grs}) tab4_sz(*v, e, b);
    b += msz::vec(LS);
    for (const SMap& m : LS) {
      e += static_cast<long long>(m.keys.size());
      b += msz::vec(m.keys) + msz::vec(m.v) + msz::umap(m.idx);
    }
    out.push_back({"l2r_tab", e, b});
    if (rlt) {
      e = 0;
      b = 0.0;
      for (const LFam* F : {&LA, &LB}) tab4_sz(F->Glt, e, b);
      out.push_back({"l2r_tail", e, b});
    }
    if (l2r_inc) {
      e = 0;
      b = 0.0;
      for (const std::vector<SMap>* T : {&LP1, &LP2, &LP3}) {
        b += msz::vec(*T);
        for (const SMap& m : *T) {
          e += static_cast<long long>(m.keys.size());
          b += msz::vec(m.keys) + msz::vec(m.v) + msz::umap(m.idx);
        }
      }
      out.push_back({"l2r_P", e, b});
      e = 0;
      b = msz::vec(lH) + msz::vec(lHdone) + msz::vec2(lArmOpen) + msz::vec(lBopen);
      for (const LH& h : lH) {
        b += msz::vec(h.c1) + msz::vec(h.aopen) + msz::vec2(h.v);
        for (const auto& v : h.v) e += static_cast<long long>(v.size());
      }
      out.push_back({"l2r_H", e, b});
    }
    e = 0;
    b = msz::vec2(spG) + msz::vec2(spCol) + msz::vec2(spRow) + msz::vec(spWmax);
    out.push_back({"spawn", static_cast<long long>(sp_row_ents), b});
    e = 0;
    b = msz::vec2(lWr) + msz::vec2(lWx) + msz::umap(lKbr) + msz::umap(lKbx) + msz::umap(lFb) +
        msz::umap(lKc) + msz::umap(lOU) + msz::umap(lOUs);
    out.push_back({"l2r_ctx", static_cast<long long>(lKbr.size() + lKbx.size()), b});
    e = 0;
    b = pmap_sz(prG, e) + pmap_sz(prGR, e) + pmap_sz(prHg, e) + pmap_sz(prGRs, e) + pmap_sz(prKls, e) +
        pmap_sz(prKrs, e);
    out.push_back({"b16_grow_prof", e, b});
    e = 0;
    b = msz::umap(prTc);
    for (const auto& kv : prTc) {
      b += msz::umap(kv.second);
      for (const auto& t : kv.second) {
        e += static_cast<long long>(t.second.ys.size());
        b += msz::vec(t.second.ys) + msz::vec(t.second.w);
      }
    }
    out.push_back({"b16_T", e, b});
    e = 0;
    b = msz::umap(prKc);
    for (const auto& kv : prKc) {
      e += static_cast<long long>(kv.second.is.size());
      b += msz::vec(kv.second.is) + msz::vec(kv.second.w);
    }
    out.push_back({"b16_kappa_K", e, b});
    e = 0;
    b = msz::umap(prUL) + msz::vec(prUR) + msz::umap(prUc) + msz::vec(prEBr) + msz::vec(prEBc) +
        msz::vec(prMBr) + msz::vec(prMBc) + msz::vec(prKap);
    for (const auto& kv : prUL) {
      e += static_cast<long long>(kv.second.w.size());
      b += msz::vec(kv.second.w) + msz::vec(kv.second.lin);
    }
    for (const PrUR& r : prUR) {
      e += static_cast<long long>(r.w.size());
      b += msz::vec(r.w) + msz::vec(r.lin);
    }
    out.push_back({"grow_sel_U", e, b});
    e = 0;
    b = msz::umap(stA) + msz::umap(stE) + msz::umap(stBc) + msz::umap(stBgc) + vmap_sz(stX, e) +
        msz::umap(stHc) + msz::umap(stOc) + profmap_sz(stP, e) + profmap_sz(stPg, e) + profmap_sz(stPls, e) +
        profmap_sz(stPrs, e) + vmap_sz(arm_cache, e) + msz::umap(stM);
    out.push_back({"b16_stem", e, b});
    e = 0;
    b = profmap_sz(chM, e) + profmap_sz(chR, e) + vmap_sz(chPaths, e) + vmap_sz(prnG, e) + vmap_sz(prnV, e) +
        msz::umap(prSc) + msz::umap(prAc);
    out.push_back({"b16_chain_naive", e, b});
  }

  // ------------------------------------------------------------------- Qg
  void fill_Qg(int j) {
    Tab4& G = Qg[j];
    const auto& pair_i = partners_left[j];

    // (a) interior loop between i·j and d·e
    for (int i : pair_i) {
      int d_hi = std::min(j - 1, i + 1 + cap);
      for (int d = i + 1; d <= d_hi; ++d) {
        int e_lo = std::max(d + MIN_PAIR, j - 1 - cap);
        for (int e = e_lo; e < j; ++e) G.acc(e, i, d, E.stem_interior(i, d, e, j) + 2.0 * pair_bonus);
      }
    }

    // (b) interior wrapping a stored Qg_{c,d,e,f}
    int f_lo = std::max(0, j - 1 - cap);
    for (int f = f_lo; f < j; ++f) {
      for (const Grp& g : Qg[f].grps) {
        int e = g.v;
        for (const Ent& s : g.ents) {
          int c = s.p, d = s.u;
          if (c == d) continue;
          double x = s.a;
          int i_lo = std::max(0, c - 1 - cap);
          for (int i = i_lo; i < c; ++i) {
            double w = E.stem_interior(i, c, f, j);
            G.acc(e, i, d, w + x + pair_bonus);
          }
        }
      }
    }

    // (f) Qgls plus unpaired run to j; seed-derived Qgls (e == f) are DP03 (c), run capped
    for (int f = 0; f < j; ++f) {
      int unpaired = j - f - 1;
      double w_ext = E.w_ap + E.w_bp + E.w_cp * unpaired + pair_bonus;
      if (dz) w_ext += E.run(f + 1, j - 1);
      for (const Grp& g : Qgls[f].grps) {
        int e = g.v;
        if (e == f && unpaired > cap) continue;
        for (const Ent& s : g.ents) {
          int i = s.p - 1;
          if (E.pairs(i, j) && E.end_ok(i, j)) G.acc(e, i, s.u, s.a + (w_ext + E.at(i, j)));
        }
      }
    }

    // (g) Qgrs at j−1, unpaired on the left
    if (j >= 1) {
      for (const Grp& g : Qgrs[j - 1].grps) {
        int e = g.v;
        for (const Ent& s : g.ents) {
          int c = s.p, d = s.u;
          double x = s.a;
          int i_lo = lead_i_lo(c);
          for (int i = i_lo; i < c; ++i) {
            if (!E.pairs(i, j) || !E.end_ok(i, j)) continue;
            double w = x + E.w_ap + E.w_bp + E.w_cp * (c - i - 1) + E.at(i, j) + pair_bonus;
            if (dz) w += E.run(i + 1, c - 1);
            G.acc(e, i, d, w);
          }
        }
      }
    }

    // (h) Qgrs + Qm on the left
    if (j >= 1) {
      for (const Grp& g : Qgrs[j - 1].grps) {
        int e = g.v;
        for (const Ent& s : g.ents) {
          int c = s.p, d = s.u;
          double x = s.a;
          if (c - 1 < 0) continue;
          const OMap& at = qmp_at(c - 1);
          qm_pair_slots(c - 1, j, [&](int r) {
            int i = at.keys[r];
            if (i < 0) return;
            double qm = at.vals[r];
            int outer = i - 1;
            if (E.pairs(outer, j) && outer < c && E.end_ok(outer, j))
              G.acc(e, outer, d, qm + x + E.w_ap + E.w_bp + E.at(outer, j) + pair_bonus);
          });
        }
      }
    }

    // (0) the pair i·j alone as a stem, Qg_{i,i,j,j} = 1
    if (single_stems) {
      for (int i : pair_i)
        if (E.end_ok(i, j)) G.acc(j, i, i, pair_bonus);
    }
  }

  void fill_Qgls(int j) {
    Tab4& out = Qgls[j];
    for (const Grp& g : Qg[j].grps) {
      int e = g.v;
      for (const Ent& s : g.ents) {
        int c = s.p, d = s.u;
        double qg = s.a;
        if (c - 1 < 0 || c == d || !E.end_ok(c, j)) continue;
        double w = E.w_bp + E.at(c, j);
        const OMap& at = qmp_at(c - 1);
        for (size_t r = 0; r < at.slots(); ++r) {
          int i = at.keys[r];
          if (i < 0) continue;
          out.acc(e, i, d, w + at.vals[r] + qg);
        }
      }
    }
    // seeds c·j (DP03 (c) once closed by (f))
    for (int c : partners_left[j]) {
      if (c - 1 < 0 || !E.end_ok(c, j)) continue;
      double w = E.w_bp + E.at(c, j) + pair_bonus;
      const OMap& at = qmp_at(c - 1);
      for (size_t r = 0; r < at.slots(); ++r) {
        int i = at.keys[r];
        if (i < 0) continue;
        out.acc(j, i, c, w + at.vals[r]);
      }
    }
  }

  void fill_Qgl(int j) {
    Tab4& out = Qgl[j];
    for (const Grp& g : Qg[j].grps) {
      int f = g.v;
      for (const Ent& s : g.ents) {
        int i = s.p, d = s.u;
        double qg = s.a;
        if (!E.end_ok(d, f)) continue;
        double w = E.w_stem + E.at(d, f);
        for (int e = d; e < f; ++e) {
          double qz = qz_log(d + 1, e);
          if (qz < HALF_NEG) continue;
          out.acc(f, i, e, qg + qz + w);
        }
      }
    }
  }

  void fill_Qgr(int j) {
    Tab4& out = Qgr[j];
    for (const Grp& g : Qgl[j].grps) {
      int f = g.v;
      for (const Ent& s : g.ents) {
        int i = s.p, d = s.u;
        double gl = s.a;
        for (int e = d + 1; e <= f; ++e) {
          double qz = qz_log(e, f - 1);
          if (qz < HALF_NEG) continue;
          out.acc(e, i, d, gl + qz);
        }
      }
    }
  }

  void fill_Qp(int j) {
    for (const Grp& g : Qgr[j].grps) {
      int e_gr = g.v;
      int f = e_gr - 1;
      if (f < 0) continue;
      Tab4& gl_right = Qgl[f];
      for (const Ent& s : g.ents) {
        int d = s.p, em1 = s.u;
        double gr = s.a;
        if (!E.end_ok(d, j)) continue;
        int e = em1 + 1;
        Grp* gl_group = gl_right.get(e);
        if (!gl_group || gl_group->ents.empty()) continue;
        for (const Ent& t : gl_group->ents) {
          int i = t.p, uu = t.u;
          if (uu != d - 1 || i > d - 1 || !E.end_ok(i, f)) continue;
          set_pair(Qp, Qp_at, i, j, t.a + gr + (E.at(i, f) + E.at(d, j)));
        }
      }
    }
  }

  void fill_Qgrs(int j) {
    Tab4& out = Qgrs[j];
    const OMap& at = qmp_at(j);
    for (size_t r = 0; r < at.slots(); ++r) {
      int p = at.keys[r];
      if (p < 0) continue;
      double qm = at.vals[r];
      int f = p - 1;
      if (f < 0) continue;
      for (const Grp& g : Qg[f].grps) {
        int e = g.v;
        for (const Ent& s : g.ents) {
          if (s.p == s.u || !E.end_ok(s.p, f)) continue;
          out.acc(e, s.p, s.u, s.a + qm + (E.w_bp + E.at(s.p, f)));
        }
      }
      // seeds c·f (DP03 (d), (e) once closed by (g), (h))
      for (int c : partners_left[f])
        if (E.end_ok(c, f)) out.acc(f, c, c, qm + (E.w_bp + E.at(c, f)) + pair_bonus);
    }
  }

  // ------------------------------------------------------------ two-index
  // Next index strictly after `after` that pairs with i, or -1. O(1) when the
  // base-type successor is a real pair; a restricted pair list may skip a few.
  void lp_build_mates() {
    if (mates_ready) return;
    mates_ready = true;
    for (int b = 0; b < 4; ++b) {
      mate_after[b].assign(static_cast<size_t>(n), -1);
      int nxt = -1;
      for (int j = n - 1; j >= 0; --j) {
        mate_after[b][static_cast<size_t>(j)] = nxt;
        const int c = E.code[j];
        const int t = c < 0 ? -1 : pair_type(b, c);
        if (t >= 0 && (allow_gu || (t != 3 && t != 5))) nxt = j;
      }
    }
  }
  int next_mate(int i, int after) const {
    if (i < 0 || after >= n - 1) return -1;
    const int b = E.code[i];
    if (b < 0) return -1;
    int j = mate_after[b][static_cast<size_t>(after)];
    while (j != -1 && !E.pairs(i, j)) j = j + 1 < n ? mate_after[b][static_cast<size_t>(j)] : -1;
    return j;
  }
  double unpaired_gap(int s, int t) const {
    if (t < s) return 0.0;
    return E.w_alpha3 * static_cast<double>(t - s + 1) + (dz ? E.run(s, t) : 0.0);
  }
  // Closing-pair left end for a first branch at `branch` (unpaired = branch - left - 1).
  int lead_i_lo(int branch) const {
    if (ml_lead_cap <= 0) return 0;
    return std::max(0, branch - ml_lead_cap);
  }
  // Qm-start form: unpaired bases from `start` up to but not including the branch.
  int lead_qm_lo(int branch) const {
    if (ml_lead_cap <= 0) return 0;
    return std::max(0, branch - (ml_lead_cap - 1));
  }
  bool lead_gap_ok(int gap) const {
    return ml_lead_cap <= 0 || gap <= ml_lead_cap - 1;
  }
  int lead_span() const { return ml_lead_cap <= 0 ? n + 1 : ml_lead_cap; }
  // Lazy left unpaired run. Off: the old window of `cap` bases (x-1-cap). On
  // (--lazy-ml-cap): the pass-2 leading run, i >= branch - ml_lead_cap.
  int ml_run_lo(int branch) const {
    return lazy_ml_cap ? lead_i_lo(branch) : std::max(0, branch - 1 - cap);
  }
  bool ml_run_ok(int gap) const {
    return lazy_ml_cap ? lead_gap_ok(gap) : gap <= cap;
  }
  void lp_beam(OMap& at) {
    const int B = tbeam[BT_QB];
    if (B <= 0 || at.live <= B) return;
    std::vector<std::pair<double, int>> ranked;
    ranked.reserve(static_cast<size_t>(at.live));
    for (size_t r = 0; r < at.slots(); ++r) {
      const int i = at.keys[r];
      if (i < 0) continue;
      const double s = pre_key_inside ? at.vals[r] : prefix_log(i) + at.vals[r];
      ranked.emplace_back(s, i);
    }
    std::sort(ranked.begin(), ranked.end(), [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
      if (a.first != b.first) return a.first > b.first;
      return a.second > b.second;
    });
    std::vector<char> keep(static_cast<size_t>(n), 0);
    const int keep_n = std::min(B, static_cast<int>(ranked.size()));
    for (int k = 0; k < keep_n; ++k) keep[static_cast<size_t>(ranked[k].second)] = 1;
    for (size_t r = 0; r < at.slots(); ++r) {
      const int i = at.keys[r];
      if (i < 0 || keep[static_cast<size_t>(i)]) continue;
      at.erase(i);
    }
    at.compact();
  }
  // Right-run class. Length >= 2 shares one class: its extension delta is
  // d3[new] - d3[old], independent of the start. Empty and one base do not.
  int lp_kind(int gap) const {
    if (gap <= 0) return 0;
    if (!dz) return 1;
    return gap == 1 ? 1 : 2;
  }
  // One next-mate step of the right run. nk is the class at qn.
  bool lp_step(int p, int q, int kind, int& qn, int& nk, double& delta) const {
    qn = next_mate(p, q);
    if (qn < 0) return false;
    const int add = qn - q;
    if (add <= 0) return false;
    delta = E.w_alpha3 * static_cast<double>(add);
    if (!dz) {
      nk = 1;
      return true;
    }
    if (kind == 0) {
      nk = add == 1 ? 1 : 2;
      delta += E.run(q, qn - 1);
    } else if (kind == 1) {
      nk = 2;
      delta += E.run(q - 1, qn - 1) - E.s1w[q - 1];
    } else {
      nk = 2;
      delta += E.d3w[qn - 1] - E.d3w[q - 1];
    }
    return true;
  }
  void mu_add(int q, int p, int kind, double w) {
    if (q < 0 || q >= n || p < 0 || w < HALF_NEG) return;
    if (lpM.empty()) lpM.assign(static_cast<size_t>(n), {});
    OMap& col = lpM[static_cast<size_t>(q)];
    const int k = p * 3 + kind;
    const double cur = col.get(k);
    col.set(k, cur < HALF_NEG ? w : log_add(cur, w));
  }
  // Beam a Multi column. The key is p*3+kind; the prefix is C[p-1].
  void lp_beam_multi(OMap& at) {
    const int B = tbeam[BT_QB];
    if (B <= 0 || at.live <= B) return;
    struct R {
      double s;
      int k, p, kind;
    };
    std::vector<R> ranked;
    ranked.reserve(static_cast<size_t>(at.live));
    for (size_t r = 0; r < at.slots(); ++r) {
      const int k = at.keys[r];
      if (k < 0) continue;
      const int p = k / 3;
      const double s = pre_key_inside ? at.vals[r] : prefix_log(p) + at.vals[r];
      ranked.push_back({s, k, p, k % 3});
    }
    std::sort(ranked.begin(), ranked.end(), [](const R& a, const R& b) {
      if (a.s != b.s) return a.s > b.s;
      if (a.p != b.p) return a.p > b.p;
      return a.kind > b.kind;
    });
    std::vector<int> drop;
    for (size_t i = static_cast<size_t>(std::min(B, static_cast<int>(ranked.size())));
         i < ranked.size(); ++i)
      drop.push_back(ranked[i].k);
    for (int k : drop) at.erase(k);
    at.compact();
  }
  // Prune H[j], seed H(j, first mate), extend only the kept hairpins, then free H[j].
  void lp_fill_qb(int j) {
    lp_build_mates();
    if (lpH.empty()) lpH.assign(static_cast<size_t>(n), {});
    if (lpM.empty()) lpM.assign(static_cast<size_t>(n), {});
    OMap& col = lpH[static_cast<size_t>(j)];
    lp_beam(col);
    const auto th = std::chrono::steady_clock::now();
    const int q0 = next_mate(j, j);
    if (q0 >= 0) {
      const double h0 = E.hairpin(j, q0);
      if (h0 > HALF_NEG) lpH[static_cast<size_t>(q0)].set(j, h0);
    }
    for (size_t r = 0; r < col.slots(); ++r) {
      const int i = col.keys[r];
      if (i < 0) continue;
      ++r_hair;
      set_pair(Qb, Qb_at, i, j, col.vals[r]);
      const int q = next_mate(i, j);
      if (q >= 0) {
        const double h = E.hairpin(i, q);
        if (h > HALF_NEG) lpH[static_cast<size_t>(q)].set(i, h);
      }
    }
    sec_hair += std::chrono::duration<double>(std::chrono::steady_clock::now() - th).count();
    col.init(0);
    const auto tc = std::chrono::steady_clock::now();
    lp_close_col(j);
    sec_mclose += std::chrono::duration<double>(std::chrono::steady_clock::now() - tc).count();
    // pair_bonus is one add on the finished pair sum (hairpin + interior + close),
    // the same place as set_pair(acc + pair_bonus) in the scan fill.
    if (pair_bonus != 0.0) {
      OMap& qb = Qb_at[j];
      for (size_t r = 0; r < qb.slots(); ++r) {
        const int i = qb.keys[r];
        if (i < 0) continue;
        const double v = qb.vals[r] + pair_bonus;
        qb.vals[r] = v;
        Qb[ix(i, j)] = v;
      }
    }
  }
  // Kept Multi(p, j) closes into Qb(p, j), then steps to the next mate.
  // The column that remains is the beam output; discarded arrivals are erased.
  void lp_close_col(int j) {
    if (lpM.empty()) return;
    OMap& col = lpM[static_cast<size_t>(j)];
    lp_beam_multi(col);
    const double base = E.w_alpha1 + E.w_alpha2;
    for (size_t r = 0; r < col.slots(); ++r) {
      const int k = col.keys[r];
      if (k < 0) continue;
      const int p = k / 3;
      const int kind = k % 3;
      const double w = col.vals[r];
      if (E.end_ok(p, j)) set_pair(Qb, Qb_at, p, j, w + base + E.at_close(p, j));
      int qn, nk;
      double delta;
      if (lp_step(p, j, kind, qn, nk, delta)) mu_add(qn, p, nk, w + delta);
    }
  }
  // From each kept pair, the outer pair sits inside the loop cap. Push it forward.
  // Same shape as LinearPartition handleP's internal/stack loop, with max_loop as the cap.
  void lp_push_interior(int j) {
    const auto ti = std::chrono::steady_clock::now();
    lp_build_mates();
    const OMap& at = Qb_at[j];
    for (size_t r = 0; r < at.slots(); ++r) {
      const int i = at.keys[r];
      if (i < 0) continue;
      const double qb = at.vals[r];
      if (qb < HALF_NEG) continue;
      const int p_lo = std::max(0, i - 1 - cap);
      for (int p = i - 1; p >= p_lo; --p) {
        const int l1 = i - p - 1;
        const int budget = cap - l1;
        if (budget < 0) break;
        for (int q = next_mate(p, j); q != -1 && q - j - 1 <= budget; q = next_mate(p, q)) {
          ++r_int;
          const double w = E.interior(p, i, j, q);
          if (w < HALF_NEG) continue;
          set_pair(Qb, Qb_at, p, q, w + qb);
        }
      }
    }
    sec_int += std::chrono::duration<double>(std::chrono::steady_clock::now() - ti).count();
  }
  // M2 that ends on a branch (lpP, no trailing run) emits Multi(p, next_mate(p, j)).
  // The right run j+1..q-1 is on that edge. The leading run is p+1..ib-1, capped
  // by ml_lead_cap. Arrivals at the same (p, q, kind) merge.
  void lp_emit_multi(int j) {
    if (lpP.empty() || j < 0 || j >= n) return;
    const auto te = std::chrono::steady_clock::now();
    lp_build_mates();
    if (lpM.empty()) lpM.assign(static_cast<size_t>(n), {});
    OMap& src = lpP[static_cast<size_t>(j)];
    std::vector<int> drop;
    for (size_t r = 0; r < src.slots(); ++r) {
      const int i = src.keys[r];
      if (i < 0) continue;
      if (!Qm2_at[j].has(i)) drop.push_back(i);
    }
    for (int i : drop) src.erase(i);
    src.compact();
    const int L = lead_span();
    for (size_t r = 0; r < src.slots(); ++r) {
      const int ib = src.keys[r];
      if (ib < 0) continue;
      const double mv = src.vals[r];
      if (mv < HALF_NEG) continue;
      const int p_lo = std::max(0, ib - L);
      for (int p = ib - 1; p >= p_lo; --p) {
        const int q = next_mate(p, j);
        if (q < 0) continue;
        ++r_m2multi;
        const int kind = lp_kind(q - j - 1);
        const double w = mv + unpaired_gap(p + 1, ib - 1) + unpaired_gap(j + 1, q - 1);
        mu_add(q, p, kind, w);
      }
    }
    sec_emit += std::chrono::duration<double>(std::chrono::steady_clock::now() - te).count();
  }

  void lp_mb_add(int q, int k, double w) {
    if (q < 0 || q >= n || w < HALF_NEG) return;
    if (lpMb.empty()) lpMb.assign(static_cast<size_t>(n), {});
    OMap& col = lpMb[static_cast<size_t>(q)];
    const double cur = col.get(k);
    col.set(k, cur < HALF_NEG ? w : log_add(cur, w));
  }
  // Outside of Multi(p, j): the closing pair, plus the jump to the next mate.
  void lp_out_multi(int j) {
    if (l2r_mode != 0 || lpM.empty() || j < 0 || j >= static_cast<int>(lpM.size())) return;
    if (lpMb.empty()) lpMb.assign(static_cast<size_t>(n), {});
    const OMap& col = lpM[static_cast<size_t>(j)];
    for (size_t r = 0; r < col.slots(); ++r) {
      const int k = col.keys[r];
      if (k < 0) continue;
      const int p = k / 3;
      const int kind = k % 3;
      if (E.end_ok(p, j)) {
        const double bb = Qb_b[ix(p, j)];
        if (bb > HALF_NEG)
          lp_mb_add(j, k, bb + E.w_alpha1 + E.w_alpha2 + E.at_close(p, j));
      }
      int qn, nk;
      double delta;
      if (!lp_step(p, j, kind, qn, nk, delta)) continue;
      const int kn = p * 3 + nk;
      if (!lpM[static_cast<size_t>(qn)].has(kn)) continue;
      const double cb = lpMb[static_cast<size_t>(qn)].get(kn);
      if (cb > HALF_NEG) lp_mb_add(j, k, cb + delta);
    }
  }
  // E3: Multi(p, q) -> M2(ib, j), q = next_mate(p, j). Qm2_b is already the
  // branch-ending (P) outside, so this does not flow into the trailing run.
  void lp_out_e3(int j) {
    if (l2r_mode != 0 || lpP.empty() || lpM.empty() || lpMb.empty()) return;
    if (j < 0 || j >= n) return;
    const int L = lead_span();
    const OMap& src = lpP[static_cast<size_t>(j)];
    for (size_t r = 0; r < src.slots(); ++r) {
      const int ib = src.keys[r];
      if (ib < 0 || !Qm2_at[j].has(ib)) continue;
      const int p_lo = std::max(0, ib - L);
      for (int p = ib - 1; p >= p_lo; --p) {
        const int q = next_mate(p, j);
        if (q < 0) continue;
        const int k = p * 3 + lp_kind(q - j - 1);
        if (q >= n || !lpM[static_cast<size_t>(q)].has(k)) continue;
        const double mb = lpMb[static_cast<size_t>(q)].get(k);
        if (mb < HALF_NEG) continue;
        const double factor = unpaired_gap(p + 1, ib - 1) + unpaired_gap(j + 1, q - 1);
        add2b(Qm2_b, ib, j, mb + factor);
      }
    }
  }

  void fill_Qb(int j) {
    // Nested pass only. --pk-l2r off is Beam-DP03: hairpin, interior, Qm2, and PKM.
    if (l2r_mode == 0 && !allow_pk) {
      lp_fill_qb(j);
      return;
    }
    // Interior spans are bounded by the loop cap, so the outer left end sits in a
    // window of width cap next to a kept inner pair. Walk those kept pairs and the
    // partners inside the window; do not walk every left partner of j.
    if (static_cast<int>(qb_acc.size()) < n) qb_acc.assign(static_cast<size_t>(n), NEG_INF);
    const std::vector<int>& pl = partners_left[j];
    for (int i : pl) {
      ++r_hair;
      qb_acc[static_cast<size_t>(i)] = E.hairpin(i, j);
    }
    const int e0 = std::max(0, j - 1 - cap);
    for (int e = e0; e < j; ++e) {
      const int budget = cap - (j - e - 1);
      const OMap& at = Qb_at[e];
      for (size_t r = 0; r < at.slots(); ++r) {
        const int d = at.keys[r];
        if (d < 0) continue;
        const double qb = at.vals[r];
        const int i_hi = std::min(d - 1, e - MIN_PAIR);
        const int i_lo = std::max(0, d - 1 - budget);
        if (i_lo > i_hi) continue;
        auto it = std::lower_bound(pl.begin(), pl.end(), i_lo);
        const auto iend = std::upper_bound(it, pl.end(), i_hi);
        for (; it != iend; ++it) {
          ++r_int;
          const int i = *it;
          const double w = E.interior(i, d, e, j);
          if (w < HALF_NEG) continue;
          double& acc = qb_acc[static_cast<size_t>(i)];
          acc = log_add(acc, w + qb);
        }
      }
    }
    for (int i : pl) {
      double acc = qb_acc[static_cast<size_t>(i)];
      qb_acc[static_cast<size_t>(i)] = NEG_INF;

      const bool closing_ok = E.end_ok(i, j);
      const double at_ij = E.at_close(i, j);

      if (closing_ok && i + 1 <= j - 1) {
        ++r_m2multi;
        double m2 = Qm2[ix(i + 1, j - 1)];
        if (m2 > HALF_NEG) acc = log_add(acc, m2 + E.w_alpha1 + E.w_alpha2 + at_ij);
      }

      if (allow_pk && closing_ok && j >= 1) {
        const OMap& pkm = PKM_at[j - 1];
        for (size_t r = 0; r < pkm.slots(); ++r) {
          int d = pkm.keys[r];
          if (d <= i) continue;
          double pk = pkm.vals[r];
          double w_pk = pk + E.w_alpha1 + E.w_beta1m + (1.0 + E.ml_pk_k) * E.w_alpha2 +
                        E.w_alpha3 * (d - i - 1) + at_ij;
          if (dz) w_pk += E.run(i + 1, d - 1);
          acc = log_add(acc, w_pk);
          if (i + 1 <= d - 1) {
            double left = Qm[ix(i + 1, d - 1)];
            if (left > HALF_NEG)
              acc = log_add(acc, left + pk + E.w_alpha1 + E.w_beta1m +
                                     (1.0 + E.ml_pk_k) * E.w_alpha2 + at_ij);
          }
        }
      }
      set_pair(Qb, Qb_at, i, j, acc + pair_bonus);
    }
  }

  void build_qlead(int t) {
    if (qlead.empty()) qlead.assign(static_cast<size_t>(n), {});
    const double w3 = E.w_alpha3;
    std::vector<LeadEnt> cand;
    if (t == 0) {
      cand.push_back({0, 0.0});
    } else {
      cand.reserve(qlead[static_cast<size_t>(t - 1)].size() + 1);
      for (const LeadEnt& e : qlead[static_cast<size_t>(t - 1)]) {
        if (t - e.x0 > E.ml_left || !lead_gap_ok(t - e.x0)) continue;
        double w = e.w + w3;
        if (dz) w += E.run(e.x0, t - 1) - E.run(e.x0, t - 2);
        cand.push_back({e.x0, w});
      }
      cand.push_back({t, 0.0});
    }
    const int B = p2_qz_beam;
    if (B > 0 && static_cast<int>(cand.size()) > B) {
      auto key = [&](const LeadEnt& e) {
        return e.w + (e.x0 > 0 ? C[static_cast<size_t>(e.x0 - 1)] : 0.0);
      };
      std::nth_element(cand.begin(), cand.begin() + B, cand.end(), [&](const LeadEnt& a, const LeadEnt& b) {
        const double ka = key(a), kb = key(b);
        if (ka != kb) return ka > kb;
        return a.x0 < b.x0;
      });
      cand.resize(static_cast<size_t>(B));
    }
    std::sort(cand.begin(), cand.end(), [](const LeadEnt& a, const LeadEnt& b) { return a.x0 < b.x0; });
    qlead[static_cast<size_t>(t)] = std::move(cand);
  }

  void add_branch(int j, int d, double obj, bool pk, double w_at) {
    const double w2 = E.w_alpha2;
    double w_single = obj + ((pk ? E.w_beta1m + E.ml_pk_k * w2 : w2) + w_at);
    int lead_lim = E.ml_left;
    if (ml_lead_cap > 0 && ml_lead_cap - 1 < lead_lim) lead_lim = ml_lead_cap - 1;
    const int lo = std::max(0, d - lead_lim);
    if (l2r_mode == 0 && !allow_pk) {
      // LinearPartition M = P: the fragment starts at the pair. Unpaired bases sit to its right.
      const auto tm = std::chrono::steady_clock::now();
      ++r_mp;
      set_Qm(d, j, w_single, false);
      sec_mp += std::chrono::duration<double>(std::chrono::steady_clock::now() - tm).count();
    } else if (static_cast<int>(qlead.size()) == n) {
      for (const LeadEnt& e : qlead[static_cast<size_t>(d)]) {
        if (e.x0 < lo || e.x0 > d) continue;
        ++r_mp;
        set_Qm(e.x0, j, w_single + e.w, false);
      }
    } else {
      const double w3 = E.w_alpha3;
      for (int i = lo; i <= d; ++i) {
        double w = w_single + w3 * (d - i);
        if (dz && i < d) w += E.run(i, d - 1);
        ++r_mp;
        set_Qm(i, j, w, false);
      }
    }
    if (d == 0) return;
    double w_more = obj + ((pk ? E.w_beta1m + E.ml_pk_k * w2 : w2) + w_at);
    const OMap& at = Qm_at[d - 1];
    const auto tx = l2r_mode == 0 ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
    for (size_t r = 0; r < at.slots(); ++r) {
      int i = at.keys[r];
      if (i < 0) continue;
      ++r_m2xp;
      const double nv = at.vals[r] + w_more;
      set_Qm(i, j, nv, !pk);
      if (!pk && nv > HALF_NEG) {
        if (lpP.empty()) lpP.assign(static_cast<size_t>(n), {});
        OMap& pc = lpP[static_cast<size_t>(j)];
        const double cur = pc.get(i);
        pc.set(i, cur < HALF_NEG ? nv : log_add(cur, nv));
      }
    }
    if (l2r_mode == 0)
      sec_m2xp += std::chrono::duration<double>(std::chrono::steady_clock::now() - tx).count();
  }

  void fill_Qm(int j) {
    const double w3 = E.w_alpha3;
    std::vector<int> surv, surv2;
    const auto t_ext = std::chrono::steady_clock::now();
    if (dz) {
      RunSt& cu = rsQm[j & 1];
      RunSt& cu2 = rsQm2[j & 1];
      cu.clear();
      cu2.clear();
      if (j > 0) {
        surv = live_keys(Qm_at[j - 1]);
        surv2 = live_keys(Qm2_at[j - 1]);
        r_mmu += static_cast<long long>(surv.size());
        r_muext += static_cast<long long>(surv2.size());
        run_ext(rsQm[(j - 1) & 1], cu, surv, j, w3);
        run_ext(rsQm2[(j - 1) & 1], cu2, surv2, j, w3);
      }
    } else if (j > 0) {
      const OMap& at = Qm_at[j - 1];
      for (size_t r = 0; r < at.slots(); ++r) {
        int i = at.keys[r];
        if (i < 0) continue;
        ++r_mmu;
        set_Qm(i, j, at.vals[r] + w3, false);
      }
      const OMap& at2 = Qm2_at[j - 1];
      for (size_t r = 0; r < at2.slots(); ++r) {
        int i = at2.keys[r];
        if (i < 0) continue;
        ++r_muext;
        double nv = at2.vals[r] + w3;
        double t = log_add(static_cast<double>(Qm2[ix(i, j)]), nv);
        Qm2[ix(i, j)] = t;
        Qm2_at[j].set(i, t);
      }
    }
    sec_qm += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_ext).count();
    const auto t_br = std::chrono::steady_clock::now();
    {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        int d = at.keys[r];
        if (d < 0 || !E.end_ok(d, j)) continue;
        add_branch(j, d, at.vals[r], false, E.at(d, j));
      }
    }
    if (allow_pk) {
      const OMap& at = Qp_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        int d = at.keys[r];
        if (d < 0) continue;
        add_branch(j, d, at.vals[r], true, 0.0);
      }
    }
    sec_branch += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_br).count();
    if (dz) {
      const auto t_fin = std::chrono::steady_clock::now();
      run_fin(&Qm, Qm_at[j], rsQm[j & 1], surv, j);
      run_fin(&Qm2, Qm2_at[j], rsQm2[j & 1], surv2, j);
      sec_qm += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_fin).count();
    }
  }

  void set_Qmp(int i, int j, double val) {
    if (val < HALF_NEG) return;
    if (dz) {
      rsQmp[j & 1].P.add(i, val);
      return;
    }
    double t = log_add(static_cast<double>(Qmp[ix(i, j)]), val);
    Qmp[ix(i, j)] = t;
    Qmp_at[j].set(i, t);
  }

  void add_branch_p(int j, int d, double obj, bool pk, double w_at) {
    const double w2 = E.w_bp;
    const double w3 = E.w_cp;
    double w_single = obj + ((pk ? E.w_beta1m + E.ml_pk_k * w2 : w2) + w_at);
    const int i_lo = lead_qm_lo(d);
    if (dz) {
      for (int i = i_lo; i <= d; ++i) set_Qmp(i, j, w_single + w3 * (d - i) + E.run(i, d - 1));
    } else {
      for (int i = i_lo; i <= d; ++i) set_Qmp(i, j, w_single + w3 * (d - i));
    }
    if (d == 0) return;
    const OMap& at = Qmp_at[d - 1];
    for (size_t r = 0; r < at.slots(); ++r) {
      int i = at.keys[r];
      if (i < 0) continue;
      set_Qmp(i, j, at.vals[r] + w_single);
    }
  }

  // fill_Qm with the stem-multiloop penalties b_p, c_p (no M2 split needed).
  void fill_Qmp(int j) {
    const double w3 = E.w_cp;
    std::vector<int> surv;
    if (dz) {
      RunSt& cu = rsQmp[j & 1];
      cu.clear();
      if (j > 0) {
        surv = live_keys(Qmp_at[j - 1]);
        run_ext(rsQmp[(j - 1) & 1], cu, surv, j, w3);
      }
    } else if (j > 0) {
      const OMap& at = Qmp_at[j - 1];
      for (size_t r = 0; r < at.slots(); ++r) {
        int i = at.keys[r];
        if (i < 0) continue;
        set_Qmp(i, j, at.vals[r] + w3);
      }
    }
    {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        int d = at.keys[r];
        if (d < 0 || !E.end_ok(d, j)) continue;
        add_branch_p(j, d, at.vals[r], false, E.at(d, j));
      }
    }
    const OMap& at = Qp_at[j];
    for (size_t r = 0; r < at.slots(); ++r) {
      int d = at.keys[r];
      if (d < 0) continue;
      add_branch_p(j, d, at.vals[r], true, 0.0);
    }
    if (dz) run_fin(&Qmp, Qmp_at[j], rsQmp[j & 1], surv, j);
  }

  void fill_Qz(int j) {
    const int pv = (j + 1) & 1, cu = j & 1;  // (j - 1) & 1 without a negative j - 1
    const double d5 = dz && j > 0 ? E.d5w[j - 1] : 0.0;
    const double s1 = dz ? E.s1w[j] : 0.0, d3 = dz ? E.d3w[j] : 0.0;
    for (int i = 0; i <= j; ++i) {
      double acc, u1 = NEG_INF, u2 = NEG_INF;
      if (dz) {
        if (i == j) {
          u1 = E.wu;  // empty region before j
        } else {
          double p0 = zP[pv][i];
          u1 = p0 > HALF_NEG ? p0 + E.wu : NEG_INF;
          double t = log_add(zU2[pv][i], zU1[pv][i] + d5);
          u2 = t > HALF_NEG ? t + E.wu : NEG_INF;
        }
        acc = NEG_INF;
      } else {
        acc = (i == j) ? E.wu : Qz[ix(i, j - 1)] + E.wu;
      }
      const OMap& atb = Qb_at[j];
      for (size_t r = 0; r < atb.slots(); ++r) {
        int d = atb.keys[r];
        if (d < 0 || d < i || !E.end_ok(d, j)) continue;
        double left = (d == i) ? 0.0 : Qz[ix(i, d - 1)];
        if (left < HALF_NEG) continue;
        acc = log_add(acc, left + atb.vals[r] + (E.w_nest + E.at(d, j)));
      }
      if (allow_pk) {
        const OMap& atp = Qp_at[j];
        for (size_t r = 0; r < atp.slots(); ++r) {
          int d = atp.keys[r];
          if (d < 0 || d < i) continue;
          double left = (d == i) ? 0.0 : Qz[ix(i, d - 1)];
          if (left < HALF_NEG) continue;
          acc = log_add(acc, left + atp.vals[r] + E.w_beta1p + E.nest_pk_k * E.w_nest);
        }
      }
      if (dz) {
        zP[cu][i] = acc;
        zU1[cu][i] = u1;
        zU2[cu][i] = u2;
        acc = run_close(acc, u1, u2, s1, d3);
      }
      Qz[ix(i, j)] = acc;
    }
  }

  // Every length. Top b2 starts at this right end, by C(k-1) + Q^z(k, j)
  // (C_hat on the scan, C_pre on the pre-pass). Dropped starts stay absent.
  void fill_Qz_kept(int j) {
    if (qz_keep.empty()) qz_keep.assign(static_cast<size_t>(n), {});
    const int B = p2_qz_beam;
    const double d5 = dz && j > 0 ? E.d5w[j - 1] : 0.0;
    const double s1 = dz ? E.s1w[j] : 0.0, d3 = dz ? E.d3w[j] : 0.0;
    struct Acc {
      double acc, u1, u2;
    };
    std::vector<int> seen(static_cast<size_t>(j + 1), 0), touch;
    std::vector<Acc> cell(static_cast<size_t>(j + 1));
    auto cell_at = [&](int k) -> Acc& {
      if (!seen[static_cast<size_t>(k)]) {
        seen[static_cast<size_t>(k)] = 1;
        cell[static_cast<size_t>(k)] = {NEG_INF, NEG_INF, NEG_INF};
        touch.push_back(k);
      }
      return cell[static_cast<size_t>(k)];
    };
    if (dz) cell_at(j).u1 = E.wu;
    else cell_at(j).acc = E.wu;
    if (j > 0) {
      for (const QkEnt& e : qz_keep[static_cast<size_t>(j - 1)]) {
        Acc& a = cell_at(e.k);
        if (dz) {
          a.u1 = e.p > HALF_NEG ? e.p + E.wu : NEG_INF;
          const double t = log_add(e.u2, e.u1 > HALF_NEG ? e.u1 + d5 : NEG_INF);
          a.u2 = t > HALF_NEG ? t + E.wu : NEG_INF;
        } else {
          a.acc = e.q + E.wu;
        }
      }
    }
    auto add_split = [&](int d, double w) {
      if (d < 0 || w < HALF_NEG) return;
      Acc& self = cell_at(d);
      self.acc = log_add(self.acc, w);
      if (d <= 0) return;
      for (const QkEnt& e : qz_keep[static_cast<size_t>(d - 1)]) {
        if (e.q < HALF_NEG || e.k > j) continue;
        Acc& a = cell_at(e.k);
        a.acc = log_add(a.acc, e.q + w);
      }
    };
    const OMap& atb = Qb_at[j];
    for (size_t r = 0; r < atb.slots(); ++r) {
      const int d = atb.keys[r];
      if (d < 0 || !E.end_ok(d, j)) continue;
      add_split(d, atb.vals[r] + (E.w_nest + E.at(d, j)));
    }
    if (allow_pk) {
      const OMap& atp = Qp_at[j];
      for (size_t r = 0; r < atp.slots(); ++r) {
        const int d = atp.keys[r];
        if (d < 0) continue;
        add_split(d, atp.vals[r] + E.w_beta1p + E.nest_pk_k * E.w_nest);
      }
    }
    struct Cand {
      int k;
      double q, p, u1, u2, key;
    };
    std::vector<Cand> cs;
    cs.reserve(touch.size());
    for (int k : touch) {
      const Acc& a = cell[static_cast<size_t>(k)];
      const double q = dz ? run_close(a.acc, a.u1, a.u2, s1, d3) : a.acc;
      if (q < HALF_NEG) continue;
      const double key = q + (k > 0 ? C[static_cast<size_t>(k - 1)] : 0.0);
      cs.push_back({k, q, dz ? a.acc : 0.0, a.u1, a.u2, key});
    }
    qz_cand += static_cast<long long>(cs.size());
    if (g_beam_stats.f) {
      std::vector<double> keys;
      keys.reserve(cs.size());
      for (const Cand& c : cs) keys.push_back(c.key);
      std::sort(keys.begin(), keys.end(), [](double a, double b) { return a > b; });
      const size_t kept_n =
          B > 0 && static_cast<int>(cs.size()) > B ? static_cast<size_t>(B) : keys.size();
      const int full = B > 0 && static_cast<int>(cs.size()) > B ? 1 : 0;
      g_beam_stats.log_prune("Qz", j, keys.size(), kept_n, BeamStats::mass_share(keys, kept_n),
                             full);
    }
    if (B > 0 && static_cast<int>(cs.size()) > B) {
      std::nth_element(cs.begin(), cs.begin() + B, cs.end(),
                       [](const Cand& a, const Cand& b) { return a.key > b.key; });
      cs.resize(static_cast<size_t>(B));
    }
    qz_kept += static_cast<long long>(cs.size());
    std::sort(cs.begin(), cs.end(), [](const Cand& a, const Cand& b) { return a.k < b.k; });
    std::vector<QkEnt>& col = qz_keep[static_cast<size_t>(j)];
    col.clear();
    col.reserve(cs.size());
    for (const Cand& c : cs) {
      Qz[ix(c.k, j)] = c.q;
      col.push_back({c.k, c.q, c.p, c.u1, c.u2});
    }
  }

  void fill_C(int j) {
    double acc, u1 = NEG_INF, u2 = NEG_INF;
    if (dz) {
      u1 = cP;
      u2 = j > 0 ? log_add(cU2, cU1 + E.d5w[j - 1]) : NEG_INF;
      acc = NEG_INF;
    } else {
      acc = (j == 0) ? 0.0 : C[j - 1];
    }
    const OMap& atb = Qb_at[j];
    for (size_t r = 0; r < atb.slots(); ++r) {
      int i = atb.keys[r];
      if (i < 0 || !E.end_ok(i, j)) continue;
      ++r_extc;
      double left = (i == 0) ? 0.0 : C[i - 1];
      if (left > HALF_NEG) acc = log_add(acc, left + atb.vals[r] + E.at(i, j));
    }
    if (allow_pk) {
      const OMap& atp = Qp_at[j];
      for (size_t r = 0; r < atp.slots(); ++r) {
        int i = atp.keys[r];
        if (i < 0) continue;
        double left = (i == 0) ? 0.0 : C[i - 1];
        if (left > HALF_NEG) acc = log_add(acc, left + atp.vals[r] + E.w_beta1);
      }
    }
    if (dz) {
      cP = acc;
      cU1 = u1;
      cU2 = u2;
      acc = run_close(acc, u1, u2, E.s1w[j], E.d3w[j]);
    }
    C[j] = acc;
    if (c_trace) {
      double qb = NEG_INF, qp = NEG_INF;
      const OMap& atb = Qb_at[j];
      for (size_t r = 0; r < atb.slots(); ++r) {
        int i = atb.keys[r];
        if (i < 0 || !E.end_ok(i, j)) continue;
        double left = (i == 0) ? 0.0 : C[i - 1];
        if (left > HALF_NEG) qb = log_add(qb, left + atb.vals[r] + E.at(i, j));
      }
      if (allow_pk) {
        const OMap& atp = Qp_at[j];
        for (size_t r = 0; r < atp.slots(); ++r) {
          int i = atp.keys[r];
          if (i < 0) continue;
          double left = (i == 0) ? 0.0 : C[i - 1];
          if (left > HALF_NEG) qp = log_add(qp, left + atp.vals[r] + E.w_beta1);
        }
      }
      std::fprintf(stderr, "[c] %d %.17g %.17g %.17g\n", j, C[j], qb, qp);
    }
  }

  // ============================================================= outside
  template <class T>
  void add2b(T& t, int i, int j, double val) {
    if (val < HALF_NEG || i < 0 || j < 0) return;
    const size_t k = ix(i, j);
    t[k] = log_add(static_cast<double>(t[k]), val);
  }

  void alloc_beta() {
    Qb_b.init(n);
    Qm_b.init(n);
    Qm2_b.init(n);
    if (E.stem_ml) Qmp_b.init(n);
    Qz_b.init(n);
    Qp_b.init(n);
    if (allow_pk) PKM_b.init(n);
    C_b.assign(n, NEG_INF);
    C_b[n - 1] = 0.0;
    if (dz) {
      for (int k = 0; k < 2; ++k) {
        rbQm[k].init(n);
        rbQm2[k].init(n);
        rbQmp[k].init(n);
        rbPKM[k].init(n);
        bz1[k].assign(n, NEG_INF);
        bz2[k].assign(n, NEG_INF);
      }
      bc1 = bc2 = NEG_INF;
    }
  }

  void outside() {
    if (!inside_done) pfunc();
    if (n == 0 || outside_done) return;
    outside_done = true;
    alloc_beta();
    if (lsp_rec)
      for (auto& v : QpS) v.assign(static_cast<size_t>(n) * n, NEG_INF);
    if (keep_closed_b && dz) {
      QmC_b.init(n);
      Qm2C_b.init(n);
    }
    for (int j = n - 1; j >= 0; --j) {
      if (allow_pk && l2r_mode == 1) l2r_back_end(j);
      outside_from_future(j);
      outside_same_j(j);
    }
  }

  void outside_from_future(int j) {
    const double w3 = E.w_alpha3;
    const double wu = E.wu;

    // PKM[j] = PKM[j-1] + α3 and + Qp[d][j]; its β is complete (Qb at j+1, PKM[j+1])
    if (allow_pk && dz) {
      const RunSt& nx = rbPKM[(j + 1) & 1];
      RunSt& cu = rbPKM[j & 1];
      cu.clear();
      const OMap& pkm = PKM_at[j];
      for (size_t r = 0; r < pkm.slots(); ++r) {
        int d = pkm.keys[r];
        if (d < 0) continue;
        double bP, b1, b2;
        run_back(PKM_b[ix(j, d)], nx.U1.get(d), nx.U2.get(d), j, w3, bP, b1, b2);
        cu.U1.put(d, b1);
        cu.U2.put(d, b2);
        if (keep_runs) run_weight(RS_PKM, d, j, bP, b1, b2);
        if (Qp_at[j].has(d)) {
          add2b(Qp_b, d, j, bP);
          qps(1, d, j, bP);
        }
      }
    } else if (allow_pk) {
      const OMap& pkm = PKM_at[j];
      for (size_t r = 0; r < pkm.slots(); ++r) {
        int d = pkm.keys[r];
        if (d < 0) continue;
        double pb = PKM_b[ix(j, d)];
        if (pb < HALF_NEG) continue;
        if (Qp_at[j].has(d)) {
          add2b(Qp_b, d, j, pb);
          qps(1, d, j, pb);
        }
        if (j > 0 && PKM_at[j - 1].has(d))
          PKM_b[ix(j - 1, d)] = log_add(static_cast<double>(PKM_b[ix(j - 1, d)]), pb + w3);
      }
    }

    if (dz) {
      double b1, b2;
      run_back(C_b[j], bc1, bc2, j, 0.0, C_b[j], b1, b2);
      bc1 = b1;
      bc2 = b2;
      const int nx = (j + 1) & 1, cu = j & 1;
      for (int i = 0; i <= j; ++i) {
        bz1[cu][i] = bz2[cu][i] = NEG_INF;
        if (Qz[ix(i, j)] > HALF_NEG) {
          double bP;
          run_back(Qz_b[ix(i, j)], bz1[nx][i], bz2[nx][i], j, wu, bP, bz1[cu][i], bz2[cu][i]);
          Qz_b[ix(i, j)] = bP;
        }
      }
    } else if (j + 1 < n) {
      C_b[j] = log_add(C_b[j], C_b[j + 1]);
      const OMap& am = Qm_at[j];
      for (size_t r = 0; r < am.slots(); ++r) {
        int i = am.keys[r];
        if (i < 0) continue;
        add2b(Qm_b, i, j, Qm_b[ix(i, j + 1)] + w3);
      }
      const OMap& am2 = Qm2_at[j];
      for (size_t r = 0; r < am2.slots(); ++r) {
        int i = am2.keys[r];
        if (i < 0) continue;
        add2b(Qm2_b, i, j, Qm2_b[ix(i, j + 1)] + w3);
      }
      if (allow_pk && E.stem_ml) {
        const OMap& amp = Qmp_at[j];
        for (size_t r = 0; r < amp.slots(); ++r) {
          int i = amp.keys[r];
          if (i < 0) continue;
          add2b(Qmp_b, i, j, Qmp_b[ix(i, j + 1)] + E.w_cp);
        }
      }
      for (int i = 0; i <= j; ++i)
        if (Qz[ix(i, j)] > HALF_NEG) add2b(Qz_b, i, j, Qz_b[ix(i, j + 1)] + wu);
    }

    // Qb interior wrap
    int jmax = std::min(n, j + 2 + cap);  // outer q with q - j - 1 <= cap
    {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        int d = at.keys[r];
        if (d < 0) continue;
        if (l2r_mode == 0) {
          const int p_lo = std::max(0, d - 1 - cap);
          for (int p = d - 1; p >= p_lo; --p) {
            const int budget = cap - (d - p - 1);
            if (budget < 0) break;
            for (int q = next_mate(p, j); q != -1 && q - j - 1 <= budget; q = next_mate(p, q)) {
              if (!Qb_at[q].has(p)) continue;
              const double w = E.interior(p, d, j, q);
              if (w < HALF_NEG) continue;
              add2b(Qb_b, d, j, Qb_b[ix(p, q)] + w);
            }
          }
          continue;
        }
        for (int q = j + 1; q < jmax; ++q) {
          const int budget = cap - (q - j - 1);
          const int i_hi = d - 1;
          const int i_lo = std::max(0, d - 1 - budget);
          if (i_lo > i_hi) continue;
          const std::vector<int>& pl = partners_left[q];
          auto it = std::lower_bound(pl.begin(), pl.end(), i_lo);
          const auto iend = std::upper_bound(it, pl.end(), i_hi);
          for (; it != iend; ++it) {
            const int i = *it;
            double w = E.interior(i, d, j, q);
            if (w < HALF_NEG) continue;
            add2b(Qb_b, d, j, Qb_b[ix(i, q)] + w);
          }
        }
      }
    }

    if (!allow_pk) return;

    // Qg (b)
    for (Grp& g : Qg[j].grps) {
      int e = g.v;
      for (Ent& s : g.ents) {
        int c = s.p, d = s.u;
        if (c == d) continue;
        for (int q = j + 1; q < jmax; ++q) {
          int i_lo = std::max(0, c - 1 - cap);
          for (int i = i_lo; i < c; ++i) {
            double w = E.stem_interior(i, c, j, q);
            if (w < HALF_NEG) continue;
            double pb = g4(Qg[q], e, i, d);
            add4(s, pb + w);
          }
        }
      }
    }

    // Qg (f)
    for (int q = j + 1; q < n; ++q) {
      int unpaired = q - j - 1;
      double w_ext = E.w_ap + E.w_bp + E.w_cp * unpaired;
      if (dz) w_ext += E.run(j + 1, q - 1);
      for (Grp& g : Qgls[j].grps) {
        int e = g.v;
        if (e == j && unpaired > cap) continue;
        for (Ent& s : g.ents) {
          int i = s.p - 1;
          if (!E.pairs(i, q) || !E.end_ok(i, q)) continue;
          double pb = g4(Qg[q], e, i, s.u);
          add4(s, pb + (w_ext + E.at(i, q)));
        }
      }
    }

    // Qg (g)(h)
    int q = j + 1;
    if (q < n) {
      for (Grp& g : Qgrs[j].grps) {
        int e = g.v;
        for (Ent& s : g.ents) {
          int c = s.p, d = s.u;
          int i_lo = lead_i_lo(c);
          for (int i = i_lo; i < c; ++i) {
            if (!E.pairs(i, q) || !E.end_ok(i, q)) continue;
            double w = E.w_ap + E.w_bp + E.w_cp * (c - i - 1) + E.at(i, q);
            if (dz) w += E.run(i + 1, c - 1);
            double pb = g4(Qg[q], e, i, d);
            add4(s, pb + w);
          }
          if (c - 1 >= 0) {
            const OMap& at = qmp_at(c - 1);
            qm_pair_slots(c - 1, q, [&](int r) {
              int ii = at.keys[r];
              if (ii < 0) return;
              double qm = at.vals[r];
              int outer = ii - 1;
              if (!E.pairs(outer, q) || outer >= c) return;
              if (!E.end_ok(outer, q)) return;
              double at_o = E.at(outer, q);
              double w = qm + E.w_ap + E.w_bp + at_o;
              double pb = g4(Qg[q], e, outer, d);
              add4(s, pb + w);
              add2b(qmp_b(), ii, c - 1, pb + s.a + E.w_ap + E.w_bp + at_o);
            });
          }
        }
      }
    }
  }

  void outside_same_j(int j) {
    const double w2 = E.w_alpha2;
    const double w3 = E.w_alpha3;
    const double a1 = E.w_alpha1;

    // C = C_left × Qb / Qp
    double cb = C_b[j];
    {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        int i = at.keys[r];
        if (i < 0 || !E.end_ok(i, j)) continue;
        double qb = at.vals[r];
        double at_ij = E.at(i, j);
        double left_a = (i == 0) ? 0.0 : C[i - 1];
        add2b(Qb_b, i, j, cb + left_a + at_ij);
        if (i > 0) C_b[i - 1] = log_add(C_b[i - 1], cb + qb + at_ij);
      }
    }
    if (allow_pk) {
      const OMap& at = Qp_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        int i = at.keys[r];
        if (i < 0) continue;
        double qp = at.vals[r];
        double left_a = (i == 0) ? 0.0 : C[i - 1];
        add2b(Qp_b, i, j, cb + left_a + E.w_beta1);
        qps(0, i, j, cb + left_a + E.w_beta1);
        if (i > 0) C_b[i - 1] = log_add(C_b[i - 1], cb + qp + E.w_beta1);
      }
    }

    // Qgrs = Qg[f] × Qm[f+1, j]
    if (allow_pk) {
      const OMap& at = qmp_at(j);
      SpMat& mpb = qmp_b();
      for (size_t r = 0; r < at.slots(); ++r) {
        int p = at.keys[r];
        if (p < 0) continue;
        double qm = at.vals[r];
        int f = p - 1;
        if (f < 0) continue;
        for (Grp& g : Qg[f].grps) {
          int e = g.v;
          for (Ent& s : g.ents) {
            if (s.p == s.u || !E.end_ok(s.p, f)) continue;
            double pb = g4(Qgrs[j], e, s.p, s.u);
            double w = E.w_bp + E.at(s.p, f);
            add4(s, pb + qm + w);
            add2b(mpb, p, j, pb + s.a + w);
          }
        }
        for (int c : partners_left[f]) {
          if (!E.end_ok(c, f)) continue;
          double pb = g4(Qgrs[j], f, c, c);
          add2b(mpb, p, j, pb + (E.w_bp + E.at(c, f)));
        }
      }
    }

    if (dz) {
      if (keep_closed_b) {
        for (size_t r = 0; r < Qm_at[j].slots(); ++r)
          if (Qm_at[j].keys[r] >= 0) QmC_b[ix(Qm_at[j].keys[r], j)] = Qm_b[ix(Qm_at[j].keys[r], j)];
        for (size_t r = 0; r < Qm2_at[j].slots(); ++r)
          if (Qm2_at[j].keys[r] >= 0)
            Qm2C_b[ix(Qm2_at[j].keys[r], j)] = Qm2_b[ix(Qm2_at[j].keys[r], j)];
      }
      run_back_at(Qm_b, Qm_at[j], j, w3, rbQm, RS_QM);
      run_back_at(Qm2_b, Qm2_at[j], j, w3, rbQm2, RS_QM2);
      if (allow_pk && E.stem_ml) run_back_at(Qmp_b, Qmp_at[j], j, E.w_cp, rbQmp, RS_QMP);
    }

    // Pass 1 E3 lands on the branch-ending outside, after the run split.
    lp_out_e3(j);

    // Qm / Qm2 from Qb, Qp at this j
    auto rev_branch = [&](int d, double obj_a, SpMat& obj_beta, bool pk,
                          double w_at) {
      double w_single_base = (pk ? E.w_beta1m + E.ml_pk_k * w2 : w2) + w_at;
      const OMap& amj = Qm_at[j];
      for (size_t r = 0; r < amj.slots(); ++r) {
        int i = amj.keys[r];
        if (i < 0 || i > d || d - i > E.ml_left || !lead_gap_ok(d - i)) continue;
        double w = w_single_base;
        if (l2r_mode == 0 && !allow_pk) {
          if (i != d) continue;
        } else {
          w += w3 * (d - i);
          if (dz) w += E.run(i, d - 1);
        }
        double mb = Qm_b[ix(i, j)];
        add2b(obj_beta, d, j, mb + w);
        if (pk) qps(1, d, j, mb + w);
      }
      if (d == 0) return;
      double w_more = w_single_base;
      const OMap& amd = Qm_at[d - 1];
      for (size_t r = 0; r < amd.slots(); ++r) {
        int i = amd.keys[r];
        if (i < 0) continue;
        double qm = amd.vals[r];
        double mb = Qm_b[ix(i, j)];
        add2b(obj_beta, d, j, mb + qm + w_more);
        if (pk) qps(1, d, j, mb + qm + w_more);
        add2b(Qm_b, i, d - 1, mb + obj_a + w_more);
        if (!pk) {
          double m2b = Qm2_b[ix(i, j)];
          add2b(obj_beta, d, j, m2b + qm + w_more);
          add2b(Qm_b, i, d - 1, m2b + obj_a + w_more);
        }
      }
    };
    {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        int d = at.keys[r];
        if (d < 0 || !E.end_ok(d, j)) continue;
        rev_branch(d, at.vals[r], Qb_b, false, E.at(d, j));
      }
    }
    if (allow_pk) {
      const OMap& at = Qp_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        int d = at.keys[r];
        if (d < 0) continue;
        rev_branch(d, at.vals[r], Qp_b, true, 0.0);
      }
    }

    // Qmp from Qb, Qp at this j
    if (allow_pk && E.stem_ml) {
      const double bp = E.w_bp, cp = E.w_cp;
      auto rev_branch_p = [&](int d, double obj_a, SpMat& obj_beta, bool pk,
                              double w_at) {
        double w_single = (pk ? E.w_beta1m + E.ml_pk_k * bp : bp) + w_at;
        const OMap& amj = Qmp_at[j];
        for (size_t r = 0; r < amj.slots(); ++r) {
          int i = amj.keys[r];
          if (i < 0 || i > d || !lead_gap_ok(d - i)) continue;
          double mb = Qmp_b[ix(i, j)];
          if (dz) {
            add2b(obj_beta, d, j, mb + w_single + cp * (d - i) + E.run(i, d - 1));
            if (pk) qps(1, d, j, mb + w_single + cp * (d - i) + E.run(i, d - 1));
            continue;
          }
          add2b(obj_beta, d, j, mb + w_single + cp * (d - i));
          if (pk) qps(1, d, j, mb + w_single + cp * (d - i));
        }
        if (d == 0) return;
        const OMap& amd = Qmp_at[d - 1];
        for (size_t r = 0; r < amd.slots(); ++r) {
          int i = amd.keys[r];
          if (i < 0) continue;
          double qm = amd.vals[r];
          double mb = Qmp_b[ix(i, j)];
          add2b(obj_beta, d, j, mb + qm + w_single);
          if (pk) qps(1, d, j, mb + qm + w_single);
          add2b(Qmp_b, i, d - 1, mb + obj_a + w_single);
        }
      };
      const OMap& atb = Qb_at[j];
      for (size_t r = 0; r < atb.slots(); ++r) {
        int d = atb.keys[r];
        if (d < 0 || !E.end_ok(d, j)) continue;
        rev_branch_p(d, atb.vals[r], Qb_b, false, E.at(d, j));
      }
      const OMap& atp = Qp_at[j];
      for (size_t r = 0; r < atp.slots(); ++r) {
        int d = atp.keys[r];
        if (d < 0) continue;
        rev_branch_p(d, atp.vals[r], Qp_b, true, 0.0);
      }
    }

    // Qz[i,j] += Qz_left × Qb/Qp
    for (int i = 0; i <= j; ++i) {
      double zb = Qz_b[ix(i, j)];
      if (zb < HALF_NEG) continue;
      const OMap& atb = Qb_at[j];
      for (size_t r = 0; r < atb.slots(); ++r) {
        int d = atb.keys[r];
        if (d < 0 || d < i || !E.end_ok(d, j)) continue;
        double qb = atb.vals[r];
        double w = E.w_nest + E.at(d, j);
        double left_a = (d == i) ? 0.0 : Qz[ix(i, d - 1)];
        add2b(Qb_b, d, j, zb + left_a + w);
        if (d > i) add2b(Qz_b, i, d - 1, zb + qb + w);
      }
      if (allow_pk) {
        const OMap& atp = Qp_at[j];
        for (size_t r = 0; r < atp.slots(); ++r) {
          int d = atp.keys[r];
          if (d < 0 || d < i) continue;
          double qp = atp.vals[r];
          double left_a = (d == i) ? 0.0 : Qz[ix(i, d - 1)];
          double w = E.w_beta1p + E.nest_pk_k * E.w_nest;
          add2b(Qp_b, d, j, zb + left_a + w);
          qps(2, d, j, zb + left_a + w);
          if (d > i) add2b(Qz_b, i, d - 1, zb + qp + w);
        }
      }
    }

    // Qb[i,j] += Qm2[i+1, j-1] × (α1+α2), plus the PK-inside terms from PKM[j-1]
    const double w_ml = a1 + E.w_alpha2;
    const double w_left = a1 + E.w_beta1m + (1.0 + E.ml_pk_k) * E.w_alpha2;
    {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        int i = at.keys[r];
        if (i < 0 || !E.end_ok(i, j)) continue;
        double at_ij = E.at_close(i, j);
        double bb = Qb_b[ix(i, j)];
        if (l2r_mode != 0 && i + 1 <= j - 1) add2b(Qm2_b, i + 1, j - 1, bb + (w_ml + at_ij));
        if (bb < HALF_NEG || !allow_pk || j < 1) continue;
        double wl = w_left + at_ij;

        const OMap& pkm = PKM_at[j - 1];
        for (size_t s = 0; s < pkm.slots(); ++s) {
          int d = pkm.keys[s];
          if (d <= i) continue;
          double pk = pkm.vals[s];
          const size_t pkix = ix(j - 1, d);
          double v = bb + wl + w3 * (d - i - 1);
          if (dz) v += E.run(i + 1, d - 1);
          PKM_b[pkix] = log_add(static_cast<double>(PKM_b[pkix]), v);
          if (i + 1 <= d - 1) {
            double left = Qm[ix(i + 1, d - 1)];
            if (left > HALF_NEG) {
              add2b(Qm_b, i + 1, d - 1, bb + pk + wl);
              PKM_b[pkix] = log_add(static_cast<double>(PKM_b[pkix]), bb + left + wl);
            }
          }
        }
      }
    }

    lp_out_multi(j);

    if (!allow_pk) return;
    if (l2r_mode == 1) l2r_back_pk(j);

    // Qp = Qgl × Qgr
    for (Grp& g : Qgr[j].grps) {
      int e_gr = g.v;
      int f = e_gr - 1;
      if (f < 0) continue;
      Tab4& gl_right = Qgl[f];
      for (Ent& s : g.ents) {
        int d = s.p, em1 = s.u;
        double gr = s.a;
        if (!E.end_ok(d, j)) continue;
        int e = em1 + 1;
        Grp* gl_group = gl_right.get(e);
        if (!gl_group) continue;
        for (Ent& t : gl_group->ents) {
          int i = t.p, uu = t.u;
          if (uu != d - 1 || i > d - 1 || !E.end_ok(i, f)) continue;
          double pb = Qp_b[ix(i, j)] + (E.at(i, f) + E.at(d, j));
          add4(s, pb + t.a);
          add4(t, pb + gr);
        }
      }
    }

    // Qgr = Qgl × Qz
    for (Grp& g : Qgl[j].grps) {
      int f = g.v;
      for (Ent& s : g.ents) {
        int i = s.p, d = s.u;
        double gl = s.a;
        for (int e = d + 1; e <= f; ++e) {
          double qz = qz_log(e, f - 1);
          if (qz < HALF_NEG) continue;
          double pb = g4(Qgr[j], e, i, d);
          add4(s, pb + qz);
          if (e <= f - 1) add2b(Qz_b, e, f - 1, pb + gl);
        }
      }
    }

    // Qgl = Qg × Qz × β2
    for (Grp& g : Qg[j].grps) {
      int f = g.v;
      for (Ent& s : g.ents) {
        int i = s.p, d = s.u;
        double qg = s.a;
        if (!E.end_ok(d, f)) continue;
        double w = E.w_stem + E.at(d, f);
        for (int e = d; e < f; ++e) {
          double qz = qz_log(d + 1, e);
          if (qz < HALF_NEG) continue;
          double pb = g4(Qgl[j], f, i, e);
          add4(s, pb + qz + w);
          if (d + 1 <= e) add2b(Qz_b, d + 1, e, pb + qg + w);
        }
      }
    }

    // Qgls = Qm × Qg
    SpMat& mpb = qmp_b();
    for (Grp& g : Qg[j].grps) {
      int e = g.v;
      for (Ent& s : g.ents) {
        int c = s.p, d = s.u;
        double qg = s.a;
        if (c - 1 < 0 || c == d || !E.end_ok(c, j)) continue;
        double w = E.w_bp + E.at(c, j);
        const OMap& at = qmp_at(c - 1);
        for (size_t r = 0; r < at.slots(); ++r) {
          int i = at.keys[r];
          if (i < 0) continue;
          double qm = at.vals[r];
          double pb = g4(Qgls[j], e, i, d);
          add4(s, pb + qm + w);
          add2b(mpb, i, c - 1, pb + qg + w);
        }
      }
    }
    for (int c : partners_left[j]) {
      if (c - 1 < 0 || !E.end_ok(c, j)) continue;
      double w = E.w_bp + E.at(c, j);
      const OMap& at = qmp_at(c - 1);
      for (size_t r = 0; r < at.slots(); ++r) {
        int i = at.keys[r];
        if (i < 0) continue;
        double pb = g4(Qgls[j], j, i, c);
        add2b(mpb, i, c - 1, pb + w);
      }
    }
  }

  // --------------------------------------------------------------- bpp
  std::vector<std::tuple<int, int, double>> collect_bpp(double cutoff) {
    std::vector<std::tuple<int, int, double>> out;
    if (n == 0 || logZ < HALF_NEG) return out;
    std::vector<double> logP(static_cast<size_t>(n) * n, NEG_INF);
    std::vector<char> seen(static_cast<size_t>(n) * n, 0);
    auto add_pair = [&](int i, int j, double lp) {
      if (lp < HALF_NEG || i < 0 || j <= i) return;
      size_t k = ix(i, j);
      logP[k] = log_add(seen[k] ? logP[k] : NEG_INF, lp);
      seen[k] = 1;
    };

    for (int j = 0; j < n; ++j) {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        int i = at.keys[r];
        if (i < 0) continue;
        add_pair(i, j, at.vals[r] + Qb_b[ix(i, j)] - logZ);
      }
    }

    if (allow_pk) {
      for (int j = 0; j < n; ++j)
        for (Grp& g : Qg[j].grps)
          for (Ent& s : g.ents) add_pair(s.p, j, s.a + s.b - logZ);
      collect_qg_inner_pairs(add_pair);
      if (l2r_mode == 1) l2r_pairs(add_pair);
    }

    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j) {
        size_t k = ix(i, j);
        if (!seen[k]) continue;
        double p = std::exp(std::min(0.0, logP[k]));
        if (p >= cutoff) out.emplace_back(i + 1, j + 1, std::min(1.0, p));
      }
    return out;
  }

  template <class AddPair>
  void collect_qg_inner_pairs(AddPair& add_pair) {
    for (int j = 0; j < n; ++j) {
      const auto& pair_i = partners_left[j];
      // (a)
      for (int i : pair_i) {
        int d_hi = std::min(j - 1, i + 1 + cap);
        for (int d = i + 1; d <= d_hi; ++d) {
          int e_lo = std::max(d + MIN_PAIR, j - 1 - cap);
          for (int e = e_lo; e < j; ++e) {
            double w = E.stem_interior(i, d, e, j);
            if (w < HALF_NEG) continue;
            add_pair(d, e, w + g4(Qg[j], e, i, d) - logZ);
          }
        }
      }
      // seed c·j with Qm on its left (Qgls group j)
      if (Grp* g = Qgls[j].get(j))
        for (const Ent& s : g->ents) add_pair(s.u, j, s.a + s.b - logZ);
      // seed c·e with Qm on its right (Qgrs entries (c, c))
      for (const Grp& g : Qgrs[j].grps)
        for (const Ent& s : g.ents)
          if (s.p == s.u) add_pair(s.p, g.v, s.a + s.b - logZ);
    }
  }
};

std::string json_num(double x) {
  char buf[40];
  std::snprintf(buf, sizeof buf, "%.17g", x);
  return buf;
}

template <class T>
std::string json_arr(const std::vector<T>& v) {
  std::string s = "[";
  for (size_t k = 0; k < v.size(); ++k) {
    if (k) s += ", ";
    if constexpr (std::is_floating_point<T>::value)
      s += json_num(v[k]);
    else
      s += std::to_string(v[k]);
  }
  return s + "]";
}

std::string json_str(const std::string& s) {
  std::string out = "\"";
  for (char c : s) {
    if (c == '"' || c == '\\') out += '\\';
    out += c;
  }
  return out + "\"";
}

// Per table (Qg, Qgl), over every stored state (run at infinite beam with the gate off for
// exact weights): w = exp(log inside + log outside - log Z), the state's posterior weight;
// "long" = outer pair p.q with q - p + 1 > long_span (a lower bound on the span of any
// pseudoknot the state is part of).  Per role (stem1: g = g1, both: g = max(g1, g2)) and
// tau: states passing g >= tau, their weight, their long weight, and how many of each
// (v, q) group's top-5 / top-20 states by w pass.
std::string DP03::gate_diag(const Gate& gd, const std::vector<double>& taus,
                            const std::string& id, int span, double delta, int long_span) {
  outside();
  const size_t T = taus.size();
  const int TOPS[2] = {5, 20};
  std::string js = "{\"id\": " + json_str(id) + ", \"n\": " + std::to_string(n) +
                   ", \"logZ\": " + json_num(logZ) + ", \"span\": " + std::to_string(span) +
                   ", \"delta\": " + json_num(delta) + ", \"long_span\": " +
                   std::to_string(long_span) + ", \"taus\": " + json_arr(taus) + ", \"tables\": {";
  const char* names[2] = {"Qg", "Qgl"};
  std::vector<Tab4>* tabs[2] = {&Qg, &Qgl};
  struct Item {
    double w, g[2];
    bool lg;
  };
  for (int ti = 0; ti < 2; ++ti) {
    long long states = 0, states_long = 0, groups = 0, top_total[2] = {0, 0};
    double w_all = 0.0, w_long = 0.0;
    std::vector<long long> pass_n[2], top_pass[2][2];
    std::vector<double> pass_w[2], pass_wl[2];
    for (int r = 0; r < 2; ++r) {
      pass_n[r].assign(T, 0);
      pass_w[r].assign(T, 0.0);
      pass_wl[r].assign(T, 0.0);
      top_pass[r][0].assign(T, 0);
      top_pass[r][1].assign(T, 0);
    }
    std::vector<Item> items;
    std::vector<size_t> order;
    for (int q = 0; q < n; ++q)
      for (const Grp& g : (*tabs[ti])[q].grps) {
        if (g.ents.empty()) continue;
        ++groups;
        items.clear();
        for (const Ent& e : g.ents) {
          Item it;
          it.w = (e.b < HALF_NEG || logZ < HALF_NEG) ? 0.0 : std::exp(e.a + e.b - logZ);
          double g1, g2;
          gd.roles(e.p, e.u, g.v, q, g1, g2);
          it.g[0] = g1;
          it.g[1] = std::max(g1, g2);
          it.lg = q - e.p + 1 > long_span;
          items.push_back(it);
          ++states;
          w_all += it.w;
          if (it.lg) {
            ++states_long;
            w_long += it.w;
          }
          for (int r = 0; r < 2; ++r)
            for (size_t k = 0; k < T; ++k)
              if (it.g[r] >= taus[k]) {
                ++pass_n[r][k];
                pass_w[r][k] += it.w;
                if (it.lg) pass_wl[r][k] += it.w;
              }
        }
        order.resize(items.size());
        for (size_t k = 0; k < order.size(); ++k) order[k] = k;
        std::stable_sort(order.begin(), order.end(),
                         [&](size_t a, size_t b) { return items[a].w > items[b].w; });
        for (int t = 0; t < 2; ++t) {
          size_t top = std::min(order.size(), static_cast<size_t>(TOPS[t]));
          top_total[t] += top;
          for (size_t m = 0; m < top; ++m)
            for (int r = 0; r < 2; ++r)
              for (size_t k = 0; k < T; ++k)
                if (items[order[m]].g[r] >= taus[k]) ++top_pass[r][t][k];
        }
      }
    if (ti) js += ", ";
    js += json_str(names[ti]) + ": {\"states\": " + std::to_string(states) +
          ", \"states_long\": " + std::to_string(states_long) + ", \"groups\": " +
          std::to_string(groups) + ", \"w\": " + json_num(w_all) + ", \"w_long\": " +
          json_num(w_long) + ", \"top5_total\": " + std::to_string(top_total[0]) +
          ", \"top20_total\": " + std::to_string(top_total[1]);
    const char* roles[2] = {"stem1", "both"};
    for (int r = 0; r < 2; ++r)
      js += ", " + json_str(roles[r]) + ": {\"pass_n\": " + json_arr(pass_n[r]) +
            ", \"pass_w\": " + json_arr(pass_w[r]) + ", \"pass_w_long\": " + json_arr(pass_wl[r]) +
            ", \"top5_pass\": " + json_arr(top_pass[r][0]) + ", \"top20_pass\": " +
            json_arr(top_pass[r][1]) + "}";
    js += "}";
  }
  return js + "}}";
}

// After outside: "Z id n logZ", then "S id table a b c d w" for every stored state of every
// beamed table with w = exp(log inside + log outside - log Z) >= min_w (keys as in
// state_key), then per table "T id table states N w W below_states K below_w V" (all
// states with w > 0; below = w < min_w, not written).  With dangles the Qm, Qm2, Qmp, PKM
// weight sums the cell's run states P, U1, U2 (needs keep_runs).
void DP03::dump_states(FILE* f, const std::string& id, double min_w) {
  outside();
  const char* sid = id.c_str();
  std::fprintf(f, "Z %s %d %.17g\n", sid, n, logZ);
  long long cnt[BT_N] = {}, below[BT_N] = {};
  double tot[BT_N] = {}, belw[BT_N] = {};
  auto emit = [&](int bt, int a, int b, int c, int d, double lw) {
    if (lw < HALF_NEG || logZ < HALF_NEG) return;
    const double w = std::exp(lw - logZ);
    ++cnt[bt];
    tot[bt] += w;
    if (w < min_w) {
      ++below[bt];
      belw[bt] += w;
      return;
    }
    std::fprintf(f, "S %s %s %d %d %d %d %.9e\n", sid, BT_NAMES[bt], a, b, c, d, w);
  };
  auto two = [&](int bt, const std::vector<OMap>& at, const SpMat& beta, int rs) {
    for (int j = 0; j < n; ++j)
      for (size_t r = 0; r < at[j].slots(); ++r) {
        int i = at[j].keys[r];
        if (i < 0) continue;
        double lw = (dz && rs >= 0) ? runW[rs][ix(i, j)] : at[j].vals[r] + beta[ix(i, j)];
        emit(bt, i, j, -1, -1, lw);
      }
  };
  two(BT_QB, Qb_at, Qb_b, -1);
  two(BT_QM, Qm_at, Qm_b, RS_QM);
  two(BT_QM2, Qm2_at, Qm2_b, RS_QM2);
  if (allow_pk && E.stem_ml) two(BT_QMP, Qmp_at, Qmp_b, RS_QMP);
  if (allow_pk) {
    two(BT_QP, Qp_at, Qp_b, -1);
    for (int j = 0; j < n; ++j)
      for (size_t r = 0; r < PKM_at[j].slots(); ++r) {
        int d = PKM_at[j].keys[r];
        if (d < 0) continue;
        double lw = dz ? runW[RS_PKM][ix(d, j)] : PKM_at[j].vals[r] + PKM_b[ix(j, d)];
        emit(BT_PKM, d, j, -1, -1, lw);
      }
    const std::pair<int, std::vector<Tab4>*> tabs[5] = {
        {BT_QG, &Qg}, {BT_QGLS, &Qgls}, {BT_QGL, &Qgl}, {BT_QGR, &Qgr}, {BT_QGRS, &Qgrs}};
    for (const auto& tb : tabs)
      for (int q = 0; q < n; ++q)
        for (const Grp& g : (*tb.second)[q].grps)
          for (const Ent& e : g.ents) emit(tb.first, e.p, e.u, g.v, q, e.a + e.b);
  }
  for (int t = 0; t < BT_N; ++t)
    std::fprintf(f, "T %s %s states %lld w %.9e below_states %lld below_w %.9e\n", sid,
                 BT_NAMES[t], cnt[t], tot[t], below[t], belw[t]);
  std::fflush(f);
}

// After outside of a step-mode run with every table at infinite beam: "Z id n logZ", then per
// table and step j with N > 0 stored states (the candidates of that pruning point) "E table
// j N K" and "C table j a b c d w r0 r1 cls" for each candidate with w >= min_w (w as in
// dump_states, keys as in state_key).  r0 = rank (0 = best) by the default pruning score
// (Qb, Qm, Qm2, Qmp, Qp, PKM: log C[i-1] + alpha, prune2's order; 4-index: four_score,
// keep_mask's order); r1 = rank by the --stem-est rank score of st with role all (Qg, Qgl,
// Qgr, partner table A built from the feasible Qgl states as in prune4_stem; infeasible
// states last; K = feasible count) and r1 = r0, K = N for the other tables.  cls (Qg, Qgl):
// stem-1 share of the outside from the tagged outside, 1 (>= 0.999) / 2 (<= 0.001) / 3
// (mixed) / 0 (w = 0); -1 for the other tables.  Qg / Qgl / Qgr lines then add "r2 r3 r4"
// (ranks as r1 of the combined score with (--stem-est-qg-r, --stem-est-consts) = (mirror,
// consistent), (partner, as-is), (partner, consistent)) and 10 tie-averaged ranks "a0 .. a9"
// of: the combined score (mirror as-is, mirror consistent, partner as-is, partner consistent),
// sigma_L (as-is, consistent), sigma_R (mirror as-is, mirror consistent, partner as-is,
// partner consistent); infeasible states and missing roles score -1e300.  Their E lines add
// the sum and the sum of squares of each a over all N candidates.
void DP03::dump_scores(FILE* f, const std::string& id, const StemEst& st, double min_w) {
  outside();
  std::fprintf(f, "Z %s %d %.17g\n", id.c_str(), n, logZ);
  std::vector<Tab4>* t4[5] = {&Qg, &Qgls, &Qgl, &Qgr, &Qgrs};
  const int b4[5] = {BT_QG, BT_QGLS, BT_QGL, BT_QGR, BT_QGRS};
  auto flat_b = [&](int t, std::vector<double>& out) {
    out.clear();
    for (Tab4& tb : *t4[t])
      for (Grp& g : tb.grps)
        for (Ent& e : g.ents) out.push_back(e.b);
  };
  std::vector<double> b0[5], b1;
  std::vector<signed char> cls[5];
  for (int t = 0; t < 5; ++t) flat_b(t, b0[t]);
  tagged_outside(1);
  for (int t : {0, 2}) {
    flat_b(t, b1);
    cls[t].resize(b1.size());
    for (size_t k = 0; k < b1.size(); ++k) {
      if (b0[t][k] < HALF_NEG) {
        cls[t][k] = 0;
        continue;
      }
      double share = b1[k] > HALF_NEG ? std::exp(b1[k] - b0[t][k]) : 0.0;
      cls[t][k] = share >= 0.999 ? 1 : share <= 0.001 ? 2 : 3;
    }
  }
  for (int t = 0; t < 5; ++t) {
    size_t k = 0;
    for (Tab4& tb : *t4[t])
      for (Grp& g : tb.grps)
        for (Ent& e : g.ents) e.b = b0[t][k++];
    b0[t].clear();
    b0[t].shrink_to_fit();
  }

  auto weight = [&](double lw) { return lw < HALF_NEG ? 0.0 : std::exp(lw - logZ); };
  auto two = [&](int bt, const std::vector<OMap>& at, const SpMat* beta, int rs,
                 bool pkm) {
    std::vector<std::pair<double, int>> ranked;
    std::vector<size_t> rank(n);
    for (int j = 0; j < n; ++j) {
      ranked.clear();
      for (size_t r = 0; r < at[j].slots(); ++r) {
        int i = at[j].keys[r];
        if (i >= 0) ranked.emplace_back(prefix_log(i) + at[j].vals[r], i);
      }
      if (ranked.empty()) continue;
      std::sort(ranked.begin(), ranked.end(),
                [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                  if (a.first != b.first) return a.first > b.first;
                  return a.second > b.second;
                });
      for (size_t k = 0; k < ranked.size(); ++k) rank[ranked[k].second] = k;
      std::fprintf(f, "E %s %d %zu %zu\n", BT_NAMES[bt], j, ranked.size(), ranked.size());
      for (size_t r = 0; r < at[j].slots(); ++r) {
        int i = at[j].keys[r];
        if (i < 0) continue;
        double lw = (dz && rs >= 0) ? runW[rs][ix(i, j)]
                                    : at[j].vals[r] + (*beta)[pkm ? ix(j, i) : ix(i, j)];
        double w = weight(lw);
        if (w < min_w || w == 0.0) continue;
        std::fprintf(f, "C %s %d %d %d -1 -1 %.9e %zu %zu -1\n", BT_NAMES[bt], j, i, j, w,
                     rank[i], rank[i]);
      }
    }
  };
  two(BT_QB, Qb_at, &Qb_b, -1, false);
  two(BT_QM, Qm_at, &Qm_b, RS_QM, false);
  two(BT_QM2, Qm2_at, &Qm2_b, RS_QM2, false);
  if (E.stem_ml) two(BT_QMP, Qmp_at, &Qmp_b, RS_QMP, false);
  two(BT_QP, Qp_at, &Qp_b, -1, false);
  two(BT_PKM, PKM_at, &PKM_b, RS_PKM, true);

  stA.clear();
  stE.clear();
  stc_q = -1;
  set_stem(&st, true, true);
  // Stem tables also get the ranks of the combined score under (qg_r, consts) = (mirror,
  // consistent), (partner, as-is), (partner, consistent) and the raw sigma_L (as-is,
  // consistent) and sigma_R (mirror as-is / consistent, partner as-is / consistent).
  static const bool VAR[3][2] = {{false, true}, {true, false}, {true, true}};
  // Tie-averaged ranks (0 = best) of the descending values.
  auto avg_ranks = [](const std::vector<double>& x) {
    std::vector<size_t> o(x.size());
    for (size_t k = 0; k < o.size(); ++k) o[k] = k;
    std::stable_sort(o.begin(), o.end(), [&](size_t a, size_t b) { return x[a] > x[b]; });
    std::vector<double> r(x.size());
    for (size_t a = 0; a < o.size();) {
      size_t b = a;
      while (b + 1 < o.size() && x[o[b + 1]] == x[o[a]]) ++b;
      for (size_t k = a; k <= b; ++k) r[o[k]] = 0.5 * static_cast<double>(a + b);
      a = b + 1;
    }
    return r;
  };
  size_t off[5] = {};
  std::vector<double> s0, s1, sv[3], raw[6];
  for (int q = 0; q < n; ++q)
    for (int t = 0; t < 5; ++t) {
      Tab4& tb = (*t4[t])[q];
      const int bt = b4[t];
      const bool stem = bt == BT_QG || bt == BT_QGL || bt == BT_QGR;
      std::unique_ptr<StemEval> ev;
      if (stem) ev.reset(new StemEval(*this, tb, bt, q));
      s0.clear();
      s1.clear();
      for (int x = 0; x < 3; ++x) sv[x].clear();
      for (int x = 0; x < 6; ++x) raw[x].clear();
      size_t k = 0, K = 0;
      for (const Grp& g : tb.grps)
        for (const Ent& e : g.ents) {
          s0.push_back(four_score(e.p, e.u, g.v, e.a));
          if (stem) {
            const bool okk = ev->ok[k];
            s1.push_back(okk ? ev->score(k) : NEG_INF);
            for (int x = 0; x < 3; ++x) {
              stem_qg_partner = VAR[x][0];
              stem_consistent = VAR[x][1];
              sv[x].push_back(okk ? ev->score(k) : NEG_INF);
            }
            for (int c = 0; c < 2; ++c) {
              stem_consistent = c == 1;
              stem_qg_partner = false;
              raw[c].push_back(okk ? ev->sigL(k) : NEG_INF);
              raw[2 + c].push_back(okk ? ev->sigR(k) : NEG_INF);
              stem_qg_partner = true;
              raw[4 + c].push_back(okk ? ev->sigR(k) : NEG_INF);
            }
            stem_qg_partner = stem_consistent = false;
            K += okk;
          }
          ++k;
        }
      const size_t N = k;
      if (bt == BT_QGL) stem_add_A(q, &ev->ok);
      if (N == 0) continue;
      std::vector<size_t> r0 = score_ranks(s0);
      std::vector<size_t> r1 = stem ? score_ranks(s1) : r0;
      std::vector<size_t> rv[3];
      std::vector<double> av[10];
      if (stem) {
        for (int x = 0; x < 3; ++x) rv[x] = score_ranks(sv[x]);
        av[0] = avg_ranks(s1);
        for (int x = 0; x < 3; ++x) av[1 + x] = avg_ranks(sv[x]);
        for (int x = 0; x < 6; ++x) av[4 + x] = avg_ranks(raw[x]);
      }
      std::fprintf(f, "E %s %d %zu %zu", BT_NAMES[bt], q, N, stem ? K : N);
      if (stem)
        for (int x = 0; x < 10; ++x) {
          double s = 0.0, s2 = 0.0;
          for (double r : av[x]) {
            s += r;
            s2 += r * r;
          }
          std::fprintf(f, " %.17g %.17g", s, s2);
        }
      std::fputc('\n', f);
      k = 0;
      for (const Grp& g : tb.grps)
        for (const Ent& e : g.ents) {
          double w = weight(e.a + e.b);
          if (w >= min_w && w > 0.0) {
            std::fprintf(f, "C %s %d %d %d %d %d %.9e %zu %zu %d", BT_NAMES[bt], q, e.p, e.u,
                         g.v, q, w, r0[k], r1[k], t == 0 || t == 2 ? cls[t][off[t] + k] : -1);
            if (stem) {
              std::fprintf(f, " %zu %zu %zu", rv[0][k], rv[1][k], rv[2][k]);
              for (int x = 0; x < 10; ++x) std::fprintf(f, " %.1f", av[x][k]);
            }
            std::fputc('\n', f);
          }
          ++k;
        }
      off[t] += N;
    }
  stem_ = nullptr;
  std::fflush(f);
}

// --dump-ext: after the inside pass of a step-mode run with every table at infinite beam (no
// outside), the Qg candidates of every step under the Result-18 scores, with the partner
// table A built from the feasible Qgl states as in prune4_stem (role all, qg-r partner,
// as-is).  "Z id n", then per step with N > 0 Qg candidates "E Qg j N K xs xb" (K feasible,
// xs / xb extension nodes visited, stack / bulge1) followed by the sum and the sum of squares
// over all N of 7 tie-averaged ranks (0 = best) of: partner (logsumexp(sigma_L, sigma_R)),
// partner + stack, partner + bulge1 (logsumexp with sigma_ext), sigma_L, sigma_R, sigma_ext
// stack, sigma_ext bulge1; infeasible states and missing terms score -1e300.  Then
// "C Qg j a b c d w rp rs rb ks kb a0 .. a6" for each candidate in w_of (the --dump-exact
// states, w as there): rp / rs / rb = keep_mask-order ranks of the three combined scores, ks /
// kb = bit m set if split selection (feasible states only) keeps it at b = 300, 1000, 3000
// (m = 0, 1, 2), stack / bulge1; a0 .. a6 as on the E line.
void DP03::dump_ext(FILE* f, const std::string& id, const StemEst& st,
                    const std::unordered_map<uint64_t, double>& w_of, int kmax) {
  stA.clear();
  stE.clear();
  stc_q = -1;
  set_stem(&st, true, true, true, false);
  stem_ext_max = kmax;
  std::fprintf(f, "Z %s %d\n", id.c_str(), n);
  auto avg_ranks = [](const std::vector<double>& x) {
    std::vector<size_t> o(x.size());
    for (size_t k = 0; k < o.size(); ++k) o[k] = k;
    std::stable_sort(o.begin(), o.end(), [&](size_t a, size_t b) { return x[a] > x[b]; });
    std::vector<double> r(x.size());
    for (size_t a = 0; a < o.size();) {
      size_t b = a;
      while (b + 1 < o.size() && x[o[b + 1]] == x[o[a]]) ++b;
      for (size_t k = a; k <= b; ++k) r[o[k]] = 0.5 * static_cast<double>(a + b);
      a = b + 1;
    }
    return r;
  };
  static const int BEAMS[3] = {300, 1000, 3000};
  std::vector<double> v[7];  // cp, cs, cb, sL, sR, xs, xb
  for (int q = 0; q < n; ++q) {
    {
      StemEval ev(*this, Qg[q], BT_QG, q);
      const size_t N = ev.ent.size();
      if (N > 0) {
        for (auto& x : v) x.assign(N, NEG_INF);
        size_t K = 0;
        for (size_t k = 0; k < N; ++k) {
          if (!ev.ok[k]) continue;
          ++K;
          v[3][k] = ev.sigL(k);
          v[4][k] = ev.sigR(k);
          v[0][k] = ev.score(k);
        }
        long long nodes[2];
        for (int m = 0; m < 2; ++m) {
          stem_ext = m + 1;
          stX.clear();
          const long long n0 = ext_nodes;
          for (size_t k = 0; k < N; ++k)
            if (ev.ok[k]) {
              const Ent& e = *ev.ent[k];
              v[5 + m][k] = stem_ext_score(e.p, e.u, ev.grp[k]->v, e.a);
            }
          nodes[m] = ext_nodes - n0;
        }
        stX.clear();
        stem_ext = 0;
        for (size_t k = 0; k < N; ++k)
          if (ev.ok[k]) {
            v[1][k] = log_add(v[0][k], v[5][k]);
            v[2][k] = log_add(v[0][k], v[6][k]);
          }
        std::vector<size_t> rk[3];
        for (int x = 0; x < 3; ++x) rk[x] = score_ranks(v[x]);
        std::vector<double> av[7];
        for (int x = 0; x < 7; ++x) av[x] = avg_ranks(v[x]);
        // split selection over the feasible states, as prune4_stem
        std::vector<unsigned char> kept[2];
        std::vector<size_t> idx;
        for (size_t k = 0; k < N; ++k)
          if (ev.ok[k]) idx.push_back(k);
        std::vector<double> sl, sr, sx, sc;
        std::vector<char> km;
        for (int m = 0; m < 2; ++m) {
          kept[m].assign(N, 0);
          sl.clear();
          sr.clear();
          sx.clear();
          sc.clear();
          for (size_t k : idx) {
            sl.push_back(v[3][k]);
            sr.push_back(v[4][k]);
            sx.push_back(v[5 + m][k]);
            sc.push_back(v[1 + m][k]);
          }
          const std::vector<double>* parts[3] = {&sl, &sr, &sx};
          for (int bi = 0; bi < 3; ++bi) {
            if (idx.size() <= static_cast<size_t>(BEAMS[bi])) {
              for (size_t k : idx) kept[m][k] |= 1 << bi;
              continue;
            }
            split_mask(parts, sc, km, BEAMS[bi]);
            for (size_t r = 0; r < idx.size(); ++r)
              if (km[r]) kept[m][idx[r]] |= 1 << bi;
          }
        }
        std::fprintf(f, "E Qg %d %zu %zu %lld %lld", q, N, K, nodes[0], nodes[1]);
        for (int x = 0; x < 7; ++x) {
          double s = 0.0, s2 = 0.0;
          for (double r : av[x]) {
            s += r;
            s2 += r * r;
          }
          std::fprintf(f, " %.17g %.17g", s, s2);
        }
        std::fputc('\n', f);
        for (size_t k = 0; k < N; ++k) {
          const Ent& e = *ev.ent[k];
          const int vv = ev.grp[k]->v;
          auto it = w_of.find(state_key(BT_QG, e.p, e.u, vv, q));
          if (it == w_of.end()) continue;
          std::fprintf(f, "C Qg %d %d %d %d %d %.9e %zu %zu %zu %d %d", q, e.p, e.u, vv, q,
                       it->second, rk[0][k], rk[1][k], rk[2][k], kept[0][k], kept[1][k]);
          for (int x = 0; x < 7; ++x) std::fprintf(f, " %.1f", av[x][k]);
          std::fputc('\n', f);
        }
      }
    }
    StemEval evl(*this, Qgl[q], BT_QGL, q);
    stem_add_A(q, &evl.ok);
  }
  stem_ = nullptr;
  std::fflush(f);
}

std::string normalize(const std::string& s) {
  std::string out;
  out.reserve(s.size());
  for (char c : s) {
    if (c == ' ' || c == '\t' || c == '\r' || c == '\n') continue;
    char u = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    out.push_back(u == 'T' ? 'U' : u);
  }
  return out;
}

void usage() {
  std::fprintf(stderr,
               "usage: dp03_cpp [--b1 N] [--b2 N] [--b4 N] [--bV N] [--bG N] "
               "[--g-cut topq|collapse|entries] [--g-sweep direct|states] "
               "[--loop-cap L] [--ml-lead-cap N|n] [--max-loop M] [--no-gu] [--rt RT] "
               "[--no-single-stems] [--wc-ends] "
               "[--score qz|prefix] [--no-pk] [--cutoff C] [--no-bpp] [--gap-beam group|step] "
               "[--params dp03|dp09|dp09st|dp09ml|dp09full] "
               "[--material dp03|rna1999|vienna2004 [--param-file PATH] [--dangles none|some|all|d2]] "
               "[--arm-score on|off --access-file PATH [--max-gap-scan N] [--max-arm-len W]] "
               "[--pair-bonus X] [--restrict-pairs] "
               "[--gate off|on [--gate-tau T] [--gate-span M] [--gate-delta D] "
               "[--gate-role stem1|both] [--gate-tables qg|qgl|qg,qgl] [--gate-reserve K]] "
               "[--gate-stats] [--gate-diag PATH [--gate-diag-taus T1,T2,...] "
               "[--gate-long-span L]] "
               "[--qzbar-old] [--pbeam-ratio a:b:c] "
               "[--beam-table T=B ...] [--dump-states PATH [--dump-min-w W]] "
               "[--dump-pruned PATH [--dump-exact PATH [--dump-pruned-min-w X]]] "
               "[--remove-state T,a,b,c,d ...] "
               "< sequences\n");
}

// "S id table a b c d w" lines of a --dump-states file, per id.
std::unordered_map<std::string, std::unordered_map<uint64_t, double>> load_states(
    const std::string& path) {
  std::unordered_map<std::string, std::unordered_map<uint64_t, double>> out;
  std::ifstream in(path);
  if (!in) {
    std::fprintf(stderr, "cannot read %s\n", path.c_str());
    std::exit(2);
  }
  std::string line, tag, id, tab;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] != 'S') continue;
    std::istringstream ls(line);
    int a, b, c, d;
    double w;
    ls >> tag >> id >> tab >> a >> b >> c >> d >> w;
    int bt = beam_table_id(tab);
    if (!ls || bt < 0) {
      std::fprintf(stderr, "%s: bad line: %s\n", path.c_str(), line.c_str());
      std::exit(2);
    }
    out[id][state_key(bt, a, b, c, d)] = w;
  }
  return out;
}

std::vector<double> parse_list(const std::string& s) {
  std::vector<double> out;
  size_t b = 0;
  while (b <= s.size()) {
    size_t e = s.find(',', b);
    if (e == std::string::npos) e = s.size();
    if (e > b) out.push_back(std::strtod(s.substr(b, e - b).c_str(), nullptr));
    b = e + 1;
  }
  return out;
}

}  // namespace

// Stem-role outside of the 4-index tables: the 4-index part of outside() with Qp_b fixed, where
// the Qp = Qgl x Qgr edge feeds only Qgl[f] (tag 1: the state ends up as the left stem, stem 1)
// or only Qgr[j] (tag 2: right stem).  Overwrites Ent::b of all 4-index tables.
void DP03::tagged_outside(int tag) {
  for (std::vector<Tab4>* tabs : {&Qg, &Qgls, &Qgl, &Qgr, &Qgrs})
    for (Tab4& t : *tabs)
      for (Grp& g : t.grps)
        for (Ent& e : g.ents) e.b = NEG_INF;
  for (int j = n - 1; j >= 0; --j) {
    const int jmax = std::min(n, j + 2 + cap);
    // Qg (b)
    for (Grp& g : Qg[j].grps) {
      int e = g.v;
      for (Ent& s : g.ents) {
        int c = s.p, d = s.u;
        if (c == d) continue;
        for (int q = j + 1; q < jmax; ++q) {
          int i_lo = std::max(0, c - 1 - cap);
          for (int i = i_lo; i < c; ++i) {
            double w = E.stem_interior(i, c, j, q);
            if (w < HALF_NEG) continue;
            add4(s, g4(Qg[q], e, i, d) + w);
          }
        }
      }
    }
    // Qg (f)
    for (int q = j + 1; q < n; ++q) {
      int unpaired = q - j - 1;
      double w_ext = E.w_ap + E.w_bp + E.w_cp * unpaired;
      if (dz) w_ext += E.run(j + 1, q - 1);
      for (Grp& g : Qgls[j].grps) {
        int e = g.v;
        if (e == j && unpaired > cap) continue;
        for (Ent& s : g.ents) {
          int i = s.p - 1;
          if (!E.pairs(i, q) || !E.end_ok(i, q)) continue;
          add4(s, g4(Qg[q], e, i, s.u) + (w_ext + E.at(i, q)));
        }
      }
    }
    // Qg (g)(h)
    if (j + 1 < n) {
      const int q = j + 1;
      for (Grp& g : Qgrs[j].grps) {
        int e = g.v;
        for (Ent& s : g.ents) {
          int c = s.p, d = s.u;
          for (int i = lead_i_lo(c); i < c; ++i) {
            if (!E.pairs(i, q) || !E.end_ok(i, q)) continue;
            double w = E.w_ap + E.w_bp + E.w_cp * (c - i - 1) + E.at(i, q);
            if (dz) w += E.run(i + 1, c - 1);
            add4(s, g4(Qg[q], e, i, d) + w);
          }
          if (c - 1 >= 0) {
            const OMap& at = qmp_at(c - 1);
            qm_pair_slots(c - 1, q, [&](int r) {
              int ii = at.keys[r];
              if (ii < 0) return;
              int outer = ii - 1;
              if (!E.pairs(outer, q) || outer >= c || !E.end_ok(outer, q)) return;
              double w = at.vals[r] + E.w_ap + E.w_bp + E.at(outer, q);
              add4(s, g4(Qg[q], e, outer, d) + w);
            });
          }
        }
      }
    }
    // Qgrs = Qg[f] x Qm[f+1, j]
    {
      const OMap& at = qmp_at(j);
      for (size_t r = 0; r < at.slots(); ++r) {
        int p = at.keys[r];
        if (p < 0 || p - 1 < 0) continue;
        double qm = at.vals[r];
        for (Grp& g : Qg[p - 1].grps) {
          int e = g.v;
          for (Ent& s : g.ents) {
            if (s.p == s.u || !E.end_ok(s.p, p - 1)) continue;
            add4(s, g4(Qgrs[j], e, s.p, s.u) + qm + (E.w_bp + E.at(s.p, p - 1)));
          }
        }
      }
    }
    // Qp = Qgl x Qgr
    for (Grp& g : Qgr[j].grps) {
      int f = g.v - 1;
      if (f < 0) continue;
      for (Ent& s : g.ents) {
        int d = s.p;
        if (!E.end_ok(d, j)) continue;
        Grp* gl_group = Qgl[f].get(s.u + 1);
        if (!gl_group) continue;
        for (Ent& t : gl_group->ents) {
          int i = t.p;
          if (t.u != d - 1 || i > d - 1 || !E.end_ok(i, f)) continue;
          double pb = Qp_b[ix(i, j)] + (E.at(i, f) + E.at(d, j));
          if (tag == 2) add4(s, pb + t.a);
          else add4(t, pb + s.a);
        }
      }
    }
    // Qgr = Qgl x Qz
    for (Grp& g : Qgl[j].grps) {
      int f = g.v;
      for (Ent& s : g.ents)
        for (int e = s.u + 1; e <= f; ++e) {
          double qz = qz_log(e, f - 1);
          if (qz < HALF_NEG) continue;
          add4(s, g4(Qgr[j], e, s.p, s.u) + qz);
        }
    }
    // Qgl = Qg x Qz x beta2
    for (Grp& g : Qg[j].grps) {
      int f = g.v;
      for (Ent& s : g.ents) {
        int d = s.u;
        if (!E.end_ok(d, f)) continue;
        double w = E.w_stem + E.at(d, f);
        for (int e = d; e < f; ++e) {
          double qz = qz_log(d + 1, e);
          if (qz < HALF_NEG) continue;
          add4(s, g4(Qgl[j], f, s.p, e) + qz + w);
        }
      }
    }
    // Qgls = Qm x Qg
    for (Grp& g : Qg[j].grps) {
      int e = g.v;
      for (Ent& s : g.ents) {
        int c = s.p, d = s.u;
        if (c - 1 < 0 || c == d || !E.end_ok(c, j)) continue;
        double w = E.w_bp + E.at(c, j);
        const OMap& at = qmp_at(c - 1);
        for (size_t r = 0; r < at.slots(); ++r) {
          int i = at.keys[r];
          if (i < 0) continue;
          add4(s, g4(Qgls[j], e, i, d) + at.vals[r] + w);
        }
      }
    }
  }
}

namespace stemdiag {

// Ranks (0-based, ties averaged) of x.
inline void avg_ranks(const std::vector<double>& x, std::vector<double>& r) {
  const size_t m = x.size();
  std::vector<size_t> o(m);
  for (size_t k = 0; k < m; ++k) o[k] = k;
  std::stable_sort(o.begin(), o.end(), [&](size_t a, size_t b) { return x[a] < x[b]; });
  r.assign(m, 0.0);
  for (size_t i = 0; i < m;) {
    size_t j = i;
    while (j + 1 < m && x[o[j + 1]] == x[o[i]]) ++j;
    double rk = 0.5 * static_cast<double>(i + j);
    for (size_t k = i; k <= j; ++k) r[o[k]] = rk;
    i = j + 1;
  }
}

// Pearson correlation; NaN if either side is constant.
inline double pearson(const std::vector<double>& a, const std::vector<double>& b) {
  const size_t m = a.size();
  double ma = 0, mb = 0;
  for (size_t k = 0; k < m; ++k) ma += a[k], mb += b[k];
  ma /= m;
  mb /= m;
  double sab = 0, saa = 0, sbb = 0;
  for (size_t k = 0; k < m; ++k) {
    sab += (a[k] - ma) * (b[k] - mb);
    saa += (a[k] - ma) * (a[k] - ma);
    sbb += (b[k] - mb) * (b[k] - mb);
  }
  if (saa <= 0 || sbb <= 0) return std::nan("");
  return sab / std::sqrt(saa * sbb);
}

// Top b of s (threshold rule, ties to the earliest), as DP03::keep_mask.
inline void top_b(const std::vector<double>& s, size_t b, std::vector<char>& keep) {
  keep.assign(s.size(), 1);
  if (s.size() <= b) return;
  std::vector<double> tmp(s);
  std::nth_element(tmp.begin(), tmp.begin() + (b - 1), tmp.end(), std::greater<double>());
  const double thr = tmp[b - 1];
  long room = static_cast<long>(b);
  for (double x : s) room -= x > thr;
  for (size_t k = 0; k < s.size(); ++k) {
    if (s[k] > thr) keep[k] = 1;
    else if (s[k] == thr && room > 0) keep[k] = 1, --room;
    else keep[k] = 0;
  }
}

inline void num(std::string& o, double x) {
  char buf[40];
  if (std::isnan(x)) std::snprintf(buf, sizeof buf, "null");
  else std::snprintf(buf, sizeof buf, "%.10g", x);
  o += buf;
}

}  // namespace stemdiag

// --pr-oracle, nested partner stems: fill_Qg on the kept pairs km (kl[j]: kept i < j) with
// interior loops <= lc per side, multiloop steps (f, g, h) if ml with the multiloop interiors
// mb (dense, both sides), seeds and run caps as the engine.
void DP03::pr_build_hbar(const std::vector<char>& km, const std::vector<std::vector<int>>& kl,
                         const std::vector<double>& mb, bool ml, int lc, std::vector<Tab4>& hb) {
  hb.assign(n, Tab4());
  std::vector<Tab4> ls(ml ? n : 0);
  auto M = [&](int i, int j) { return mb[ix(i, j)]; };
  for (int j = 0; j < n; ++j) {
    Tab4& G = hb[j];
    for (int i : kl[j]) {
      const int d_hi = std::min(j - 1, i + 1 + lc);
      for (int d = i + 1; d <= d_hi; ++d)
        for (int e = std::max(d + MIN_PAIR, j - 1 - lc); e < j; ++e)
          if (km[ix(d, e)]) G.acc(e, i, d, E.stem_interior(i, d, e, j) + 2.0 * pair_bonus);
    }
    for (int f = std::max(0, j - 1 - lc); f < j; ++f)
      for (const Grp& g : hb[f].grps)
        for (const Ent& s : g.ents) {
          if (s.p == s.u) continue;
          for (int i = std::max(0, s.p - 1 - lc); i < s.p; ++i)
            if (km[ix(i, j)]) G.acc(g.v, i, s.u, E.stem_interior(i, s.p, f, j) + s.a + pair_bonus);
        }
    if (ml) {
      for (int f = 0; f < j; ++f) {
        const int unpaired = j - f - 1;
        double w_ext = E.w_ap + E.w_bp + E.w_cp * unpaired + pair_bonus;
        if (dz) w_ext += E.run(f + 1, j - 1);
        for (const Grp& g : ls[f].grps) {
          if (g.v == f && unpaired > cap) continue;
          for (const Ent& s : g.ents) {
            const int i = s.p - 1;
            if (i >= 0 && km[ix(i, j)] && E.end_ok(i, j))
              G.acc(g.v, i, s.u, s.a + (w_ext + E.at(i, j)));
          }
        }
      }
      if (j >= 1) {
        Tab4 rs;
        const int jj = j - 1;
        for (int p = 1; p <= jj; ++p) {
          const double qm = M(p, jj);
          if (qm < HALF_NEG) continue;
          const int f = p - 1;
          for (const Grp& g : hb[f].grps)
            for (const Ent& s : g.ents)
              if (s.p != s.u && E.end_ok(s.p, f))
                rs.acc(g.v, s.p, s.u, s.a + qm + (E.w_bp + E.at(s.p, f)));
          for (int c : kl[f])
            if (E.end_ok(c, f)) rs.acc(f, c, c, qm + (E.w_bp + E.at(c, f)) + pair_bonus);
        }
        for (const Grp& g : rs.grps)
          for (const Ent& s : g.ents) {
            const int c = s.p, d = s.u;
            for (int i = lead_i_lo(c); i < c; ++i) {
              if (!km[ix(i, j)] || !E.end_ok(i, j)) continue;
              double w = s.a + E.w_ap + E.w_bp + E.w_cp * (c - i - 1) + E.at(i, j) + pair_bonus;
              if (dz) w += E.run(i + 1, c - 1);
              G.acc(g.v, i, d, w);
            }
            for (int i = 1; i < c; ++i) {
              const double qm = M(i, c - 1);
              if (qm < HALF_NEG) continue;
              const int o = i - 1;
              if (km[ix(o, j)] && E.end_ok(o, j))
                G.acc(g.v, o, d, qm + s.a + E.w_ap + E.w_bp + E.at(o, j) + pair_bonus);
            }
          }
      }
    }
    if (single_stems)
      for (int i : kl[j])
        if (E.end_ok(i, j)) G.acc(j, i, i, pair_bonus);
    if (ml) {
      Tab4& out = ls[j];
      for (const Grp& g : G.grps)
        for (const Ent& s : g.ents) {
          const int c = s.p;
          if (c - 1 < 1 || c == s.u || !E.end_ok(c, j)) continue;
          const double w = E.w_bp + E.at(c, j) + s.a;
          for (int i = 1; i <= c - 1; ++i) {
            const double qm = M(i, c - 1);
            if (qm > HALF_NEG) out.acc(g.v, i, s.u, w + qm);
          }
        }
      for (int c : kl[j]) {
        if (c - 1 < 1 || !E.end_ok(c, j)) continue;
        const double w = E.w_bp + E.at(c, j) + pair_bonus;
        for (int i = 1; i <= c - 1; ++i) {
          const double qm = M(i, c - 1);
          if (qm > HALF_NEG) out.acc(j, i, c, w + qm);
        }
      }
    }
  }
}

// --pr-oracle: the principled Out_hat as a 4-index outside pass (tagged_outside's order).
// Qg: growth onto kept pairs (interior <= lc per side; single-pair states too, as the
// inner=true roots), Qgls / Qgrs closings onto kept pairs (if ml), the stop through Qgl.
// Qgl's L role: the partner's Qgr from exact Qgl (hx) or from hb (nested H), times the future
// filler Qz (zx) or Qzbar; R role: the exact Qgl alpha of the partner.  Ctx: Qp_b (cx) or
// pr_ctx.  The right multiloop branch of a growth step: exact Qm (mx) or Mbar.  Ent.b holds
// the result; own[0] / own[1]: Qgls / Qgrs over every closing (their Out_hat), flat in
// off[1] / off[4] order.
void DP03::pr_hybrid_outside(const PrVar& v, const StemEst& st, const std::vector<char>& km,
                             const std::vector<std::vector<int>>& kr, std::vector<Tab4>* hb,
                             std::vector<double> own[2], const std::vector<size_t> off[5]) {
  for (std::vector<Tab4>* tabs : {&Qg, &Qgls, &Qgl, &Qgr, &Qgrs})
    for (Tab4& t : *tabs)
      for (Grp& g : t.grps)
        for (Ent& e : g.ents) e.b = NEG_INF;
  std::fill(own[0].begin(), own[0].end(), NEG_INF);
  std::fill(own[1].begin(), own[1].end(), NEG_INF);
  const int lc = std::min(v.lc, cap);
  auto ctx = [&](int i, int j) { return v.cx ? Qp_b[ix(i, j)] : pr_ctx(i, j); };
  Tab4 gl, gr;
  for (int j = n - 1; j >= 0; --j) {
    // Qg growth
    for (Grp& g : Qg[j].grps) {
      const int e = g.v;
      for (Ent& s : g.ents) {
        const int c = s.p, d = s.u;
        for (int i = std::max(0, c - 1 - lc); i < c; ++i) {
          const std::vector<int>& r = kr[i];
          for (auto it = std::upper_bound(r.begin(), r.end(), j);
               it != r.end() && *it <= j + 1 + lc; ++it) {
            const double w = E.stem_interior(i, c, j, *it);
            if (w < HALF_NEG) continue;
            add4(s, g4(Qg[*it], e, i, d) + w + pair_bonus);
          }
        }
      }
    }
    // Qgls closings (f)
    {
      size_t k = off[1][j];
      for (Grp& g : Qgls[j].grps) {
        const int e = g.v;
        for (Ent& s : g.ents) {
          double& ow = own[0][k++];
          const int i = s.p - 1;
          if (i < 0) continue;
          const std::vector<int>& pr = partners_right[i];
          for (auto it = std::upper_bound(pr.begin(), pr.end(), j); it != pr.end(); ++it) {
            const int q = *it, unpaired = q - j - 1;
            if (e == j && unpaired > cap) break;
            if (!E.end_ok(i, q)) continue;
            double w = E.w_ap + E.w_bp + E.w_cp * unpaired + pair_bonus + E.at(i, q);
            if (dz) w += E.run(j + 1, q - 1);
            const double t = g4(Qg[q], e, i, s.u);
            if (t < HALF_NEG) continue;
            ow = log_add(ow, t + w);
            if (v.ml && km[ix(i, q)]) add4(s, t + w);
          }
        }
      }
    }
    // Qgrs closings (g)(h) at j + 1
    if (j + 1 < n) {
      const int q = j + 1;
      size_t k = off[4][j];
      for (Grp& g : Qgrs[j].grps) {
        const int e = g.v;
        for (Ent& s : g.ents) {
          double& ow = own[1][k++];
          const int c = s.p, d = s.u;
          auto close = [&](int i, double w) {
            const double t = g4(Qg[q], e, i, d);
            if (t < HALF_NEG) return;
            ow = log_add(ow, t + w);
            if (v.ml && km[ix(i, q)]) add4(s, t + w);
          };
          for (int i = lead_i_lo(c); i < c; ++i) {
            if (!E.pairs(i, q) || !E.end_ok(i, q)) continue;
            double w = E.w_ap + E.w_bp + E.w_cp * (c - i - 1) + E.at(i, q) + pair_bonus;
            if (dz) w += E.run(i + 1, c - 1);
            close(i, w);
          }
          if (c - 1 >= 0) {
            const OMap& at = qmp_at(c - 1);
            for (size_t r = 0; r < at.slots(); ++r) {
              const int o = at.keys[r] - 1;
              if (at.keys[r] < 0 || o < 0 || o >= c || !E.pairs(o, q) || !E.end_ok(o, q)) continue;
              close(o, at.vals[r] + E.w_ap + E.w_bp + E.at(o, q) + pair_bonus);
            }
          }
        }
      }
    }
    // Qg[p-1] -> Qgrs[j] through the right branch [p, j]
    if (v.ml) {
      auto push = [&](int p, double qm) {
        const int f = p - 1;
        for (Grp& g : Qg[f].grps) {
          Grp* rg = Qgrs[j].get(g.v);
          if (!rg) continue;
          for (Ent& s : g.ents) {
            if (!E.end_ok(s.p, f)) continue;
            auto it = rg->idx.find(key2(s.p, s.u));
            if (it == rg->idx.end()) continue;
            add4(s, rg->ents[it->second].b + qm + (E.w_bp + E.at(s.p, f)));
          }
        }
      };
      if (v.mx) {
        const OMap& at = qmp_at(j);
        for (size_t r = 0; r < at.slots(); ++r)
          if (at.keys[r] >= 1) push(at.keys[r], at.vals[r]);
      } else {
        for (int p = 1; p <= j; ++p) {
          const double qm = st.PrMbar(p, j);
          if (qm > HALF_NEG) push(p, qm);
        }
      }
    }
    // Qp, R role: the partner L-stem's exact Qgl alpha
    for (Grp& g : Qgr[j].grps) {
      const int f = g.v - 1;
      if (f < 0) continue;
      for (Ent& s : g.ents) {
        const int d = s.p;
        if (!E.end_ok(d, j)) continue;
        Grp* yg = Qgl[f].get(s.u + 1);
        if (!yg) continue;
        for (const Ent& t : yg->ents) {
          if (t.u != d - 1 || t.p > d - 1 || !E.end_ok(t.p, f)) continue;
          add4(s, t.a + ctx(t.p, j) + (E.at(t.p, f) + E.at(d, j)));
        }
      }
    }
    // Qp, L role: the partner R-stem's (hybrid) Qgr alpha
    Tab4* GR = &Qgr[j];
    if (!(v.hx && v.zx)) {
      gl = Tab4();
      gr = Tab4();
      const Tab4* L = &Qgl[j];
      if (!v.hx) {
        for (const Grp& g : (*hb)[j].grps) {
          const int f = g.v;
          for (const Ent& s : g.ents) {
            const int d = s.u;
            if (!E.end_ok(d, f)) continue;
            const double w = E.w_stem + E.at(d, f) + s.a;
            for (int e = d; e < f; ++e) {
              const double qz = qz_log(d + 1, e);
              if (qz > HALF_NEG) gl.acc(f, s.p, e, w + qz);
            }
          }
        }
        L = &gl;
      }
      for (const Grp& g : L->grps) {
        const int f = g.v;
        for (const Ent& s : g.ents)
          for (int e = s.u + 1; e <= f; ++e) {
            const double z = v.zx ? qz_log(e, f - 1) : st.Qzbar(e, f - 1);
            if (z > HALF_NEG) gr.acc(e, s.p, s.u, s.a + z);
          }
      }
      GR = &gr;
    }
    for (const Grp& g : GR->grps) {
      const int f = g.v - 1;
      if (f < 0) continue;
      for (const Ent& s : g.ents) {
        const int d = s.p;
        if (!E.end_ok(d, j)) continue;
        Grp* yg = Qgl[f].get(s.u + 1);
        if (!yg) continue;
        for (Ent& t : yg->ents) {
          if (t.u != d - 1 || t.p > d - 1 || !E.end_ok(t.p, f)) continue;
          add4(t, s.a + ctx(t.p, j) + (E.at(t.p, f) + E.at(d, j)));
        }
      }
    }
    // Qgr = Qgl x Qz
    for (Grp& g : Qgl[j].grps) {
      const int f = g.v;
      for (Ent& s : g.ents)
        for (int e = s.u + 1; e <= f; ++e) {
          const double qz = qz_log(e, f - 1);
          if (qz < HALF_NEG) continue;
          add4(s, g4(Qgr[j], e, s.p, s.u) + qz);
        }
    }
    // Qgl = Qg x Qz x beta2
    for (Grp& g : Qg[j].grps) {
      const int f = g.v;
      for (Ent& s : g.ents) {
        const int d = s.u;
        if (!E.end_ok(d, f)) continue;
        const double w = E.w_stem + E.at(d, f);
        for (int e = d; e < f; ++e) {
          const double qz = qz_log(d + 1, e);
          if (qz < HALF_NEG) continue;
          add4(s, g4(Qgl[j], f, s.p, e) + qz + w);
        }
      }
    }
    // Qgls = Qm x Qg (single-pair states: their seeds)
    if (v.ml)
      for (Grp& g : Qg[j].grps) {
        Grp* lg = Qgls[j].get(g.v);
        if (!lg) continue;
        for (Ent& s : g.ents) {
          const int c = s.p;
          if (c - 1 < 0 || !E.end_ok(c, j)) continue;
          const double w = E.w_bp + E.at(c, j);
          const OMap& at = qmp_at(c - 1);
          for (size_t r = 0; r < at.slots(); ++r) {
            if (at.keys[r] < 0) continue;
            auto it = lg->idx.find(key2(at.keys[r], s.u));
            if (it != lg->idx.end()) add4(s, lg->ents[it->second].b + at.vals[r] + w);
          }
        }
      }
  }
}

// --pr-oracle PATH (--beam inf, --stem-est-principled on, after the outside pass).  Variant
// sets (comma list): step1 (every piece exact, kept = nested P > 0: ex_B ml on lc 30, ex_noml
// ml off lc 30, ex_A ml off lc 2), step2 (one piece exact, the others nested, floors 1e-7 and
// 0, B growth), step2A (the same, A growth), check (all nested on the run's kept pairs,
// pr_loop / pr_ml, against pr_score4; run with --pr-arms-k 0 --pr-cut 0 and large
// --pr-ml-span / --pr-ml-close).  For every state, score - alpha = Ent.b (Qg, Qgl, Qgr) or
// own (Qgls, Qgrs) against the exact beta; w = exp(alpha + beta - logZ).  Lines:
//   V id var hb_sec pass_sec stats_sec n_kept
//   R id var tab q N W kept1000 oracle1000 rho nboth sumd maxd wsumd nmiss wmiss nspur
//   D id var Qg2 nboth sumd maxd nmiss    (Qg states with two or more pairs)
//   C id tab N nmis maxd    (check)
void DP03::pr_oracle_diag(StemEst& st, FILE* f, const std::string& id, const std::string& sets) {
  outside();
  stem_ = &st;
  const char* sid = id.c_str();
  const char* TN[5] = {"Qg", "Qgls", "Qgl", "Qgr", "Qgrs"};
  std::vector<Tab4>* T5[5] = {&Qg, &Qgls, &Qgl, &Qgr, &Qgrs};
  std::vector<size_t> off[5];
  std::vector<double> bx[5];
  for (int t = 0; t < 5; ++t) {
    off[t].assign(n + 1, 0);
    for (int q = 0; q < n; ++q) {
      size_t c = 0;
      for (const Grp& g : (*T5[t])[q].grps) c += g.ents.size();
      off[t][q + 1] = off[t][q] + c;
    }
    bx[t].resize(off[t][n]);
    size_t k = 0;
    for (int q = 0; q < n; ++q)
      for (const Grp& g : (*T5[t])[q].grps)
        for (const Ent& e : g.ents) bx[t][k++] = e.b;
  }
  std::vector<double> own[2];
  own[0].assign(off[1][n], NEG_INF);
  own[1].assign(off[4][n], NEG_INF);
  const size_t nn = static_cast<size_t>(n) * n;
  std::vector<double> mbx(nn, NEG_INF);
  for (int j = 0; j < n; ++j) {
    const OMap& at = qmp_at(j);
    for (size_t r = 0; r < at.slots(); ++r)
      if (at.keys[r] >= 0) mbx[ix(at.keys[r], j)] = at.vals[r];
  }

  std::vector<PrVar> vars;
  bool check = false;
  auto has = [&](const char* s) {
    return ("," + sets + ",").find(std::string(",") + s + ",") != std::string::npos;
  };
  if (has("step1")) {
    vars.push_back({"ex_B", true, true, true, true, 1, true, 30});
    vars.push_back({"ex_noml", true, true, true, true, 1, false, 30});
    vars.push_back({"ex_A", true, true, true, true, 1, false, 2});
  }
  for (int a = 0; a < 2; ++a) {
    if (!has(a == 0 ? "step2" : "step2A")) continue;
    const bool ml = a == 0;
    const int lc = a == 0 ? 30 : 2;
    const std::string g = a == 0 ? "" : "A";
    for (int fi = 0; fi < 2; ++fi) {
      const std::string s = g + (fi == 0 ? "_f7" : "_f0");
      vars.push_back({"nest" + s, false, false, false, false, fi, ml, lc});
      vars.push_back({"Z" + s, false, true, false, false, fi, ml, lc});
      vars.push_back({"C" + s, false, false, false, true, fi, ml, lc});
      vars.push_back({"M" + s, false, false, true, false, fi, ml, lc});
      vars.push_back({"H" + s, true, false, false, false, fi, ml, lc});
    }
  }
  if (has("check")) {
    check = true;
    vars.push_back({"check", false, false, false, false, -1, pr_ml, pr_loop});
  }

  std::vector<char> km;
  std::vector<std::vector<int>> kr, kl;
  int km_fi = -2;
  auto set_mask = [&](int fi) {
    if (fi == km_fi) return;
    km_fi = fi;
    if (fi < 0) {
      km = st.prKept;
    } else {
      const double lf = fi == 0 ? std::log(1e-7) : NEG_INF;
      km.assign(nn, 0);
      for (size_t k = 0; k < nn; ++k)
        if (st.prLP[k] > HALF_NEG && !(st.prLP[k] < lf)) km[k] = 1;
    }
    kr.assign(n, {});
    kl.assign(n, {});
    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j)
        if (km[ix(i, j)]) kr[i].push_back(j), kl[j].push_back(i);
  };
  std::vector<Tab4> hb;
  std::tuple<int, bool, bool, int> hb_key{-9, false, false, 0};

  std::vector<double> sc, w, rs, rw;
  std::vector<char> keep;
  for (const PrVar& v : vars) {
    const auto t0 = std::chrono::steady_clock::now();
    set_mask(v.fi);
    long long nk = 0;
    for (char c : km) nk += c;
    if (!v.hx) {
      const std::tuple<int, bool, bool, int> key{v.fi, v.mx, v.ml, std::min(v.lc, cap)};
      if (key != hb_key) {
        hb_key = key;
        pr_build_hbar(km, kl, v.mx ? mbx : st.prMbar, v.ml, std::min(v.lc, cap), hb);
      }
    }
    const auto t1 = std::chrono::steady_clock::now();
    pr_hybrid_outside(v, st, km, kr, v.hx ? nullptr : &hb, own, off);
    const auto t2 = std::chrono::steady_clock::now();
    {
      long long nb = 0, nm = 0;
      double sd = 0, md = 0;
      size_t k = 0;
      for (int q = 0; q < n; ++q)
        for (const Grp& g : Qg[q].grps)
          for (const Ent& e : g.ents) {
            const double b0 = bx[0][k++];
            if (e.p == e.u || b0 < HALF_NEG) continue;
            if (e.b < HALF_NEG) {
              ++nm;
              continue;
            }
            const double d = std::fabs(e.b - b0);
            ++nb;
            sd += d;
            md = std::max(md, d);
          }
      std::fprintf(f, "D %s %s Qg2 %lld %.9e %.9e %lld\n", sid, v.name.c_str(), nb, sd, md, nm);
    }
    for (int t = 0; t < 5; ++t)
      for (int q = 0; q < n; ++q) {
        const size_t lo = off[t][q], N = off[t][q + 1] - lo;
        if (N == 0) continue;
        sc.resize(N);
        w.resize(N);
        double W = 0, sumd = 0, maxd = 0, wsumd = 0, wmiss = 0;
        long long nboth = 0, nmiss = 0, nspur = 0;
        size_t k = 0;
        for (const Grp& g : (*T5[t])[q].grps)
          for (const Ent& e : g.ents) {
            const double b0 = bx[t][lo + k];
            const double b1 = t == 1 ? own[0][lo + k] : t == 4 ? own[1][lo + k] : e.b;
            const double wk = (b0 > HALF_NEG && logZ > HALF_NEG) ? std::exp(e.a + b0 - logZ) : 0.0;
            w[k] = wk;
            W += wk;
            sc[k] = b1 > HALF_NEG ? e.a + b1 : NEG_INF;
            if (b0 > HALF_NEG) {
              if (b1 > HALF_NEG) {
                const double d = std::fabs(b1 - b0);
                ++nboth;
                sumd += d;
                wsumd += wk * d;
                maxd = std::max(maxd, d);
              } else {
                ++nmiss;
                wmiss += wk;
              }
            } else if (b1 > HALF_NEG) {
              ++nspur;
            }
            ++k;
          }
        double kept = W, orc = W, rho = std::nan("");
        if (N > 1000) {
          stemdiag::top_b(sc, 1000, keep);
          kept = 0;
          for (size_t m = 0; m < N; ++m)
            if (keep[m]) kept += w[m];
          stemdiag::top_b(w, 1000, keep);
          orc = 0;
          for (size_t m = 0; m < N; ++m)
            if (keep[m]) orc += w[m];
        }
        if (N >= 3) {
          stemdiag::avg_ranks(sc, rs);
          stemdiag::avg_ranks(w, rw);
          rho = stemdiag::pearson(rs, rw);
        }
        std::fprintf(f, "R %s %s %s %d %zu %.9e %.9e %.9e %.6f %lld %.9e %.9e %.9e %lld %.9e %lld\n",
                     sid, v.name.c_str(), TN[t], q, N, W, kept, orc, rho, nboth, sumd, maxd, wsumd,
                     nmiss, wmiss, nspur);
      }
    const auto t3 = std::chrono::steady_clock::now();
    auto sec = [](auto a, auto b) { return std::chrono::duration<double>(b - a).count(); };
    std::fprintf(f, "V %s %s %.2f %.2f %.2f %lld\n", sid, v.name.c_str(), sec(t0, t1), sec(t1, t2),
                 sec(t2, t3), nk);
    std::fflush(f);
    std::fprintf(stderr, "[pr-oracle] %s %s kept %lld hb %.1fs pass %.1fs stats %.1fs\n", sid,
                 v.name.c_str(), nk, sec(t0, t1), sec(t1, t2), sec(t2, t3));

    if (check && v.name == "check") {
      stA.clear();
      stE.clear();
      prKc.clear();
      prSc.clear();
      pr_sc_n = 0;
      stc_q = -1;
      const int BT[5] = {BT_QG, BT_QGLS, BT_QGL, BT_QGR, BT_QGRS};
      long long cn[5] = {}, cm[5] = {};
      double cd[5] = {};
      for (int q = 0; q < n; ++q) {
        stem_cache_step(q);
        for (int t = 0; t < 5; ++t) {
          size_t k = off[t][q];
          for (const Grp& g : (*T5[t])[q].grps)
            for (const Ent& e : g.ents) {
              const double b1 = t == 1 ? own[0][k] : t == 4 ? own[1][k] : e.b;
              ++k;
              const double s = pr_score4(BT[t], e, g.v, q);
              const double x = b1 > HALF_NEG ? e.a + b1 : NEG_INF;
              ++cn[t];
              if (s < HALF_NEG && x < HALF_NEG) continue;
              const double d = (s < HALF_NEG || x < HALF_NEG) ? INFINITY : std::fabs(s - x);
              if (d > 1e-9) ++cm[t];
              cd[t] = std::max(cd[t], d);
            }
        }
        stem_add_A(q);
      }
      for (int t = 0; t < 5; ++t) {
        std::fprintf(f, "C %s %s %lld %lld %.3e\n", sid, TN[t], cn[t], cm[t], cd[t]);
        std::fprintf(stderr, "[pr-oracle check] %s %s states %lld mismatch %lld max_absdiff %.3e\n",
                     sid, TN[t], cn[t], cm[t], cd[t]);
      }
    }
  }
  // score: pr_score4 as configured (--pr-abl etc.) on every state against its exact w, per
  // table and step (beams 100, 300, 1000; see the header)
  if (has("score")) {
    const auto t0 = std::chrono::steady_clock::now();
    stA.clear();
    stE.clear();
    prKc.clear();
    prSc.clear();
    pr_sc_n = 0;
    stc_q = -1;
    const int BT[5] = {BT_QG, BT_QGLS, BT_QGL, BT_QGR, BT_QGRS};
    for (int q = 0; q < n; ++q) {
      stem_cache_step(q);
      for (int t = 0; t < 5; ++t) {
        const size_t lo = off[t][q], N = off[t][q + 1] - lo;
        if (N == 0) continue;
        sc.resize(N);
        w.resize(N);
        double W = 0, wmiss = 0;
        long long nmiss = 0;
        size_t k = 0;
        for (const Grp& g : (*T5[t])[q].grps)
          for (const Ent& e : g.ents) {
            const double b0 = bx[t][lo + k];
            const double wk = (b0 > HALF_NEG && logZ > HALF_NEG) ? std::exp(e.a + b0 - logZ) : 0.0;
            const double s = pr_score4(BT[t], e, g.v, q);
            w[k] = wk;
            sc[k] = s;
            W += wk;
            if (s < HALF_NEG && wk > 0) ++nmiss, wmiss += wk;
            ++k;
          }
        double kb[3] = {W, W, W}, ob[3] = {W, W, W};
        const size_t BB[3] = {100, 300, 1000};
        for (int bi = 0; bi < 3; ++bi) {
          if (N <= BB[bi]) continue;
          stemdiag::top_b(sc, BB[bi], keep);
          kb[bi] = 0;
          for (size_t m = 0; m < N; ++m)
            if (keep[m]) kb[bi] += w[m];
          stemdiag::top_b(w, BB[bi], keep);
          ob[bi] = 0;
          for (size_t m = 0; m < N; ++m)
            if (keep[m]) ob[bi] += w[m];
        }
        double rho = std::nan("");
        if (N >= 3) {
          stemdiag::avg_ranks(sc, rs);
          stemdiag::avg_ranks(w, rw);
          rho = stemdiag::pearson(rs, rw);
        }
        std::fprintf(f, "A %s %s %d %zu %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.6f %lld %.9e\n", sid,
                     TN[t], q, N, W, kb[0], ob[0], kb[1], ob[1], kb[2], ob[2], rho, nmiss, wmiss);
      }
      stem_add_A(q);
    }
    std::fprintf(stderr, "[pr-oracle score] %s abl %d pk %d sep_evals %lld pp_mean %.3f ppk_mean %.3f %.1fs\n",
                 sid, pr_abl, pr_abl_pk, pr_sep_evals,
                 pr_sep_evals ? static_cast<double>(pr_sep_pp) / pr_sep_evals : 0.0,
                 pr_sep_evals ? static_cast<double>(pr_sep_ppk) / pr_sep_evals : 0.0,
                 std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
  }
  for (int t = 0; t < 5; ++t) {
    size_t k = 0;
    for (int q = 0; q < n; ++q)
      for (Grp& g : (*T5[t])[q].grps)
        for (Ent& e : g.ents) e.b = bx[t][k++];
  }
}

// --stem-diag (run at --beam inf after the outside pass): for every Qg / Qgl state, the exact
// posterior weight w = inside x outside / Z, its role (stem1: >= 99.9% of w through the left-stem
// Qp edge, stem2: <= 0.1%, mixed), and the pruning scores
//   s1 = log C[p-1] + log X + log Qz[u+1][v-1]                         (current)
//   s1H_L = s1 + log H_L(u, v, q) - log C0(q) + beta1 + 2 beta2      (stem 1, loops <= L)
//   s2 = log B(p, u+1, v) + log X + log S(q+1) + at(p, q) + beta1     (stem 2, Qgl only)
//   s_L = logsumexp(s1H_L, s2)                                        (Qgl only)
//   s2H_L = log S(q+1) + log X + log Qz[u+1][v-1] + log H2_L(p, u, v) - log S(p) + beta1
//           + 2 beta2  (Qg only: mirrored stem-2 estimate of --stem-est-role all)
//   sH_L = logsumexp(s1H_L, s2H_L)                                    (Qg only)
// with A, B from the exact Qgl table.  Per table, score, role class (all, stem1, stem2, mixed):
// per-step Spearman with log w (sum of rho x states, sum of states; live states only), set
// recall@b of the per-step top b by w (steps with > b states) and mass recall@b (all steps),
// b = 100, 300, 1000; per L the H = 0 states and the states the engine would remove.
std::string DP03::stem_diag(const std::vector<const StemEst*>& ests, const std::string& id) {
  using namespace stemdiag;
  auto t0 = std::chrono::steady_clock::now();
  outside();
  const int NL = static_cast<int>(ests.size());
  const size_t BEAMS[3] = {100, 300, 1000};
  std::vector<Tab4>* tabs[2] = {&Qg, &Qgl};
  auto flat_b = [&](int t, std::vector<double>& out) {
    out.clear();
    for (Tab4& tb : *tabs[t])
      for (Grp& g : tb.grps)
        for (Ent& e : g.ents) out.push_back(e.b);
  };
  std::vector<double> b0[2], b1[2], b2[2];
  for (int t = 0; t < 2; ++t) flat_b(t, b0[t]);
  tagged_outside(1);
  for (int t = 0; t < 2; ++t) flat_b(t, b1[t]);
  tagged_outside(2);
  for (int t = 0; t < 2; ++t) flat_b(t, b2[t]);
  for (int t = 0; t < 2; ++t) {
    size_t k = 0;
    for (Tab4& tb : *tabs[t])
      for (Grp& g : tb.grps)
        for (Ent& e : g.ents) e.b = b0[t][k++];
  }
  double tag_err = 0.0;
  long long tag_bad = 0;
  for (int t = 0; t < 2; ++t)
    for (size_t k = 0; k < b0[t].size(); ++k) {
      double s = log_add(b1[t][k], b2[t][k]);
      if (b0[t][k] > HALF_NEG && s > HALF_NEG) tag_err = std::max(tag_err, std::fabs(s - b0[t][k]));
      else if ((b0[t][k] > HALF_NEG) != (s > HALF_NEG)) ++tag_bad;
    }
  stA.clear();
  stE.clear();
  stc_q = -1;
  first_f = stemidx::MinTree(n + 1);
  for (int q = 0; q < n; ++q) stem_add_A(q);

  const double w_pk = E.w_beta1 + 2.0 * E.w_stem;
  const std::vector<double>& S = ests[0]->S;
  std::string out = "{\"id\":\"" + id + "\",\"n\":" + std::to_string(n) + ",\"logZ\":";
  num(out, logZ);
  out += ",\"loops\":[";
  for (int l = 0; l < NL; ++l) out += (l ? "," : "") + std::to_string(ests[l]->loop);
  out += "],\"stem\":[";
  for (int l = 0; l < NL; ++l)
    out += std::string(l ? "," : "") + "{\"entries\":" + std::to_string(ests[l]->entries) +
           ",\"dropped\":" + std::to_string(ests[l]->dropped) +
           ",\"nodes\":" + std::to_string(ests[l]->index.nodes()) +
           ",\"bytes\":" + std::to_string(ests[l]->bytes()) + "}";
  out += "],\"tag_err\":";
  num(out, tag_err);
  out += ",\"tag_bad\":" + std::to_string(tag_bad) + ",\"tables\":{";

  for (int t = 0; t < 2; ++t) {
    const bool gl = t == 1;
    std::vector<std::string> names = {"s1"};
    for (int l = 0; l < NL; ++l) names.push_back("s1H_L" + std::to_string(ests[l]->loop));
    if (gl) {
      names.push_back("s2");
      for (int l = 0; l < NL; ++l) names.push_back("s_L" + std::to_string(ests[l]->loop));
    } else {
      for (int l = 0; l < NL; ++l) names.push_back("s2H_L" + std::to_string(ests[l]->loop));
      for (int l = 0; l < NL; ++l) names.push_back("sH_L" + std::to_string(ests[l]->loop));
    }
    const int NS = static_cast<int>(names.size());
    // accumulators [score][class]
    std::vector<double> sp_sum(NS * 4, 0.0), sp_w(NS * 4, 0.0);
    std::vector<double> hit(NS * 4 * 3, 0.0), tot(NS * 4 * 3, 0.0), mk(NS * 4 * 3, 0.0),
        mt(NS * 4 * 3, 0.0);
    double cnt[4] = {}, mass[4] = {};
    long long states = 0, live_states = 0, dead_states = 0;
    // per L: H = 0 states, their live count / mass, with stem-1 outside, removable (+ live)
    std::vector<double> h0n(NL, 0), h0live(NL, 0), h0mass(NL, 0), h0b1(NL, 0), rmn(NL, 0),
        rmlive(NL, 0), rmmass(NL, 0);
    size_t off = 0;
    std::vector<double> w, sc_flat, rw, rs, sub_w, sub_s;
    std::vector<int> cls;
    std::vector<char> ks, ke;
    for (int q = 0; q < n; ++q) {
      Tab4& tb = (*tabs[t])[q];
      size_t N = 0;
      for (Grp& g : tb.grps) N += g.ents.size();
      if (N == 0) continue;
      w.assign(N, NEG_INF);
      cls.assign(N, -1);
      sc_flat.assign(static_cast<size_t>(NS) * N, NEG_INF);
      size_t k = 0;
      for (Grp& g : tb.grps)
        for (Ent& e : g.ents) {
          const size_t fk = off + k;
          const int p = e.p, u = e.u, v = g.v;
          if (b0[t][fk] > HALF_NEG) {
            w[k] = e.a + b0[t][fk] - logZ;
            double r1 = b1[t][fk] - log_add(b1[t][fk], b2[t][fk]);
            double share = b1[t][fk] > HALF_NEG ? std::exp(r1) : 0.0;
            cls[k] = share >= 0.999 ? 1 : share <= 0.001 ? 2 : 3;
          }
          const double s1 = four_score(p, u, v, e.a);
          sc_flat[k] = s1;
          double bval = NEG_INF, s2 = NEG_INF;
          if (gl) {
            if (E.end_ok(p, q)) bval = stem_B(p, u + 1, v);
            if (bval > HALF_NEG) s2 = bval + e.a + S[q + 1] + E.at(p, q) + E.w_beta1;
            sc_flat[static_cast<size_t>(1 + NL) * N + k] = s2;
          }
          for (int l = 0; l < NL; ++l) {
            double h = ests[l]->H(u, v, q);
            double s1h = h < HALF_NEG ? NEG_INF : s1 + h - ests[l]->C0[q] + w_pk;
            sc_flat[static_cast<size_t>(1 + l) * N + k] = s1h;
            if (gl) {
              sc_flat[static_cast<size_t>(2 + NL + l) * N + k] = log_add(s1h, s2);
            } else {
              double h2 = ests[l]->H2(p, u, v), s2h = NEG_INF;
              if (h2 > HALF_NEG) {
                s2h = S[q + 1] + e.a;
                if (score_qz) s2h += qz_log(u + 1, v - 1);
                s2h = s2h + h2 - S[p] + w_pk;
              }
              sc_flat[static_cast<size_t>(1 + NL + l) * N + k] = s2h;
              sc_flat[static_cast<size_t>(1 + 2 * NL + l) * N + k] = log_add(s1h, s2h);
            }
            const bool live = w[k] > HALF_NEG;
            if (h < HALF_NEG) {
              h0n[l] += 1;
              if (live) h0live[l] += 1, h0mass[l] += std::exp(w[k]);
              if (b1[t][fk] > HALF_NEG) h0b1[l] += 1;
            }
            bool rm = gl ? (!E.end_ok(p, q) || (h < HALF_NEG && bval < HALF_NEG))
                         : (h < HALF_NEG && first_f.query(u + 1, v) > v - 1);
            if (rm) {
              rmn[l] += 1;
              if (live) rmlive[l] += 1, rmmass[l] += std::exp(w[k]);
            }
          }
          ++k;
        }
      off += N;
      states += N;
      // live counts / mass per class (0 = all)
      for (size_t i = 0; i < N; ++i) {
        if (cls[i] < 0) {
          ++dead_states;
          continue;
        }
        ++live_states;
        double ew = std::exp(w[i]);
        cnt[0] += 1, mass[0] += ew;
        cnt[cls[i]] += 1, mass[cls[i]] += ew;
      }
      // Spearman per score / class
      for (int c = 0; c < 4; ++c) {
        std::vector<size_t> idx;
        for (size_t i = 0; i < N; ++i)
          if (cls[i] >= 0 && (c == 0 || cls[i] == c)) idx.push_back(i);
        if (idx.size() < 3) continue;
        sub_w.resize(idx.size());
        for (size_t m = 0; m < idx.size(); ++m) sub_w[m] = w[idx[m]];
        avg_ranks(sub_w, rw);
        for (int s = 0; s < NS; ++s) {
          sub_s.resize(idx.size());
          for (size_t m = 0; m < idx.size(); ++m) sub_s[m] = sc_flat[static_cast<size_t>(s) * N + idx[m]];
          avg_ranks(sub_s, rs);
          double rho = pearson(rw, rs);
          if (std::isnan(rho)) continue;
          sp_sum[s * 4 + c] += rho * idx.size();
          sp_w[s * 4 + c] += idx.size();
        }
      }
      // recall
      std::vector<size_t> live;
      std::vector<double> lw;
      for (size_t i = 0; i < N; ++i)
        if (cls[i] >= 0) live.push_back(i), lw.push_back(w[i]);
      for (int bi = 0; bi < 3; ++bi) {
        const size_t b = BEAMS[bi];
        std::vector<char> inE(N, 0);
        top_b(lw, b, ke);
        for (size_t m = 0; m < live.size(); ++m) inE[live[m]] = ke[m];
        for (int s = 0; s < NS; ++s) {
          std::vector<double> col(sc_flat.begin() + static_cast<size_t>(s) * N,
                                  sc_flat.begin() + static_cast<size_t>(s + 1) * N);
          top_b(col, b, ks);
          for (size_t i : live) {
            const double ew = std::exp(w[i]);
            for (int c : {0, cls[i]}) {
              const size_t a = (static_cast<size_t>(s) * 4 + c) * 3 + bi;
              mt[a] += ew;
              if (ks[i]) mk[a] += ew;
              if (N > b && inE[i]) {
                tot[a] += 1;
                if (ks[i]) hit[a] += 1;
              }
            }
          }
        }
      }
    }
    auto arr = [&](const std::vector<double>& v) {
      std::string o = "[";
      for (size_t i = 0; i < v.size(); ++i) {
        if (i) o += ",";
        num(o, v[i]);
      }
      return o + "]";
    };
    out += std::string(t ? "," : "") + "\"" + (gl ? "Qgl" : "Qg") + "\":{\"scores\":[";
    for (int s = 0; s < NS; ++s) out += std::string(s ? "," : "") + "\"" + names[s] + "\"";
    out += "],\"states\":" + std::to_string(states) + ",\"live\":" + std::to_string(live_states) +
           ",\"dead\":" + std::to_string(dead_states) + ",\"count\":" +
           arr(std::vector<double>(cnt, cnt + 4)) + ",\"mass\":" +
           arr(std::vector<double>(mass, mass + 4)) + ",\"sp_sum\":" + arr(sp_sum) +
           ",\"sp_w\":" + arr(sp_w) + ",\"hit\":" + arr(hit) + ",\"tot\":" + arr(tot) +
           ",\"mkept\":" + arr(mk) + ",\"mtot\":" + arr(mt) + ",\"h0_n\":" + arr(h0n) +
           ",\"h0_live\":" + arr(h0live) + ",\"h0_mass\":" + arr(h0mass) + ",\"h0_b1\":" +
           arr(h0b1) + ",\"rm_n\":" + arr(rmn) + ",\"rm_live\":" + arr(rmlive) +
           ",\"rm_mass\":" + arr(rmmass) + "}";
  }
  out += "},\"diag_sec\":";
  num(out, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
  out += "}";
  return out;
}

// --l2r-p1 time split: process CPU seconds and peak RSS so far.
static double cpu_sec() { return static_cast<double>(std::clock()) / CLOCKS_PER_SEC; }
static double ru_cpu() {
  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru) != 0) return 0.0;
  return static_cast<double>(ru.ru_utime.tv_sec) + ru.ru_utime.tv_usec * 1e-6 +
         static_cast<double>(ru.ru_stime.tv_sec) + ru.ru_stime.tv_usec * 1e-6;
}
static double hwm_mb() {
  std::ifstream f("/proc/self/status");
  std::string k;
  while (f >> k)
    if (k == "VmHWM:") {
      double kb = 0.0;
      f >> kb;
      return kb / 1024.0;
    }
  return 0.0;
}

int main(int argc, char** argv) {
  std::setvbuf(stderr, nullptr, _IOLBF, 0);
  std::signal(SIGXCPU, on_cpu_stop);
  std::signal(SIGUSR1, on_cpu_stop);
  int beam = 0;
  int loop_cap = MAX_LOOP;
  // Temporary per-pass loop caps (--loop-cap-p1 / -p2 / -scan). -1 = use --loop-cap.
  // 0 / off / n = uncapped (NO_CAP). Also drives max_loop for that pass.
  constexpr int NO_LOOP_CAP = 1000000;
  int loop_cap_p1 = -1, loop_cap_p2 = -1, loop_cap_scan = -1;
  bool score_qz = true;
  bool allow_pk = true;
  double cutoff = 1e-6;
  bool do_bpp = true;
  int pk_l2r = 0, l2r_beam = -1, l2r_spawn_k = 0;
  bool l2r_spawn_cube = false, l2r_spawn_static = false;  // --l2r-spawn all|cube|cube-static|inc
  bool l2r_spawn_inc = false;
  int l2r_pbeam[3] = {1000, 1000, 1000};  // --l2r-pbeam b1,b2,b3 (P1, P2, P3; 0 or inf: none)
  int l2r_p_hk = 64;                      // --l2r-p-hk K (H entries per suffix)
  int l2r_p3_split[2] = {0, 0};           // --l2r-p3-split Ba,bh
  bool l2r_p_raw = false;                 // --l2r-p-aopen q|raw
  bool l2r_p3_hex = false;                // --l2r-p3-hex on
  bool l2r_p3_lin = false;                // --l2r-p3-lin on
  bool l2r_inc_check = false;             // --l2r-inc-check on
  bool l2r_ninf_mass = false, l2r_keep_ninf = false;  // --l2r-ninf-mass on, --l2r-keep-ninf on
  bool l2r_qz_audit = false, ntr_z_set = false;
  bool l2r_spawn_wk = false;                // --l2r-spawn-w suf|k
  int l2r_spawn_pool = 128, l2r_spawn_maxpop = 0;  // --l2r-spawn-pool, --l2r-spawn-maxpop (0: 4 pool)
  std::string l2r_tbeam;  // --l2r-tbeam NAME=B[,NAME=B..], B = 0 or inf: no beam
  std::string l2r_odump, l2r_oload, l2r_sub, l2r_trace, l2r_trace_key;
  bool l2r_perf = false;
  bool k_memo = false, k_bound = false, k_exact_ub = false, k_exact_check = false, k_exact = false;
  bool rlt = false, next_pair = false, grammar_simple = false, c_trace = false, pass1_only = false;
  int rlt_beam = 0;
  std::string sample_prof;
  bool l2r_ctx = false, l2r_kx = false;  // --l2r-ctx on: span context (needs Bml); --l2r-kx on: exact K
  bool l2r_ctx_p1 = false;               // --l2r-ctx p1: pass-1 context ratio lambda(a')
  double l2r_ctx_tau = 0.1;              // --l2r-ctx-tau
  int l2r_ctx_rounds = 1;                // --l2r-ctx-rounds: pass-1 rounds that recompute lambda
  bool l2r_split = false;                // --l2r-split on: consumer split of the outsides (diag)
  int l2r_p1 = 0;                        // --l2r-p1 R: rounds of the beamed pass-1 sweep
  int l2r_p1_beam = 250;                 // --l2r-p1-beam B
  bool l2r_p1_own = true;                // --l2r-p1-qb pass1|nested
  double pk_db[3] = {0.0, 0.0, 0.0};     // --pk-dbeta1 / -1m / -1p: kcal/mol added after --params
  bool pk_l2r_R = false, l2r_diag = false;
  bool arm_score = false;
  std::string access_file;
  int max_gap_scan = 200;
  int max_arm_len = 0;
  bool gap_step = false;
  bool loop_cap_set = false;
  bool restrict_pairs = false;
  std::vector<std::pair<int, int>> restrict_list;
  GateCfg gcfg;
  bool gate_stats_out = false;
  std::string gate_diag_path;
  std::vector<double> diag_taus = {0.0, 1e-4, 1e-3, 1e-2, 0.05, 0.1};
  int long_span = 100;
  int table_beam[BT_N];
  std::fill(table_beam, table_beam + BT_N, -1);
  std::string states_path, pruned_path, exact_path, scores_path, ext_path, kept_path, beam_stats_path;
  double dump_min_w = 0.0, pruned_min_w = 0.0;
  std::unordered_set<uint64_t> remove_keys;
  bool stem_rank = false, stem_both = false, stem_qg2 = false;
  bool stem_qg_partner = false, stem_consistent = false, stem_qg_r_set = false,
       stem_consts_set = false, stem_bg_check = false;
  int stem_ext = 0, stem_ext_max = 30;
  bool stem_split = false, stem_ext_set = false, stem_ext_check = false;
  int grow_lg = 3, grow_nodes = 2000, ml_span = 30;
  double grow_cut = 1e-8;
  int stem_helpers = 0;
  bool grow_set = false, helpers_set = false;
  bool stem_profile = false, profile_set = false, ml_tail = false, prof_check = false;
  bool bg_reorder = false, ro_check = false;
  double prof_eps = 0.0;
  bool stem_chain = false, chain_set = false, chain_ml = true;
  double chain_floor = 1e-7, chain_cut = 1e-8;
  int chain_span = 150, chain_loop = 30;
  bool stem_hr_av = false, stem_hr_phi = false, stem_hr_set = false, stem_hr_check = false;
  double stem_hr_gamma = 1.0;
  bool stem_pr = false, stem_pr_nested = false, pr_set = false, pr_check = false;
  int pr_span = 300, pr_ml_span = 150, pr_ml_close = 30, pr_arms_k = 64;
  double pr_floor = 1e-7, pr_cut = 1e-8;
  std::string pr_oracle_path, pr_oracle_set = "step1";
  bool pr_ratio = false, pr_gapfree = false, pr_cut_roots = false;
  // --pr-arms exact|sep, --pr-arm-z qz|unpaired, --pr-partner-floor (sep default 1e-20; off:
  // the kept pairs), --pr-test-qz unpaired, --pr-sep-check; sep defaults --pr-cut-by gapfree
  int pr_arms_mode = 0, pr_arm_z = 0;
  bool pr_qz_unp = false, pr_sep_check = false, pr_cut_by_set = false;
  double pr_partner_floor = -1.0;
  double pr_arms_eps = 1e-3;
  int pr_arms_max = 256;
  // --pr-grow-select floor|bound, --pr-grow-eps, --pr-grow-B, --pr-grow-check, --pr-score-dump
  bool pr_gsel = false, pr_gcheck = false;
  bool pr_aq_sep = false;  // --pr-armq-sep on|off
  size_t pr_aq_max = 128;  // --pr-aq-max: Step A entries kept
  bool pr_arm_persist = false, pr_arm_gaponly = false;  // --pr-arm-persist, --pr-arm-gaponly
  int pr_pspan = 0;                                     // --pr-partner-span W
  int pr_pmax = 0;         // --pr-partner-maxpairs K (0 = off)
  bool pr_pcheck = false;  // --pr-partner-check
  bool pr_pml = true, pr_pmlcheck = false;  // --pr-partner-ml, --pr-partner-ml-check
  int pr_fast = 0;                          // --pr-fast off|acc|ml
  int pr_abl = 0;                           // --pr-abl bml,q,p,xs,ks,sep
  int pr_abl_pk = 0;                        // --pr-abl-pk
  double pr_veps = 1e-8;                                // --pr-v-eps (0: no trimming)
  bool pr_vsuffix = false;                              // --pr-v-trim global|suffix
  int pr_vcap = 0;                                      // --bV (was --pr-v-cap); default 256
  int p2_qz_beam = 0, p2_qz_key = 1, p2_v_mode = 0, p2_v_bp = 0, p2_v_cap = 0, p2_v_trim = 0;
  bool bV_set = false;
  int p2_g_mode = 0, p2_g_bp = 0, p2_g_cap = 0, nest_beam = -1, pre_beam = 0, pre_key = 0;
  bool b1_set = false, rlt_beam_set = false, sweep_set = false, vsweep_set = false;
  int b2_width = 0;
  int b2_at = -1, pbeam_explicit_at = -1, pbeam_ratio_at = -1, pbr[3] = {1, 2, 2};
  bool user_pbeam = false, user_pratio = false;
  bool sw_check = false, p2_sweep = false, p2_vsweep = false;
  int p2_v_floor = 0;  // 0 off, 1 smallest kept mu, 2 log-weight 0
  bool p2_v_cut_fb = false;
  bool lazy_ml_cap = false;
  bool p2_g_suffix = false, p2_g_eager = false, p2_g_fb = false, profile = false, dump_qb_on = false;
  bool qzbar_old = false;  // --qzbar-old: dense Qz→prQz like Oct-8 incperf (no sparse Qz_pre)
  int pr_vmlk = 0;                                      // --pr-v-mlk K (0: all)
  bool pr_vmlsteps = false;                             // --pr-v-mlsteps on
  bool pr_vxs = false;                                  // --pr-v-key xsys|xs
  bool pr_vstream = false;                              // --pr-v-stream off|on
  bool l2r_arm_fb = false;                              // --l2r-arm-fb on|off
  // --l2r-lean on (--pk-l2r on with a beam): skip the pre-pass structures only the b16
  // scoring reads (beamed nested run with its stem index, S, C0, nested outside; chain;
  // prKept; prQbIn unless pass 1 or the V diagnostics need it)
  bool l2r_lean = false;
  bool pr_v_ratio = false;  // --pr-v-ratio on|off
  std::string pr_v_diag;    // --pr-v-diag PREFIX
  double pr_v_sample = 0.0; // --pr-v-sample P
  double pr_geps = 1e-6;
  int pr_gB = 512;
  std::string pr_dump_path;
  int stem_loop = 2;
  std::string stem_diag_path;
  std::vector<int> stem_diag_loops = {0, 2, 6};

  // --preset default|accurate expands in place (options after it override; before are
  // overridden). Energy model (--material, --dangles, --param-file) is not part of it.
  // Paper beams are --b1/--b2/--b4/--bV/--bG; pbeam is derived as (b2,b4,b4) unless
  // --l2r-pbeam or --pbeam-ratio is given. accurate = paper values; default = smaller b4/bG.
  // Scripts should pass --preset only (no extra beam flags).

  for (int a = 1; a < argc; ++a) {
    if (std::strcmp(argv[a], "--l2r-pbeam") == 0) user_pbeam = true;
    if (std::strcmp(argv[a], "--pbeam-ratio") == 0) user_pratio = true;
  }
  std::vector<std::string> args_s;
  for (int a = 0; a < argc; ++a) {
    if (a > 0 && std::string(argv[a]) == "--preset") {
      const std::string p = a + 1 < argc ? argv[++a] : "";
      if (p != "default" && p != "accurate") {
        std::fprintf(stderr, "--preset must be 'default' or 'accurate'\n");
        return 2;
      }
      const bool acc = (p == "accurate");
      // Paper preset beams (accurate) / smaller default. pbeam = (b2,b4,b4) after parse.
      args_s.push_back("--b1");
      args_s.push_back("200");
      args_s.push_back("--b2");
      args_s.push_back(acc ? "1000" : "200");
      args_s.push_back("--b4");
      args_s.push_back(acc ? "3000" : "1000");
      args_s.push_back("--bV");
      args_s.push_back(acc ? "256" : "512");
      args_s.push_back("--bG");
      args_s.push_back(acc ? "32" : "16");
      for (const char* x :
           {"--gap-beam", "step", "--loop-cap", "30", "--max-loop", "30", "--wc-ends", "--rt", "0.6163207755",
            "--cutoff", "1e-5", "--stem-est", "rank", "--stem-est-role", "all", "--stem-est-loop", "0",
            "--stem-est-qg-r", "partner", "--stem-est-qg-ext", "loop", "--stem-est-grow-loop", "2",
            "--stem-est-helpers", "grow", "--stem-est-profile", "on", "--stem-est-chain", "on",
            "--stem-est-principled", "on", "--pr-partner-floor", "1e-20", "--pr-grow-select", "bound",
            "--pr-grow-eps", "1e-6", "--pr-cut", "0", "--pr-armq-sep", "on", "--pr-fast", "ml", "--pr-abl", "sepv",
            "--grow-rec", "on", "--grow-ml-c", "16", "--grow-trunc", "16", "--grow-trunc-mode", "side",
            "--pk-l2r", "on", "--l2r-spawn", "inc",
            "--l2r-p3-hex", "on",
            "--l2r-p3-lin", "on", "--pr-ml-span", "100000", "--grow-ml-ub", "on", "--pr-v-key", "xs",
            "--pr-v-eps", "0", "--pr-v-trim", "suffix", "--pr-v-stream", "on",
            "--pr-v-mlk", "0", "--grow-node-trunc-r", "-1",
            "--grow-node-trunc-z", "bar", "--grow-r-from-v", "on", "--malloc-trim", "on",
            "--sweep", "on", "--v-sweep", "on", "--grammar", "full", "--ml-lead-cap", "30"})
        args_s.push_back(x);
    } else if (a > 0 && std::string(argv[a]) == "--stem-est-preset") {
      const std::string p = a + 1 < argc ? argv[++a] : "";
      if (p != "v1") {
        std::fprintf(stderr, "--stem-est-preset must be 'v1'\n");
        return 2;
      }
      for (const char* x : {"--stem-est", "rank", "--stem-est-role", "all", "--stem-est-loop", "0",
                            "--stem-est-qg-r", "partner", "--stem-est-qg-ext", "loop",
                            "--stem-est-grow-loop", "2", "--stem-est-helpers", "grow",
                            "--stem-est-profile", "on", "--stem-est-chain", "on",
                            "--stem-est-chain-loop", "2", "--stem-est-chain-ml", "off"})
        args_s.push_back(x);
    } else {
      args_s.push_back(argv[a]);
    }
  }
  if (static_cast<int>(args_s.size()) != argc) {
    std::string line = "[preset] expanded:";
    for (std::size_t k = 1; k < args_s.size(); ++k) line += " " + args_s[k];
    std::fprintf(stderr, "%s\n", line.c_str());
  }
  std::vector<char*> args_v;
  for (std::string& s : args_s) args_v.push_back(s.data());
  args_v.push_back(nullptr);
  argc = static_cast<int>(args_s.size());
  argv = args_v.data();

  for (int a = 1; a < argc; ++a) {
    std::string arg = argv[a];
    auto next = [&]() -> const char* {
      if (a + 1 >= argc) {
        usage();
        std::exit(2);
      }
      return argv[++a];
    };
    if (arg == "--beam" || arg.rfind("--beam=", 0) == 0) {
      std::fprintf(stderr, "--beam is removed; use --b4\n");
      return 2;
    }
    else if (arg == "--b4") {
      std::string b = next();
      beam = (b == "inf" || b == "infinity" || b == "none") ? 0 : std::atoi(b.c_str());
    }
    else if (arg == "--loop-cap") {
      loop_cap = std::atoi(next());
      loop_cap_set = true;
    }
    else if (arg == "--loop-cap-p1" || arg == "--loop-cap-p2" || arg == "--loop-cap-scan") {
      // Temporary: per-pass --loop-cap / --max-loop. "off"/"n"/0 = uncapped.
      auto parse_pc = [](const std::string& v) -> int {
        if (v == "n" || v == "off" || v == "0") return NO_LOOP_CAP;
        return std::atoi(v.c_str());
      };
      const int v = parse_pc(next());
      if (arg == "--loop-cap-p1") loop_cap_p1 = v;
      else if (arg == "--loop-cap-p2") loop_cap_p2 = v;
      else loop_cap_scan = v;
    }
    else if (arg == "--ml-lead-cap") {
      std::string v = next();
      if (v == "n" || v == "off") ml_lead_cap = 0;
      else ml_lead_cap = std::atoi(v.c_str());
    }
    else if (arg == "--max-loop") max_loop = std::atoi(next());
    else if (arg == "--no-gu") allow_gu = false;
    else if (arg == "--no-single-stems") single_stems = false;
    else if (arg == "--wc-ends") wc_ends = true;
    else if (arg == "--pair-bonus") pair_bonus = std::atof(next());
    else if (arg == "--restrict-pairs") restrict_pairs = true;
    else if (arg == "--rt") RT = std::atof(next());
    else if (arg == "--score") {
      std::string s = next();
      if (s == "qz") score_qz = true;
      else if (s == "prefix") score_qz = false;
      else {
        std::fprintf(stderr, "score must be 'qz' or 'prefix'\n");
        return 2;
      }
    } else if (arg == "--no-pk") allow_pk = false;
    else if (arg == "--cutoff") cutoff = std::atof(next());
    else if (arg == "--no-bpp") do_bpp = false;
    else if (arg == "--arm-score") {
      std::string s = next();
      if (s == "on") arm_score = true;
      else if (s == "off") arm_score = false;
      else {
        std::fprintf(stderr, "--arm-score must be 'on' or 'off'\n");
        return 2;
      }
    } else if (arg == "--access-file") access_file = next();
    else if (arg == "--max-gap-scan") max_gap_scan = std::atoi(next());
    else if (arg == "--max-arm-len") max_arm_len = std::atoi(next());
    else if (arg == "--gap-beam") {
      std::string m = next();
      if (m == "group") gap_step = false;
      else if (m == "step") gap_step = true;
      else {
        std::fprintf(stderr, "--gap-beam must be 'group' or 'step'\n");
        return 2;
      }
    }
    else if (arg == "--params") {
      std::string m = next();
      if (m == "dp03") PK = DP03_PARAMS;
      else if (m == "dp09") PK = DP09_PARAMS;
      else if (m == "dp09st") PK = DP09ST_PARAMS;
      else if (m == "dp09ml") PK = DP09ML_PARAMS;
      else if (m == "dp09full") PK = DP09FULL_PARAMS;
      else {
        std::fprintf(stderr, "--params must be dp03, dp09, dp09st, dp09ml or dp09full\n");
        return 2;
      }
    }
    else if (arg == "--material") {
      std::string m = next();
      mat_rna1999 = m == "rna1999";
      mat_v04 = m == "vienna2004";
      if (m != "dp03" && !mat_rna1999 && !mat_v04) {
        std::fprintf(stderr, "--material must be dp03, rna1999 or vienna2004\n");
        return 2;
      }
    }
    else if (arg == "--param-file") param_file = next();
    else if (arg == "--dangles") {
      std::string m = next();
      if (m == "none") dangles_mode = 0;
      else if (m == "some") dangles_mode = 1;
      else if (m == "all") dangles_mode = 2;
      else if (m == "d2") dangles_mode = 3;
      else {
        std::fprintf(stderr, "--dangles must be none, some, all or d2\n");
        return 2;
      }
    }
    else if (arg == "--gate") {
      std::string s = next();
      if (s == "on") gcfg.on = true;
      else if (s == "off") gcfg.on = false;
      else {
        std::fprintf(stderr, "--gate must be 'on' or 'off'\n");
        return 2;
      }
    }
    else if (arg == "--gate-tau") gcfg.tau = std::strtod(next(), nullptr);
    else if (arg == "--gate-span") gcfg.span = std::atoi(next());
    else if (arg == "--gate-delta") gcfg.delta = std::strtod(next(), nullptr);
    else if (arg == "--gate-role") {
      std::string s = next();
      if (s == "stem1") gcfg.both = false;
      else if (s == "both") gcfg.both = true;
      else {
        std::fprintf(stderr, "--gate-role must be 'stem1' or 'both'\n");
        return 2;
      }
    }
    else if (arg == "--gate-tables") {
      std::string s = next();
      if (s == "qgl") gcfg.tables[GATE_QG] = false, gcfg.tables[GATE_QGL] = true;
      else if (s == "qg") gcfg.tables[GATE_QG] = true, gcfg.tables[GATE_QGL] = false;
      else if (s == "qg,qgl" || s == "qgl,qg") gcfg.tables[GATE_QG] = gcfg.tables[GATE_QGL] = true;
      else {
        std::fprintf(stderr, "--gate-tables must be qg, qgl or qg,qgl\n");
        return 2;
      }
    }
    else if (arg == "--gate-reserve") gcfg.reserve = std::atoi(next());
    else if (arg == "--gate-stats") gate_stats_out = true;
    else if (arg == "--stem-est") {
      std::string s = next();
      if (s == "rank") stem_rank = true;
      else if (s == "off") stem_rank = false;
      else {
        std::fprintf(stderr, "--stem-est must be 'off' or 'rank'\n");
        return 2;
      }
    }
    else if (arg == "--stem-est-loop") stem_loop = std::atoi(next());
    else if (arg == "--stem-est-role") {
      std::string s = next();
      stem_both = s == "both" || s == "all";
      stem_qg2 = s == "all";
      if (s != "stem1" && s != "both" && s != "all") {
        std::fprintf(stderr, "--stem-est-role must be 'stem1', 'both' or 'all'\n");
        return 2;
      }
    }
    else if (arg == "--stem-est-qg-r") {
      std::string s = next();
      stem_qg_r_set = true;
      stem_qg_partner = s == "partner";
      if (s != "mirror" && s != "partner") {
        std::fprintf(stderr, "--stem-est-qg-r must be 'mirror' or 'partner'\n");
        return 2;
      }
    }
    else if (arg == "--stem-est-consts") {
      std::string s = next();
      stem_consts_set = true;
      stem_consistent = s == "consistent";
      if (s != "as-is" && s != "consistent") {
        std::fprintf(stderr, "--stem-est-consts must be 'as-is' or 'consistent'\n");
        return 2;
      }
    }
    else if (arg == "--stem-est-qg-ext") {
      std::string s = next();
      stem_ext_set = true;
      stem_ext = s == "stack" ? 1 : s == "bulge1" ? 2 : s == "loop" ? 3 : 0;
      if (s != "off" && s != "stack" && s != "bulge1" && s != "loop") {
        std::fprintf(stderr, "--stem-est-qg-ext must be 'off', 'stack', 'bulge1' or 'loop'\n");
        return 2;
      }
    }
    else if (arg == "--stem-est-grow-loop") {
      grow_set = true;
      grow_lg = std::atoi(next());
    }
    else if (arg == "--stem-est-grow-cut") {
      grow_set = true;
      grow_cut = std::atof(next());
    }
    else if (arg == "--stem-est-grow-nodes") {
      grow_set = true;
      grow_nodes = std::atoi(next());
    }
    else if (arg == "--stem-est-helpers") {
      std::string s = next();
      helpers_set = true;
      stem_helpers = s == "grow" ? 1 : s == "nested" ? 2 : 0;
      if (s != "off" && s != "grow" && s != "nested") {
        std::fprintf(stderr, "--stem-est-helpers must be 'off', 'grow' or 'nested'\n");
        return 2;
      }
    }
    else if (arg == "--stem-est-ml-span") {
      helpers_set = true;
      ml_span = std::atoi(next());
    }
    else if (arg == "--stem-est-profile") {
      std::string s = next();
      stem_profile = s == "on";
      if (s != "off" && s != "on") {
        std::fprintf(stderr, "--stem-est-profile must be 'off' or 'on'\n");
        return 2;
      }
    }
    else if (arg == "--stem-est-profile-eps") {
      profile_set = true;
      prof_eps = std::strtod(next(), nullptr);
    }
    else if (arg == "--stem-est-ml-tail") {
      std::string s = next();
      profile_set = true;
      ml_tail = s == "nested";
      if (s != "off" && s != "nested") {
        std::fprintf(stderr, "--stem-est-ml-tail must be 'off' or 'nested'\n");
        return 2;
      }
    }
    else if (arg == "--stem-profile-check") prof_check = true;
    else if (arg == "--stem-est-bg-reorder") {
      std::string s = next();
      bg_reorder = s == "on";
      if (s != "off" && s != "on") {
        std::fprintf(stderr, "--stem-est-bg-reorder must be 'off' or 'on'\n");
        return 2;
      }
    }
    else if (arg == "--stem-bg-reorder-check") ro_check = true;
    else if (arg == "--stem-est-chain" || arg == "--stem-est-chain-ml") {
      std::string s = next();
      if (arg == "--stem-est-chain-ml") chain_set = true;
      (arg == "--stem-est-chain" ? stem_chain : chain_ml) = s == "on";
      if (s != "off" && s != "on") {
        std::fprintf(stderr, "%s must be 'off' or 'on'\n", arg.c_str());
        return 2;
      }
    }
    else if (arg == "--stem-est-chain-floor") {
      chain_set = true;
      chain_floor = std::strtod(next(), nullptr);
    }
    else if (arg == "--stem-est-chain-span") {
      chain_set = true;
      chain_span = std::atoi(next());
    }
    else if (arg == "--stem-est-chain-loop") {
      chain_set = true;
      chain_loop = std::atoi(next());
    }
    else if (arg == "--stem-est-chain-cut") {
      chain_set = true;
      chain_cut = std::strtod(next(), nullptr);
    }
    else if (arg == "--stem-est-hr") {
      std::string s = next();
      stem_hr_set = true;
      stem_hr_av = s == "av";
      stem_hr_phi = s == "phi";
      if (s != "orig" && s != "av" && s != "phi") {
        std::fprintf(stderr, "--stem-est-hr must be 'orig', 'av' or 'phi'\n");
        return 2;
      }
    }
    else if (arg == "--stem-est-hr-gamma") {
      stem_hr_set = true;
      stem_hr_gamma = std::strtod(next(), nullptr);
    }
    else if (arg == "--stem-hr-check") stem_hr_check = true;
    else if (arg == "--stem-est-principled" || arg == "--stem-est-principled-nested") {
      std::string s = next();
      (arg == "--stem-est-principled" ? stem_pr : stem_pr_nested) = s == "on";
      if (s != "off" && s != "on") {
        std::fprintf(stderr, "%s must be 'off' or 'on'\n", arg.c_str());
        return 2;
      }
    }
    else if (arg == "--pr-floor" || arg == "--pr-cut") {
      pr_set = true;
      (arg == "--pr-floor" ? pr_floor : pr_cut) = std::strtod(next(), nullptr);
    }
    else if (arg == "--pr-span" || arg == "--pr-ml-span" || arg == "--pr-ml-close" ||
             arg == "--pr-arms-k") {
      pr_set = true;
      (arg == "--pr-span" ? pr_span : arg == "--pr-ml-span" ? pr_ml_span
                                    : arg == "--pr-ml-close" ? pr_ml_close : pr_arms_k) =
          std::atoi(next());
    }
    else if (arg == "--stem-pr-check") pr_check = true;
    else if (arg == "--pr-arms-rank" || arg == "--pr-cut-by" || arg == "--pr-cut-roots") {
      pr_set = true;
      const std::string s = next();
      const char* a = arg == "--pr-arms-rank" ? "ratio" : arg == "--pr-cut-by" ? "gapfree" : "on";
      const char* b = arg == "--pr-arms-rank" ? "phi" : arg == "--pr-cut-by" ? "weight" : "off";
      if (s != a && s != b) {
        std::fprintf(stderr, "%s must be '%s' or '%s'\n", arg.c_str(), b, a);
        return 2;
      }
      (arg == "--pr-arms-rank" ? pr_ratio : arg == "--pr-cut-by" ? pr_gapfree : pr_cut_roots) = s == a;
      if (arg == "--pr-cut-by") pr_cut_by_set = true;
    }
    else if (arg == "--pr-arms" || arg == "--pr-arm-z" || arg == "--pr-test-qz") {
      pr_set = true;
      const std::string s = next();
      const char* a = arg == "--pr-arms" ? "sep" : "unpaired";
      const char* b = arg == "--pr-arms" ? "exact" : arg == "--pr-arm-z" ? "qz" : "off";
      if (s != a && s != b) {
        std::fprintf(stderr, "%s must be '%s' or '%s'\n", arg.c_str(), b, a);
        return 2;
      }
      if (arg == "--pr-arms") pr_arms_mode = s == a ? 1 : 0;
      else if (arg == "--pr-arm-z") pr_arm_z = s == a ? 1 : 0;
      else pr_qz_unp = s == a;
    }
    else if (arg == "--pr-partner-floor") {
      pr_set = true;
      pr_partner_floor = std::strtod(next(), nullptr);
      if (!(pr_partner_floor >= 0.0 && pr_partner_floor < 1.0)) {
        std::fprintf(stderr, "--pr-partner-floor must be in [0, 1)\n");
        return 2;
      }
    }
    else if (arg == "--pr-sep-check") {
      pr_set = true;
      pr_sep_check = true;
    }
    else if (arg == "--pr-grow-select") {
      pr_set = true;
      const std::string s = next();
      if (s != "floor" && s != "bound") {
        std::fprintf(stderr, "--pr-grow-select must be 'floor' or 'bound'\n");
        return 2;
      }
      pr_gsel = s == "bound";
    }
    else if (arg == "--pr-grow-eps") {
      pr_set = true;
      pr_geps = std::strtod(next(), nullptr);
    }
    else if (arg == "--pr-grow-B") {
      pr_set = true;
      pr_gB = std::atoi(next());
    }
    else if (arg == "--pr-grow-check") {
      pr_set = true;
      pr_gcheck = true;
    }
    else if (arg == "--grow-rec") {
      const std::string s = next();
      if (s != "on" && s != "off") {
        std::fprintf(stderr, "--grow-rec must be on or off\n");
        return 2;
      }
      g_grow_rec = s == "on";
    }
    else if (arg == "--grow-ml-L" || arg == "--grow-ml-R") {
      const std::string s = next();
      const int m = s == "all" ? 2 : s == "depth1" ? 1 : s == "none" ? 0 : -1;
      if (m < 0) {
        std::fprintf(stderr, "%s must be all, depth1 or none\n", arg.c_str());
        return 2;
      }
      (arg == "--grow-ml-L" ? g_ml_L : g_ml_R) = m;
    }
    else if (arg == "--grow-ml-c") {
      g_ml_c = std::atoi(next());
    }
    else if (arg == "--grow-ml-ub") {
      g_ml_ub = std::strcmp(next(), "on") == 0;
    }
    else if (arg == "--grow-ml-diag") {
      g_ml_diag = true;
    }
    else if (arg == "--grow-rec-check") {
      g_rec_check = true;
    }
    else if (arg == "--pr-U-nobml") {
      g_U_nobml = true;
    }
    else if (arg == "--grow-trunc") {
      g_tr_K = std::atoi(next());
    }
    else if (arg == "--grow-trunc-eps") {
      g_tr_eps = std::strtod(next(), nullptr);
    }
    else if (arg == "--grow-trunc-mode") {
      const std::string s = next();
      if (s != "side" && s != "pool") {
        std::fprintf(stderr, "--grow-trunc-mode must be side or pool\n");
        return 2;
      }
      g_tr_pool = s == "pool";
    }
    else if (arg == "--pk-l2r") {
      const std::string s = next();
      if (s != "on" && s != "off" && s != "check") {
        std::fprintf(stderr, "--pk-l2r must be on, off or check\n");
        return 2;
      }
      pk_l2r = s == "on" ? 1 : s == "check" ? 2 : 0;
    }
    else if (arg == "--pk-l2r-checkR") {
      pk_l2r_R = true;
    }
    else if (arg == "--l2r-beam") {
      const std::string b = next();
      l2r_beam = (b == "inf" || b == "none") ? 0 : std::atoi(b.c_str());
    }
    else if (arg == "--l2r-diag") {
      l2r_diag = true;
    }
    else if (arg == "--l2r-tbeam") {
      l2r_tbeam = next();
    }
    else if (arg == "--l2r-spawn-k") {
      const std::string b = next();
      l2r_spawn_k = b == "all" ? 0 : std::atoi(b.c_str());
    }
    else if (arg == "--l2r-spawn") {
      const std::string v = next();
      if (v != "all" && v != "cube" && v != "cube-static" && v != "inc") {
        std::fprintf(stderr, "--l2r-spawn must be all, cube, cube-static or inc\n");
        return 2;
      }
      l2r_spawn_cube = v == "cube" || v == "cube-static";
      l2r_spawn_static = v == "cube-static";
      l2r_spawn_inc = v == "inc";
    }
    else if (arg == "--l2r-p3-split") {
      if (std::sscanf(next(), "%d,%d", &l2r_p3_split[0], &l2r_p3_split[1]) != 2 || l2r_p3_split[0] < 1 ||
          l2r_p3_split[1] < 1) {
        std::fprintf(stderr, "--l2r-p3-split needs Ba,bh >= 1\n");
        return 2;
      }
    }
    else if (arg == "--l2r-pbeam") {
      std::stringstream ss(next());
      std::string item;
      int k = 0;
      while (std::getline(ss, item, ',') && k < 3)
        l2r_pbeam[k++] = (item == "inf" || item == "none") ? 0 : std::atoi(item.c_str());
      if (k != 3) {
        std::fprintf(stderr, "--l2r-pbeam needs b1,b2,b3\n");
        return 2;
      }
      pbeam_explicit_at = a;
    }
    else if (arg == "--l2r-p-hk") {
      l2r_p_hk = std::atoi(next());
      if (l2r_p_hk < 1) {
        std::fprintf(stderr, "--l2r-p-hk must be >= 1\n");
        return 2;
      }
    }
    else if (arg == "--l2r-p-aopen") {
      const std::string v = next();
      if (v != "q" && v != "raw") {
        std::fprintf(stderr, "--l2r-p-aopen must be q or raw\n");
        return 2;
      }
      l2r_p_raw = v == "raw";
    }
    else if (arg == "--l2r-p3-hex") {
      l2r_p3_hex = std::string(next()) == "on";
    }
    else if (arg == "--l2r-p3-lin") {
      l2r_p3_lin = std::string(next()) == "on";
    }
    else if (arg == "--l2r-inc-check") {
      l2r_inc_check = std::string(next()) == "on";
    }
    else if (arg == "--l2r-ninf-mass") {
      l2r_ninf_mass = std::string(next()) == "on";
    }
    else if (arg == "--l2r-keep-ninf") {
      l2r_keep_ninf = std::string(next()) == "on";
    }
    else if (arg == "--l2r-qz-audit") {
      l2r_qz_audit = std::string(next()) == "on";
    }
    else if (arg == "--grow-node-trunc-z") {
      const std::string v = next();
      if (v != "scan" && v != "bar") {
        std::fprintf(stderr, "--grow-node-trunc-z must be scan or bar\n");
        return 2;
      }
      ntr_z_set = true;
      g_ntr_zbar = v == "bar";
    }
    else if (arg == "--l2r-spawn-w") {
      l2r_spawn_wk = std::string(next()) != "suf";
    }
    else if (arg == "--l2r-spawn-pool") {
      l2r_spawn_pool = std::atoi(next());
    }
    else if (arg == "--l2r-spawn-maxpop") {
      l2r_spawn_maxpop = std::atoi(next());
    }
    else if (arg == "--l2r-trace") {
      l2r_trace = next();
    }
    else if (arg == "--l2r-perf") {
      l2r_perf = std::strcmp(next(), "on") == 0;
    }
    else if (arg == "--l2r-trace-key") {
      l2r_trace_key = next();
    }
    else if (arg == "--l2r-oracle-dump") {
      l2r_odump = next();
    }
    else if (arg == "--l2r-oracle") {
      l2r_oload = next();
    }
    else if (arg == "--l2r-sub") {
      l2r_sub = next();
    }
    else if (arg == "--l2r-ctx") {
      const std::string s = next();
      if (s != "off" && s != "on" && s != "p1") {
        std::fprintf(stderr, "--l2r-ctx must be off, on or p1\n");
        return 2;
      }
      if (s != "off") {
        std::fprintf(stderr, "[unsupported] --l2r-ctx %s is off; it builds Bml\n", s.c_str());
        return 2;
      }
      l2r_ctx = false;
      l2r_ctx_p1 = false;
    }
    else if (arg == "--l2r-ctx-tau") {
      l2r_ctx_tau = std::atof(next());
    }
    else if (arg == "--l2r-ctx-rounds") {
      l2r_ctx_rounds = std::max(1, std::atoi(next()));
    }
    else if (arg == "--l2r-kx") {
      l2r_kx = std::string(next()) == "on";
    }
    else if (arg == "--l2r-p1") {
      l2r_p1 = std::atoi(next());
    }
    else if (arg == "--l2r-p1-qb") {
      l2r_p1_own = std::string(next()) == "pass1";
    }
    else if (arg == "--l2r-p1-beam") {
      l2r_p1_beam = std::atoi(next());
    }
    else if (arg == "--l2r-split") {
      l2r_split = std::string(next()) == "on";
    }
    else if (arg == "--pk-dbeta1" || arg == "--pk-dbeta1m" || arg == "--pk-dbeta1p") {
      (arg == "--pk-dbeta1" ? pk_db[0] : arg == "--pk-dbeta1m" ? pk_db[1] : pk_db[2]) = std::atof(next());
    }
    else if (arg == "--pr-score-dump") {
      pr_set = true;
      pr_dump_path = next();
    }
    else if (arg == "--pr-armq-sep") {
      pr_set = true;
      pr_aq_sep = std::string(next()) == "on";
    }
    else if (arg == "--pr-v-ratio") {
      pr_set = true;
      pr_v_ratio = std::string(next()) == "on";
    }
    else if (arg == "--pr-v-diag") {
      pr_set = true;
      pr_v_diag = next();
    }
    else if (arg == "--pr-v-ratio2") {
      pr_set = true;
      g_v2ratio = std::string(next()) == "on";
      if (g_v2ratio) pr_v_ratio = true;
    }
    else if (arg == "--pr-v2-dump") {
      pr_set = true;
      g_v2dump = next();
    }
    else if (arg == "--pr-static-only") {
      pr_set = true;
      g_static_only = true;
    }
    else if (arg == "--pr-v-sample") {
      pr_set = true;
      pr_v_sample = std::strtod(next(), nullptr);
    }
    else if (arg == "--pr-aq-max") {
      pr_set = true;
      pr_aq_max = static_cast<size_t>(std::max(1, std::atoi(next())));
    }
    else if (arg == "--pr-arm-persist" || arg == "--pr-arm-gaponly") {
      pr_set = true;
      const std::string s = next();
      if (s != "on" && s != "off") {
        std::fprintf(stderr, "%s on|off\n", arg.c_str());
        return 2;
      }
      (arg == "--pr-arm-persist" ? pr_arm_persist : pr_arm_gaponly) = s == "on";
    }
    else if (arg == "--pr-partner-maxpairs") {
      pr_set = true;
      const std::string s = next();
      pr_pmax = s == "off" ? 0 : std::atoi(s.c_str());
      if (pr_pmax < 0 || pr_pmax > 255) {
        std::fprintf(stderr, "--pr-partner-maxpairs must be off or 1..255\n");
        return 2;
      }
    }
    else if (arg == "--pr-partner-check") {
      pr_set = true;
      pr_pcheck = true;
    }
    else if (arg == "--pr-partner-ml") {
      pr_set = true;
      const std::string s = next();
      if (s != "on" && s != "off") {
        std::fprintf(stderr, "--pr-partner-ml must be on or off\n");
        return 2;
      }
      pr_pml = s == "on";
    }
    else if (arg == "--pr-fast") {
      pr_set = true;
      const std::string s = next();
      if (s != "off" && s != "acc" && s != "ml") {
        std::fprintf(stderr, "--pr-fast must be off, acc or ml\n");
        return 2;
      }
      pr_fast = s == "off" ? 0 : s == "acc" ? 1 : 2;
    }
    else if (arg == "--pr-abl") {
      pr_set = true;
      std::stringstream ss(next());
      std::string t;
      while (std::getline(ss, t, ',')) {
        if (t == "bml") pr_abl |= DP03::PA_BML;
        else if (t == "q") pr_abl |= DP03::PA_Q;
        else if (t == "p") pr_abl |= DP03::PA_P;
        else if (t == "xs") pr_abl |= DP03::PA_XS;
        else if (t == "ks") pr_abl |= DP03::PA_KS;
        else if (t == "xb") pr_abl |= DP03::PA_XB;
        else if (t == "qr") pr_abl |= DP03::PA_QR;
        else if (t == "sep") pr_abl |= DP03::PA_SEP | DP03::PA_BML | DP03::PA_Q;
        else if (t == "sepv") pr_abl |= DP03::PA_SEP | DP03::PA_BML;
        else if (t != "none") {
          std::fprintf(stderr, "--pr-abl: comma list of bml, q, p, xs, xb, ks, qr, sep, sepv, none\n");
          return 2;
        }
      }
    }
    else if (arg == "--pr-abl-pk") {
      pr_set = true;
      pr_abl_pk = std::stoi(next());
    }
    else if (arg == "--pr-partner-ml-check") {
      pr_set = true;
      pr_pmlcheck = true;
    }
    else if (arg == "--pr-partner-span") {
      pr_set = true;
      pr_pspan = std::atoi(next());
    }
    else if (arg == "--pr-v-eps") {
      pr_set = true;
      pr_veps = std::atof(next());
    }
    else if (arg == "--pr-v-trim") {
      pr_set = true;
      const std::string s = next();
      if (s != "global" && s != "suffix") {
        std::fprintf(stderr, "--pr-v-trim must be global or suffix\n");
        return 2;
      }
      pr_vsuffix = s == "suffix";
    }
    else if (arg == "--pr-v-cap" || arg.rfind("--pr-v-cap=", 0) == 0) {
      std::fprintf(stderr, "--pr-v-cap is removed; use --bV\n");
      return 2;
    }
    else if (arg == "--bV") {
      pr_set = true;
      bV_set = true;
      const int n = std::max(0, std::atoi(next()));
      pr_vcap = n;
      p2_v_cap = n;  // one flag: xs-keyed cap and C_pre left-end cap
    }
    else if (arg == "--pr-v-mlsteps") {
      pr_vmlsteps = std::string(next()) == "on";
    }
    else if (arg == "--pr-v-mlk") {
      pr_vmlk = std::max(0, std::atoi(next()));
    }
    else if (arg == "--pr-v-key") {
      pr_set = true;
      const std::string s = next();
      if (s != "xsys" && s != "xs") {
        std::fprintf(stderr, "--pr-v-key must be xsys or xs\n");
        return 2;
      }
      pr_vxs = s == "xs";
    }
    else if (arg == "--pr-v-stream") {
      pr_set = true;
      pr_vstream = std::string(next()) == "on";
    }
    else if (arg == "--grow-node-trunc" || arg.rfind("--grow-node-trunc=", 0) == 0) {
      // Exact --grow-node-trunc only (not -r / -z).
      std::fprintf(stderr, "--grow-node-trunc is removed; use --bG\n");
      return 2;
    }
    else if (arg == "--bG") {
      g_ntr_K = std::atoi(next());
    }
    else if (arg == "--g-cut") {
      const std::string s = next();
      if (s == "topq") g_cut = G_CUT_TOPQ;
      else if (s == "collapse") g_cut = G_CUT_COLLAPSE;
      else if (s == "entries") g_cut = G_CUT_ENTRIES;
      else {
        std::fprintf(stderr, "--g-cut must be topq, collapse, or entries\n");
        return 2;
      }
    }
    else if (arg == "--g-sweep") {
      const std::string s = next();
      if (s == "direct") g_sweep = G_SWEEP_DIRECT;
      else if (s == "states") g_sweep = G_SWEEP_STATES;
      else {
        std::fprintf(stderr, "--g-sweep must be direct or states\n");
        return 2;
      }
    }
    else if (arg == "--grow-node-trunc-r") {
      g_ntr_KR = std::max(-1, std::atoi(next()));
    }
    else if (arg == "--grow-r-from-v") {
      g_rfv = std::string(next()) == "on";
    }
    else if (arg == "--malloc-trim") {
      g_mtrim = std::string(next()) == "on";
    }
    else if (arg == "--l2r-arm-fb") {
      l2r_arm_fb = std::string(next()) == "on";
    }
    else if (arg == "--l2r-lean") {
      l2r_lean = std::string(next()) == "on";
    }
    else if (arg == "--mem-log") {
      g_mem_log = std::string(next()) == "on";
    }
    else if (arg == "--pr-arms-eps") {
      pr_set = true;
      pr_arms_eps = std::strtod(next(), nullptr);
    }
    else if (arg == "--pr-arms-max") {
      pr_set = true;
      pr_arms_max = std::atoi(next());
    }
    else if (arg == "--pr-oracle") pr_oracle_path = next();
    else if (arg == "--pr-oracle-set") pr_oracle_set = next();
    else if (arg == "--stem-est-ext-max") {
      stem_ext_set = true;
      stem_ext_max = std::atoi(next());
    }
    else if (arg == "--stem-est-qg-combine") {
      std::string s = next();
      stem_ext_set = true;
      stem_split = s == "split";
      if (s != "sum" && s != "split") {
        std::fprintf(stderr, "--stem-est-qg-combine must be 'sum' or 'split'\n");
        return 2;
      }
    }
    else if (arg == "--stem-ext-check") stem_ext_check = true;
    else if (arg == "--stem-bg-check") stem_bg_check = true;
    else if (arg == "--stem-diag") stem_diag_path = next();
    else if (arg == "--stem-diag-loops") {
      stem_diag_loops.clear();
      for (double x : parse_list(next())) stem_diag_loops.push_back(static_cast<int>(x));
    }
    else if (arg == "--gate-diag") gate_diag_path = next();
    else if (arg == "--gate-diag-taus") diag_taus = parse_list(next());
    else if (arg == "--gate-long-span") long_span = std::atoi(next());
    else if (arg == "--beam-table") {
      std::string s = next();
      size_t eq = s.find('=');
      int bt = eq == std::string::npos ? -1 : beam_table_id(s.substr(0, eq));
      std::string b = eq == std::string::npos ? "" : s.substr(eq + 1);
      if (bt < 0 || b.empty()) {
        std::fprintf(stderr, "--beam-table needs T=B with T in Qb,Qm,Qm2,Qmp,PKM,Qp,Qg,Qgls,"
                             "Qgl,Qgr,Qgrs and B a beam or inf\n");
        return 2;
      }
      table_beam[bt] = (b == "inf" || b == "infinity" || b == "none") ? 0
                                                                      : std::max(0, std::atoi(b.c_str()));
    }
    else if (arg == "--dump-states") states_path = next();
    else if (arg == "--dump-scores") scores_path = next();
    else if (arg == "--beam-stats") beam_stats_path = next();
    else if (arg == "--dump-ext") ext_path = next();
    else if (arg == "--dump-min-w") dump_min_w = std::strtod(next(), nullptr);
    else if (arg == "--dump-pruned") pruned_path = next();
    else if (arg == "--dump-kept") kept_path = next();
    else if (arg == "--dump-exact") exact_path = next();
    else if (arg == "--dump-pruned-min-w") pruned_min_w = std::strtod(next(), nullptr);
    else if (arg == "--remove-state") {
      std::string s = next();
      size_t c = s.find(',');
      int bt = c == std::string::npos ? -1 : beam_table_id(s.substr(0, c));
      int v[4] = {-1, -1, -1, -1};
      if (bt >= 0) {
        std::vector<double> xs = parse_list(s.substr(c + 1));
        for (size_t k = 0; k < xs.size() && k < 4; ++k) v[k] = static_cast<int>(xs[k]);
      }
      if (bt < 0) {
        std::fprintf(stderr, "--remove-state needs T,a,b,c,d (keys as in --dump-states)\n");
        return 2;
      }
      remove_keys.insert(state_key(bt, v[0], v[1], v[2], v[3]));
    }
    else if (arg == "--qz-beam" || arg.rfind("--qz-beam=", 0) == 0) {
      std::fprintf(stderr, "--qz-beam is removed; Q^z uses b2 at every length\n");
      return 2;
    }
    else if (arg == "--qz-key" || arg.rfind("--qz-key=", 0) == 0) {
      std::fprintf(stderr, "--qz-key is removed; Q^z is ranked by C(k-1)+Q^z\n");
      return 2;
    }
    else if (arg == "--v-pairs" || arg.rfind("--v-pairs=", 0) == 0) {
      const std::string s = arg == "--v-pairs" ? next() : arg.substr(10);
      if (s == "P1" || s == "p1") {
        p2_v_mode = 1;
        p2_v_bp = 0;
      } else if ((s.rfind("own:", 0) == 0 || s.rfind("inside:", 0) == 0) &&
                 std::atoi(s.c_str() + (s[0] == 'o' ? 4 : 7)) > 0) {
        p2_v_mode = 2;
        p2_v_bp = std::atoi(s.c_str() + (s[0] == 'o' ? 4 : 7));
      } else if (s.rfind("prefix:", 0) == 0 && std::atoi(s.c_str() + 7) > 0) {
        p2_v_mode = 3;
        p2_v_bp = std::atoi(s.c_str() + 7);
      } else {
        std::fprintf(stderr, "--v-pairs must be P1, inside:b_P, prefix:b_P, or own:b_P\n");
        return 2;
      }
    }
    else if (arg == "--v-cap" || arg.rfind("--v-cap=", 0) == 0) {
      std::fprintf(stderr, "--v-cap is removed; use --bV\n");
      return 2;
    }
    else if (arg == "--v-trim" || arg.rfind("--v-trim=", 0) == 0) {
      const std::string s = arg == "--v-trim" ? next() : arg.substr(9);
      if (s == "both") p2_v_trim = 0;
      else if (s == "rel") p2_v_trim = 1;
      else if (s == "cap") p2_v_trim = 2;
      else {
        std::fprintf(stderr, "--v-trim must be rel, cap, or both\n");
        return 2;
      }
    }
    else if (arg == "--g-pairs" || arg.rfind("--g-pairs=", 0) == 0) {
      const std::string s = arg == "--g-pairs" ? next() : arg.substr(10);
      if (s == "P1" || s == "p1") {
        p2_g_mode = 1;
        p2_g_bp = 0;
      } else if ((s.rfind("own:", 0) == 0 || s.rfind("inside:", 0) == 0) &&
                 std::atoi(s.c_str() + (s[0] == 'o' ? 4 : 7)) > 0) {
        p2_g_mode = 2;
        p2_g_bp = std::atoi(s.c_str() + (s[0] == 'o' ? 4 : 7));
      } else if (s.rfind("prefix:", 0) == 0 && std::atoi(s.c_str() + 7) > 0) {
        p2_g_mode = 3;
        p2_g_bp = std::atoi(s.c_str() + 7);
      } else {
        std::fprintf(stderr, "--g-pairs must be P1, inside:b_P, prefix:b_P, or own:b_P\n");
        return 2;
      }
    }
    else if (arg == "--g-mode" || arg.rfind("--g-mode=", 0) == 0) {
      const std::string s = arg == "--g-mode" ? next() : arg.substr(9);
      if (s == "eager") p2_g_eager = true;
      else if (s == "lazy") p2_g_eager = false;
      else {
        std::fprintf(stderr, "--g-mode must be eager or lazy\n");
        return 2;
      }
    }
    else if (arg == "--g-cap" || arg.rfind("--g-cap=", 0) == 0) {
      const std::string s = arg == "--g-cap" ? next() : arg.substr(8);
      p2_g_cap = std::max(0, std::atoi(s.c_str()));
      if (p2_g_cap > 0) {
        g_ntr_K = p2_g_cap;
        p2_g_suffix = true;
      }
    }
    else if (arg == "--g-fallback" || arg.rfind("--g-fallback=", 0) == 0) {
      const std::string s = arg == "--g-fallback" ? next() : arg.substr(13);
      if (s == "onestep" || s == "on") p2_g_fb = true;
      else if (s == "off") p2_g_fb = false;
      else {
        std::fprintf(stderr, "--g-fallback must be onestep or off\n");
        return 2;
      }
    }
    else if (arg == "--profile") {
      profile = true;
    }
    else if (arg == "--g-mass" || arg.rfind("--g-mass=", 0) == 0) {
      const std::string s = arg == "--g-mass" ? next() : arg.substr(9);
      g_mass_f = std::fopen(s.c_str(), "wb");
      if (!g_mass_f) {
        std::fprintf(stderr, "--g-mass cannot open %s\n", s.c_str());
        return 2;
      }
    }
    else if (arg == "--g-key" || arg.rfind("--g-key=", 0) == 0) {
      const std::string s = arg == "--g-key" ? next() : arg.substr(8);
      if (s != "suffix") {
        std::fprintf(stderr, "--g-key must be suffix\n");
        return 2;
      }
      p2_g_suffix = true;
    }
    else if (arg == "--sweep-check") sw_check = true;
    else if (arg == "--sweep") {
      sweep_set = true;
      if (a + 1 < argc && (std::string(argv[a + 1]) == "off" || std::string(argv[a + 1]) == "on"))
        p2_sweep = std::string(argv[++a]) == "on";
      else
        p2_sweep = true;
    }
    else if (arg == "--v-sweep") {
      vsweep_set = true;
      if (a + 1 < argc && (std::string(argv[a + 1]) == "off" || std::string(argv[a + 1]) == "on"))
        p2_vsweep = std::string(argv[++a]) == "on";
      else
        p2_vsweep = true;
    }
    else if (arg == "--k-memo") k_memo = true;
    else if (arg == "--k-exact") k_exact = true;
    else if (arg == "--k-bound") k_bound = true;
    else if (arg == "--k-exact-ub") k_exact_ub = true;
    else if (arg == "--k-exact-check") k_exact_check = true;
    else if (arg == "--rlt") rlt = true;
    else if (arg == "--pass1-only") pass1_only = true;
    else if (arg == "--c-trace") c_trace = true;
    else if (arg == "--grammar") {
      if (a + 1 >= argc) {
        std::fprintf(stderr, "--grammar needs full or simple\n");
        return 2;
      }
      const std::string g = argv[++a];
      if (g == "simple") grammar_simple = true;
      else if (g == "full") grammar_simple = false;
      else {
        std::fprintf(stderr, "--grammar must be full or simple\n");
        return 2;
      }
    }
    else if (arg == "--next-pair") next_pair = true;
    else if (arg == "--rlt-beam") {
      rlt = true;
      rlt_beam = std::max(0, std::atoi(next()));
      rlt_beam_set = true;
    }
    else if (arg == "--sample-prof") sample_prof = next();
    else if (arg == "--v-floor" || arg.rfind("--v-floor=", 0) == 0) {
      const std::string s = arg == "--v-floor" ? std::string("min") : arg.substr(10);
      if (s == "min") p2_v_floor = 1;
      else if (s == "0") p2_v_floor = 2;
      else {
        std::fprintf(stderr, "--v-floor must be min or 0\n");
        return 2;
      }
    }
    else if (arg == "--v-cut-fb") {
      if (a + 1 < argc && (std::string(argv[a + 1]) == "off" || std::string(argv[a + 1]) == "on"))
        p2_v_cut_fb = std::string(argv[++a]) == "on";
      else
        p2_v_cut_fb = true;
    }
    else if (arg == "--mu-dump") {
      g_mu_dump = next();
    }
    else if (arg == "--g-eq") {
      g_eq_path = next();
    }
    else if (arg == "--miss-sample") {
      g_miss_cap = std::max(0, std::atoi(next()));
      g_miss_path = next();
    }
    else if (arg == "--v-query") {
      const std::string in_path = next();
      g_vquery_out = next();
      if (a + 1 < argc && std::string(argv[a + 1]).rfind("--", 0) != 0) g_vtrace_path = next();
      std::ifstream vin(in_path);
      std::string line;
      while (std::getline(vin, line)) {
        if (line.empty() || line[0] == 's') continue;
        std::istringstream ls(line);
        std::string tag;
        int st, x, y, b;
        if (!(ls >> tag >> st >> x >> y >> b)) continue;
        int kind = tag == "K" ? 1 : tag == "H" ? 2 : 0;
        g_vquery.push_back({kind, st, x, y, b});
      }
    }
    else if (arg == "--lazy-ml-cap") {
      if (a + 1 < argc && (std::string(argv[a + 1]) == "off" || std::string(argv[a + 1]) == "on"))
        lazy_ml_cap = std::string(argv[++a]) == "on";
      else
        lazy_ml_cap = true;
    }
    else if (arg == "--nest-beam" || arg.rfind("--nest-beam=", 0) == 0) {
      const std::string s = arg == "--nest-beam" ? next() : arg.substr(12);
      nest_beam = std::max(0, std::atoi(s.c_str()));
      b1_set = true;
    }
    else if (arg == "--b1") {
      nest_beam = std::max(0, std::atoi(next()));
      b1_set = true;
    }
    else if (arg == "--b2") {
      const int nbeam = std::max(0, std::atoi(next()));
      table_beam[BT_QB] = table_beam[BT_QM] = table_beam[BT_QM2] = table_beam[BT_QMP] =
          table_beam[BT_QP] = nbeam;
      b2_at = a;
    }
    else if (arg == "--qzbar-old") {
      qzbar_old = true;
    }
    else if (arg == "--pbeam-ratio") {
      int ra, rb, rc;
      if (std::sscanf(next(), "%d:%d:%d", &ra, &rb, &rc) != 3 || ra < 0 || rb < 0 || rc < 0) {
        std::fprintf(stderr, "--pbeam-ratio needs a:b:c of non-negative integers\n");
        return 2;
      }
      pbr[0] = ra;
      pbr[1] = rb;
      pbr[2] = rc;
      pbeam_ratio_at = a;
    }
    else if (arg == "--pre-beam" || arg.rfind("--pre-beam=", 0) == 0) {
      const std::string s = arg == "--pre-beam" ? next() : arg.substr(11);
      pre_beam = std::max(0, std::atoi(s.c_str()));
    }
    else if (arg == "--pre-key" || arg.rfind("--pre-key=", 0) == 0) {
      const std::string s = arg == "--pre-key" ? next() : arg.substr(10);
      if (s == "prefix") pre_key = 0;
      else if (s == "inside") pre_key = 1;
      else {
        std::fprintf(stderr, "--pre-key must be prefix or inside\n");
        return 2;
      }
    }
    else if (arg == "--dump-qb") dump_qb_on = true;
    else if (arg == "-h" || arg == "--help") {
      usage();
      return 0;
    } else {
      usage();
      return 2;
    }
  }

  // --pk-l2r on defaults: --l2r-spawn inc; pbeam derived below as (b2,b4,b4) unless
  // --l2r-pbeam / --pbeam-ratio given; --l2r-p3-hex/lin on; --pr-ml-span 100000; with
  // --grow-rec on (ml all/all) --bG 16 --grow-node-trunc-r -1; --pr-v-key xs when safe;
  // then --grow-r-from-v on if its requirements hold.
  if (pk_l2r == 1) {
    auto given = [&](const char* f) {
      for (int k = 1; k < argc; ++k)
        if (std::strcmp(argv[k], f) == 0) return true;
      return false;
    };
    if (!given("--l2r-spawn")) {
      l2r_spawn_inc = true;
      l2r_spawn_cube = l2r_spawn_static = false;
    }
    // pbeam filled after b2/b4 are known (see below).
    if (!given("--l2r-p3-hex")) l2r_p3_hex = true;
    if (!given("--l2r-p3-lin")) l2r_p3_lin = true;
    if (!given("--pr-ml-span")) pr_ml_span = 100000;
    const bool grow_ok = g_grow_rec && g_ml_L == 2 && g_ml_R == 2 && !g_ml_diag && !g_rec_check;
    if (grow_ok && !given("--bG")) g_ntr_K = 16;
    if (p2_g_cap > 0) g_ntr_K = p2_g_cap;
    if (grow_ok && !given("--grow-node-trunc-r")) g_ntr_KR = -1;
    if (!given("--pr-v-key") && pr_fast > 0 && !l2r_ctx && !l2r_ctx_p1 && l2r_p1 == 0 && pr_pmax == 0 &&
        !pr_v_ratio)
      pr_vxs = true;
    if (!given("--grow-r-from-v")) {
      if (g_ntr_K > 0 && g_ntr_KR < 0 && pr_vxs && pr_fast >= 2 && !l2r_ctx_p1 && !l2r_diag)
        g_rfv = true;
      else
        std::fprintf(stderr, "[defaults] --grow-r-from-v left off: needs --bG K > 0, "
                             "--grow-node-trunc-r -1, --pr-v-key xs, --pr-fast ml\n");
    }
  }
  if (pbeam_ratio_at > pbeam_explicit_at) {
    const int b4 = beam > 0 ? beam : 0;
    for (int k = 0; k < 3; ++k) l2r_pbeam[k] = pbr[k] * b4;
  }
  // b1 and b2 default to 200, independent of b4. rlt-beam defaults to b4.
  // Sweep and the V sweep are on. No length gate. Q^z and Q^z_pre use b2.
  if (!b1_set) nest_beam = 200;
  if (b2_at < 0) {
    const int b2tabs[] = {BT_QB, BT_QM, BT_QM2, BT_QMP, BT_QP};
    for (int t : b2tabs)
      if (table_beam[t] < 0) table_beam[t] = 200;
  }
  b2_width = table_beam[BT_QM] >= 0 ? table_beam[BT_QM] : 200;
  // pbeam = (b2, b4, b4) unless --l2r-pbeam or --pbeam-ratio was set.
  if (!user_pbeam && !user_pratio) {
    l2r_pbeam[0] = b2_width;
    l2r_pbeam[1] = beam;
    l2r_pbeam[2] = beam;
  }
  // --bV default 256 (replaces --v-cap / --pr-v-cap; 0 = uncapped).
  if (!bV_set) {
    pr_vcap = 256;
    p2_v_cap = 256;
  }
  std::fprintf(stderr, "[arity] b2 %d b4 %d pbeam %d,%d,%d\n", b2_width, beam, l2r_pbeam[0],
               l2r_pbeam[1], l2r_pbeam[2]);
  if (!rlt_beam_set) rlt_beam = beam;
  // Default sweep ON: sw_as_rn with q'-aware intermediate cap (= lazy trunc key).
  // Pass `--sweep off` for rec_node/pr_steps (hard14/OLD path).
  if (!sweep_set) p2_sweep = true;
  if (!vsweep_set) p2_vsweep = true;

  if (arm_score && access_file.empty()) {
    std::fprintf(stderr, "--arm-score on needs --access-file\n");
    return 2;
  }
  if (arm_score && !score_qz) {
    std::fprintf(stderr, "--arm-score on needs --score qz\n");
    return 2;
  }
  const bool stem_diag = !stem_diag_path.empty();
  if ((stem_rank || stem_diag) && (!allow_pk || arm_score || gcfg.on)) {
    std::fprintf(stderr, "--stem-est rank and --stem-diag need pseudoknots, --arm-score off and "
                         "--gate off\n");
    return 2;
  }
  const bool dump_scores = !scores_path.empty();
  const bool dump_ext = !ext_path.empty();
  bool all_inf = beam <= 0;
  for (int t = 0; t < BT_N; ++t) all_inf &= table_beam[t] <= 0;
  if (dump_ext && (!all_inf || !gap_step || stem_rank || stem_diag || dump_scores || !allow_pk ||
                   arm_score || gcfg.on || exact_path.empty())) {
    std::fprintf(stderr, "--dump-ext runs at --beam inf (every table) in step mode with "
                         "pseudoknots and --dump-exact, --stem-est off, --arm-score off, "
                         "--gate off, no --stem-diag / --dump-scores\n");
    return 2;
  }
  if (dump_scores && (!all_inf || !gap_step || stem_rank || stem_diag || !allow_pk ||
                      arm_score || gcfg.on)) {
    std::fprintf(stderr, "--dump-scores runs at --beam inf (every table) in step mode with "
                         "pseudoknots, --stem-est off, --arm-score off, --gate off, no "
                         "--stem-diag\n");
    return 2;
  }
  if (stem_diag && (stem_rank || beam > 0 || !states_path.empty() || !gate_diag_path.empty() ||
                    !do_bpp)) {
    std::fprintf(stderr, "--stem-diag runs at --beam inf with --stem-est off, the outside pass "
                         "on, and without --dump-states / --gate-diag\n");
    return 2;
  }
  for (int t = 0; t < BT_N; ++t)
    if (stem_diag && table_beam[t] > 0) {
      std::fprintf(stderr, "--stem-diag needs every table at infinite beam\n");
      return 2;
    }
  if (stem_loop < 0) {
    std::fprintf(stderr, "--stem-est-loop must be >= 0\n");
    return 2;
  }
  if (!pr_oracle_path.empty() && (!stem_pr || !all_inf || !gap_step || !allow_pk)) {
    std::fprintf(stderr, "--pr-oracle needs --stem-est-principled on, --beam inf (every "
                         "table), step mode and pseudoknots\n");
    return 2;
  }
  if ((stem_qg_r_set || stem_consts_set) && !stem_rank) {
    std::fprintf(stderr, "--stem-est-qg-r and --stem-est-consts need --stem-est rank\n");
    return 2;
  }
  if (stem_qg_partner && !stem_qg2) {
    std::fprintf(stderr, "--stem-est-qg-r partner needs --stem-est-role all\n");
    return 2;
  }
  if ((stem_ext_set && !stem_rank) || (stem_ext && !stem_qg_partner)) {
    std::fprintf(stderr, "--stem-est-qg-ext needs --stem-est rank, --stem-est-role all and "
                         "--stem-est-qg-r partner\n");
    return 2;
  }
  if (stem_ext_max < 1) {
    std::fprintf(stderr, "--stem-est-ext-max must be >= 1\n");
    return 2;
  }
  if ((stem_split || stem_ext_check) && !stem_ext) {
    std::fprintf(stderr, "--stem-est-qg-combine split and --stem-ext-check need --stem-est-qg-ext\n");
    return 2;
  }
  if (grow_set && stem_ext != 3) {
    std::fprintf(stderr, "--stem-est-grow-loop / -cut / -nodes need --stem-est-qg-ext loop\n");
    return 2;
  }
  if (grow_lg < 0 || grow_cut < 0.0 || grow_cut > 1.0 || grow_nodes < 0 || ml_span < 0) {
    std::fprintf(stderr, "--stem-est-grow-loop, -nodes, --stem-est-ml-span must be >= 0 and "
                         "--stem-est-grow-cut in [0, 1]\n");
    return 2;
  }
  if (helpers_set && !stem_ext) {
    std::fprintf(stderr, "--stem-est-helpers and --stem-est-ml-span need --stem-est-qg-ext\n");
    return 2;
  }
  if (stem_profile && (stem_ext != 3 || stem_helpers != 1 || stem_split || stem_ext_check)) {
    std::fprintf(stderr, "--stem-est-profile on needs --stem-est-qg-ext loop, --stem-est-helpers "
                         "grow, --stem-est-qg-combine sum and no --stem-ext-check\n");
    return 2;
  }
  if ((profile_set || prof_check) && !stem_profile) {
    std::fprintf(stderr, "--stem-est-profile-eps, --stem-est-ml-tail and --stem-profile-check "
                         "need --stem-est-profile on\n");
    return 2;
  }
  if (!(prof_eps >= 0.0 && prof_eps < 1.0)) {
    std::fprintf(stderr, "--stem-est-profile-eps must be in [0, 1)\n");
    return 2;
  }
  if ((bg_reorder && !stem_profile) || (ro_check && !bg_reorder)) {
    std::fprintf(stderr, "--stem-est-bg-reorder on needs --stem-est-profile on; "
                         "--stem-bg-reorder-check needs --stem-est-bg-reorder on\n");
    return 2;
  }
  if ((chain_set && !stem_chain) || (stem_chain && (!stem_profile || ml_tail))) {
    std::fprintf(stderr, "--stem-est-chain-* need --stem-est-chain on, which needs "
                         "--stem-est-profile on and --stem-est-ml-tail off\n");
    return 2;
  }
  if (!(chain_floor >= 0.0 && chain_floor < 1.0) || !(chain_cut >= 0.0 && chain_cut < 1.0) ||
      chain_span < 0 || chain_loop < 0) {
    std::fprintf(stderr, "--stem-est-chain-floor / -cut must be in [0, 1), -span / -loop >= 0\n");
    return 2;
  }
  if ((stem_hr_set || stem_hr_check) && !stem_rank) {
    std::fprintf(stderr, "--stem-est-hr, --stem-est-hr-gamma and --stem-hr-check need "
                         "--stem-est rank\n");
    return 2;
  }
  if (g_grow_rec && (!stem_pr || !(pr_abl & DP03::PA_SEP))) {
    std::fprintf(stderr, "--grow-rec on needs --stem-est-principled on and --pr-abl sep or sepv\n");
    return 2;
  }
  if (g_grow_rec && (g_ml_c < 0 || pr_check || pr_gcheck)) {
    std::fprintf(stderr, "--grow-rec on: --grow-ml-c >= 0, no --stem-pr-check / --pr-grow-check\n");
    return 2;
  }
  if (g_rec_check && (!g_grow_rec || g_ml_L != 2 || g_ml_R != 2 || g_ml_c != 0 ||
                      (pr_gsel && (pr_geps > 0.0 || pr_gB > 0)))) {
    std::fprintf(stderr, "--grow-rec-check needs --grow-rec on, ml all/all, no --grow-ml-c and "
                         "unselected profiles (--pr-grow-eps 0 --pr-grow-B 0)\n");
    return 2;
  }
  if (g_tr_K < 0 || !(g_tr_eps >= 0.0 && g_tr_eps < 1.0) ||
      (g_tr_K > 0 && (!g_grow_rec || !pr_gsel || g_ml_diag || g_rec_check || pr_abl_pk > 0))) {
    std::fprintf(stderr, "--grow-trunc K >= 0 needs --grow-rec on and --pr-grow-select bound, "
                         "no --grow-ml-diag / --grow-rec-check / --pr-abl-pk; eps in [0, 1)\n");
    return 2;
  }
  if (g_ntr_K < 0 || (g_ntr_K > 0 && (!g_grow_rec || pk_l2r != 1 || g_ml_L != 2 || g_ml_R != 2 ||
                                      g_ml_diag || g_rec_check))) {
    std::fprintf(stderr, "--bG K >= 0 needs --grow-rec on, --pk-l2r on, ml all/all, "
                         "no --grow-ml-diag / --grow-rec-check\n");
    return 2;
  }
  if (pr_vxs && (pk_l2r != 1 || l2r_ctx || l2r_ctx_p1 || l2r_p1 > 0 || pr_fast <= 0 || pr_pmax > 0 ||
                 pr_v_ratio)) {
    std::fprintf(stderr, "--pr-v-key xs needs --pk-l2r on and --pr-fast, no --l2r-ctx / --l2r-p1 / "
                         "--pr-partner-maxpairs / --pr-v-ratio (they read per-ys arm entries)\n");
    return 2;
  }
  if (g_rfv && (g_ntr_K <= 0 || g_ntr_KR >= 0 || !pr_vxs || pr_fast < 2 || l2r_ctx_p1 || l2r_diag)) {
    std::fprintf(stderr, "--grow-r-from-v on needs --bG K > 0, --grow-node-trunc-r -1, "
                         "--pr-v-key xs, --pr-fast ml; no --l2r-ctx p1 / --l2r-diag\n");
    return 2;
  }
  if (pr_vstream && (pr_fast <= 0 || pr_pmax > 0 || pr_v_ratio || pr_pmlcheck)) {
    std::fprintf(stderr, "--pr-v-stream on needs --pr-fast, no --pr-partner-maxpairs / --pr-v-ratio / "
                         "--pr-partner-ml-check\n");
    return 2;
  }
  g_skip_bml = (g_grow_rec || g_U_nobml) && (pr_abl & DP03::PA_BML) && !pr_check &&
               !(pk_l2r == 1 && l2r_ctx);
  if (stem_pr && !g_skip_bml) {
    std::fprintf(stderr, "[unsupported] principled score with Bml is off\n");
    g_skip_bml = true;
  }
  PK.beta1 += pk_db[0];
  PK.beta1m += pk_db[1];
  PK.beta1p += pk_db[2];
  if (!std::isfinite(stem_hr_gamma) || stem_hr_gamma < 0.0 ||
      (!stem_hr_phi && !(stem_hr_gamma > 0.0))) {
    std::fprintf(stderr, "--stem-est-hr-gamma must be finite and > 0 (>= 0 with phi)\n");
    return 2;
  }
  if (stem_hr_check && !stem_hr_av && !stem_hr_phi) {
    std::fprintf(stderr, "--stem-hr-check needs --stem-est-hr av or phi\n");
    return 2;
  }
  if ((stem_pr || stem_pr_nested) && (!stem_rank || !stem_qg2 || !stem_both)) {
    std::fprintf(stderr, "--stem-est-principled(-nested) on needs --stem-est rank and "
                         "--stem-est-role all\n");
    return 2;
  }
  if (stem_pr && (stem_split || stem_ext_check || stem_hr_av || stem_hr_phi || stem_hr_check ||
                  ro_check || prof_check)) {
    std::fprintf(stderr, "--stem-est-principled on excludes --stem-est-qg-combine split, "
                         "--stem-est-hr av / phi and the other stem checks\n");
    return 2;
  }
  if ((pr_set || pr_check) && !stem_pr && !stem_pr_nested) {
    std::fprintf(stderr, "--pr-* and --stem-pr-check need --stem-est-principled on\n");
    return 2;
  }
  if (pr_arms_mode == 1) {
    if (!pr_cut_by_set) pr_gapfree = true;
    if (pr_partner_floor < 0.0) pr_partner_floor = 1e-20;
  }
  if ((pr_arm_z == 1 || pr_sep_check) && pr_arms_mode != 1) {
    std::fprintf(stderr, "--pr-arm-z and --pr-sep-check need --pr-arms sep\n");
    return 2;
  }
  if (!(pr_geps >= 0.0 && pr_geps < 1.0) || pr_gB < 0 || (pr_gsel && pr_cut_roots)) {
    std::fprintf(stderr, "--pr-grow-eps must be in [0, 1), --pr-grow-B >= 0 (0 = no cap); "
                         "--pr-grow-select bound excludes --pr-cut-roots on\n");
    return 2;
  }
  if ((pr_arm_persist || pr_arm_gaponly) &&
      (pr_arms_k <= 0 || pr_ratio || pr_arms_mode == 1 || pr_gcheck)) {
    std::fprintf(stderr, "--pr-arm-persist / --pr-arm-gaponly need --pr-arms-k > 0, --pr-arms-rank "
                         "phi, no --pr-arms sep and no --pr-grow-check\n");
    return 2;
  }
  if ((pr_arm_persist && !pr_aq_sep) || pr_pspan < 0 || !(pr_veps >= 0.0 && pr_veps < 1.0)) {
    std::fprintf(stderr, "--pr-arm-persist needs --pr-armq-sep on; --pr-partner-span >= 0; "
                         "--pr-v-eps in [0, 1)\n");
    return 2;
  }
  if (pr_gsel) pr_floor = 0.0;
  FILE* pr_dump_f = nullptr;
  if (!pr_dump_path.empty() && !(pr_dump_f = std::fopen(pr_dump_path.c_str(), "w"))) {
    std::fprintf(stderr, "cannot write %s\n", pr_dump_path.c_str());
    return 2;
  }
  if (!(pr_arms_eps >= 0.0) || pr_arms_max < 1 || (pr_cut_roots && !pr_gapfree)) {
    std::fprintf(stderr, "--pr-arms-eps must be >= 0, --pr-arms-max >= 1; --pr-cut-roots on "
                         "needs --pr-cut-by gapfree\n");
    return 2;
  }
  if ((p2_qz_beam > 0 || p2_v_mode > 0 || p2_g_mode > 0 || p2_v_cap > 0 || p2_g_suffix ||
       p2_g_eager || p2_g_fb || p2_g_cap > 0) &&
      !stem_pr) {
    std::fprintf(stderr, "--qz-beam / --v-pairs / --v-cap / --g-pairs / --g-key need "
                         "--stem-est-principled on\n");
    return 2;
  }
  if (!(pr_floor >= 0.0 && pr_floor < 1.0) || !(pr_cut >= 0.0 && pr_cut < 1.0) || pr_span < 1 ||
      pr_ml_span < 0 || pr_ml_close < 0 || pr_arms_k < 0) {
    std::fprintf(stderr, "--pr-floor / --pr-cut must be in [0, 1), --pr-span >= 1, "
                         "--pr-ml-span / --pr-ml-close / --pr-arms-k >= 0\n");
    return 2;
  }
  if (stem_bg_check && (!stem_rank || !all_inf || !do_bpp)) {
    std::fprintf(stderr, "--stem-bg-check needs --stem-est rank, --beam inf (every table) and "
                         "the outside pass\n");
    return 2;
  }
  FILE* stem_diag_f = nullptr;
  if (stem_diag && !(stem_diag_f = std::fopen(stem_diag_path.c_str(), "w"))) {
    std::fprintf(stderr, "cannot write %s\n", stem_diag_path.c_str());
    return 2;
  }
  const bool gate_diag = !gate_diag_path.empty();
  const bool need_access = arm_score || gcfg.on || gate_diag;
  if ((gcfg.on || gate_diag) && access_file.empty()) {
    std::fprintf(stderr, "--gate on and --gate-diag need --access-file\n");
    return 2;
  }
  if ((gcfg.on || gate_diag) && !allow_pk) {
    std::fprintf(stderr, "--gate on and --gate-diag need pseudoknots (no --no-pk)\n");
    return 2;
  }
  AccessIndex access_idx;
  if (need_access) access_idx = index_access_file(access_file);
  FILE* diag_f = nullptr;
  if (gate_diag) {
    diag_f = std::fopen(gate_diag_path.c_str(), "w");
    if (!diag_f) {
      std::fprintf(stderr, "cannot write %s\n", gate_diag_path.c_str());
      return 2;
    }
  }
  if ((dangles_mode == 1 || dangles_mode == 2) && !mat_rna1999) {
    std::fprintf(stderr, "--dangles some|all needs --material rna1999\n");
    return 2;
  }
  if (dangles_mode == 3 && !mat_v04) {
    std::fprintf(stderr, "--dangles d2 needs --material vienna2004\n");
    return 2;
  }
  if (mat_rna1999 || mat_v04) {
    if (param_file.empty()) {
      char buf[4096];
      ssize_t len = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
      std::string exe = len > 0 ? std::string(buf, len) : std::string(argv[0]);
      std::string dir = exe.substr(0, exe.find_last_of('/') + 1);
      param_file = dir + (mat_v04 ? "../dp03/params/vienna2004.par" : "../dp03/params/rna1999.dG");
    }
    if (mat_v04) load_v04(param_file);
    else load_rna1999(param_file);
  }
  if (mat_v04 && !loop_cap_set) loop_cap = V04_MAX_LOOP;
  // Resolve per-pass caps (temporary --loop-cap-p1/p2/scan). Default: --loop-cap.
  if (loop_cap_p1 < 0) loop_cap_p1 = loop_cap;
  if (loop_cap_p2 < 0) loop_cap_p2 = loop_cap;
  if (loop_cap_scan < 0) loop_cap_scan = loop_cap;
  auto apply_max_loop = [&](int lc) {
    // Keep energy max_loop in lockstep with the DP window for this pass.
    max_loop = (lc >= NO_LOOP_CAP) ? NO_LOOP_CAP : lc;
  };
  if (loop_cap_p1 != loop_cap || loop_cap_p2 != loop_cap || loop_cap_scan != loop_cap)
    std::fprintf(stderr, "[loop-cap-pass] p1 %d p2 %d scan %d (base %d)\n", loop_cap_p1,
                 loop_cap_p2, loop_cap_scan, loop_cap);
  if (!pruned_path.empty() && gcfg.on) {
    std::fprintf(stderr, "--dump-pruned does not log the gated pruning (needs --gate off)\n");
    return 2;
  }
  if (!exact_path.empty() && pruned_path.empty() && !dump_ext) {
    std::fprintf(stderr, "--dump-exact needs --dump-pruned or --dump-ext\n");
    return 2;
  }
  FILE* states_f = nullptr;
  if (!states_path.empty() && !(states_f = std::fopen(states_path.c_str(), "w"))) {
    std::fprintf(stderr, "cannot write %s\n", states_path.c_str());
    return 2;
  }
  FILE* oracle_f = nullptr;
  if (!pr_oracle_path.empty() && !(oracle_f = std::fopen(pr_oracle_path.c_str(), "w"))) {
    std::fprintf(stderr, "cannot write %s\n", pr_oracle_path.c_str());
    return 2;
  }
  FILE* pruned_f = nullptr;
  if (!pruned_path.empty() && !(pruned_f = std::fopen(pruned_path.c_str(), "w"))) {
    std::fprintf(stderr, "cannot write %s\n", pruned_path.c_str());
    return 2;
  }
  if (!kept_path.empty() && pruned_path.empty()) {
    std::fprintf(stderr, "--dump-kept needs --dump-pruned\n");
    return 2;
  }
  FILE* kept_f = nullptr;
  if (!kept_path.empty() && !(kept_f = std::fopen(kept_path.c_str(), "w"))) {
    std::fprintf(stderr, "cannot write %s\n", kept_path.c_str());
    return 2;
  }
  FILE* scores_f = nullptr;
  if (!beam_stats_path.empty()) {
    g_beam_stats.open(beam_stats_path);
    if (!g_beam_stats.f) {
      std::fprintf(stderr, "--beam-stats: cannot open %s\n", beam_stats_path.c_str());
      return 2;
    }
  }
  if (dump_scores && !(scores_f = std::fopen(scores_path.c_str(), "w"))) {
    std::fprintf(stderr, "cannot write %s\n", scores_path.c_str());
    return 2;
  }
  FILE* ext_f = nullptr;
  if (dump_ext && !(ext_f = std::fopen(ext_path.c_str(), "w"))) {
    std::fprintf(stderr, "cannot write %s\n", ext_path.c_str());
    return 2;
  }
  std::unordered_map<std::string, std::unordered_map<uint64_t, double>> exact_states;
  if (!exact_path.empty()) exact_states = load_states(exact_path);
  const std::unordered_map<uint64_t, double> no_states;

  // Every allocated table, printed before any DP. A table with 2 or more indices and no
  // positive beam aborts the run.
  auto beam_audit = [&](int n) -> int {
    int bad = 0;
    auto row = [&](const char* stage, const char* name, int idx, const char* src, int w) {
      const int active = w > 0;
      if (idx == 2) {
        std::fprintf(stderr, "[beams] %s %s idx %d beam %s %d active %d binds %d\n", stage, name,
                     idx, src, w, active, active && w < n);
      } else {
        std::fprintf(stderr, "[beams] %s %s idx %d beam %s %d active %d\n", stage, name, idx, src,
                     w, active);
      }
      // Width 0 keeps every state (beam off). A negative width is an unset table.
      if (idx >= 2 && w < 0) {
        std::fprintf(stderr, "[beams] abort: %s %s has %d indices and no active beam\n", stage,
                     name, idx);
        bad = 1;
      }
    };
    auto tw = [&](int t) { return table_beam[t] >= 0 ? table_beam[t] : beam; };
    auto nest_tables = [&](const char* stage, int b1) {
      row(stage, "Qb", 2, "b1", b1);
      row(stage, "Qm", 2, "b1", b1);
      row(stage, "Qm2", 2, "b1", b1);
      row(stage, "Qp", 2, "b1", b1);
      row(stage, "Qmp", 2, "b1", b1);
      row(stage, "Qz", 2, "b2", b2_width);
      row(stage, "C", 1, "none", 0);
    };
    const int bV = p2_v_cap > 0 ? p2_v_cap : pr_vcap;
    const char* bVsrc = "bV";
    std::fprintf(stderr, "[beams] n %d b1 %d b2 %d b4 %d bV %d bG %d pbeam %d,%d,%d\n", n,
                 nest_beam, b2_width, beam, bV, g_ntr_K, l2r_pbeam[0], l2r_pbeam[1],
                 l2r_pbeam[2]);
    if (stem_rank) nest_tables("stem", nest_beam);
    if (stem_pr) {
      nest_tables("pass1", pre_beam > 0 ? pre_beam : nest_beam);
      row("pass2", "V", 3, bVsrc, bV);
      if (g_ntr_K > 0)
        row("pass2", "G.L", 3, "b_G", g_ntr_K);
      else
        std::fprintf(stderr, "[beams] pass2 G.L idx 3 beam b_G 0 active 0\n");
      if (g_rfv)
        row("pass2", "G.R", 3, bVsrc, bV);
      else if (g_ntr_KR > 0)
        row("pass2", "G.R", 3, "b_G", g_ntr_KR);
      else if (g_ntr_K > 0)
        row("pass2", "G.R", 3, "b_G", g_ntr_K);
      else
        std::fprintf(stderr, "[beams] pass2 G.R idx 3 beam none 0 active 0\n");
      row("pass2", "Qz_pre", 2, "b2", b2_width);
    }
    row("scan", "Qb", 2, "b2", tw(BT_QB));
    row("scan", "Qm", 2, "b2", tw(BT_QM));
    row("scan", "Qm2", 2, "b2", tw(BT_QM2));
    row("scan", "Qp", 2, "b2", tw(BT_QP));
    row("scan", "Qmp", 2, "b2", tw(BT_QMP));
    row("scan", "Qz", 2, "b2", b2_width);
    row("scan", "C", 1, "none", 0);
    if (allow_pk) {
      row("scan", "PKM", 2, "b4", tw(BT_PKM));
      for (const char* nm : {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us", "Uls", "Ulss", "Urs", "S"})
        row("scan", nm, 4, "b4", beam);
      for (int t : {BT_QG, BT_QGLS, BT_QGL, BT_QGR, BT_QGRS})
        row("scan", BT_NAMES[t], 4, "b4", tw(t));
      if (rlt) {
        row("scan", "Rlt", 4, "rlt", rlt_beam);
        row("scan", "Ult", 4, "rlt", rlt_beam);
      }
    }
    if (pk_l2r == 1) {
      row("scan", "P1", 2, "b2", l2r_pbeam[0]);
      row("scan", "P2", 3, grammar_simple ? "b4" : "2b4", l2r_pbeam[1]);
      row("scan", "P3", 4, grammar_simple ? "b4" : "2b4", l2r_pbeam[2]);
      row("scan", "H", 3, "lp_hk", l2r_p_hk);
      if (l2r_p1 > 0) row("p1sweep", "states", 4, "p1b", l2r_p1_beam);
      if (grammar_simple) {
        for (const char* nm : {"R^M", "R^M2", "R^Multi", "U^M", "U^M2", "U^Multi"})
          row("scan", nm, 4, "b4", beam);
      }
    }
    auto dev = [&](const char* name, int on, const std::string& what) {
      std::fprintf(stderr, "[deviation] %s active %d beam 0 %s\n", name, on, what.c_str());
    };
    dev("max_gap_scan", arm_score ? 1 : 0,
        "arm-score window " + std::to_string(max_gap_scan) + "; unused when --arm-score is off");
    dev("max_arm_len", max_arm_len > 0 ? 1 : 0,
        "longest arm window; 0 means the sequence length");
    dev("gap_step", gap_step ? 1 : 0,
        "groups the b4 beam per step j instead of per (v,j); not a separate count");
    dev("gate", gcfg.on ? 1 : 0,
        "drops Qg/Qgl candidates by an access threshold; not a count by key");
    dev("arm_score", arm_score ? 1 : 0,
        "replaces the arm weight with an access score; not a count by key");
    dev("access", (arm_score || gcfg.on) ? 1 : 0,
        "access file for arm_score and the gate; not loaded in the L2R default");
    dev("loop_cap", 1, "energy-model max unpaired loop (" + std::to_string(loop_cap) +
                           "); not a count by key");
    dev("ml_lead_cap", ml_lead_cap > 0 ? 1 : 0,
        "leading unpaired run before the first multiloop branch is at most " +
            std::to_string(ml_lead_cap > 0 ? ml_lead_cap - 1 : n) +
            " bases; not a beam");
    dev("pair_bonus", pair_bonus != 0.0 ? 1 : 0,
        "added to every pair weight (value " + std::to_string(pair_bonus) + ")");
    dev("lK", k_exact ? 0 : 1, "K uses the ratio lK, not a bound");
    dev("V_miss", p2_vsweep ? 0 : 1,
        "a V miss falls back to the old V builder; with --v-sweep a miss is -inf");
    dev("v-floor", p2_v_floor ? 1 : 0,
        "fills a missing V entry with the smallest kept mu or log-weight 0");
    dev("v-cut-fb", p2_v_cut_fb ? 1 : 0,
        "a cut V entry (full map) is scored kappa_min / C_pre(b-1); ranking only");
    if (bad) {
      std::fprintf(stderr, "[beams] abort: a table with 2 or more indices has no active beam\n");
      return 2;
    }
    return 0;
  };

  std::string line, id;
  while (std::getline(std::cin, line)) {
    if (!line.empty() && line[0] == '>') {
      id = line.substr(1);
      continue;
    }
    std::string seq = normalize(line);
    if (seq.empty()) continue;
    if (restrict_pairs) {
      std::string pl;
      if (!std::getline(std::cin, pl)) {
        std::fprintf(stderr, "--restrict-pairs: missing pair line after the sequence\n");
        return 2;
      }
      restrict_list.clear();
      std::istringstream ps(pl);
      int a, b;
      while (ps >> a >> b) restrict_list.emplace_back(std::min(a, b) - 1, std::max(a, b) - 1);
      allowed_pairs = &restrict_list;
    }
    if (beam_audit(static_cast<int>(seq.size()))) return 2;

    auto t0 = std::chrono::steady_clock::now();
    const double cpu0 = ru_cpu();
    double pass1_s = 0, pass2_s = 0, scan_s = 0;
    double cpu_pass1 = 0, cpu_pass2 = 0;
    AccessRec access;
    if (need_access) access = load_access(access_idx, seq);
    std::unique_ptr<Gate> gate;
    if (gcfg.on || gate_diag)
      gate.reset(new Gate(static_cast<int>(seq.size()), access, gcfg.span, gcfg.delta, gcfg.both));
    apply_max_loop(loop_cap_scan);
    DP03 dp(seq, beam, allow_pk, loop_cap_scan, score_qz, arm_score, arm_score ? &access : nullptr,
            max_gap_scan, max_arm_len, gap_step, &gcfg, gcfg.on ? gate.get() : nullptr,
            table_beam, states_f != nullptr || dump_scores);
    dp.qz_beam_fill = !qzbar_old;
    dp.p2_qz_beam = b2_width;
    dp.qzbar_old = qzbar_old;
    if (qzbar_old)
      std::fprintf(stderr, "[qzbar-old] dense Qz→prQz (no beamed Qz_pre)\n");
    const std::string sid = id.empty() ? "-" : id;
    g_mem_sid = sid;
    if (g_mem_log) {
      std::vector<MemItem> it;
      dp.mem_items(it, "", true, false);
      mem_log("dp_alloc", it);
    }
    // Stem estimate, principled pre-pass and chain all use b1. Q^z uses b2.
    StemEst stem;
    std::vector<StemEst> diag_ests;
    double stem_sec = 0.0, chain_prep_sec = 0.0;
    // b1 is the nested beam at every length. The n<=1000 lean skip stays off once b1 is set.
    const bool lean = l2r_lean && pk_l2r == 1 && (l2r_beam < 0 ? beam : l2r_beam) > 0 &&
                      stem_rank && stem_pr && g_grow_rec && !stem_diag && !dump_scores &&
                      !dump_ext && !stem_hr_av && !stem_hr_phi && !stem_bg_check && nest_beam < 0;
    if (l2r_lean && !lean)
      std::fprintf(stderr, "[l2r-lean] %s not applicable to this run's flags: full pre-pass\n",
                   sid.c_str());
    stem.lean = lean;
    if (stem_rank || stem_diag || dump_scores || dump_ext) {
      auto ts = std::chrono::steady_clock::now();
      std::unique_ptr<DP03> nest;
      if (!lean) {
        const int scan_nb = nest_beam;
        apply_max_loop(loop_cap_p1);
        nest.reset(new DP03(seq, scan_nb, false, loop_cap_p1, score_qz, false, nullptr,
                            max_gap_scan, max_arm_len, true, nullptr, nullptr, nullptr, false));
        nest->chain_qmp = stem_chain;
        nest->qz_beam_fill = !qzbar_old;
        nest->p2_qz_beam = b2_width;
        nest->qzbar_old = qzbar_old;
        std::fprintf(stderr, "[nest-beam] b1 %d\n", scan_nb);
        if (pass1_only && !sample_prof.empty()) sample_start(sample_prof.c_str());
        const auto t_stem = std::chrono::steady_clock::now();
        g_beam_stats.set_pass("stem");
        nest->pfunc();
        pass1_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_stem).count();
        nest->log_nest("stem");
        if (g_mem_log) {
          std::vector<MemItem> it;
          nest->mem_items(it, "nest_", true, false);
          mem_log("nest_pfunc", it);
        }
      }
      if (stem_rank) {
        if (!lean) {
          const auto t_se = std::chrono::steady_clock::now();
          nest->build_stem_est(stem_loop, stem, stem_qg2);
          const double se_dn = static_cast<double>(seq.size());
          std::fprintf(stderr,
                       "[stem-est] sec %.6f entries %zu per_nt %.3f idx_mb %.1f t_mb %.1f hwm_mb %.1f\n",
                       std::chrono::duration<double>(std::chrono::steady_clock::now() - t_se).count(),
                       stem.entries, se_dn > 0 ? static_cast<double>(stem.entries) / se_dn : 0.0,
                       (12.0 * static_cast<double>(stem.entries + stem.entries2 + stem.entries_av) +
                        16.0 * static_cast<double>(stem.index.nodes() + stem.index2.nodes() +
                                                   stem.index_av.nodes())) /
                           1048576.0,
                       8.0 * se_dn * se_dn / 1048576.0, hwm_mb());
          if (g_mem_log) {
            std::vector<MemItem> it;
            stem_mem_items(stem, it);
            mem_log("stem_est", it);
          }
          if (stem_helpers == 2 || ml_tail) {
            nest->build_nested_beta(stem);
            if (g_mem_log) {
              std::vector<MemItem> it;
              nest->mem_items(it, "nest_", true, false);
              stem_mem_items(stem, it);
              mem_log("nested_beta", it);
            }
          }
        }
        dp.set_stem(&stem, stem_both, stem_qg2, stem_qg_partner, stem_consistent);
        dp.set_stem_ext(stem_ext, stem_ext_max, stem_split, stem_ext_check);
        dp.set_stem_grow(grow_lg, grow_cut, grow_nodes, stem_helpers, ml_span);
        if (stem_profile) dp.set_stem_profile(true, prof_eps, ml_tail, prof_check);
        if (bg_reorder) dp.set_stem_bg_reorder(true, ro_check);
        // exact nested-only run, shared by the chain and the principled pre-pass when b1 is 0
        std::unique_ptr<DP03> nexact;
        auto exact_nested = [&]() -> DP03& {
          if (!nexact) {
            apply_max_loop(loop_cap_p1);
            nexact.reset(new DP03(seq, 0, false, loop_cap_p1, score_qz, false, nullptr, max_gap_scan,
                                  max_arm_len, true, nullptr, nullptr, nullptr, false));
            nexact->chain_qmp = true;
            nexact->keep_closed_b = stem_pr_nested;
            nexact->prof_on = profile;
            nexact->qz_beam_fill = !qzbar_old;
            nexact->p2_qz_beam = b2_width;
            nexact->qzbar_old = qzbar_old;
            const auto t_ex = std::chrono::steady_clock::now();
            g_beam_stats.set_pass("pass1");
            nexact->pfunc();
            pass1_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_ex).count();
            nexact->log_nest("pass1");
            if (g_mem_log) {
              std::vector<MemItem> it;
              nexact->mem_items(it, "exact_", true, false);
              mem_log("exact_pfunc", it);
            }
          }
          return *nexact;
        };
        if (stem_pr || stem_pr_nested) {
          // Nested pre-pass at b1. --pre-beam replaces b1 when it is set.
          const auto tp = std::chrono::steady_clock::now();
          std::unique_ptr<DP03> nbeam;
          const int nb = pre_beam > 0 ? pre_beam : nest_beam;
          const bool use_beam = nb > 0;
          std::fprintf(stderr, "[nest-beam] pre %d\n", use_beam ? nb : 0);
          if (use_beam) {
            apply_max_loop(loop_cap_p1);
            nbeam.reset(new DP03(seq, nb, false, loop_cap_p1, score_qz, false, nullptr, max_gap_scan,
                                 max_arm_len, true, nullptr, nullptr, nullptr, false));
            nbeam->chain_qmp = true;
            nbeam->keep_closed_b = stem_pr_nested;
            nbeam->prof_on = profile;
            nbeam->qz_beam_fill = !qzbar_old;
            nbeam->p2_qz_beam = b2_width;
            nbeam->qzbar_old = qzbar_old;
            nbeam->pre_key_inside = pre_key == 1;
            nbeam->p2_sweep = p2_sweep;
            const auto t_pre = std::chrono::steady_clock::now();
            g_beam_stats.set_pass("pass1");
            nbeam->pfunc();
            pass1_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_pre).count();
            nbeam->log_nest("pass1");
            if (pass1_only) {
              if (!sample_prof.empty()) sample_stop();
              std::fprintf(stderr,
                           "[pass1-only] pass1_s %.6f stem_idx %zu pr_Mbar 0 pr_QbIn 0 grow_roots 0 "
                           "hwm_mb %.1f\n",
                           pass1_s, stem.entries, hwm_mb());
              nbeam->log_pr_side();
              return 0;
            }
          }
          DP03& np = nbeam ? *nbeam : exact_nested();
          apply_max_loop(loop_cap_p2);
          np.set_loop_cap(loop_cap_p2);
          auto stamp_p2 = [&](DP03& o) {
            o.qz_beam_fill = !qzbar_old;
            o.p2_qz_beam = b2_width;
            o.qzbar_old = qzbar_old;
            o.p2_qz_key = 1;
            o.p2_v_mode = p2_v_mode;
            o.p2_v_bp = p2_v_bp;
            o.p2_v_cap = p2_v_cap;
            o.p2_v_trim = p2_v_trim;
            o.p2_g_mode = p2_g_mode;
            o.p2_g_bp = p2_g_bp;
            o.p2_g_suffix = p2_g_suffix;
            o.p2_g_eager = p2_g_eager;
            o.p2_g_fb = p2_g_fb;
            o.p2_g_cap = p2_g_cap;
            o.sw_check = sw_check;
            o.p2_sweep = p2_sweep;
            o.p2_vsweep = p2_vsweep;
            o.p2_v_floor = p2_v_floor;
            o.p2_v_cut_fb = p2_v_cut_fb;
            o.lazy_ml_cap = lazy_ml_cap;
            o.prof_on = profile;
          };
          stamp_p2(np);
          stamp_p2(dp);
          stem.want_lp = oracle_f != nullptr || !pr_v_diag.empty() || !g_v2dump.empty() || g_v2ratio;
          stem.skip_kept = lean && oracle_f == nullptr && !pr_check;
          stem.skip_qbin = lean && l2r_p1 == 0 && !l2r_ctx_p1 && !pr_v_ratio &&
                           pr_v_diag.empty() && pr_v_sample <= 0.0 && g_v2dump.empty();
          cpu_pass1 = ru_cpu() - cpu0;
          const double hwm_pass1 = hwm_mb();
          const int nseq = static_cast<int>(seq.size());
          std::fprintf(stderr,
                       "[stage] pass1 n %d cpu %.3f cpu_s_nt %.6f hwm_mb %.1f hwm_mb_nt %.4f\n",
                       nseq, cpu_pass1, nseq > 0 ? cpu_pass1 / nseq : 0.0, hwm_pass1,
                       nseq > 0 ? hwm_pass1 / nseq : 0.0);
          const auto t_p2 = std::chrono::steady_clock::now();
          const double cpu_p2_0 = ru_cpu();
          g_beam_stats.set_pass("pass2");
          np.build_principled(stem, pr_floor, stem_pr_nested, pr_check, pr_partner_floor);
          if (g_sweep == G_SWEEP_STATES && nbeam) {
            dp.g_pushes = nbeam->g_pushes;
            dp.g_pre_n = nbeam->g_pre_n;
            dp.g_pre_sum = nbeam->g_pre_sum;
            dp.g_pre_max = nbeam->g_pre_max;
            dp.g_post_n = nbeam->g_post_n;
            dp.g_post_sum = nbeam->g_post_sum;
            dp.g_post_max = nbeam->g_post_max;
            dp.g_merge_sec = nbeam->g_merge_sec;
            dp.g_cut_sec = nbeam->g_cut_sec;
          }
          pass2_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_p2).count();
          cpu_pass2 = ru_cpu() - cpu_p2_0;
          const double hwm_pass2 = hwm_mb();
          std::fprintf(stderr,
                       "[stage] pass2 n %d cpu %.3f cpu_s_nt %.6f hwm_mb %.1f hwm_mb_nt %.4f\n",
                       nseq, cpu_pass2, nseq > 0 ? cpu_pass2 / nseq : 0.0, hwm_pass2,
                       nseq > 0 ? hwm_pass2 / nseq : 0.0);
          stem.prof_qz_pre = np.prof_qz_sec;
          if (g_mem_log) {
            std::vector<MemItem> it;
            np.mem_items(it, nbeam ? "beam_" : "exact_", true, false);
            stem_mem_items(stem, it);
            mem_log("principled", it);
          }
          if (stem_pr) {
            dp.set_stem_principled(true, chain_loop, chain_ml, pr_ml_span, pr_ml_close, pr_cut,
                                   pr_arms_k, pr_check);
            dp.set_stem_pr_select(pr_ratio, pr_arms_eps, pr_arms_max, pr_gapfree, pr_cut_roots);
            dp.set_stem_pr_sep(pr_arms_mode, pr_arm_z, pr_qz_unp, pr_sep_check);
            dp.set_stem_pr_grow(pr_gsel, pr_geps, pr_gB, pr_gcheck, pr_dump_f);
            dp.set_stem_pr_armq_sep(pr_aq_sep, pr_aq_max);
            dp.set_stem_pr_arms(pr_arm_persist, pr_arm_gaponly, pr_pspan, pr_veps);
            dp.pr_vsuffix = pr_vsuffix;
            dp.pr_vcap = pr_vcap;
            dp.pr_vmlk = pr_vmlk;
            dp.pr_vmlsteps = pr_vmlsteps;
            dp.pr_vxs = pr_vxs;
            dp.pr_vstream = pr_vstream;
            dp.set_stem_pr_partner_cap(pr_pmax, pr_pcheck);
            dp.set_stem_pr_partner_ml(pr_pml, pr_pmlcheck);
            dp.set_stem_pr_v(pr_v_ratio, pr_v_diag, pr_v_sample);
            dp.pr_fast = pr_fast;
            dp.pr_abl = pr_abl;
            dp.pr_abl_pk = pr_abl_pk;
            if (g_grow_rec) {
              dp.rnL = g_ml_diag ? 3 : g_ml_L + 1;
              dp.rnR = g_ml_diag ? 3 : g_ml_R + 1;
              dp.rec_noml = dp.rnL == 1 && dp.rnR == 1;
            }
          }
          if (stem_pr_nested) dp.set_stem_pr_nested(true);
          stem.pr_prep_sec =
              std::chrono::duration<double>(std::chrono::steady_clock::now() - tp).count();
        }
        if ((stem_chain && !lean) || stem_hr_av || stem_hr_phi) {
          const auto tc = std::chrono::steady_clock::now();
          // Chain and hairpin tables use the same b1 nested DP.
          if (!lean && nest) {
            if (stem_chain) nest->build_chain(stem, chain_floor);
            if (stem_hr_av) nest->build_hr_av(stem, stem_hr_gamma);
            if (stem_hr_phi) nest->build_hr_phi(stem, stem_hr_gamma);
          }
          if (stem_chain) dp.set_stem_chain(true, chain_span, chain_cut, chain_ml, chain_loop);
          if (stem_hr_av || stem_hr_phi) dp.set_stem_hr(true, stem_hr_av, stem_hr_check);
          chain_prep_sec =
              std::chrono::duration<double>(std::chrono::steady_clock::now() - tc).count();
          if (g_mem_log) {
            std::vector<MemItem> it;
            stem_mem_items(stem, it);
            mem_log("chain", it);
          }
        } else if (stem_chain) {
          dp.set_stem_chain(true, chain_span, chain_cut, chain_ml, chain_loop);
        }
      }
      if (dump_scores || dump_ext) nest->build_stem_est(stem_loop, stem, true);
      diag_ests.resize(stem_diag ? stem_diag_loops.size() : 0);
      for (size_t k = 0; k < diag_ests.size(); ++k)
        nest->build_stem_est(stem_diag_loops[k], diag_ests[k], true);
      stem_sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - ts).count();
    }
    if (g_mem_log) {
      std::vector<MemItem> it;
      stem_mem_items(stem, it);
      mem_log("stem_done", it);
    }
    DropLog drops;
    if (pruned_f) {
      drops.f = pruned_f;
      drops.id = sid;
      drops.n = std::max(1, static_cast<int>(seq.size()));
      drops.min_w = pruned_min_w;
      drops.kf = kept_f;
      if (!exact_path.empty()) {
        auto it = exact_states.find(sid);
        drops.exact = it == exact_states.end() ? &no_states : &it->second;
        if (it == exact_states.end())
          std::fprintf(stderr, "--dump-exact: no states for %s\n", sid.c_str());
      }
      dp.drops = &drops;
    }
    if (!remove_keys.empty()) dp.remove_keys = &remove_keys;
    if (pk_l2r != 0) {
      if (pk_l2r == 2 && beam != 0) {
        std::fprintf(stderr, "--pk-l2r check needs --beam inf\n");
        return 2;
      }
      dp.l2r_mode = pk_l2r;
      dp.l2r_check_R = pk_l2r_R;
      dp.l2r_beam = pk_l2r == 1 ? (l2r_beam < 0 ? beam : l2r_beam) : 0;
      dp.l2r_spawn_k = dp.l2r_beam > 0 ? l2r_spawn_k : 0;
      dp.sp_cube = l2r_spawn_cube && dp.l2r_beam > 0;
      dp.l2r_arm_fb = l2r_arm_fb && dp.l2r_beam > 0;
      if (dp.sp_cube && l2r_spawn_k <= 0) {
        std::fprintf(stderr, "--l2r-spawn cube needs --l2r-spawn-k k > 0\n");
        return 2;
      }
      dp.sp_pool = l2r_spawn_pool;
      dp.sp_wk = l2r_spawn_wk;
      dp.sp_static = l2r_spawn_static;
      dp.sp_maxpop = l2r_spawn_maxpop;
      if (l2r_spawn_inc && pk_l2r != 1) {
        std::fprintf(stderr, "--l2r-spawn inc needs --pk-l2r on\n");
        return 2;
      }
      if (l2r_inc_check && (!l2r_spawn_inc || dp.l2r_beam > 0)) {
        std::fprintf(stderr, "--l2r-inc-check needs --l2r-spawn inc and no l2r beam\n");
        return 2;
      }
      if (l2r_spawn_inc && g_ntr_K > 0 && !g_ntr_zbar) {
        if (ntr_z_set) {
          std::fprintf(stderr, "--l2r-spawn inc ranks growth nodes from step a2 on: needs --grow-node-trunc-z bar\n");
          return 2;
        }
        g_ntr_zbar = true;
      }
      dp.l2r_inc = l2r_spawn_inc;
      dp.l2r_ninf_mass = l2r_ninf_mass;
      dp.l2r_keep_ninf = l2r_keep_ninf;
      dp.l2r_inc_check = l2r_inc_check;
      dp.qz_audit_req = l2r_qz_audit;
      dp.lp_hk = l2r_p_hk;
      dp.lp_aopen_raw = l2r_p_raw;
      dp.lp_hex = l2r_p3_hex;
      dp.lp_lin = l2r_p3_lin;
      for (int k = 0; k < 3; ++k) dp.l2r_tb[DP03::PT_P1 + k] = dp.l2r_beam > 0 ? l2r_pbeam[k] : 0;
      if (l2r_p3_split[0] > 0 && dp.l2r_beam > 0) {
        dp.lp_split_a = l2r_p3_split[0];
        dp.lp_split_h = l2r_p3_split[1];
        dp.l2r_tb[DP03::PT_P3] = l2r_p3_split[0] * l2r_p3_split[1];
      }
      {
        static const char* pn[DP03::PT_N] = {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us",
                                            "Uls", "Ulss", "Urs", "S", "P1", "P2", "P3"};
        std::stringstream ss(l2r_tbeam);
        std::string item;
        while (std::getline(ss, item, ',')) {
          const size_t eq = item.find('=');
          int t = 0;
          while (t < DP03::PT_N && (eq == std::string::npos || item.substr(0, eq) != pn[t])) ++t;
          if (t == DP03::PT_N) {
            std::fprintf(stderr, "--l2r-tbeam: bad item %s\n", item.c_str());
            return 2;
          }
          const std::string v = item.substr(eq + 1);
          dp.l2r_tb[t] = v == "inf" ? 0 : std::atoi(v.c_str());
        }
      }
      if (dp.l2r_beam > 0 && !(stem_pr && g_grow_rec)) {
        std::fprintf(stderr, "--pk-l2r on with a beam needs --stem-est-principled on and --grow-rec on\n");
        return 2;
      }
      dp.l2r_diag = l2r_diag && dp.l2r_beam > 0;
      dp.rlt = rlt;
      dp.next_pair = next_pair;
      dp.rlt_beam = rlt_beam;
      dp.grammar_simple = grammar_simple;
      dp.c_trace = c_trace;
      if (grammar_simple && !l2r_spawn_inc) {
        std::fprintf(stderr, "--grammar simple needs --l2r-spawn inc\n");
        return 2;
      }
      dp.l2r_init();
      if (g_mem_log) {
        std::vector<MemItem> it;
        dp.mem_items(it, "", true, true);
        stem_mem_items(stem, it);
        mem_log("l2r_init", it);
      }
      // --l2r-p1 R: R rounds of a beamed pseudoknot pass (a copy of this run, beams capped at
      // --l2r-p1-beam) scored with the current Sbar; its stored Qp items with the pre-pass
      // nested Qb (--l2r-p1-qb nested) or its own Qb (pass1) give Sbar_pk, which replaces Sbar
      // for the next round and the final pass.
      // --l2r-ctx p1: every pass-1 round also gives lambda(a') from its Qp items' inside and
      // outside (applied to the next round's scores and the final pass).
      if (l2r_ctx_p1 && dp.l2r_beam > 0) l2r_p1 = std::max(l2r_p1, l2r_ctx_rounds);
      if (l2r_p1 > 0 && dp.l2r_beam > 0) {
        if (!dp.l2r_stem_is(&stem)) {
          std::fprintf(stderr, "--l2r-p1: the run's stem estimator is not main's\n");
          return 2;
        }
        const std::vector<double> S0 = stem.prS;
        for (int r = 1; r <= l2r_p1; ++r) {
          const double c0 = cpu_sec();
          std::vector<double> sn, sp;
          double z1 = NEG_INF, c1 = c0;
          long long nqp = 0, nqb = 0;
          std::string prof;
          std::vector<double> lam;
          std::vector<int> lcnt;
          double c3 = 0.0;
          {
            DP03 p1(dp);
            p1.l2r_cap_beam(l2r_p1_beam);
            p1.l2r_diag = false;
            p1.drops = nullptr;
            p1.mem_scan = false;
            p1.l2r_yw_rebuild();
            p1.l2r_init_scores();
            z1 = p1.pfunc();
            c1 = cpu_sec();
            prof = p1.l2r_prof(c1 - c0);
            dp.l2r_sfx_pass(p1, sn, sp, nqp, nqb);
            if (l2r_ctx_p1) {
              const double c4 = cpu_sec();
              dp.l2r_ctx_lambda(p1, l2r_p1_own ? sp : sn, l2r_ctx_tau, lam, lcnt);
              c3 = cpu_sec() - c4;
            }
            if (g_mem_log) {
              std::vector<MemItem> it;
              p1.mem_items(it, "p1_", true, true);
              mem_log(("p1_round" + std::to_string(r)).c_str(), it);
            }
          }
          if (l2r_ctx_p1) {
            std::vector<double> lv;
            int hi = 0, lo = 0;
            std::string pos;
            char b2[64];
            for (size_t a = 0; a < lam.size(); ++a) {
              if (lcnt[a] == 0) continue;
              lv.push_back(lam[a]);
              if (lam[a] > std::log(2.0)) ++hi;
              if (lam[a] < std::log(0.5)) ++lo;
              std::snprintf(b2, sizeof b2, " %zu:%.3f", a, lam[a]);
              pos += b2;
            }
            std::sort(lv.begin(), lv.end());
            char b[400];
            std::snprintf(b, sizeof b,
                          "[l2r-ctx] %s round %d tau %.3g a_items %zu lambda_gt2 %d lambda_lt0.5 %d "
                          "loglam_min %.4f loglam_med %.4f loglam_max %.4f ctx_cpu %.3f | a:loglam",
                          sid.c_str(), r, l2r_ctx_tau, lv.size(), hi, lo, lv.empty() ? 0.0 : lv.front(),
                          lv.empty() ? 0.0 : lv[lv.size() / 2], lv.empty() ? 0.0 : lv.back(), c3);
            std::fprintf(stderr, "%s%s\n", b, pos.c_str());
            dp.lLam = std::move(lam);
          }
          const double c2 = cpu_sec();
          std::vector<double>& spk = l2r_p1_own ? sp : sn;
          std::string ln = "[l2r-p1] " + sid + " round " + std::to_string(r) + " beam " +
                           std::to_string(l2r_p1_beam) + " qb " + (l2r_p1_own ? "pass1" : "nested");
          char b[400];
          std::snprintf(b, sizeof b,
                        " logZ_pass1 %.9f logS_pk1_nestedqb %.9f logS_pk1_pass1qb %.9f "
                        "logS_nested1 %.9f qp_items %lld qb_items %lld p1_cpu %.3f sfx_cpu %.3f "
                        "hwm_mb %.0f | %s | y:log(Sbar_pk[y+1]/Sbar_nested[y+1])",
                        z1, sn[0], sp[0], S0[0], nqp, nqb, c1 - c0, c2 - c1, hwm_mb(), prof.c_str());
          ln += b;
          for (size_t y = 0; y + 1 < spk.size(); ++y) {
            std::snprintf(b, sizeof b, " %zu:%.2f", y, spk[y + 1] - S0[y + 1]);
            ln += b;
          }
          std::fprintf(stderr, "%s\n", ln.c_str());
          stem.prS = std::move(spk);
        }
        dp.l2r_yw_rebuild();
      }
      dp.lpf = l2r_perf;
      dp.k_memo = k_memo || k_exact_ub;
      dp.k_bound = k_bound;
      dp.k_exact_ub = k_exact_ub;
      dp.k_exact_check = k_exact_check;
      if (k_exact) dp.sx_k_exact = true;
      dp.rlt = rlt;
      dp.next_pair = next_pair;
      dp.rlt_beam = rlt_beam;
      dp.grammar_simple = grammar_simple;
      dp.c_trace = c_trace;
      dp.l2r_init_scores();
      if (p2_g_eager) dp.growth_eager();
      if (g_mem_log) {
        std::vector<MemItem> it;
        dp.mem_items(it, "", true, true);
        stem_mem_items(stem, it);
        mem_log("init_scores", it);
        dp.mem_scan = true;
      }
      dp.l2r_split = l2r_split && dp.l2r_diag;
      dp.lsp_rec = dp.l2r_split;
      dp.l2r_ctx = l2r_ctx && dp.l2r_beam > 0;
      dp.l2r_ctx_init();
      dp.lsub_a = l2r_kx || l2r_sub.find('a') != std::string::npos;
      dp.lsub_c = l2r_sub.find('c') != std::string::npos;
      dp.lsub_d = l2r_sub.find('d') != std::string::npos;
      dp.lsub_e = l2r_sub.find('e') != std::string::npos;
      if ((dp.lsub_c || dp.lsub_d || dp.lsub_e) && !dp.l2r_oracle_load(l2r_oload)) {
        std::fprintf(stderr, "--l2r-sub c/d/e needs a readable --l2r-oracle for this sequence\n");
        return 1;
      }
      if (!l2r_trace.empty()) {
        if (std::sscanf(l2r_trace_key.c_str(), "%d,%d,%d,%d,%d", &dp.ltr_ap, &dp.ltr_xs, &dp.ltr_xp, &dp.ltr_a2,
                        &dp.ltr_c1) < 3 ||
            !(dp.ltr_f = std::fopen(l2r_trace.c_str(), "w"))) {
          std::fprintf(stderr, "--l2r-trace needs a writable file and --l2r-trace-key a',xs,x'\n");
          return 1;
        }
      }
    }
    double logz = NEG_INF;
    bool partial = false;
    const double cpu_p2 = cpu_sec(), hwm_p2 = hwm_mb();
    if (!sample_prof.empty()) sample_start(sample_prof.c_str());
    const auto t_scan = std::chrono::steady_clock::now();
    g_beam_stats.set_pass("scan");
    apply_max_loop(loop_cap_scan);
    dp.set_loop_cap(loop_cap_scan);
    try {
      logz = dp.pfunc();
    } catch (const CpuStop&) {
      partial = true;
    }
    scan_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_scan).count();
    {
      const double cpu_scan = ru_cpu() - cpu0 - cpu_pass1 - cpu_pass2;
      const double hwm_scan = hwm_mb();
      const int nseq = static_cast<int>(seq.size());
      const double ntd = nseq > 0 ? static_cast<double>(nseq) : 1.0;
      const double idx_mb =
          (12.0 * static_cast<double>(stem.entries + stem.entries2 + stem.entries_av) +
           16.0 * static_cast<double>(stem.index.nodes() + stem.index2.nodes() +
                                      stem.index_av.nodes())) /
          1048576.0;
      std::fprintf(stderr,
                   "[stage] scan n %d cpu %.3f cpu_s_nt %.6f hwm_mb %.1f hwm_mb_nt %.4f\n",
                   nseq, cpu_scan, cpu_scan / ntd, hwm_scan, hwm_scan / ntd);
      std::fprintf(stderr,
                   "[tables] n %d stem_idx %zu stem_idx_mb %.1f pr_Mbar %zu pr_Mbar_mb %.1f "
                   "pr_QbIn %zu pr_QbIn_mb %.1f n2_mb %.1f\n",
                   nseq, stem.entries, idx_mb, stem.prMbar.size(),
                   static_cast<double>(stem.prMbar.size()) * 8.0 / 1048576.0, stem.prQbIn.size(),
                   static_cast<double>(stem.prQbIn.size()) * 8.0 / 1048576.0,
                   8.0 * ntd * ntd / 1048576.0);
    }
    if (!sample_prof.empty()) sample_stop();
    if (g_mem_log) {
      std::vector<MemItem> it;
      dp.mem_items(it, "", true, pk_l2r != 0);
      mem_log("scan_done", it);
    }
    if (pk_l2r == 1 && dp.l2r_beam > 0) {
      const double fin = cpu_sec() - cpu_p2;
      std::fprintf(stderr, "[l2r-time] %s before_final_cpu %.3f final_cpu %.3f hwm_before_mb %.0f "
                   "hwm_mb %.0f | final %s\n", sid.c_str(), cpu_p2, fin, hwm_p2, hwm_mb(),
                   dp.l2r_prof(fin + dp.l2r_isec).c_str());
      dp.profile_dump(sid.c_str(), fin, stem.prof_qz_pre);
    }
    if (g_mass_f) g_mass_flush();
    if (b2_width > 0 || p2_qz_beam > 0 || p2_v_mode > 0 || p2_g_mode > 0 || p2_v_cap > 0 ||
        p2_g_suffix || p2_g_eager || p2_g_fb || p2_sweep || p2_vsweep) {
      std::fprintf(stderr,
                   "[pass2] %s qz %lld qz_mb %.3f qz_s %.4f v_pairs %lld v_ent %lld v_s %.4f "
                   "g_L %lld prep_s %.4f\n",
                   sid.c_str(), stem.p2_qz_n,
                   static_cast<double>(stem.p2_qz_n * static_cast<long long>(sizeof(StemEst::QzEnt))) /
                       1048576.0,
                   stem.p2_qz_sec, dp.v_pairs, dp.pr_V_entries, dp.pr_static_sec, dp.ntr_kept[0],
                   stem.pr_prep_sec);
    }
    if (pk_l2r != 0) {
      const std::vector<long long> sc = dp.l2r_states();
      std::fprintf(stderr,
                   "[l2r] %s n %d mode %s states R %lld Rs %lld Rls %lld Rlss %lld Rrs %lld U %lld "
                   "Us %lld Uls %lld Ulss %lld Urs %lld S %lld tr R %lld Rs %lld Rls %lld Rrs %lld "
                   "U %lld Us %lld Uls %lld Urs %lld S %lld check_s %.3f",
                   sid.c_str(), static_cast<int>(seq.size()), pk_l2r == 1 ? "on" : "check", sc[0],
                   sc[1], sc[2], sc[3], sc[4], sc[5], sc[6], sc[7], sc[8], sc[9], sc[10],
                   dp.l2r_tr[DP03::LT_AG], dp.l2r_tr[DP03::LT_AGS], dp.l2r_tr[DP03::LT_ALS],
                   dp.l2r_tr[DP03::LT_ARS], dp.l2r_tr[DP03::LT_BG], dp.l2r_tr[DP03::LT_BGS],
                   dp.l2r_tr[DP03::LT_BLS], dp.l2r_tr[DP03::LT_BRS], dp.l2r_tr[DP03::LT_S],
                   dp.l2r_sec[5]);
      if (pk_l2r == 2)
        std::fprintf(stderr, " checkR %lld mis %lld max_rel %.3e checkQp %lld mis %lld max_rel %.3e",
                     dp.lc_R_n, dp.lc_R_mis, dp.lc_R_rel, dp.lc_Qp_n, dp.lc_Qp_mis, dp.lc_Qp_rel);
      std::fprintf(stderr, "\n");
      const double nn = static_cast<double>(seq.size());
      auto per = [&](long long x) { return nn > 0.0 ? static_cast<double>(x) / nn : 0.0; };
      const long long r_ml = dp.ext_n[1] + dp.ext_n[2] + dp.ext_n[3];
      const long long u_ml = dp.ext_n[5] + dp.ext_n[6] + dp.ext_n[7];
      std::fprintf(stderr,
                   "[ext] n %d R_int %.3f R_ml %.3f R_f %.3f R_g %.3f R_h %.3f "
                   "U_int %.3f U_ml %.3f U_f %.3f U_g %.3f U_h %.3f\n",
                   static_cast<int>(seq.size()), per(dp.ext_n[0]), per(r_ml), per(dp.ext_n[1]),
                   per(dp.ext_n[2]), per(dp.ext_n[3]), per(dp.ext_n[4]), per(u_ml), per(dp.ext_n[5]),
                   per(dp.ext_n[6]), per(dp.ext_n[7]));
      std::fprintf(stderr,
                   "[it] n %d R_int %.3f R_h %.3f R_ls %.3f R_rs %.3f "
                   "U_int %.3f U_h %.3f U_ls %.3f U_rs %.3f\n",
                   static_cast<int>(seq.size()), per(dp.it_int[0]), per(dp.it_h[0]),
                   per(dp.it_ls[0]), per(dp.it_rs[0]), per(dp.it_int[1]), per(dp.it_h[1]),
                   per(dp.it_ls[1]), per(dp.it_rs[1]));
      std::fprintf(stderr,
                   "[qp] n %d beam %d pre_per %.3f kept_per %.3f pre_max %d kept_max %d bind %lld\n",
                   static_cast<int>(seq.size()), dp.qp_beam, per(dp.qp_pre_n),
                   per(dp.qp_kept_n), dp.qp_pre_max, dp.qp_kept_max, dp.qp_bind_n);
      std::fprintf(stderr,
                   "[qm] n %d beam %d pre_per %.3f kept_per %.3f pre_max %d kept_max %d bind %lld\n",
                   static_cast<int>(seq.size()), dp.qm_beam, per(dp.qm_pre_n), per(dp.qm_kept_n),
                   dp.qm_pre_max, dp.qm_kept_max, dp.qm_bind_n);
      long long scored = 0;
      for (int t = DP03::PT_R; t < DP03::PT_N; ++t) scored += dp.l2r_pre[t];
      const double wall_s =
          std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      const double ntd = nn > 0.0 ? nn : 1.0;
      std::fprintf(stderr,
                   "[pos] n %d per_nt interior %.3f h %.3f ls %.3f rs %.3f tail %.3f "
                   "S %.3f P1 %.3f P2 %.3f P3 %.3f Qz %.3f Qz_pre %.3f K %.3f scored %.3f\n",
                   static_cast<int>(seq.size()),
                   per(dp.it_int[0] + dp.it_int[1]), per(dp.it_h[0] + dp.it_h[1]),
                   per(dp.it_ls[0] + dp.it_ls[1]), per(dp.it_rs[0] + dp.it_rs[1]), per(dp.rlt_ext),
                   per(dp.l2r_tr[DP03::LT_S]), per(dp.l2r_pre[DP03::PT_P1]),
                   per(dp.l2r_pre[DP03::PT_P2]), per(dp.l2r_pre[DP03::PT_P3]), per(dp.qz_cand),
                   per(stem.p2_qz_n), per(dp.sx_k_n), per(scored));
      std::fprintf(stderr,
                   "[vmiss] n %d floor %d cut_fb %d look %lld miss %lld look_U %lld miss_U %lld "
                   "look_R %lld miss_R %lld look_S %lld miss_S %lld key0_U %lld part_U %lld "
                   "key0_R %lld fb %lld\n",
                   static_cast<int>(seq.size()), dp.p2_v_floor, dp.p2_v_cut_fb ? 1 : 0, dp.rfv_look,
                   dp.rfv_miss, dp.rfv_look_u, dp.rfv_miss_u, dp.rfv_look_r, dp.rfv_miss_r,
                   dp.rfv_look_s, dp.rfv_miss_s, dp.rfv_key0_u, dp.rfv_part_u, dp.rfv_key0_r,
                   dp.vfb_n);
      {
        static const char* sn[DP03::ST_N] = {"P1", "P2", "P3", "R", "S", "U"};
        for (int s = 0; s < DP03::ST_N; ++s)
          std::fprintf(stderr,
                       "[vmiss3] %s V %lld %lld %lld K %lld %lld %lld H %lld %lld %lld\n",
                       sn[s], dp.vm_i[s], dp.vm_ii[s], dp.vm_iii[s], dp.km_i[s], dp.km_ii[s],
                       dp.km_iii[s], dp.hm_i[s], dp.hm_ii[s], dp.hm_iii[s]);
      }
      if (g_miss_cap > 0 && !g_miss_path.empty()) {
        std::ofstream mf(g_miss_path);
        mf << "seen V " << g_miss_seen[0] << " K " << g_miss_seen[1] << " H " << g_miss_seen[2]
           << "\n";
        static const char kn[3] = {'V', 'K', 'H'};
        for (int k = 0; k < 3; ++k)
          for (const MissRec& r : g_miss[k])
            mf << kn[k] << " " << r.st << " " << r.x << " " << r.y << " " << r.b << "\n";
      }
      std::fprintf(stderr,
                   "[times] n %d pass1 %.3f pass2 %.3f scan %.3f wall %.3f "
                   "pass1_s_nt %.6f pass2_s_nt %.6f scan_s_nt %.6f wall_s_nt %.6f\n",
                   static_cast<int>(seq.size()), pass1_s, pass2_s, scan_s - dp.rfv_sec, wall_s,
                   pass1_s / ntd, pass2_s / ntd, (scan_s - dp.rfv_sec) / ntd, wall_s / ntd);
      std::fprintf(stderr,
                   "[vfb] sec %.4f calls %lld  Qz cand_per %.3f kept_per %.3f reads_per %.3f\n",
                   dp.rfv_sec, dp.rfv_od, per(dp.qz_cand), per(dp.qz_kept), per(dp.qz_reads));
      if (!dp.col_cnt.empty()) {
        const int nc = static_cast<int>(dp.col_cnt.size());
        const int a = (nc * 4) / 5;
        auto tail = [&](auto get) {
          long long prev = a > 0 ? get(dp.col_cnt[static_cast<size_t>(a - 1)]) : 0;
          long long s = 0;
          for (int j = a; j < nc; ++j) {
            const long long v = get(dp.col_cnt[static_cast<size_t>(j)]);
            s += v - prev;
            prev = v;
          }
          return nc > a ? static_cast<double>(s) / (nc - a) : 0.0;
        };
        std::fprintf(stderr,
                     "[pos20] n %d from %d int %.3f h %.3f ls %.3f rs %.3f tail %.3f "
                     "edge %.3f qz_reads %.3f K %.3f qm_kept %.3f\n",
                     nc, a, tail([](const DP03::ColCnt& c) { return c.iti; }),
                     tail([](const DP03::ColCnt& c) { return c.ith; }),
                     tail([](const DP03::ColCnt& c) { return c.itls; }),
                     tail([](const DP03::ColCnt& c) { return c.itrs; }),
                     tail([](const DP03::ColCnt& c) { return c.tail; }),
                     tail([](const DP03::ColCnt& c) { return c.edge; }),
                     tail([](const DP03::ColCnt& c) { return c.qz; }),
                     tail([](const DP03::ColCnt& c) { return c.k; }),
                     tail([](const DP03::ColCnt& c) { return c.qm; }));
      }
      {
        std::vector<MemItem> mit;
        dp.mem_items(mit, "", true, pk_l2r != 0);
        stem_mem_items(stem, mit);
        double tot = 0;
        for (const MemItem& m : mit) tot += m.bytes;
        std::fprintf(stderr, "[mem-nt] n %d total_mb %.3f total_mb_nt %.6f",
                     static_cast<int>(seq.size()), tot / 1048576.0, tot / 1048576.0 / ntd);
        for (const MemItem& m : mit) {
          if (m.bytes < 1024.0 && m.ents == 0) continue;
          std::fprintf(stderr, " %s %lld %.3f %.6f", m.name.c_str(), m.ents, m.bytes / 1048576.0,
                       m.bytes / 1048576.0 / ntd);
        }
        std::fprintf(stderr, "\n");
      }
      dp.log_nest("scan");
      if (dp.rlt) {
        std::fprintf(stderr, "[rlt] n %d beam %d kept_per %.3f max %d\n", static_cast<int>(seq.size()),
                     dp.rlt_beam, per(dp.rlt_kept), dp.rlt_max);
        std::fprintf(stderr, "[rlt-fam] R %lld max %d U %lld max %d\n", dp.rlt_kept_f[0], dp.rlt_max_f[0],
                     dp.rlt_kept_f[1], dp.rlt_max_f[1]);
        const int rn = static_cast<int>(dp.rlt_pos[0].size());
        const int un = static_cast<int>(dp.rlt_pos[1].size());
        const int pn = std::max(rn, un);
        for (int j = 0; j < pn; ++j) {
          const int rv = j < rn ? dp.rlt_pos[0][j] : 0;
          const int uv = j < un ? dp.rlt_pos[1][j] : 0;
          if (rv || uv)
            std::fprintf(stderr, "[rlt-pos] %d %d %d\n", j, rv, uv);
        }
      }
      if (dp.l2r_inc) {
        long long st[3] = {};
        for (int k = 0; k < 3; ++k)
          for (const auto& m : k == 0 ? dp.LP1 : k == 1 ? dp.LP2 : dp.LP3)
            st[k] += static_cast<long long>(m.keys.size());
        std::string b;
        char buf[256];
        for (int k = 0; k < 3; ++k) {
          const int t = DP03::PT_P1 + k;
          std::snprintf(buf, sizeof buf,
                        " P%d beam %d pre %lld kept %lld max %lld stored %lld tr %lld build_s %.3f "
                        "score_s %.3f miss_arm %lld miss_H %lld miss_K %lld miss_open %lld",
                        k + 1, dp.ltb(t), dp.l2r_pre[t], dp.l2r_kept[t], dp.l2r_kmax[t], st[k],
                        dp.l2r_tr[DP03::LT_P1 + k], dp.l2r_bsec[t], dp.l2r_psec[t], dp.lp_miss[k][0],
                        dp.lp_miss[k][1], dp.lp_miss[k][2], dp.lp_miss[k][3]);
          b += buf;
        }
        std::fprintf(stderr,
                     "[l2r-inc] %s n %d hk %d aopen %s%s seed_checks %lld seeds %lld Rs_pre %lld "
                     "Rs_kept %lld seed_s %.3f H_lists %lld H_roots %lld H_ents %lld H_s %.3f "
                     "check %lld mis %lld max_rel %.3e p3_split %d,%d a2_groups %lld a2_kept %lld "
                     "p3_hex %d hex_builds %lld hex_terms %lld hex_ents %lld hex_s %.3f p3_lin %d\n",
                     sid.c_str(), static_cast<int>(seq.size()), dp.lp_hk,
                     dp.lp_aopen_raw ? "raw" : "q", b.c_str(), dp.lp_seed_checks, dp.lp_seeds,
                     dp.l2r_pre[DP03::PT_RS], dp.l2r_kept[DP03::PT_RS], dp.l2r_bsec[DP03::PT_RS],
                     dp.lp_h_lists, dp.lp_h_roots, dp.lp_h_ents, dp.lp_h_sec, dp.lc_I_n, dp.lc_I_mis,
                     dp.lc_I_rel, dp.lp_split_a, dp.lp_split_h, dp.lp_split_grp, dp.lp_split_grp_kept,
                     dp.lp_hex ? 1 : 0, dp.lp_hx_builds, dp.lp_hx_terms, dp.lp_hx_ents, dp.lp_hx_sec,
                     dp.lp_lin_ok() ? 1 : 0);
      }
      if (dp.lpf) {
        std::string b;
        char buf[160];
        auto I = [&](const char* k, long long v) {
          std::snprintf(buf, sizeof buf, " %s %lld", k, v);
          b += buf;
        };
        auto D = [&](const char* k, double v) {
          std::snprintf(buf, sizeof buf, " %s %.4f", k, v);
          b += buf;
        };
        static const char* pn[DP03::PT_N] = {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us",
                                            "Uls", "Ulss", "Urs", "S", "P1", "P2", "P3"};
        std::string k;
        for (int t = 0; t < DP03::PT_N; ++t) {
          if (dp.lpf_npos[t] == 0) continue;
          k = std::string("a_") + pn[t] + "_";
          I((k + "beam").c_str(), dp.ltb(t));
          I((k + "pos").c_str(), dp.lpf_npos[t]);
          I((k + "pre").c_str(), dp.l2r_pre[t]);
          I((k + "premax").c_str(), dp.lpf_premax[t]);
          I((k + "kept").c_str(), dp.l2r_kept[t]);
          I((k + "keptmax").c_str(), dp.l2r_kmax[t]);
          I((k + "full").c_str(), dp.lpf_full[t]);
          I((k + "scored").c_str(), dp.lpf_nsc[t]);
          D((k + "build_s").c_str(), dp.l2r_bsec[t]);
          D((k + "prune_s").c_str(), dp.lpf_prune_sec[t]);
        }
        static const char* fn[3] = {"S", "P1", "P3"};
        for (int f = 0; f < 3; ++f) {
          if (dp.lpf_fx_pos[f] == 0) continue;
          k = std::string("c_") + fn[f] + "_";
          I((k + "pos").c_str(), dp.lpf_fx_pos[f]);
          I((k + "items").c_str(), dp.lpf_fx_items[f]);
          I((k + "pair").c_str(), dp.lpf_fx_pair[f]);
          I((k + "pairmax").c_str(), dp.lpf_fx_pairmax[f]);
          I((k + "unp").c_str(), dp.lpf_fx_unp[f]);
          D((k + "s").c_str(), dp.lpf_fx_sec[f]);
        }
        for (int t = 0; t < 3; ++t) {
          k = std::string("d_P") + std::to_string(t + 1) + "_";
          I((k + "calls").c_str(), dp.lpf_sc_calls[t]);
          I((k + "hsteps").c_str(), dp.lpf_sc_hsteps[t]);
          I((k + "dot").c_str(), dp.lpf_sc_dot[t]);
          I((k + "kcalls").c_str(), dp.lpf_sc_kcalls[t]);
          I((k + "ksteps").c_str(), dp.lpf_sc_ksteps[t]);
          I((k + "osteps").c_str(), dp.lpf_sc_osteps[t]);
          D((k + "h_s").c_str(), dp.lpf_sc_hsec[t]);
          D((k + "k_s").c_str(), dp.lpf_sc_ksec[t]);
          D((k + "o_s").c_str(), dp.lpf_sc_osec[t]);
        }
        I("d_pspan", dp.pr_pspan);
        I("d_aq_sep", dp.pr_aq_sep ? 1 : 0);
        I("e_lists", dp.lp_h_lists);
        I("e_roots", dp.lp_h_roots);
        I("e_newroots", dp.lpf_h_newroots);
        I("e_lam", dp.lpf_h_lam);
        I("e_merge", dp.lpf_h_merge);
        I("e_trunc", dp.lpf_h_trunc);
        I("e_trunc_in", dp.lpf_h_trunc_in);
        I("e_ents", dp.lp_h_ents);
        I("e_helpers", dp.lpf_h_helpers);
        I("e_clears", dp.lpf_h_clears);
        D("e_s", dp.lp_h_sec);
        D("e_root_s", dp.lpf_h_root_sec);
        D("e_help_s", dp.lpf_h_help_sec);
        D("e_merge_s", dp.lpf_h_merge_sec);
        I("f_hg_clears", dp.lpf_hg_clears);
        D("f_hg_s", dp.lpf_hg_sec);
        D("f_mtrim_s", dp.lpf_mtrim_sec);
        I("f_help_built", dp.lpf_help_built);
        I("f_help_distinct", static_cast<long long>(dp.lpf_help_keys.size()));
        I("f_help_rebuilt", dp.lpf_help_rebuilt);
        D("f_help_s", dp.lpf_help_sec);
        I("f_ml_nodes", dp.lpf_ml_nodes);
        I("f_ml_cand", dp.lpf_ml_cand);
        I("f_ml_max", dp.lpf_ml_max);
        I("f_ub_nodes", dp.ml_ub_nodes);
        I("f_ub_cand", dp.ml_ub_cand);
        I("f_ub_built", dp.ml_ub_built);
        I("f_ub_keys", static_cast<long long>(dp.prRHub.size()));
        I("f_racc_dense", dp.lpf_racc_dense);
        I("f_racc_fin_n", dp.lpf_racc_n);
        I("f_racc_fin", dp.lpf_racc_fin);
        I("f_rn_roots", dp.rn_roots);
        I("f_rn_nodes", dp.rn_nodes);
        I("f_stc_n", dp.lpf_stc_n);
        I("f_stc_buckets", dp.lpf_stc_buckets);
        I("f_stc_scan", dp.lpf_stc_scan);
        D("f_stc_s", dp.lpf_stc_sec);
        I("v_pairs", dp.vi_pairs);
        D("v_pairs_s", dp.vi_pairs_s);
        I("v_nodes_out", dp.vi_nodes[0]);
        I("v_nodes_in", dp.vi_nodes[1]);
        I("v_nb_int", dp.vi_nb_int);
        I("v_nb_ml", dp.vi_nb_ml);
        I("v_ml_nodes", dp.vi_ml_nodes);
        I("v_mlc", dp.vi_mlc);
        I("v_mlc_max", dp.vi_mlc_max);
        I("v_help_built", dp.pr_hbuilt[1]);
        I("v_help_terms", dp.pr_hterms[1]);
        I("v_help_ents", dp.pr_hentries[1]);
        I("v_add_int", dp.vi_add_int);
        I("v_add_ml", dp.vi_add_ml);
        I("v_add_help", dp.vi_help_adds);
        I("v_ent_fin", dp.vi_ent_fin);
        I("v_ent_kept", dp.vi_ent_kept);
        I("v_ats", dp.vi_ats);
        D("v_build_s", dp.vi_vbuild_s);
        D("v_steps_s", dp.vi_steps_s);
        D("v_terms_s", dp.vi_terms_s);
        D("v_help_s", dp.vi_help_s);
        D("v_add_int_s", dp.vi_add_int_s);
        D("v_add_ml_s", dp.vi_add_ml_s);
        D("v_fin_s", dp.vi_fin_s);
        D("v_ats_s", dp.vi_ats_s);
        D("v_arm_s", dp.vi_arm_s);
        D("v_rest_s", dp.vi_rest_s);
        D("v_init_s", dp.l2r_isec);
        D("v_prep_s", stem.pr_prep_sec);
        std::fprintf(stderr, "[l2r-perf] %s n %d%s\n", sid.c_str(), static_cast<int>(seq.size()), b.c_str());
      }
      if (dp.qz_audit || g_ntr_zbar)
        std::fprintf(stderr,
                     "[l2r-qz] %s n %d ntr_z %s ntr_first_j %d ntr_arm_skip %lld audit %d qz_reads %lld "
                     "qz_unfilled %lld qz_unfilled_pr %lld\n",
                     sid.c_str(), static_cast<int>(seq.size()), g_ntr_zbar ? "bar" : "scan", dp.ntr_first_j,
                     dp.ntr_arm_skip, static_cast<int>(dp.qz_audit), dp.qz_reads, dp.qz_unf, dp.qz_unf_pr);
      // per table: states before / after the beam (summed over steps), max kept per step,
      // build and score+prune seconds
      static const char* pn[DP03::PT_N] = {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us",
                                          "Uls", "Ulss", "Urs", "S", "P1", "P2", "P3"};
      std::string b;
      char buf[256];
      double bs = 0.0, ps = 0.0;
      for (int t = 0; t < (dp.l2r_inc ? DP03::PT_N : DP03::PT_P1); ++t) {
        std::snprintf(buf, sizeof buf, " %s %lld/%lld/%lld/%.3f/%.3f", pn[t], dp.l2r_pre[t],
                      dp.l2r_kept[t], dp.l2r_kmax[t], dp.l2r_bsec[t], dp.l2r_psec[t]);
        b += buf;
        bs += dp.l2r_bsec[t];
        ps += dp.l2r_psec[t];
      }
      std::fprintf(stderr,
                   "[l2r-beam] %s n %d beam %d spawn_k %d steps %d tables(pre/kept/max/build_s/"
                   "score_s)%s close_s %.3f build_s %.3f score_s %.3f init_s %.3f spawn events "
                   "%lld cand %lld kept %lld arms %d\n",
                   sid.c_str(), static_cast<int>(seq.size()), dp.l2r_beam, dp.l2r_spawn_k,
                   static_cast<int>(seq.size()), b.c_str(), dp.l2r_csec, bs, ps, dp.l2r_isec,
                   dp.sp_events, dp.sp_cand, dp.sp_kept, dp.l2r_beam > 0 ? dp.l2r_arms() : 0);
      dp.scan_ct_dump();
      if (dp.l2r_beam > 0 && !partial)
        std::fprintf(stderr,
                     "[l2r-spawn] %s n %d mode %s w %s k %d pool %d maxpop %d events %lld cand %lld "
                     "pops %lld kept %lld row_ents %zu row_mb %.3f col_ents %zu col_max %zu "
                     "cols_s %.3f rows_s %.3f heap_s %.3f rank_s %.3f spawn_s %.3f rs_mass %.9g\n",
                     sid.c_str(), static_cast<int>(seq.size()),
                     dp.l2r_inc ? "inc" : !dp.sp_cube ? "all" : dp.sp_static ? "cube-static" : "cube",
                     dp.sp_wk ? "k" : "suf",
                     dp.l2r_spawn_k, dp.sp_pool, dp.sp_maxpop > 0 ? dp.sp_maxpop : 4 * dp.sp_pool,
                     dp.sp_events, dp.sp_cand, dp.sp_pops, dp.sp_kept, dp.sp_row_ents,
                     dp.sp_row_ents * sizeof(int) / 1048576.0, dp.sp_col_ents, dp.sp_col_max,
                     dp.sp_col_sec, dp.sp_row_sec,
                     dp.sp_heap_sec, dp.sp_rank_sec, dp.l2r_bsec[DP03::PT_RS] + dp.l2r_psec[DP03::PT_RS],
                     dp.l2r_rs_mass());
      if (dp.l2r_beam > 0 && !partial && dp.l2r_ninf_mass)
        std::fprintf(stderr, "[l2r-ninf] %s n %d keep_ninf %d%s\n",
                     sid.c_str(), static_cast<int>(seq.size()), static_cast<int>(dp.l2r_keep_ninf),
                     dp.l2r_ninf_report().c_str());
      std::fprintf(stderr,
                   "[l2r-v] %s n %d trim %s eps %.3g cap %d v_before %lld v_after %lld v_capped %lld "
                   "arms %d arm_entries %zu arm_groups %zu arm_fb %d fb_hit %lld fb_floor %lld\n",
                   sid.c_str(), static_cast<int>(seq.size()), dp.pr_vsuffix ? "suffix" : "global",
                   dp.pr_veps, dp.pr_vcap, dp.v_before, dp.v_after, dp.v_capped, dp.l2r_arms(),
                   dp.l2r_arm_ents(), dp.l2r_arm_groups(), static_cast<int>(dp.l2r_arm_fb), dp.l2r_fb_hit,
                   dp.l2r_fb_floor);
      if (dp.pr_vmlk > 0)
        std::fprintf(stderr,
                     "[l2r-vk] %s n %d k %d nodes %lld cand %lld drop %lld resid_max %.4g resid_mean %.4g "
                     "mass_s %.3f init_s %.3f\n",
                     sid.c_str(), static_cast<int>(seq.size()), dp.pr_vmlk, dp.vk_nodes, dp.vk_cand,
                     dp.vk_drop, dp.vk_rmax, dp.vk_nodes > 0 ? dp.vk_rsum / dp.vk_nodes : 0.0,
                     dp.vk_mass_s, dp.l2r_isec);
      std::fprintf(stderr,
                   "[l2r-mem] %s n %d v_key %s v_stream %d v_pairs %lld v_entries %lld v_per_pair %.1f "
                   "v_live_max %lld v_sec %.3f arm_entries %zu node_trunc %d ntr_nodes %lld "
                   "ntr_L_in %lld ntr_L_kept %lld ntr_R_in %lld ntr_R_kept %lld ntr_R_K %d ntr_R_miss %lld "
                   "malloc_trim %d\n",
                   sid.c_str(), static_cast<int>(seq.size()), dp.pr_vxs ? "xs" : "xsys",
                   static_cast<int>(dp.pr_vstream), dp.v_pairs, dp.pr_V_entries,
                   dp.v_pairs > 0 ? static_cast<double>(dp.pr_V_entries) / dp.v_pairs : 0.0, dp.v_live_max,
                   dp.pr_static_sec, dp.l2r_arm_ents(), g_ntr_K, dp.ntr_nodes, dp.ntr_in[0], dp.ntr_kept[0],
                   dp.ntr_in[1], dp.ntr_kept[1], g_ntr_KR > 0 ? g_ntr_KR : g_ntr_K, dp.ntr_rmiss,
                   static_cast<int>(g_mtrim));
      if (dp.sw_as_rn_calls > 0 || dp.sw_fill_sec > 0.0) {
        const double mean_in =
            dp.sw_as_rn_calls > 0
                ? static_cast<double>(dp.sw_as_rn_ent_in) / static_cast<double>(dp.sw_as_rn_calls)
                : 0.0;
        const double mean_out =
            dp.sw_as_rn_calls > 0
                ? static_cast<double>(dp.sw_as_rn_ent_out) / static_cast<double>(dp.sw_as_rn_calls)
                : 0.0;
        const double mean_ga =
            dp.sw_ga_nodes > 0 ? static_cast<double>(dp.sw_ga_ent_sum) / static_cast<double>(dp.sw_ga_nodes)
                              : 0.0;
        const double work_per_out =
            mean_out > 0.0 ? mean_in / mean_out : 0.0;
        const char* gcut =
            g_cut == G_CUT_COLLAPSE ? "collapse" : (g_cut == G_CUT_ENTRIES ? "entries" : "topq");
        std::fprintf(stderr,
                     "[g-read] n %d bG %d g_cut %s ntr_init_sec %.6f sw_fill_sec %.6f "
                     "sw_as_rn_calls %lld sw_as_rn_sec %.6f ms_per_read %.4f ent_in_mean %.1f "
                     "ent_out_mean %.1f work_per_out %.1f ga_nodes %lld ga_ent_mean %.1f "
                     "ga_ent_max %lld rn_keys %lld grow_sec %.6f\n",
                     static_cast<int>(seq.size()), g_ntr_K, gcut, dp.sw_ntr_init_sec, dp.sw_fill_sec,
                     dp.sw_as_rn_calls, dp.sw_as_rn_sec,
                     dp.sw_as_rn_calls > 0 ? 1e3 * dp.sw_as_rn_sec / static_cast<double>(dp.sw_as_rn_calls)
                                          : 0.0,
                     mean_in, mean_out, work_per_out, dp.sw_ga_nodes, mean_ga, dp.sw_ga_ent_max,
                     dp.rn_keys, dp.pr_grow_sec);
        auto meanll = [](long long s, long long n) {
          return n > 0 ? static_cast<double>(s) / static_cast<double>(n) : 0.0;
        };
        std::fprintf(stderr,
                     "[vg-fill] V pushes %lld pre_mean %.1f pre_max %lld post_mean %.1f post_max %lld "
                     "merge_sec %.6f cut_sec %.6f | G pushes %lld pre_mean %.1f pre_max %lld "
                     "post_mean %.1f post_max %lld merge_sec %.6f cut_sec %.6f\n",
                     dp.v_pushes, meanll(dp.v_pre_sum, dp.v_pre_n), dp.v_pre_max,
                     meanll(dp.v_post_sum, dp.v_post_n), dp.v_post_max, dp.v_merge_sec, dp.v_cut_sec,
                     dp.g_pushes, meanll(dp.g_pre_sum, dp.g_pre_n), dp.g_pre_max,
                     meanll(dp.g_post_sum, dp.g_post_n), dp.g_post_max, dp.g_merge_sec, dp.g_cut_sec);
        auto yo_mean = [](long long s, long long n) {
          return n > 0 ? static_cast<double>(s) / static_cast<double>(n) : 0.0;
        };
        std::fprintf(stderr,
                     "[ml-yo] qm_ix %d S %d path %s pairs %lld yo_mean %.2f yo_max %lld | "
                     "scan_path pairs %lld yo_mean %.2f yo_max %lld\n",
                     dp.ml_yo_qm_ix ? 1 : 0, dp.ml_yo_S, dp.ml_yo_qm_ix ? "qm_ix" : "scan",
                     dp.ml_yo_qm_n, yo_mean(dp.ml_yo_qm_sum, dp.ml_yo_qm_n), dp.ml_yo_qm_max,
                     dp.ml_yo_sc_n, yo_mean(dp.ml_yo_sc_sum, dp.ml_yo_sc_n), dp.ml_yo_sc_max);
        double vf_med = 0.0;
        if (!dp.arm_Vf_samples.empty()) {
          auto v = dp.arm_Vf_samples;
          std::nth_element(v.begin(), v.begin() + static_cast<long>(v.size() / 2), v.end());
          vf_med = v[v.size() / 2];
        }
        std::fprintf(stderr,
                     "[kl-m] read_sec %.6f kl_in_read_sec %.6f arms %lld Vf_mean %.1f Vf_med %.1f "
                     "Vf_max %lld Kq_mean %.1f Kq_max %lld\n",
                     dp.sw_as_rn_sec, dp.kl_in_read_sec, dp.arm_read_n,
                     meanll(dp.arm_Vf_sum, dp.arm_read_n), vf_med, dp.arm_Vf_max,
                     meanll(dp.arm_Kq_sum, dp.arm_read_n), dp.arm_Kq_max);
      }
      if (g_mem_log) {
        static const char* pn[DP03::PT_N] = {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us",
                                             "Uls", "Ulss", "Urs", "S", "P1", "P2", "P3"};
        std::string tb;
        char buf[160];
        for (int t = 0; t < DP03::PT_N; ++t) {
          std::snprintf(buf, sizeof buf, " %s %lld/%lld/%lld", pn[t], dp.l2r_nsc[t], dp.l2r_ninf[t],
                        dp.l2r_noarm[t]);
          tb += buf;
        }
        long long re = 0, ve = 0;
        const double rb = dp.rside_bytes(re), vb = dp.rfv_bytes(ve);
        auto mean = [](long long s, long long c) { return c > 0 ? static_cast<double>(s) / c : 0.0; };
        std::fprintf(stderr,
                     "[l2r-mem2] %s n %d v_cap_pairs %lld v_cap_entries %lld v_cap_frac_mean %.3g "
                     "v_cap_frac_max %.3g v_cap_vfrac_max %.3g v_cap_vfrac_tot %.3g R_node %.1f/%lld "
                     "R_root %.1f/%lld R_key %.1f/%lld R_entries %lld R_mb %.1f rfv %d V_keep %lld "
                     "V_keep_mb %.1f rfv_look %lld rfv_miss %lld rfv_seed %lld rfv_od %lld rfv_od_ent %lld "
                     "rfv_child_miss %lld rfv_keys %lld rfv_key_ent %lld | scored/ninf/noarm%s\n",
                     sid.c_str(), static_cast<int>(seq.size()), dp.v_cap_pairs, dp.v_capped,
                     dp.v_cap_pairs > 0 ? dp.v_cap_fsum / dp.v_cap_pairs : 0.0,
                     dp.v_cap_fmax, dp.v_cap_fwmax,
                     dp.v_cap_pairs > 0 ? std::exp(dp.v_cap_drop - dp.v_cap_all) : 0.0,
                     mean(dp.rz_sum[0], dp.rz_n[0]), dp.rz_max[0], mean(dp.rz_sum[1], dp.rz_n[1]),
                     dp.rz_max[1], mean(dp.rz_sum[2], dp.rz_n[2]), dp.rz_max[2], re, rb / 1048576.0,
                     static_cast<int>(g_rfv), ve, vb / 1048576.0, dp.rfv_look, dp.rfv_miss,
                     dp.rfv_seed, dp.rfv_od, dp.rfv_od_ent, dp.rfv_child_miss, dp.rfv_keys,
                     dp.rfv_key_ent, tb.c_str());
      }
    }
    dp.pr_v_finish(sid);
    if (pruned_f) {
      for (int t = 0; t < BT_N; ++t) {
        std::fprintf(pruned_f, "P %s %s events %lld dropped %lld matched %lld w %.9e pos",
                     sid.c_str(), BT_NAMES[t], drops.events[t], drops.dropped[t],
                     drops.matched[t], drops.w[t]);
        for (int k = 0; k < 10; ++k) std::fprintf(pruned_f, " %.6e", drops.wpos[t][k]);
        std::fprintf(pruned_f, "\n");
      }
      std::fflush(pruned_f);
      if (kept_f) std::fflush(kept_f);
    }
    if (gcfg.on || gate_stats_out) {
      const GateStats& gs = dp.gate_stats;
      std::string body;
      const char* names[2] = {"Qg", "Qgl"};
      for (int t = 0; t < 2; ++t)
        body += std::string(" ") + names[t] + " cand " + std::to_string(gs.cand[t]) + " elig " +
                std::to_string(gs.elig[t]) + " ranked " + std::to_string(gs.ranked[t]) +
                " kept " + std::to_string(gs.kept[t]) + " reserved " +
                std::to_string(gs.reserved[t]);
      std::fprintf(stderr, "[gate] %s n %d on %d tau %.17g span %d%s\n",
                   id.empty() ? "-" : id.c_str(), static_cast<int>(seq.size()), gcfg.on ? 1 : 0,
                   gcfg.tau, gcfg.span, body.c_str());
    }
    if (arm_score) {
      const ArmStats& st = dp.arm_stats;
      std::string hist;
      for (int k = 0; k < LK_HIST_BINS; ++k)
        hist += (k ? "," : "") + std::to_string(st.lk_hist[k]);
      std::fprintf(stderr,
                   "[arm] %s n %d scored %lld arm_wins %lld frac %.6f mean_log_kappa %.6f "
                   "max_log_kappa %.6f windows %lld floor_windows %lld chosen_floor %lld "
                   "lk_hist %s\n",
                   id.empty() ? "-" : id.c_str(), static_cast<int>(seq.size()), st.scored,
                   st.arm_wins, st.scored ? static_cast<double>(st.arm_wins) / st.scored : 0.0,
                   st.arm_wins ? st.sum_log_kappa / st.arm_wins : 0.0, st.max_log_kappa,
                   st.windows, st.floor_windows, st.chosen_floor, hist.c_str());
    }
    if (stem_rank)
      std::fprintf(stderr,
                   "[stem] %s n %d loop %d role %s entries %zu dropped %zu nodes %zu bytes %zu "
                   "removed_qg %lld removed_qgl %lld removed_qgr %lld sec %.6f\n",
                   sid.c_str(), static_cast<int>(seq.size()), stem.loop,
                   stem_qg2 ? "all" : stem_both ? "both" : "stem1", stem.entries, stem.dropped,
                   stem.index.nodes() + (stem.mirror ? stem.index2.nodes() : 0),
                   stem.bytes(), dp.stem_removed[0], dp.stem_removed[1], dp.stem_removed[2],
                   stem_sec);
    if (stem_rank && (stem_hr_av || stem_hr_phi))
      std::fprintf(stderr,
                   "[stem-hr] %s n %d mode %s gamma %.6g entries %zu dropped %zu nodes %zu "
                   "qg_states %lld av_only %lld orig_only %lld mean_delta %.9g "
                   "min_delta %.9g max_delta %.9g check_states %lld check_mismatch %lld "
                   "check_max_rel %.3e\n",
                   sid.c_str(), static_cast<int>(seq.size()), stem_hr_av ? "av" : "phi",
                   stem_hr_gamma, stem.entries_av,
                   stem.dropped_av, stem.index_av.nodes(), dp.hr_n, dp.hr_av_only,
                   dp.hr_orig_only, dp.hr_n ? dp.hr_dsum / dp.hr_n : 0.0, dp.hr_dmin,
                   dp.hr_dmax, dp.hrc_states, dp.hrc_mis, dp.hrc_rel);
    if (stem_rank && (stem_qg_r_set || stem_consts_set)) {
      // Pseudoknot factors each role adds on top of alpha and the outside estimate (H_R,
      // the mirrored H_L, or the partner tables); at_in = at(u, v), at_out = at(p, q).
      const char* qg_l = stem_consistent ? "beta1+2beta2+at_in+at_out" : "beta1+2beta2";
      const char* qg_r = stem_qg_partner
                             ? (stem_consistent ? "Bg[beta2+at_in]+beta1+at_out" : "Bg[beta2+at_in]")
                             : (stem_consistent ? "beta1+2beta2+at_in+at_out" : "beta1+2beta2");
      const char* gl_l = stem_consistent ? "beta1+beta2+at_out" : "beta1+2beta2";
      std::fprintf(stderr,
                   "[stem-consts] %s n %d qg_r %s consts %s Qg.L %s Qg.R %s Qgl.L %s Qgl.R %s "
                   "Qgr.R %s qg_score_sec %.6f\n",
                   sid.c_str(), static_cast<int>(seq.size()),
                   stem_qg_partner ? "partner" : "mirror", stem_consistent ? "consistent" : "as-is",
                   qg_l, stem_qg2 ? qg_r : "-", gl_l, stem_both ? "beta1+at_out" : "-",
                   stem_both ? "beta1+at_out" : "-", dp.stem_qg_sec);
    }
    if (stem_rank && stem_ext) {
      // loop / helpers only: growth and helper counts (Result 19)
      char grow[1024] = "";
      if (stem_ext == 3 || stem_helpers) {
        const long long hs = dp.hlp_states[0] + dp.hlp_states[1];
        const double lay = static_cast<double>(dp.grow_layer);
        std::snprintf(
            grow, sizeof grow,
            " lg %d cut %.3g cap %d helpers %s ml_span %d ranked_qg %lld ranked_qgls %lld "
            "ranked_qgrs %lld lists %lld list_nodes %lld nodes_per_list %.4f layer_nodes %lld "
            "cut_drop %lld cap_drop %lld cut_frac %.6f cap_frac %.6f helper_closings %lld "
            "closings_per_helper %.4f helper_nodes %lld helper_finite %lld helper_sec %.6f",
            stem_ext == 3 ? grow_lg : -1, grow_cut, grow_nodes,
            stem_helpers == 1 ? "grow" : stem_helpers == 2 ? "nested" : "off",
            ml_span, dp.ext_states, dp.hlp_states[0], dp.hlp_states[1], dp.grow_lists,
            dp.grow_list_nodes,
            dp.grow_lists ? static_cast<double>(dp.grow_list_nodes) / dp.grow_lists : 0.0,
            dp.grow_layer, dp.grow_cut_drop, dp.grow_cap_drop,
            lay > 0 ? dp.grow_cut_drop / lay : 0.0, lay > 0 ? dp.grow_cap_drop / lay : 0.0,
            dp.hlp_closings, hs ? static_cast<double>(dp.hlp_closings) / hs : 0.0,
            dp.hlp_nodes, dp.hlp_finite, dp.hlp_sec);
      }
      // --stem-est-profile (Result 21): keys per step, mean key-profile support, entries
      // evaluated per scored state, Bg lookups / raw computations / partner entries visited,
      // profile build and per-state scoring time
      char prof[1024] = "";
      if (stem_profile) {
        const long long keys = dp.pf_keys[0] + dp.pf_keys[1] + dp.pf_keys[2];
        std::snprintf(
            prof, sizeof prof,
            " profile on eps %.3g tail %s steps %lld keys_qg %lld keys_qgls %lld keys_qgrs %lld "
            "keys_per_step %.4f pair_profiles %lld mean_ly %.4f mean_rx %.4f trunc_drop %lld "
            "scored %lld evals_per_state %.4f tail_closings %lld bg_lookups %lld bg_misses %lld "
            "bg_visits %lld bg_visits_per_miss %.4f prof_build_sec %.6f prof_score_sec %.6f",
            prof_eps, ml_tail ? "nested" : "off", dp.pf_steps, dp.pf_keys[0], dp.pf_keys[1],
            dp.pf_keys[2], dp.pf_steps ? static_cast<double>(keys) / dp.pf_steps : 0.0,
            dp.pf_pairs, keys ? static_cast<double>(dp.pf_ly) / keys : 0.0,
            keys ? static_cast<double>(dp.pf_rx) / keys : 0.0, dp.pf_trunc, dp.pf_states,
            dp.pf_states ? static_cast<double>(dp.pf_evals) / dp.pf_states : 0.0, dp.pf_tail_cl,
            dp.bg_lookups, dp.bg_misses, dp.bg_visits,
            dp.bg_misses ? static_cast<double>(dp.bg_visits) / dp.bg_misses : 0.0,
            dp.pf_build_sec, dp.pf_total_sec - dp.pf_build_sec);
        if (bg_reorder) {
          const size_t len = std::strlen(prof);
          std::snprintf(prof + len, sizeof prof - len,
                        " bg_reorder on m_rows %lld m_build_visits %lld m_eval_visits %lld",
                        dp.ro_rows, dp.ro_build, dp.ro_eval);
        }
      }
      // --stem-est-chain (Result 24): kept pairs, chain nodes, steps by kind, root profiles and
      // their mean support, entries cut; chain time (part of prof_build_sec), nested pass time
      char chn[512] = "";
      if (stem_chain)
        std::snprintf(
            chn, sizeof chn,
            " chain on floor %.3g span %d loop %d cut %.3g ml %s kept %zu chain_nodes %lld "
            "edges_b %lld edges_f %lld edges_g %lld edges_h %lld roots %lld root_ly %.4f "
            "root_rx %.4f chain_cut_drop %lld chain_sec %.6f chain_prep_sec %.6f",
            chain_floor, chain_span, chain_loop, chain_cut, chain_ml ? "on" : "off", stem.n_kept,
            dp.ch_nodes, dp.ch_edges[0], dp.ch_edges[1], dp.ch_edges[2], dp.ch_edges[3],
            dp.ch_roots, dp.ch_roots ? static_cast<double>(dp.ch_ly) / dp.ch_roots : 0.0,
            dp.ch_roots ? static_cast<double>(dp.ch_rx) / dp.ch_roots : 0.0, dp.ch_cut_drop,
            dp.ch_sec, chain_prep_sec);
      // --stem-est-principled: ranked states and scoring time per table, key-profile entries
      // (G endpoints) per state, arms per T, T / K / Stop builds and hits, per-q' arm
      // structures, static sizes, pre-pass kept pairs and time, Bml check
      std::string prs;
      if (stem_pr) {
        long long st_all = 0;
        for (int t = 0; t < 6; ++t) st_all += dp.pr_states[t];
        const long long st4 = st_all - dp.pr_states[5];
        char b[2048];
        std::snprintf(
            b, sizeof b,
            " principled on floor %.3g span %d ml_span %d ml_close %d cut %.3g arms_k %d "
            "kept %zu states_qg %lld states_qgl %lld states_qgr %lld states_qgls %lld "
            "states_qgrs %lld states_qp %lld sec_qg %.6f sec_qgl %.6f sec_qgr %.6f "
            "sec_qgls %.6f sec_qgrs %.6f sec_qp %.6f endpoints_per_state %.4f "
            "arms_per_T %.4f arms_skipped %lld ys_per_T %.4f T_built %lld T_hits %lld "
            "K_built %lld K_hits %lld i_per_K %.4f stop_built %lld stop_hits %lld "
            "armq_built %lld V_entries %lld arms %lld arm_entries %lld G_nodes %lld "
            "G_entries %lld roots %lld keys %lld key_entries %lld cut_drop %lld "
            "static_sec %.6f pr_prep_sec %.6f bml_check_cells %lld bml_check_mis %lld "
            "bml_check_rel %.3e",
            pr_floor, 0, pr_ml_span, pr_ml_close, pr_cut, pr_arms_k, stem.pr_n_kept,
            dp.pr_states[0], dp.pr_states[1], dp.pr_states[2], dp.pr_states[3],
            dp.pr_states[4], dp.pr_states[5], dp.pr_sec[0], dp.pr_sec[1], dp.pr_sec[2],
            dp.pr_sec[3], dp.pr_sec[4], dp.pr_sec[5],
            st4 ? static_cast<double>(dp.pr_evals) / st4 : 0.0,
            dp.pr_T_built ? static_cast<double>(dp.pr_T_arms) / dp.pr_T_built : 0.0,
            dp.pr_T_skip,
            dp.pr_T_built ? static_cast<double>(dp.pr_T_ys) / dp.pr_T_built : 0.0, dp.pr_T_built,
            dp.pr_T_hits, dp.pr_K_built, dp.pr_K_hits,
            dp.pr_K_built ? static_cast<double>(dp.pr_K_is) / dp.pr_K_built : 0.0,
            dp.pr_slr_miss, dp.pr_slr_hits, dp.pr_aq_built, dp.pr_V_entries, dp.pr_arms,
            dp.pr_arm_entries, dp.pr_G_nodes, dp.pr_G_entries, dp.pr_roots, dp.pr_keys,
            dp.pr_key_entries, dp.pr_cut_drop, dp.pr_static_sec, stem.pr_prep_sec,
            stem.bml_check_cells, stem.bml_check_mis, stem.bml_check_rel);
        prs = b;
        if (stem.pr_partner) {
          std::snprintf(b, sizeof b, " partner_floor %.3g partner_pairs %zu", pr_partner_floor,
                        stem.pr_n_partner);
          prs += b;
        }
        if (pr_qz_unp) prs += " test_qz unpaired";
        if (pr_pmax > 0) {
          std::snprintf(b, sizeof b,
                        " partner_maxpairs %d V_pairs %lld V_layer_entries %lld V_map_per_pair %.4f"
                        " V_layers_per_pair %.4f pcheck_n %lld pcheck_key_mis %lld pcheck_maxd %.3e"
                        " penum_n %lld penum_key_mis %lld penum_maxd %.3e",
                        pr_pmax, dp.pr_V_pairs, dp.pr_V_layer_entries,
                        dp.pr_V_pairs ? static_cast<double>(dp.pr_V_entries) / dp.pr_V_pairs : 0.0,
                        dp.pr_V_pairs ? static_cast<double>(dp.pr_V_layer_entries) / dp.pr_V_pairs : 0.0,
                        dp.pc_n, dp.pc_key_mis, dp.pc_maxd, dp.pe_n, dp.pe_key_mis, dp.pe_maxd);
          prs += b;
        }
        if (!pr_pml) {
          std::snprintf(b, sizeof b, " partner_ml off pmlcheck_n %lld pmlcheck_key_mis %lld"
                        " pmlcheck_maxd %.3e", dp.pm_n, dp.pm_key_mis, dp.pm_maxd);
          prs += b;
        }
        {
          auto avg = [](long long a, long long c) { return c ? static_cast<double>(a) / c : 0.0; };
          auto shr = [](double a, long long c) { return c ? a / c : 0.0; };
          std::snprintf(
              b, sizeof b,
              " grow_select %s grow_eps %.3g grow_B %d grow_sec %.6f sel_sec %.6f sel_profiles "
              "%lld sel_before %.4f sel_after %.4f sel_iota_share %.6f sel_zero %lld "
              "node_sel_before %.4f node_sel_after %.4f node_iota_share %.6f node_zero %lld "
              "uL_cols %lld uR_entries %lld U_evals %lld",
              pr_gsel ? "bound" : "floor", pr_geps, pr_gB, dp.pr_grow_sec, dp.pr_sel_sec,
              dp.gs_prof[1], avg(dp.gs_before[1], dp.gs_prof[1]),
              avg(dp.gs_after[1], dp.gs_prof[1]), shr(dp.gs_share[1], dp.gs_prof[1]),
              dp.gs_zero[1], avg(dp.gs_before[0], dp.gs_prof[0]),
              avg(dp.gs_after[0], dp.gs_prof[0]), shr(dp.gs_share[0], dp.gs_prof[0]),
              dp.gs_zero[0], dp.gs_cols, dp.gs_rent, dp.gs_U);
          prs += b;
          std::snprintf(b, sizeof b,
                        " grow_built %lld grow_int %lld grow_ml %lld fut_built %lld fut_int %lld "
                        "fut_ml %lld partner_sec %.6f aq_sep %d aq_sec %.6f aq_max %zu "
                        "T_sec %.6f K_sec %.6f SLR_sec %.6f W_sec %.6f W_built %lld W_ent %lld",
                        dp.st_nodes[0], dp.st_kind[0][0], dp.st_kind[0][1], dp.st_nodes[1],
                        dp.st_kind[1][0], dp.st_kind[1][1], dp.pr_partner_sec,
                        static_cast<int>(dp.pr_aq_sep), dp.pr_aq_sec, dp.pr_aq_max, dp.pr_sec_T,
                        dp.pr_sec_K, dp.pr_sec_SLR, dp.pr_sec_W, dp.pr_W_built, dp.pr_W_ent);
          prs += b;
          const double pidx_mb = dp.pidx_mb;
          std::snprintf(b, sizeof b,
                        " arm_persist %d arm_gaponly %d partner_span %d v_eps %.3g pidx_sec %.6f "
                        "pidx_nodes %lld pidx_mb %.3f pidx_upd %lld aq_arms %lld "
                        "aq_arms_per_build %.4f gap_epochs %lld gap_miss %lld v_before %lld "
                        "v_after %lld",
                        static_cast<int>(pr_arm_persist), static_cast<int>(pr_arm_gaponly), pr_pspan,
                        pr_veps, dp.pidx_sec, dp.pidx_nodes, pidx_mb, dp.pidx_upd, dp.aq_arms,
                        avg(dp.aq_arms, dp.pr_aq_built), dp.gap_epochs, dp.gap_miss, dp.v_before,
                        dp.v_after);
          prs += b;
          std::snprintf(b, sizeof b,
                        " fast %d hg_built %lld hg_entries %lld hg_terms %lld hg_clears %lld "
                        "hv_built %lld hv_entries %lld hv_terms %lld node_terms_g %lld "
                        "node_terms_v %lld",
                        pr_fast, dp.pr_hbuilt[0], dp.pr_hentries[0], dp.pr_hterms[0],
                        dp.pr_hg_clears, dp.pr_hbuilt[1], dp.pr_hentries[1], dp.pr_hterms[1],
                        dp.pr_nterms[0], dp.pr_nterms[1]);
          prs += b;
        }
        if (pr_arms_mode == 1) {
          const long long steps = dp.pr_sep_steps_all();
          std::snprintf(
              b, sizeof b,
              " arms sep arm_z %s sep_queries %lld sep_cols %lld sep_groups %lld "
              "groups_per_step %.4f arms_per_group %.4f queries_per_group %.4f sep_sec %.6f",
              pr_arm_z == 1 ? "unpaired" : "qz", dp.pr_sep_q, dp.pr_sep_cols, dp.pr_sep_groups,
              steps ? static_cast<double>(dp.pr_sep_groups) / steps : 0.0,
              dp.pr_sep_groups ? static_cast<double>(dp.pr_sep_garms) / dp.pr_sep_groups : 0.0,
              dp.pr_sep_groups ? static_cast<double>(dp.pr_sep_q) / dp.pr_sep_groups : 0.0,
              dp.pr_sep_sec);
          prs += b;
        }
      }
      std::fprintf(stderr,
                   "[stem-ext] %s n %d ext %s kmax %d combine %s states %lld nodes %lld "
                   "nodes_per_state %.4f finite %lld ext_sec %.6f qg_score_sec %.6f%s%s%s%s\n",
                   sid.c_str(), static_cast<int>(seq.size()),
                   stem_ext == 1 ? "stack" : stem_ext == 2 ? "bulge1" : stem_ext == 3 ? "loop" : "off",
                   stem_ext_max, stem_split ? "split" : "sum", dp.ext_states, dp.ext_nodes,
                   dp.ext_states ? static_cast<double>(dp.ext_nodes) / dp.ext_states : 0.0,
                   dp.ext_finite, dp.ext_sec, dp.stem_qg_sec, grow, prof, chn, prs.c_str());
      if (stem_pr) {
        auto qs = [](std::vector<int> v) {
          char t[96];
          if (v.empty()) {
            std::snprintf(t, sizeof t, "0/0/0");
            return std::string(t);
          }
          std::sort(v.begin(), v.end());
          double s = 0;
          for (int x : v) s += x;
          std::snprintf(t, sizeof t, "%.1f/%d/%d", s / static_cast<double>(v.size()),
                        v[static_cast<size_t>(0.9 * static_cast<double>(v.size() - 1))], v.back());
          return std::string(t);
        };
        std::string g = g_grow_rec ? "rec" : "prof";
        char b[1024];
        std::snprintf(b, sizeof b,
                      " ml_L %d ml_R %d ml_c %d diag %d U_nobml %d bml_built %d grow_sec %.6f "
                      "sel_sec %.6f pr_prep_sec %.6f partner_sec %.6f mem_nodes_mb %.3f "
                      "mem_roots_mb %.3f mem_help_peak_mb %.3f mem_step_peak_mb %.3f "
                      "mem_total_peak_mb %.3f help_entries_peak %lld help_clears %lld",
                      g_ml_L, g_ml_R, g_ml_c, static_cast<int>(g_ml_diag),
                      static_cast<int>(g_U_nobml), static_cast<int>(!stem.prBml.empty()),
                      dp.pr_grow_sec, dp.pr_sel_sec, stem.pr_prep_sec, dp.pr_partner_sec,
                      dp.gm_nodes / 1048576.0, dp.gm_roots / 1048576.0, dp.gm_help_peak / 1048576.0,
                      dp.gm_step_peak / 1048576.0, dp.gm_total_peak / 1048576.0,
                      dp.gm_help_entries_peak, g_grow_rec ? dp.rh_clears : dp.pr_hg_clears);
        g += b;
        const char* nm[3] = {"node", "root", "key"};
        for (int w = 0; w < 3; ++w) {
          g += std::string(" ") + nm[w] + "_n " + std::to_string(dp.gz_after[w].size());
          if (!g_grow_rec) g += std::string(" ") + nm[w] + "_before " + qs(dp.gz_before[w]);
          g += std::string(" ") + nm[w] + (g_grow_rec ? "_elems " : "_after ") + qs(dp.gz_after[w]);
          g += std::string(" ") + nm[w] + "_dp " + qs(dp.gz_dp[w]);
          g += std::string(" ") + nm[w] + "_dq " + qs(dp.gz_dq[w]);
          if (!g_grow_rec) {
            g += std::string(" ") + nm[w] + "_dp_before " + qs(dp.gz_dpb[w]);
            g += std::string(" ") + nm[w] + "_dq_before " + qs(dp.gz_dqb[w]);
          }
        }
        if (g_grow_rec) {
          std::snprintf(b, sizeof b,
                        " rn_nodes %lld rn_roots %lld rn_keys %lld rn_terms %lld rh_built %lld "
                        "rh_terms %lld mlc_drop_left %lld mlc_drop_right %lld",
                        dp.rn_nodes, dp.rn_roots, dp.rn_keys, dp.rn_terms, dp.rh_built, dp.rh_terms,
                        dp.rn_mlc_drop_l, dp.rn_mlc_drop_r);
          g += b;
        }
        if (g_tr_K > 0) {
          std::snprintf(b, sizeof b, " trunc_K %d trunc_eps %g trunc_mode %s", g_tr_K, g_tr_eps,
                        g_tr_pool ? "pool" : "side");
          g += b;
          const char* tn[2] = {"root", "key"};
          for (int w = 0; w < 2; ++w) {
            const double nn = static_cast<double>(dp.tr_before[w].size());
            g += std::string(" ") + tn[w] + "_tr_n " + std::to_string(dp.tr_before[w].size());
            g += std::string(" ") + tn[w] + "_tr_before " + qs(dp.tr_before[w]);
            g += std::string(" ") + tn[w] + "_tr_after " + qs(dp.tr_after[w]);
            g += std::string(" ") + tn[w] + "_tr_L " + qs(dp.tr_kL[w]);
            g += std::string(" ") + tn[w] + "_tr_R " + qs(dp.tr_kR[w]);
            std::snprintf(b, sizeof b, " %s_tr_share %.8f %s_tr_zero %.0f", tn[w],
                          nn > 0 ? dp.tr_share[w] / nn : 1.0, tn[w], dp.tr_zero[w]);
            g += b;
          }
        }
        if (g_rec_check) {
          std::snprintf(b, sizeof b, " check_states %lld check_max_dscore %.3e check_entries %lld "
                        "check_max_dvec %.3e", dp.rc_n, dp.rc_max_score, dp.rc_vec_n, dp.rc_max_vec);
          g += b;
        }
        if (g_ml_diag && dp.rd_n > 0) {
          const char* cn[DP03::RD_CFG] = {"all_all", "all_none", "depth1_none", "depth1_depth1",
                                          "none_none"};
          std::snprintf(b, sizeof b, " diag_states %lld", dp.rd_n);
          g += b;
          for (int c = 0; c < DP03::RD_CFG; ++c) {
            // kept-fraction quantiles from the -log10 histogram (bin lower edges)
            auto q = [&](double frac) {
              const long long goal = static_cast<long long>(std::ceil(frac * static_cast<double>(dp.rd_n)));
              long long cum = 0;
              for (int k = DP03::RD_BINS; k >= 0; --k) {
                cum += dp.rd_hist[c][k];
                if (cum >= goal) return k == DP03::RD_BINS ? 0.0 : std::pow(10.0, -k / 100.0);
              }
              return 1.0;
            };
            std::snprintf(b, sizeof b,
                          " kept_%s_mean %.6f kept_%s_p50 %.4g kept_%s_p10 %.4g kept_%s_p1 %.4g "
                          "kept_%s_min %.3g",
                          cn[c], dp.rd_sum[c] / static_cast<double>(dp.rd_n), cn[c], q(0.5), cn[c],
                          q(0.1), cn[c], q(0.01), cn[c], dp.rd_min[c]);
            g += b;
          }
          for (int side = 0; side < 2; ++side)
            for (int m = 0; m < 3; ++m) {
              std::snprintf(b, sizeof b, " share_%s%d %.6f", side ? "R" : "L", m,
                            dp.rd_layer_sum[side][m] / static_cast<double>(dp.rd_n));
              g += b;
            }
        }
        std::fprintf(stderr, "[growth] %s n %d mode %s\n", sid.c_str(), static_cast<int>(seq.size()),
                     g.c_str());
      }
    }
    if (partial) {
      std::fprintf(stderr, "[partial] %s n %d stopped_at_step %d phase %s\n", sid.c_str(),
                   static_cast<int>(seq.size()), dp.stop_j, dp.stop_phase == 1 ? "vbuild" : "scan");
      std::fflush(stderr);
      return 3;
    }
    if (pr_gcheck)
      std::fprintf(stderr,
                   "[pr-grow-check] %s n %d endpoints %lld viol_L %lld viol_R %lld "
                   "max_log_ratio_L %.6e max_log_ratio_R %.6e sel_above_all %lld\n",
                   sid.c_str(), static_cast<int>(seq.size()), dp.gc_n, dp.gc_vL, dp.gc_vR,
                   dp.gc_maxL, dp.gc_maxR, dp.gc_sel_gt);
    if (pr_check)
      std::fprintf(stderr,
                   "[pr-check] %s n %d states %lld/%lld/%lld/%lld/%lld/%lld mismatch "
                   "%lld/%lld/%lld/%lld/%lld/%lld max_rel %.3e/%.3e/%.3e/%.3e/%.3e/%.3e "
                   "bml_cells %lld bml_mismatch %lld bml_max_rel %.3e\n",
                   sid.c_str(), static_cast<int>(seq.size()), dp.prc_states[0], dp.prc_states[1],
                   dp.prc_states[2], dp.prc_states[3], dp.prc_states[4], dp.prc_states[5],
                   dp.prc_mis[0], dp.prc_mis[1], dp.prc_mis[2], dp.prc_mis[3], dp.prc_mis[4],
                   dp.prc_mis[5], dp.prc_rel[0], dp.prc_rel[1], dp.prc_rel[2], dp.prc_rel[3],
                   dp.prc_rel[4], dp.prc_rel[5], stem.bml_check_cells, stem.bml_check_mis,
                   stem.bml_check_rel);
    if (ro_check)
      std::fprintf(stderr,
                   "[bg-reorder-check] %s n %d states %lld mismatch %lld max_rel %.3e events %lld "
                   "flip_events %lld flip_states %lld\n",
                   sid.c_str(), static_cast<int>(seq.size()), dp.roc_states, dp.roc_mis,
                   dp.roc_rel, dp.roc_events, dp.roc_flip_events, dp.roc_flip_states);
    if (prof_check)
      std::fprintf(stderr,
                   "[profile-check] %s n %d qg_states %lld qg_mismatch %lld qg_max_rel %.3e "
                   "helper_states %lld helper_mismatch %lld helper_max_rel %.3e\n",
                   sid.c_str(), static_cast<int>(seq.size()), dp.pc_states[0], dp.pc_mis[0],
                   dp.pc_rel[0], dp.pc_states[1], dp.pc_mis[1], dp.pc_rel[1]);
    if (stem_ext_check) {
      char hc[256] = "";
      if (stem_helpers)
        std::snprintf(hc, sizeof hc, " helper_states %lld helper_mismatch %lld helper_max_rel %.3e",
                      dp.hc_states, dp.hc_mis, dp.hc_rel);
      std::fprintf(stderr, "[ext-check] %s n %d states %lld mismatch %lld max_rel %.3e%s\n",
                   sid.c_str(), static_cast<int>(seq.size()), dp.xc_states, dp.xc_mis, dp.xc_rel,
                   hc);
    }
    std::vector<std::tuple<int, int, double>> P;
    if (do_bpp) P = dp.bpp(cutoff);
    if (g_mem_log) {
      std::vector<MemItem> it;
      dp.mem_items(it, "", true, pk_l2r != 0);
      mem_log("outside_done", it);
    }
    if (dp.l2r_diag) {
      dp.outside();
      if (dp.l2r_split) dp.l2r_split_run();
      std::fprintf(stderr, "[l2r-diag] %s n %d logZ %.12f%s\n", sid.c_str(), static_cast<int>(seq.size()),
                   logz, dp.l2r_diag_report({200, 500, 1000, 2000, 4000}).c_str());
    }
    if (!l2r_odump.empty() && pk_l2r == 1) {
      if (!dp.l2r_diag) dp.outside();
      dp.l2r_oracle_dump(l2r_odump);
    }
    if (dp.ltr_f) {
      dp.outside();
      dp.ltr_final();
    }
    if (stem_bg_check) std::fprintf(stderr, "%s\n", dp.stem_bg_check(sid, logz).c_str());
    double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (stem_diag) {
      std::vector<const StemEst*> ests;
      for (const StemEst& s : diag_ests) ests.push_back(&s);
      std::string js = dp.stem_diag(ests, sid);
      std::fprintf(stem_diag_f, "%s\n", js.c_str());
      std::fflush(stem_diag_f);
    }
    if (gate_diag) {
      std::string js = dp.gate_diag(*gate, diag_taus, id.empty() ? "-" : id, gcfg.span,
                                    gcfg.delta, long_span);
      std::fprintf(diag_f, "%s\n", js.c_str());
      std::fflush(diag_f);
    }
    if (states_f) dp.dump_states(states_f, sid, dump_min_w);
    if (oracle_f) dp.pr_oracle_diag(stem, oracle_f, sid, pr_oracle_set);
    if (pr_sep_check) dp.pr_sep_report(stderr, sid);
    if (scores_f) dp.dump_scores(scores_f, sid, stem, dump_min_w);
    if (ext_f) {
      auto it = exact_states.find(sid);
      if (it == exact_states.end()) std::fprintf(stderr, "--dump-exact: no states for %s\n", sid.c_str());
      dp.dump_ext(ext_f, sid, stem, it == exact_states.end() ? no_states : it->second,
                  stem_ext_max);
    }

    if (const char* p2dump = std::getenv("DUMP_P2_ALL")) {
      dp.outside();
      dp.dump_p2_all(p2dump);
    }
    if (!id.empty()) std::printf(">%s\n", id.c_str());
    std::printf("logZ %.17g\n", logz);
    std::printf("sec %.6f\n", sec);
    std::printf("pairs %zu\n", P.size());
    for (const auto& t : P)
      std::printf("%d %d %.17g\n", std::get<0>(t), std::get<1>(t), std::get<2>(t));
    if (dump_qb_on) dp.dump_qb(stdout);
    std::fflush(stdout);
    id.clear();
  }
  if (diag_f) std::fclose(diag_f);
  if (stem_diag_f) std::fclose(stem_diag_f);
  g_beam_stats.close();
  return 0;
}
