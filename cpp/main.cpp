// Entry point of the unity build. main.cpp #includes the parts below in the order of the
// original single file (tag paper-v1, cpp/dp03_l2r_merged.cpp), so the compiler sees the same
// code: common/ energy/ outside class DP03, then the class body (pass1/ pass2/ pass3/ exact/
// diag/), then main() (cli/). Map from files to the paper: docs/ARCHITECTURE.md.
//
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
// --g-sweep direct|states (default states): how G is built for growth reads.
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

#include "common/globals.hpp"
namespace {
#include "common/logspace.hpp"
#include "energy/parameters.hpp"
#include "energy/loop_energy.hpp"
#include "common/tables.hpp"
#include "pass3/arm_gate.hpp"
#include "common/beam_tables.hpp"
#include "pass3/stem_est_tables.hpp"
#include "common/memlog.hpp"
// ------------------------------------------------------------------- DP03
class DP03 {
#include "pass3/scan_driver.hpp"
#include "pass1/prepass_products.hpp"
#include "pass2/g_states.hpp"
#include "pass3/rank_options.hpp"
#include "common/dp03_state.hpp"
#include "pass3/rank_stem_variants.hpp"
#include "pass2/growth_lazy.hpp"
#include "pass2/v_build.hpp"
#include "pass2/reads_v.hpp"
#include "diag/sweep_checks.hpp"
#include "pass2/g_direct.hpp"
#include "diag/g_eq.hpp"
#include "pass2/v_static.hpp"
#include "diag/v_diag.hpp"
#include "pass2/arms.hpp"
#include "pass2/reads_k.hpp"
#include "pass2/g_accum.hpp"
#include "pass2/sweep.hpp"
#include "pass2/reads_eval.hpp"
#include "pass3/rank_score4.hpp"
#include "pass3/rank_prune.hpp"
#include "pass3/scan_inside.hpp"
#include "pass3/scan_scores.hpp"
#include "pass3/scan_helpers.hpp"
#include "diag/scan_trace.hpp"
#include "pass3/scan_prune.hpp"
#include "pass3/scan_steps_outside.hpp"
#include "diag/memory.hpp"
#include "exact/pk_tables.hpp"
#include "pass1/inside.hpp"
#include "pass2/qz.hpp"
#include "pass3/outside.hpp"
};
#include "diag/dumps.hpp"
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
#include "diag/stem_diag.hpp"
#include "common/timers.hpp"
int main(int argc, char** argv) {
#include "cli/options.hpp"
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
      args_s.push_back(acc ? "256" : "16");
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

#include "cli/parse_args.hpp"
#include "cli/defaults_validate.hpp"
  std::string line, id;
  while (std::getline(std::cin, line)) {
#include "cli/run_setup.hpp"
#include "cli/run_passes.hpp"
#include "cli/run_report.hpp"
  }
  if (diag_f) std::fclose(diag_f);
  if (stem_diag_f) std::fclose(stem_diag_f);
  g_beam_stats.close();
  return 0;
}
