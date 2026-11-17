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
//   dp03_cpp [--beam B] [--loop-cap L] [--max-loop M] [--rt RT] [--wc-ends] [--cutoff C]
//            [--no-bpp] [--gap-beam group|step]
//            [--material dp03|rna1999|vienna2004 [--param-file PATH] [--dangles none|some|all|d2]]
//            < sequences
//
// --material: nested-loop energies; dp03 = the compact table below (default), rna1999 =
//   NUPACK 3's Turner 1999 file (default path ../dp03/params/rna1999.dG next to this
//   binary) with terminal AU/GU penalties; --dangles some|all adds NUPACK's dangles
//   (ExplDangle: d5 + d3 per unpaired run, "some" = the smaller one for a single base);
//   vienna2004 = LinearPartition-V's Turner 2004 model (../dp03/params/vienna2004.par),
//   --loop-cap 70 by default, --dangles none|d2. Pseudoknot penalties are the dp03 set
//   (PkParams below); the nested loops use the Turner 2004 tables.
// --beam inf (or 0) = infinite beam: nothing is pruned.
// --gap-beam group|step: 4-index tables keep `beam` states per (v, j) group (default)
//   or per step j over all groups.
//   Build with -DDP03_DEBUG (make debug) to assert kappa >= 1 and G >= F(gap).
// --b2 N: beam of the scan's two-index columns Qb, Qm, Qm2, Qp, Qmp, and of Q^z / Q^z_pre
//   Unset: 200.
// --b4 N: four-index scan beam (was --beam). Unset: 1000 under --pk-l2r on defaults, else 0.
// --bV N: V left-end cap per map. Sets both p2_v_cap and pr_vcap.
// --bG N: growth node trunc K. Unset: 16.
// --g-sweep direct|states|lazy (default lazy; states if unset and lazy cannot run, i.e. no
//   --pk-l2r on, --sweep off or --g-trace): how G is built for growth reads. states = paper-v1.
//   direct = sw_fill (pair→pair via PrMbar); states = sw_hyper_push do_g on Qm/Qm2/pair (E1–E9),
//   cut every state to bG entries by C_pre*G*S_pre; scan still applies Chat at read.
//   lazy = the same edges recorded in pass 2, pulled from the scan with Chat summed in at every
//   node (pass2/g_lazy.hpp); --bG-node K (default 32), --bG-read K_read (default K),
//   --bG-node-key ul|sbar (default ul).
// --b1 N: the nested pre-pass and the stem-estimate run at every length.
//   Unset: 200, independent of --b4. 0 keeps every state.
// Removed (error if passed): --beam.
// --stem-est-qg-r mirror|partner (default mirror; partner needs role all): Qg's stem-2 score
//   from the mirrored index H_L, or from the partner tables: log Bg(d, u, v) + alpha +
//   log S[q+1] (DP03::stem_Bg). Only the ranking changes; the removal of infeasible states
//   does not.
// --stem-est-qg-ext off|stack|bulge1 (default off; needs role all and qg-r partner): Qg's
//   score also sums sigma_ext, the stem-1 and stem-2 scores of the state grown outward by
//   case-(b) wraps (DP03::stem_ext_score); at most 30 pairs. "[stem-ext]" stderr line.
//   --stem-est-qg-ext loop: growth steps with up to --stem-est-grow-loop Lg (default 3)
//   unpaired bases per side (case-(b) wraps); per layer, nodes with log G < log(1e-8) + the
//   layer max are dropped, and at most 2000 nodes are kept per outer pair.
//   --stem-est-helpers off|grow|nested (default off; needs
//   --stem-est-qg-ext): Qgls / Qgrs are ranked by their closings into a Qg (cases f / g / h;
//   f over j <= q + 1 + 30) times, for grow, the stop-or-grow
//   outside O of the closed stem (DP03::helper_sigma) or, for nested, the outside weight of
//   the closing pair in the nested-only run (DP03::helper_sigma_nested; one lookup per state).
// --stem-est-principled off|on (default off; needs --stem-est rank, --stem-est-role all):
//   every pseudoknot table (Qg, Qgl, Qgr, Qgls, Qgrs, Qp) ranked by alpha + log Out_hat,
//   Out_hat following DP03's edges outward from the state: exact scan tables at positions
//   <= q, a nested-only pre-pass beyond q (exact for n <= 1000, else the nested beam,
//   every pair span), outside the whole pseudoknot Ctx(i, j) = e^beta1 C[i-1]
//   Sbar[j+1] + e^beta1m Bml[i][j] (Bml: the pseudoknot as one multiloop branch, from the
//   nested outside; a pseudoknot inside another one's gap is not represented).  Stems grow
//   and future stems are built by DP03's stem steps over the nested pairs with probability
//   >= 1e-7: interior steps (the chain-loop cap), multiloop steps, multiloop parts of span
//   <= --pr-ml-span S (default 150), for the growth and the future stems.
//   Qgls closings up to 30 unpaired (the engine's cap for seeds), profiles cut at --pr-cut c
//   (default 1e-8) times their max, 64 arms per gap by their gap-free future factor.
//   Exact removal is unchanged.
//   --pr-grow-select floor|bound (default floor): bound grows over every pair (pair floor 0)
//   and keeps the endpoints p'.q' of every growth node, root and key profile by
//   iota = w (U_L + U_R), U the gap-free bound of the step on the endpoint's L and R parts
//   relative to the gap's Qz[u+1][v-1] (DP03::pr_U), largest first until 1 - --pr-grow-eps
//   (default 1e-6) of the total, at most 512; nodes when built, roots and keys at the step
//   that evaluates them; replaces --pr-cut on the growth profiles.
// --pr-abl LIST (comma list; ablation of sigma_pr's couplings): bml Ctx without its multiloop term; q the partner of Stop_L beyond the state's q,
//   T(q; u, v) / Qz[q+1][q'] for T(q'; u, v); p Stop_R with K(p; u, v) of the state's p for
//   K(p'; u, v); xs every arm of the gap with Qz[u+1][v-1] / Qz[xs][x'] for its two fillers;
//   ks the same bound for K's fillers, Qz[u+1][v-1] / Qz[e][f]; sep (implies bml, q) the
//   separable evaluation at DP03::pr_eval_sep, Stop_R over all left ends p'; sepv (implies bml) the same with Stop_L over the right ends q' of the
//   profile (exact under bml); xb the xs bound on the selected arms only; qr Qz[0][q'] / Qz[0][q] for
//   q's Qz[q+1][q'].
// stdin: one sequence per line; an optional ">id" line names the next one.
// stdout per sequence: [>id] / logZ / sec / pairs K / K lines "i j p" (1-based).

#include "common/globals.hpp"
// Third-party header, so at global scope: fast_map.hpp opens ankerl::unordered_dense.
#include "common/fast_map.hpp"
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
#include "pass2/g_lazy.hpp"
#include "pass3/rank_options.hpp"
#include "common/dp03_state.hpp"
#include "pass3/rank_stem_variants.hpp"
#include "pass2/growth_lazy.hpp"
#include "pass2/v_build.hpp"
#include "pass2/reads_v.hpp"
#include "pass2/g_direct.hpp"
#include "pass2/v_static.hpp"
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
               "usage: dp03 [--b1 N] [--b2 N] [--b4 N] [--bV N] [--bG N] < sequences\n"
               "  --b1  beam of the nested pre-pass (pass I)          default 200\n"
               "  --b2  beam of the two-index states (pass III)       default 1000\n"
               "  --b4  beam of the four-index pseudoknot states (pass III)  default 3000\n"
               "  --bV  beam of V (pass II)                           default 256\n"
               "  --bG  beam of G (pass II)                           default 32\n"
               "Input: one sequence per line, optional \">name\" lines. Output: logZ and pair "
               "probabilities on stdout; diagnostics on stderr.\n");
}

// Paper setting: what --preset accurate expands to, with --bG 32 (lazy G node beam; the
// public default) instead of the historical states-mode --bG 256, plus the energy model and
// lazy G key. Always prepended before user arguments so a bare run matches the paper.
void append_paper_defaults(std::vector<std::string>& args_s) {
  args_s.push_back("--b1");
  args_s.push_back("200");
  args_s.push_back("--b2");
  args_s.push_back("1000");
  args_s.push_back("--b4");
  args_s.push_back("3000");
  args_s.push_back("--bV");
  args_s.push_back("256");
  args_s.push_back("--bG");
  args_s.push_back("32");
  for (const char* x :
       {"--gap-beam", "step", "--loop-cap", "30", "--max-loop", "30", "--wc-ends", "--rt", "0.6163207755",
        "--cutoff", "1e-5", "--stem-est", "rank", "--stem-est-role", "all", "--stem-est-loop", "0",
        "--stem-est-qg-r", "partner", "--stem-est-qg-ext", "loop", "--stem-est-grow-loop", "2",
        "--stem-est-helpers", "grow",
        "--stem-est-principled", "on", "--pr-partner-floor", "1e-20", "--pr-grow-select", "bound",
        "--pr-grow-eps", "1e-6", "--pr-cut", "0", "--pr-armq-sep", "on", "--pr-fast", "ml", "--pr-abl", "sepv",
        "--grow-ml-c", "16", "--grow-trunc", "16", "--grow-trunc-mode", "side",
        "--pk-l2r", "on", "--l2r-spawn", "inc",
        "--l2r-p3-hex", "on",
        "--l2r-p3-lin", "on", "--pr-ml-span", "100000", "--grow-ml-ub", "on", "--pr-v-key", "xs",
        "--pr-v-eps", "0", "--pr-v-trim", "suffix", "--pr-v-stream", "on",
        "--pr-v-mlk", "0", "--grow-node-trunc-r", "-1",
        "--grow-node-trunc-z", "bar", "--grow-r-from-v", "on", "--malloc-trim", "on",
        "--sweep", "on", "--v-sweep", "on", "--grammar", "full", "--ml-lead-cap", "30",
        "--material", "rna1999", "--dangles", "some", "--g-sweep", "lazy", "--bG-node-key", "ul"})
    args_s.push_back(x);
}

}  // namespace
#include "common/timers.hpp"
int main(int argc, char** argv) {
#include "cli/options.hpp"
  // Always prepend the paper setting (accurate expansion + rna1999/some + lazy G). User
  // arguments follow and override. --preset default|accurate remains parseable for old
  // scripts (options after it override; before are overridden). Paper beams are
  // --b1/--b2/--b4/--bV/--bG; pbeam is derived as (b2,b4,b4).

  for (int a = 1; a < argc; ++a) {
    if (std::string(argv[a]) == "-h" || std::string(argv[a]) == "--help") {
      usage();
      return 0;
    }
  }

  std::vector<std::string> args_s;
  args_s.push_back(argv[0] ? argv[0] : "dp03");
  const std::size_t defaults_begin = args_s.size();
  append_paper_defaults(args_s);
  const std::size_t defaults_end = args_s.size();
  for (int a = 1; a < argc; ++a) {
    if (std::string(argv[a]) == "--preset") {
      const std::string p = a + 1 < argc ? argv[++a] : "";
      if (p != "default" && p != "accurate") {
        std::fprintf(stderr, "--preset must be 'default' or 'accurate'\n");
        return 2;
      }
      const bool acc = (p == "accurate");
      // Historical preset beams (accurate keeps --bG 256 for --g-sweep states scripts).
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
            "--stem-est-helpers", "grow",
            "--stem-est-principled", "on", "--pr-partner-floor", "1e-20", "--pr-grow-select", "bound",
            "--pr-grow-eps", "1e-6", "--pr-cut", "0", "--pr-armq-sep", "on", "--pr-fast", "ml", "--pr-abl", "sepv",
            "--grow-ml-c", "16", "--grow-trunc", "16", "--grow-trunc-mode", "side",
            "--pk-l2r", "on", "--l2r-spawn", "inc",
            "--l2r-p3-hex", "on",
            "--l2r-p3-lin", "on", "--pr-ml-span", "100000", "--grow-ml-ub", "on", "--pr-v-key", "xs",
            "--pr-v-eps", "0", "--pr-v-trim", "suffix", "--pr-v-stream", "on",
            "--pr-v-mlk", "0", "--grow-node-trunc-r", "-1",
            "--grow-node-trunc-z", "bar", "--grow-r-from-v", "on", "--malloc-trim", "on",
            "--sweep", "on", "--v-sweep", "on", "--grammar", "full", "--ml-lead-cap", "30"})
        args_s.push_back(x);
    } else {
      args_s.push_back(argv[a]);
    }
  }
  {
    std::string line = "[defaults]";
    for (std::size_t k = defaults_begin; k < defaults_end; ++k) line += " " + args_s[k];
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
  return 0;
}
