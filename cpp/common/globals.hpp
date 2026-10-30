// Global switches set by the CLI (pass-2 G build and cut: g_sweep, g_cut, g_ntr_K = b_G;
// diagnostic traces and dumps), the standard-library includes, and the two index headers
// (gate_index.hpp, stem_est.hpp).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

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
static int g_sweep = G_SWEEP_STATES;
// --out-id: compare reconstructed outsides to pass-1 Qb_b. Stop term uses E.at,
// the exterior-loop pair factor pass 1 adds into Qb_b (this check only; production seed is the same E.at).
static bool g_out_id = false;
// --g-cut-at all|pair (--g-sweep states): cut G at every state, or at pair states only.
static bool g_cut_pair_only = false;
// --g-trace P,Q[,X,Y] (--g-sweep states, 0-based): log each state whose G cut drops entry
// (P,Q); with X,Y only states (i,j) with P <= i <= X, Y <= j <= Q (between the two pairs).
static int g_trace_p = -1, g_trace_q = -1, g_trace_x = -1, g_trace_y = -1;
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

