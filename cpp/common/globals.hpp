// Global switches set by the CLI (pass-2 G build and cut: g_sweep, g_ntr_K = b_G;
// diagnostic traces), the standard-library includes, and the index header stem_est.hpp.
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
// g_L / g_R come from layered marginal vectors (DP03::RN).
// --grow-ml-L / -R: highest multiloop layer m used by the L / R stop (0 none, 1 depth1, 2 all);
// --grow-ml-c C: top C choices per multiloop half-step (0: all).
static int g_ml_L = 2, g_ml_R = 2, g_ml_c = 0;
// --grow-ml-ub on|off (default on): the --grow-ml-c cut of a node's multiloop terms builds only the
// helpers whose tot bound can reach the top C (same terms kept, same order)
static bool g_ml_ub = true;
// g_skip_bml: build_principled leaves Bml empty.
static bool g_skip_bml = false;
// --grow-trunc K: at root / key evaluation, the g_L entries ranked by
// g_L[q'] u_L(q') and the g_R entries by g_R[p'] u_R(p') (step-q pr_U halves without Bml),
// largest first until 1 - g_tr_eps of the total, at most K; --grow-trunc-mode side:
// per side (K per side), pool: one ranking over both sides (K in total).
static int g_tr_K = 0, g_tr_pool = 0;
static double g_tr_eps = 1e-6;
// --bG K (ml all/all): every growth node / root / helper / key
// keeps one layer-summed list per side, its top K entries by --grow-trunc's ranking at build
// time, built from its terms' lists; the full layer vectors are never stored.
static int g_ntr_K = 0;
enum { G_SWEEP_DIRECT = 0, G_SWEEP_STATES = 1, G_SWEEP_LAZY = 2 };
static int g_sweep = G_SWEEP_LAZY;
// --g-sweep lazy: --bG-node K keeps K entries q' per node (0: all), ranked by
// H[q'] u_L(q') (--bG-node-key ul, default) or H[q'] Sbar_pre(q'+1) (sbar); --bG-read K_read
// sets b_G (g_ntr_K) for the read, default K. Ther.ther needs ul and K >= 24.
static int g_lz_K = 32, g_lz_Kread = -1, g_lz_key = 1;
static const char* g_sweep_name() {
  return g_sweep == G_SWEEP_DIRECT ? "direct" : g_sweep == G_SWEEP_LAZY ? "lazy" : "states";
}
// --out-id: compare reconstructed outsides to pass-1 Qb_b. Stop term uses E.at,
// the exterior-loop pair factor pass 1 adds into Qb_b (this check only; production seed is the same E.at).
static bool g_out_id = false;
// --g-trace P,Q[,X,Y] (--g-sweep states, 0-based): log each state whose G cut drops entry
// (P,Q); with X,Y only states (i,j) with P <= i <= X, Y <= j <= Q (between the two pairs).
static int g_trace_p = -1, g_trace_q = -1, g_trace_x = -1, g_trace_y = -1;
// --grow-node-trunc-r KR: K of the R side (0: K; -1: R not truncated, stored dense in R[0], the
// layers summed)
static int g_ntr_KR = 0;
// --grow-r-from-v on (with --bG K --grow-node-trunc-r -1, --pr-v-key xs): the
// U-type scores read V (kept after the build; seed roots through the arm entries) in place of
// the R side of the growth vectors, which is then not built (only its total, the --grow-ml-c rank).
static bool g_rfv = false;
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

#include "stem_est.hpp"

