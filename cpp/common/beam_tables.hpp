// Beam bookkeeping shared by all passes: table ids BT_* (Qb, Qm, Qm2, Qmp, PKM, Qp, Qg, Qgls,
// Qgl, Qgr, ...), state_key, per-table beam statistics (--beam-stats) and the drop log
// (--dump-pruned).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

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

