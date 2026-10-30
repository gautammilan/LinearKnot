// StemEst: the pass-1 products the pass-3 ranking reads (S, C0, stem index, chain and
// principled data: prS, prQz, PrMbar, ...). Built by pass1/prepass_products.hpp.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

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

