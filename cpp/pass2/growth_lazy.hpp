// Old lazy growth path (--sweep off): pr_qz, distinct-key accumulators (dacc), profile
// finishing (pr_finish_*), pr_steps (stem steps x.y -> x'.y'), split / helper mass and
// pr_node_acc; partner caps.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

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
    return qz_log(a, b);
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
  int pr_abl = 0, pr_abl_p = -1, pr_abl_q = -1;
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
    if (!ml || (scan && rec_noml) || !E.end_ok(x, y)) return;
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
    if ((scan && rec_noml) || !E.end_ok(x, y)) return;
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
    size_t n_int = 0;
    if (pr_fast >= 2) {
      for (const ChainNbr& t : nb)
        if (t.t == 0) terms.emplace_back(child(t.x, t.y), t.w);
      n_int = terms.size();
      MlSplit s;
      pr_split(x, y, scan, s);
      if (s.on) {
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
    double ref = self_w;
    for (const auto& [q, w] : terms)
      if (q->mx > HALF_NEG) ref = std::max(ref, w + q->mx);
    dacc_begin(ref > HALF_NEG ? ref : 0.0);
    if (self_w > HALF_NEG) dacc_add1(!scan && pr_vxs ? key2(x, 0) : key2(x, y), self_w);
    for (size_t i = 0; i < n_int; ++i) dacc_add(*terms[i].first, terms[i].second);
    for (size_t i = n_int; i < terms.size(); ++i) dacc_add(*terms[i].first, terms[i].second);
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

 public:

 private:
  static int kcl_xs(uint64_t k) { return static_cast<int>(k >> 40); }
  static int kcl_ys(uint64_t k) { return static_cast<int>((k >> 8) & 0xffffffffu); }

