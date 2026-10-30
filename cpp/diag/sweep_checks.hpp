// Diagnostics of the pass-2 sweep (--sweep-check): sw_verify (outside sweep of the pass-0
// future-stem profile) and sw_growth_check (Gamma vs the lazy growth).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

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

