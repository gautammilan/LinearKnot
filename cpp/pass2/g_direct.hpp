// Pass 2, direct G (--g-sweep direct, not the paper default): sw_fill builds Gamma pair -> pair
// through PrMbar with right-gap limit rf.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

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
    // Cap a pair's map before it is pushed: the top b_G q' by Ghat + ntrUL, all (p', q') kept for those q'.
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
      // Rank by q' with Ghat[q'] = sum_p' G*C_pre (+ ntrUL).
      std::unordered_map<int, double> byq;
      byq.reserve(m.size());
      for (const auto& [k, w] : m) {
        if (w < HALF_NEG) continue;
        const int p = k_hi(k), q = k_lo(k);
        double g = w;
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
        if (q >= 0 && q < static_cast<int>(ntrUL.size()) && ntrUL[q] > HALF_NEG) {
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
      // Retain all raw (p',q') mass for the kept q'.
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

