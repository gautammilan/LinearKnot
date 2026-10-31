// Pass 2 / scan, G through the multiloop states pulled on demand (--g-sweep lazy). With Chat
// the exterior prefix of the reading DP (the scan's C), every node s keeps one weight per
// right end q':
//   H_s[q'] = sum_p' Chat(p'-1) G_s[p', q'] = stop'_s[q'] + sum_{e: t -> s} w_e H_t[q'],
//   stop'_{x.y}[y] = Chat(x-1) at(x, y)   (end_ok),
// over the hypergraph sw_hyper_push records (StemEst::lz*), so the edges are those of
// --g-sweep states. A node is built when the first root that reaches it is read, memoized,
// and cut to its top K entries by H_s[q'] u_L(q') (--bG-node-key sbar: H_s[q'] Sbar_pre(q'+1)).
// Uncut with Chat = C_pre, sum_q' H_{x.y}[q'] Sbar_pre(q'+1) = Qb_b(x, y) (lz_out_id).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  struct LzMemo {
    std::vector<std::vector<std::pair<int, double>>> h;  // per node, q' ascending
    std::vector<char> state;                              // 0 new, 1 on the stack, 2 built
    long long built = 0, pulled = 0, ents = 0, late = 0;
    double sec = 0.0;
    void reset(size_t nodes) {
      h.assign(nodes, {});
      state.assign(nodes, 0);
      built = pulled = ents = late = 0;
      sec = 0.0;
    }
  };
  LzMemo lzm;
  std::vector<double> lzAcc;
  std::vector<int> lzStamp, lzTouch;
  int lzGen = 0;
  struct LzFrame {
    int u;
    size_t k;
  };
  std::vector<LzFrame> lzStack;
  std::vector<std::pair<double, int>> lzRank;

  // Chat(k) = Ch[k] must be final for k < cfin; a stop read at or past cfin is counted in late.
  // Parent chains along helices are O(n) deep: explicit stack, post-order.
  const std::vector<std::pair<int, double>>& lz_get(const StemEst& st, LzMemo& m, int id,
                                                    const std::vector<double>& Ch, int K, int key,
                                                    int cfin) {
    if (m.state[static_cast<size_t>(id)] == 2) return m.h[static_cast<size_t>(id)];
    const auto t0 = std::chrono::steady_clock::now();
    if (static_cast<int>(lzAcc.size()) < n) {
      lzAcc.assign(static_cast<size_t>(n), NEG_INF);
      lzStamp.assign(static_cast<size_t>(n), 0);
      lzGen = 0;
    }
    lzStack.clear();
    lzStack.push_back({id, 0});
    m.state[static_cast<size_t>(id)] = 1;
    while (!lzStack.empty()) {
      const int u = lzStack.back().u;
      const size_t b = st.lzOff[static_cast<size_t>(u)], e = st.lzOff[static_cast<size_t>(u) + 1];
      bool down = false;
      while (b + lzStack.back().k < e) {
        const int p = st.lzPar[b + lzStack.back().k++];
        char& sp = m.state[static_cast<size_t>(p)];
        if (sp == 2) continue;
        if (sp == 1) throw std::runtime_error("--g-sweep lazy: cycle in the G hypergraph");
        sp = 1;
        lzStack.push_back({p, 0});
        down = true;
        break;
      }
      if (down) continue;
      if (++lzGen == 0) {
        std::fill(lzStamp.begin(), lzStamp.end(), 0);
        lzGen = 1;
      }
      lzTouch.clear();
      auto add = [&](int q, double v) {
        if (lzStamp[static_cast<size_t>(q)] != lzGen) {
          lzStamp[static_cast<size_t>(q)] = lzGen;
          lzAcc[static_cast<size_t>(q)] = v;
          lzTouch.push_back(q);
        } else {
          lzAcc[static_cast<size_t>(q)] = log_add(lzAcc[static_cast<size_t>(q)], v);
        }
      };
      for (size_t r = b; r < e; ++r) {
        const double w = st.lzW[r];
        for (const auto& [q, v] : m.h[static_cast<size_t>(st.lzPar[r])]) add(q, w + v);
      }
      const int x = st.lzStopX[static_cast<size_t>(u)];
      if (x >= 0) {
        if (x > 0 && x - 1 >= cfin) ++m.late;
        const double c = x > 0 ? Ch[static_cast<size_t>(x - 1)] : 0.0;
        if (c > HALF_NEG) add(st.lzStopY[static_cast<size_t>(u)], c + st.lzStopW[static_cast<size_t>(u)]);
      }
      std::vector<std::pair<int, double>>& out = m.h[static_cast<size_t>(u)];
      out.clear();
      if (K > 0 && lzTouch.size() > static_cast<size_t>(K)) {
        lzRank.clear();
        for (int q : lzTouch) {
          const double v = lzAcc[static_cast<size_t>(q)];
          double s = key == 1 ? (static_cast<size_t>(q) < ntrUL.size() ? ntrUL[static_cast<size_t>(q)] : NEG_INF)
                              : st.prS[static_cast<size_t>(q) + 1];
          s = s > HALF_NEG ? v + s : NEG_INF;
          lzRank.emplace_back(s, q);
        }
        std::nth_element(lzRank.begin(), lzRank.begin() + K, lzRank.end(),
                         [](const std::pair<double, int>& a, const std::pair<double, int>& c) {
                           if (a.first != c.first) return a.first > c.first;
                           return a.second < c.second;
                         });
        out.reserve(static_cast<size_t>(K));
        for (int t = 0; t < K; ++t) out.emplace_back(lzRank[static_cast<size_t>(t)].second,
                                                     lzAcc[static_cast<size_t>(lzRank[static_cast<size_t>(t)].second)]);
      } else {
        out.reserve(lzTouch.size());
        for (int q : lzTouch) out.emplace_back(q, lzAcc[static_cast<size_t>(q)]);
      }
      std::sort(out.begin(), out.end());
      ++m.built;
      m.pulled += static_cast<long long>(e - b);
      m.ents += static_cast<long long>(out.size());
      m.state[static_cast<size_t>(u)] = 2;
      lzStack.pop_back();
    }
    m.sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return m.h[static_cast<size_t>(id)];
  }

  // --out-id under --g-sweep lazy (on the nested pre-pass): uncut, Chat = C_pre,
  // sum_q' H_{x.y}[q'] Sbar_pre(q'+1) against Qb_b(x, y) for every kept pair with end_ok.
  void lz_out_id(const StemEst& st) {
    LzMemo m;
    m.reset(st.lzStopX.size());
    auto rel = [](double a, double b) {
      if (a < HALF_NEG && b < HALF_NEG) return 0.0;
      if (a < HALF_NEG || b < HALF_NEG) return 1e300;
      return std::fabs(std::expm1(a - b));
    };
    std::vector<double> err;
    for (int j = 0; j < n; ++j) {
      for (size_t r = 0; r < Qb_at[j].slots(); ++r) {
        const int i = Qb_at[j].keys[r];
        if (i < 0 || !E.end_ok(i, j)) continue;
        const double ref = Qb_b[ix(i, j)];
        double got = NEG_INF;
        auto it = st.lzPair.find(key2(i, j));
        if (it != st.lzPair.end()) {
          for (const auto& [q, v] : lz_get(st, m, it->second, C, 0, 0, n + 1)) {
            const double sb = st.prS[static_cast<size_t>(q) + 1];
            if (sb > HALF_NEG) got = log_add(got, v + sb);
          }
        }
        if (got < HALF_NEG && ref < HALF_NEG) continue;
        err.push_back(rel(got, ref));
      }
    }
    double mx = 0, med = 0;
    if (!err.empty()) {
      mx = *std::max_element(err.begin(), err.end());
      std::nth_element(err.begin(), err.begin() + static_cast<long>(err.size() / 2), err.end());
      med = err[err.size() / 2];
    }
    std::fprintf(stderr, "[out-id] G-lazy pairs %zu max_rel %.3e med_rel %.3e built %lld pulled %lld\n",
                 err.size(), mx, med, m.built, m.pulled);
  }
