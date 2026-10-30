// Stem-estimate ranking variants (--stem-est-qg-ext, -grow-loop, -helpers, -chain, -hr):
// stem_B, stem_Bg, growth nodes, helper sigma, chains. Score only; never changes values.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  // ------------------------------------------------------- stem estimate
  // log B(d, e, e2) = log sum over f <= e2 - 1 of A(d, e, f) Qz[f+1][e2-1].
  double stem_B(int d, int e, int e2) const {
    auto it = stA.find(key2(d, e));
    double acc = NEG_INF;
    if (it == stA.end()) return acc;
    const StemA& a = it->second;
    for (size_t k = 0; k < a.fs.size(); ++k) {
      if (a.fs[k] > e2 - 1) break;
      acc = log_add(acc, a.vals[k] + qz_log(a.fs[k] + 1, e2 - 1));
    }
    return acc;
  }

  static uint64_t key3(int a, int b, int c) {
    return (static_cast<uint64_t>(a) << 42) | (static_cast<uint64_t>(b) << 21) |
           static_cast<uint64_t>(c);
  }

  void stem_cache_step(int q) {
    if (q == stc_q) return;
    stc_q = q;
    std::unique_ptr<SecTimer> ptm(lpf ? new SecTimer(lpf_stc_sec) : nullptr);
    if (lpf) {
      ++lpf_stc_n;
      lpf_stc_buckets += static_cast<long long>(stBc.bucket_count() + stBgc.bucket_count() + stX.bucket_count() +
                                                stHc.bucket_count() + stOc.bucket_count());
      if (pr_on)
        lpf_stc_scan += static_cast<long long>(prTc.size() + prAQ.size() + prKls.bucket_count() +
                                               prKrs.bucket_count() + prRKls.bucket_count() +
                                               prRKrs.bucket_count() + prRT.bucket_count() + prAc.bucket_count());
    }
    stBc.clear();
    stBgc.clear();
    stX.clear();
    stHc.clear();
    stOc.clear();
    if (stem_profile) {
      stP.clear();
      stPg.clear();
      stPls.clear();
      stPrs.clear();
      stM.clear();
    }
    if (pr_on) pr_step(q);
  }

  double stem_B_cached(int d, int e, int v) {
    const uint64_t key = key3(d, e, v);
    auto it = stBc.find(key);
    if (it != stBc.end()) return it->second;
    double b = stem_B(d, e, v);
    stBc.emplace(key, b);
    return b;
  }

  // Partner weight of a Qg state (d, u, v) as the right stem (DP03._stem_Bg): the engine's
  // Qg -> Qgl edge (beta2 + at(u, v), end_ok(u, v), filler Qz[u+1][e-1]) times
  // B(d, e, v) (Qgl -> Qgr filler, Qp with the surviving left stems), over the e of stE[d]
  // in [u+1, v-1], ascending.
  double stem_Bg_raw(int d, int u, int v) {
    if (!E.end_ok(u, v)) return NEG_INF;
    auto it = stE.find(d);
    if (it == stE.end()) return NEG_INF;
    const std::vector<int>& es = it->second;
    double acc = NEG_INF;
    for (auto k = std::lower_bound(es.begin(), es.end(), u + 1); k != es.end() && *k <= v - 1; ++k) {
      const int e = *k;
      ++bg_visits;
      double b = stem_B_cached(d, e, v);
      if (b < HALF_NEG) continue;
      double qz = qz_log(u + 1, e - 1);
      if (qz < HALF_NEG) continue;
      acc = log_add(acc, qz + b);
    }
    if (acc < HALF_NEG) return NEG_INF;
    return acc + (E.w_stem + E.at(u, v));
  }

  double stem_Bg(int d, int u, int v) {
    ++bg_lookups;
    const uint64_t key = key3(d, u, v);
    auto it = stBgc.find(key);
    if (it != stBgc.end()) return it->second;
    ++bg_misses;
    double bg = stem_Bg_raw(d, u, v);
    stBgc.emplace(key, bg);
    return bg;
  }

  // L-stem factor of a final right end q: log H(u, v, q) - log C0[q]; with hr_av the index
  // is st->index_av (av or phi entries), and with hr_sq (av) the constant is log S[q+1]
  // (h + hr_c(q) equals h - C0[q] bit for bit otherwise).
  bool hr_av = false, hr_sq = false, hr_check = false;
  double hr_h(int u, int v, int q) const { return hr_av ? stem_->Hav(u, v, q) : stem_->H(u, v, q); }
  double hr_c(int q) const { return hr_sq ? stem_->S[q + 1] : -stem_->C0[q]; }

  // --stem-hr-check: Hav(u, v, q) by listing every partner stem (innermost x'.y', x' in
  // [u+1, v-1], y' > q, stacked outward while the pairs can form) and summing
  // E_out E_in prod w_stack prod phi_gamma in long double.
  double hr_naive(int u, int v, int q) const {
    if (stem_->hr_phi) return hr_naive_phi(u, v, q);
    std::vector<long double> terms;
    for (int x = u + 1; x <= v - 1; ++x)
      for (int y = q + 1; y < n; ++y) {
        if (y - x < MIN_PAIR || !E.pairs(x, y) || !E.end_ok(x, y)) continue;
        const double lp0 = stem_->lphi[ix(x, y)];
        if (lp0 < HALF_NEG) continue;
        long double w = static_cast<long double>(E.at(x, y)) + lp0;
        for (int i = x, j = y;;) {
          if (E.end_ok(i, j)) terms.push_back(w + E.at(i, j));
          const int i0 = i - 1, j0 = j + 1;
          if (i0 < 0 || j0 > n - 1 || j0 - i0 < MIN_PAIR || !E.pairs(i0, j0)) break;
          const double s = E.stem_interior(i0, i, j, j0), lp = stem_->lphi[ix(i0, j0)];
          if (s < HALF_NEG || lp < HALF_NEG) break;
          w += static_cast<long double>(s) + lp;
          i = i0;
          j = j0;
        }
      }
    if (terms.empty()) return NEG_INF;
    long double m = *std::max_element(terms.begin(), terms.end()), sum = 0.0L;
    for (long double t : terms) sum += std::exp(t - m);
    return static_cast<double>(m + std::log(sum));
  }

  // --stem-hr-check with phi: H_R'(u, v, q) by listing every partner stem of the original
  // table (innermost x'.y' an end pair, x' in [u+1, v-1], y' > q; outward steps with at most
  // the table's loop unpaired bases per side through pairs that can form; outermost pair an
  // end pair) and summing E_in C0[y'-1] prod phi_gamma prod w E_out S[j_1+1] in long double.
  double hr_naive_phi(int u, int v, int q) const {
    const int lc = std::min(stem_->loop, cap);
    auto partner = [&](int i, int j) {
      if (i < 0 || j > n - 1) return false;
      const auto& pr = partners_right[i];
      return std::find(pr.begin(), pr.end(), j) != pr.end();
    };
    std::vector<long double> terms;
    std::vector<std::tuple<int, int, long double>> stack;
    for (int x = u + 1; x <= v - 1; ++x)
      for (int y : partners_right[x]) {
        if (y <= q || !E.end_ok(x, y)) continue;
        const double lp = stem_->lphi[ix(x, y)];
        if (lp < HALF_NEG) continue;
        stack.emplace_back(x, y, static_cast<long double>(E.at(x, y)) + stem_->C0[y - 1] + lp);
        while (!stack.empty()) {
          auto [i, j, w] = stack.back();
          stack.pop_back();
          if (E.end_ok(i, j)) terms.push_back(w + E.at(i, j) + stem_->S[j + 1]);
          for (int i0 = std::max(0, i - 1 - lc); i0 < i; ++i0)
            for (int j0 = j + 1; j0 <= std::min(n - 1, j + 1 + lc); ++j0) {
              if (!partner(i0, j0)) continue;
              const double s = E.stem_interior(i0, i, j, j0), lp0 = stem_->lphi[ix(i0, j0)];
              if (s < HALF_NEG || lp0 < HALF_NEG) continue;
              stack.emplace_back(i0, j0, w + s + lp0);
            }
        }
      }
    if (terms.empty()) return NEG_INF;
    long double m = *std::max_element(terms.begin(), terms.end()), sum = 0.0L;
    for (long double t : terms) sum += std::exp(t - m);
    return static_cast<double>(m + std::log(sum));
  }

  double stem_H_cached(int u, int v, int q2) {
    const uint64_t key = key3(u, v, q2);
    auto it = stHc.find(key);
    if (it != stHc.end()) return it->second;
    double h = hr_h(u, v, q2);
    stHc.emplace(key, h);
    return h;
  }

  // Extension nodes of a Qg state with outer pair p.stc_q (DP03._stem_ext_nodes).
  const std::vector<ExtNode>& stem_ext_nodes(int p) { return stem_grow_nodes(p, stc_q); }

  // Growth nodes of outer pair p.q (DP03._stem_grow_nodes): the outer pairs (d', q') reached
  // by 1..stem_ext_max case-(b) wraps (weight stem_interior + pair_bonus each; stack: d'-1,
  // q'+1; bulge1 also d'-2, q'+1 and d'-1, q'+2; loop: d'-1-a, q'+1+b for a, b <=
  // min(grow_lg, cap), a outer), with log W summed over the paths (layer by layer; nodes in
  // first-reached order).  loop: after each layer, nodes with log W < log(grow_cut) + the
  // layer max are dropped, then new nodes beyond grow_nodes in total are dropped (lowest W
  // first, ties by order); dropped nodes are neither kept nor expanded.  Only end_ok nodes
  // are kept: every exit of a Qg state needs end_ok of its outer pair.
  const std::vector<ExtNode>& stem_grow_nodes(int p, int q) {
    const uint64_t key = key2(p, q);
    auto it = stX.find(key);
    if (it != stX.end()) return it->second;
    std::vector<ExtNode> out;
    if (stem_ext == 3) {
      const int lg = std::min(grow_lg, cap);
      const double lcut = grow_cut > 0.0 ? std::log(grow_cut) : NEG_INF;
      std::vector<ExtNode> cur{{p, q, 0.0}}, nxt, tot;
      std::unordered_map<uint64_t, size_t> nidx, tidx;
      std::vector<size_t> fresh, order;
      std::vector<char> keep;
      for (int t = 1; t <= stem_ext_max && !cur.empty(); ++t) {
        nxt.clear();
        nidx.clear();
        for (const ExtNode& c : cur)
          for (int a = 0; a <= lg; ++a)
            for (int b = 0; b <= lg; ++b) {
              const int i = c.d - 1 - a, j = c.q + 1 + b;
              if (i < 0 || j > n - 1) continue;
              const double s = E.stem_interior(i, c.d, c.q, j);
              if (s < HALF_NEG) continue;
              const double val = c.w + (s + pair_bonus);
              auto r = nidx.emplace(key2(i, j), nxt.size());
              if (r.second) nxt.push_back({i, j, val});
              else nxt[r.first->second].w = log_add(nxt[r.first->second].w, val);
            }
        grow_layer += static_cast<long long>(nxt.size());
        keep.assign(nxt.size(), 1);
        if (grow_cut > 0.0 && !nxt.empty()) {
          double mx = NEG_INF;
          for (const ExtNode& x : nxt) mx = std::max(mx, x.w);
          const double thr = lcut + mx;
          for (size_t k = 0; k < nxt.size(); ++k)
            if (nxt[k].w < thr) {
              keep[k] = 0;
              ++grow_cut_drop;
            }
        }
        if (grow_nodes > 0) {
          fresh.clear();
          for (size_t k = 0; k < nxt.size(); ++k)
            if (keep[k] && !tidx.count(key2(nxt[k].d, nxt[k].q))) fresh.push_back(k);
          const size_t room = static_cast<size_t>(grow_nodes) > tot.size()
                                  ? static_cast<size_t>(grow_nodes) - tot.size() : 0;
          if (fresh.size() > room) {
            order.resize(fresh.size());
            for (size_t m = 0; m < fresh.size(); ++m) order[m] = m;
            std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
              return nxt[fresh[a]].w > nxt[fresh[b]].w;
            });
            for (size_t m = room; m < order.size(); ++m) keep[fresh[order[m]]] = 0;
            grow_cap_drop += static_cast<long long>(fresh.size() - room);
          }
        }
        cur.clear();
        for (size_t k = 0; k < nxt.size(); ++k) {
          if (!keep[k]) continue;
          const ExtNode& x = nxt[k];
          cur.push_back(x);
          auto r = tidx.emplace(key2(x.d, x.q), tot.size());
          if (r.second) tot.push_back(x);
          else tot[r.first->second].w = log_add(tot[r.first->second].w, x.w);
        }
      }
      grow_list_nodes += static_cast<long long>(tot.size());
      for (const ExtNode& x : tot)
        if (E.end_ok(x.d, x.q)) out.push_back(x);
    } else if (stem_ext == 1) {
      double w = 0.0;
      for (int t = 1; t <= stem_ext_max; ++t) {
        const int i = p - t, j = q + t;
        if (i < 0 || j > n - 1) break;
        const double s = E.stem_interior(i, i + 1, j - 1, j);
        if (s < HALF_NEG) break;
        w = w + (s + pair_bonus);
        ++grow_list_nodes;
        if (E.end_ok(i, j)) out.push_back({i, j, w});
      }
    } else {
      static const int MV[3][2] = {{1, 1}, {2, 1}, {1, 2}};
      std::vector<ExtNode> cur{{p, q, 0.0}}, nxt, tot;
      std::unordered_map<uint64_t, size_t> nidx, tidx;
      for (int t = 1; t <= stem_ext_max && !cur.empty(); ++t) {
        nxt.clear();
        nidx.clear();
        for (const ExtNode& c : cur)
          for (const auto& m : MV) {
            const int i = c.d - m[0], j = c.q + m[1];
            if (i < 0 || j > n - 1) continue;
            const double s = E.stem_interior(i, c.d, c.q, j);
            if (s < HALF_NEG) continue;
            const double val = c.w + (s + pair_bonus);
            auto r = nidx.emplace(key2(i, j), nxt.size());
            if (r.second) nxt.push_back({i, j, val});
            else nxt[r.first->second].w = log_add(nxt[r.first->second].w, val);
          }
        for (const ExtNode& x : nxt) {
          auto r = tidx.emplace(key2(x.d, x.q), tot.size());
          if (r.second) tot.push_back(x);
          else tot[r.first->second].w = log_add(tot[r.first->second].w, x.w);
        }
        cur.swap(nxt);
      }
      grow_list_nodes += static_cast<long long>(tot.size());
      for (const ExtNode& x : tot)
        if (E.end_ok(x.d, x.q)) out.push_back(x);
    }
    ++grow_lists;
    return stX.emplace(key, std::move(out)).first->second;
  }

  // sigma_ext of the Qg state (p, u, v, stc_q) with log weight a (DP03._stem_ext_score): log
  // sum over its extension nodes (d', q', W) of the stem-1 score four_score(d', u, v, a) + W +
  // H(u, v, q') - C0[q'] + w_pk and the stem-2 score Bg(d', u, v) + a + W + S[q'+1].
  double stem_ext_score(int p, int u, int v, double a) {
    const std::vector<ExtNode>& xs = stem_ext_nodes(p);
    const double w_pk = E.w_beta1 + 2.0 * E.w_stem;
    double acc = NEG_INF;
    for (const ExtNode& x : xs) {
      double sl = NEG_INF, sr = NEG_INF;
      const double h = stem_H_cached(u, v, x.q);
      if (h > HALF_NEG) {
        sl = four_score(x.d, u, v, a) + x.w + h + hr_c(x.q) + w_pk;
        if (stem_consistent) sl += (E.at(u, v) + E.at(x.d, x.q));
      }
      const double bg = stem_Bg(x.d, u, v);
      if (bg > HALF_NEG) {
        sr = bg + a + x.w + stem_->S[x.q + 1];
        if (stem_consistent) sr += (E.at(x.d, x.q) + E.w_beta1);
      }
      acc = log_add(acc, log_add(sl, sr));
    }
    ext_nodes += static_cast<long long>(xs.size());
    return acc;
  }

  // --stem-ext-check: sigma_ext by walking every wrap path separately (no shared path sums;
  // H and Bg recomputed), every term summed in long double.
  void ext_paths(int d, int q, int t, long double w, std::vector<std::tuple<int, int, long double>>& out) {
    if (t > 0 && E.end_ok(d, q)) out.emplace_back(d, q, w);
    if (t == stem_ext_max) return;
    std::vector<std::pair<int, int>> mv;
    if (stem_ext == 3) {
      const int lg = std::min(grow_lg, cap);
      for (int a = 0; a <= lg; ++a)
        for (int b = 0; b <= lg; ++b) mv.emplace_back(a + 1, b + 1);
    } else {
      mv = {{1, 1}, {2, 1}, {1, 2}};
      if (stem_ext == 1) mv.resize(1);
    }
    for (const auto& m : mv) {
      const int i = d - m.first, j = q + m.second;
      if (i < 0 || j > n - 1) continue;
      const double s = E.stem_interior(i, d, q, j);
      if (s < HALF_NEG) continue;
      ext_paths(i, j, t + 1, w + static_cast<long double>(s) + pair_bonus, out);
    }
  }

  double stem_ext_naive(int p, int u, int v, double a) {
    std::vector<std::tuple<int, int, long double>> paths;
    ext_paths(p, stc_q, 0, 0.0L, paths);
    const long double w_pk = static_cast<long double>(E.w_beta1) + 2.0L * E.w_stem;
    std::vector<long double> terms;
    for (const auto& [d, q, w] : paths) {
      const double h = hr_h(u, v, q);
      if (h > HALF_NEG) {
        long double s = static_cast<long double>(prefix_log(d)) + a + w + h + hr_c(q) + w_pk;
        if (score_qz) s += qz_log(u + 1, v - 1);
        if (stem_consistent) s += static_cast<long double>(E.at(u, v)) + E.at(d, q);
        terms.push_back(s);
      }
      const double bg = stem_Bg_raw(d, u, v);
      if (bg > HALF_NEG) {
        long double s = static_cast<long double>(bg) + a + w + stem_->S[q + 1];
        if (stem_consistent) s += static_cast<long double>(E.at(d, q)) + E.w_beta1;
        terms.push_back(s);
      }
    }
    if (terms.empty()) return NEG_INF;
    long double m = *std::max_element(terms.begin(), terms.end()), sum = 0.0L;
    for (long double t : terms) sum += std::exp(t - m);
    return static_cast<double>(m + std::log(sum));
  }

  // ------------------------------------------------- helper scores (Qgls, Qgrs)
  // Stop term of a stem with final outer pair d.q2 and inner pair (u, v), without alpha
  // (DP03._helper_stop): logsumexp of C[d-1] Qz[u+1][v-1] H(u, v, q2) / C0[q2] w_pk and
  // Bg(d, u, v) S[q2+1] (consistent constants as stem_ext_score).
  double helper_stop(int d, int q2, int u, int v) {
    double sl = NEG_INF, sr = NEG_INF;
    const double h = stem_H_cached(u, v, q2);
    if (h > HALF_NEG) {
      double s = prefix_log(d);
      if (score_qz) s += qz_log(u + 1, v - 1);
      sl = s + h + hr_c(q2) + (E.w_beta1 + 2.0 * E.w_stem);
      if (stem_consistent) sl += (E.at(u, v) + E.at(d, q2));
    }
    const double bg = stem_Bg(d, u, v);
    if (bg > HALF_NEG) {
      sr = bg + stem_->S[q2 + 1];
      if (stem_consistent) sr += (E.at(d, q2) + E.w_beta1);
    }
    return log_add(sl, sr);
  }

  // O(x, y; u, v) (DP03._helper_O): the stem closed at x.y stops there or grows on (growth
  // nodes of x.y with their path weights), cached per step.
  double helper_O(int x, int y, int u, int v) {
    const uint64_t key = (static_cast<uint64_t>(x) << 48) | (static_cast<uint64_t>(y) << 32) |
                         (static_cast<uint64_t>(u) << 16) | static_cast<uint64_t>(v);
    auto it = stOc.find(key);
    if (it != stOc.end()) return it->second;
    double acc = helper_stop(x, y, u, v);
    const std::vector<ExtNode>& xs = stem_grow_nodes(x, y);
    for (const ExtNode& nd : xs) acc = log_add(acc, nd.w + helper_stop(nd.d, nd.q, u, v));
    hlp_nodes += static_cast<long long>(xs.size());
    stOc.emplace(key, acc);
    return acc;
  }

  // A closing of a helper state into a Qg with outer pair x.y and closing weight w (alpha
  // not included); filled by helper_closings.
  struct Closing {
    int x, y;
    double w;
  };

  // Closings of Qgls[q] state (i, u, v): case (f) into outer pair (i-1, j), j = q+1 ..
  // min(n-1, q+1+stem_ml_span) (seed group v == q: run capped as in fill_Qg).
  void closings_ls(int i, int v, int q, std::vector<Closing>& out) const {
    out.clear();
    const int x = i - 1;
    if (x < 0) return;
    const int j_hi = std::min(n - 1, q + 1 + stem_ml_span);
    for (int j = q + 1; j <= j_hi; ++j) {
      const int unpaired = j - q - 1;
      if (v == q && unpaired > cap) continue;
      if (!E.pairs(x, j) || !E.end_ok(x, j)) continue;
      double w_ext = E.w_ap + E.w_bp + E.w_cp * unpaired + pair_bonus;
      if (dz) w_ext += E.run(q + 1, j - 1);
      out.push_back({x, j, w_ext + E.at(x, j)});
    }
  }

  // Closings of Qgrs[q] state (c, u, v) into outer pair (i, q+1): case (g) (unpaired run on
  // the left, i ascending), then case (h) (Qm[i+1][c-1] on the left, i ascending).
  void closings_rs(int c, int q, std::vector<Closing>& out) const {
    out.clear();
    const int j = q + 1;
    if (j > n - 1) return;
    for (int i = lead_i_lo(c); i < c; ++i) {
      if (!E.pairs(i, j) || !E.end_ok(i, j)) continue;
      double w = E.w_ap + E.w_bp + E.w_cp * (c - i - 1) + E.at(i, j) + pair_bonus;
      if (dz) w += E.run(i + 1, c - 1);
      out.push_back({i, j, w});
    }
    if (c - 1 < 0) return;
    const OMap& at = qmp_at(c - 1);
    std::vector<std::pair<int, double>> qm;
    qm_pair_slots(c - 1, j, [&](int r) {
      if (at.keys[r] >= 0) qm.emplace_back(at.keys[r], at.vals[r]);
    });
    std::sort(qm.begin(), qm.end());
    for (const auto& [i, val] : qm) {
      const int outer = i - 1;
      if (outer < 0 || outer >= c || !E.pairs(outer, j) || !E.end_ok(outer, j)) continue;
      out.push_back({outer, j, val + E.w_ap + E.w_bp + E.at(outer, j) + pair_bonus});
    }
  }

  // sigma of a helper state with inner pair (u, v) and log weight a: logsumexp over its
  // closings of w + a + O(x, y; u, v).
  double helper_sigma(const std::vector<Closing>& cl, int u, int v, double a) {
    double acc = NEG_INF;
    for (const Closing& c : cl) acc = log_add(acc, c.w + a + helper_O(c.x, c.y, u, v));
    return acc;
  }

  // --stem-ext-check: the helper score by walking every closing and growth path separately
  // (H and Bg recomputed), every term summed in long double.
  double helper_naive(const std::vector<Closing>& cl, int u, int v, double a) {
    std::vector<long double> terms;
    std::vector<std::tuple<int, int, long double>> paths;
    const long double w_pk = static_cast<long double>(E.w_beta1) + 2.0L * E.w_stem;
    for (const Closing& c : cl) {
      paths.clear();
      paths.emplace_back(c.x, c.y, 0.0L);
      ext_paths(c.x, c.y, 0, 0.0L, paths);
      for (const auto& [d, q, w] : paths) {
        const long double base = static_cast<long double>(c.w) + a + w;
        const double h = hr_h(u, v, q);
        if (h > HALF_NEG) {
          long double s = base + prefix_log(d) + h + hr_c(q) + w_pk;
          if (score_qz) s += qz_log(u + 1, v - 1);
          if (stem_consistent) s += static_cast<long double>(E.at(u, v)) + E.at(d, q);
          terms.push_back(s);
        }
        const double bg = stem_Bg_raw(d, u, v);
        if (bg > HALF_NEG) {
          long double s = base + bg + stem_->S[q + 1];
          if (stem_consistent) s += static_cast<long double>(E.at(d, q)) + E.w_beta1;
          terms.push_back(s);
        }
      }
    }
    if (terms.empty()) return NEG_INF;
    long double m = *std::max_element(terms.begin(), terms.end()), sum = 0.0L;
    for (long double t : terms) sum += std::exp(t - m);
    return static_cast<double>(m + std::log(sum));
  }

  // --stem-est-helpers nested: log sum over the closings (x, y) of w + beta_nested(x, y), in
  // closing order.  The closings depend only on (i, v == q) for Qgls and on c for Qgrs, so the
  // sum is cached per step; returns it and sets ncl to the number of closings.
  double helper_nested_sum(bool ls, int p, int v, int q, int& ncl) {
    if (hn_q != q) {
      hn_q = q;
      hnL.assign(2 * static_cast<size_t>(n), {NEG_INF, -1});
      hnR.assign(static_cast<size_t>(n), {NEG_INF, -1});
    }
    std::pair<double, int>& slot = ls ? hnL[2 * static_cast<size_t>(p) + (v == q)] : hnR[p];
    if (slot.second < 0) {
      std::vector<Closing> cl;
      if (ls) closings_ls(p, v, q, cl);
      else closings_rs(p, q, cl);
      double acc = NEG_INF;
      for (const Closing& c : cl) acc = log_add(acc, c.w + stem_->B(c.x, c.y));
      slot = {acc, static_cast<int>(cl.size())};
    }
    ncl = slot.second;
    return slot.first;
  }

  // sigma of a helper state (nested): alpha + its closing sum + the gap filler, as four_score.
  double helper_sigma_nested(bool ls, const Ent& e, int v, int q, int& ncl) {
    const double s = helper_nested_sum(ls, e.p, v, q, ncl);
    if (s < HALF_NEG) return NEG_INF;
    double r = e.a + s;
    if (score_qz) r += qz_log(e.u + 1, v - 1);
    return r;
  }

  // --stem-ext-check (nested): the same score without the per-step cache, in long double.
  double helper_nested_naive(bool ls, const Ent& e, int v, int q) {
    std::vector<Closing> cl;
    if (ls) closings_ls(e.p, v, q, cl);
    else closings_rs(e.p, q, cl);
    std::vector<long double> terms;
    for (const Closing& c : cl) {
      const double b = stem_->B(c.x, c.y);
      if (b < HALF_NEG) continue;
      long double t = static_cast<long double>(e.a) + c.w + b;
      if (score_qz) t += qz_log(e.u + 1, v - 1);
      terms.push_back(t);
    }
    if (terms.empty()) return NEG_INF;
    long double m = *std::max_element(terms.begin(), terms.end()), sum = 0.0L;
    for (long double t : terms) sum += std::exp(t - m);
    return static_cast<double>(m + std::log(sum));
  }

  // ------------------------------------------------ profiles (--stem-est-profile)
  // Sparse profile accumulator: entries in first-reached order, each summed in log space in the
  // order added; finish() sorts them by index (DP03._ProfAcc).
  struct ProfAcc {
    std::unordered_map<int, size_t> li, ri;
    std::vector<std::pair<int, double>> l, r;
    static void add(std::unordered_map<int, size_t>& ix, std::vector<std::pair<int, double>>& v,
                    int key, double w) {
      auto res = ix.emplace(key, v.size());
      if (res.second) v.emplace_back(key, w);
      else v[res.first->second].second = log_add(v[res.first->second].second, w);
    }
    void addL(int y, double w) { add(li, l, y, w); }
    void addR(int x, double w) { add(ri, r, x, w); }
    void finish(Prof& p) {
      auto by_ix = [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
        return a.first < b.first;
      };
      std::sort(l.begin(), l.end(), by_ix);
      std::sort(r.begin(), r.end(), by_ix);
      for (const auto& e : l) {
        p.ly.push_back(e.first);
        p.lw.push_back(e.second);
      }
      for (const auto& e : r) {
        p.rx.push_back(e.first);
        p.rw.push_back(e.second);
      }
    }
  };

  void pf_note_step() {
    if (pf_q != stc_q) {
      pf_q = stc_q;
      ++pf_steps;
    }
  }

  // Profile of the stem closed at outer pair x.y (DP03._pair_prof): the stop-here node (x, y)
  // with G = 1, then its growth nodes in list order.  Node (d', q', G) adds
  // G C[d'-1] w_pk / C0[q'] to L at q' and G S[q'+1] to R at d' (consistent: at(d', q') to
  // both and beta1 to R).  Untruncated; cached per step.
  const Prof& pair_prof(int x, int y) {
    if (stem_chain) return chain_root(x, y);
    const uint64_t key = key2(x, y);
    auto it = stP.find(key);
    if (it != stP.end()) return it->second;
    const double w_pk = E.w_beta1 + 2.0 * E.w_stem;
    ProfAcc acc;
    auto add = [&](int d, int q2, double g) {
      double cl = ((g + prefix_log(d)) + hr_c(q2)) + w_pk;
      double cr = g + stem_->S[q2 + 1];
      if (stem_consistent) {
        cl += E.at(d, q2);
        cr += (E.at(d, q2) + E.w_beta1);
      }
      acc.addL(q2, cl);
      acc.addR(d, cr);
    };
    add(x, y, 0.0);
    for (const ExtNode& nd : stem_grow_nodes(x, y)) add(nd.d, nd.q, nd.w);
    Prof p;
    acc.finish(p);
    ++pf_pairs;
    return stP.emplace(key, std::move(p)).first->second;
  }

  // ------------------------------------------------ arm chains (--stem-est-chain)
  // A step from pair x.y out to a kept pair x'.y' enclosing it, with log weight w; t = 0 (b)
  // interior, 1 (f), 2 (g), 3 (h) multiloop.
  struct ChainNbr {
    int x, y, t;
    double w;
  };
  void count_steps(int which, const std::vector<ChainNbr>& nb) {
    ++st_nodes[which];
    for (const ChainNbr& s : nb) ++st_kind[which][s.t == 0 ? 0 : 1];
  }
  bool stem_chain = false, chain_ml = true;
  int chain_span = 150, chain_loop = 30;
  double chain_cut = 1e-8;
  // chain node profiles Prof*(x, y) and root profiles, by key2; valid for the whole run
  std::unordered_map<uint64_t, Prof> chM, chR;

  // Steps out of x.y (DP03._chain_nbrs): (b) x' = x-1 down to x-1-Lc, y' ascending in
  // (y, y+1+Lc], Lc = min(chain_loop, cap); with chain_ml and end_ok(x, y): (f) x' = i-1 over
  // the keys i of qmp_at(x-1) ascending (x - i <= S), y' in (y, y+1+ml_span]; (g) x' ascending
  // in [x-1-cap, x), y' with mbar(y+1, y'-1) and y'-y-1 <= S; (h) x' as (f), y' as (g).
  // x'.y' kept (and end_ok for f-h).
  void chain_nbrs(int x, int y, std::vector<ChainNbr>& out) {
    out.clear();
    const StemEst& st = *stem_;
    const int lc = std::min(chain_loop, cap);
    for (int xo = x - 1; xo >= std::max(0, x - 1 - lc); --xo) {
      const std::vector<int>& kr = st.kept_right[xo];
      for (auto it = std::upper_bound(kr.begin(), kr.end(), y);
           it != kr.end() && *it <= y + 1 + lc; ++it) {
        const double s = E.stem_interior(xo, x, y, *it);
        if (s < HALF_NEG) continue;
        out.push_back({xo, *it, 0, s + pair_bonus});
      }
    }
    if (!chain_ml || !E.end_ok(x, y)) return;
    const double br = E.w_bp + E.at(x, y);
    std::vector<std::pair<int, double>> qm;
    if (x >= 1) {
      const OMap& at = qmp_at(x - 1);
      for (size_t r = 0; r < at.slots(); ++r) {
        const int i = at.keys[r];
        if (i >= 1 && x - i <= chain_span) qm.emplace_back(i - 1, at.vals[r]);
      }
      std::sort(qm.begin(), qm.end());
    }
    auto close = [&](int xo, int yo) { return (E.w_ap + E.w_bp + E.at(xo, yo)) + pair_bonus; };
    for (const auto& [xo, qv] : qm) {
      const std::vector<int>& kr = st.kept_right[xo];
      for (auto it = std::upper_bound(kr.begin(), kr.end(), y);
           it != kr.end() && *it <= y + 1 + stem_ml_span; ++it) {
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
      const std::vector<int>& kr = st.kept_right[xo];
      for (auto it = std::upper_bound(kr.begin(), kr.end(), y + 1);
           it != kr.end() && *it - y - 1 <= chain_span; ++it) {
        const int yo = *it;
        const double mb = st.Mbar(y + 1, yo - 1);
        if (mb < HALF_NEG || !E.end_ok(xo, yo)) continue;
        out.push_back({xo, yo, 2, ((mb + br) + rl) + close(xo, yo)});
      }
    }
    for (const auto& [xo, qv] : qm) {
      const std::vector<int>& kr = st.kept_right[xo];
      for (auto it = std::upper_bound(kr.begin(), kr.end(), y + 1);
           it != kr.end() && *it - y - 1 <= chain_span; ++it) {
        const int yo = *it;
        const double mb = st.Mbar(y + 1, yo - 1);
        if (mb < HALF_NEG || !E.end_ok(xo, yo)) continue;
        out.push_back({xo, yo, 3, ((qv + mb) + br) + close(xo, yo)});
      }
    }
  }

  // Stop-here entries of node x.y with G = 1 (as pair_prof's add).
  void chain_stop(ProfAcc& acc, int x, int y) const {
    double cl = (prefix_log(x) + hr_c(y)) + (E.w_beta1 + 2.0 * E.w_stem);
    double cr = stem_->S[y + 1];
    if (stem_consistent) {
      cl += E.at(x, y);
      cr += (E.at(x, y) + E.w_beta1);
    }
    acc.addL(y, cl);
    acc.addR(x, cr);
  }

  // Adds w + Prof*(x', y') of every step in nb (step order; all in chM); counts the edges.
  void chain_merge(ProfAcc& acc, const std::vector<ChainNbr>& nb) {
    for (const ChainNbr& t : nb) {
      const Prof& q = chM.find(key2(t.x, t.y))->second;
      for (size_t k = 0; k < q.ly.size(); ++k) acc.addL(q.ly[k], t.w + q.lw[k]);
      for (size_t k = 0; k < q.rx.size(); ++k) acc.addR(q.rx[k], t.w + q.rw[k]);
      ++ch_edges[t.t];
    }
  }

  // Drop entries below log(chain_cut) + the part's max (order kept).
  long long chain_trunc(std::vector<int>& ix, std::vector<double>& w) const {
    if (chain_cut <= 0.0 || w.empty()) return 0;
    const double thr = *std::max_element(w.begin(), w.end()) + std::log(chain_cut);
    size_t o = 0;
    for (size_t k = 0; k < w.size(); ++k)
      if (w[k] >= thr) {
        ix[o] = ix[k];
        w[o] = w[k];
        ++o;
      }
    const long long d = static_cast<long long>(w.size() - o);
    ix.resize(o);
    w.resize(o);
    return d;
  }

  // Prof*(x, y) = [end_ok(x, y)] stop(x, y) + sum over the steps of w Prof*(x', y'), cut;
  // built bottom-up with an explicit stack (steps strictly enclose, so no cycles).
  void chain_build(int x, int y) {
    if (chM.count(key2(x, y))) return;
    std::vector<std::pair<int, int>> stk{{x, y}};
    std::vector<ChainNbr> nb;
    while (!stk.empty()) {
      const auto [a, b] = stk.back();
      if (chM.count(key2(a, b))) {
        stk.pop_back();
        continue;
      }
      chain_nbrs(a, b, nb);
      bool ready = true;
      for (const ChainNbr& t : nb)
        if (!chM.count(key2(t.x, t.y))) {
          stk.emplace_back(t.x, t.y);
          ready = false;
        }
      if (!ready) continue;
      ProfAcc acc;
      if (E.end_ok(a, b)) chain_stop(acc, a, b);
      chain_merge(acc, nb);
      Prof p;
      acc.finish(p);
      ch_cut_drop += chain_trunc(p.ly, p.lw) + chain_trunc(p.rx, p.rw);
      ++ch_nodes;
      chM.emplace(key2(a, b), std::move(p));
      stk.pop_back();
    }
  }

  // Root profile of the stem closed at x.y (DP03._chain_root): stop(x, y) unconditionally plus
  // the steps out of x.y; untruncated, cached for the run.
  const Prof& chain_root(int x, int y) {
    const uint64_t key = key2(x, y);
    auto it = chR.find(key);
    if (it != chR.end()) return it->second;
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<ChainNbr> nb;
    chain_nbrs(x, y, nb);
    for (const ChainNbr& t : nb) chain_build(t.x, t.y);
    ProfAcc acc;
    chain_stop(acc, x, y);
    chain_merge(acc, nb);
    Prof p;
    acc.finish(p);
    ++ch_roots;
    ch_ly += static_cast<long long>(p.ly.size());
    ch_rx += static_cast<long long>(p.rx.size());
    ch_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return chR.emplace(key, std::move(p)).first->second;
  }

  // Every chain path from x.y (root: stop unconditional; later nodes need end_ok) with its
  // weight summed in long double.
  void chain_paths(int x, int y, bool root, long double w,
                   std::vector<std::tuple<int, int, long double>>& out) {
    if (root || E.end_ok(x, y)) out.emplace_back(x, y, w);
    std::vector<ChainNbr> nb;
    chain_nbrs(x, y, nb);
    for (const ChainNbr& t : nb) chain_paths(t.x, t.y, false, w + t.w, out);
  }

  // --stem-profile-check with chain: the score of a state with inner pair (u, v) and log
  // weight a closed by cl, every chain path's stop terms summed separately (as helper_naive);
  // the path list of each root is enumerated once per run.
  std::unordered_map<uint64_t, std::vector<std::tuple<int, int, long double>>> chPaths;
  double chain_naive(const std::vector<Closing>& cl, int u, int v, double a) {
    std::vector<long double> terms;
    const long double w_pk = static_cast<long double>(E.w_beta1) + 2.0L * E.w_stem;
    for (const Closing& c : cl) {
      auto it = chPaths.find(key2(c.x, c.y));
      if (it == chPaths.end()) {
        it = chPaths.emplace(key2(c.x, c.y), std::vector<std::tuple<int, int, long double>>())
                 .first;
        chain_paths(c.x, c.y, true, 0.0L, it->second);
      }
      for (const auto& [d, q, w] : it->second) {
        const long double base = static_cast<long double>(c.w) + a + w;
        const double h = hr_h(u, v, q);
        if (h > HALF_NEG) {
          long double s = base + prefix_log(d) + h + hr_c(q) + w_pk;
          if (score_qz) s += qz_log(u + 1, v - 1);
          if (stem_consistent) s += static_cast<long double>(E.at(u, v)) + E.at(d, q);
          terms.push_back(s);
        }
        const double bg = stem_Bg_raw(d, u, v);
        if (bg > HALF_NEG) {
          long double s = base + bg + stem_->S[q + 1];
          if (stem_consistent) s += static_cast<long double>(E.at(d, q)) + E.w_beta1;
          terms.push_back(s);
        }
      }
    }
    if (terms.empty()) return NEG_INF;
    long double m = *std::max_element(terms.begin(), terms.end()), sum = 0.0L;
    for (long double t : terms) sum += std::exp(t - m);
    return static_cast<double>(m + std::log(sum));
  }

  // ------------------------------------------ principled score (--stem-est-principled)
  // sigma(X) = alpha_X + log Out_hat(X); Out_hat follows DP03's edges outward from X: exact
  // scan tables at positions <= q, the nested pre-pass (StemEst pr*) beyond q, and outside
  // the whole pseudoknot Ctx(i, j) = e^beta1 C[i-1] Sbar[j+1] + e^beta1m Bml[i][j].
  // Stem steps from a pair outward (pr_steps) run over the pre-pass kept pairs.  The growth G
  // of the scored stem takes Qm from the scan on its left; future stems H (profiles over
  // their outer ends) take Mbar on both sides.
  bool pr_on = false, pr_check = false, pr_nested = false, pr_ml = true, pr_built = false;
  int pr_loop = 30, pr_ml_span = 150, pr_ml_close = 30, pr_arms_k = 64;
  bool pr_ratio = false, pr_gapfree = false, pr_cut_roots = false;
  double pr_arms_eps = 1e-3;
  int pr_arms_max = 256;
  double pr_cut = 1e-8;
  // --pr-arms sep: the exterior term of the L-stem arm part over all arms, the gap filler made
  // separable: Qz[u+1][v-1] S_q'(u, v), S = sum over arms xs >= u+1 (two: = u+1), x' <= v-1
  // of e^phi / Zarm (arm_z 0: Qz[xs][x'], 1: e^(wu (x' - xs + 1))).  The multiloop-context
  // term (Bml(p', ys), not one factor per arm) keeps the exact-placement arms of pr_select.
  // qz_unp: the score's gap filler unpaired-only (test); sep_check: every distinct query
  // against brute force and exact placement.
  int pr_arms_mode = 0, pr_arm_z = 0;
  bool pr_qz_unp = false, pr_sep_check = false;
  bool pr_gsel = false, pr_gcheck = false;
  double pr_geps = 1e-6;
  int pr_gB = 512;
  FILE* pr_dump = nullptr;

 public:
  // Growth instrumentation, per profile kind (0 node, 1 root, 2 key): entries before / after
  // selection, distinct p' and q' after (and before) selection; under --grow-rec the nonzero
  // q' of g_L and p' of g_R (all built layers).  Memory in bytes of the stored growth data:
  // nodes and roots (kept for the run), helpers (cleared at PR_HG_MAX), the per-step root /
  // key caches; peaks over the run.
  bool rec_noml = false;
  std::vector<int> gz_before[3], gz_after[3], gz_dp[3], gz_dq[3], gz_dpb[3], gz_dqb[3];
  int gs_who = 0;
  double gm_nodes = 0, gm_roots = 0, gm_help = 0, gm_help_peak = 0, gm_step = 0, gm_step_peak = 0,
         gm_total_peak = 0;
  long long gm_help_entries_peak = 0;
