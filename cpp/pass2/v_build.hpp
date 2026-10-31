// Pass 2, V^pre: layered profile finish and the b_V cut (pr_vtrim: keep the top-b_V entries
// xs.ys by V Sbar_pre[ys+1], --pr-v-trim suffix), V dumps, V suffix / group builds and
// pr_static_ats (static arm table).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.


  // --pr-v-eps: keep the entries xs.ys (layered: xs.ys.c) of the map of a pair with left end x
  // by V Sbar[ys+1] / Qzbar[xs][x], largest first, until 1 - eps of their total (key order kept).
  // --pr-v-trim suffix: the same weights, but xs.ys is dropped only if below eps times the sum
  // over the entries with xs' >= xs (the lookups take xs > a2, so only those always come with
  // it); then at most pr_vcap entries, by weight over that suffix sum.
  void pr_vtrim_suffix(PProf& p, int x, bool layered) {
    const StemEst& st = *stem_;
    const size_t m = p.k.size();
    std::vector<double> io(m);
    std::vector<int> xsv(m);
    for (size_t i = 0; i < m; ++i) {
      const int xs = layered ? kcl_xs(p.k[i]) : k_hi(p.k[i]);
      const int ys = layered ? kcl_ys(p.k[i]) : k_lo(p.k[i]);
      const double z = st.Qzbar(xs, x);
      io[i] = p.w[i] + (pr_vxs ? 0.0 : st.prS[ys + 1]) - (z > HALF_NEG ? z : 0.0);
      xsv[i] = xs;
    }
    std::vector<size_t> ord(m);
    for (size_t i = 0; i < m; ++i) ord[i] = i;
    std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return xsv[a] > xsv[b]; });
    std::vector<double> rel(m, NEG_INF);
    const double le = pr_veps > 0.0 ? std::log(pr_veps) : NEG_INF;
    double suf = NEG_INF;
    for (size_t a = 0; a < m;) {
      size_t b = a;
      while (b < m && xsv[ord[b]] == xsv[ord[a]]) suf = log_add(suf, io[ord[b++]]);
      for (size_t c = a; c < b; ++c) rel[ord[c]] = io[ord[c]] - suf;
      a = b;
    }
    std::vector<char> keep(m, 0);
    std::vector<size_t> kept;
    for (size_t i = 0; i < m; ++i)
      if (io[i] < HALF_NEG ? false : rel[i] >= le) kept.push_back(i);
    if (pr_vcap > 0 && kept.size() > static_cast<size_t>(pr_vcap)) {
      v_capped += static_cast<long long>(kept.size()) - pr_vcap;
      std::stable_sort(kept.begin(), kept.end(), [&](size_t a, size_t b) { return rel[a] > rel[b]; });
      double ai = NEG_INF, di = NEG_INF, aw = NEG_INF, dw = NEG_INF;
      for (size_t c = 0; c < kept.size(); ++c) {
        ai = log_add(ai, io[kept[c]]);
        aw = log_add(aw, p.w[kept[c]]);
        if (c >= static_cast<size_t>(pr_vcap)) {
          di = log_add(di, io[kept[c]]);
          dw = log_add(dw, p.w[kept[c]]);
        }
      }
      ++v_cap_pairs;
      const double fi = std::exp(di - ai), fw = std::exp(dw - aw);
      v_cap_fsum += fi;
      v_cap_fmax = std::max(v_cap_fmax, fi);
      v_cap_fwmax = std::max(v_cap_fwmax, fw);
      v_cap_all = log_add(v_cap_all, aw);
      v_cap_drop = log_add(v_cap_drop, dw);
      kept.resize(static_cast<size_t>(pr_vcap));
    }
    for (size_t i : kept) keep[i] = 1;
    size_t o = 0;
    for (size_t i = 0; i < m; ++i)
      if (keep[i]) {
        p.k[o] = p.k[i];
        p.w[o] = p.w[i];
        ++o;
      }
    p.k.resize(o);
    p.w.resize(o);
  }
  // --bV n_V: keep the top n_V left ends b of one V map, ranked by C_pre(b-1) * V.
  void pr_vtrim_cpre(PProf& p, bool layered) {
    if (p2_v_cap <= 0 || p.k.empty() || !stem_ || stem_->prC.empty()) return;
    const size_t m = p.k.size();
    std::vector<int> xs(m);
    for (size_t i = 0; i < m; ++i) xs[i] = layered ? kcl_xs(p.k[i]) : k_hi(p.k[i]);
    std::vector<size_t> ord(m);
    for (size_t i = 0; i < m; ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return xs[a] < xs[b]; });
    std::vector<std::pair<int, double>> grp;
    for (size_t a = 0; a < m;) {
      size_t b = a;
      double s = NEG_INF;
      const int x0 = xs[ord[a]];
      while (b < m && xs[ord[b]] == x0) s = log_add(s, p.w[ord[b++]]);
      const double c = x0 > 0 ? stem_->prC[static_cast<size_t>(x0 - 1)] : 0.0;
      grp.emplace_back(x0, s + c);
      a = b;
    }
    if (static_cast<int>(grp.size()) <= p2_v_cap) return;
    std::nth_element(grp.begin(), grp.begin() + p2_v_cap, grp.end(),
                     [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                       return a.second > b.second;
                     });
    std::vector<char> keepx(static_cast<size_t>(n), 0);
    for (int i = 0; i < p2_v_cap; ++i)
      if (grp[static_cast<size_t>(i)].first >= 0 && grp[static_cast<size_t>(i)].first < n)
        keepx[static_cast<size_t>(grp[static_cast<size_t>(i)].first)] = 1;
    size_t o = 0;
    for (size_t i = 0; i < m; ++i)
      if (xs[i] >= 0 && xs[i] < n && keepx[static_cast<size_t>(xs[i])]) {
        p.k[o] = p.k[i];
        p.w[o] = p.w[i];
        ++o;
      }
    p.k.resize(o);
    p.w.resize(o);
  }
  void pr_vtrim(PProf& p, int x, bool layered = false) {
    v_before += static_cast<long long>(p.k.size());
    if (pr_vsuffix && p.k.size() > 1) {
    pr_vtrim_suffix(p, x, layered);
  } else if (pr_veps > 0.0 && p.k.size() > 1) {
    const StemEst& st = *stem_;
    const size_t m = p.k.size();
    std::vector<double> io(m);
    double tot = NEG_INF;
    for (size_t i = 0; i < m; ++i) {
      const int xs = layered ? kcl_xs(p.k[i]) : k_hi(p.k[i]);
      const int ys = layered ? kcl_ys(p.k[i]) : k_lo(p.k[i]);
      const double z = st.Qzbar(xs, x);
      io[i] = p.w[i] + (pr_vxs ? 0.0 : st.prS[ys + 1]) - (z > HALF_NEG ? z : 0.0);
      tot = log_add(tot, io[i]);
    }
    if (tot > HALF_NEG) {
      std::vector<size_t> ord(m);
      for (size_t i = 0; i < m; ++i) ord[i] = i;
      std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return io[a] > io[b]; });
      const double goal = tot + std::log1p(-pr_veps);
      std::vector<char> keep(m, 0);
      double s = NEG_INF;
      for (size_t i : ord) {
        keep[i] = 1;
        s = log_add(s, io[i]);
        if (s >= goal) break;
      }
      size_t o = 0;
      for (size_t i = 0; i < m; ++i)
        if (keep[i]) {
          p.k[o] = p.k[i];
          p.w[o] = p.w[i];
          ++o;
        }
      p.k.resize(o);
      p.w.resize(o);
    }
  }
    if (p2_v_cap > 0) pr_vtrim_cpre(p, layered);
    v_after += static_cast<long long>(p.k.size());
  }


  void pr_static_ats(bool count, const std::vector<std::pair<int, int>>& pairs,
                     std::vector<AT>& ats, std::vector<std::vector<AT>>* bk = nullptr) {
    const StemEst& st = *stem_;
    std::unordered_map<uint64_t, PProf> V;
    V.reserve(pairs.size() * 2 + 1);
    std::vector<ChainNbr> nb;
    if (pr_fast > 0) {
      // y descending: every step's x'.y' has y' > y, so it is built first, and the helpers
      // of one y are dropped when y changes.
      std::vector<std::pair<int, int>> ord(pairs);
      std::stable_sort(ord.begin(), ord.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
        return a.second > b.second;
      });
      std::unordered_map<uint64_t, PProf> H;
      auto child = [&](int a, int b) -> const PProf* { return &V.find(key2(a, b))->second; };
      // pass 0: the stop at x.y (end_ok); pass 1: the pair alone (single_stems)
      auto node = [&](int x, int y, bool inner, std::unordered_map<uint64_t, PProf>& hc, PProf& p) {
        std::chrono::steady_clock::time_point t0;
        pr_steps(x, y, inner, false, nb, pr_fast < 2);
        if (count) count_steps(1, nb);
        double self = !inner ? (E.end_ok(x, y) ? E.at(x, y) : NEG_INF)
                             : (single_stems ? E.at(x, y) : NEG_INF);
        if (pr_vxs && self > HALF_NEG) self += E.w_beta1 + st.prS[y + 1];
        pr_node_acc(x, y, inner, false, nb, self, child, hc);
        const long long d = pr_finish_fut_d(p, x);
        pr_vtrim(p, x);
        if (count) pr_cut_drop += d;
      };
      auto put_v = [&](int x, int y, PProf& p) {
        if (count) {
          pr_V_entries += static_cast<long long>(p.k.size());
          ++v_pairs;
        }
        pr_linify(p);
        V.emplace(key2(x, y), std::move(p));
      };
      auto put_ats = [&](int x, int y, const PProf& p) {
        const double base = (E.at(x, y) + E.w_stem) + pair_bonus;
        for (size_t k = 0; k < p.k.size(); ++k) {
          const AT a{k_hi(p.k[k]), x, y, k_lo(p.k[k]), base + p.w[k]};
          if (bk) (*bk)[a.xs].push_back(a);
          else ats.push_back(a);
        }
      };
      if (pr_vstream) {
        // V(x', y') is read only by pairs with y' - y - 1 <= R (pr_steps / pr_helper reach)
        const int lc = std::min(pr_loop, cap);
        const int R = std::max({lc, cap, pr_ml_span});
        std::vector<std::vector<int>> byY(n);
        std::unordered_map<uint64_t, PProf> H1;
        long long live = 0;
        int yfree = n - 1, ylast = -1;
        for (const auto& [x, y] : ord) {
          cpu_stop_check(1);
          if (y != ylast) {
            H.clear();
            H1.clear();
            ylast = y;
            for (; yfree > y + 1 + R; --yfree) {
              for (int a : byY[yfree]) {
                auto it = V.find(key2(a, yfree));
                if (g_rfv) rfv_put(it->first, it->second);
                live -= static_cast<long long>(it->second.k.size());
                V.erase(it);
              }
              std::vector<int>().swap(byY[yfree]);
            }
          }
          PProf p;
          node(x, y, false, H, p);
          live += static_cast<long long>(p.k.size());
          v_live_max = std::max(v_live_max, live);
          put_v(x, y, p);
          byY[y].push_back(x);
          if (E.end_ok(x, y)) {
            PProf q;
            node(x, y, true, H1, q);
            put_ats(x, y, q);
          }
        }
        if (g_rfv) rfv_keep(V);
        return;
      }
      for (int pass = 0; pass < 2; ++pass) {
        H.clear();
        int ylast = -1;
        for (const auto& [x, y] : ord) {
          cpu_stop_check(1);
          if (pass == 1 && !E.end_ok(x, y)) continue;
          if (y != ylast) {
            H.clear();
            ylast = y;
          }
          PProf p;
          node(x, y, pass == 1, H, p);
          if (pass == 0) put_v(x, y, p);
          else put_ats(x, y, p);
        }
        if (pass == 0) v_live_max = pr_V_entries;
      }
      if (g_rfv) rfv_keep(V);
      return;
    }
    {
      auto merge = [&](PAcc& acc) {
        for (const ChainNbr& t : nb) {
          const PProf& q = V.find(key2(t.x, t.y))->second;
          for (size_t k = 0; k < q.k.size(); ++k) acc.add(q.k[k], t.w + q.w[k]);
        }
      };
      for (const auto& [x, y] : pairs) {
        cpu_stop_check(1);
        pr_steps(x, y, false, false, nb);
        if (count) count_steps(1, nb);
        PAcc acc;
        if (E.end_ok(x, y)) acc.add(key2(x, y), E.at(x, y));
        merge(acc);
        PProf p;
        const long long d = pr_finish_fut(acc, p, x);
        pr_vtrim(p, x);
        if (count) {
          pr_cut_drop += d;
          pr_V_entries += static_cast<long long>(p.k.size());
        }
        V.emplace(key2(x, y), std::move(p));
      }
      for (const auto& [x, y] : pairs) {
        cpu_stop_check(1);
        if (!E.end_ok(x, y)) continue;
        pr_steps(x, y, true, false, nb);
        if (count) count_steps(1, nb);
        PAcc acc;
        if (single_stems) acc.add(key2(x, y), E.at(x, y));
        merge(acc);
        PProf p;
        const long long d = pr_finish_fut(acc, p, x);
        pr_vtrim(p, x);
        if (count) pr_cut_drop += d;
        const double base = (E.at(x, y) + E.w_stem) + pair_bonus;
        for (size_t k = 0; k < p.k.size(); ++k)
          ats.push_back({k_hi(p.k[k]), x, y, k_lo(p.k[k]), base + p.w[k]});
      }
    }
  }

  // --grow-r-from-v: V kept after the build as flat arrays (xs ascending, float log V) per pair
  // (prVkI: key2(x, y) -> offset, length), filled as the build frees each pair's profile;
  // prVod: V of pairs outside V's pair set (key3(x, y, inner)), built from the kept V of their
  // steps with the growth vectors' stop at x.y; prVodH: their helpers; prVch: the kept V of a
  // node's steps as profiles while it is built; prVkey: per step, the U keys' sums over their
  // closings (index order).  Counters: lookups, lookups with no entry, seed lookups through the
  // arm entries, on-demand nodes and entries, steps onto pairs not kept, keys and entries.
 public:
  struct VSpan {
    const int32_t* xs = nullptr;
    const float* lw = nullptr;
    uint32_t len = 0;
  };
  struct VOwn {
    std::vector<int32_t> xs;
    std::vector<float> lw;
    VSpan span() const { return {xs.data(), lw.data(), static_cast<uint32_t>(xs.size())}; }
  };
  std::unordered_map<uint64_t, std::pair<uint32_t, uint32_t>> prVkI;
  std::vector<int32_t> prVkX;
  std::vector<float> prVkW;
  std::unordered_map<uint64_t, VOwn> prVod;
  std::unordered_map<uint64_t, PProf> prVodH, prVch;
  std::unordered_map<int, std::vector<std::pair<int, double>>> prVkey[2];
  std::vector<double> rfvAcc;
  std::vector<int> rfvT;
  long long rfv_look = 0, rfv_miss = 0, rfv_seed = 0, rfv_od = 0, rfv_od_ent = 0,
            rfv_child_miss = 0, rfv_keys = 0, rfv_key_ent = 0;
  double rfv_sec = 0;
  // v-sweep lookups split by U keys (rfv_key) and R/K terms (rfv_root). key0 is a state
  // scored in a column over its beam whose key is -inf because the V read missed.
  long long rfv_look_u = 0, rfv_miss_u = 0, rfv_look_r = 0, rfv_miss_r = 0;
  long long rfv_key0_u = 0, rfv_part_u = 0, rfv_key0_r = 0;
  long long rfv_look_s = 0, rfv_miss_s = 0, vfb_n = 0;
  // Miss kind: (i) map full (b_V), entry cut; (ii) map exists, not full; (iii) no map.
  enum { ST_P1 = 0, ST_P2, ST_P3, ST_R, ST_S, ST_U, ST_N };
  long long vm_i[ST_N] = {}, vm_ii[ST_N] = {}, vm_iii[ST_N] = {};
  long long km_i[ST_N] = {}, km_ii[ST_N] = {}, km_iii[ST_N] = {};
  long long hm_i[ST_N] = {}, hm_ii[ST_N] = {}, hm_iii[ST_N] = {};
