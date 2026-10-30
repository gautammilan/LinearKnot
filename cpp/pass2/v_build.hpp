// Pass 2, V^pre: layered profile finish and the b_V cut (pr_vtrim: keep the top-b_V entries
// xs.ys by V Sbar_pre[ys+1], --pr-v-trim suffix), V dumps, V suffix / group builds and
// pr_static_ats (static arm table).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  // Layered profile of a future stem with left end x: entries sorted, cut by the per-(xs, ys)
  // sum over c (pr_finish_fut's rule); returns the dropped (xs, ys); groups: kept (xs, ys).
  long long pr_finish_layers(PAcc& acc, PProf& p, int x, long long& groups) const {
    std::vector<std::pair<uint64_t, double>>& v = acc.v;
    std::sort(v.begin(), v.end(), [](const std::pair<uint64_t, double>& a,
                                     const std::pair<uint64_t, double>& b) {
      return a.first < b.first;
    });
    const StemEst& st = *stem_;
    auto adj = [&](uint64_t k) {
      return pr_gapfree ? st.prS[kcl_ys(k) + 1] - E.wu * (x - kcl_xs(k)) : 0.0;
    };
    std::vector<size_t> g0;
    std::vector<double> col;
    for (size_t i = 0; i < v.size(); ++i) {
      if (i == 0 || (v[i].first >> 8) != (v[i - 1].first >> 8)) {
        g0.push_back(i);
        col.push_back(v[i].second);
      } else {
        col.back() = log_add(col.back(), v[i].second);
      }
    }
    g0.push_back(v.size());
    double thr = NEG_INF;
    if (pr_cut > 0.0 && !col.empty()) {
      double mx = NEG_INF;
      for (size_t g = 0; g < col.size(); ++g) mx = std::max(mx, col[g] + adj(v[g0[g]].first));
      thr = mx + std::log(pr_cut);
    }
    long long d = 0;
    for (size_t g = 0; g < col.size(); ++g) {
      if (col[g] + adj(v[g0[g]].first) >= thr) {
        ++groups;
        for (size_t i = g0[g]; i < g0[g + 1]; ++i) {
          p.k.push_back(v[i].first);
          p.w.push_back(v[i].second);
        }
      } else {
        ++d;
      }
    }
    return d;
  }

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
  // --v-cap n_V: keep the top n_V left ends b of one V map, ranked by C_pre(b-1) * V.
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
    if (!(p2_v_cap > 0 && p2_v_trim == 2)) {
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
    }
    if (p2_v_cap > 0 && p2_v_trim != 1) pr_vtrim_cpre(p, layered);
    v_after += static_cast<long long>(p.k.size());
  }

  // --pr-v-diag / --pr-v-ratio: Vr(xs; x', y') = sum over partner pairs xs.ys, ys > y', of
  // Qb(xs, ys) / Qb(x', y') with the entry's end factors (base(x', y') at(xs, ys) e^beta1
  // Sbar[ys+1]); xs = x': the pair alone.  prVsuf[xs][k]: log suffix sum over the k-th and
  // later partners ys of xs of Qb(xs, ys) at(xs, ys) e^beta1 Sbar[ys+1] (end_ok).
  bool pr_vratio = false;
  std::string pr_vdiag;
  // --pr-v2-dump PATH: per full-support group (xs; x', y') log Vr and log Vr2, where Vr2 takes
  // per partner ys min(Qb(xs,ys)/Qb(x',y'), Qout(x',y')/Qout(xs,ys)) = the inside ratio times
  // min(1, P(x'.y') / P(xs.ys)); then exits.
  std::string pr_v2dump;
  void pr_v2_dump(const std::vector<std::pair<int, int>>& pairs) {
    const StemEst& st = *stem_;
    const auto& KR = pr_KRp();
    FILE* f = std::fopen(pr_v2dump.c_str(), "w");
    if (!f) throw std::runtime_error("cannot open --pr-v2-dump file");
    std::fprintf(f, "xs\txp\typ\tlogVr\tlogVr2\n");
    long long nin = 0, nout = 0, groups = 0;
    for (const auto& [x, y] : pairs) {
      if (!E.end_ok(x, y)) continue;
      const double qi = st.QbIn(x, y);
      const double lpi = st.prLP[static_cast<size_t>(x) * n + y];
      if (qi < HALF_NEG || lpi < HALF_NEG) continue;
      const double base = (E.at(x, y) + E.w_stem) + pair_bonus;
      for (int xs = 0; xs < x; ++xs) {
        double a = NEG_INF, b = NEG_INF;
        for (auto it = std::upper_bound(KR[xs].begin(), KR[xs].end(), y); it != KR[xs].end(); ++it) {
          const int ys = *it;
          const double qb = st.QbIn(xs, ys);
          if (qb < HALF_NEG || !E.end_ok(xs, ys)) continue;
          const double t = qb + E.at(xs, ys) + (E.w_beta1 + st.prS[ys + 1]);
          const double g = std::min(0.0, lpi - st.prLP[static_cast<size_t>(xs) * n + ys]);
          ++(g < 0.0 ? nout : nin);
          a = log_add(a, t);
          b = log_add(b, t + g);
        }
        if (a < HALF_NEG) continue;
        ++groups;
        std::fprintf(f, "%d\t%d\t%d\t%.9g\t%.9g\n", xs, x, y, base - qi + a, base - qi + b);
      }
    }
    std::fclose(f);
    std::fprintf(stderr, "[v2] groups %lld terms_inside %lld terms_outside %lld\n", groups, nin, nout);
  }
  double pr_vsample_p = 0.0;
  long long pr_vsample_max = 10000;
  std::vector<std::vector<double>> prVsuf, prVsufO;
  std::vector<double> prVrg, prVint, prVlen, prVcon, prVshare;
  // --pr-v-diag, per group: Vr2 (per partner min of the inside and outside ratios), Vsep
  // (min of the two sums), Vo (the outside-ratio sum), LB (sum of the per-partner lower bounds
  // Qb(xs,ys)/Qb(x',y') max(0, P(x'.y') + P(xs.ys) - 1) / P(xs.ys)).
  std::vector<double> prVr2g, prVsepg, prVog, prVlbg;
  std::vector<long long> prVuse;
  std::unordered_map<uint64_t, int> prArmOf;
  FILE* pr_vsel = nullptr;
  // --pr-v-sample: PREFIX.arms.tsv.gz, per sampled lookup and arm (xs, x') of the gap, the
  // partner value with the gap fillers under exact V, Vr, Vr2, Vsep (log; -inf: none)
  FILE* pr_varms = nullptr;
  long long pr_vsel_n = 0, pr_vT_n = 0, pr_vr_groups = 0, pr_vr_extra = 0;
  double pr_vr_sec = 0.0, pr_vint_sec = 0.0;
  uint64_t pr_vrng = 0x9E3779B97F4A7C15ull;
  double pr_vrand() {
    pr_vrng ^= pr_vrng << 13;
    pr_vrng ^= pr_vrng >> 7;
    pr_vrng ^= pr_vrng << 17;
    return static_cast<double>(pr_vrng >> 11) * (1.0 / 9007199254740992.0);
  }
  const std::vector<std::vector<int>>& pr_KRp() const {
    const StemEst& st = *stem_;
    return st.pr_partner ? st.prKeptRightP : st.prKeptRight;
  }
  void pr_vsuf_build() {
    if (!prVsuf.empty()) return;
    const StemEst& st = *stem_;
    const auto& KR = pr_KRp();
    prVsuf.assign(n, {});
    prVsufO.assign(n, {});
    for (int xs = 0; xs < n; ++xs) {
      const std::vector<int>& kr = KR[xs];
      std::vector<double>& s = prVsuf[xs];
      std::vector<double>& so = prVsufO[xs];
      s.assign(kr.size() + 1, NEG_INF);
      so.assign(kr.size() + 1, NEG_INF);
      for (int k = static_cast<int>(kr.size()) - 1; k >= 0; --k) {
        const int ys = kr[k];
        const double qb = st.QbIn(xs, ys);
        double t = NEG_INF;
        if (qb > HALF_NEG && E.end_ok(xs, ys)) t = qb + E.at(xs, ys) + (E.w_beta1 + st.prS[ys + 1]);
        s[k] = log_add(s[k + 1], t);
        const double lps = st.prLP.empty() ? NEG_INF : st.prLP[static_cast<size_t>(xs) * n + ys];
        so[k] = log_add(so[k + 1], t > HALF_NEG && lps > HALF_NEG ? t - lps : NEG_INF);
      }
    }
  }
  // kind 0: Vsep = min(Vr, Vo); 1: Vr2; 2: LB; 3: Vo = P(x'.y') sum_ys t / P(xs.ys) / Qb(x',y').
  // xs = x': the pair alone (exact).  Same units as pr_vr_group.
  double pr_vb_group(int xs, int xp, int yp, int kind) const {
    const StemEst& st = *stem_;
    if (xs == xp || st.prLP.empty()) return pr_vr_group(xs, xp, yp);
    const double base = (E.at(xp, yp) + E.w_stem) + pair_bonus;
    const double qi = st.QbIn(xp, yp);
    const double lpi = st.prLP[static_cast<size_t>(xp) * n + yp];
    if (qi < HALF_NEG || lpi < HALF_NEG) return NEG_INF;
    const std::vector<int>& kr = pr_KRp()[xs];
    const size_t k0 = std::upper_bound(kr.begin(), kr.end(), yp) - kr.begin();
    if (kind == 0 || kind == 3) {
      const double so = prVsufO[xs][k0];
      const double vo = so > HALF_NEG ? base - qi + lpi + so : NEG_INF;
      if (kind == 3) return vo;
      const double s = prVsuf[xs][k0];
      const double vr = s > HALF_NEG ? base - qi + s : NEG_INF;
      return std::min(vr, vo);
    }
    const uint64_t key = key3(xs, xp, yp);
    auto& memo = prVbc[kind - 1];
    if (auto it = memo.find(key); it != memo.end()) return it->second;
    const double pi = std::exp(lpi);
    double acc = NEG_INF;
    for (size_t k = k0; k < kr.size(); ++k) {
      const int ys = kr[k];
      const double qb = st.QbIn(xs, ys);
      const double lps = st.prLP[static_cast<size_t>(xs) * n + ys];
      if (qb < HALF_NEG || lps < HALF_NEG || !E.end_ok(xs, ys)) continue;
      const double t = qb + E.at(xs, ys) + (E.w_beta1 + st.prS[ys + 1]);
      if (kind == 1) {
        acc = log_add(acc, t + std::min(0.0, lpi - lps));
      } else {
        const double j = pi + std::expm1(lps);  // P(x'.y') + P(xs.ys) - 1
        if (j > 0.0) acc = log_add(acc, t + std::log(j) - lps);
      }
    }
    const double out = acc > HALF_NEG ? base - qi + acc : NEG_INF;
    memo.emplace(key, out);
    return out;
  }
  mutable std::unordered_map<uint64_t, double> prVbc[2];
  // log Vr(xs; x', y') (the V group value, prYw's units); x'.y' end_ok.
  double pr_vr_group(int xs, int xp, int yp) const {
    const StemEst& st = *stem_;
    const double base = (E.at(xp, yp) + E.w_stem) + pair_bonus;
    if (xs == xp)
      return single_stems ? base + E.at(xp, yp) + (E.w_beta1 + st.prS[yp + 1]) : NEG_INF;
    const double qi = st.QbIn(xp, yp);
    if (qi < HALF_NEG) return NEG_INF;
    const std::vector<int>& kr = pr_KRp()[xs];
    const size_t k = std::upper_bound(kr.begin(), kr.end(), yp) - kr.begin();
    const double s = prVsuf[xs][k];
    return s > HALF_NEG ? base - qi + s : NEG_INF;
  }

  void pr_static_ats(int K, bool count, const std::vector<std::pair<int, int>>& pairs,
                     std::vector<AT>& ats, std::vector<std::vector<AT>>* bk = nullptr) {
    const StemEst& st = *stem_;
    if (pr_vratio && K == 0) {
      const auto& KR = pr_KRp();
      for (const auto& [x, y] : pairs) {
        if (!E.end_ok(x, y)) continue;
        const double base = (E.at(x, y) + E.w_stem) + pair_bonus;
        if (single_stems) ats.push_back({x, x, y, y, base + E.at(x, y)});
        const double qi = st.QbIn(x, y);
        if (qi < HALF_NEG) continue;
        for (int xs = 0; xs < x; ++xs)
          for (auto it = std::upper_bound(KR[xs].begin(), KR[xs].end(), y); it != KR[xs].end(); ++it) {
            const int ys = *it;
            const double qb = st.QbIn(xs, ys);
            if (qb < HALF_NEG || !E.end_ok(xs, ys)) continue;
            const double g =
                g_v2ratio ? std::min(0.0, st.prLP[static_cast<size_t>(x) * n + y] -
                                              st.prLP[static_cast<size_t>(xs) * n + ys])
                          : 0.0;
            ats.push_back({xs, x, y, ys, base + E.at(xs, ys) + qb - qi + g});
            if (count) ++pr_V_entries;
          }
      }
      return;
    }
    std::unordered_map<uint64_t, PProf> V;
    V.reserve(pairs.size() * 2 + 1);
    std::vector<ChainNbr> nb;
    if (K == 0 && pr_fast > 0) {
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
        if (lpf) t0 = std::chrono::steady_clock::now();
        pr_steps(x, y, inner, false, nb, pr_fast < 2 || pr_vmlsteps);
        if (count) count_steps(1, nb);
        if (lpf) {
          ++vi_nodes[inner ? 1 : 0];
          for (const ChainNbr& t : nb) ++(t.t == 0 ? vi_nb_int : vi_nb_ml);
          vi_lap(t0, vi_steps_s);
        }
        double self = !inner ? (E.end_ok(x, y) ? E.at(x, y) : NEG_INF)
                             : (single_stems ? E.at(x, y) : NEG_INF);
        if (pr_vxs && self > HALF_NEG) self += E.w_beta1 + st.prS[y + 1];
        pr_node_acc(x, y, inner, false, nb, self, child, hc);
        if (lpf) t0 = std::chrono::steady_clock::now();
        const long long d = pr_finish_fut_d(p, x);
        if (lpf) vi_ent_fin += static_cast<long long>(p.k.size());
        pr_vtrim(p, x);
        if (lpf) {
          vi_ent_kept += static_cast<long long>(p.k.size());
          vi_lap(t0, vi_fin_s);
        }
        if (count) pr_cut_drop += d;
      };
      auto put_v = [&](int x, int y, PProf& p) {
        if (sw_check) sw_ref.emplace(key2(x, y), p);
        if (count) {
          pr_V_entries += static_cast<long long>(p.k.size());
          ++v_pairs;
        }
        pr_linify(p);
        V.emplace(key2(x, y), std::move(p));
      };
      auto put_ats = [&](int x, int y, const PProf& p) {
        std::unique_ptr<SecTimer> ptm(lpf ? new SecTimer(vi_ats_s) : nullptr);
        if (lpf) vi_ats += static_cast<long long>(p.k.size());
        const double base = (E.at(x, y) + E.w_stem) + pair_bonus;
        for (size_t k = 0; k < p.k.size(); ++k) {
          const AT a{k_hi(p.k[k]), x, y, k_lo(p.k[k]), base + p.w[k]};
          if (bk) (*bk)[a.xs].push_back(a);
          else ats.push_back(a);
        }
      };
      auto ats_sz = [&](long long& e) {
        if (!bk) {
          e = static_cast<long long>(ats.size());
          return msz::vec(ats);
        }
        double b = 0.0;
        e = 0;
        for (const std::vector<AT>& v : *bk) {
          e += static_cast<long long>(v.size());
          b += msz::vec(v);
        }
        return b;
      };
      if (pr_vstream) {
        // V(x', y') is read only by pairs with y' - y - 1 <= R (pr_steps / pr_helper reach)
        const int lc = std::min(pr_loop, cap);
        const int R = pr_ml && pr_partner_ml ? std::max({lc, cap, pr_ml_span}) : lc;
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
        if (g_mem_log) {
          long long ve = 0, he = 0, ae = 0;
          const double vb = pmap_sz(V, ve), hb = pmap_sz(H, he) + pmap_sz(H1, he);
          const double ab = ats_sz(ae);
          mem_log("V_stream", {{"V_prof", ve, vb},
                               {"V_help", he, hb},
                               {"V_ats", ae, ab},
                               {"V_live_max", v_live_max, 0.0},
                               {"V_untrimmed", v_before, 0.0},
                               {"V_trimmed", v_after, 0.0}});
        }
        vquery_dump(V);
        mu_dump_lazy(V);
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
        if (g_mem_log) {
          long long ve = 0, he = 0, ae = 0;
          const double vb = pmap_sz(V, ve), hb = pmap_sz(H, he);
          const double ab = ats_sz(ae);
          mem_log(pass == 0 ? "V_pass0" : "V_pass1",
                  {{"V_prof", ve, vb},
                   {"V_help", he, hb},
                   {"V_ats", ae, ab},
                   {"V_untrimmed", v_before, 0.0},
                   {"V_trimmed", v_after, 0.0}});
        }
      }
      vquery_dump(V);
      if (g_rfv) rfv_keep(V);
      return;
    }
    if (K == 0) {
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
      return;
    }
    const uint64_t cmax = static_cast<uint64_t>(K - 1);
    for (const auto& [x, y] : pairs) {
      cpu_stop_check(1);
      pr_steps(x, y, false, false, nb);
      if (count) count_steps(1, nb);
      PAcc acc;
      if (cmax >= 1) {
        if (E.end_ok(x, y)) acc.add(kcl(x, y, 1), E.at(x, y));
        for (const ChainNbr& t : nb) {
          const PProf& q = V.find(key2(t.x, t.y))->second;
          for (size_t k = 0; k < q.k.size(); ++k)
            if ((q.k[k] & 0xffu) < cmax) acc.add(q.k[k] + 1, t.w + q.w[k]);
        }
      }
      PProf p;
      long long groups = 0;
      const long long d = pr_finish_layers(acc, p, x, groups);
      pr_vtrim(p, x, true);
      if (count) {
        pr_cut_drop += d;
        pr_V_entries += groups;
        pr_V_layer_entries += static_cast<long long>(p.k.size());
        ++pr_V_pairs;
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
      for (const ChainNbr& t : nb) {
        const PProf& q = V.find(key2(t.x, t.y))->second;
        for (size_t k = 0; k < q.k.size(); ++k)
          acc.add(key2(kcl_xs(q.k[k]), kcl_ys(q.k[k])), t.w + q.w[k]);
      }
      PProf p;
      const long long d = pr_finish_fut(acc, p, x);
      pr_vtrim(p, x);
      if (count) pr_cut_drop += d;
      const double base = (E.at(x, y) + E.w_stem) + pair_bonus;
      for (size_t k = 0; k < p.k.size(); ++k)
        ats.push_back({k_hi(p.k[k]), x, y, k_lo(p.k[k]), base + p.w[k]});
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
