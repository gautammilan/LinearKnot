// Pass-3 ranking scores during the scan: arms, lambda, K_L lookups (lK, lKx), the S score and
// lp_score (linear P3).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  void l2r_init_scores() {
    if (l2r_beam <= 0) return;
    SecTimer t(l2r_isec);
    pr_build_static();
    std::unique_ptr<SecTimer> rtm(lpf ? new SecTimer(vi_rest_s) : nullptr);
    pr_sf_build();
    for (int r = 0; r < prA; ++r) lArm.emplace(key2(prArmXs[r], prArmXp[r]), r);
    lYwS.assign(prYy.size(), NEG_INF);
    for (int r = 0; r < prA; ++r) {
      double acc = NEG_INF;
      for (int k = prArmY0[r + 1] - 1; k >= prArmY0[r]; --k) {
        acc = log_add(acc, prYw[k]);
        lYwS[k] = acc;
      }
    }
    if (sp_cube && sp_static) l2r_spawn_cols();
    if (l2r_inc) lp_init_scores();
  }
  // --g-mode eager: G_L for every pair in P, before the scan, prefix frozen at C_pre.
  // y descending so an enclosing pair (larger right end) is built first.
  void growth_eager() {
    if (!p2_g_eager || !stem_ || l2r_beam <= 0) return;
    const std::vector<std::vector<int>>& KR = pr_pair_list(true);
    std::vector<std::pair<int, int>> pairs;
    for (int x = 0; x < n; ++x)
      for (int y : KR[x]) pairs.emplace_back(x, y);
    std::sort(pairs.begin(), pairs.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
      if (a.second != b.second) return a.second > b.second;
      return a.first > b.first;
    });
    for (const auto& [x, y] : pairs) {
      rec_root(x, y, false);
      rec_root(x, y, true);
    }
  }
  // Posterior mass of the spawned Rs states, summed over the steps (after the outside pass).
  double l2r_rs_mass() {
    outside();
    double m = 0.0;
    for (int j = 0; j < n; ++j)
      for (const Grp& g : LA.Gs[j].grps)
        for (const Ent& e : g.ents)
          if (e.a > HALF_NEG && e.b > HALF_NEG) m += std::exp(e.a + e.b - logZ);
    return m;
  }
  int l2r_arms() const { return prA; }
  size_t l2r_arm_ents() const { return pr_vxs ? prYy.size() : prEy.size(); }
  size_t l2r_arm_groups() const { return prYy.size(); }
  int larm(int xs, int xp) const {
    auto it = lArm.find(key2(xs, xp));
    return it == lArm.end() ? -1 : it->second;
  }
  // --l2r-arm-fb on: an arm xs..x' not in V scores as the nearest V arm xs''..x', xs < xs'' <= x',
  // times per extra outer pair x the best stack x.y onto (x+1).(y-1) over x's partners (with the
  // terminal penalties); with no such arm, or a score of -inf, the floor kFbFloor (finite).
  bool l2r_arm_fb = false;
  static constexpr double kFbFloor = -1e4;
  long long l2r_fb_hit = 0, l2r_fb_floor = 0;
  std::unordered_map<uint64_t, std::pair<int, double>> lFb;
  std::vector<double> lStk;
  double lstack_best(int x) {
    if (lStk.empty()) lStk.assign(n, std::nan(""));
    if (!std::isnan(lStk[x])) return lStk[x];
    double b = NEG_INF;
    for (int y : partners_right[x]) {
      if (x + 1 >= y - 1 || !E.pairs(x + 1, y - 1)) continue;
      const double s = E.stem_interior(x, x + 1, y - 1, y);
      if (s < HALF_NEG) continue;
      b = std::max(b, s + pair_bonus + E.at(x, y) - E.at(x + 1, y - 1));
    }
    return lStk[x] = b;
  }
  std::pair<int, double> larm_fb(int xs, int xp) {
    const int r = larm(xs, xp);
    if (r >= 0 || !l2r_arm_fb) return {r, 0.0};
    const uint64_t key = key2(xs, xp);
    auto it = lFb.find(key);
    if (it != lFb.end()) return it->second;
    std::pair<int, double> res{-1, 0.0};
    double extra = 0.0;
    for (int x = xs; x < xp; ++x) {
      const double s = lstack_best(x);
      if (s < HALF_NEG) break;
      extra += s;
      const int r2 = larm(x + 1, xp);
      if (r2 >= 0) {
        res = {r2, extra};
        break;
      }
    }
    return lFb.emplace(key, res).first->second;
  }
  // ---- D2 oracle substitutions (--l2r-sub a,c,d; (b) untrimmed V = --pr-v-eps 0)
  // a: K_q'(arm) = sum_{y' > q'} Qzbar(q'+1, y'-1) E_in V (not separable), cached per
  //    (q', arm) for the current step.  c: + rho(a'), the exact run's context correction of
  //    C(a'-1) (S and the U tables).  d: the exact outside of U / Us states of the exact run
  //    (U, Us; S through its B seeds), with rho divided out unless c.
  bool lsub_a = false, lsub_c = false, lsub_d = false;
  std::vector<double> lRho;
  std::unordered_map<uint64_t, double> lOU, lOUs;
  std::unordered_map<uint64_t, double> lKc;
  int lKc_q = -1;
  static uint64_t key4(int a, int b, int c, int d) {
    return (static_cast<uint64_t>(a) << 48) | (static_cast<uint64_t>(b) << 32) |
           (static_cast<uint64_t>(c) << 16) | static_cast<uint64_t>(d);
  }
  double lrho(int a) const { return lRho.empty() ? 0.0 : lRho[a]; }
  // --l2r-ctx p1: log lambda(a') from a finished pass 1 (empty: off).  Scores only: g_L's stop
  // term C(p-1) lambda(p), and the S and U-type scores times lambda(a').
  std::vector<double> lLam;
  double llam(int a) const { return lLam.empty() ? 0.0 : lLam[a]; }
  // lambda(a') = [sum_ys I1 O1 + tau W] / [sum_ys I1 E + tau W] over src's Qp items (a', ys),
  // I1 / O1 their inside / outside, E = e^beta1 C1(a'-1) Sbar_pk(ys+1) (C1: src's exterior
  // prefix), W = mean of I1 E at a'; 1 where src has no item at a'.  src.outside() is run here.
  void l2r_ctx_lambda(DP03& src, const std::vector<double>& spk, double tau,
                      std::vector<double>& lam, std::vector<int>& cnt) const {
    src.outside();
    std::vector<double> num(n, NEG_INF), den(n, NEG_INF);
    cnt.assign(n, 0);
    for (int ys = 0; ys < n; ++ys) {
      const OMap& at = src.Qp_at[ys];
      for (size_t r = 0; r < at.slots(); ++r) {
        const int a = at.keys[r];
        if (a < 0) continue;
        const double I = at.vals[r];
        num[a] = log_add(num[a], I + src.Qp_b[ix(a, ys)]);
        den[a] = log_add(den[a], I + (E.w_beta1 + src.prefix_log(a) + spk[ys + 1]));
        ++cnt[a];
      }
    }
    lam.assign(n, 0.0);
    const double lt = std::log(tau);
    for (int a = 0; a < n; ++a) {
      if (cnt[a] == 0 || den[a] < HALF_NEG) continue;
      const double tw = tau > 0.0 ? lt + den[a] - std::log(static_cast<double>(cnt[a])) : NEG_INF;
      lam[a] = log_add(num[a], tw) - log_add(den[a], tw);
    }
  }
  double lKx_sum(int r, int qp) const {
    const StemEst& st = *stem_;
    const int hi = prArmY0[r + 1];
    double acc = NEG_INF;
    for (int k = static_cast<int>(std::upper_bound(prYy.begin() + prArmY0[r], prYy.begin() + hi, qp) -
                                  prYy.begin());
         k < hi; ++k) {
      const double z = st.Qzbar(qp + 1, prYy[k] - 1);
      if (z > HALF_NEG && prYw[k] > HALF_NEG) acc = log_add(acc, z + prYw[k]);
    }
    return acc;
  }
  double lKx(int r, int qp) {
    if (lKc_q != cur_j) {
      lKc.clear();
      lKc_q = cur_j;
    }
    const uint64_t key = key2(qp, r);
    auto it = lKc.find(key);
    if (it != lKc.end()) return it->second;
    const double acc = lKx_sum(r, qp);
    lKc.emplace(key, acc);
    return acc;
  }
  // Persistent exact-sum memo for --k-exact-ub. Same key as --k-memo, value is lKx.
  double lKx_memo(int r, int qp) {
    ++lKx_calls;
    const uint64_t key = key2(qp, r);
    auto it = lKxM.find(key);
    if (it != lKxM.end()) {
      ++lKx_hits;
      return it->second;
    }
    const double v = lKx_sum(r, qp);
    lKxM.emplace(key, v);
    return v;
  }
  double lK(int r, int qp) {
    if (r < 0 || qp + 1 > n - 1) return NEG_INF;
    if (sx_k_exact) return lKx_memo(r, qp);
    if (lsub_a) return lKx(r, qp);
    if (!k_memo) return arm_phi(r, qp, pr_aq_sep ? stem_->Qzbar(0, qp) : 0.0);
    ++lKm_calls;
    const uint64_t key = key2(qp, r);
    auto it = lKm.find(key);
    if (it != lKm.end()) {
      ++lKm_hits;
      return it->second;
    }
    const double v = arm_phi(r, qp, pr_aq_sep ? stem_->Qzbar(0, qp) : 0.0);
    lKm.emplace(key, v);
    return v;
  }
  // --l2r-ctx on: the nested pre-pass's span context pr_ctx(c, ys) = e^beta1 C[c-1] Sbar[ys+1] +
  // e^beta1m Bml(c, ys) in place of the scores' exterior e^beta1 C[c-1] Sbar[ys+1], as one log
  // factor per state: kbar(c) = log sum_ys W(ys) e^(pr_ctx(c, ys) - C[c-1]) -
  // log sum_ys W(ys) e^(beta1 + Sbar[ys+1]), W the pre-pass V weight by ys of the state's arm
  // (S and R types; c = a' for S, the A pair's left end for R) or of every arm with the state's
  // xs (U types, c = a').  Cached per (c, arm) / (c, xs).
  bool l2r_ctx = false;
  std::vector<std::vector<std::pair<int, double>>> lWr, lWx;
  std::unordered_map<uint64_t, double> lKbr, lKbx;
  void l2r_ctx_init() {
    if (!l2r_ctx) return;
    if (pr_vxs) throw std::runtime_error("--pr-v-key xs: --l2r-ctx needs per-ys arm entries");
    auto add = [](std::vector<std::pair<int, double>>& v, int ys, double w) {
      for (auto& [y, s] : v)
        if (y == ys) {
          s = log_add(s, w);
          return;
        }
      v.emplace_back(ys, w);
    };
    lWr.assign(prA, {});
    lWx.assign(n, {});
    for (int r = 0; r < prA; ++r)
      for (int k = prArmE0[r]; k < prArmE0[r + 1]; ++k) {
        add(lWr[r], prEys[k], prEw[k]);
        add(lWx[prArmXs[r]], prEys[k], prEw[k]);
      }
  }
  double lkbar(const std::vector<std::pair<int, double>>& w, int c) const {
    const StemEst& st = *stem_;
    const double pc = prefix_log(c);
    double num = NEG_INF, den = NEG_INF;
    for (const auto& [ys, lw] : w) {
      if (ys <= c) continue;
      const double ex = E.w_beta1 + st.prS[ys + 1];
      den = log_add(den, lw + ex);
      num = log_add(num, lw + log_add(ex, E.w_beta1m + st.Bml(c, ys) - pc));
    }
    return den > HALF_NEG ? num - den : 0.0;
  }
  double lkb_arm(int c, int r) {
    if (!l2r_ctx || r < 0) return 0.0;
    const uint64_t key = key2(c, r);
    auto it = lKbr.find(key);
    return it != lKbr.end() ? it->second : lKbr.emplace(key, lkbar(lWr[r], c)).first->second;
  }
  double lkb_xs(int c, int xs) {
    if (!l2r_ctx) return 0.0;
    const uint64_t key = key2(c, xs);
    auto it = lKbx.find(key);
    return it != lKbx.end() ? it->second : lKbx.emplace(key, lkbar(lWx[xs], c)).first->second;
  }

  double scoreS_k(int r, int j) {
    ++rfv_look_s;
    const double k = lK(r, j);
    if (k < HALF_NEG) {
      ++rfv_miss_s;
      const int mk = k_map_kind(r, j);
      note_kind(vm_i, vm_ii, vm_iii, ST_S, mk);
      note_kind(km_i, km_ii, km_iii, ST_S, mk);
      if (mk == 2) miss_k(r, j, ST_S);
    }
    return k;
  }
  double l2r_scoreS(uint64_t key, double f, int j) {
    if (f < HALF_NEG) return NEG_INF;
    const int a = k3a(key), xs = k3b(key), xp = k3c(key);
    if (!lsub_d && l2r_arm_fb) {
      const auto [r, fx] = larm_fb(xs, xp);
      if (r >= 0 && r != larm(xs, xp)) ++l2r_fb_hit;
      const double k = scoreS_k(r, j);
      if (k < HALF_NEG) {
        ++l2r_fb_floor;
        return f + prefix_log(a) + kFbFloor;
      }
      return f + (k + fx) + prefix_log(a) + (lsub_c ? lrho(a) : lkb_arm(a, r)) + llam(a);
    }
    if (!lsub_d && !lLam.empty()) {
      const int r = larm(xs, xp);
      return f + scoreS_k(r, j) + prefix_log(a) + (lsub_c ? lrho(a) : lkb_arm(a, r)) + lLam[a];
    }
    if (!lsub_d) {
      const int r = larm(xs, xp);
      if (r < 0) ++l2r_noarm[PT_S];
      return f + scoreS_k(r, j) + prefix_log(a) + (lsub_c ? lrho(a) : lkb_arm(a, r));
    }
    const StemEst& st = *stem_;
    double acc = NEG_INF;
    for (int y : partners_right[xp]) {
      if (y <= j || !E.end_ok(xp, y)) continue;
      auto it = lOUs.find(key4(a, xs, xp, y));
      if (it == lOUs.end()) continue;
      const double fill = lsub_a ? st.Qzbar(j + 1, y - 1) : st.Qzbar(0, y - 1) - st.Qzbar(0, j);
      if (!(fill > HALF_NEG) || !(fill < 1e300)) continue;
      acc = log_add(acc, fill + (E.w_stem + E.at(xp, y) + pair_bonus) + it->second);
    }
    return acc > HALF_NEG ? f + acc - (lsub_c ? 0.0 : lrho(a)) : NEG_INF;
  }

  // ---- --l2r-spawn inc scores (log; values exact, scores only rank the states of one step j).
  // Pz(i) = Qzbar(0, i), the pre-pass prefix gap filler: Qz(i, j) ~ Pz(j) / Pz(i-1), an
  // overestimate, used only as a ratio whose 1 / Pz(j) is common to the states at j.
  //   H(a2, j)[q'] = sum_{c1 > j, a2.c1} Pz(c1-1) E_A(a2, c1) Lambda_L(a2, c1)[q'], Lambda_L the
  //   truncated L vector of rec_root(a2, c1, inner) (C[p'-1] folded in); per a2 a suffix list
  //   over c1, each suffix kept as its top lp_hk entries by H uL(q') (ntrUL; q' with uL = 0
  //   dropped).  Built at the scan's step a2 (the roots read C left of a2 only).
  //   P3: F - Pz(j) + log sum_q' H(a2, j)[q'] K_q'(arm) + ctx(a2, arm)
  //   P2: F - Pz(j) + Aopen(a2, j) + ArmOpen(xs, j) + ctx(a2, xs)
  //   P1: F - Pz(j) + Aopen(a2, j) + Bopen(j)
  //   Aopen = log sum_q' H[q'] / Pz(q') (--l2r-p-aopen q; raw: sum_q' H[q']); ArmOpen(xs, j) =
  //   log sum_{V arms xs..x', x' >= j} W(xs, x') / Pz(x'), W = sum_y' Pz(y'-1) E_in V (prYw);
  //   Bopen(j) = log sum_{xs > j} Pz(xs-1) ArmOpen(xs, xs).
  // A score of -inf (no floor) is counted in lp_miss[table][reason]: 0 arm not in V, 1 H(a2, j)
  // empty, 2 no kept q' with K_q'(arm) > 0 (every y' <= q'), 3 ArmOpen / Bopen empty.
  int lp_hk = 16;
  bool lp_aopen_raw = false;
  struct LH {
    std::vector<int> c1;                                 // ascending
    std::vector<std::vector<std::pair<int, double>>> v;  // suffix c1[k..], q' ascending
    std::vector<double> aopen;
    std::vector<const RN*> rt;                           // rec_root(a2, c1[k], true)
    std::vector<std::vector<double>> l;                  // --l2r-p3-lin: e^(v - lpz(q') - lm)
    std::vector<double> lm;
  };
  // --l2r-p3-lin on (separable arms, no span, no sub a): the P3 dot product in linear space,
  // H as e^(v - lpz(q') - m) with one m per list, K as e^(prSF - prSF[first of arm]).
  bool lp_lin = false;
  std::vector<double> lSFl, lSFm;
  bool lp_lin_ok() const { return lp_lin && pr_aq_sep && pr_pspan <= 0 && !lsub_a; }
  void lp_lin_of(const std::vector<std::pair<int, double>>& v, std::vector<double>& l, double& m) const {
    l.resize(v.size());
    m = NEG_INF;
    for (size_t i = 0; i < v.size(); ++i) {
      const int q = v[i].first;
      const double z = q + 1 <= n - 1 ? lpz(q) : NEG_INF;
      l[i] = z > HALF_NEG && v[i].second > HALF_NEG ? v[i].second - z : NEG_INF;
      m = std::max(m, l[i]);
    }
    for (double& x : l) x = x > HALF_NEG ? std::exp(x - m) : 0.0;
  }
  std::vector<LH> lH;
  std::vector<char> lHdone;
  std::vector<std::vector<std::pair<int, double>>> lArmOpen;  // per xs: (x', suffix), x' ascending
  std::vector<double> lBopen;
  long long lp_h_ents = 0, lp_h_roots = 0, lp_h_lists = 0;
  double lp_h_sec = 0.0;
  long long lp_miss[3][4] = {};
  double lpz(int i) const { return stem_->Qzbar(0, i); }
  std::vector<std::pair<int, double>> lamL(const RN& r) {
    if (g_ntr_K > 0) return r.T[0];
    if (g_tr_K > 0) return rec_trunc(r, 0).L;
    std::vector<double> acc(n, NEG_INF);
    for (int m = 0; m <= g_ml_L && m < rnL; ++m) {
      const RV& a = r.L[m];
      for (size_t i = 0; i < a.v.size(); ++i)
        if (a.v[i] > 0.0) acc[a.lo + i] = log_add(acc[a.lo + i], a.mx + std::log(a.v[i]));
    }
    std::vector<std::pair<int, double>> o;
    for (int q = 0; q < n; ++q)
      if (acc[q] > HALF_NEG) o.emplace_back(q, acc[q]);
    return o;
  }
  const LH& lh_get(int a2) {
    if (lHdone[a2]) return lH[a2];
    if (a2 > cur_j + 1) {
      std::fprintf(stderr, "[l2r] abort: H(%d) built at scan %d (needs C left of a2)\n", a2, cur_j);
      std::abort();
    }
    SecTimer t(lp_h_sec);
    ProfReg reg(*this, 3);
    lHdone[a2] = 1;
    ntr_init();
    LH& h = lH[a2];
    for (int c1 : partners_right[a2])
      if (c1 >= a2 + 2 && E.end_ok(a2, c1)) h.c1.push_back(c1);
    const size_t m = h.c1.size();
    h.v.resize(m);
    h.aopen.assign(m, NEG_INF);
    h.rt.assign(m, nullptr);
    std::vector<std::pair<int, double>> cum, mrg;
    std::vector<std::pair<double, int>> rk;
    const long long hb0 = lpf_help_built, hc0 = lpf_hg_clears, rr0 = rn_roots;
    const double hs0 = lpf_help_sec;
    for (size_t k = m; k-- > 0;) {
      const int c1 = h.c1[k];
      const double base = lpz(c1 - 1) + (pair_bonus + E.w_stem + E.at(a2, c1));
      std::chrono::steady_clock::time_point ht0;
      if (lpf) ht0 = std::chrono::steady_clock::now();
      h.rt[k] = &rec_root(a2, c1, true);
      std::vector<std::pair<int, double>> add = lamL(*h.rt[k]);
      if (lpf) {
        const auto ht1 = std::chrono::steady_clock::now();
        lpf_h_root_sec += std::chrono::duration<double>(ht1 - ht0).count();
        ht0 = ht1;
        lpf_h_lam += static_cast<long long>(add.size());
        lpf_h_merge += static_cast<long long>(cum.size() + add.size());
        if (cum.size() + add.size() > static_cast<size_t>(lp_hk)) {
          ++lpf_h_trunc;
          lpf_h_trunc_in += static_cast<long long>(cum.size() + add.size());
        }
      }
      ++lp_h_roots;
      mrg.clear();
      size_t x = 0, y = 0;
      while (x < cum.size() || y < add.size()) {
        if (y < add.size() && ntrUL[add[y].first] < HALF_NEG) {
          ++y;
          continue;
        }
        if (y == add.size() || (x < cum.size() && cum[x].first < add[y].first)) {
          mrg.push_back(cum[x++]);
        } else if (x == cum.size() || add[y].first < cum[x].first) {
          mrg.emplace_back(add[y].first, base + add[y].second);
          ++y;
        } else {
          mrg.emplace_back(cum[x].first, log_add(cum[x].second, base + add[y].second));
          ++x;
          ++y;
        }
      }
      if (mrg.size() > static_cast<size_t>(lp_hk)) {
        rk.clear();
        for (size_t i = 0; i < mrg.size(); ++i) rk.emplace_back(mrg[i].second + ntrUL[mrg[i].first], static_cast<int>(i));
        std::nth_element(rk.begin(), rk.begin() + lp_hk, rk.end(),
                         [](const std::pair<double, int>& p, const std::pair<double, int>& q) {
                           return p.first > q.first || (p.first == q.first && p.second < q.second);
                         });
        std::vector<char> keep(mrg.size(), 0);
        for (int i = 0; i < lp_hk; ++i) keep[rk[i].second] = 1;
        size_t w = 0;
        for (size_t i = 0; i < mrg.size(); ++i)
          if (keep[i]) mrg[w++] = mrg[i];
        mrg.resize(w);
      }
      cum.swap(mrg);
      h.v[k] = cum;
      if (lp_lin_ok()) {
        h.l.resize(m);
        h.lm.resize(m);
        lp_lin_of(cum, h.l[k], h.lm[k]);
      }
      double ao = NEG_INF;
      for (const auto& [qp, hv] : cum) ao = log_add(ao, lp_aopen_raw ? hv : hv - lpz(qp));
      h.aopen[k] = ao;
      lp_h_ents += static_cast<long long>(cum.size());
      if (lpf) lpf_h_merge_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - ht0).count();
    }
    if (lpf) {
      lpf_h_newroots += rn_roots - rr0;
      lpf_h_helpers += lpf_help_built - hb0;
      lpf_h_clears += lpf_hg_clears - hc0;
      lpf_h_help_sec += lpf_help_sec - hs0;
    }
    ++lp_h_lists;
    return h;
  }
  // index of the first c1 > j in H(a2), or -1
  int lh_at(const LH& h, int j) const {
    const auto it = std::upper_bound(h.c1.begin(), h.c1.end(), j);
    return it == h.c1.end() ? -1 : static_cast<int>(it - h.c1.begin());
  }
  // --l2r-p3-hex on: the P3 score's continuation with the pre-pass filler in place of the prefix
  // ratio, H_ex(a2, j)[q'] = sum_{c1 > j} Qzbar(j+1, c1-1) E_A(a2, c1) Lambda_L(a2, c1)[q'] (the
  // backward filler recursion's value), top lp_hk entries as H ranks them; built per (a2, j) at
  // step j, kept for the current j only.
  bool lp_hex = false;
  int lhx_j = -1;
  struct LX {
    std::vector<std::pair<int, double>> v;
    std::vector<double> l;
    double m = NEG_INF;
  };
  std::unordered_map<int, LX> lHx;
  std::vector<double> lhxA;
  std::vector<int> lhxT;
  long long lp_hx_builds = 0, lp_hx_terms = 0, lp_hx_ents = 0;
  double lp_hx_sec = 0.0;
  const LX& lh_ex(const LH& h, int a2, int j) {
    if (j != lhx_j) {
      lHx.clear();
      lhx_j = j;
    }
    auto it = lHx.find(a2);
    if (it != lHx.end()) return it->second;
    SecTimer tm(lp_hx_sec);
    if (lhxA.empty()) lhxA.assign(n, NEG_INF);
    lhxT.clear();
    const int k0 = lh_at(h, j);
    for (int k = k0 < 0 ? static_cast<int>(h.c1.size()) : k0; k < static_cast<int>(h.c1.size()); ++k) {
      const int c1 = h.c1[k];
      const double z = stem_->Qzbar(j + 1, c1 - 1);
      if (z < HALF_NEG) continue;
      const double base = z + (pair_bonus + E.w_stem + E.at(a2, c1));
      std::vector<std::pair<int, double>> own;
      const std::vector<std::pair<int, double>>& lam = g_ntr_K > 0 ? h.rt[k]->T[0] : (own = lamL(*h.rt[k]));
      for (const auto& [q, g] : lam) {
        if (ntrUL[q] < HALF_NEG) continue;
        if (lhxA[q] < HALF_NEG) lhxT.push_back(q);
        lhxA[q] = log_add(lhxA[q], base + g);
        ++lp_hx_terms;
      }
    }
    std::sort(lhxT.begin(), lhxT.end());
    std::vector<std::pair<int, double>> o;
    o.reserve(lhxT.size());
    for (int q : lhxT) {
      o.emplace_back(q, lhxA[q]);
      lhxA[q] = NEG_INF;
    }
    if (o.size() > static_cast<size_t>(lp_hk)) {
      std::vector<std::pair<double, int>> rk;
      for (size_t i = 0; i < o.size(); ++i) rk.emplace_back(o[i].second + ntrUL[o[i].first], static_cast<int>(i));
      std::nth_element(rk.begin(), rk.begin() + lp_hk, rk.end(),
                       [](const std::pair<double, int>& p, const std::pair<double, int>& q) {
                         return p.first > q.first || (p.first == q.first && p.second < q.second);
                       });
      std::vector<char> keep(o.size(), 0);
      for (int i = 0; i < lp_hk; ++i) keep[rk[i].second] = 1;
      size_t w = 0;
      for (size_t i = 0; i < o.size(); ++i)
        if (keep[i]) o[w++] = o[i];
      o.resize(w);
    }
    ++lp_hx_builds;
    lp_hx_ents += static_cast<long long>(o.size());
    LX x;
    if (lp_lin_ok()) lp_lin_of(o, x.l, x.m);
    x.v = std::move(o);
    return lHx.emplace(a2, std::move(x)).first->second;
  }
  double l_armopen(int xs, int j) const {
    if (xs >= n) return NEG_INF;
    const auto& v = lArmOpen[xs];
    const auto it = std::lower_bound(v.begin(), v.end(), j,
                                     [](const std::pair<int, double>& e, int x) { return e.first < x; });
    return it == v.end() ? NEG_INF : it->second;
  }
  void lp_init_scores() {
    if (lp_lin_ok()) {
      lSFl.assign(prSF.size(), 0.0);
      lSFm.assign(prA, NEG_INF);
      for (int r = 0; r < prA; ++r) {
        if (prArmY0[r] >= prArmY0[r + 1]) continue;
        const double m = prSF[prArmY0[r]];
        if (m < HALF_NEG) continue;
        lSFm[r] = m;
        for (int k = prArmY0[r]; k < prArmY0[r + 1]; ++k)
          if (prSF[k] > HALF_NEG) lSFl[k] = std::exp(prSF[k] - m);
      }
    }
    lH.assign(n, LH{});
    lHdone.assign(n, 0);
    lArmOpen.assign(n, {});
    for (int xs = 0; xs < n; ++xs) {
      auto& v = lArmOpen[xs];
      for (int r = prXsFirst[xs]; r < prXsFirst[xs + 1]; ++r) {
        double w = NEG_INF;
        for (int k = prArmY0[r]; k < prArmY0[r + 1]; ++k) {
          const double z = lpz(prYy[k] - 1);
          if (z > HALF_NEG) w = log_add(w, z + prYw[k]);
        }
        if (w > HALF_NEG) v.emplace_back(prArmXp[r], w - lpz(prArmXp[r]));
      }
      for (size_t k = v.size(); k-- > 1;) v[k - 1].second = log_add(v[k - 1].second, v[k].second);
    }
    lBopen.assign(n + 1, NEG_INF);
    for (int j = n - 2; j >= 0; --j) lBopen[j] = log_add(lBopen[j + 1], lpz(j) + l_armopen(j + 1, j + 1));
  }
  std::chrono::steady_clock::time_point lpf_t0;
  void lpf_lap(double& acc) {
    const auto t1 = std::chrono::steady_clock::now();
    acc += std::chrono::duration<double>(t1 - lpf_t0).count();
    lpf_t0 = t1;
  }
  double lp_score(int pt, uint64_t key, double f, int j) {
    if (f < HALF_NEG) return NEG_INF;
    const int t = pt - PT_P1, a2 = k3a(key), xs = k3b(key);
    if (lpf) {
      ++lpf_sc_calls[t];
      lpf_t0 = std::chrono::steady_clock::now();
    }
    int arm = -1;
    if (pt == PT_P3) {
      arm = larm(xs, k3c(key));
      if (lpf) {
        ++lpf_sc_osteps[t];
        lpf_lap(lpf_sc_osec[t]);
      }
      if (arm < 0) {
        ++lp_miss[t][0];
        return NEG_INF;
      }
    }
    const LH& h = lh_get(a2);
    const int k = lh_at(h, j);
    if (lpf) {
      lpf_sc_hsteps[t] += 1 + lpf_log2(h.c1.size());
      lpf_lap(lpf_sc_hsec[t]);
    }
    if (k < 0 || h.v[k].empty()) {
      ++lp_miss[t][1];
      const int hk = h_map_kind(a2, j);
      note_kind(hm_i, hm_ii, hm_iii, t, hk);
      if (hk == 2 && (t == ST_P3 || t == ST_S)) miss_take(2, t, a2, j, -1);
      return NEG_INF;
    }
    if (pt == PT_P3) {
      const LX* hx = lp_hex ? &lh_ex(h, a2, j) : nullptr;
      const std::vector<std::pair<int, double>>& hl = hx ? hx->v : h.v[k];
      double s = NEG_INF;
      if (lp_lin_ok()) {
        const std::vector<double>& ll = hx ? hx->l : h.l[k];
        const int hi = prArmY0[arm + 1];
        const int* yy = prYy.data();
        int kk = prArmY0[arm];
        double a = 0.0;
        for (size_t i = 0; i < hl.size(); ++i) {
          if (ll[i] == 0.0) continue;
          kk = static_cast<int>(std::upper_bound(yy + kk, yy + hi, hl[i].first) - yy);
          if (kk >= hi) break;
          a += ll[i] * lSFl[kk];
        }
        if (a > 0.0) s = (hx ? hx->m : h.lm[k]) + lSFm[arm] + std::log(a);
      } else {
        for (const auto& [qp, hv] : hl) {
          const double kq = lK(arm, qp);
          if (kq > HALF_NEG) s = log_add(s, hv + kq);
        }
      }
      if (lpf) {
        const long long m = static_cast<long long>(hl.size());
        lpf_sc_dot[t] += m;
        lpf_sc_kcalls[t] += m;
        lpf_sc_ksteps[t] += m * lpf_log2(static_cast<size_t>(prArmY0[arm + 1] - prArmY0[arm]));
        lpf_lap(lpf_sc_ksec[t]);
      }
      if (s < HALF_NEG) {
        ++lp_miss[t][2];
        const int qp = hl.empty() ? j : hl.front().first;
        const int mk = k_map_kind(arm, qp);
        note_kind(km_i, km_ii, km_iii, t, mk);
        if (mk == 2) miss_k(arm, qp, t);
        return NEG_INF;
      }
      if (lp_hex) return f + s + lkb_arm(a2, arm);
      return f - lpz(j) + s + lkb_arm(a2, arm);
    }
    if (lpf) {
      lpf_sc_osteps[t] += pt == PT_P2 && xs < n ? lpf_log2(lArmOpen[xs].size()) : 1;
      lpf_sc_dot[t] += 1;
    }
    const double o = pt == PT_P2 ? l_armopen(xs, j) : lBopen[j];
    if (lpf) lpf_lap(lpf_sc_osec[t]);
    if (o < HALF_NEG) {
      ++lp_miss[t][3];
      return NEG_INF;
    }
    return f - lpz(j) + h.aopen[k] + o + (pt == PT_P2 ? lkb_xs(a2, xs) : 0.0);
  }

  // Pre-pass suffix recursion (build_principled's prS) with qp[y] = (j, log Qp(y, j)) items as
  // one more paired element, closed as fill_C does: Qp e^beta1 Sbar[j+1].  qb: (j, log Qb(y, j))
  // items in place of the pre-pass nested Qb.
  using QpItems = std::vector<std::vector<std::pair<int, double>>>;
