// Pass-3 ranking scores during the scan: arms, lambda, K_L lookup (lK), the S score and
// lp_score (linear P3).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  void l2r_init_scores() {
    if (l2r_beam <= 0) return;
    SecTimer t(l2r_isec);
    pr_build_static();
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
    if (l2r_inc) lp_init_scores();
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
  double lK(int r, int qp) {
    if (r < 0 || qp + 1 > n - 1) return NEG_INF;
    return arm_phi(r, qp, pr_aq_sep ? stem_->Qzbar(0, qp) : 0.0);
  }

  double scoreS_k(int r, int j) {
    ++rfv_look_s;
    const double k = lK(r, j);
    if (k < HALF_NEG) {
      ++rfv_miss_s;
      const int mk = k_map_kind(r, j);
      note_kind(vm_i, vm_ii, vm_iii, ST_S, mk);
      note_kind(km_i, km_ii, km_iii, ST_S, mk);
    }
    return k;
  }
  double l2r_scoreS(uint64_t key, double f, int j) {
    if (f < HALF_NEG) return NEG_INF;
    const int a = k3a(key), xs = k3b(key), xp = k3c(key);
    const int r = larm(xs, xp);
    if (r < 0) ++l2r_noarm[PT_S];
    return f + scoreS_k(r, j) + prefix_log(a);
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
  //   Aopen = log sum_q' H[q'] / Pz(q'); ArmOpen(xs, j) =
  //   log sum_{V arms xs..x', x' >= j} W(xs, x') / Pz(x'), W = sum_y' Pz(y'-1) E_in V (prYw);
  //   Bopen(j) = log sum_{xs > j} Pz(xs-1) ArmOpen(xs, xs).
  // A score of -inf (no floor) is counted in lp_miss[table][reason]: 0 arm not in V, 1 H(a2, j)
  // empty, 2 no kept q' with K_q'(arm) > 0 (every y' <= q'), 3 ArmOpen / Bopen empty.
  int lp_hk = 16;
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
  bool lp_lin_ok() const { return lp_lin && pr_aq_sep && pr_pspan <= 0; }
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
    for (size_t k = m; k-- > 0;) {
      const int c1 = h.c1[k];
      const double base = lpz(c1 - 1) + (pair_bonus + E.w_stem + E.at(a2, c1));
      h.rt[k] = &rec_root(a2, c1, true);
      std::vector<std::pair<int, double>> add = lamL(*h.rt[k]);
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
      for (const auto& [qp, hv] : cum) ao = log_add(ao, hv - lpz(qp));
      h.aopen[k] = ao;
      lp_h_ents += static_cast<long long>(cum.size());
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
  double lp_score(int pt, uint64_t key, double f, int j) {
    if (f < HALF_NEG) return NEG_INF;
    const int t = pt - PT_P1, a2 = k3a(key), xs = k3b(key);
    int arm = -1;
    if (pt == PT_P3) {
      arm = larm(xs, k3c(key));
      if (arm < 0) {
        ++lp_miss[t][0];
        return NEG_INF;
      }
    }
    const LH& h = lh_get(a2);
    const int k = lh_at(h, j);
    if (k < 0 || h.v[k].empty()) {
      ++lp_miss[t][1];
      const int hk = h_map_kind(a2, j);
      note_kind(hm_i, hm_ii, hm_iii, t, hk);
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
      if (s < HALF_NEG) {
        ++lp_miss[t][2];
        const int qp = hl.empty() ? j : hl.front().first;
        const int mk = k_map_kind(arm, qp);
        note_kind(km_i, km_ii, km_iii, t, mk);
        return NEG_INF;
      }
      if (lp_hex) return f + s;
      return f - lpz(j) + s;
    }
    const double o = pt == PT_P2 ? l_armopen(xs, j) : lBopen[j];
    if (o < HALF_NEG) {
      ++lp_miss[t][3];
      return NEG_INF;
    }
    return f - lpz(j) + h.aopen[k] + o;
  }

  // Pre-pass suffix recursion (build_principled's prS) with qp[y] = (j, log Qp(y, j)) items as
  // one more paired element, closed as fill_C does: Qp e^beta1 Sbar[j+1].  qb: (j, log Qb(y, j))
  // items in place of the pre-pass nested Qb.
  using QpItems = std::vector<std::vector<std::pair<int, double>>>;
