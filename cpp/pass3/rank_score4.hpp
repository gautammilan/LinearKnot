// Pass-3 ranking score of a 4-index state (pr_score4) and the naive reference scores (prn_*,
// pr_naive4, --pr-check).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  // sigma of a state of Qg / Qgl / Qgr / Qgls / Qgrs (outer p, inner u, group v) at step q.
  double pr_score4(int bt, const Ent& e, int v, int q) {
    pr_build_static();
    const int p = e.p, u = e.u;
    pr_abl_p = p;
    pr_abl_q = q;
    if (g_grow_rec && (bt == BT_QG || bt == BT_QGLS || bt == BT_QGRS)) {
      const RN& r = bt == BT_QG ? rec_root(p, q, p == u)
                                : rec_key(bt == BT_QGLS, p, bt == BT_QGLS && v == q, q);
      const double s = pr_eval_rec(r, u, v, bt == BT_QG ? 0 : 1);
      if (g_rec_check) {
        const PProf& kp = bt == BT_QG ? pr_root_eval(p, q, p == u)
                                      : pr_key(bt == BT_QGLS, p, bt == BT_QGLS && v == q, q);
        rec_check(r, kp, s, pr_eval(kp, u, v));
      }
      return e.a + s;
    }
    if (bt == BT_QG) return e.a + pr_eval(pr_root_eval(p, q, p == u), u, v);
    if (bt == BT_QGLS) return e.a + pr_eval(pr_key(true, p, v == q, q), u, v);
    if (bt == BT_QGRS) return e.a + pr_eval(pr_key(false, p, false, q), u, v);
    if (!E.end_ok(p, q)) return NEG_INF;
    const StemEst& st = *stem_;
    double s = NEG_INF;
    if (bt == BT_QGL) {
      if (pr_arms_mode == 1) {
        const double se = pr_sep_ext(q, u, v, true);
        const PrTV& t = pr_T(q, u, v, true);
        double bm = NEG_INF;
        for (size_t k = 0; k < t.ys.size(); ++k) bm = log_add(bm, t.w[k] + st.Bml(p, t.ys[k]));
        s = log_add(prefix_log(p) + se, E.w_beta1m + bm);
      } else {
        const PrTV& t = pr_T(q, u, v, true);
        if (!t.ys.empty()) {
          double bm = NEG_INF;
          if (!(pr_abl & PA_BML))
            for (size_t k = 0; k < t.ys.size(); ++k) bm = log_add(bm, t.w[k] + st.Bml(p, t.ys[k]));
          s = log_add(prefix_log(p) + t.ext, E.w_beta1m + bm);
        }
      }
      auto ia = stA.find(key2(p, u + 1));
      if (ia != stA.end())
        for (int f : ia->second.fs) {
          if (f > v - 1) break;
          const double gr = pr_qz(f + 1, v - 1);
          if (gr < HALF_NEG) continue;
          const double c = pr_actx(p, u + 1, f, q);
          if (c > HALF_NEG) s = log_add(s, gr + c);
        }
    } else {
      s = pr_actx(p, u + 1, v - 1, q);
    }
    if (s < HALF_NEG) return NEG_INF;
    return e.a + E.at(p, q) + s;
  }

  static int pr_tab(int bt) {
    return bt == BT_QG ? 0 : bt == BT_QGL ? 1 : bt == BT_QGR ? 2 : bt == BT_QGLS ? 3 : bt == BT_QGRS ? 4 : 5;
  }

  // At a new step: per-step caches cleared, q'-keyed ones below q dropped.
  void pr_step(int q) {
    if (!prSepTouched.empty()) {
      ++pr_sep_steps;
      prSepTouched.clear();
    }
    prKls.clear();
    prKrs.clear();
    prRKls.clear();
    prRKrs.clear();
    prVkey[0].clear();
    prVkey[1].clear();
    prRT.clear();
    gm_step = 0;
    prAc.clear();
    pr_q = q;
    if (pr_gsel || pr_gcheck) {
      prGRs.clear();
      prUc.clear();
      for (auto it = prUL.begin(); it != prUL.end();) it = it->first < q ? prUL.erase(it) : std::next(it);
    }
    for (auto it = prTc.begin(); it != prTc.end();) it = it->first < q ? prTc.erase(it) : std::next(it);
    for (auto it = prAQ.begin(); it != prAQ.end();) it = it->first < q ? prAQ.erase(it) : std::next(it);
  }

  // --stem-pr-check: Out_hat of a ranked state with every DP03 path (growth, future stem,
  // partner, closing, context) as its own term, long double.
  struct PrPath {
    int x, y;
    long double w;
  };
  std::unordered_map<uint64_t, std::vector<PrPath>> prnG, prnV;

  static double lse(const std::vector<long double>& t) {
    if (t.empty()) return NEG_INF;
    long double m = *std::max_element(t.begin(), t.end()), s = 0.0L;
    for (long double x : t) s += std::exp(x - m);
    return static_cast<double>(m + std::log(s));
  }

  // Weight of the step x.y -> xo.yo, every step kind tested directly.
  double prn_step(int x, int y, int xo, int yo, bool inner, bool scan) const {
    const StemEst& st = *stem_;
    if (st.skip_kept) StemEst::lean_abort("prKept");
    if (!(xo < x && yo > y) || !st.prKept[ix(xo, yo)]) return NEG_INF;
    std::vector<long double> t;
    const int lc = std::min(pr_loop, cap);
    if (x - xo - 1 <= lc && yo - y - 1 <= lc) {
      const double s = E.stem_interior(xo, x, y, yo);
      if (s > HALF_NEG) t.push_back(static_cast<long double>(s) + pair_bonus);
    }
    if (pr_ml && (scan || pr_partner_ml) && E.end_ok(x, y) && E.end_ok(xo, yo)) {
      const long double br = static_cast<long double>(E.w_bp) + E.at(x, y);
      const long double cl = static_cast<long double>(E.w_ap) + E.w_bp + E.at(xo, yo) + pair_bonus;
      double lq = NEG_INF, rq = NEG_INF;
      if (xo + 1 <= x - 1 && x - (xo + 1) <= pr_ml_span)
        lq = scan ? qmp_at(x - 1).get(xo + 1) : st.PrMbar(xo + 1, x - 1);
      if (y + 1 <= yo - 1 && yo - y - 1 <= pr_ml_span) rq = st.PrMbar(y + 1, yo - 1);
      if (lq > HALF_NEG && yo - y - 1 <= (inner ? cap : pr_ml_span)) {
        long double w = lq + br + static_cast<long double>(E.w_cp) * (yo - y - 1) + cl;
        if (dz) w += E.run(y + 1, yo - 1);
        t.push_back(w);
      }
      if (ml_run_ok(x - xo - 1) && rq > HALF_NEG) {
        long double w = rq + br + static_cast<long double>(E.w_cp) * (x - xo - 1) + cl;
        if (dz) w += E.run(xo + 1, x - 1);
        t.push_back(w);
      }
      if (lq > HALF_NEG && rq > HALF_NEG) t.push_back(lq + rq + br + cl);
    }
    return lse(t);
  }

  void prn_paths(int x, int y, bool inner, bool scan, bool start_ok, long double w,
                 std::vector<PrPath>& out) const {
    if (start_ok && E.end_ok(x, y)) out.push_back({x, y, w});
    for (int xo = 0; xo < x; ++xo)
      for (int yo : stem_->prKeptRight[xo]) {
        if (yo <= y) continue;
        const double s = prn_step(x, y, xo, yo, inner, scan);
        if (s > HALF_NEG) prn_paths(xo, yo, false, scan, true, w + s, out);
      }
  }

  const std::vector<PrPath>& prn_G(int x, int y, bool inner) {
    const uint64_t key = key3(x, y, inner ? 1 : 0);
    auto it = prnG.find(key);
    if (it != prnG.end()) return it->second;
    std::vector<PrPath> v;
    prn_paths(x, y, inner, true, true, 0.0L, v);
    return prnG.emplace(key, std::move(v)).first->second;
  }

  const std::vector<PrPath>& prn_V(int x, int y) {
    const uint64_t key = key2(x, y);
    auto it = prnV.find(key);
    if (it != prnV.end()) return it->second;
    std::vector<PrPath> v;
    prn_paths(x, y, true, false, single_stems, 0.0L, v);
    return prnV.emplace(key, std::move(v)).first->second;
  }

  long double prn_ctx(int i, int j) const {
    const StemEst& st = *stem_;
    std::vector<long double> t;
    if (st.prS[j + 1] > HALF_NEG)
      t.push_back(static_cast<long double>(E.w_beta1) + prefix_log(i) + st.prS[j + 1]);
    if (st.Bml(i, j) > HALF_NEG) t.push_back(static_cast<long double>(E.w_beta1m) + st.Bml(i, j));
    return lse(t);
  }

  // Future R-stem terms: for the L-stem's outer pair pp.qp and inner gap (u, v), every kept
  // inner pair x'.y' (x' <= v-1, y' > qp) and stem path to xs.ys (xs >= u+1; xs == u+1 only
  // if exact_xs), with base added.
  void prn_future(int pp, int qp, int u, int v, bool exact_xs, long double base,
                  std::vector<long double>& t) {
    const StemEst& st = *stem_;
    for (int xp = u + 1; xp <= v - 1; ++xp)
      for (int yp : st.prKeptRight[xp]) {
        if (yp <= qp || !E.end_ok(xp, yp)) continue;
        const double gr = qz_log(xp + 1, v - 1), zb = st.Qzbar(qp + 1, yp - 1);
        if (gr < HALF_NEG || zb < HALF_NEG) continue;
        for (const PrPath& h : prn_V(xp, yp)) {
          if (h.x < u + 1 || (exact_xs && h.x != u + 1)) continue;
          const double gl = exact_xs ? 0.0 : qz_log(u + 1, h.x - 1);
          if (gl < HALF_NEG) continue;
          const long double c = prn_ctx(pp, h.y);
          if (c < HALF_NEG) continue;
          t.push_back(base + gl + gr + zb + E.w_stem + pair_bonus + E.at(xp, yp) + h.w +
                      E.at(h.x, h.y) + c);
        }
      }
  }

  // Surviving L-stems Qgl[i][pp-1][e][f] (e in [e_lo, e_hi], f <= v-1) closing a pseudoknot
  // ending at qp, read from the Qgl tables, with base added.
  void prn_partner(int pp, int qp, int u, int v, bool exact_e, long double base,
                   std::vector<long double>& t) {
    for (int f = u + 1; f <= v - 1; ++f)
      for (const Grp& g : Qgl[f].grps) {
        const int e = g.v;
        if (e < u + 1 || e > v - 1 || (exact_e && e != u + 1)) continue;
        const double gl = exact_e ? 0.0 : qz_log(u + 1, e - 1), gr = qz_log(f + 1, v - 1);
        if (gl < HALF_NEG || gr < HALF_NEG) continue;
        for (const Ent& s : g.ents) {
          if (s.u != pp - 1 || !E.end_ok(s.p, f)) continue;
          const long double c = prn_ctx(s.p, qp);
          if (c < HALF_NEG) continue;
          t.push_back(base + gl + s.a + E.at(s.p, f) + gr + c);
        }
      }
  }

  double prn_outqg(int x, int y, bool inner, int u, int v, long double w0) {
    std::vector<long double> t;
    if (!E.end_ok(u, v)) return NEG_INF;
    const long double own = static_cast<long double>(E.w_stem) + E.at(u, v);
    for (const PrPath& g : prn_G(x, y, inner)) {
      const long double base = w0 + g.w + E.at(g.x, g.y) + own;
      prn_future(g.x, g.y, u, v, false, base, t);
      prn_partner(g.x, g.y, u, v, false, base, t);
    }
    return lse(t);
  }

  double pr_naive4(int bt, const Ent& e, int v, int q) {
    const int p = e.p, u = e.u;
    if (bt == BT_QG) return prn_outqg(p, q, p == u, u, v, e.a);
    if (bt == BT_QGLS || bt == BT_QGRS) {
      std::vector<Closing> cl;
      if (bt == BT_QGLS) pr_ls_closings(p, v == q, q, cl);
      else closings_rs(p, q, cl);
      std::vector<long double> t;
      for (const Closing& c : cl) {
        const double o = prn_outqg(c.x, c.y, false, u, v, static_cast<long double>(e.a) + c.w);
        if (o > HALF_NEG) t.push_back(o);
      }
      return lse(t);
    }
    if (!E.end_ok(p, q)) return NEG_INF;
    std::vector<long double> t;
    const long double base = static_cast<long double>(e.a) + E.at(p, q);
    if (bt == BT_QGL) {
      prn_future(p, q, u, v, true, base, t);
      prn_partner(p, q, u, v, true, base, t);
    } else {
      for (const Grp& g : Qgl[v - 1].grps) {
        if (g.v != u + 1) continue;
        for (const Ent& s : g.ents) {
          if (s.u != p - 1 || !E.end_ok(s.p, v - 1)) continue;
          const long double c = prn_ctx(s.p, q);
          if (c > HALF_NEG) t.push_back(base + s.a + E.at(s.p, v - 1) + c);
        }
      }
    }
    return lse(t);
  }

  void prc_cmp(int t, double s, double ref) {
    ++prc_states[t];
    if ((ref > HALF_NEG) != (s > HALF_NEG)) ++prc_mis[t];
    else if (ref > HALF_NEG) prc_rel[t] = std::max(prc_rel[t], std::fabs(std::expm1(s - ref)));
  }

