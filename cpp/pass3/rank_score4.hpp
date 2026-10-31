// Pass-3 ranking score of a 4-index state (pr_score4).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  // sigma of a state of Qg / Qgl / Qgr / Qgls / Qgrs (outer p, inner u, group v) at step q.
  double pr_score4(int bt, const Ent& e, int v, int q) {
    pr_build_static();
    const int p = e.p, u = e.u;
    pr_abl_p = p;
    pr_abl_q = q;
    if ((bt == BT_QG || bt == BT_QGLS || bt == BT_QGRS)) {
      const RN& r = bt == BT_QG ? rec_root(p, q, p == u)
                                : rec_key(bt == BT_QGLS, p, bt == BT_QGLS && v == q, q);
      const double s = pr_eval_rec(r, u, v, bt == BT_QG ? 0 : 1);
      return e.a + s;
    }
    if (bt == BT_QG) return e.a + pr_eval(pr_root_eval(p, q, p == u), u, v);
    if (bt == BT_QGLS) return e.a + pr_eval(pr_key(true, p, v == q, q), u, v);
    if (bt == BT_QGRS) return e.a + pr_eval(pr_key(false, p, false, q), u, v);
    if (!E.end_ok(p, q)) return NEG_INF;
    const StemEst& st = *stem_;
    double s = NEG_INF;
    if (bt == BT_QGL) {
      {
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
    if (pr_gsel) {
      prGRs.clear();
      prUc.clear();
      for (auto it = prUL.begin(); it != prUL.end();) it = it->first < q ? prUL.erase(it) : std::next(it);
    }
    for (auto it = prTc.begin(); it != prTc.end();) it = it->first < q ? prTc.erase(it) : std::next(it);
    for (auto it = prAQ.begin(); it != prAQ.end();) it = it->first < q ? prAQ.erase(it) : std::next(it);
  }













