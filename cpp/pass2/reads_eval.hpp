// Pass-2 reads Ghat for pass 3: rec_root / rec_key, Lambda dumps (DUMP_P2_ALL), truncation
// (rec_trunc, tau), pr_eval_rec / pr_eval_sep / pr_eval, and their checks.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  const RN& rec_root(int x, int y, bool inner) {
    const uint64_t key = key3(x, y, inner ? 1 : 0);
    auto it = prRR.find(key);
    const bool dump_focus = std::getenv("DUMP_ROOT") &&
                            ((x == 4 && y == 35) || (x == 9 && y == 37) || (x == 9 && y == 41));
    auto dump_t0 = [&](const RN& placed) {
      if (!dump_focus) return;
      std::fprintf(stderr,
                   "DUMP_ROOT_CREATE\tx\t%d\ty\t%d\tinner\t%d\tp2_sweep\t%d\tK\t%d\tntr\t%d\n", x, y,
                   inner ? 1 : 0, p2_sweep ? 1 : 0, g_ntr_K, ntrUL.empty() ? 0 : 1);
      if (!placed.T[0].empty()) {
        for (const auto& e : placed.T[0])
          std::fprintf(stderr, "DUMP_ROOT_T0\tx\t%d\ty\t%d\tq\t%d\tg\t%.17g\n", x, y, e.first,
                       e.second);
        std::fprintf(stderr, "DUMP_ROOT_T0META\tx\t%d\ty\t%d\tn\t%zu\tsrc\tT0\n", x, y,
                     placed.T[0].size());
      } else {
        auto lam = lamL(placed);
        for (const auto& e : lam)
          std::fprintf(stderr, "DUMP_ROOT_T0\tx\t%d\ty\t%d\tq\t%d\tg\t%.17g\n", x, y, e.first,
                       e.second);
        std::fprintf(stderr, "DUMP_ROOT_T0META\tx\t%d\ty\t%d\tn\t%zu\tsrc\tlamL\n", x, y, lam.size());
      }
    };
    if (it != prRR.end()) {
      if (dump_focus)
        std::fprintf(stderr, "DUMP_ROOT_HIT\tx\t%d\ty\t%d\tinner\t%d\tp2_sweep\t%d\n", x, y,
                     inner ? 1 : 0, p2_sweep ? 1 : 0);
      return it->second;
    }
    if (p2_sweep) {
      PrGrowTimer gt(*this);
      RN r = sw_as_rn(x, y, inner);
      ++rn_roots;
      auto& placed = prRR.emplace(key, std::move(r)).first->second;
      dump_t0(placed);
      return placed;
    }
    PrGrowTimer gt(*this);
    std::vector<ChainNbr> nb;
    pr_steps(x, y, inner, !geq_pre, nb);
    if (nb.empty() && p2_g_fb && p2_g_mode > 0 && !g_in_set(x, y)) g_one_step(x, y, nb);
    count_steps(0, nb);
    for (const ChainNbr& t : nb) rec_gbuild(t.x, t.y);
    RN r;
    rec_node(x, y, inner, nb, r);
    ++rn_roots;
    rec_stat(1, r);
    gm_roots += rbytes(r);
    gm_update();
    auto& placed = prRR.emplace(key, std::move(r)).first->second;
    dump_t0(placed);
    return placed;
  }
  const RN& rec_key(bool ls, int p, bool seed, int q) {
    const auto tick = []() { return std::chrono::steady_clock::now(); };
    auto t0 = tick();
    auto charge = [&](std::chrono::steady_clock::time_point t1) {
      sx_rkey_sec += std::chrono::duration<double>(t1 - t0).count();
    };
    std::unordered_map<int, RN>& m = ls ? prRKls : prRKrs;
    const int key = ls ? 2 * p + (seed ? 1 : 0) : p;
    auto it = m.find(key);
    if (it != m.end()) {
      charge(tick());
      return it->second;
    }
    PrGrowTimer gt(*this);
    std::vector<Closing> cl;
    if (ls) pr_ls_closings(p, seed, q, cl);
    else closings_rs(p, q, cl);
    charge(tick());
    std::vector<RTerm> t;
    for (const Closing& c : cl) t.push_back({&rec_root(c.x, c.y, false), c.w, 0});
    t0 = tick();
    rTerms = std::move(t);
    racc_begin();
    racc_refs();
    racc_add_terms();
    RN r;
    racc_finish(r);
    ++rn_keys;
    rec_stat(2, r);
    gm_step += rbytes(r);
    gm_update();
    auto placed = m.emplace(key, std::move(r)).first;
    charge(tick());
    return placed->second;
  }

  // DUMP_P2_ALL=prefix: write Lambda_L/R (roots+nodes+keys), V-key lists, Qz_pre for offline
  // sweep-vs-lazy comparison (Ther |d| gap diagnosis).
  void dump_p2_all(const char* prefix) const {
    if (!prefix || !*prefix) return;
    auto dec3 = [](uint64_t k, int& a, int& b, int& c) {
      a = static_cast<int>(k >> 42);
      b = static_cast<int>((k >> 21) & 0x1fffff);
      c = static_cast<int>(k & 0x1fffff);
    };
    auto wmap = [&](const char* tag, const std::unordered_map<uint64_t, RN>& M, int side) {
      std::string path = std::string(prefix) + "." + tag + ".tsv";
      std::ofstream out(path);
      out << "key\tx\ty\tinner\ti\tg\n";
      for (const auto& [k, r] : M) {
        int x, y, inn;
        dec3(k, x, y, inn);
        if (g_ntr_K > 0) {
          for (const auto& e : r.T[side])
            out << k << '\t' << x << '\t' << y << '\t' << inn << '\t' << e.first << '\t'
                << std::setprecision(17) << e.second << '\n';
        } else {
          const RV& a = side ? r.R[0] : r.L[0];
          for (size_t i = 0; i < a.v.size(); ++i) {
            if (!(a.v[i] > 0.0)) continue;
            out << k << '\t' << x << '\t' << y << '\t' << inn << '\t' << (a.lo + (int)i) << '\t'
                << std::setprecision(17) << (a.mx + std::log(a.v[i])) << '\n';
          }
        }
      }
      std::fprintf(stderr, "[dump-p2] %s n=%zu\n", path.c_str(), M.size());
    };
    auto wkey = [&](const char* tag, const std::unordered_map<int, RN>& M, int side) {
      std::string path = std::string(prefix) + "." + tag + ".tsv";
      std::ofstream out(path);
      out << "key\ti\tg\n";
      for (const auto& [k, r] : M) {
        if (g_ntr_K > 0) {
          for (const auto& e : r.T[side])
            out << k << '\t' << e.first << '\t' << std::setprecision(17) << e.second << '\n';
        }
      }
      std::fprintf(stderr, "[dump-p2] %s n=%zu\n", path.c_str(), M.size());
    };
    wmap("Lambda_L_root", prRR, 0);
    wmap("Lambda_R_root", prRR, 1);
    wmap("Lambda_L_node", prRG, 0);
    wmap("Lambda_R_node", prRG, 1);
    wkey("K_L_ls", prRKls, 0);
    wkey("K_R_ls", prRKls, 1);
    wkey("K_L_rs", prRKrs, 0);
    wkey("K_R_rs", prRKrs, 1);
    for (int s = 0; s < 2; ++s) {
      std::string path = std::string(prefix) + ".V_s" + std::to_string(s) + ".tsv";
      std::ofstream out(path);
      out << "xs\ty\tw\n";
      for (const auto& [xs, vec] : prVkey[s])
        for (const auto& e : vec)
          out << xs << '\t' << e.first << '\t' << std::setprecision(17) << e.second << '\n';
      std::fprintf(stderr, "[dump-p2] %s keys=%zu\n", path.c_str(), prVkey[s].size());
    }
    if (stem_) {
      std::string path = std::string(prefix) + ".Qz_pre.tsv";
      std::ofstream out(path);
      out << "i\tj\tz\n";
      long long nq = 0;
      if (!stem_->qz_col.empty()) {
        for (int j = 0; j < n; ++j)
          for (const auto& e : stem_->qz_col[static_cast<size_t>(j)]) {
            out << e.k << '\t' << j << '\t' << std::setprecision(17) << e.q << '\n';
            ++nq;
          }
      } else if (!stem_->prQz.empty()) {
        for (int i = 0; i < n; ++i)
          for (int j = i; j < n; ++j) {
            const double z = stem_->prQz[static_cast<size_t>(i) * n + j];
            if (z < HALF_NEG) continue;
            out << i << '\t' << j << '\t' << std::setprecision(17) << z << '\n';
            ++nq;
          }
      }
      std::fprintf(stderr, "[dump-p2] %s n=%lld\n", path.c_str(), nq);
    }
  }

  // --grow-trunc: the kept (index, log g) of g_L (q') and g_R (p') of an evaluated root / key,
  // summed over the used layers; per step (u_L, u_R move with q).
  struct RT {
    std::vector<std::pair<int, double>> L, R;
  };
  std::unordered_map<const RN*, RT> prRT;
  std::vector<int> tr_before[2], tr_after[2], tr_kL[2], tr_kR[2];
  double tr_share[2] = {}, tr_zero[2] = {};
  double tr_uL(int qp) {
    if (qp + 1 > n - 1) return NEG_INF;
    const PrUL& u = pr_ul(qp);
    return u.lmax > HALF_NEG ? u.phi : NEG_INF;
  }
  double tr_uR(int pp) {
    if (pp < 1) return NEG_INF;
    const PrUR& r = pr_ur(pp);
    return r.lmax > HALF_NEG ? r.r1 : NEG_INF;
  }
  const RT& rec_trunc(const RN& r, int who) {
    auto it = prRT.find(&r);
    if (it != prRT.end()) return it->second;
    if (g_ntr_K > 0) throw std::runtime_error("--grow-node-trunc: no layer vectors for rec_trunc");
    const auto t0 = std::chrono::steady_clock::now();
    pr_ub_init();
    pr_ub_sync();
    // entries: side, index, log g, log score
    struct E4 {
      int side, i;
      double g, s;
    };
    std::vector<E4> e[2];
    for (int side = 0; side < 2; ++side) {
      const RV* a = side ? r.R : r.L;
      const int nl = std::min(side ? g_ml_R + 1 : g_ml_L + 1, side ? rnR : rnL);
      int lo = n, hi = -1;
      for (int m = 0; m < nl; ++m)
        if (!a[m].v.empty()) {
          lo = std::min(lo, a[m].lo);
          hi = std::max(hi, a[m].lo + static_cast<int>(a[m].v.size()) - 1);
        }
      for (int i = lo; i <= hi; ++i) {
        double g = NEG_INF;
        for (int m = 0; m < nl; ++m) {
          const int k = i - a[m].lo;
          if (k >= 0 && k < static_cast<int>(a[m].v.size()) && a[m].v[k] > 0.0)
            g = log_add(g, a[m].mx + std::log(a[m].v[k]));
        }
        if (g < HALF_NEG) continue;
        // --g-key suffix: rank g_L(q') by G_L * S̄_pre(q'+1). Beam-DP03 uses this
        // path (--grow-trunc); LR-DP03 ranks in racc_finish_trunc.
        double u;
        if (p2_g_suffix && side == 0 && stem_)
          u = stem_->prS[static_cast<size_t>(i + 1)];
        else
          u = side ? tr_uR(i) : tr_uL(i);
        e[side].push_back({side, i, g, u > HALF_NEG ? g + u : NEG_INF});
      }
    }
    RT t;
    const size_t before = e[0].size() + e[1].size();
    double share_num = NEG_INF, share_den = NEG_INF;
    auto pick = [&](std::vector<E4>& v) {
      double tot = NEG_INF;
      for (const E4& x : v) tot = log_add(tot, x.s);
      const size_t cap = std::min(v.size(), static_cast<size_t>(g_tr_K));
      std::vector<size_t> ord(v.size());
      for (size_t k = 0; k < v.size(); ++k) ord[k] = k;
      double kept = NEG_INF;
      if (tot < HALF_NEG) {
        tr_zero[who] += 1.0;
        std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return v[a].g > v[b].g; });
        for (size_t c = 0; c < cap; ++c) (v[ord[c]].side ? t.R : t.L).emplace_back(v[ord[c]].i, v[ord[c]].g);
      } else {
        std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return v[a].s > v[b].s; });
        const double target = g_tr_eps > 0.0 ? tot + std::log1p(-g_tr_eps) : tot;
        size_t c = 0;
        for (size_t k : ord) {
          if (c >= cap || v[k].s < HALF_NEG) break;
          (v[k].side ? t.R : t.L).emplace_back(v[k].i, v[k].g);
          ++c;
          kept = log_add(kept, v[k].s);
          if (g_tr_eps > 0.0 && kept >= target) break;
        }
        share_num = log_add(share_num, kept);
        share_den = log_add(share_den, tot);
      }
    };
    if (g_tr_pool) {
      std::vector<E4> all = e[0];
      all.insert(all.end(), e[1].begin(), e[1].end());
      pick(all);
    } else {
      pick(e[0]);
      pick(e[1]);
    }
    tr_before[who].push_back(static_cast<int>(before));
    tr_after[who].push_back(static_cast<int>(t.L.size() + t.R.size()));
    tr_kL[who].push_back(static_cast<int>(t.L.size()));
    tr_kR[who].push_back(static_cast<int>(t.R.size()));
    tr_share[who] += share_den > HALF_NEG ? std::exp(share_num - share_den) : 1.0;
    if (p2_g_suffix || p2_g_mode > 0) ntr_kept[0] += static_cast<long long>(t.L.size());
    gm_step += 16.0 * static_cast<double>(t.L.capacity() + t.R.capacity()) + sizeof(RT) + 32;
    gm_update();
    pr_sel_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return prRT.emplace(&r, std::move(t)).first->second;
  }

  // Score of a stem with inner pair u.v from the RN r of its outer pairs: per layer m,
  //   SL[m] = sum_q' L[m][q'] tau_q'(u, v)     (sep: sum_q' L[m][q'] / Qz[q+1][q'] times tau_q)
  //   SR[m] = sum_p' R[m][p'] kappa_p'(u, v)   (--pr-abl-pk K: the top K p' by the used g_R)
  // and sigma = sum_{m <= ml_L} SL[m] + sum_{m <= ml_R} SR[m] (pr_eval_sep with the layer sums
  // for g_L, g_R).  --grow-ml-diag: every layer evaluated, mass kept of each configuration.
  std::vector<double> rTau;
  std::vector<int> rTauStamp;
  int rTauSt = 0;
  static constexpr int RD_CFG = 5, RD_BINS = 2000;
  static constexpr int RD_L[RD_CFG] = {2, 2, 1, 1, 0}, RD_R[RD_CFG] = {2, 0, 0, 1, 0};
  long long rd_n = 0, rd_hist[RD_CFG][RD_BINS + 1] = {};
  double rd_sum[RD_CFG] = {}, rd_layer_sum[2][3] = {};
  double rd_min[RD_CFG] = {1, 1, 1, 1, 1};
  double rec_tau(int qp, int u, int v) {
    if (rTauStamp.size() < static_cast<size_t>(n)) {
      rTauStamp.assign(n, 0);
      rTau.assign(n, NEG_INF);
    }
    if (rTauStamp[qp] == rTauSt) return rTau[qp];
    rTauStamp[qp] = rTauSt;
    PrPartnerTimer pt(*this);
    return rTau[qp] = pr_T(qp, u, v, false).ext;
  }
  double pr_eval_rec(const RN& r, int u, int v, int who) {
    if (g_ntr_K > 0) throw std::runtime_error("--grow-node-trunc: pr_eval_rec needs the layer vectors");
    if (!E.end_ok(u, v)) return NEG_INF;
    ++rTauSt;
    ++pr_sep_evals;
    const int q = pr_abl_q;
    const bool vecL = !(pr_abl & PA_Q);
    if (g_tr_K > 0) {
      const RT& t = rec_trunc(r, who);
      double sl = NEG_INF, sr = NEG_INF;
      for (const auto& [qp, g] : t.L) {
        if (vecL) {
          const double ta = rec_tau(qp, u, v);
          if (ta > HALF_NEG) sl = log_add(sl, g + ta);
        } else {
          const double z = pr_abl_qz(q, qp);
          if (z > HALF_NEG) sl = log_add(sl, g - z);
        }
      }
      if (!vecL && sl > HALF_NEG) {
        const double ta = rec_tau(q, u, v);
        sl = ta > HALF_NEG ? sl + ta : NEG_INF;
      }
      for (const auto& [pp, g] : t.R) {
        PrPartnerTimer pt(*this);
        const double kx = pr_K(pp, u, v).ext;
        if (kx > HALF_NEG) sr = log_add(sr, g + kx);
      }
      const double s = log_add(sl, sr);
      if (s < HALF_NEG) return NEG_INF;
      return s + (E.w_stem + E.at(u, v));
    }
    const int nL = g_ml_diag ? rnL : g_ml_L + 1, nR = g_ml_diag ? rnR : g_ml_R + 1;
    double SL[3] = {NEG_INF, NEG_INF, NEG_INF}, SR[3] = {NEG_INF, NEG_INF, NEG_INF};
    for (int m = 0; m < nL; ++m) {
      const RV& a = r.L[m];
      double acc = NEG_INF;
      for (size_t i = 0; i < a.v.size(); ++i) {
        if (a.v[i] <= 0.0) continue;
        const int qp = a.lo + static_cast<int>(i);
        const double g = a.mx + std::log(a.v[i]);
        if (vecL) {
          const double t = rec_tau(qp, u, v);
          if (t > HALF_NEG) acc = log_add(acc, g + t);
        } else {
          const double z = pr_abl_qz(q, qp);
          if (z > HALF_NEG) acc = log_add(acc, g - z);
        }
      }
      if (!vecL && acc > HALF_NEG) {
        const double t = rec_tau(q, u, v);
        acc = t > HALF_NEG ? acc + t : NEG_INF;
      }
      SL[m] = acc;
    }
    if (pr_abl_pk > 0) {
      // top K p' by the configured g_R, then every layer over those p'
      std::vector<std::pair<double, int>> gr;
      std::vector<double> g(n, NEG_INF);
      for (int m = 0; m <= g_ml_R && m < rnR; ++m) {
        const RV& a = r.R[m];
        for (size_t i = 0; i < a.v.size(); ++i)
          if (a.v[i] > 0.0) g[a.lo + i] = log_add(g[a.lo + i], a.mx + std::log(a.v[i]));
      }
      for (int p = 0; p < n; ++p)
        if (g[p] > HALF_NEG) gr.emplace_back(g[p], p);
      std::stable_sort(gr.begin(), gr.end(), [](const auto& x, const auto& y) { return x.first > y.first; });
      if (gr.size() > static_cast<size_t>(pr_abl_pk)) gr.resize(static_cast<size_t>(pr_abl_pk));
      for (const auto& [gw, pp] : gr) {
        PrPartnerTimer pt(*this);
        const double kx = pr_K(pp, u, v).ext;
        if (kx < HALF_NEG) continue;
        for (int m = 0; m < nR; ++m) {
          const RV& a = r.R[m];
          const int i = pp - a.lo;
          if (i < 0 || i >= static_cast<int>(a.v.size()) || a.v[i] <= 0.0) continue;
          SR[m] = log_add(SR[m], a.mx + std::log(a.v[i]) + kx);
        }
      }
    } else {
      for (int m = 0; m < nR; ++m) {
        const RV& a = r.R[m];
        double acc = NEG_INF;
        for (size_t i = 0; i < a.v.size(); ++i) {
          if (a.v[i] <= 0.0) continue;
          PrPartnerTimer pt(*this);
          const double kx = pr_K(a.lo + static_cast<int>(i), u, v).ext;
          if (kx > HALF_NEG) acc = log_add(acc, a.mx + std::log(a.v[i]) + kx);
        }
        SR[m] = acc;
      }
    }
    double s = NEG_INF;
    for (int m = 0; m <= g_ml_L && m < nL; ++m) s = log_add(s, SL[m]);
    for (int m = 0; m <= g_ml_R && m < nR; ++m) s = log_add(s, SR[m]);
    if (g_ml_diag) {
      double all = NEG_INF;
      for (int m = 0; m < 3; ++m) all = log_add(all, log_add(SL[m], SR[m]));
      if (all > HALF_NEG) {
        ++rd_n;
        for (int m = 0; m < 3; ++m) {
          if (SL[m] > HALF_NEG) rd_layer_sum[0][m] += std::exp(SL[m] - all);
          if (SR[m] > HALF_NEG) rd_layer_sum[1][m] += std::exp(SR[m] - all);
        }
        for (int c = 0; c < RD_CFG; ++c) {
          double sc = NEG_INF;
          for (int m = 0; m <= RD_L[c]; ++m) sc = log_add(sc, SL[m]);
          for (int m = 0; m <= RD_R[c]; ++m) sc = log_add(sc, SR[m]);
          const double f = sc > HALF_NEG ? std::min(1.0, std::exp(sc - all)) : 0.0;
          rd_sum[c] += f;
          rd_min[c] = std::min(rd_min[c], f);
          // bins of -log10 f, width 0.01; the last bin holds f < 1e-20 and f = 0
          const double d = f > 0.0 ? -std::log10(f) : 1e9;
          const int b = std::min(RD_BINS, static_cast<int>(d * 100.0));
          ++rd_hist[c][b];
        }
      }
    }
    if (s < HALF_NEG) return NEG_INF;
    return s + (E.w_stem + E.at(u, v));
  }
  // --grow-rec-check: recursion score and summed vectors against the profile of the state.
  long long rc_n = 0, rc_vec_n = 0;
  double rc_max_score = 0.0, rc_max_vec = 0.0;
  void rec_check(const RN& r, const PProf& kp, double s_rec, double s_prof) {
    ++rc_n;
    if (s_rec > HALF_NEG || s_prof > HALF_NEG) {
      const double d = (s_rec > HALF_NEG && s_prof > HALF_NEG) ? std::fabs(s_rec - s_prof) : INFINITY;
      rc_max_score = std::max(rc_max_score, d);
    }
    const StemEst& st = *stem_;
    std::vector<double> L(n, NEG_INF), R(n, NEG_INF);
    for (size_t k = 0; k < kp.k.size(); ++k) {
      const int pp = k_hi(kp.k[k]), qp = k_lo(kp.k[k]);
      L[qp] = log_add(L[qp], kp.w[k] + prefix_log(pp));
      R[pp] = log_add(R[pp], kp.w[k] + E.w_beta1 + st.prS[qp + 1]);
    }
    auto cmp = [&](const RV* a, const std::vector<double>& ref) {
      std::vector<double> x(n, NEG_INF);
      double mx = NEG_INF;
      for (int m = 0; m < 3; ++m)
        for (size_t i = 0; i < a[m].v.size(); ++i)
          if (a[m].v[i] > 0.0) x[a[m].lo + i] = log_add(x[a[m].lo + i], a[m].mx + std::log(a[m].v[i]));
      for (int i = 0; i < n; ++i) mx = std::max(mx, std::max(x[i], ref[i]));
      if (mx < HALF_NEG) return;
      for (int i = 0; i < n; ++i) {
        if (std::max(x[i], ref[i]) < mx - 27.6) continue;  // entries below 1e-12 of the max
        ++rc_vec_n;
        const double d = (x[i] > HALF_NEG && ref[i] > HALF_NEG) ? std::fabs(x[i] - ref[i]) : INFINITY;
        rc_max_vec = std::max(rc_max_vec, d);
      }
    };
    cmp(r.L, L);
    cmp(r.R, R);
  }

  // log Out of a stem with inner pair u.v whose outer pairs p'.q' carry the profile kp:
  // [end_ok(u, v)] e^beta2 at(u, v) sum kp(p', q') (Stop_L + Stop_R)(p', q').
  // --pr-abl sep (with bml, q): the separable form of pr_eval,
  //   Gamma_L tau_q(u, v) + sum over the top --pr-abl-pk p' (by g_R; all if 0) g_R(p') kappa_p'(u, v),
  //   Gamma_L = sum kp(p', q') C[p'-1] / Qz[q+1][q'],  g_R(p') = sum_q' kp(p', q') e^beta1 Sbar[q'+1],
  //   tau_q = T(q; u, v).ext,  kappa_p' = K(p'; u, v).ext.
  // Without q (sepv: bml only, exact): Stop_L as sum_q' g_L(q') tau_q'(u, v) with
  // g_L(q') = sum_p' kp(p', q') C[p'-1].
  std::vector<std::pair<double, int>> prSepR, prSepL;
  double pr_eval_sep(const PProf& kp, int u, int v) {
    const StemEst& st = *stem_;
    const int q = pr_abl_q;
    const bool vecL = !(pr_abl & PA_Q);
    double gl = NEG_INF;
    prSepR.clear();
    prSepL.clear();
    for (size_t k = 0; k < kp.k.size(); ++k) {
      const int pp = k_hi(kp.k[k]), qp = k_lo(kp.k[k]);
      if (vecL) {
        prSepL.emplace_back(kp.w[k] + prefix_log(pp), qp);
      } else {
        const double z = pr_abl_qz(q, qp);
        if (z > HALF_NEG) gl = log_add(gl, kp.w[k] + prefix_log(pp) - z);
      }
      prSepR.emplace_back(kp.w[k] + (E.w_beta1 + st.prS[qp + 1]), pp);
    }
    double sv = NEG_INF;
    if (vecL) {
      std::sort(prSepL.begin(), prSepL.end(), [](const auto& a, const auto& b) { return a.second < b.second; });
      for (size_t k = 0; k < prSepL.size();) {
        double g = NEG_INF;
        const int qp = prSepL[k].second;
        for (; k < prSepL.size() && prSepL[k].second == qp; ++k) g = log_add(g, prSepL[k].first);
        if (g < HALF_NEG) continue;
        PrPartnerTimer pt(*this);
        const PrTV& t = pr_T(qp, u, v, false);
        if (t.ext > HALF_NEG) sv = log_add(sv, g + t.ext);
      }
    }
    std::sort(prSepR.begin(), prSepR.end(), [](const auto& a, const auto& b) { return a.second < b.second; });
    size_t m = 0;
    for (size_t k = 0; k < prSepR.size(); ++k) {
      if (m > 0 && prSepR[m - 1].second == prSepR[k].second)
        prSepR[m - 1].first = log_add(prSepR[m - 1].first, prSepR[k].first);
      else prSepR[m++] = prSepR[k];
    }
    prSepR.resize(m);
    std::sort(prSepR.begin(), prSepR.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
    const size_t K = pr_abl_pk > 0 ? std::min(m, static_cast<size_t>(pr_abl_pk)) : m;
    ++pr_sep_evals;
    pr_sep_pp += static_cast<long long>(m);
    pr_sep_ppk += static_cast<long long>(K);
    double s = sv;
    if (gl > HALF_NEG) {
      PrPartnerTimer pt(*this);
      const PrTV& t = pr_T(q, u, v, false);
      if (t.ext > HALF_NEG) s = gl + t.ext;
    }
    for (size_t k = 0; k < K; ++k) {
      if (prSepR[k].first < HALF_NEG) continue;
      PrPartnerTimer pt(*this);
      const PrKV& kv = pr_K(prSepR[k].second, u, v);
      if (kv.ext > HALF_NEG) s = log_add(s, prSepR[k].first + kv.ext);
    }
    if (s < HALF_NEG) return NEG_INF;
    return s + (E.w_stem + E.at(u, v));
  }

  double pr_eval(const PProf& kp, int u, int v) {
    if (!E.end_ok(u, v)) return NEG_INF;
    if (pr_abl & PA_SEP) return pr_eval_sep(kp, u, v);
    if (pr_sc_n >= PR_SC_MAX) {
      prSc.clear();
      pr_sc_n = 0;
    }
    std::unordered_map<uint64_t, double>& b = prSc[key2(u, v)];
    double acc = NEG_INF;
    for (size_t k = 0; k < kp.k.size(); ++k) {
      const double s = pr_SLR(b, k_hi(kp.k[k]), k_lo(kp.k[k]), u, v);
      if (s > HALF_NEG) acc = log_add(acc, kp.w[k] + s);
    }
    pr_evals += static_cast<long long>(kp.k.size());
    if (pr_gcheck)
      for (size_t k = 0; k < kp.k.size(); ++k) pr_gcheck_ep(k_hi(kp.k[k]), k_lo(kp.k[k]), u, v);
    if (acc < HALF_NEG) return NEG_INF;
    return acc + (E.w_stem + E.at(u, v));
  }

  // --pr-grow-check: Stop_L over every arm in the gap (exact placement) and as scored (the
  // selected arms), Stop_R as scored, against Qz[u+1][v-1] U_L and Qz[u+1][v-1] U_R.
  void pr_gcheck_ep(int pp, int qp, int u, int v) {
    const double g = pr_qz(u + 1, v - 1);
    if (g < HALF_NEG) return;
    const auto [ul, ur] = pr_U(pp, qp);
    double sl_all = NEG_INF, sl_sel = NEG_INF, sr = NEG_INF;
    if (qp + 1 <= n - 1 && u + 1 <= v - 1) {
      ArmQ& a = pr_armq(qp);
      acc_begin();
      for (int pos = prXpFirst[u + 1]; pos < prXpFirst[v]; ++pos) {
        const int r = prByXp[pos];
        if (prArmXs[r] < u + 1 || a.phi[r] < HALF_NEG) continue;
        const double gl = pr_qz(u + 1, prArmXs[r] - 1), gr = pr_qz(prArmXp[r] + 1, v - 1);
        if (gl < HALF_NEG || gr < HALF_NEG) continue;
        for (const auto& [ys, w] : pr_W(a, r, qp)) acc_add(ys, (gl + gr) + w);
      }
      for (int ys : prTouched) sl_all = log_add(sl_all, prAcc[ys] + pr_ctx(pp, ys));
      const PrTV& t = pr_T(qp, u, v, false);
      for (size_t k = 0; k < t.ys.size(); ++k) sl_sel = log_add(sl_sel, t.w[k] + pr_ctx(pp, t.ys[k]));
    }
    const PrKV& kv = pr_K(pp, u, v);
    for (size_t k = 0; k < kv.is.size(); ++k) sr = log_add(sr, kv.w[k] + pr_ctx(kv.is[k], qp));
    const double slack = std::log1p(1e-12);
    ++gc_n;
    if (sl_all > HALF_NEG) {
      const double d = ul > HALF_NEG ? sl_all - (g + ul) : INFINITY;
      gc_maxL = std::max(gc_maxL, d);
      if (d > slack) ++gc_vL;
    }
    if (sr > HALF_NEG) {
      const double d = ur > HALF_NEG ? sr - (g + ur) : INFINITY;
      gc_maxR = std::max(gc_maxR, d);
      if (d > slack) ++gc_vR;
    }
    if (sl_sel > HALF_NEG && !(sl_all > HALF_NEG && sl_sel <= sl_all + slack)) ++gc_sel_gt;
  }

  // log sum_i alpha(Qgl[i][d-1][e][f]) at(i, f) Ctx(i, q) over the surviving Qgl; per step.
  std::unordered_map<uint64_t, double> prAc;
  double pr_actx(int d, int e, int f, int q) {
    const uint64_t key = key3(d, e, f);
    auto it = prAc.find(key);
    if (it != prAc.end()) return it->second;
    const StemEst& st = *stem_;
    double s = NEG_INF;
    auto ia = stA.find(key2(d, e));
    if (ia != stA.end()) {
      const StemA& sa = ia->second;
      const auto kf = std::upper_bound(sa.fs.begin(), sa.fs.end(), f) - sa.fs.begin() - 1;
      if (kf >= 0 && sa.fs[kf] == f && sa.vals[kf] > HALF_NEG) {
        double bm = NEG_INF;
        if (!(pr_abl & PA_BML))
          for (const auto& [i, lw] : sa.il[kf]) bm = log_add(bm, lw + st.Bml(i, q));
        s = log_add((E.w_beta1 + st.prS[q + 1]) + sa.vals[kf], E.w_beta1m + bm);
      }
    }
    prAc.emplace(key, s);
    return s;
  }

