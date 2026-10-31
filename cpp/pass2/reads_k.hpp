// Principled-B reads (lazy path): W, T, K_L (pr_K), context, S_L / S_R (pr_SLR), U_L / U_R / U,
// growth selection (pr_gselect by iota = w (U_L + U_R)), growth build, roots and keys. The paper
// preset runs only pr_ub_init and pr_ls_closings here; its K_L read is lK in pass3/scan_scores.hpp.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  const std::vector<std::pair<int, double>>& pr_W(ArmQ& a, int r, int qp) {
    auto it = a.W.find(r);
    if (it != a.W.end()) return it->second;
    if (pr_vxs) throw std::runtime_error("--pr-v-key xs: pr_W needs per-ys arm entries");
    SecTimer tw(pr_sec_W);
    ++pr_W_built;
    const StemEst& st = *stem_;
    std::vector<std::pair<int, double>> w;
    const int hi = prArmE0[r + 1];
    for (int k = static_cast<int>(std::upper_bound(prEy.begin() + prArmE0[r], prEy.begin() + hi, qp) -
                                  prEy.begin());
         k < hi; ++k) {
      const double z = st.Qzbar(qp + 1, prEy[k] - 1);
      if (z > HALF_NEG) w.emplace_back(prEys[k], z + prEw[k]);
      ++pr_W_ent;
    }
    std::stable_sort(w.begin(), w.end(), [](const std::pair<int, double>& x,
                                            const std::pair<int, double>& y) {
      return x.first < y.first;
    });
    std::vector<std::pair<int, double>> m;
    for (const auto& e : w) {
      if (!m.empty() && m.back().first == e.first) m.back().second = log_add(m.back().second, e.second);
      else m.push_back(e);
    }
    return a.W.emplace(r, std::move(m)).first->second;
  }

  // Dense accumulator over [0, n) (stamped).
  std::vector<double> prAcc;
  std::vector<int> prStamp, prTouched;
  int pr_stamp = 0;
  void acc_begin() {
    if (prAcc.size() != static_cast<size_t>(n)) {
      prAcc.assign(n, NEG_INF);
      prStamp.assign(n, 0);
    }
    ++pr_stamp;
    prTouched.clear();
  }
  void acc_add(int i, double w) {
    if (prStamp[i] != pr_stamp) {
      prStamp[i] = pr_stamp;
      prAcc[i] = NEG_INF;
      prTouched.push_back(i);
    }
    prAcc[i] = log_add(prAcc[i], w);
  }

  // T(q'; u, v)(ys) = sum over the selected arms (xs, x') of Qz[u+1][xs-1] Qz[x'+1][v-1]
  // W(arm)(ys) (two: xs = u+1, no left filler); ext = log sum_ys T e^beta1 Sbar[ys+1].
  // Cached by q' until the scan passes q' (no term depends on the step).
  struct PrTV {
    std::vector<int> ys;
    std::vector<double> w;
    double ext = NEG_INF;
  };
  std::unordered_map<int, std::unordered_map<uint64_t, PrTV>> prTc;
  std::vector<int> prSel;

  const PrTV& pr_T(int qp, int u, int v, bool two) {
    std::unordered_map<uint64_t, PrTV>& bucket = prTc[qp];
    const uint64_t key = key3(u, v, two ? 1 : 0);
    auto it = bucket.find(key);
    if (it != bucket.end()) {
      ++pr_T_hits;
      return it->second;
    }
    ++pr_T_built;
    PrPartnerTimer pt(*this);
    SecTimer tt(pr_sec_T);
    const StemEst& st = *stem_;
    PrTV t;
    if (qp + 1 <= n - 1) {
      ArmQ& a = pr_armq(qp);
      const bool xs = pr_abl & (PA_XS | PA_XB);
      if (pr_abl & PA_XS) {
        const int k0 = pr_arms_k;
        const bool r0 = pr_ratio;
        pr_arms_k = 0;
        pr_ratio = false;
        pr_select(a, u, v, two, prSel);
        pr_arms_k = k0;
        pr_ratio = r0;
      } else {
        pr_select(a, u, v, two, prSel);
      }
      pr_T_arms += static_cast<long long>(prSel.size());
      const double gz = xs ? pr_qz(u + 1, v - 1) : 0.0;
      acc_begin();
      for (int r : prSel) {
        double g;
        if (xs) {
          const double za = pr_qz(prArmXs[r], prArmXp[r]);
          if (gz < HALF_NEG || za < HALF_NEG) continue;
          g = gz - za;
        } else {
          const double gl = two ? 0.0 : pr_qz(u + 1, prArmXs[r] - 1);
          const double gr = pr_qz(prArmXp[r] + 1, v - 1);
          if (gl < HALF_NEG || gr < HALF_NEG) continue;
          g = gl + gr;
        }
        for (const auto& [ys, w] : pr_W(a, r, qp)) acc_add(ys, g + w);
      }
      std::sort(prTouched.begin(), prTouched.end());
      for (int ys : prTouched) {
        t.ys.push_back(ys);
        t.w.push_back(prAcc[ys]);
        t.ext = log_add(t.ext, prAcc[ys] + (E.w_beta1 + st.prS[ys + 1]));
      }
      pr_T_ys += static_cast<long long>(t.ys.size());
    }
    return bucket.emplace(key, std::move(t)).first->second;
  }

  // K(p'; u, v)(i) = sum over the surviving L-stems Qgl[i][p'-1][e][f], u+1 <= e <= f <= v-1,
  // of Qz[u+1][e-1] alpha at(i, f) Qz[f+1][v-1]; ext = log sum_i K C[i-1].  Final once the
  // scan passes v - 1, so kept across steps (cleared when it outgrows PR_KC_MAX).
  struct PrKV {
    std::vector<int> is;
    std::vector<double> w;
    double ext = NEG_INF;
  };
  static constexpr size_t PR_KC_MAX = size_t(1) << 21;
  std::unordered_map<uint64_t, PrKV> prKc;

  const PrKV& pr_K(int pp, int u, int v) {
    const uint64_t key = key3(pp, u, v);
    auto it = prKc.find(key);
    if (it != prKc.end()) {
      ++pr_K_hits;
      return it->second;
    }
    if (prKc.size() >= PR_KC_MAX) prKc.clear();
    ++pr_K_built;
    PrPartnerTimer pt(*this);
    SecTimer tk(pr_sec_K);
    PrKV kv;
    auto ie = stE.find(pp);
    if (ie != stE.end()) {
      acc_begin();
      const std::vector<int>& es = ie->second;
      for (auto ke = std::lower_bound(es.begin(), es.end(), u + 1); ke != es.end() && *ke <= v - 1;
           ++ke) {
        const int e = *ke;
        const double gl = pr_qz(u + 1, e - 1);
        if (gl < HALF_NEG) continue;
        auto ia = stA.find(key2(pp, e));
        if (ia == stA.end()) continue;
        const StemA& sa = ia->second;
        for (size_t kf = 0; kf < sa.fs.size(); ++kf) {
          const int f = sa.fs[kf];
          if (f > v - 1) break;
          const double gr = pr_qz(f + 1, v - 1);
          if (gr < HALF_NEG) continue;
          const double g = (pr_abl & PA_KS) ? pr_qz(u + 1, v - 1) - pr_qz(e, f) : gl + gr;
        if (!(g > HALF_NEG)) continue;
        for (const auto& [i, lw] : sa.il[kf]) acc_add(i, g + lw);
        }
      }
      std::sort(prTouched.begin(), prTouched.end());
      for (int i : prTouched) {
        kv.is.push_back(i);
        kv.w.push_back(prAcc[i]);
        kv.ext = log_add(kv.ext, prAcc[i] + prefix_log(i));
      }
      pr_K_is += static_cast<long long>(kv.is.size());
    }
    return prKc.emplace(key, std::move(kv)).first->second;
  }

  double pr_ctx(int i, int j) const {
    const StemEst& st = *stem_;
    return log_add(E.w_beta1 + prefix_log(i) + st.prS[j + 1], E.w_beta1m + st.Bml(i, j));
  }

  // Stop_L + Stop_R of an outer pair p'.q' for a stem with inner pair u.v, without that
  // stem's own beta2 + at(u, v) and the growth factor:
  //   Stop_L = sum_ys T(q'; u, v)(ys) Ctx(p', ys),  Stop_R = sum_i K(p'; u, v)(i) Ctx(i, q').
  // Final once built: cached by (u, v), then (p', q'); cleared when it outgrows PR_SC_MAX.
  static constexpr size_t PR_SC_MAX = size_t(1) << 22;
  std::unordered_map<uint64_t, std::unordered_map<uint64_t, double>> prSc;
  size_t pr_sc_n = 0;

  // --pr-abl q: T(q; u, v) / Qz[q+1][q'] for T(q'; u, v); p: K(p; u, v) for K(p'; u, v)
  // (p, q: the scored state's); bml: Ctx without its multiloop term.
  double pr_SLR_abl(int pp, int qp, int u, int v) {
    const StemEst& st = *stem_;
    const bool bml = !(pr_abl & PA_BML);
    const int tq = (pr_abl & PA_Q) && pr_abl_q >= 0 ? pr_abl_q : qp;
    const double tz = -pr_abl_qz(tq, qp);
    double sl = NEG_INF;
    if (tz > HALF_NEG && -tz > HALF_NEG) {
      const PrTV& t = pr_T(tq, u, v, false);
      if (!t.ys.empty()) {
        double bm = NEG_INF;
        if (bml)
          for (size_t k = 0; k < t.ys.size(); ++k) bm = log_add(bm, t.w[k] + st.Bml(pp, t.ys[k]));
        sl = tz + log_add(prefix_log(pp) + t.ext, E.w_beta1m + bm);
      }
    }
    const int kp = (pr_abl & PA_P) && pr_abl_p >= 0 ? pr_abl_p : pp;
    double sr = NEG_INF;
    const PrKV& kv = pr_K(kp, u, v);
    if (!kv.is.empty()) {
      double bm = NEG_INF;
      if (bml)
        for (size_t k = 0; k < kv.is.size(); ++k) bm = log_add(bm, kv.w[k] + st.Bml(kv.is[k], qp));
      sr = log_add((E.w_beta1 + st.prS[qp + 1]) + kv.ext, E.w_beta1m + bm);
    }
    return log_add(sl, sr);
  }

  double pr_SLR(std::unordered_map<uint64_t, double>& b, int pp, int qp, int u, int v) {
    if (pr_abl & (PA_Q | PA_P)) {
      ++pr_slr_miss;
      PrPartnerTimer pt(*this);
      SecTimer ts(pr_sec_SLR);
      return pr_SLR_abl(pp, qp, u, v);
    }
    const uint64_t key = key2(pp, qp);
    auto it = b.find(key);
    if (it != b.end()) {
      ++pr_slr_hits;
      return it->second;
    }
    ++pr_slr_miss;
    PrPartnerTimer pt(*this);
    SecTimer ts(pr_sec_SLR);
    if (pr_abl & PA_BML) {
      const double s = pr_SLR_abl(pp, qp, u, v);
      b.emplace(key, s);
      ++pr_sc_n;
      return s;
    }
    const StemEst& st = *stem_;
    double sl = NEG_INF;
    {
      const PrTV& t = pr_T(qp, u, v, false);
      if (!t.ys.empty()) {
        double bm = NEG_INF;
        for (size_t k = 0; k < t.ys.size(); ++k) bm = log_add(bm, t.w[k] + st.Bml(pp, t.ys[k]));
        sl = log_add(prefix_log(pp) + t.ext, E.w_beta1m + bm);
      }
    }
    double sr = NEG_INF;
    const PrKV& kv = pr_K(pp, u, v);
    if (!kv.is.empty()) {
      double bm = NEG_INF;
      for (size_t k = 0; k < kv.is.size(); ++k) bm = log_add(bm, kv.w[k] + st.Bml(kv.is[k], qp));
      sr = log_add((E.w_beta1 + st.prS[qp + 1]) + kv.ext, E.w_beta1m + bm);
    }
    const double s = log_add(sl, sr);
    b.emplace(key, s);
    ++pr_sc_n;
    return s;
  }

  // --pr-grow-select bound.  At step q every gap lies in [0, q-1], every arm (xs, x') with
  // x' <= q-1 has its final Qz[xs][x'], and every partner entry with f <= q-1 is in stA.  An
  // arm inside the gap has Qz[u+1][xs-1] Qz[x'+1][v-1] <= Qz[u+1][v-1] / Qz[xs][x'] (a filler
  // of [u+1, v-1] contains every concatenation), so for every gap of the step
  //   Stop_L(p', q'; u, v) <= Qz[u+1][v-1] U_L(p', q'),
  //     U_L = sum over those arms and their ys of W(arm)(ys) Ctx(p', ys) / Qz[xs][x'],
  //   Stop_R(p', q'; u, v) <= Qz[u+1][v-1] U_R(p', q'),
  //     U_R = sum over those entries (p', e, f; i) of alpha at(i, f) Ctx(i, q') / Qz[e][f]
  // (the selected arms are a subset of all arms).  With dangles a cut b | b+1 can split an
  // unpaired run, adding d3(b) + d5(b+1), so each cut of the arm (xs - 1 | xs, x' | x' + 1)
  // or partner entry multiplies its 1 / Qz by e^kappa(b), kappa(b) = max(0, d3w[b] + d5w[b+1]).
  // Per q': the ys-vector of sum W / Qz, columns added as the scan passes them; per p': the
  // i-vector of sum alpha at / Qz.  The Bml parts of Ctx are dot products with exp(Bml) scaled
  // per row (U_L) or per column (U_R), in long double (the product of the two scalings stays
  // in range).
  struct PrPend {
    int d, e, f, i;
    double w;
  };
  std::vector<PrPend> prPend;
  size_t prPendDone = 0;
  int pr_q = -1, pr_sync_q = -1;
  struct PrUL {
    int known = -1;
    bool dirty = false;
    std::vector<double> w;
    std::vector<long double> lin;
    double phi = NEG_INF, lmax = NEG_INF;
    double acc = NEG_INF;  // --pr-v-key xs: phi accumulated directly (w / lin unused)
  };
  std::unordered_map<int, PrUL> prUL;
  struct PrUR {
    bool dirty = false;
    std::vector<double> w;
    std::vector<long double> lin;
    double r1 = NEG_INF, lmax = NEG_INF;
  };
  std::vector<PrUR> prUR;
  std::vector<long double> prEBr, prEBc;
  std::vector<double> prMBr, prMBc, prKap;
  std::unordered_map<uint64_t, std::pair<double, double>> prUc;
  // log 1 / Qz[a][b] with the cut factors of both ends (a >= 1)
  double pr_iqz(int a, int b) const {
    const double z = pr_qz(a, b);
    if (z < HALF_NEG) return NEG_INF;
    return (prKap[a - 1] + prKap[b]) - z;
  }

  void pr_ub_init() {
    if (!prMBr.empty()) return;
    prKap.assign(n + 1, 0.0);
    if (dz)
      for (int b = 0; b + 1 < n; ++b) prKap[b] = std::max(0.0, E.d3w[b] + E.d5w[b + 1]);
    prMBr.assign(n, NEG_INF);
    prMBc.assign(n, NEG_INF);
    prUR.assign(n, PrUR());
  }

  // Partner entries with f <= q-1 into their p' vectors (entries arrive in f order).
  void pr_ub_sync() {
    if (pr_sync_q == pr_q) return;
    pr_sync_q = pr_q;
    while (prPendDone < prPend.size() && prPend[prPendDone].f <= pr_q - 1) {
      const PrPend& e = prPend[prPendDone++];
      const double iz = e.e >= 1 ? pr_iqz(e.e, e.f) : NEG_INF;
      if (iz < HALF_NEG) continue;
      PrUR& r = prUR[e.d];
      if (r.w.empty()) r.w.assign(n, NEG_INF);
      r.w[e.i] = log_add(r.w[e.i], e.w + iz);
      r.dirty = true;
      ++gs_rent;
    }
  }

  PrUL& pr_ul(int qp) {
    const StemEst& st = *stem_;
    PrUL& u = prUL[qp];
    if (u.w.empty() && !pr_vxs) u.w.assign(n, NEG_INF);
    for (int c = u.known + 1; c <= pr_q - 1; ++c) {
      for (int pos = prXpFirst[c]; pos < prXpFirst[c + 1]; ++pos) {
        const int r = prByXp[pos];
        const double iz = prArmXs[r] >= 1 ? pr_iqz(prArmXs[r], c) : NEG_INF;
        if (iz < HALF_NEG) continue;
        if (pr_vxs) {
          const int hy = prArmY0[r + 1];
          for (int k = static_cast<int>(std::upper_bound(prYy.begin() + prArmY0[r], prYy.begin() + hy, qp) -
                                        prYy.begin());
               k < hy; ++k) {
            const double zb = st.Qzbar(qp + 1, prYy[k] - 1);
            if (zb < HALF_NEG) continue;
            u.acc = log_add(u.acc, (zb + prYw[k]) + iz);
            u.dirty = true;
          }
          continue;
        }
        const int hi = prArmE0[r + 1];
        for (int k = static_cast<int>(std::upper_bound(prEy.begin() + prArmE0[r], prEy.begin() + hi, qp) -
                                      prEy.begin());
             k < hi; ++k) {
          const double zb = st.Qzbar(qp + 1, prEy[k] - 1);
          if (zb < HALF_NEG) continue;
          double& t = u.w[prEys[k]];
          t = log_add(t, (zb + prEw[k]) + iz);
          u.dirty = true;
        }
      }
      ++gs_cols;
    }
    u.known = std::max(u.known, pr_q - 1);
    if (u.dirty && pr_vxs) {
      u.dirty = false;
      u.phi = u.lmax = u.acc;
    } else if (u.dirty) {
      u.dirty = false;
      u.phi = NEG_INF;
      u.lmax = NEG_INF;
      for (int y = qp + 1; y < n; ++y)
        if (u.w[y] > HALF_NEG) {
          u.phi = log_add(u.phi, u.w[y] + (E.w_beta1 + st.prS[y + 1]));
          u.lmax = std::max(u.lmax, u.w[y]);
        }
      u.lin.assign(n, 0.0L);
      for (int y = qp + 1; y < n; ++y)
        if (u.w[y] > HALF_NEG) u.lin[y] = std::exp(static_cast<long double>(u.w[y]) - u.lmax);
    }
    return u;
  }

  PrUR& pr_ur(int d) {
    PrUR& r = prUR[d];
    if (r.dirty) {
      r.dirty = false;
      r.r1 = NEG_INF;
      r.lmax = NEG_INF;
      for (int i = 0; i < d; ++i)
        if (r.w[i] > HALF_NEG) {
          r.r1 = log_add(r.r1, r.w[i] + prefix_log(i));
          r.lmax = std::max(r.lmax, r.w[i]);
        }
      r.lin.assign(n, 0.0L);
      for (int i = 0; i < d; ++i)
        if (r.w[i] > HALF_NEG) r.lin[i] = std::exp(static_cast<long double>(r.w[i]) - r.lmax);
    }
    return r;
  }

  // (log U_L, log U_R) of the endpoint p'.q' at the current step (cached per step).
  std::pair<double, double> pr_U(int pp, int qp) {
    const uint64_t key = key2(pp, qp);
    auto it = prUc.find(key);
    if (it != prUc.end()) return it->second;
    if (pr_vxs) throw std::runtime_error("--pr-v-key xs: pr_U needs per-ys arm entries");
    pr_ub_init();
    pr_ub_sync();
    ++gs_U;
    const StemEst& st = *stem_;
    double ul = NEG_INF, ur = NEG_INF;
    if (qp + 1 <= n - 1) {
      const PrUL& u = pr_ul(qp);
      if (u.lmax > HALF_NEG) {
        ul = prefix_log(pp) + u.phi;
        {
          long double dot = 0.0L;
          const long double* eb = &prEBr[ix(pp, 0)];
          for (int y = qp + 1; y < n; ++y) dot += u.lin[y] * eb[y];
          if (dot > 0.0L)
            ul = log_add(ul, (E.w_beta1m + u.lmax) + (prMBr[pp] + static_cast<double>(std::log(dot))));
        }
      }
    }
    if (pp >= 1) {
      const PrUR& r = pr_ur(pp);
      if (r.lmax > HALF_NEG) {
        ur = (E.w_beta1 + st.prS[qp + 1]) + r.r1;
        {
          long double dot = 0.0L;
          const long double* eb = &prEBc[ix(qp, 0)];
          for (int i = 0; i < pp; ++i) dot += r.lin[i] * eb[i];
          if (dot > 0.0L)
            ur = log_add(ur, (E.w_beta1m + r.lmax) + (prMBc[qp] + static_cast<double>(std::log(dot))));
        }
      }
    }
    return prUc.emplace(key, std::make_pair(ul, ur)).first->second;
  }

  // Keeps the entries of `all` (key order) with the largest iota = w (U_L + U_R) until their
  // sum reaches (1 - eps) of the total, at most B (0: no cap); a zero total keeps the B
  // largest by w.  kind 0: growth node, 1: evaluated root / key.
  void pr_gselect(PProf& all, PProf& out, int kind) {
    const auto t0 = std::chrono::steady_clock::now();
    const size_t m = all.k.size();
    ++gs_prof[kind];
    gs_before[kind] += static_cast<long long>(m);
    out.k.clear();
    out.w.clear();
    const int who = gs_who;
    {
      const auto [dp, dq] = gz_distinct(all.k);
      gz_before[who].push_back(static_cast<int>(m));
      gz_dpb[who].push_back(dp);
      gz_dqb[who].push_back(dq);
    }
    auto rec_after = [&]() {
      const auto [dp, dq] = gz_distinct(out.k);
      gz_after[who].push_back(static_cast<int>(out.k.size()));
      gz_dp[who].push_back(dp);
      gz_dq[who].push_back(dq);
    };
    std::vector<double> io(m);
    double tot = NEG_INF;
    for (size_t k = 0; k < m; ++k) {
      const auto [ul, ur] = pr_U(k_hi(all.k[k]), k_lo(all.k[k]));
      const double u = log_add(ul, ur);
      io[k] = u > HALF_NEG ? all.w[k] + u : NEG_INF;
      tot = log_add(tot, io[k]);
    }
    std::vector<size_t> ord(m);
    for (size_t k = 0; k < m; ++k) ord[k] = k;
    std::vector<char> keep(m, 0);
    const size_t cap = pr_gB > 0 ? std::min(m, static_cast<size_t>(pr_gB)) : m;
    if (tot < HALF_NEG) {
      ++gs_zero[kind];
      std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return all.w[a] > all.w[b]; });
      for (size_t c = 0; c < cap; ++c) keep[ord[c]] = 1;
      gs_share[kind] += 1.0;
    } else {
      std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return io[a] > io[b]; });
      const double target = pr_geps > 0.0 ? tot + std::log1p(-pr_geps) : tot;
      double kept = NEG_INF;
      size_t c = 0;
      for (size_t k : ord) {
        if (c >= cap || io[k] < HALF_NEG) break;
        keep[k] = 1;
        ++c;
        kept = log_add(kept, io[k]);
        if (pr_geps > 0.0 && kept >= target) break;
      }
      gs_share[kind] += std::exp(kept - tot);
    }
    for (size_t k = 0; k < m; ++k)
      if (keep[k]) {
        out.k.push_back(all.k[k]);
        out.w.push_back(all.w[k]);
      }
    gs_after[kind] += static_cast<long long>(out.k.size());
    rec_after();
    pr_sel_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }

  // Wall time of the outermost growth / root / key build, less the selection inside it.
  int pr_gdepth = 0;
  struct PrGrowTimer {
    DP03& d;
    std::chrono::steady_clock::time_point t0;
    double s0;
    explicit PrGrowTimer(DP03& x) : d(x), t0(std::chrono::steady_clock::now()), s0(x.pr_sel_sec) {
      ++d.pr_gdepth;
    }
    ~PrGrowTimer() {
      if (--d.pr_gdepth == 0) {
        const double dt =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() -
            (d.pr_sel_sec - s0);
        d.pr_grow_sec += dt;
      }
    }
  };

  // Growth profiles (scan left, final once built): G(x, y) = [end_ok] at(x, y) at key (x, y)
  // plus the steps' w G(x', y'), cut (bound: selected); roots (x, y, inner) the same,
  // untruncated (bound: selected per step when a Qg state evaluates them, pr_root_eval).
  std::unordered_map<uint64_t, PProf> prG, prGR;
  // --pr-fast 2: growth helpers HF / HM over the G profiles (final, so the helpers are too);
  // dropped whenever they hold more than PR_HG_MAX entries (only between node builds).
  std::unordered_map<uint64_t, PProf> prHg;
  static constexpr long long PR_HG_MAX = 8000000;
  void pr_hg_trim() {
    if (pr_hentries[0] - pr_hg_ent <= PR_HG_MAX) return;
    prHg.clear();
    pr_hg_ent = pr_hentries[0];
    ++pr_hg_clears;
    gm_help = 0;
  }
  auto pr_gchild() {
    return [this](int a, int b) -> const PProf* { return &prG.find(key2(a, b))->second; };
  }

  void pr_gbuild(int x, int y) {
    if (prG.count(key2(x, y))) return;
    std::vector<std::pair<int, int>> stk{{x, y}};
    std::vector<ChainNbr> nb;
    while (!stk.empty()) {
      const auto [a, b] = stk.back();
      if (prG.count(key2(a, b))) {
        stk.pop_back();
        continue;
      }
      pr_steps(a, b, false, true, nb);
      bool ready = true;
      for (const ChainNbr& t : nb)
        if (!prG.count(key2(t.x, t.y))) {
          stk.emplace_back(t.x, t.y);
          ready = false;
        }
      if (!ready) continue;
      count_steps(0, nb);
      PProf p;
      if (pr_fast > 0) {
        pr_hg_trim();
        pr_node_acc(a, b, false, true, nb, E.end_ok(a, b) ? E.at(a, b) : NEG_INF, pr_gchild(), prHg);
        if (pr_gsel) {
          PProf all;
          dacc_finish(all, 0.0);
          gs_who = 0;
          pr_gselect(all, p, 0);
        } else {
          pr_cut_drop += pr_finish_out_d(p, pr_cut);
        }
        pr_linify(p);
        ++pr_G_nodes;
        pr_G_entries += static_cast<long long>(p.k.size());
        gm_nodes += pbytes(p);
        gm_update();
        prG.emplace(key2(a, b), std::move(p));
        stk.pop_back();
        continue;
      }
      PAcc acc;
      if (E.end_ok(a, b)) acc.add(key2(a, b), E.at(a, b));
      for (const ChainNbr& t : nb) {
        const PProf& q = prG.find(key2(t.x, t.y))->second;
        for (size_t k = 0; k < q.k.size(); ++k) acc.add(q.k[k], t.w + q.w[k]);
      }
      if (pr_gsel) {
        PProf all;
        acc.finish(all, 0.0);
        gs_who = 0;
        pr_gselect(all, p, 0);
      } else {
        pr_cut_drop += pr_finish_out(acc, p, pr_cut);
      }
      ++pr_G_nodes;
      pr_G_entries += static_cast<long long>(p.k.size());
      gm_nodes += pbytes(p);
      gm_update();
      prG.emplace(key2(a, b), std::move(p));
      stk.pop_back();
    }
  }

  const PProf& pr_root(int x, int y, bool inner) {
    const uint64_t key = key3(x, y, inner ? 1 : 0);
    auto it = prGR.find(key);
    if (it != prGR.end()) return it->second;
    PrGrowTimer gt(*this);
    std::vector<ChainNbr> nb;
    pr_steps(x, y, inner, true, nb);
    count_steps(0, nb);
    for (const ChainNbr& t : nb) pr_gbuild(t.x, t.y);
    PProf p;
    if (pr_fast > 0) {
      pr_hg_trim();
      pr_node_acc(x, y, inner, true, nb, E.end_ok(x, y) ? E.at(x, y) : NEG_INF, pr_gchild(), prHg);
      dacc_finish(p, 0.0);
      pr_linify(p);
      ++pr_roots;
      gm_roots += pbytes(p);
      gm_update();
      return prGR.emplace(key, std::move(p)).first->second;
    }
    PAcc acc;
    if (E.end_ok(x, y)) acc.add(key2(x, y), E.at(x, y));
    for (const ChainNbr& t : nb) {
      const PProf& q = prG.find(key2(t.x, t.y))->second;
      for (size_t k = 0; k < q.k.size(); ++k) acc.add(q.k[k], t.w + q.w[k]);
    }
    acc.finish(p, 0.0);
    ++pr_roots;
    gm_roots += pbytes(p);
    gm_update();
    return prGR.emplace(key, std::move(p)).first->second;
  }

  // The root profile a Qg state at this step evaluates (bound: selected, cached per step).
  std::unordered_map<uint64_t, PProf> prGRs;
  const PProf& pr_root_eval(int x, int y, bool inner) {
    if (!pr_gsel) return pr_root(x, y, inner);
    const uint64_t key = key3(x, y, inner ? 1 : 0);
    auto it = prGRs.find(key);
    if (it != prGRs.end()) return it->second;
    PrGrowTimer gt(*this);
    PProf all = pr_root(x, y, inner), out;
    gs_who = 1;
    pr_gselect(all, out, 1);
    gm_step += pbytes(out);
    gm_update();
    return prGRs.emplace(key, std::move(out)).first->second;
  }

  // Key profiles of the helper states at step q: Qgls (i, seed) closes by (f) into
  // (i-1, j), j - q - 1 <= cap (seed) or pr_ml_close; Qgrs c by closings_rs.  Per step.
  std::unordered_map<int, PProf> prKls, prKrs;

  void pr_ls_closings(int i, bool seed, int q, std::vector<Closing>& out) const {
    out.clear();
    const int x = i - 1;
    if (x < 0) return;
    const int j_hi = std::min(n - 1, q + 1 + pr_ml_close);
    for (int j = q + 1; j <= j_hi; ++j) {
      const int unpaired = j - q - 1;
      if (seed && unpaired > cap) continue;
      if (!E.pairs(x, j) || !E.end_ok(x, j)) continue;
      double w_ext = E.w_ap + E.w_bp + E.w_cp * unpaired + pair_bonus;
      if (dz) w_ext += E.run(q + 1, j - 1);
      out.push_back({x, j, w_ext + E.at(x, j)});
    }
  }

  const PProf& pr_key(bool ls, int p, bool seed, int q) {
    std::unordered_map<int, PProf>& m = ls ? prKls : prKrs;
    const int key = ls ? 2 * p + (seed ? 1 : 0) : p;
    auto it = m.find(key);
    if (it != m.end()) return it->second;
    PrGrowTimer gt(*this);
    std::vector<Closing> cl;
    if (ls) pr_ls_closings(p, seed, q, cl);
    else closings_rs(p, q, cl);
    PProf out;
    if (pr_fast > 0) {
      std::vector<std::pair<const PProf*, double>> rs;
      rs.reserve(cl.size());
      for (const Closing& c : cl) rs.emplace_back(&pr_root(c.x, c.y, false), c.w);
      double ref = NEG_INF;
      for (const auto& [r, w] : rs)
        if (r->mx > HALF_NEG) ref = std::max(ref, w + r->mx);
      dacc_begin(ref > HALF_NEG ? ref : 0.0);
      for (const auto& [r, w] : rs) dacc_add(*r, w);
      if (pr_gsel) {
        PProf all;
        dacc_finish(all, 0.0);
        gs_who = 2;
        pr_gselect(all, out, 1);
      } else {
        dacc_finish(out, 0.0);
      }
      ++pr_keys;
      pr_key_entries += static_cast<long long>(out.k.size());
      gm_step += pbytes(out);
      gm_update();
      return m.emplace(key, std::move(out)).first->second;
    }
    PAcc acc;
    for (const Closing& c : cl) {
      const PProf& r = pr_root(c.x, c.y, false);
      for (size_t k = 0; k < r.k.size(); ++k) acc.add(r.k[k], c.w + r.w[k]);
    }
    if (pr_gsel) {
      PProf all;
      acc.finish(all, 0.0);
      gs_who = 2;
      pr_gselect(all, out, 1);
    } else {
      acc.finish(out, 0.0);
    }
    ++pr_keys;
    pr_key_entries += static_cast<long long>(out.k.size());
    gm_step += pbytes(out);
    gm_update();
    return m.emplace(key, std::move(out)).first->second;
  }

 public:
  // ------------------------------------------------ layered marginal vectors
  // RN of a growth node / root / key / helper: L[m][q'] = sum_p' G_m[p', q'] C[p'-1] and
  // R[m][p'] = sum_q' G_m[p', q'] e^beta1 Sbar[q'+1], G_m the growth weight with m multiloop
  // steps (m = 2: two or more).  Same steps, tables and caps as the profiles (pr_steps,
  // pr_split, pr_helper's ranges); interior steps keep m, multiloop steps add one.  Layers
  // m >= rnL (L) or >= rnR (R) are not built.  Values: e^mx v[i - lo].
  struct RV {
    int lo = 0;
    double mx = NEG_INF;
    std::vector<double> v;
  };
  struct RN {
    RV L[3], R[3];
    double tot = NEG_INF;  // log of the summed built layers (L and R), the --grow-ml-c rank
    // --bG: kept (index, log value) of L (0) and R (1), layers summed, index order;
    // rest: log of the R mass dropped in this node's subtree (bounds any missing R entry)
    std::vector<std::pair<int, double>> T[2];
    double rest = NEG_INF;
    double rsum = NEG_INF;  // --grow-r-from-v: log of the R side's sum (R itself not built)
  };
