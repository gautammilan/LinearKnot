// Scan scores read from the growth vectors (lsc_L, lsc_R) and the l2r beam selection
// (l2r_prune4).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  double lsc_L(const RN& r, int arm, int who) {
    const auto t0 = std::chrono::steady_clock::now();
    auto done = [&](double s) {
      const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      sx_k_sec += dt;
      kl_in_read_sec += dt;
      return s;
    };
    if (arm < 0) return done(NEG_INF);
    double s = NEG_INF;
    auto note = [&](int m, int list) {
      ++sx_k_n;
      ++sx_k_src[sx_ksrc];
      if (m >= 0 && m < 8) ++sx_m_h[m];
      if (list >= 0) ++sx_list_h[std::min(list, 64)];
    };
    auto note_arm = [&]() {
      if (arm < 0) return;
      if (arm_seen.empty()) arm_seen.assign(static_cast<size_t>(std::max(prA, 1)), 0);
      if (arm >= prA || arm_seen[static_cast<size_t>(arm)]) return;
      arm_seen[static_cast<size_t>(arm)] = 1;
      ++arm_read_n;
      // |{f : V(c,f)[b] kept}| for arm (b,c)=(xs,xp)
      int nf = 0;
      if (stem_ && !stem_->prMu.empty()) {
        const int b = prArmXs[arm], c = prArmXp[arm];
        for (int f : pr_pair_list(true)[c]) {
          auto it = stem_->prMu.find(key2(c, f));
          if (it == stem_->prMu.end()) continue;
          if (it->second.find(key2(b, 0)) != it->second.end()) ++nf;
        }
      }
      arm_Vf_sum += nf;
      if (nf > arm_Vf_max) arm_Vf_max = nf;
      arm_Vf_samples.push_back(nf);
      // |{q' : K_L(arm,q') nonzero}|
      int n_kq = 0;
      for (int qp = 0; qp + 1 <= n - 1; ++qp)
        if (lK(arm, qp) > HALF_NEG) ++n_kq;
      arm_Kq_sum += n_kq;
      if (n_kq > arm_Kq_max) arm_Kq_max = n_kq;
    };
    if (g_ntr_K > 0) {
      const int len = static_cast<int>(r.T[0].size());
      for (const auto& [qp, g] : r.T[0]) {
        note(-1, len);
        const double k = lK(arm, qp);
        if (k > HALF_NEG) s = log_add(s, g + k);
      }
      note_arm();
      return done(s);
    }
    if (g_tr_K > 0) {
      const auto& L = rec_trunc(r, who).L;
      const int len = static_cast<int>(L.size());
      for (const auto& [qp, g] : L) {
        note(-1, len);
        const double k = lK(arm, qp);
        if (k > HALF_NEG) s = log_add(s, g + k);
      }
      return done(s);
    }
    for (int m = 0; m <= g_ml_L && m < rnL; ++m) {
      const RV& a = r.L[m];
      for (size_t i = 0; i < a.v.size(); ++i) {
        if (a.v[i] <= 0.0) continue;
        note(m, -1);
        sx_m_sum += m;
        const double k = lK(arm, a.lo + static_cast<int>(i));
        if (k > HALF_NEG) s = log_add(s, a.mx + std::log(a.v[i]) + k);
      }
    }
    return done(s);
  }
  double lsc_R(const RN& r, int xs) const {
    if (g_ntr_K > 0 && g_ntr_KR < 0) {
      const RV& a = r.R[0];
      const int k = xs - a.lo;
      return k >= 0 && k < static_cast<int>(a.v.size()) && a.v[k] > 0.0 ? a.mx + std::log(a.v[k]) : NEG_INF;
    }
    if (g_ntr_K > 0) {
      const auto& t = r.T[1];
      auto it = std::lower_bound(t.begin(), t.end(), xs,
                                 [](const std::pair<int, double>& e, int v) { return e.first < v; });
      if (it != t.end() && it->first == xs) return it->second;
      ++ntr_rmiss;
      return r.rest;
    }
    double s = NEG_INF;
    for (int m = 0; m <= g_ml_R && m < rnR; ++m) {
      const RV& a = r.R[m];
      const int k = xs - a.lo;
      if (k >= 0 && k < static_cast<int>(a.v.size()) && a.v[k] > 0.0)
        s = log_add(s, a.mx + std::log(a.v[k]));
    }
    return s;
  }

  // Top l2r_beam entries of t by score(group, entry).
  template <class Score>
