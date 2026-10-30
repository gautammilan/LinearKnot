// Scan diagnostics: oracle dump / load (exact run rho), lsc_L / lsc_R, ldiag, ltr_* traces,
// prune statistics.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  // Exact run (beam inf, after the outside pass): rho(a') = log sum_ys Qp Qp_out -
  // log sum_ys Qp e^(beta1) C(a'-1) Sbar(ys+1); the outside of every U and Us state.
  void l2r_oracle_dump(const std::string& path) {
    const StemEst& st = *stem_;
    std::vector<double> num(n, NEG_INF), den(n, NEG_INF), dcx(n, NEG_INF);
    for (int ys = 0; ys < n; ++ys) {
      const OMap& at = Qp_at[ys];
      for (size_t r = 0; r < at.slots(); ++r) {
        const int a = at.keys[r];
        if (a < 0) continue;
        num[a] = log_add(num[a], at.vals[r] + Qp_b[ix(a, ys)]);
        den[a] = log_add(den[a], at.vals[r] + E.w_beta1 + prefix_log(a) + st.prS[ys + 1]);
        dcx[a] = log_add(dcx[a], at.vals[r] + pr_ctx(a, ys));
      }
    }
    std::string ln = "[l2r-rho] a:rho_ext/rho_ctx/log_mass";
    for (int a = 0; a < n; ++a) {
      if (num[a] < HALF_NEG) continue;
      char b[96];
      std::snprintf(b, sizeof b, " %d:%.2f/%.2f/%.1f", a, num[a] - den[a], num[a] - dcx[a],
                    num[a] - logZ);
      ln += b;
    }
    std::fprintf(stderr, "%s\n", ln.c_str());
    // exterior suffix with pseudoknots (C_b[y]) against the pre-pass's nested Sbar[y+1]
    ln = "[l2r-suffix] y:log(C_b[y]/Sbar[y+1])";
    for (int y = 0; y < n; ++y) {
      char b[48];
      std::snprintf(b, sizeof b, " %d:%.2f", y, C_b[y] - st.prS[y + 1]);
      ln += b;
    }
    std::fprintf(stderr, "%s\n", ln.c_str());
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) {
      std::fprintf(stderr, "--l2r-oracle-dump: cannot write %s\n", path.c_str());
      return;
    }
    const int32_t nn = n;
    std::fwrite(&nn, sizeof nn, 1, f);
    for (int a = 0; a < n; ++a) {
      const double r = num[a] > HALF_NEG && den[a] > HALF_NEG ? num[a] - den[a] : 0.0;
      std::fwrite(&r, sizeof r, 1, f);
    }
    for (int seed = 0; seed < 2; ++seed) {
      std::vector<std::pair<uint64_t, double>> rec;
      for (int j = 0; j < n; ++j)
        for (const Grp& g : (seed ? LB.Gs[j] : LB.G[j]).grps)
          for (const Ent& e : g.ents)
            if (e.b > HALF_NEG) rec.emplace_back(key4(e.u, g.v, e.p, j), e.b);
      const int64_t c = static_cast<int64_t>(rec.size());
      std::fwrite(&c, sizeof c, 1, f);
      for (const auto& [k, b] : rec) {
        std::fwrite(&k, sizeof k, 1, f);
        std::fwrite(&b, sizeof b, 1, f);
      }
    }
    std::vector<std::pair<uint64_t, double>> rec;
    for (int j = 0; j < n; ++j)
      for (const Grp& g : LA.Gs[j].grps)
        for (const Ent& e : g.ents)
          if (e.b > HALF_NEG) rec.emplace_back(key4(e.u, g.v, e.p, j), e.b);
    const int64_t c = static_cast<int64_t>(rec.size());
    std::fwrite(&c, sizeof c, 1, f);
    for (const auto& [k, b] : rec) {
      std::fwrite(&k, sizeof k, 1, f);
      std::fwrite(&b, sizeof b, 1, f);
    }
    std::fclose(f);
  }
  bool l2r_oracle_load(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return false;
    int32_t nn = 0;
    bool ok = std::fread(&nn, sizeof nn, 1, f) == 1 && nn == n;
    lRho.assign(n, 0.0);
    if (ok) ok = std::fread(lRho.data(), sizeof(double), n, f) == static_cast<size_t>(n);
    for (int seed = 0; ok && seed < 2; ++seed) {
      int64_t c = 0;
      ok = std::fread(&c, sizeof c, 1, f) == 1;
      std::unordered_map<uint64_t, double>& m = seed ? lOUs : lOU;
      m.reserve(static_cast<size_t>(c));
      for (int64_t i = 0; ok && i < c; ++i) {
        uint64_t k;
        double b;
        ok = std::fread(&k, sizeof k, 1, f) == 1 && std::fread(&b, sizeof b, 1, f) == 1;
        m.emplace(k, b);
      }
    }
    int64_t c = 0;
    if (ok && std::fread(&c, sizeof c, 1, f) == 1) {
      lORs.reserve(static_cast<size_t>(c));
      for (int64_t i = 0; ok && i < c; ++i) {
        uint64_t k;
        double b;
        ok = std::fread(&k, sizeof k, 1, f) == 1 && std::fread(&b, sizeof b, 1, f) == 1;
        lORs.emplace(k, b);
      }
    } else if (lsub_e) {
      ok = false;
    }
    std::fclose(f);
    return ok;
  }
  // --l2r-sub e (oracle diagnostic): A seeds ranked by alpha + their exact outside (exact run's dump)
  bool lsub_e = false;
  std::unordered_map<uint64_t, double> lORs;
  double lors(int u, int v, int p, int j) const {
    auto it = lORs.find(key4(u, v, p, j));
    return it == lORs.end() ? NEG_INF : it->second;
  }
  // log sum_{y' > j} E_in V of arm r
  double lS(int r, int j) const {
    const int hi = prArmY0[r + 1];
    const int k = static_cast<int>(std::upper_bound(prYy.begin() + prArmY0[r], prYy.begin() + hi, j) -
                                   prYy.begin());
    return k < hi ? lYwS[k] : NEG_INF;
  }
  double lsc_L(const RN& r, int arm, int who) {
    const auto t0 = std::chrono::steady_clock::now();
    auto done = [&](double s) {
      const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      if (sx_k_exact) sx_kx_sec += dt;
      else sx_k_sec += dt;
      kl_in_read_sec += dt;
      return s;
    };
    if (arm < 0) return done(NEG_INF);
    double s = NEG_INF;
    auto note = [&](int m, int list) {
      if (sx_k_exact) return;
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

  // --l2r-diag: no pruning; every state's score kept per (table, step) in entry order, and
  // after the outside pass the posterior mass a top-b cut by score drops against the cut by
  // the posterior itself (oracle).
  bool l2r_diag = false;
  std::vector<std::vector<double>> ldS[PT_N];
  void ldiag_put(int pt, int j, std::vector<double>&& sc) {
    if (ldS[pt].empty()) ldS[pt].resize(n);
    ldS[pt][j] = std::move(sc);
  }

  // --l2r-trace FILE a',xs,x' (0-based): every state of one pseudoknot path at its step (A
  // tables with arm (xs, x'), B tables with (xs, a'), S with (a', xs, x')); scores only where
  // the beam prunes anyway, so the run is unchanged.  After the outside pass: a, b, posterior.
  FILE* ltr_f = nullptr;
  int ltr_ap = -1, ltr_xs = -1, ltr_xp = -1;
  static const char* ltr_name(int pt) {
    static const char* pn[PT_N] = {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us", "Uls", "Ulss",
                                   "Urs", "S", "P1", "P2", "P3"};
    return pn[pt];
  }
  bool ltr_match4(int pt, const Grp& g, const Ent& e) const {
    return pt < PT_U ? (g.v == ltr_xp && e.u == ltr_xs) : (g.v == ltr_xs && e.u == ltr_ap);
  }
  bool ltr_match3(uint64_t key) const {
    return k3a(key) == ltr_ap && k3b(key) == ltr_xs && k3c(key) == ltr_xp;
  }
  void ltr_line(int pt, int j, int p, int u, int v, double a, const std::vector<double>* sc,
                const std::vector<char>* keep, size_t k, size_t tot, int b) {
    std::fprintf(ltr_f, "scan %s j %d p %d u %d v %d a %.6f", ltr_name(pt), j, p, u, v, a);
    if (sc) {
      const double s = (*sc)[k];
      long long rank = 1;
      for (double x : *sc) rank += x > s;
      std::vector<double> c(*sc);
      std::nth_element(c.begin(), c.begin() + (b - 1), c.end(), std::greater<double>());
      std::fprintf(ltr_f, " s %.6f rank %lld tot %zu b %d cut %.6f kept %d\n", s, rank, tot, b,
                   c[b - 1], static_cast<int>((*keep)[k]));
    } else {
      std::fprintf(ltr_f, " tot %zu b %d noprune\n", tot, b);
    }
  }
  void ltr_step4(const Tab4& t, int pt, int j, const std::vector<double>* sc,
                 const std::vector<char>* keep, size_t tot, int b) {
    size_t k = 0;
    for (const Grp& g : t.grps)
      for (const Ent& e : g.ents) {
        if (ltr_match4(pt, g, e)) ltr_line(pt, j, e.p, e.u, g.v, e.a, sc, keep, k, tot, b);
        ++k;
      }
  }
  void ltr_final() {
    if (!ltr_f) return;
    std::fprintf(ltr_f, "arm %d,%d idx %d logZ %.9f\n", ltr_xs, ltr_xp, larm(ltr_xs, ltr_xp), logZ);
    for (int j = 0; j < n; ++j) {
      for (int pt = 0; pt < PT_S; ++pt)
        for (const Grp& g : l2r_tab(pt, j).grps)
          for (const Ent& e : g.ents)
            if (ltr_match4(pt, g, e))
              std::fprintf(ltr_f, "final %s j %d p %d u %d v %d a %.6f b %.6f post %.6e\n",
                           ltr_name(pt), j, e.p, e.u, g.v, e.a, e.b,
                           e.a > HALF_NEG && e.b > HALF_NEG ? std::exp(e.a + e.b - logZ) : 0.0);
      if (j == ltr_c1)
        for (const Grp& g : LA.Gs[j].grps)
          for (const Ent& e : g.ents)
            if (e.p == ltr_a2)
              std::fprintf(ltr_f, "final-seed j %d p %d u %d v %d arm %d a %.6f b %.6f post %.6e\n", j, e.p, e.u,
                           g.v, larm(e.u, g.v), e.a, e.b,
                           e.a > HALF_NEG && e.b > HALF_NEG ? std::exp(e.a + e.b - logZ) : 0.0);
      const SMap& m = LS[j];
      for (size_t k = 0; k < m.keys.size(); ++k)
        if (ltr_match3(m.keys[k]))
          std::fprintf(ltr_f, "final S j %d p %d u %d v %d a %.6f b %.6f post %.6e\n", j,
                       k3a(m.keys[k]), k3b(m.keys[k]), k3c(m.keys[k]), m.v[k].F, m.v[k].bF,
                       m.v[k].F > HALF_NEG && m.v[k].bF > HALF_NEG
                           ? std::exp(m.v[k].F + m.v[k].bF - logZ) : 0.0);
    }
    std::fclose(ltr_f);
    ltr_f = nullptr;
  }

  // --l2r-perf on: spawn-stage counters and wall-clock timers, reported on stderr ([l2r-perf]);
  // nothing they touch feeds a value or a pruning decision.
  bool lpf = false;
  long long lpf_npos[PT_N] = {}, lpf_full[PT_N] = {}, lpf_nsc[PT_N] = {}, lpf_premax[PT_N] = {};
  double lpf_prune_sec[PT_N] = {};
  // filler (l2r_fext) per table 0 S, 1 P1, 2 P3: items d.j, state x item visits, unpaired visits
  long long lpf_fx_pos[3] = {}, lpf_fx_items[3] = {}, lpf_fx_pair[3] = {}, lpf_fx_unp[3] = {},
            lpf_fx_pairmax[3] = {};
  double lpf_fx_sec[3] = {};
  // lp_score per table 0 P1, 1 P2, 2 P3
  long long lpf_sc_calls[3] = {}, lpf_sc_hsteps[3] = {}, lpf_sc_dot[3] = {}, lpf_sc_kcalls[3] = {},
            lpf_sc_ksteps[3] = {}, lpf_sc_osteps[3] = {};
  double lpf_sc_hsec[3] = {}, lpf_sc_ksec[3] = {}, lpf_sc_osec[3] = {};
  // H construction
  long long lpf_h_merge = 0, lpf_h_trunc = 0, lpf_h_trunc_in = 0, lpf_h_newroots = 0, lpf_h_lam = 0,
            lpf_h_helpers = 0, lpf_h_clears = 0;
  double lpf_h_root_sec = 0.0, lpf_h_merge_sec = 0.0, lpf_h_help_sec = 0.0;
  // growth helpers and per-position overheads
  long long lpf_hg_clears = 0, lpf_help_built = 0, lpf_help_rebuilt = 0, lpf_racc_dense = 0,
            lpf_racc_fin = 0, lpf_racc_n = 0;
  double lpf_hg_sec = 0.0, lpf_mtrim_sec = 0.0, lpf_help_sec = 0.0;
  long long lpf_ml_cand = 0, lpf_ml_max = 0, lpf_ml_nodes = 0, lpf_stc_n = 0, lpf_stc_buckets = 0,
            lpf_stc_scan = 0;
  double lpf_stc_sec = 0.0;
  // init (V build): nodes (outer, inner), pr_steps' interior / multiloop steps, multiloop
  // candidates (helper terms) per node, entries added (interior children, helpers, inside the
  // helpers), entries after the finish / after the trim, arm entries, partner pairs
  long long vi_nodes[2] = {}, vi_nb_int = 0, vi_nb_ml = 0, vi_mlc = 0, vi_mlc_max = 0, vi_ml_nodes = 0,
            vi_add_int = 0, vi_add_ml = 0, vi_help_adds = 0, vi_ent_fin = 0, vi_ent_kept = 0, vi_ats = 0,
            vi_pairs = 0;
  double vi_steps_s = 0.0, vi_terms_s = 0.0, vi_add_int_s = 0.0, vi_add_ml_s = 0.0, vi_help_s = 0.0,
         vi_fin_s = 0.0, vi_ats_s = 0.0, vi_vbuild_s = 0.0, vi_arm_s = 0.0, vi_rest_s = 0.0, vi_pairs_s = 0.0;
  static void vi_lap(std::chrono::steady_clock::time_point& t0, double& acc) {
    const auto t1 = std::chrono::steady_clock::now();
    acc += std::chrono::duration<double>(t1 - t0).count();
    t0 = t1;
  }
  std::unordered_set<uint64_t> lpf_help_keys;
  static long long lpf_log2(size_t m) {
    long long s = 0;
    while (m > 0) {
      ++s;
      m >>= 1;
    }
    return s;
  }
  void lpf_prune_stat(int pt, size_t pre, size_t kept, int b, bool scored) {
    ++lpf_npos[pt];
    lpf_premax[pt] = std::max(lpf_premax[pt], static_cast<long long>(pre));
    if (b > 0 && kept >= static_cast<size_t>(b)) ++lpf_full[pt];
    if (scored) lpf_nsc[pt] += static_cast<long long>(pre);
  }

  // Top l2r_beam entries of t by score(group, entry).
  template <class Score>
