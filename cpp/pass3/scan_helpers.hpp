// Scan helpers: pseudoknot suffix Sbar_pk, items, per-pass profile and counters, beam caps,
// table views, b-collection; scan_ct_dump / profile_dump.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  std::vector<double> l2r_pk_suffix(const QpItems& qp, const QpItems* qb = nullptr) const {
    const StemEst& st = *stem_;
    std::vector<double> S2(n + 1, NEG_INF);
    S2[n] = 0.0;
    auto paired = [&](int y) {
      double acc = NEG_INF;
      if (qb) {
        for (const auto& it : (*qb)[y])
          if (E.end_ok(y, it.first)) acc = log_add(acc, it.second + E.at(y, it.first) + S2[it.first + 1]);
      } else {
        for (int j = y + 1; j < n; ++j) {
          const double q = st.QbIn(y, j);
          if (q > HALF_NEG && E.end_ok(y, j)) acc = log_add(acc, q + E.at(y, j) + S2[j + 1]);
        }
      }
      for (const auto& it : qp[y]) acc = log_add(acc, it.second + E.w_beta1 + S2[it.first + 1]);
      return acc;
    };
    if (dz) {
      std::vector<double> P(n + 1, NEG_INF);
      P[n] = 0.0;
      for (int y = n - 1; y >= 0; --y) {
        P[y] = paired(y);
        double s = P[y];
        for (int t = y; t < n; ++t) s = log_add(s, E.run(y, t) + P[t + 1]);
        S2[y] = s;
      }
    } else {
      for (int y = n - 1; y >= 0; --y) S2[y] = log_add(S2[y + 1], paired(y));
    }
    return S2;
  }
  // --l2r-p1 R: Sbar_pk from a finished beamed pseudoknot pass src: its stored Qp items with
  // the pre-pass nested Qb (Sn) or with src's own Qb items (Sp; Sp[0] = src's log Z, fill_C).
  static QpItems l2r_items(const std::vector<OMap>& at, int n, long long& cnt) {
    QpItems it(n);
    cnt = 0;
    for (int j = 0; j < n; ++j)
      for (size_t r = 0; r < at[j].slots(); ++r)
        if (at[j].keys[r] >= 0 && at[j].vals[r] > HALF_NEG) {
          it[at[j].keys[r]].emplace_back(j, at[j].vals[r]);
          ++cnt;
        }
    return it;
  }
  void l2r_sfx_pass(const DP03& src, std::vector<double>& Sn, std::vector<double>& Sp,
                    long long& qp_items, long long& qb_items) const {
    const QpItems qp = l2r_items(src.Qp_at, n, qp_items);
    const QpItems qb = l2r_items(src.Qb_at, n, qb_items);
    Sn = l2r_pk_suffix(qp);
    Sp = l2r_pk_suffix(qp, &qb);
  }
  // One pass's time split (s): init (l2r_init_scores: arm entries + prYw + prSF), growth and
  // selection (counted where they run), Rs = spawn (build + score), all tables' score and build
  // (Rs included), close; other = pass CPU - init - build - score - close (nested tables).
  std::string l2r_prof(double pass_cpu) const {
    double bs = 0.0, ps = 0.0;
    for (int t = 0; t < PT_N; ++t) {
      bs += l2r_bsec[t];
      ps += l2r_psec[t];
    }
    char b[512];
    std::snprintf(b, sizeof b,
                  "cpu %.2f init %.2f static %.2f grow %.2f sel %.2f spawn %.2f (cand %lld) "
                  "score %.2f build %.2f close %.2f other %.2f",
                  pass_cpu, l2r_isec, pr_static_sec, pr_grow_sec, pr_sel_sec,
                  l2r_bsec[PT_RS] + l2r_psec[PT_RS], sp_cand, ps, bs, l2r_csec,
                  pass_cpu - l2r_isec - bs - ps - l2r_csec);
    return b;
  }
  static int pct_p99(const long long* h, int nbin) {
    long long tot = 0;
    for (int i = 0; i < nbin; ++i) tot += h[i];
    if (tot <= 0) return -1;
    const long long need = (tot * 99 + 99) / 100;
    long long c = 0;
    for (int i = 0; i < nbin; ++i) {
      c += h[i];
      if (c >= need) return i;
    }
    return nbin - 1;
  }
  void scan_ct_dump() const {
    auto grp = [&](int a, int b, long long& pre, long long& kept, long long& mx, long long& bind,
                   long long& pos, int& beam) {
      pre = kept = mx = bind = pos = 0;
      beam = ltb(a);
      for (int t = a; t <= b; ++t) {
        pre += l2r_pre[t];
        kept += l2r_kept[t];
        mx = std::max(mx, l2r_kmax[t]);
        bind += l2r_bind[t];
        pos += l2r_pos[t];
      }
    };
    auto one = [&](const char* name, int a, int b) {
      long long pre, kept, mx, bind, pos;
      int beam;
      grp(a, b, pre, kept, mx, bind, pos, beam);
      const double per = pos > 0 ? static_cast<double>(pre) / static_cast<double>(pos) : 0.0;
      std::fprintf(stderr,
                   " %s beam %d pre %lld kept %lld max %lld bind %lld pos %lld cand_per_pos %.3f",
                   name, beam, pre, kept, mx, bind, pos, per);
    };
    long long m_n = 0, m_sum = 0, list_n = 0, list_sum = 0;
    for (int i = 0; i < 8; ++i) {
      m_n += sx_m_h[i];
      m_sum += static_cast<long long>(i) * sx_m_h[i];
    }
    for (int i = 0; i < 65; ++i) {
      list_n += sx_list_h[i];
      list_sum += static_cast<long long>(i) * sx_list_h[i];
    }
    const double k_per = sx_rkey_n > 0 ? static_cast<double>(sx_k_src[2]) / static_cast<double>(sx_rkey_n)
                                       : 0.0;
    long long rpre = 0;
    for (int t = PT_R; t <= PT_RRS; ++t) rpre += l2r_pre[t];
    const double rkey_per = rpre > 0 ? static_cast<double>(sx_rkey_n) / static_cast<double>(rpre) : 0.0;
    std::fprintf(stderr, "[scan-ct] n %d", n);
    one("P1", PT_P1, PT_P1);
    one("P2", PT_P2, PT_P2);
    one("P3", PT_P3, PT_P3);
    one("R", PT_R, PT_RRS);
    one("S", PT_S, PT_S);
    one("U", PT_U, PT_URS);
    std::fprintf(stderr,
                 " rkey_calls %lld r_cand %lld rkey_per_cand %.6f k_rkey %lld k_per_rkey %.3f "
                 "k_root %lld k_spawn %lld k_n %lld m_n %lld m_mean %.3f m_p99 %d "
                 "qlist_n %lld qlist_mean %.3f qlist_p99 %d\n",
                 sx_rkey_n, rpre, rkey_per, sx_k_src[2], k_per, sx_k_src[1], sx_k_src[3], sx_k_n,
                 m_n, m_n > 0 ? static_cast<double>(m_sum) / static_cast<double>(m_n) : 0.0,
                 pct_p99(sx_m_h, 8), list_n,
                 list_n > 0 ? static_cast<double>(list_sum) / static_cast<double>(list_n) : 0.0,
                 pct_p99(sx_list_h, 65));
    if (k_memo) {
      const double bytes =
          static_cast<double>(lKm.size()) * static_cast<double>(sizeof(uint64_t) + sizeof(double)) +
          static_cast<double>(lKm.bucket_count()) * static_cast<double>(sizeof(void*));
      std::fprintf(stderr,
                   "[k-memo] calls %lld hits %lld hit_rate %.6f entries %zu bytes %.0f k_s %.4f\n",
                   lKm_calls, lKm_hits,
                   lKm_calls > 0 ? static_cast<double>(lKm_hits) / static_cast<double>(lKm_calls) : 0.0,
                   lKm.size(), bytes, sx_k_sec);
    }
    if (k_exact_ub || sx_k_exact) {
      const double bytes = static_cast<double>(lKxM.size()) *
                               static_cast<double>(sizeof(uint64_t) + sizeof(double)) +
                           static_cast<double>(lKxM.bucket_count()) * static_cast<double>(sizeof(void*));
      const double per = n > 0 ? static_cast<double>(sx_kx_sigma) / static_cast<double>(n) : 0.0;
      std::fprintf(stderr,
                   "[k-exact] sigma %lld per_pos %.3f kx_calls %lld kx_hits %lld hit_rate %.6f "
                   "entries %zu bytes %.0f kx_s %.4f cols %lld mis %lld\n",
                   sx_kx_sigma, per, lKx_calls, lKx_hits,
                   lKx_calls > 0 ? static_cast<double>(lKx_hits) / static_cast<double>(lKx_calls) : 0.0,
                   lKxM.size(), bytes, sx_kx_sec, k_exact_cols, k_exact_mis);
    }
  }
  void profile_dump(const char* sid, double wall, double qz_pre) const {
    if (!prof_on) return;
    auto bucket = [&](int a, int b) {
      double s = 0.0;
      for (int t = a; t <= b; ++t) s += l2r_bsec[t] + l2r_psec[t];
      return s;
    };
    double rscore = 0.0;
    for (int t = PT_R; t <= PT_RRS; ++t) rscore += l2r_psec[t];
    const double rkey = std::max(0.0, rscore - prof_grow_reg[1]);
    long long pairs = 0;
    if (p2_g_mode > 0 && stem_)
      for (const std::vector<int>& kr : stem_->prKeptRightG)
        pairs += static_cast<long long>(kr.size());
    std::fprintf(stderr,
                 "[profile] %s wall %.3f qz_pre %.4f qz_scan %.4f v_s %.4f grow_s %.4f "
                 "grow_R %.4f grow_spawn %.4f grow_H %.4f grow_other %.4f rkey_s %.4f "
                 "P1 %.4f P2 %.4f P3 %.4f R %.4f U %.4f S %.4f close %.4f init %.4f "
                 "vectors %lld roots %lld closure %lld steps %lld steps_ml %lld help_steps %lld "
                 "pairs_P %lld\n",
                 sid, wall, qz_pre, prof_qz_sec, pr_static_sec, pr_grow_sec, prof_grow_reg[1],
                 prof_grow_reg[2], prof_grow_reg[3], prof_grow_reg[0], rkey, bucket(PT_P1, PT_P1),
                 bucket(PT_P2, PT_P2), bucket(PT_P3, PT_P3), bucket(PT_R, PT_RRS),
                 bucket(PT_U, PT_URS), bucket(PT_S, PT_S), l2r_csec, l2r_isec, rn_nodes + rn_roots,
                 rn_roots, static_cast<long long>(prRG.size()), st_kind[0][0] + st_kind[0][1],
                 st_kind[0][1], prof_help_steps, pairs);
    // Exclusive scan split. R and P3 above still include the growth built inside them.
    // Growth charged past a bucket is not removed twice; `over` is that excess and `rest`
    // is whatever of the scan none of these timers owns. v_s and init sit before the scan.
    const double P1 = bucket(PT_P1, PT_P1);
    const double P2 = bucket(PT_P2, PT_P2);
    const double P3b = bucket(PT_P3, PT_P3);
    const double Rb = bucket(PT_R, PT_RRS);
    const double Ub = bucket(PT_U, PT_URS);
    const double Sb = bucket(PT_S, PT_S);
    const double gR = std::min(prof_grow_reg[1], Rb);
    const double gH = std::min(prof_grow_reg[3], P3b);
    const double over = (prof_grow_reg[1] - gR) + (prof_grow_reg[3] - gH);
    const double rest = wall - prof_qz_sec - pr_grow_sec - P1 - P2 - (P3b - gH) - (Rb - gR) - Ub -
                        Sb - l2r_csec;
    std::fprintf(stderr,
                 "[profile-x] %s scan %.3f qz %.4f grow %.4f P1 %.4f P2 %.4f P3 %.4f R %.4f "
                 "U %.4f S %.4f close %.4f rest %.4f over %.4f\n",
                 sid, wall, prof_qz_sec, pr_grow_sec, P1, P2, P3b - gH, Rb - gR, Ub, Sb, l2r_csec,
                 rest, over);
    const double spawn = l2r_bsec[PT_RS];
    const double ext = l2r_bsec[PT_R] + l2r_bsec[PT_U] + l2r_bsec[PT_S];
    const double nest = l2r_bsec[PT_P1] + l2r_bsec[PT_P2] + l2r_bsec[PT_P3];
    const double km = sx_k_n > 0 ? static_cast<double>(sx_m_sum) / static_cast<double>(sx_k_n) : 0.0;
    std::fprintf(stderr,
                 "[profile-s] %s spawn %.4f ext %.4f nest %.4f k %.4f k_n %lld k_m %.3f rkey %.4f "
                 "per_nt %.6f\n",
                 sid, spawn, ext, nest, sx_k_sec, sx_k_n, km, sx_rkey_sec,
                 n > 0 ? wall / static_cast<double>(n) : 0.0);
  }
  // The pass-1 copy: every beam (nested tables and l2r) capped at B.
  void l2r_cap_beam(int B) {
    beam = beam > 0 ? std::min(beam, B) : B;
    for (int t = 0; t < BT_N; ++t)
      if (tbeam[t] <= 0 || tbeam[t] > B) tbeam[t] = B;
    l2r_beam = B;
  }
  // After st.prS is replaced: prYw (per arm and y', log sum_ys w e^beta1 Sbar[ys+1]) again.
  void l2r_yw_rebuild() {
    if (pr_vxs) throw std::runtime_error("--pr-v-key xs: no per-ys arm entries to rebuild prYw from");
    const StemEst& st = *stem_;
    std::fill(prYw.begin(), prYw.end(), NEG_INF);
    for (int r = 0; r < prA; ++r) {
      int g = prArmY0[r];
      for (int k = prArmE0[r]; k < prArmE0[r + 1]; ++k) {
        if (k > prArmE0[r] && prEy[k] != prEy[k - 1]) ++g;
        prYw[g] = log_add(prYw[g], prEw[k] + (E.w_beta1 + st.prS[prEys[k] + 1]));
      }
    }
    prSF.clear();
  }
  bool l2r_stem_is(const StemEst* s) const { return s == stem_; }
  int l2r_arm_entries() const { return static_cast<int>(l2r_arm_ents()); }

  // --l2r-split on (with --l2r-diag): Qp_b by consumer, QpS[0] exterior (beta1), [1] multiloop
  // (beta1m: Qm / Qmp branch, PK inside a pair via PKM), [2] gap (beta1p: Qz gap, S gap).  The
  // l2r states' outsides are linear in the Qp_b read at close B, so an outside pass reading
  // QpS[c] there gives every state's share by consumer (lspB[c]); a last normal pass restores b.
  bool l2r_split = false;
  bool lsp_rec = false;
  int lsp_use = -1;
  std::vector<double> QpS[3];
  std::vector<std::vector<double>> lspB[3][PT_N];
  void qps(int c, int i, int j, double v) {
    if (lsp_rec) add2b(QpS[c], i, j, v);
  }
  SMap& l2r_smap(int pt, int j) {
    return pt == PT_S ? LS[j] : pt == PT_P1 ? LP1[j] : pt == PT_P2 ? LP2[j] : LP3[j];
  }
  Tab4& l2r_tab(int pt, int j) {
    LFam& F = pt < PT_U ? LA : LB;
    const int r = pt - (pt < PT_U ? PT_R : PT_U);
    return r == 0 ? F.G[j] : r == 1 ? F.Gs[j] : r == 2 ? F.Gls[j] : r == 3 ? F.Glss[j] : F.Grs[j];
  }
  void l2r_reset_b() {
    for (int j = 0; j < n; ++j) {
      for (int pt = 0; pt < PT_S; ++pt)
        for (Grp& g : l2r_tab(pt, j).grps)
          for (Ent& e : g.ents) e.b = NEG_INF;
      for (SSt& s : LS[j].v) s.bP = s.bU1 = s.bU2 = s.bc = s.bF = NEG_INF;
      if (l2r_inc)
        for (std::vector<SMap>* T : {&LP1, &LP2, &LP3})
          for (SSt& s : (*T)[j].v) s.bP = s.bU1 = s.bU2 = s.bc = s.bF = NEG_INF;
    }
  }
  void l2r_collect_b(int c) {
    for (int pt = 0; pt < PT_N; ++pt) {
      lspB[c][pt].assign(n, {});
      for (int j = 0; j < n; ++j) {
        std::vector<double>& o = lspB[c][pt][j];
        if (pt >= PT_S) {
          if (pt > PT_S && !l2r_inc) continue;
          for (const SSt& s : l2r_smap(pt, j).v) o.push_back(s.bF);
        } else {
          for (const Grp& g : l2r_tab(pt, j).grps)
            for (const Ent& e : g.ents) o.push_back(e.b);
        }
      }
    }
  }
  void l2r_split_run() {
    lsp_rec = false;
    for (int c = 0; c < 3; ++c) {
      l2r_reset_b();
      outside_done = false;
      lsp_use = c;
      outside();
      l2r_collect_b(c);
    }
    lsp_use = -1;
    l2r_reset_b();
    outside_done = false;
    outside();
  }

