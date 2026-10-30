// Pass 1 products for pass 3 (from the nested-only pre-pass at beam b1): build_stem_est
// (S, C0), nested outside beta, chain data, and build_principled (prS, prQz, PrMbar: the
// nested outside Sbar_pre / C_pre that the pass-2 reads multiply).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  // --stem-est tables from this finished nested-only run (StemEst above).
  void build_stem_est(int L, StemEst& st, bool mirror = false) const {
    if (!mates_ready) const_cast<DP03*>(this)->lp_build_mates();
    st.n = n;
    st.loop = L;
    st.C0 = C;
    std::vector<double>& S = st.S;
    S.assign(n + 1, NEG_INF);
    S[n] = 0.0;
    if (dz) {
      std::vector<double> P(n + 1, NEG_INF);
      P[n] = 0.0;
      for (int y = n - 1; y >= 0; --y) {
        double acc = NEG_INF;
        for (int j = y + 1; j < n; ++j) {
          double qb = Qb[ix(y, j)];
          if (qb > HALF_NEG && E.end_ok(y, j)) acc = log_add(acc, qb + E.at(y, j) + S[j + 1]);
        }
        P[y] = acc;
        double s = P[y];
        for (int t = y; t < n; ++t) s = log_add(s, E.run(y, t) + P[t + 1]);
        S[y] = s;
      }
    } else {
      for (int y = n - 1; y >= 0; --y) {
        double acc = S[y + 1];
        for (int j = y + 1; j < n; ++j) {
          double qb = Qb[ix(y, j)];
          if (qb > HALF_NEG && E.end_ok(y, j)) acc = log_add(acc, qb + E.at(y, j) + S[j + 1]);
        }
        S[y] = acc;
      }
    }
    const int lc = std::min(L, cap);
    std::vector<double> T(static_cast<size_t>(n) * n, NEG_INF);
    std::vector<stemidx::Entry> ents;
    for (int x = 0; x < n; ++x)
      for (int y = next_mate(x, x); y >= 0; y = next_mate(x, y)) {
        double acc = E.end_ok(x, y) ? E.at(x, y) + S[y + 1] : NEG_INF;
        for (int x0 = std::max(0, x - 1 - lc); x0 < x; ++x0)
          for (int y0 = y + 1; y0 <= std::min(n - 1, y + 1 + lc); ++y0) {
            double t0 = T[ix(x0, y0)];
            if (t0 < HALF_NEG) continue;
            double w = E.stem_interior(x0, x, y, y0);
            if (w < HALF_NEG) continue;
            acc = log_add(acc, w + t0);
          }
        T[ix(x, y)] = acc;
        if (acc > HALF_NEG && E.end_ok(x, y)) ents.push_back({x, y, (E.at(x, y) + acc) + C[y - 1]});
      }
    st.dropped = stem_sparse(ents);
    st.entries = ents.size();
    st.index = stemidx::PersistentSum(n, ents);
    st.mirror = mirror;
    if (!mirror) return;
    std::fill(T.begin(), T.end(), NEG_INF);
    std::vector<stemidx::Entry> ents2;
    for (int x = 0; x < n; ++x) {
      const double pre = x > 0 ? C[x - 1] : 0.0;
      for (int y = next_mate(x, x); y >= 0; y = next_mate(x, y)) {
        double acc = E.end_ok(x, y) ? E.at(x, y) + pre : NEG_INF;
        for (int x0 = std::max(0, x - 1 - lc); x0 < x; ++x0)
          for (int y0 = y + 1; y0 <= std::min(n - 1, y + 1 + lc); ++y0) {
            double t0 = T[ix(x0, y0)];
            if (t0 < HALF_NEG) continue;
            double w = E.stem_interior(x0, x, y, y0);
            if (w < HALF_NEG) continue;
            acc = log_add(acc, w + t0);
          }
        T[ix(x, y)] = acc;
        if (acc > HALF_NEG && E.end_ok(x, y))
          ents2.push_back({n - 1 - y, n - 1 - x, (E.at(x, y) + acc) + S[x + 1]});
      }
    }
    st.dropped2 = stem_sparse(ents2);
    st.entries2 = ents2.size();
    st.index2 = stemidx::PersistentSum(n, ents2);
  }

  // st.beta from this nested-only run (after build_stem_est, which sets st.S): runs the outside
  // pass and keeps the Qb outside value of every pair.
  void build_nested_beta(StemEst& st) {
    if (!mates_ready) lp_build_mates();
    outside();
    st.beta.assign(static_cast<size_t>(n) * n, NEG_INF);
    for (int x = 0; x < n; ++x)
      for (int y = next_mate(x, x); y >= 0; y = next_mate(x, y)) {
        double b = Qb_b[ix(x, y)];
        if (b < HALF_NEG && E.end_ok(x, y)) b = prefix_log(x) + E.at(x, y) + st.S[y + 1];
        st.beta[ix(x, y)] = b;
      }
  }

  // Nested-only run: also fill Qmp (stem_ml), the chain's mbar.  Set before pfunc().
  bool chain_qmp = false;

  // st's --stem-est-chain data from this nested-only run: runs the outside pass; kept pairs
  // have (Qb + Qb outside) - logZ >= log(floor).
  void build_chain(StemEst& st, double floor) {
    outside();
    const size_t nn = static_cast<size_t>(n) * n;
    st.chain = true;
    st.kept.assign(nn, 0);
    st.kept_right.assign(n, {});
    st.mbar.assign(nn, NEG_INF);
    st.n_kept = 0;
    const double lf = floor > 0.0 ? std::log(floor) : NEG_INF;
    for (int j = 0; j < n; ++j) {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        const int i = at.keys[r];
        if (i < 0) continue;
        const double b = Qb_b[ix(i, j)];
        if (at.vals[r] < HALF_NEG || b < HALF_NEG) continue;
        if (at.vals[r] + b - logZ < lf) continue;
        st.kept[ix(i, j)] = 1;
        st.kept_right[i].push_back(j);
        ++st.n_kept;
      }
    }
    for (std::vector<int>& kr : st.kept_right) std::sort(kr.begin(), kr.end());
    for (int j = 0; j < n; ++j) {
      const OMap& at = qmp_at(j);
      for (size_t r = 0; r < at.slots(); ++r)
        if (at.keys[r] >= 0) st.mbar[ix(at.keys[r], j)] = at.vals[r];
    }
  }

  // Closed-value outsides of Qm / Qm2 with dangles (before run_back_at); set before outside().
  bool keep_closed_b = false;
  SpMat QmC_b, Qm2C_b;

  // st's --stem-est-principled data from this nested-only run (chain_qmp set; runs the outside
  // pass): prS, prQz, prMbar, kept pairs (P >= floor), prBml; nested: the outsides of Qb / Qm /
  // Qm2 (closed values).  Bml[d][e] = (sum of
  //   (a) Qm_b[i][e] w3^(d-i) run(i, d-1)          over cells (i, e), d - i <= ml_left,
  //   (b) Qm_b[i][e] Qm[i][d-1]                     over cells (i, d-1),
  // times alpha2^ml_pk_k: rev_branch's pseudoknot branch) + (c) the PKM closing of every pair
  // i.j' (i < d, j' > e): Qb_b[i][j'] alpha1 alpha2^(1+ml_pk_k) at_close(i, j') w3^(j'-1-e)
  // run(e+1, j'-1) [w3^(d-i-1) run(i+1, d-1) + Qm[i+1][d-1]] (fill_Qb's PKM terms); beta1m
  // is left to Ctx.  check: every cell recomputed directly (n <= 200).
  void build_principled(StemEst& st, double floor, bool nested, bool check,
                        double partner_floor = -1.0) {
    outside();
    const size_t nn = static_cast<size_t>(n) * n;
    st.n = n;
    st.pr = true;
    std::vector<double>& S = st.prS;
    S.assign(n + 1, NEG_INF);
    S[n] = 0.0;
    if (dz) {
      std::vector<double> P(n + 1, NEG_INF);
      P[n] = 0.0;
      for (int y = n - 1; y >= 0; --y) {
        double acc = NEG_INF;
        for (int j = y + 1; j < n; ++j) {
          double qb = Qb[ix(y, j)];
          if (qb > HALF_NEG && E.end_ok(y, j)) acc = log_add(acc, qb + E.at(y, j) + S[j + 1]);
        }
        P[y] = acc;
        double s = P[y];
        for (int t = y; t < n; ++t) s = log_add(s, E.run(y, t) + P[t + 1]);
        S[y] = s;
      }
    } else {
      for (int y = n - 1; y >= 0; --y) {
        double acc = S[y + 1];
        for (int j = y + 1; j < n; ++j) {
          double qb = Qb[ix(y, j)];
          if (qb > HALF_NEG && E.end_ok(y, j)) acc = log_add(acc, qb + E.at(y, j) + S[j + 1]);
        }
        S[y] = acc;
      }
    }
    // --qzbar-old: always dense-copy nested Qz into prQz (snap/incperf), even if qz_beam_fill.
    if (!qz_beam_fill || qzbar_old) {
      st.prQz.assign(nn, NEG_INF);
      for (int i = 0; i < n; ++i)
        for (int j = i; j < n; ++j) st.prQz[ix(i, j)] = Qz[ix(i, j)];
      if (qzbar_old) {
        st.qz_on = 0;
        st.qz_col.clear();
        st.qz_col.shrink_to_fit();
      }
    }
    if (!st.skip_qbin) st.prQbIn = Qb.dense();
    st.prMbar.assign(nn, NEG_INF);
    if (p2_sweep) {
      st.prMbarAt.assign(n, {});
      st.prMbarRight.assign(n, {});
    }
    for (int j = 0; j < n; ++j) {
      const OMap& at = qmp_at(j);
      for (size_t r = 0; r < at.slots(); ++r) {
        const int i = at.keys[r];
        if (i < 0) continue;
        st.prMbar[ix(i, j)] = at.vals[r];
        if (p2_sweep) {
          st.prMbarAt[j].push_back(i);
          st.prMbarRight[i].push_back(j);
        }
      }
    }
    if (p2_sweep) {
      st.pr_mbar_n = 0;
      st.pr_mbar_max = 0;
      for (int j = 0; j < n; ++j) {
        st.pr_mbar_n += static_cast<long long>(st.prMbarAt[j].size());
        st.pr_mbar_max = std::max(st.pr_mbar_max, static_cast<int>(st.prMbarAt[j].size()));
        std::sort(st.prMbarRight[j].begin(), st.prMbarRight[j].end());
      }
      if (beam > 0) {
        auto stat = [&](const std::vector<OMap>& tabs) {
          long long c = 0;
          int m = 0;
          for (int j = 0; j < n; ++j) {
            int k = 0;
            for (size_t r = 0; r < tabs[j].slots(); ++r)
              if (tabs[j].keys[r] >= 0) ++k;
            c += k;
            m = std::max(m, k);
          }
          return std::pair<long long, int>(c, m);
        };
        const auto qb = stat(Qb_at), qm = stat(Qm_at), q2 = stat(Qm2_at);
        std::fprintf(stderr,
                     "[pre-beam] b %d key %s Qb %lld max %d Qm %lld max %d Qm2 %lld max %d\n",
                     beam, pre_key_inside ? "inside" : "prefix", qb.first, qb.second, qm.first,
                     qm.second, q2.first, q2.second);
      }
    }
    if (prof_on && !(p2_sweep && beam > 0)) {
      int qm_max = 0;
      for (int j = 0; j < n; ++j) {
        int k = 0;
        const OMap& at = qmp_at(j);
        for (size_t r = 0; r < at.slots(); ++r)
          if (at.keys[r] >= 0) ++k;
        qm_max = std::max(qm_max, k);
      }
      std::fprintf(stderr, "[qm-col] max %d\n", qm_max);
    }
    if (!st.skip_kept) st.prKept.assign(nn, 0);
    st.prKeptRight.assign(n, {});
    st.pr_n_kept = 0;
    const double lf = floor > 0.0 ? std::log(floor) : NEG_INF;
    st.pr_partner = partner_floor >= 0.0;
    const double lfp = partner_floor > 0.0 ? std::log(partner_floor) : NEG_INF;
    st.prKeptRightP.assign(st.pr_partner ? n : 0, {});
    st.pr_n_partner = 0;
    if (st.want_lp) st.prLP.assign(nn, NEG_INF);
    for (int j = 0; j < n; ++j) {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        const int i = at.keys[r];
        if (i < 0) continue;
        const double b = Qb_b[ix(i, j)];
        if (at.vals[r] < HALF_NEG || b < HALF_NEG) continue;
        if (st.want_lp) st.prLP[ix(i, j)] = at.vals[r] + b - logZ;
        if (st.pr_partner && at.vals[r] + b - logZ >= lfp) {
          st.prKeptRightP[i].push_back(j);
          ++st.pr_n_partner;
        }
        if (at.vals[r] + b - logZ < lf) continue;
        if (!st.skip_kept) st.prKept[ix(i, j)] = 1;
        st.prKeptRight[i].push_back(j);
        ++st.pr_n_kept;
      }
    }
    for (std::vector<int>& kr : st.prKeptRight) std::sort(kr.begin(), kr.end());
    for (std::vector<int>& kr : st.prKeptRightP) std::sort(kr.begin(), kr.end());
    st.prBml.clear();
    if (g_skip_bml) {
      if (nested) {
        st.pr_nested = true;
        st.prQb_b = Qb_b.dense();
        st.prQm_b = dz ? QmC_b.dense() : Qm_b.dense();
        st.prQm2_b = dz ? Qm2C_b.dense() : Qm2_b.dense();
      }
      pass2_apply(st);
      return;
    }

    const double w2 = E.w_alpha2, w3 = E.w_alpha3, kw = E.ml_pk_k * w2;
    std::vector<double> Bab(nn, NEG_INF);
    for (int e = 0; e < n; ++e) {
      const OMap& amj = Qm_at[e];
      for (size_t r = 0; r < amj.slots(); ++r) {
        const int i = amj.keys[r];
        if (i < 0) continue;
        const double mb = Qm_b[ix(i, e)];
        if (mb < HALF_NEG) continue;
        for (int d = i; d <= e && d - i <= E.ml_left && lead_gap_ok(d - i); ++d) {
          double w = w3 * (d - i);
          if (dz) w += E.run(i, d - 1);
          Bab[ix(d, e)] = log_add(Bab[ix(d, e)], mb + w);
        }
      }
    }
    for (int d = 1; d < n; ++d) {
      const OMap& amd = Qm_at[d - 1];
      for (size_t r = 0; r < amd.slots(); ++r) {
        const int i = amd.keys[r];
        if (i < 0) continue;
        const double qm = amd.vals[r];
        for (int e = d; e < n; ++e) {
          const double mb = Qm_b[ix(i, e)];
          if (mb > HALF_NEG) Bab[ix(d, e)] = log_add(Bab[ix(d, e)], mb + qm);
        }
      }
    }
    // (c): R[i][e] = sum over closings i.j' > e of their weight up to the PKM run end e
    std::vector<double> R(nn, NEG_INF);
    const double wc = E.w_alpha1 + (1.0 + E.ml_pk_k) * w2;
    for (int j = 0; j < n; ++j) {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        const int i = at.keys[r];
        if (i < 0 || !E.end_ok(i, j)) continue;
        const double bb = Qb_b[ix(i, j)];
        if (bb < HALF_NEG) continue;
        const double base = bb + wc + E.at_close(i, j);
        for (int e = i + 1; e <= j - 1; ++e) {
          double w = base + w3 * (j - 1 - e);
          if (dz) w += E.run(e + 1, j - 1);
          R[ix(i, e)] = log_add(R[ix(i, e)], w);
        }
      }
    }
    st.prBml.assign(nn, NEG_INF);
    std::vector<double> Lc(n, NEG_INF);
    for (int d = 1; d < n; ++d) {
      for (int i = 0; i < d; ++i) {
        double l = w3 * (d - i - 1);
        if (dz) l += E.run(i + 1, d - 1);
        if (i + 1 <= d - 1 && Qm[ix(i + 1, d - 1)] > HALF_NEG) l = log_add(l, Qm[ix(i + 1, d - 1)]);
        Lc[i] = l;
      }
      for (int e = d; e < n; ++e) {
        double c = NEG_INF;
        for (int i = 0; i < d; ++i) {
          const double rv = R[ix(i, e)];
          if (rv > HALF_NEG) c = log_add(c, rv + Lc[i]);
        }
        const double ab = Bab[ix(d, e)];
        st.prBml[ix(d, e)] = log_add(ab > HALF_NEG ? ab + kw : NEG_INF, c);
      }
    }
    if (nested) {
      st.pr_nested = true;
      st.prQb_b = Qb_b.dense();
      st.prQm_b = dz ? QmC_b.dense() : Qm_b.dense();
      st.prQm2_b = dz ? Qm2C_b.dense() : Qm2_b.dense();
    }
    if (check && n <= 200) bml_check(st);
    pass2_apply(st);
  }

  // Pass 2 limits. Every flag defaults off, and then this function returns without touching tables.
  // Q^z_pre(k,l): DP03 Q^z split on the last base, Q^p dropped. Kept starts ranked by
  // C_pre(k-1) * Q^z_pre (key1) or that times Pr(k-1 is a pair right end) (key2).
  // Pair sets: P1 = Q^b states that survived pass 1; own:b_P = top b_P left ends per right
  // end by the inside Q^b.
  // Outside push of gamma (growth) and mu (V) along the nested Qb / Qm / Qm2 hyperedges.
  // Pseudoknots are off on this pre-pass. Right ends run n-1..0. At j the order is the
  // reverse of the inside: extension from j+1, interior from kept outer pairs, the exterior
  // edge, the Qm/Qm2 branch split over kept states, then the pair pushes into Qm2.
  // Each map is cut before it pushes. Qm and Qm2 maps are freed at the end of j.
  struct LeadEnt {
    int x0;
    double w;
  };
