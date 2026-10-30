// Pass-3 ranking options: setters for --stem-est* / --pr-* (stem score, profile, chain,
// principled, arms, V options), Bml check, hr variants, cpu_stop_check, background reorder
// and --stem-bg-check.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  // build_principled's Bml cell by cell from the three edge lists (long double).
  void bml_check(StemEst& st) const {
    const long double w2 = E.w_alpha2, w3 = E.w_alpha3;
    for (int d = 0; d < n; ++d)
      for (int e = d; e < n; ++e) {
        std::vector<long double> t;
        for (int i = std::max(0, d - std::min(E.ml_left, n)); i <= d; ++i) {
          if (!lead_gap_ok(d - i)) continue;
          if (!Qm_at[e].has(i) || Qm_b[ix(i, e)] < HALF_NEG) continue;
          long double w = Qm_b[ix(i, e)] + E.ml_pk_k * w2 + w3 * (d - i);
          if (dz) w += E.run(i, d - 1);
          t.push_back(w);
        }
        if (d >= 1)
          for (int i = 0; i <= d - 1; ++i) {
            if (!Qm_at[d - 1].has(i) || Qm_b[ix(i, e)] < HALF_NEG) continue;
            t.push_back(static_cast<long double>(Qm_b[ix(i, e)]) + Qm[ix(i, d - 1)] +
                        E.ml_pk_k * w2);
          }
        for (int i = 0; i < d; ++i)
          for (int j = e + 1; j < n; ++j) {
            if (!Qb_at[j].has(i) || !E.end_ok(i, j) || Qb_b[ix(i, j)] < HALF_NEG) continue;
            long double base = static_cast<long double>(Qb_b[ix(i, j)]) + E.w_alpha1 +
                               (1.0L + E.ml_pk_k) * w2 + E.at_close(i, j) + w3 * (j - 1 - e);
            if (dz) base += E.run(e + 1, j - 1);
            long double l1 = w3 * (d - i - 1);
            if (dz) l1 += E.run(i + 1, d - 1);
            t.push_back(base + l1);
            if (i + 1 <= d - 1 && Qm[ix(i + 1, d - 1)] > HALF_NEG)
              t.push_back(base + Qm[ix(i + 1, d - 1)]);
          }
        const double got = st.Bml(d, e);
        ++st.bml_check_cells;
        if (t.empty()) {
          if (got > HALF_NEG) ++st.bml_check_mis;
          continue;
        }
        long double m = *std::max_element(t.begin(), t.end()), s = 0.0L;
        for (long double x : t) s += std::exp(x - m);
        const double ref = static_cast<double>(m + std::log(s));
        if (got < HALF_NEG) ++st.bml_check_mis;
        else st.bml_check_rel = std::max(st.bml_check_rel, std::fabs(std::expm1(got - ref)));
      }
  }

  // st's --stem-est-hr av data from this nested-only run (runs the outside pass):
  // phi(i, j) = (pu(i) + P(i, j)) (pu(j) + P(i, j)), pu(k) = max(0, 1 - sum_j P(k, j)),
  // raised to gamma; Tav(x, y) = phi [E_out(x, y) + w_stack(x-1, y+1; x, y) Tav(x-1, y+1)]
  // over the pairs that can form, outermost first; entries Uav(x', y') = E_in Tav for the
  // end_ok pairs (E_out = E_in = at if end_ok, else 0).
  void build_hr_av(StemEst& st, double gamma) {
    nested_lphi(st, gamma);
    const size_t nn = static_cast<size_t>(n) * n;
    std::vector<double> T(nn, NEG_INF);
    std::vector<stemidx::Entry> ents;
    for (int x = 0; x < n; ++x)
      for (int y : partners_right[x]) {
        const double lp = st.lphi[ix(x, y)];
        if (lp < HALF_NEG) continue;
        double acc = E.end_ok(x, y) ? E.at(x, y) : NEG_INF;
        if (x >= 1 && y + 1 <= n - 1) {
          const double t0 = T[ix(x - 1, y + 1)];
          if (t0 > HALF_NEG) {
            const double w = E.stem_interior(x - 1, x, y, y + 1);
            if (w > HALF_NEG) acc = log_add(acc, w + t0);
          }
        }
        if (acc < HALF_NEG) continue;
        T[ix(x, y)] = lp + acc;
        if (E.end_ok(x, y)) ents.push_back({x, y, E.at(x, y) + T[ix(x, y)]});
      }
    st.dropped_av = stem_sparse(ents);
    st.entries_av = ents.size();
    st.index_av = stemidx::PersistentSum(n, ents);
    st.hr_av = true;
  }

  // st's --stem-est-hr phi data (after build_stem_est on the stem-table run, whose S, C0 and
  // loop st keeps): the original stem table with phi_gamma(x, y) multiplied in once per pair,
  //   T'(x, y) = phi_gamma(x, y) [E_out(x, y) S[y+1] + sum w(x0, y0 | x, y) T'(x0, y0)],
  // entries E_in(x', y') T'(x', y') C0[y'-1], same loop order and index as build_stem_est.
  // gamma = 0 adds nothing, so the entries equal the original ones bit for bit.
  void build_hr_phi(StemEst& st, double gamma) {
    nested_lphi(st, gamma);
    const bool one = gamma == 0.0;
    if (one)
      for (int x = 0; x < n; ++x)
        for (int y : partners_right[x]) st.lphi[ix(x, y)] = 0.0;
    const std::vector<double>& S = st.S;
    const std::vector<double>& C0 = st.C0;
    const int lc = std::min(st.loop, cap);
    std::vector<double> T(static_cast<size_t>(n) * n, NEG_INF);
    std::vector<stemidx::Entry> ents;
    for (int x = 0; x < n; ++x)
      for (int y : partners_right[x]) {
        double acc = E.end_ok(x, y) ? E.at(x, y) + S[y + 1] : NEG_INF;
        for (int x0 = std::max(0, x - 1 - lc); x0 < x; ++x0)
          for (int y0 = y + 1; y0 <= std::min(n - 1, y + 1 + lc); ++y0) {
            double t0 = T[ix(x0, y0)];
            if (t0 < HALF_NEG) continue;
            double w = E.stem_interior(x0, x, y, y0);
            if (w < HALF_NEG) continue;
            acc = log_add(acc, w + t0);
          }
        const double t = one ? acc : st.lphi[ix(x, y)] + acc;
        T[ix(x, y)] = t;
        if (t > HALF_NEG && E.end_ok(x, y)) ents.push_back({x, y, (E.at(x, y) + t) + C0[y - 1]});
      }
    st.dropped_av = stem_sparse(ents);
    st.entries_av = ents.size();
    st.index_av = stemidx::PersistentSum(n, ents);
    st.hr_av = true;
    st.hr_phi = true;
  }

  // st.lphi = gamma (log(pu(x) + P(x, y)) + log(pu(y) + P(x, y))) over the pairs that can form
  // (NEG_INF where a factor is 0), from this nested-only run's outside pass.
  void nested_lphi(StemEst& st, double gamma) {
    outside();
    const size_t nn = static_cast<size_t>(n) * n;
    std::vector<double> P(nn, 0.0), sp(n, 0.0);
    for (int j = 0; j < n; ++j) {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        const int i = at.keys[r];
        if (i < 0) continue;
        const double b = Qb_b[ix(i, j)];
        if (at.vals[r] < HALF_NEG || b < HALF_NEG) continue;
        const double p = std::exp(at.vals[r] + b - logZ);
        P[ix(i, j)] = p;
        sp[i] += p;
        sp[j] += p;
      }
    }
    st.lphi.assign(nn, NEG_INF);
    for (int x = 0; x < n; ++x)
      for (int y : partners_right[x]) {
        const double p = P[ix(x, y)];
        const double a = std::max(0.0, 1.0 - sp[x]) + p, b = std::max(0.0, 1.0 - sp[y]) + p;
        if (a > 0.0 && b > 0.0) st.lphi[ix(x, y)] = gamma * (std::log(a) + std::log(b));
      }
  }

  // --stem-est rank (role both: Qgl / Qgr also scored as stem 2).
  // qg2 (role all): Qg also scored as stem 2 (needs st->mirror).  qg_partner
  // (--stem-est-qg-r partner): Qg's stem-2 score from the partner tables (stem_Bg) instead of
  // the mirrored index.  consistent (--stem-est-consts consistent): every role adds beta1, the
  // closing at(p, q) and only the stem terms not yet in alpha (StemEval).
  void set_loop_cap(int lc) { cap = lc; }
  void set_stem(const StemEst* st, bool both, bool qg2 = false, bool qg_partner = false,
                bool consistent = false) {
    stem_ = st;
    stem_both = both;
    stem_qg2 = qg2;
    stem_qg_partner = qg_partner;
    stem_consistent = consistent;
    first_f = stemidx::MinTree(n + 1);
  }
  // --stem-est-qg-ext (mode 1 stack, 2 bulge1; needs the partner Qg score), K_max, split.
  void set_stem_ext(int mode, int kmax, bool split, bool check) {
    stem_ext = mode;
    stem_ext_max = kmax;
    stem_split = split;
    ext_check = check;
  }
  // --stem-est-qg-ext loop parameters (Lg, cut c, node cap N) and --stem-est-helpers grow
  // helpers: 0 off, 1 grow, 2 nested
  void set_stem_grow(int lg, double cut, int nodes, int helpers, int ml_span) {
    grow_lg = lg;
    grow_cut = cut;
    grow_nodes = nodes;
    stem_helpers = helpers;
    stem_ml_span = ml_span;
  }
  // --stem-est-profile on, --stem-est-profile-eps, --stem-est-ml-tail nested,
  // --stem-profile-check
  void set_stem_profile(bool on, double eps, bool tail, bool check) {
    stem_profile = on;
    prof_eps = eps;
    ml_tail = tail;
    prof_check = check;
  }
  // --stem-est-chain on (needs set_stem_profile on and st->chain), -span, -cut, -ml, -loop
  void set_stem_chain(bool on, int span, double cut, bool ml, int loop) {
    stem_chain = on;
    chain_span = span;
    chain_cut = cut;
    chain_ml = ml;
    chain_loop = loop;
  }
  // --stem-est-hr av / phi (needs st->hr_av): every L-stem score uses Hav(u, v, q) S[q+1]
  // (av) or Hav(u, v, q) / C0[q] (phi) in place of H(u, v, q) / C0[q]; feasibility
  // (StemEval::ok) keeps H.  check: --stem-hr-check.
  void set_stem_hr(bool av, bool sq, bool check) {
    hr_av = av;
    hr_sq = sq;
    hr_check = check;
  }
  // --stem-est-principled on (needs st->pr): every pseudoknot table ranked by
  // alpha + log Out_hat; growth / future stems use the interior loop cap loop and multiloop
  // steps (ml) of span <= ml_span; Qgls closings run <= ml_close; profiles cut at cut x max;
  // arms_k arms per gap (0 = all); check: --stem-pr-check.
  void set_stem_principled(bool on, int loop, bool ml, int ml_span, int ml_close, double cut,
                           int arms_k, bool check) {
    pr_on = on;
    pr_loop = loop;
    pr_ml = ml;
    pr_ml_span = ml_span;
    pr_ml_close = ml_close;
    pr_cut = cut;
    pr_arms_k = arms_k;
    pr_check = check;
  }
  // --pr-arms-rank ratio: arms by phi / Z_arm (Z_arm = wu per arm base x'-xs+1), added in
  // decreasing order until the rest of the gap's total < arms_eps x the selected sum, at most
  // arms_max (replaces --pr-arms-k).  --pr-cut-by gapfree: profile entries cut by their
  // gap-free contribution (growth: x C[p'-1] Sbar[q'+1]; future stems: x Sbar[ys+1] /
  // e^(wu (x - xs))); cut_roots: the same cut on the root and key profiles.
  void set_stem_pr_select(bool ratio, double arms_eps, int arms_max, bool gapfree,
                          bool cut_roots) {
    pr_ratio = ratio;
    pr_arms_eps = arms_eps;
    pr_arms_max = arms_max;
    pr_gapfree = gapfree;
    pr_cut_roots = cut_roots;
  }
  void set_stem_pr_sep(int mode, int arm_z, bool qz_unp, bool check) {
    pr_arms_mode = mode;
    pr_arm_z = arm_z;
    pr_qz_unp = qz_unp;
    pr_sep_check = check;
  }
  // --pr-grow-select bound: growth node, root and key profiles keep their endpoints by
  // iota = w U(p', q') (largest first, until 1 - eps of the total, at most B; B = 0: no cap;
  // eps = 0 and B = 0: everything); check: --pr-grow-check; dump: --pr-score-dump.
  // --pr-armq-sep on: phi(arm) at q' with Qzbar[q'+1][y'-1] replaced by the upper bound
  // Qzbar[0][y'-1] / Qzbar[0][q'], from per-arm suffix sums built once.
  void set_stem_pr_armq_sep(bool on, size_t aq_max) {
    pr_aq_sep = on;
    pr_aq_max = aq_max;
  }
  // --pr-arm-persist on: the top-k arm rankings of every q' from one persistent argmax tree per
  // order (DP03::pidx_build).  --pr-arm-gaponly on: Step A over the arms whose x' lies in the
  // union of the gaps of the states being scored (DP03::pr_gap_begin).  --pr-partner-span W:
  // only partner entries with y' - q' <= W (0: all).  --pr-v-eps e: V maps and arm maps keep
  // their entries by V(xs, ys) Sbar[ys+1] / Qzbar[xs][x] until 1 - e of the total (default
  // 1e-8; 0: all).
  void set_stem_pr_arms(bool persist, bool gaponly, int pspan, double veps) {
    pr_arm_persist = persist;
    pr_arm_gaponly = gaponly;
    pr_pspan = pspan;
    pr_veps = veps;
  }
  // --pr-v-ratio on: arm entries from Vr; --pr-v-diag PREFIX: PREFIX.ventries.tsv.gz and, with
  // --pr-v-sample P, PREFIX.sel.tsv (each L-lookup sampled with probability P, at most 10000).
  void set_stem_pr_v(bool ratio, const std::string& diag, double sample_p) {
    pr_vratio = ratio;
    pr_vdiag = diag;
    pr_v2dump = g_v2dump;
    pr_vsample_p = diag.empty() ? 0.0 : sample_p;
    if (pr_vsample_p > 0.0) {
      pr_vsel = std::fopen((diag + ".sel.tsv").c_str(), "w");
      std::fprintf(pr_vsel, "qp\tu\tv\ttwo\tnarmsV\tnarmsVr\toverlap\tvr_in_v\tshare_vr\tshare_vtop\tshare_eng\t"
                            "overlap_r2\tr2_in_v\tshare_r2\toverlap_sep\tsep_in_v\tshare_sep\n");
      pr_varms = popen(("gzip -1 > '" + diag + ".arms.tsv.gz'").c_str(), "w");
      std::fprintf(pr_varms, "sample\txs\txp\tlogV\tlogVr\tlogVr2\tlogVsep\n");
    }
  }
  void pr_v_finish(const std::string& sid) {
    pr_vdump(sid);
    if (pr_vsel) std::fclose(pr_vsel);
    pr_vsel = nullptr;
    if (pr_varms) pclose(pr_varms);
    pr_varms = nullptr;
  }
  void set_stem_pr_grow(bool bound, double eps, int B, bool check, FILE* dump) {
    pr_gsel = bound;
    pr_geps = eps;
    pr_gB = B;
    pr_gcheck = check;
    pr_dump = dump;
  }
  // sep: distinct queries, columns of the range-sum tables built, (step, q') groups queried,
  // their finite arms, steps with a query, table time; check: queries / arm terms compared,
  // max relative error of S against brute force, max |log sep - log exact|, mismatches
  // (finite against infinite), histograms of the log ratio separable / exact per arm term
  // (weighted by its share of the query's exact sum) and per query, on [-10, 10] in 0.25 bins
  // plus two tails
  void pr_sep_report(FILE* f, const std::string& sid) const {
    std::string ht, hq;
    char b[64];
    double wt = 0.0;
    for (int k = 0; k < PSC_NB; ++k) {
      std::snprintf(b, sizeof b, "%s%.6g", k ? "," : "", psc_ht[k]);
      ht += b;
      std::snprintf(b, sizeof b, "%s%.6g", k ? "," : "", psc_hq[k]);
      hq += b;
      wt += psc_ht[k];
    }
    std::fprintf(f,
                 "[pr-sep-check] %s queries %lld terms %lld tab_rel %.3e tab_mis %lld "
                 "sep_max %.6e sep_mis %lld q_mean %.6e q_abs %.6e term_wmean %.6e "
                 "term_wabs %.6e\n[pr-sep-hist] %s term %s\n[pr-sep-hist] %s query %s\n",
                 sid.c_str(), psc_q, psc_terms, psc_tab_rel, psc_tab_mis, psc_sep_max, psc_sep_mis,
                 psc_q ? psc_d_sum / psc_q : 0.0, psc_q ? psc_d_abs / psc_q : 0.0,
                 wt ? psc_lr_w / wt : 0.0, wt ? psc_lr_wabs / wt : 0.0, sid.c_str(), ht.c_str(),
                 sid.c_str(), hq.c_str());
  }
  long long pr_sep_steps_all() const { return pr_sep_steps + (prSepTouched.empty() ? 0 : 1); }
  static constexpr int PSC_NB = 82;
  long long pr_sep_q = 0, pr_sep_cols = 0, pr_sep_groups = 0, pr_sep_garms = 0,
            pr_sep_steps = 0, psc_q = 0, psc_terms = 0, psc_tab_mis = 0, psc_sep_mis = 0;
  double pr_sep_sec = 0.0, psc_tab_rel = 0.0, psc_sep_max = 0.0, psc_lr_w = 0.0,
         psc_lr_wabs = 0.0, psc_d_sum = 0.0, psc_d_abs = 0.0;
  double psc_ht[PSC_NB] = {}, psc_hq[PSC_NB] = {};
  static int psc_bin(double x) {
    if (x < -10.0) return 0;
    if (x >= 10.0) return PSC_NB - 1;
    return 1 + std::min(PSC_NB - 3, static_cast<int>((x + 10.0) / 0.25));
  }
  // --stem-est-principled-nested on (needs st->pr_nested): Qb / Qm / Qm2 ranked by
  // alpha + their nested outside.
  void set_stem_pr_nested(bool on) { pr_nested = on; }
  // principled: per table (Qg, Qgl, Qgr, Qgls, Qgrs, Qp) ranked states and scoring time;
  // key-profile entries evaluated; T builds (cache hits), arms selected / skipped (outside
  // the gap) over them, ys entries of the T's; K builds (hits) and their entries; SL+SR
  // misses (hits); per-q' arm structures built; static part: V entries, arms, arm entries,
  // G profiles and entries, entries cut; check: states compared per table, finite / infinite
  // mismatches, max relative error
  long long pr_states[6] = {}, pr_evals = 0, pr_T_built = 0, pr_T_hits = 0, pr_T_arms = 0,
            pr_T_skip = 0, pr_T_ys = 0, pr_K_built = 0, pr_K_hits = 0, pr_K_is = 0,
            pr_slr_miss = 0, pr_slr_hits = 0, pr_aq_built = 0, pr_V_entries = 0, pr_arms = 0,
            pr_arm_entries = 0, pr_G_nodes = 0, pr_G_entries = 0, pr_cut_drop = 0,
            pr_roots = 0, pr_keys = 0, pr_key_entries = 0;
  double pr_sec[6] = {}, pr_static_sec = 0.0;
  // grow select: profiles selected (0 growth nodes, 1 evaluated roots / keys), their entries
  // before / after, summed kept iota shares, profiles with zero total; U_L columns added,
  // partner entries added, U evaluations; growth and selection time (all modes: growth).
  // check: endpoints, bound violations (L, R), max log(L / (Qz U_L)), max log(R / (Qz U_R)),
  // endpoints whose selected-arm L exceeds the all-arm L
  long long gs_prof[2] = {}, gs_before[2] = {}, gs_after[2] = {}, gs_zero[2] = {}, gs_cols = 0,
            gs_rent = 0, gs_U = 0, gc_n = 0, gc_vL = 0, gc_vR = 0, gc_sel_gt = 0;
  double gs_share[2] = {}, pr_grow_sec = 0.0, pr_sel_sec = 0.0, gc_maxL = NEG_INF,
         gc_maxR = NEG_INF;
  // g-read profile: one-time fill vs per-root sw_as_rn aggregation.
  double sw_ntr_init_sec = 0.0, sw_fill_sec = 0.0, sw_as_rn_sec = 0.0;
  long long sw_as_rn_calls = 0, sw_as_rn_ent_in = 0, sw_as_rn_ent_out = 0;
  long long sw_ga_nodes = 0, sw_ga_ent_sum = 0, sw_ga_ent_max = 0;
  // V/G fill split (pass-2 V sweep vs scan sw_fill).
  double v_merge_sec = 0.0, v_cut_sec = 0.0, g_merge_sec = 0.0, g_cut_sec = 0.0;
  long long v_pushes = 0, g_pushes = 0;
  long long v_pre_n = 0, v_pre_sum = 0, v_pre_max = 0, v_post_n = 0, v_post_sum = 0, v_post_max = 0;
  long long g_pre_n = 0, g_pre_sum = 0, g_pre_max = 0, g_post_n = 0, g_post_sum = 0, g_post_max = 0;
  // sw_fill multiloop linearity: right-gap yo visits attributed to each pair at column y.
  long long ml_yo_qm_n = 0, ml_yo_qm_sum = 0, ml_yo_qm_max = 0;
  long long ml_yo_sc_n = 0, ml_yo_sc_sum = 0, ml_yo_sc_max = 0;
  int ml_yo_S = 0;
  bool ml_yo_qm_ix = false;
  // K_L inside growth reads (lsc_L) and per-arm (b,c) V/K stats.
  double kl_in_read_sec = 0.0;
  long long arm_read_n = 0;
  long long arm_Vf_sum = 0, arm_Vf_max = 0;
  std::vector<int> arm_Vf_samples;  // for median
  long long arm_Kq_sum = 0, arm_Kq_max = 0;
  std::vector<char> arm_seen;
  // steps of the built nodes by kind: [0] growth nodes + roots, [1] future stems;
  // [.][0] interior, [.][1] multiloop.  Partner terms: wall time of the outermost T / K /
  // Stop build.  stop_j: the step at which a CPU-limit stop was taken (-1: none).
  long long st_kind[2][2] = {}, st_nodes[2] = {};
  double pr_partner_sec = 0.0;
  bool pr_aq_sep = false;
  double pr_aq_sec = 0.0;
  std::vector<double> prSF;
  // inclusive seconds of the pr_T / pr_K / pr_SLR / pr_W misses; W misses and entries read
  double pr_sec_T = 0.0, pr_sec_K = 0.0, pr_sec_SLR = 0.0, pr_sec_W = 0.0;
  long long pr_W_built = 0, pr_W_ent = 0;
  size_t pr_aq_max = 128;  // --pr-aq-max
  bool pr_arm_persist = false, pr_arm_gaponly = false;
  int pr_pspan = 0;
  double pr_veps = 0.0;
  // --pr-v-trim suffix: drop xs.ys only against the entries with xs' >= xs; --pr-v-cap: entries
  // per key kept (largest share of their suffix first), 0 = no cap
  bool pr_vsuffix = false;
  int pr_vcap = 0;
  // Pass 2 limits. All default off: 0 leaves the current tables and ranking keys.
  int p2_qz_beam = 0, p2_qz_key = 1, p2_v_mode = 0, p2_v_bp = 0, p2_v_cap = 0, p2_v_trim = 0;
  int p2_g_mode = 0, p2_g_bp = 0, p2_g_cap = 0;
  bool sw_check = false, p2_sweep = false, p2_vsweep = false, geq_pre = false;
  int p2_v_floor = 0;  // 0 off, 1 smallest kept mu, 2 log-weight 0
  bool p2_v_cut_fb = false;  // --v-cut-fb: kind-(i) miss = kappa_min / C_pre(b-1)
  bool lazy_ml_cap = false;  // --lazy-ml-cap: lazy leading run uses ml_lead_cap
  bool p2_g_suffix = false, p2_g_eager = false, p2_g_fb = false;
  // --pre-key inside: beam rank is the inside weight. Default is C(i-1) + inside.
  bool pre_key_inside = false;
  // --qz-beam on the nested pre-pass: fill Q^z_pre beamed instead of the dense table.
  // The scan's own Q^z (part of Z) is never beamed.
  bool qz_beam_fill = false, qzbar_old = false, prof_on = false;
  long long qz_cand = 0, qz_kept = 0;
  // Pass-1 rule iterations (per nested pfunc). Printed as [p1rule] per nucleotide.
  long long r_hair = 0, r_int = 0, r_m2xp = 0, r_m2multi = 0, r_muext = 0, r_mmu = 0, r_mp = 0,
            r_extc = 0;
  double sec_qb = 0, sec_qm = 0, sec_branch = 0, sec_qz = 0, sec_c = 0;
  double sec_hair = 0, sec_int = 0, sec_m2xp = 0, sec_emit = 0, sec_mclose = 0, sec_mp = 0;
  double prof_qz_sec = 0.0;
  int prof_region = 0;             // 0 other, 1 R prune, 2 spawn rank, 3 H / P3
  double prof_grow_reg[4] = {};
  double sx_k_sec = 0, sx_rkey_sec = 0;
  long long sx_k_n = 0, sx_m_sum = 0;
  // 0 other, 1 R root (rec_root), 2 R key (rec_key), 3 spawn rank
  int sx_ksrc = 0;
  long long sx_k_src[4] = {};
  long long sx_rkey_n = 0;
  long long sx_m_h[8] = {};
  // query-weighted length of the node-trunc list that the K call walks (no layer m on that path)
  long long sx_list_h[65] = {};
  // --k-memo: lK(arm, q') = arm_phi, a function of the pass-2 arm tables only
  bool k_memo = false;
  std::unordered_map<uint64_t, double> lKm;
  long long lKm_calls = 0, lKm_hits = 0;
  // --k-bound: compare lK against the exact gap sum once, after the arms exist.
  // --k-exact-ub: R beams rank by the exact gap sum, in upper-bound order, and stop early.
  bool k_bound = false, k_exact_ub = false, k_exact_check = false, sx_k_exact = false;
  std::unordered_map<uint64_t, double> lKxM;
  long long lKx_calls = 0, lKx_hits = 0, sx_kx_sigma = 0, k_exact_cols = 0, k_exact_mis = 0;
  double sx_kx_sec = 0.0;
  long long prof_help_steps = 0;
  struct QkEnt {
    int k;
    double q, p, u1, u2;
  };
  std::vector<std::vector<QkEnt>> qz_keep;
  long long v_capped = 0;
  // cap: pairs capped; per capped pair the dropped share of the trim weights (sum, max) and of
  // V itself (max); log sums of V over the capped pairs' entries and over the dropped ones
  long long v_cap_pairs = 0;
  double v_cap_fsum = 0.0, v_cap_fmax = 0.0, v_cap_fwmax = 0.0, v_cap_all = NEG_INF,
         v_cap_drop = NEG_INF;
  // --pr-v-mlk K (0 = all): a V node (scan = false, pr_fast 2) keeps its K multiloop helper
  // terms of largest mass w + log sum(helper) and scales them by (kept + dropped) / kept, so
  // the node's mass is unchanged
  int pr_vmlk = 0;
  // --pr-v-mlsteps on: the V nodes enumerate their multiloop steps under pr_fast 2 as well (only
  // the step counts use them; the helpers carry those terms)
  bool pr_vmlsteps = false;
  long long vk_nodes = 0, vk_cand = 0, vk_drop = 0;
  double vk_rmax = 0.0, vk_rsum = 0.0, vk_mass_s = 0.0;
  // --pr-v-key xs: V keyed by xs alone (key2(xs, 0)), the value the sum over ys of
  // V e^beta1 Sbar[ys+1]; arm entries one per (xs, x', y') (prYw), no per-ys entries
  // (prEy / prEys / prEw empty, prArmE0 = prArmY0).
  bool pr_vxs = false;
  // --pr-v-stream on: pass 0 and pass 1 in one y-descending sweep; V(x', y') freed once
  // y' - y - 1 exceeds every step's reach; arm entries bucketed by xs, each bucket freed as
  // it is appended to the arm tables.  v_live_max: most V entries alive at once.
  bool pr_vstream = false;
  long long v_live_max = 0, v_pairs = 0;
  // persistent index: build seconds, nodes, leaf updates; gaponly: Step A arms indexed, union
  // epochs; V trimming: entries before / after
  double pidx_sec = 0.0, pidx_mb = 0.0;
  long long pidx_nodes = 0, pidx_upd = 0, aq_arms = 0, gap_epochs = 0, gap_miss = 0, v_before = 0,
            v_after = 0;
  struct SecTimer {
    double& acc;
    std::chrono::steady_clock::time_point t0;
    explicit SecTimer(double& a) : acc(a), t0(std::chrono::steady_clock::now()) {}
    ~SecTimer() { acc += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); }
  };
  // stop_phase: where the stop was taken (0 scan, 1 the one-time V build); cur_j: current step.
  int pr_pdepth = 0, stop_j = -1, stop_phase = 0, cur_j = -1;
  void cpu_stop_check(int phase) {
    if (!g_cpu_stop) return;
    stop_j = cur_j;
    stop_phase = phase;
    throw CpuStop{};
  }
  struct PrPartnerTimer {
    DP03& d;
    std::chrono::steady_clock::time_point t0;
    explicit PrPartnerTimer(DP03& x) : d(x), t0(std::chrono::steady_clock::now()) { ++d.pr_pdepth; }
    ~PrPartnerTimer() {
      if (--d.pr_pdepth == 0)
        d.pr_partner_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }
  };
  long long prc_states[6] = {}, prc_mis[6] = {};
  double prc_rel[6] = {};
  // av / phi: ranked Qg states with both L-stem factors finite, the sum / min / max of
  // (log Hav + hr_c(q)) - (log H - log C0[q]) over them, states finite under the new mode only
  // / orig only; check: Qg / Qgl states compared with hr_naive, finite / infinite mismatches,
  // max relative error
  long long hr_n = 0, hr_av_only = 0, hr_orig_only = 0, hrc_states = 0, hrc_mis = 0;
  double hr_dsum = 0.0, hr_dmin = 0.0, hr_dmax = 0.0, hrc_rel = 0.0;
  // chain: chain nodes built, edges by step (b interior, f, g, h multiloop), root profiles,
  // their L / R entries, entries dropped by the cut; chain time (inside pair_prof)
  long long ch_nodes = 0, ch_edges[4] = {}, ch_roots = 0, ch_ly = 0, ch_rx = 0, ch_cut_drop = 0;
  double ch_sec = 0.0;
  // profiles: steps with a profile, keys built (Qg pairs, Qgls (i, seed), Qgrs c), pair
  // profiles built, L / R entries of the key profiles, entries dropped by eps, entries
  // evaluated by the ranked states, tail closings; Bg lookups, raw computations and the
  // partner entries those visit; build / total scoring time
  long long pf_steps = 0, pf_keys[3] = {}, pf_pairs = 0, pf_ly = 0, pf_rx = 0, pf_trunc = 0,
            pf_evals = 0, pf_tail_cl = 0, pf_states = 0;
  long long bg_lookups = 0, bg_misses = 0, bg_visits = 0;
  double pf_build_sec = 0.0, pf_total_sec = 0.0;
  // --stem-est-bg-reorder on: M rows built (key, v), the (x', e) entries they visit, the e
  // entries visited by the scored states; check: states compared, finite / infinite
  // mismatches, max relative error, ranking events / states whose keep decision differs
  bool bg_reorder = false, ro_check = false;
  long long ro_rows = 0, ro_build = 0, ro_eval = 0, roc_states = 0, roc_mis = 0, roc_events = 0,
            roc_flip_events = 0, roc_flip_states = 0;
  double roc_rel = 0.0, ro_off = NEG_INF;
  void set_stem_bg_reorder(bool on, bool check) {
    bg_reorder = on;
    ro_check = check;
  }
  // --stem-profile-check: Qg / helper states compared, finite / infinite mismatches, max
  // relative error against the Result-19 score
  long long pc_states[2] = {}, pc_mis[2] = {};
  double pc_rel[2] = {};
  long long stem_removed[3] = {};  // Qg, Qgl, Qgr
  double stem_qg_sec = 0.0;  // time in the Qg stem scoring (feasibility + ranking)
  // sigma_ext: ranked Qg states, extension nodes over them, states with a finite sigma_ext,
  // time; --stem-ext-check: states compared, finite / infinite mismatches, max relative error
  long long ext_states = 0, ext_nodes = 0, ext_finite = 0;
  double ext_sec = 0.0;
  long long xc_states = 0, xc_mis = 0;
  double xc_rel = 0.0;
  // growth node lists built (one per outer pair and step), their nodes, layer nodes reached
  // and dropped by the cut / the node cap; helpers: ranked Qgls / Qgrs states, closings and
  // growth nodes over them, finite scores, time; --stem-ext-check on the helper scores
  long long grow_lists = 0, grow_list_nodes = 0, grow_layer = 0, grow_cut_drop = 0,
            grow_cap_drop = 0;
  long long hlp_states[2] = {}, hlp_closings = 0, hlp_nodes = 0, hlp_finite = 0;
  double hlp_sec = 0.0;
  long long hc_states = 0, hc_mis = 0;
  double hc_rel = 0.0;

  // --stem-bg-check (after the outside pass, every table at infinite beam): stem_Bg of every
  // stored Qg state against a brute-force sum over the stored Qgl states, and the posterior
  // weight of the Qg states a d-dependent removal (H_R = 0 and Bg(d, u, v) = 0) would drop.
  std::string stem_bg_check(const std::string& id, double logZ) {
    struct L { int i, e, f; double val; };
    std::unordered_map<int, std::vector<L>> by_d;
    for (int f = 0; f < n; ++f)
      for (const Grp& g : Qgl[f].grps)
        for (const Ent& s : g.ents)
          if (E.end_ok(s.p, f))
            by_d[s.u + 1].push_back({s.p, g.v, f, prefix_log(s.p) + s.a + E.at(s.p, f)});
    stBc.clear();
    stBgc.clear();
    stc_q = -1;
    long long states = 0, zero_mis = 0, rule_n = 0, ff_only = 0, ff_viol = 0;
    double max_rel = 0.0, rule_w = 0.0;
    std::vector<long double> terms;
    for (int q = 0; q < n; ++q)
      for (const Grp& g : Qg[q].grps)
        for (const Ent& s : g.ents) {
          const int d = s.p, u = s.u, v = g.v;
          ++states;
          terms.clear();
          auto it = by_d.find(d);
          if (E.end_ok(u, v) && it != by_d.end())
            for (const L& x : it->second) {
              if (x.e < u + 1 || x.e > v - 1 || x.f > v - 1) continue;
              const double t = x.val + qz_log(x.f + 1, v - 1) + qz_log(u + 1, x.e - 1);
              if (t > HALF_NEG) terms.push_back(t);
            }
          double brute = NEG_INF;
          if (!terms.empty()) {
            long double m = *std::max_element(terms.begin(), terms.end()), sum = 0.0L;
            for (long double t : terms) sum += std::exp(t - m);
            brute = static_cast<double>(m + std::log(sum) +
                                        (static_cast<long double>(E.w_stem) + E.at(u, v)));
          }
          const double bg = stem_Bg_raw(d, u, v);
          if ((bg > HALF_NEG) != (brute > HALF_NEG)) ++zero_mis;
          else if (bg > HALF_NEG) max_rel = std::max(max_rel, std::fabs(std::expm1(bg - brute)));
          const bool h0 = stem_->H(u, v, q) < HALF_NEG;
          const bool ff = first_f.query(u + 1, v) <= v - 1;
          if (bg > HALF_NEG && !ff) ++ff_viol;
          if (h0 && bg < HALF_NEG) {
            ++rule_n;
            if (s.b > HALF_NEG) rule_w += std::exp(s.a + s.b - logZ);
            if (ff) ++ff_only;
          }
        }
    char buf[512];
    std::snprintf(buf, sizeof buf,
                  "[bg-check] %s n %d states %lld max_rel %.3e zero_mismatch %lld ff_violations "
                  "%lld rule_drop %lld rule_drop_ff_ok %lld rule_drop_w %.6e",
                  id.c_str(), n, states, max_rel, zero_mis, ff_viol, rule_n, ff_only, rule_w);
    return buf;
  }

  // --stem-diag: one JSON line of stem-estimate diagnostics (definition at stem_diag).
  std::string stem_diag(const std::vector<const StemEst*>& ests, const std::string& id);
  void tagged_outside(int tag);

