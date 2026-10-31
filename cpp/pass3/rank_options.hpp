// Pass-3 ranking options: setters for --stem-est* / --pr-* (stem score, profile, chain,
// principled, arms, V options), their counters and timers, and cpu_stop_check.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.


  // --stem-est rank (role both: Qgl / Qgr also scored as stem 2).
  // qg2 (role all): Qg also scored as stem 2 (needs st->mirror).  qg_partner
  // (--stem-est-qg-r partner): Qg's stem-2 score from the partner tables (stem_Bg) instead of
  // the mirrored index.
  void set_loop_cap(int lc) { cap = lc; }
  void set_stem(const StemEst* st, bool both, bool qg2, bool qg_partner) {
    stem_ = st;
    stem_both = both;
    stem_qg2 = qg2;
    stem_qg_partner = qg_partner;
    first_f = stemidx::MinTree(n + 1);
  }
  // --stem-est-qg-ext (mode 1 stack, 2 bulge1, 3 loop; needs the partner Qg score).
  void set_stem_ext(int mode) { stem_ext = mode; }
  // --stem-est-grow-loop Lg and --stem-est-helpers (0 off, 1 grow, 2 nested).
  void set_stem_grow(int lg, int helpers) {
    grow_lg = lg;
    stem_helpers = helpers;
  }
  void set_stem_profile(bool on) { stem_profile = on; }
  // Chain scores (needs set_stem_profile on and st->chain).
  void set_stem_chain(bool on) { stem_chain = on; }
  // --stem-est-principled on (needs st->pr): every pseudoknot table ranked by
  // alpha + log Out_hat; multiloop steps of span <= ml_span; profiles cut at cut x max.
  void set_stem_principled(bool on, int ml_span, double cut) {
    pr_on = on;
    pr_ml_span = ml_span;
    pr_cut = cut;
  }
  // gapfree: profile entries cut by their gap-free contribution (growth: x C[p'-1] Sbar[q'+1];
  // future stems: x Sbar[ys+1] / e^(wu (x - xs))).
  void set_stem_pr_select(bool gapfree) { pr_gapfree = gapfree; }
  // --pr-grow-select bound: growth node, root and key profiles keep their endpoints by
  // iota = w U(p', q') (largest first, until 1 - eps of the total).
  // --pr-armq-sep on: phi(arm) at q' with Qzbar[q'+1][y'-1] replaced by the upper bound
  // Qzbar[0][y'-1] / Qzbar[0][q'], from per-arm suffix sums built once.
  void set_stem_pr_armq_sep(bool on) { pr_aq_sep = on; }
  // --pr-v-eps e: V maps and arm maps keep their entries by V(xs, ys) Sbar[ys+1] / Qzbar[xs][x]
  // until 1 - e of the total (default 1e-8; 0: all).
  void set_stem_pr_arms(double veps) { pr_veps = veps; }
  void set_stem_pr_grow(bool bound, double eps) {
    pr_gsel = bound;
    pr_geps = eps;
  }
  long long pr_sep_steps = 0;
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
  long long gs_prof[2] = {}, gs_before[2] = {}, gs_after[2] = {}, gs_zero[2] = {}, gs_cols = 0, gs_rent = 0, gs_U = 0;
  double gs_share[2] = {}, pr_grow_sec = 0.0, pr_sel_sec = 0.0;
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
  size_t pr_aq_max = 128;
  bool pr_arm_persist = false, pr_arm_gaponly = false;
  int pr_pspan = 0;
  double pr_veps = 0.0;
  // --pr-v-trim suffix: drop xs.ys only against the entries with xs' >= xs; --bV: entries
  // per key kept (largest share of their suffix first), 0 = no cap
  bool pr_vsuffix = false;
  int pr_vcap = 0;
  // Pass 2 limits. All default off: 0 leaves the current tables and ranking keys.
  int p2_qz_beam = 0, p2_qz_key = 1, p2_v_cap = 0;
  bool p2_sweep = false, p2_vsweep = false;
  // Nested pre-pass: fill Q^z_pre beamed instead of the dense table. The scan's own Q^z
  // (part of Z) is never beamed.
  bool qz_beam_fill = false;
  long long qz_cand = 0, qz_kept = 0;
  // Pass-1 rule iterations (per nested pfunc). Printed as [p1rule] per nucleotide.
  long long r_hair = 0, r_int = 0, r_m2xp = 0, r_m2multi = 0, r_muext = 0, r_mmu = 0, r_mp = 0,
            r_extc = 0;
  double sec_qb = 0, sec_qm = 0, sec_branch = 0, sec_qz = 0, sec_c = 0;
  // Pass 3 outside + the pair table built from it (bpp), always on (two clock reads per run).
  double sec_outside = 0, sec_bpp = 0;
  double sec_hair = 0, sec_int = 0, sec_m2xp = 0, sec_emit = 0, sec_mclose = 0, sec_mp = 0;
  double prof_qz_sec = 0.0;
  double sx_k_sec = 0, sx_rkey_sec = 0;
  long long sx_k_n = 0, sx_m_sum = 0;
  // 0 other, 1 R root (rec_root), 2 R key (rec_key), 3 spawn rank
  int sx_ksrc = 0;
  long long sx_k_src[4] = {};
  long long sx_rkey_n = 0;
  long long sx_m_h[8] = {};
  // query-weighted length of the node-trunc list that the K call walks (no layer m on that path)
  long long sx_list_h[65] = {};
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
  // The section clocks are a compile-time option (-DLK_FINE_TIMERS): two clock reads per
  // section cost more than the sections themselves on the per-column paths, so without the
  // macro every *_sec below stays 0 and no clock is read.
  struct SecTimer {
#ifdef LK_FINE_TIMERS
    double& acc;
    std::chrono::steady_clock::time_point t0;
    explicit SecTimer(double& a) : acc(a), t0(std::chrono::steady_clock::now()) {}
    ~SecTimer() { acc += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); }
#else
    explicit SecTimer(double&) {}
#endif
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
#ifdef LK_FINE_TIMERS
    std::chrono::steady_clock::time_point t0;
    explicit PrPartnerTimer(DP03& x) : d(x), t0(std::chrono::steady_clock::now()) { ++d.pr_pdepth; }
    ~PrPartnerTimer() {
      if (--d.pr_pdepth == 0)
        d.pr_partner_sec +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }
#else
    explicit PrPartnerTimer(DP03& x) : d(x) { ++d.pr_pdepth; }
    ~PrPartnerTimer() { --d.pr_pdepth; }
#endif
  };
  // av / phi: ranked Qg states with both L-stem factors finite, the sum / min / max of
  // (log Hav + hr_c(q)) - (log H - log C0[q]) over them, states finite under the new mode only
  // / orig only; check: Qg / Qgl states compared with hr_naive, finite / infinite mismatches,
  // max relative error
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
  long long stem_removed[3] = {};  // Qg, Qgl, Qgr
  double stem_qg_sec = 0.0;  // time in the Qg stem scoring (feasibility + ranking)
  // sigma_ext: ranked Qg states, extension nodes over them, states with a finite sigma_ext,
  // time
  long long ext_states = 0, ext_nodes = 0, ext_finite = 0;
  double ext_sec = 0.0;
  // growth node lists built (one per outer pair and step), their nodes, layer nodes reached
  // and dropped by the cut / the node cap; helpers: ranked Qgls / Qgrs states, closings and
  // growth nodes over them, finite scores, time
  long long grow_lists = 0, grow_list_nodes = 0, grow_layer = 0, grow_cut_drop = 0,
            grow_cap_drop = 0;
  long long hlp_states[2] = {}, hlp_closings = 0, hlp_nodes = 0, hlp_finite = 0;
  double hlp_sec = 0.0;

