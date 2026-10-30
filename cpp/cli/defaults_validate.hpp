// main(), part 4: derived defaults (--pk-l2r on, b1 / b2 / pbeam / b_V, sweep on), flag
// consistency checks, output files, beam audit.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  // --pk-l2r on defaults: --l2r-spawn inc; pbeam derived below as (b2,b4,b4) unless
  // --l2r-pbeam / --pbeam-ratio given; --l2r-p3-hex/lin on; --pr-ml-span 100000; with
  // --grow-rec on (ml all/all) --bG 16 --grow-node-trunc-r -1; --pr-v-key xs when safe;
  // then --grow-r-from-v on if its requirements hold.
  if (pk_l2r == 1) {
    auto given = [&](const char* f) {
      for (int k = 1; k < argc; ++k)
        if (std::strcmp(argv[k], f) == 0) return true;
      return false;
    };
    if (!given("--l2r-spawn")) {
      l2r_spawn_inc = true;
      l2r_spawn_cube = l2r_spawn_static = false;
    }
    // pbeam filled after b2/b4 are known (see below).
    if (!given("--l2r-p3-hex")) l2r_p3_hex = true;
    if (!given("--l2r-p3-lin")) l2r_p3_lin = true;
    if (!given("--pr-ml-span")) pr_ml_span = 100000;
    const bool grow_ok = g_grow_rec && g_ml_L == 2 && g_ml_R == 2 && !g_ml_diag && !g_rec_check;
    if (grow_ok && !given("--bG")) g_ntr_K = 16;
    if (p2_g_cap > 0) g_ntr_K = p2_g_cap;
    if (grow_ok && !given("--grow-node-trunc-r")) g_ntr_KR = -1;
    if (!given("--pr-v-key") && pr_fast > 0 && !l2r_ctx && !l2r_ctx_p1 && l2r_p1 == 0 && pr_pmax == 0 &&
        !pr_v_ratio)
      pr_vxs = true;
    if (!given("--grow-r-from-v")) {
      if (g_ntr_K > 0 && g_ntr_KR < 0 && pr_vxs && pr_fast >= 2 && !l2r_ctx_p1 && !l2r_diag)
        g_rfv = true;
      else
        std::fprintf(stderr, "[defaults] --grow-r-from-v left off: needs --bG K > 0, "
                             "--grow-node-trunc-r -1, --pr-v-key xs, --pr-fast ml\n");
    }
  }
  if (pbeam_ratio_at > pbeam_explicit_at) {
    const int b4 = beam > 0 ? beam : 0;
    for (int k = 0; k < 3; ++k) l2r_pbeam[k] = pbr[k] * b4;
  }
  // b1 and b2 default to 200, independent of b4. rlt-beam defaults to b4.
  // Sweep and the V sweep are on. No length gate. Q^z and Q^z_pre use b2.
  if (!b1_set) nest_beam = 200;
  if (b2_at < 0) {
    const int b2tabs[] = {BT_QB, BT_QM, BT_QM2, BT_QMP, BT_QP};
    for (int t : b2tabs)
      if (table_beam[t] < 0) table_beam[t] = 200;
  }
  b2_width = table_beam[BT_QM] >= 0 ? table_beam[BT_QM] : 200;
  // pbeam = (b2, b4, b4) unless --l2r-pbeam or --pbeam-ratio was set.
  if (!user_pbeam && !user_pratio) {
    l2r_pbeam[0] = b2_width;
    l2r_pbeam[1] = beam;
    l2r_pbeam[2] = beam;
  }
  // --bV default 256 (replaces --v-cap / --pr-v-cap; 0 = uncapped).
  if (!bV_set) {
    pr_vcap = 256;
    p2_v_cap = 256;
  }
  std::fprintf(stderr, "[arity] b2 %d b4 %d pbeam %d,%d,%d\n", b2_width, beam, l2r_pbeam[0],
               l2r_pbeam[1], l2r_pbeam[2]);
  if (!rlt_beam_set) rlt_beam = beam;
  // Default sweep ON: sw_as_rn with q'-aware intermediate cap (= lazy trunc key).
  // Pass `--sweep off` for rec_node/pr_steps (hard14/OLD path).
  if (!sweep_set) p2_sweep = true;
  if (!vsweep_set) p2_vsweep = true;

  if (arm_score && access_file.empty()) {
    std::fprintf(stderr, "--arm-score on needs --access-file\n");
    return 2;
  }
  if (arm_score && !score_qz) {
    std::fprintf(stderr, "--arm-score on needs --score qz\n");
    return 2;
  }
  const bool stem_diag = !stem_diag_path.empty();
  if ((stem_rank || stem_diag) && (!allow_pk || arm_score || gcfg.on)) {
    std::fprintf(stderr, "--stem-est rank and --stem-diag need pseudoknots, --arm-score off and "
                         "--gate off\n");
    return 2;
  }
  const bool dump_scores = !scores_path.empty();
  const bool dump_ext = !ext_path.empty();
  bool all_inf = beam <= 0;
  for (int t = 0; t < BT_N; ++t) all_inf &= table_beam[t] <= 0;
  if (dump_ext && (!all_inf || !gap_step || stem_rank || stem_diag || dump_scores || !allow_pk ||
                   arm_score || gcfg.on || exact_path.empty())) {
    std::fprintf(stderr, "--dump-ext runs at --beam inf (every table) in step mode with "
                         "pseudoknots and --dump-exact, --stem-est off, --arm-score off, "
                         "--gate off, no --stem-diag / --dump-scores\n");
    return 2;
  }
  if (dump_scores && (!all_inf || !gap_step || stem_rank || stem_diag || !allow_pk ||
                      arm_score || gcfg.on)) {
    std::fprintf(stderr, "--dump-scores runs at --beam inf (every table) in step mode with "
                         "pseudoknots, --stem-est off, --arm-score off, --gate off, no "
                         "--stem-diag\n");
    return 2;
  }
  if (stem_diag && (stem_rank || beam > 0 || !states_path.empty() || !gate_diag_path.empty() ||
                    !do_bpp)) {
    std::fprintf(stderr, "--stem-diag runs at --beam inf with --stem-est off, the outside pass "
                         "on, and without --dump-states / --gate-diag\n");
    return 2;
  }
  for (int t = 0; t < BT_N; ++t)
    if (stem_diag && table_beam[t] > 0) {
      std::fprintf(stderr, "--stem-diag needs every table at infinite beam\n");
      return 2;
    }
  if (stem_loop < 0) {
    std::fprintf(stderr, "--stem-est-loop must be >= 0\n");
    return 2;
  }
  if (!pr_oracle_path.empty() && (!stem_pr || !all_inf || !gap_step || !allow_pk)) {
    std::fprintf(stderr, "--pr-oracle needs --stem-est-principled on, --beam inf (every "
                         "table), step mode and pseudoknots\n");
    return 2;
  }
  if ((stem_qg_r_set || stem_consts_set) && !stem_rank) {
    std::fprintf(stderr, "--stem-est-qg-r and --stem-est-consts need --stem-est rank\n");
    return 2;
  }
  if (stem_qg_partner && !stem_qg2) {
    std::fprintf(stderr, "--stem-est-qg-r partner needs --stem-est-role all\n");
    return 2;
  }
  if ((stem_ext_set && !stem_rank) || (stem_ext && !stem_qg_partner)) {
    std::fprintf(stderr, "--stem-est-qg-ext needs --stem-est rank, --stem-est-role all and "
                         "--stem-est-qg-r partner\n");
    return 2;
  }
  if (stem_ext_max < 1) {
    std::fprintf(stderr, "--stem-est-ext-max must be >= 1\n");
    return 2;
  }
  if ((stem_split || stem_ext_check) && !stem_ext) {
    std::fprintf(stderr, "--stem-est-qg-combine split and --stem-ext-check need --stem-est-qg-ext\n");
    return 2;
  }
  if (grow_set && stem_ext != 3) {
    std::fprintf(stderr, "--stem-est-grow-loop / -cut / -nodes need --stem-est-qg-ext loop\n");
    return 2;
  }
  if (grow_lg < 0 || grow_cut < 0.0 || grow_cut > 1.0 || grow_nodes < 0 || ml_span < 0) {
    std::fprintf(stderr, "--stem-est-grow-loop, -nodes, --stem-est-ml-span must be >= 0 and "
                         "--stem-est-grow-cut in [0, 1]\n");
    return 2;
  }
  if (helpers_set && !stem_ext) {
    std::fprintf(stderr, "--stem-est-helpers and --stem-est-ml-span need --stem-est-qg-ext\n");
    return 2;
  }
  if (stem_profile && (stem_ext != 3 || stem_helpers != 1 || stem_split || stem_ext_check)) {
    std::fprintf(stderr, "--stem-est-profile on needs --stem-est-qg-ext loop, --stem-est-helpers "
                         "grow, --stem-est-qg-combine sum and no --stem-ext-check\n");
    return 2;
  }
  if ((profile_set || prof_check) && !stem_profile) {
    std::fprintf(stderr, "--stem-est-profile-eps, --stem-est-ml-tail and --stem-profile-check "
                         "need --stem-est-profile on\n");
    return 2;
  }
  if (!(prof_eps >= 0.0 && prof_eps < 1.0)) {
    std::fprintf(stderr, "--stem-est-profile-eps must be in [0, 1)\n");
    return 2;
  }
  if ((bg_reorder && !stem_profile) || (ro_check && !bg_reorder)) {
    std::fprintf(stderr, "--stem-est-bg-reorder on needs --stem-est-profile on; "
                         "--stem-bg-reorder-check needs --stem-est-bg-reorder on\n");
    return 2;
  }
  if ((chain_set && !stem_chain) || (stem_chain && (!stem_profile || ml_tail))) {
    std::fprintf(stderr, "--stem-est-chain-* need --stem-est-chain on, which needs "
                         "--stem-est-profile on and --stem-est-ml-tail off\n");
    return 2;
  }
  if (!(chain_floor >= 0.0 && chain_floor < 1.0) || !(chain_cut >= 0.0 && chain_cut < 1.0) ||
      chain_span < 0 || chain_loop < 0) {
    std::fprintf(stderr, "--stem-est-chain-floor / -cut must be in [0, 1), -span / -loop >= 0\n");
    return 2;
  }
  if ((stem_hr_set || stem_hr_check) && !stem_rank) {
    std::fprintf(stderr, "--stem-est-hr, --stem-est-hr-gamma and --stem-hr-check need "
                         "--stem-est rank\n");
    return 2;
  }
  if (g_grow_rec && (!stem_pr || !(pr_abl & DP03::PA_SEP))) {
    std::fprintf(stderr, "--grow-rec on needs --stem-est-principled on and --pr-abl sep or sepv\n");
    return 2;
  }
  if (g_grow_rec && (g_ml_c < 0 || pr_check || pr_gcheck)) {
    std::fprintf(stderr, "--grow-rec on: --grow-ml-c >= 0, no --stem-pr-check / --pr-grow-check\n");
    return 2;
  }
  if (g_rec_check && (!g_grow_rec || g_ml_L != 2 || g_ml_R != 2 || g_ml_c != 0 ||
                      (pr_gsel && (pr_geps > 0.0 || pr_gB > 0)))) {
    std::fprintf(stderr, "--grow-rec-check needs --grow-rec on, ml all/all, no --grow-ml-c and "
                         "unselected profiles (--pr-grow-eps 0 --pr-grow-B 0)\n");
    return 2;
  }
  if (g_tr_K < 0 || !(g_tr_eps >= 0.0 && g_tr_eps < 1.0) ||
      (g_tr_K > 0 && (!g_grow_rec || !pr_gsel || g_ml_diag || g_rec_check || pr_abl_pk > 0))) {
    std::fprintf(stderr, "--grow-trunc K >= 0 needs --grow-rec on and --pr-grow-select bound, "
                         "no --grow-ml-diag / --grow-rec-check / --pr-abl-pk; eps in [0, 1)\n");
    return 2;
  }
  if (g_ntr_K < 0 || (g_ntr_K > 0 && (!g_grow_rec || pk_l2r != 1 || g_ml_L != 2 || g_ml_R != 2 ||
                                      g_ml_diag || g_rec_check))) {
    std::fprintf(stderr, "--bG K >= 0 needs --grow-rec on, --pk-l2r on, ml all/all, "
                         "no --grow-ml-diag / --grow-rec-check\n");
    return 2;
  }
  if (pr_vxs && (pk_l2r != 1 || l2r_ctx || l2r_ctx_p1 || l2r_p1 > 0 || pr_fast <= 0 || pr_pmax > 0 ||
                 pr_v_ratio)) {
    std::fprintf(stderr, "--pr-v-key xs needs --pk-l2r on and --pr-fast, no --l2r-ctx / --l2r-p1 / "
                         "--pr-partner-maxpairs / --pr-v-ratio (they read per-ys arm entries)\n");
    return 2;
  }
  if (g_rfv && (g_ntr_K <= 0 || g_ntr_KR >= 0 || !pr_vxs || pr_fast < 2 || l2r_ctx_p1 || l2r_diag)) {
    std::fprintf(stderr, "--grow-r-from-v on needs --bG K > 0, --grow-node-trunc-r -1, "
                         "--pr-v-key xs, --pr-fast ml; no --l2r-ctx p1 / --l2r-diag\n");
    return 2;
  }
  if (pr_vstream && (pr_fast <= 0 || pr_pmax > 0 || pr_v_ratio || pr_pmlcheck)) {
    std::fprintf(stderr, "--pr-v-stream on needs --pr-fast, no --pr-partner-maxpairs / --pr-v-ratio / "
                         "--pr-partner-ml-check\n");
    return 2;
  }
  g_skip_bml = (g_grow_rec || g_U_nobml) && (pr_abl & DP03::PA_BML) && !pr_check &&
               !(pk_l2r == 1 && l2r_ctx);
  if (stem_pr && !g_skip_bml) {
    std::fprintf(stderr, "[unsupported] principled score with Bml is off\n");
    g_skip_bml = true;
  }
  PK.beta1 += pk_db[0];
  PK.beta1m += pk_db[1];
  PK.beta1p += pk_db[2];
  if (!std::isfinite(stem_hr_gamma) || stem_hr_gamma < 0.0 ||
      (!stem_hr_phi && !(stem_hr_gamma > 0.0))) {
    std::fprintf(stderr, "--stem-est-hr-gamma must be finite and > 0 (>= 0 with phi)\n");
    return 2;
  }
  if (stem_hr_check && !stem_hr_av && !stem_hr_phi) {
    std::fprintf(stderr, "--stem-hr-check needs --stem-est-hr av or phi\n");
    return 2;
  }
  if ((stem_pr || stem_pr_nested) && (!stem_rank || !stem_qg2 || !stem_both)) {
    std::fprintf(stderr, "--stem-est-principled(-nested) on needs --stem-est rank and "
                         "--stem-est-role all\n");
    return 2;
  }
  if (stem_pr && (stem_split || stem_ext_check || stem_hr_av || stem_hr_phi || stem_hr_check ||
                  ro_check || prof_check)) {
    std::fprintf(stderr, "--stem-est-principled on excludes --stem-est-qg-combine split, "
                         "--stem-est-hr av / phi and the other stem checks\n");
    return 2;
  }
  if ((pr_set || pr_check) && !stem_pr && !stem_pr_nested) {
    std::fprintf(stderr, "--pr-* and --stem-pr-check need --stem-est-principled on\n");
    return 2;
  }
  if (pr_arms_mode == 1) {
    if (!pr_cut_by_set) pr_gapfree = true;
    if (pr_partner_floor < 0.0) pr_partner_floor = 1e-20;
  }
  if ((pr_arm_z == 1 || pr_sep_check) && pr_arms_mode != 1) {
    std::fprintf(stderr, "--pr-arm-z and --pr-sep-check need --pr-arms sep\n");
    return 2;
  }
  if (!(pr_geps >= 0.0 && pr_geps < 1.0) || pr_gB < 0 || (pr_gsel && pr_cut_roots)) {
    std::fprintf(stderr, "--pr-grow-eps must be in [0, 1), --pr-grow-B >= 0 (0 = no cap); "
                         "--pr-grow-select bound excludes --pr-cut-roots on\n");
    return 2;
  }
  if ((pr_arm_persist || pr_arm_gaponly) &&
      (pr_arms_k <= 0 || pr_ratio || pr_arms_mode == 1 || pr_gcheck)) {
    std::fprintf(stderr, "--pr-arm-persist / --pr-arm-gaponly need --pr-arms-k > 0, --pr-arms-rank "
                         "phi, no --pr-arms sep and no --pr-grow-check\n");
    return 2;
  }
  if ((pr_arm_persist && !pr_aq_sep) || pr_pspan < 0 || !(pr_veps >= 0.0 && pr_veps < 1.0)) {
    std::fprintf(stderr, "--pr-arm-persist needs --pr-armq-sep on; --pr-partner-span >= 0; "
                         "--pr-v-eps in [0, 1)\n");
    return 2;
  }
  if (pr_gsel) pr_floor = 0.0;
  FILE* pr_dump_f = nullptr;
  if (!pr_dump_path.empty() && !(pr_dump_f = std::fopen(pr_dump_path.c_str(), "w"))) {
    std::fprintf(stderr, "cannot write %s\n", pr_dump_path.c_str());
    return 2;
  }
  if (!(pr_arms_eps >= 0.0) || pr_arms_max < 1 || (pr_cut_roots && !pr_gapfree)) {
    std::fprintf(stderr, "--pr-arms-eps must be >= 0, --pr-arms-max >= 1; --pr-cut-roots on "
                         "needs --pr-cut-by gapfree\n");
    return 2;
  }
  if ((p2_qz_beam > 0 || p2_v_mode > 0 || p2_g_mode > 0 || p2_v_cap > 0 || p2_g_suffix ||
       p2_g_eager || p2_g_fb || p2_g_cap > 0) &&
      !stem_pr) {
    std::fprintf(stderr, "--qz-beam / --v-pairs / --v-cap / --g-pairs / --g-key need "
                         "--stem-est-principled on\n");
    return 2;
  }
  if (!(pr_floor >= 0.0 && pr_floor < 1.0) || !(pr_cut >= 0.0 && pr_cut < 1.0) || pr_span < 1 ||
      pr_ml_span < 0 || pr_ml_close < 0 || pr_arms_k < 0) {
    std::fprintf(stderr, "--pr-floor / --pr-cut must be in [0, 1), --pr-span >= 1, "
                         "--pr-ml-span / --pr-ml-close / --pr-arms-k >= 0\n");
    return 2;
  }
  if (stem_bg_check && (!stem_rank || !all_inf || !do_bpp)) {
    std::fprintf(stderr, "--stem-bg-check needs --stem-est rank, --beam inf (every table) and "
                         "the outside pass\n");
    return 2;
  }
  FILE* stem_diag_f = nullptr;
  if (stem_diag && !(stem_diag_f = std::fopen(stem_diag_path.c_str(), "w"))) {
    std::fprintf(stderr, "cannot write %s\n", stem_diag_path.c_str());
    return 2;
  }
  const bool gate_diag = !gate_diag_path.empty();
  const bool need_access = arm_score || gcfg.on || gate_diag;
  if ((gcfg.on || gate_diag) && access_file.empty()) {
    std::fprintf(stderr, "--gate on and --gate-diag need --access-file\n");
    return 2;
  }
  if ((gcfg.on || gate_diag) && !allow_pk) {
    std::fprintf(stderr, "--gate on and --gate-diag need pseudoknots (no --no-pk)\n");
    return 2;
  }
  AccessIndex access_idx;
  if (need_access) access_idx = index_access_file(access_file);
  FILE* diag_f = nullptr;
  if (gate_diag) {
    diag_f = std::fopen(gate_diag_path.c_str(), "w");
    if (!diag_f) {
      std::fprintf(stderr, "cannot write %s\n", gate_diag_path.c_str());
      return 2;
    }
  }
  if ((dangles_mode == 1 || dangles_mode == 2) && !mat_rna1999) {
    std::fprintf(stderr, "--dangles some|all needs --material rna1999\n");
    return 2;
  }
  if (dangles_mode == 3 && !mat_v04) {
    std::fprintf(stderr, "--dangles d2 needs --material vienna2004\n");
    return 2;
  }
  if (mat_rna1999 || mat_v04) {
    if (param_file.empty()) {
      char buf[4096];
      ssize_t len = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
      std::string exe = len > 0 ? std::string(buf, len) : std::string(argv[0]);
      std::string dir = exe.substr(0, exe.find_last_of('/') + 1);
      param_file = dir + (mat_v04 ? "../dp03/params/vienna2004.par" : "../dp03/params/rna1999.dG");
    }
    if (mat_v04) load_v04(param_file);
    else load_rna1999(param_file);
  }
  if (mat_v04 && !loop_cap_set) loop_cap = V04_MAX_LOOP;
  // Resolve per-pass caps (temporary --loop-cap-p1/p2/scan). Default: --loop-cap.
  if (loop_cap_p1 < 0) loop_cap_p1 = loop_cap;
  if (loop_cap_p2 < 0) loop_cap_p2 = loop_cap;
  if (loop_cap_scan < 0) loop_cap_scan = loop_cap;
  auto apply_max_loop = [&](int lc) {
    // Keep energy max_loop in lockstep with the DP window for this pass.
    max_loop = (lc >= NO_LOOP_CAP) ? NO_LOOP_CAP : lc;
  };
  if (loop_cap_p1 != loop_cap || loop_cap_p2 != loop_cap || loop_cap_scan != loop_cap)
    std::fprintf(stderr, "[loop-cap-pass] p1 %d p2 %d scan %d (base %d)\n", loop_cap_p1,
                 loop_cap_p2, loop_cap_scan, loop_cap);
  if (!pruned_path.empty() && gcfg.on) {
    std::fprintf(stderr, "--dump-pruned does not log the gated pruning (needs --gate off)\n");
    return 2;
  }
  if (!exact_path.empty() && pruned_path.empty() && !dump_ext) {
    std::fprintf(stderr, "--dump-exact needs --dump-pruned or --dump-ext\n");
    return 2;
  }
  FILE* states_f = nullptr;
  if (!states_path.empty() && !(states_f = std::fopen(states_path.c_str(), "w"))) {
    std::fprintf(stderr, "cannot write %s\n", states_path.c_str());
    return 2;
  }
  FILE* oracle_f = nullptr;
  if (!pr_oracle_path.empty() && !(oracle_f = std::fopen(pr_oracle_path.c_str(), "w"))) {
    std::fprintf(stderr, "cannot write %s\n", pr_oracle_path.c_str());
    return 2;
  }
  FILE* pruned_f = nullptr;
  if (!pruned_path.empty() && !(pruned_f = std::fopen(pruned_path.c_str(), "w"))) {
    std::fprintf(stderr, "cannot write %s\n", pruned_path.c_str());
    return 2;
  }
  if (!kept_path.empty() && pruned_path.empty()) {
    std::fprintf(stderr, "--dump-kept needs --dump-pruned\n");
    return 2;
  }
  FILE* kept_f = nullptr;
  if (!kept_path.empty() && !(kept_f = std::fopen(kept_path.c_str(), "w"))) {
    std::fprintf(stderr, "cannot write %s\n", kept_path.c_str());
    return 2;
  }
  FILE* scores_f = nullptr;
  if (!beam_stats_path.empty()) {
    g_beam_stats.open(beam_stats_path);
    if (!g_beam_stats.f) {
      std::fprintf(stderr, "--beam-stats: cannot open %s\n", beam_stats_path.c_str());
      return 2;
    }
  }
  if (dump_scores && !(scores_f = std::fopen(scores_path.c_str(), "w"))) {
    std::fprintf(stderr, "cannot write %s\n", scores_path.c_str());
    return 2;
  }
  FILE* ext_f = nullptr;
  if (dump_ext && !(ext_f = std::fopen(ext_path.c_str(), "w"))) {
    std::fprintf(stderr, "cannot write %s\n", ext_path.c_str());
    return 2;
  }
  std::unordered_map<std::string, std::unordered_map<uint64_t, double>> exact_states;
  if (!exact_path.empty()) exact_states = load_states(exact_path);
  const std::unordered_map<uint64_t, double> no_states;

  // Every allocated table, printed before any DP. A table with 2 or more indices and no
  // positive beam aborts the run.
  auto beam_audit = [&](int n) -> int {
    int bad = 0;
    auto row = [&](const char* stage, const char* name, int idx, const char* src, int w) {
      const int active = w > 0;
      if (idx == 2) {
        std::fprintf(stderr, "[beams] %s %s idx %d beam %s %d active %d binds %d\n", stage, name,
                     idx, src, w, active, active && w < n);
      } else {
        std::fprintf(stderr, "[beams] %s %s idx %d beam %s %d active %d\n", stage, name, idx, src,
                     w, active);
      }
      // Width 0 keeps every state (beam off). A negative width is an unset table.
      if (idx >= 2 && w < 0) {
        std::fprintf(stderr, "[beams] abort: %s %s has %d indices and no active beam\n", stage,
                     name, idx);
        bad = 1;
      }
    };
    auto tw = [&](int t) { return table_beam[t] >= 0 ? table_beam[t] : beam; };
    auto nest_tables = [&](const char* stage, int b1) {
      row(stage, "Qb", 2, "b1", b1);
      row(stage, "Qm", 2, "b1", b1);
      row(stage, "Qm2", 2, "b1", b1);
      row(stage, "Qp", 2, "b1", b1);
      row(stage, "Qmp", 2, "b1", b1);
      row(stage, "Qz", 2, "b2", b2_width);
      row(stage, "C", 1, "none", 0);
    };
    const int bV = p2_v_cap > 0 ? p2_v_cap : pr_vcap;
    const char* bVsrc = "bV";
    std::fprintf(stderr, "[beams] n %d b1 %d b2 %d b4 %d bV %d bG %d pbeam %d,%d,%d\n", n,
                 nest_beam, b2_width, beam, bV, g_ntr_K, l2r_pbeam[0], l2r_pbeam[1],
                 l2r_pbeam[2]);
    if (stem_rank) nest_tables("stem", nest_beam);
    if (stem_pr) {
      nest_tables("pass1", pre_beam > 0 ? pre_beam : nest_beam);
      row("pass2", "V", 3, bVsrc, bV);
      if (g_ntr_K > 0)
        row("pass2", "G.L", 3, "b_G", g_ntr_K);
      else
        std::fprintf(stderr, "[beams] pass2 G.L idx 3 beam b_G 0 active 0\n");
      if (g_rfv)
        row("pass2", "G.R", 3, bVsrc, bV);
      else if (g_ntr_KR > 0)
        row("pass2", "G.R", 3, "b_G", g_ntr_KR);
      else if (g_ntr_K > 0)
        row("pass2", "G.R", 3, "b_G", g_ntr_K);
      else
        std::fprintf(stderr, "[beams] pass2 G.R idx 3 beam none 0 active 0\n");
      row("pass2", "Qz_pre", 2, "b2", b2_width);
    }
    row("scan", "Qb", 2, "b2", tw(BT_QB));
    row("scan", "Qm", 2, "b2", tw(BT_QM));
    row("scan", "Qm2", 2, "b2", tw(BT_QM2));
    row("scan", "Qp", 2, "b2", tw(BT_QP));
    row("scan", "Qmp", 2, "b2", tw(BT_QMP));
    row("scan", "Qz", 2, "b2", b2_width);
    row("scan", "C", 1, "none", 0);
    if (allow_pk) {
      row("scan", "PKM", 2, "b4", tw(BT_PKM));
      for (const char* nm : {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us", "Uls", "Ulss", "Urs", "S"})
        row("scan", nm, 4, "b4", beam);
      for (int t : {BT_QG, BT_QGLS, BT_QGL, BT_QGR, BT_QGRS})
        row("scan", BT_NAMES[t], 4, "b4", tw(t));
      if (rlt) {
        row("scan", "Rlt", 4, "rlt", rlt_beam);
        row("scan", "Ult", 4, "rlt", rlt_beam);
      }
    }
    if (pk_l2r == 1) {
      row("scan", "P1", 2, "b2", l2r_pbeam[0]);
      row("scan", "P2", 3, grammar_simple ? "b4" : "2b4", l2r_pbeam[1]);
      row("scan", "P3", 4, grammar_simple ? "b4" : "2b4", l2r_pbeam[2]);
      row("scan", "H", 3, "lp_hk", l2r_p_hk);
      if (l2r_p1 > 0) row("p1sweep", "states", 4, "p1b", l2r_p1_beam);
      if (grammar_simple) {
        for (const char* nm : {"R^M", "R^M2", "R^Multi", "U^M", "U^M2", "U^Multi"})
          row("scan", nm, 4, "b4", beam);
      }
    }
    auto dev = [&](const char* name, int on, const std::string& what) {
      std::fprintf(stderr, "[deviation] %s active %d beam 0 %s\n", name, on, what.c_str());
    };
    dev("max_gap_scan", arm_score ? 1 : 0,
        "arm-score window " + std::to_string(max_gap_scan) + "; unused when --arm-score is off");
    dev("max_arm_len", max_arm_len > 0 ? 1 : 0,
        "longest arm window; 0 means the sequence length");
    dev("gap_step", gap_step ? 1 : 0,
        "groups the b4 beam per step j instead of per (v,j); not a separate count");
    dev("gate", gcfg.on ? 1 : 0,
        "drops Qg/Qgl candidates by an access threshold; not a count by key");
    dev("arm_score", arm_score ? 1 : 0,
        "replaces the arm weight with an access score; not a count by key");
    dev("access", (arm_score || gcfg.on) ? 1 : 0,
        "access file for arm_score and the gate; not loaded in the L2R default");
    dev("loop_cap", 1, "energy-model max unpaired loop (" + std::to_string(loop_cap) +
                           "); not a count by key");
    dev("ml_lead_cap", ml_lead_cap > 0 ? 1 : 0,
        "leading unpaired run before the first multiloop branch is at most " +
            std::to_string(ml_lead_cap > 0 ? ml_lead_cap - 1 : n) +
            " bases; not a beam");
    dev("pair_bonus", pair_bonus != 0.0 ? 1 : 0,
        "added to every pair weight (value " + std::to_string(pair_bonus) + ")");
    dev("lK", k_exact ? 0 : 1, "K uses the ratio lK, not a bound");
    dev("V_miss", p2_vsweep ? 0 : 1,
        "a V miss falls back to the old V builder; with --v-sweep a miss is -inf");
    dev("v-floor", p2_v_floor ? 1 : 0,
        "fills a missing V entry with the smallest kept mu or log-weight 0");
    dev("v-cut-fb", p2_v_cut_fb ? 1 : 0,
        "a cut V entry (full map) is scored kappa_min / C_pre(b-1); ranking only");
    if (bad) {
      std::fprintf(stderr, "[beams] abort: a table with 2 or more indices has no active beam\n");
      return 2;
    }
    return 0;
  };

