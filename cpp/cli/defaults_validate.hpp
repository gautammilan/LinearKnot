// main(), part 4: derived defaults (--pk-l2r on, b1 / b2 / pbeam / b_V, sweep on), flag
// consistency checks, output files, beam audit.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  // --pk-l2r on defaults: --l2r-spawn inc; pbeam derived below as (b2,b4,b4);
  // --l2r-p3-hex/lin on; --pr-ml-span 100000; --bG 16 --grow-node-trunc-r -1; --pr-v-key xs when safe;
  // then --grow-r-from-v on if its requirements hold.
  if (pk_l2r == 1) {
    auto given = [&](const char* f) {
      for (int k = 1; k < argc; ++k)
        if (std::strcmp(argv[k], f) == 0) return true;
      return false;
    };
    if (!given("--l2r-spawn")) {
      l2r_spawn_inc = true;
      l2r_spawn_cube = false;
    }
    // pbeam filled after b2/b4 are known (see below).
    if (!given("--l2r-p3-hex")) l2r_p3_hex = true;
    if (!given("--l2r-p3-lin")) l2r_p3_lin = true;
    if (!given("--pr-ml-span")) pr_ml_span = 100000;
    if (!given("--bG")) g_ntr_K = 16;
    if (!given("--grow-node-trunc-r")) g_ntr_KR = -1;
    if (!given("--pr-v-key") && (pr_fast > 0))
      pr_vxs = true;
    if (!given("--grow-r-from-v")) {
      if ((g_ntr_K > 0) && (g_ntr_KR < 0) && pr_vxs && (pr_fast >= 2))
        g_rfv = true;
      else
        std::fprintf(stderr, "[defaults] --grow-r-from-v left off: needs --bG K > 0, "
                             "--grow-node-trunc-r -1, --pr-v-key xs, --pr-fast ml\n");
    }
  }
  // b1 and b2 default to 200, independent of b4.
  // Sweep and the V sweep are on. No length gate. Q^z and Q^z_pre use b2.
  if (!b1_set) nest_beam = 200;
  if (b2_at < 0) {
    const int b2tabs[] = {BT_QB, BT_QM, BT_QM2, BT_QMP, BT_QP};
    for (int t : b2tabs)
      if (table_beam[t] < 0) table_beam[t] = 200;
  }
  b2_width = table_beam[BT_QM] >= 0 ? table_beam[BT_QM] : 200;
  // pbeam = (b2, b4, b4).
  l2r_pbeam[0] = b2_width;
  l2r_pbeam[1] = beam;
  l2r_pbeam[2] = beam;
  // --bV default 256 (0 = uncapped).
  if (!bV_set) {
    pr_vcap = 256;
    p2_v_cap = 256;
  }
  std::fprintf(stderr, "[arity] b2 %d b4 %d pbeam %d,%d,%d\n", b2_width, beam, l2r_pbeam[0],
               l2r_pbeam[1], l2r_pbeam[2]);
  // Default sweep ON: sw_as_rn with q'-aware intermediate cap (= lazy trunc key).
  // Pass `--sweep off` for rec_node/pr_steps (hard14/OLD path).
  if (!sweep_set) p2_sweep = true;
  if (!vsweep_set) p2_vsweep = true;

  if (stem_loop < 0) {
    std::fprintf(stderr, "--stem-est-loop must be >= 0\n");
    return 2;
  }
  if (stem_qg_r_set && !stem_rank) {
    std::fprintf(stderr, "--stem-est-qg-r needs --stem-est rank\n");
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
  if (grow_set && stem_ext != 3) {
    std::fprintf(stderr, "--stem-est-grow-loop needs --stem-est-qg-ext loop\n");
    return 2;
  }
  if (grow_lg < 0) {
    std::fprintf(stderr, "--stem-est-grow-loop must be >= 0\n");
    return 2;
  }
  if (helpers_set && !stem_ext) {
    std::fprintf(stderr, "--stem-est-helpers needs --stem-est-qg-ext\n");
    return 2;
  }
  if (stem_ext != 3 || stem_helpers != 1) {
    std::fprintf(stderr, "the stem profile needs --stem-est-qg-ext loop and --stem-est-helpers grow\n");
    return 2;
  }
  if (!stem_pr || !(pr_abl & DP03::PA_SEP)) {
    std::fprintf(stderr, "pass-2 growth needs --stem-est-principled on and --pr-abl sep or sepv\n");
    return 2;
  }
  if (g_ml_c < 0) {
    std::fprintf(stderr, "--grow-ml-c must be >= 0\n");
    return 2;
  }
  if (g_tr_K < 0 || (g_tr_K > 0 && !pr_gsel)) {
    std::fprintf(stderr, "--grow-trunc K >= 0 needs --pr-grow-select bound\n");
    return 2;
  }
  if (g_ntr_K < 0 || (g_ntr_K > 0 && pk_l2r != 1)) {
    std::fprintf(stderr, "--bG K >= 0 needs --pk-l2r on\n");
    return 2;
  }
  if (pr_vxs && (pk_l2r != 1 || pr_fast <= 0)) {
    std::fprintf(stderr, "--pr-v-key xs needs --pk-l2r on and --pr-fast\n");
    return 2;
  }
  if (g_rfv && (g_ntr_K <= 0 || g_ntr_KR >= 0 || !pr_vxs || pr_fast < 2)) {
    std::fprintf(stderr, "--grow-r-from-v on needs --bG K > 0, --grow-node-trunc-r -1, "
                         "--pr-v-key xs, --pr-fast ml\n");
    return 2;
  }
  if (pr_vstream && pr_fast <= 0) {
    std::fprintf(stderr, "--pr-v-stream on needs --pr-fast\n");
    return 2;
  }
  g_skip_bml = (pr_abl & DP03::PA_BML) != 0;
  if (stem_pr && !g_skip_bml) {
    std::fprintf(stderr, "[unsupported] principled score with Bml is off\n");
    g_skip_bml = true;
  }
  if (stem_pr && (!stem_rank || !stem_qg2 || !stem_both)) {
    std::fprintf(stderr, "--stem-est-principled on needs --stem-est rank and --stem-est-role all\n");
    return 2;
  }
  if (pr_set && !stem_pr) {
    std::fprintf(stderr, "--pr-* need --stem-est-principled on\n");
    return 2;
  }
  if (!(pr_geps >= 0.0 && pr_geps < 1.0)) {
    std::fprintf(stderr, "--pr-grow-eps must be in [0, 1)\n");
    return 2;
  }
  if (!(pr_veps >= 0.0 && pr_veps < 1.0)) {
    std::fprintf(stderr, "--pr-v-eps must be in [0, 1)\n");
    return 2;
  }
  if (pr_gsel) pr_floor = 0.0;
  if (p2_v_cap > 0 && !stem_pr) {
    std::fprintf(stderr, "--bV needs --stem-est-principled on\n");
    return 2;
  }
  if (!(pr_floor >= 0.0 && pr_floor < 1.0) || !(pr_cut >= 0.0 && pr_cut < 1.0) || pr_ml_span < 0) {
    std::fprintf(stderr, "--pr-cut must be in [0, 1), --pr-ml-span >= 0\n");
    return 2;
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
  // Per-pass caps default to --loop-cap.
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
      nest_tables("pass1", nest_beam);
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
    row("scan", "PKM", 2, "b4", tw(BT_PKM));
    for (const char* nm : {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us", "Uls", "Ulss", "Urs", "S"})
      row("scan", nm, 4, "b4", beam);
    for (int t : {BT_QG, BT_QGLS, BT_QGL, BT_QGR, BT_QGRS})
      row("scan", BT_NAMES[t], 4, "b4", tw(t));
    if (pk_l2r == 1) {
      row("scan", "P1", 2, "b2", l2r_pbeam[0]);
      row("scan", "P2", 3, grammar_simple ? "b4" : "2b4", l2r_pbeam[1]);
      row("scan", "P3", 4, grammar_simple ? "b4" : "2b4", l2r_pbeam[2]);
      row("scan", "H", 3, "lp_hk", l2r_p_hk);
      if (grammar_simple) {
        for (const char* nm : {"R^M", "R^M2", "R^Multi", "U^M", "U^M2", "U^Multi"})
          row("scan", nm, 4, "b4", beam);
      }
    }
    auto dev = [&](const char* name, int on, const std::string& what) {
      std::fprintf(stderr, "[deviation] %s active %d beam 0 %s\n", name, on, what.c_str());
    };
    dev("gap_step", gap_step ? 1 : 0,
        "groups the b4 beam per step j instead of per (v,j); not a separate count");
    dev("loop_cap", 1, "energy-model max unpaired loop (" + std::to_string(loop_cap) +
                           "); not a count by key");
    dev("ml_lead_cap", ml_lead_cap > 0 ? 1 : 0,
        "leading unpaired run before the first multiloop branch is at most " +
            std::to_string(ml_lead_cap > 0 ? ml_lead_cap - 1 : n) +
            " bases; not a beam");
    dev("pair_bonus", pair_bonus != 0.0 ? 1 : 0,
        "added to every pair weight (value " + std::to_string(pair_bonus) + ")");
    dev("lK", 1, "K uses the ratio lK, not a bound");
    dev("V_miss", p2_vsweep ? 0 : 1,
        "a V miss falls back to the old V builder; with --v-sweep a miss is -inf");
    if (bad) {
      std::fprintf(stderr, "[beams] abort: a table with 2 or more indices has no active beam\n");
      return 2;
    }
    return 0;
  };

