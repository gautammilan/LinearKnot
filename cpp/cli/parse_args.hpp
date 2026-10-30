// main(), part 3: the argument loop (every --flag).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  for (int a = 1; a < argc; ++a) {
    std::string arg = argv[a];
    auto next = [&]() -> const char* {
      if (a + 1 >= argc) {
        usage();
        std::exit(2);
      }
      return argv[++a];
    };
    if (arg == "--beam" || arg.rfind("--beam=", 0) == 0) {
      std::fprintf(stderr, "--beam is removed; use --b4\n");
      return 2;
    }
    else if (arg == "--b4") {
      std::string b = next();
      beam = (b == "inf" || b == "infinity" || b == "none") ? 0 : std::atoi(b.c_str());
    }
    else if (arg == "--loop-cap") {
      loop_cap = std::atoi(next());
      loop_cap_set = true;
    }
    else if (arg == "--loop-cap-p1" || arg == "--loop-cap-p2" || arg == "--loop-cap-scan") {
      // Temporary: per-pass --loop-cap / --max-loop. "off"/"n"/0 = uncapped.
      auto parse_pc = [](const std::string& v) -> int {
        if (v == "n" || v == "off" || v == "0") return NO_LOOP_CAP;
        return std::atoi(v.c_str());
      };
      const int v = parse_pc(next());
      if (arg == "--loop-cap-p1") loop_cap_p1 = v;
      else if (arg == "--loop-cap-p2") loop_cap_p2 = v;
      else loop_cap_scan = v;
    }
    else if (arg == "--ml-lead-cap") {
      std::string v = next();
      if (v == "n" || v == "off") ml_lead_cap = 0;
      else ml_lead_cap = std::atoi(v.c_str());
    }
    else if (arg == "--max-loop") max_loop = std::atoi(next());
    else if (arg == "--no-gu") allow_gu = false;
    else if (arg == "--no-single-stems") single_stems = false;
    else if (arg == "--wc-ends") wc_ends = true;
    else if (arg == "--pair-bonus") pair_bonus = std::atof(next());
    else if (arg == "--restrict-pairs") restrict_pairs = true;
    else if (arg == "--rt") RT = std::atof(next());
    else if (arg == "--score") {
      std::string s = next();
      if (s == "qz") score_qz = true;
      else if (s == "prefix") score_qz = false;
      else {
        std::fprintf(stderr, "score must be 'qz' or 'prefix'\n");
        return 2;
      }
    } else if (arg == "--no-pk") allow_pk = false;
    else if (arg == "--cutoff") cutoff = std::atof(next());
    else if (arg == "--no-bpp") do_bpp = false;
    else if (arg == "--arm-score") {
      std::string s = next();
      if (s == "on") arm_score = true;
      else if (s == "off") arm_score = false;
      else {
        std::fprintf(stderr, "--arm-score must be 'on' or 'off'\n");
        return 2;
      }
    } else if (arg == "--access-file") access_file = next();
    else if (arg == "--max-gap-scan") max_gap_scan = std::atoi(next());
    else if (arg == "--max-arm-len") max_arm_len = std::atoi(next());
    else if (arg == "--gap-beam") {
      std::string m = next();
      if (m == "group") gap_step = false;
      else if (m == "step") gap_step = true;
      else {
        std::fprintf(stderr, "--gap-beam must be 'group' or 'step'\n");
        return 2;
      }
    }
    else if (arg == "--params") {
      std::string m = next();
      if (m == "dp03") PK = DP03_PARAMS;
      else if (m == "dp09") PK = DP09_PARAMS;
      else if (m == "dp09st") PK = DP09ST_PARAMS;
      else if (m == "dp09ml") PK = DP09ML_PARAMS;
      else if (m == "dp09full") PK = DP09FULL_PARAMS;
      else {
        std::fprintf(stderr, "--params must be dp03, dp09, dp09st, dp09ml or dp09full\n");
        return 2;
      }
    }
    else if (arg == "--material") {
      std::string m = next();
      mat_rna1999 = m == "rna1999";
      mat_v04 = m == "vienna2004";
      if (m != "dp03" && !mat_rna1999 && !mat_v04) {
        std::fprintf(stderr, "--material must be dp03, rna1999 or vienna2004\n");
        return 2;
      }
    }
    else if (arg == "--param-file") param_file = next();
    else if (arg == "--dangles") {
      std::string m = next();
      if (m == "none") dangles_mode = 0;
      else if (m == "some") dangles_mode = 1;
      else if (m == "all") dangles_mode = 2;
      else if (m == "d2") dangles_mode = 3;
      else {
        std::fprintf(stderr, "--dangles must be none, some, all or d2\n");
        return 2;
      }
    }
    else if (arg == "--gate") {
      std::string s = next();
      if (s == "on") gcfg.on = true;
      else if (s == "off") gcfg.on = false;
      else {
        std::fprintf(stderr, "--gate must be 'on' or 'off'\n");
        return 2;
      }
    }
    else if (arg == "--gate-tau") gcfg.tau = std::strtod(next(), nullptr);
    else if (arg == "--gate-span") gcfg.span = std::atoi(next());
    else if (arg == "--gate-delta") gcfg.delta = std::strtod(next(), nullptr);
    else if (arg == "--gate-role") {
      std::string s = next();
      if (s == "stem1") gcfg.both = false;
      else if (s == "both") gcfg.both = true;
      else {
        std::fprintf(stderr, "--gate-role must be 'stem1' or 'both'\n");
        return 2;
      }
    }
    else if (arg == "--gate-tables") {
      std::string s = next();
      if (s == "qgl") gcfg.tables[GATE_QG] = false, gcfg.tables[GATE_QGL] = true;
      else if (s == "qg") gcfg.tables[GATE_QG] = true, gcfg.tables[GATE_QGL] = false;
      else if (s == "qg,qgl" || s == "qgl,qg") gcfg.tables[GATE_QG] = gcfg.tables[GATE_QGL] = true;
      else {
        std::fprintf(stderr, "--gate-tables must be qg, qgl or qg,qgl\n");
        return 2;
      }
    }
    else if (arg == "--gate-reserve") gcfg.reserve = std::atoi(next());
    else if (arg == "--gate-stats") gate_stats_out = true;
    else if (arg == "--stem-est") {
      std::string s = next();
      if (s == "rank") stem_rank = true;
      else if (s == "off") stem_rank = false;
      else {
        std::fprintf(stderr, "--stem-est must be 'off' or 'rank'\n");
        return 2;
      }
    }
    else if (arg == "--stem-est-loop") stem_loop = std::atoi(next());
    else if (arg == "--stem-est-role") {
      std::string s = next();
      stem_both = s == "both" || s == "all";
      stem_qg2 = s == "all";
      if (s != "stem1" && s != "both" && s != "all") {
        std::fprintf(stderr, "--stem-est-role must be 'stem1', 'both' or 'all'\n");
        return 2;
      }
    }
    else if (arg == "--stem-est-qg-r") {
      std::string s = next();
      stem_qg_r_set = true;
      stem_qg_partner = s == "partner";
      if (s != "mirror" && s != "partner") {
        std::fprintf(stderr, "--stem-est-qg-r must be 'mirror' or 'partner'\n");
        return 2;
      }
    }
    else if (arg == "--stem-est-consts") {
      std::string s = next();
      stem_consts_set = true;
      stem_consistent = s == "consistent";
      if (s != "as-is" && s != "consistent") {
        std::fprintf(stderr, "--stem-est-consts must be 'as-is' or 'consistent'\n");
        return 2;
      }
    }
    else if (arg == "--stem-est-qg-ext") {
      std::string s = next();
      stem_ext_set = true;
      stem_ext = s == "stack" ? 1 : s == "bulge1" ? 2 : s == "loop" ? 3 : 0;
      if (s != "off" && s != "stack" && s != "bulge1" && s != "loop") {
        std::fprintf(stderr, "--stem-est-qg-ext must be 'off', 'stack', 'bulge1' or 'loop'\n");
        return 2;
      }
    }
    else if (arg == "--stem-est-grow-loop") {
      grow_set = true;
      grow_lg = std::atoi(next());
    }
    else if (arg == "--stem-est-grow-cut") {
      grow_set = true;
      grow_cut = std::atof(next());
    }
    else if (arg == "--stem-est-grow-nodes") {
      grow_set = true;
      grow_nodes = std::atoi(next());
    }
    else if (arg == "--stem-est-helpers") {
      std::string s = next();
      helpers_set = true;
      stem_helpers = s == "grow" ? 1 : s == "nested" ? 2 : 0;
      if (s != "off" && s != "grow" && s != "nested") {
        std::fprintf(stderr, "--stem-est-helpers must be 'off', 'grow' or 'nested'\n");
        return 2;
      }
    }
    else if (arg == "--stem-est-ml-span") {
      helpers_set = true;
      ml_span = std::atoi(next());
    }
    else if (arg == "--stem-est-profile") {
      std::string s = next();
      stem_profile = s == "on";
      if (s != "off" && s != "on") {
        std::fprintf(stderr, "--stem-est-profile must be 'off' or 'on'\n");
        return 2;
      }
    }
    else if (arg == "--stem-est-profile-eps") {
      profile_set = true;
      prof_eps = std::strtod(next(), nullptr);
    }
    else if (arg == "--stem-est-ml-tail") {
      std::string s = next();
      profile_set = true;
      ml_tail = s == "nested";
      if (s != "off" && s != "nested") {
        std::fprintf(stderr, "--stem-est-ml-tail must be 'off' or 'nested'\n");
        return 2;
      }
    }
    else if (arg == "--stem-profile-check") prof_check = true;
    else if (arg == "--stem-est-bg-reorder") {
      std::string s = next();
      bg_reorder = s == "on";
      if (s != "off" && s != "on") {
        std::fprintf(stderr, "--stem-est-bg-reorder must be 'off' or 'on'\n");
        return 2;
      }
    }
    else if (arg == "--stem-bg-reorder-check") ro_check = true;
    else if (arg == "--stem-est-chain" || arg == "--stem-est-chain-ml") {
      std::string s = next();
      if (arg == "--stem-est-chain-ml") chain_set = true;
      (arg == "--stem-est-chain" ? stem_chain : chain_ml) = s == "on";
      if (s != "off" && s != "on") {
        std::fprintf(stderr, "%s must be 'off' or 'on'\n", arg.c_str());
        return 2;
      }
    }
    else if (arg == "--stem-est-chain-floor") {
      chain_set = true;
      chain_floor = std::strtod(next(), nullptr);
    }
    else if (arg == "--stem-est-chain-span") {
      chain_set = true;
      chain_span = std::atoi(next());
    }
    else if (arg == "--stem-est-chain-loop") {
      chain_set = true;
      chain_loop = std::atoi(next());
    }
    else if (arg == "--stem-est-chain-cut") {
      chain_set = true;
      chain_cut = std::strtod(next(), nullptr);
    }
    else if (arg == "--stem-est-hr") {
      std::string s = next();
      stem_hr_set = true;
      stem_hr_av = s == "av";
      stem_hr_phi = s == "phi";
      if (s != "orig" && s != "av" && s != "phi") {
        std::fprintf(stderr, "--stem-est-hr must be 'orig', 'av' or 'phi'\n");
        return 2;
      }
    }
    else if (arg == "--stem-est-hr-gamma") {
      stem_hr_set = true;
      stem_hr_gamma = std::strtod(next(), nullptr);
    }
    else if (arg == "--stem-hr-check") stem_hr_check = true;
    else if (arg == "--stem-est-principled" || arg == "--stem-est-principled-nested") {
      std::string s = next();
      (arg == "--stem-est-principled" ? stem_pr : stem_pr_nested) = s == "on";
      if (s != "off" && s != "on") {
        std::fprintf(stderr, "%s must be 'off' or 'on'\n", arg.c_str());
        return 2;
      }
    }
    else if (arg == "--pr-floor" || arg == "--pr-cut") {
      pr_set = true;
      (arg == "--pr-floor" ? pr_floor : pr_cut) = std::strtod(next(), nullptr);
    }
    else if (arg == "--pr-span" || arg == "--pr-ml-span" || arg == "--pr-ml-close" ||
             arg == "--pr-arms-k") {
      pr_set = true;
      (arg == "--pr-span" ? pr_span : arg == "--pr-ml-span" ? pr_ml_span
                                    : arg == "--pr-ml-close" ? pr_ml_close : pr_arms_k) =
          std::atoi(next());
    }
    else if (arg == "--stem-pr-check") pr_check = true;
    else if (arg == "--pr-arms-rank" || arg == "--pr-cut-by" || arg == "--pr-cut-roots") {
      pr_set = true;
      const std::string s = next();
      const char* a = arg == "--pr-arms-rank" ? "ratio" : arg == "--pr-cut-by" ? "gapfree" : "on";
      const char* b = arg == "--pr-arms-rank" ? "phi" : arg == "--pr-cut-by" ? "weight" : "off";
      if (s != a && s != b) {
        std::fprintf(stderr, "%s must be '%s' or '%s'\n", arg.c_str(), b, a);
        return 2;
      }
      (arg == "--pr-arms-rank" ? pr_ratio : arg == "--pr-cut-by" ? pr_gapfree : pr_cut_roots) = s == a;
      if (arg == "--pr-cut-by") pr_cut_by_set = true;
    }
    else if (arg == "--pr-arms" || arg == "--pr-arm-z" || arg == "--pr-test-qz") {
      pr_set = true;
      const std::string s = next();
      const char* a = arg == "--pr-arms" ? "sep" : "unpaired";
      const char* b = arg == "--pr-arms" ? "exact" : arg == "--pr-arm-z" ? "qz" : "off";
      if (s != a && s != b) {
        std::fprintf(stderr, "%s must be '%s' or '%s'\n", arg.c_str(), b, a);
        return 2;
      }
      if (arg == "--pr-arms") pr_arms_mode = s == a ? 1 : 0;
      else if (arg == "--pr-arm-z") pr_arm_z = s == a ? 1 : 0;
      else pr_qz_unp = s == a;
    }
    else if (arg == "--pr-partner-floor") {
      pr_set = true;
      pr_partner_floor = std::strtod(next(), nullptr);
      if (!(pr_partner_floor >= 0.0 && pr_partner_floor < 1.0)) {
        std::fprintf(stderr, "--pr-partner-floor must be in [0, 1)\n");
        return 2;
      }
    }
    else if (arg == "--pr-sep-check") {
      pr_set = true;
      pr_sep_check = true;
    }
    else if (arg == "--pr-grow-select") {
      pr_set = true;
      const std::string s = next();
      if (s != "floor" && s != "bound") {
        std::fprintf(stderr, "--pr-grow-select must be 'floor' or 'bound'\n");
        return 2;
      }
      pr_gsel = s == "bound";
    }
    else if (arg == "--pr-grow-eps") {
      pr_set = true;
      pr_geps = std::strtod(next(), nullptr);
    }
    else if (arg == "--pr-grow-B") {
      pr_set = true;
      pr_gB = std::atoi(next());
    }
    else if (arg == "--pr-grow-check") {
      pr_set = true;
      pr_gcheck = true;
    }
    else if (arg == "--grow-rec") {
      const std::string s = next();
      if (s != "on" && s != "off") {
        std::fprintf(stderr, "--grow-rec must be on or off\n");
        return 2;
      }
      g_grow_rec = s == "on";
    }
    else if (arg == "--grow-ml-L" || arg == "--grow-ml-R") {
      const std::string s = next();
      const int m = s == "all" ? 2 : s == "depth1" ? 1 : s == "none" ? 0 : -1;
      if (m < 0) {
        std::fprintf(stderr, "%s must be all, depth1 or none\n", arg.c_str());
        return 2;
      }
      (arg == "--grow-ml-L" ? g_ml_L : g_ml_R) = m;
    }
    else if (arg == "--grow-ml-c") {
      g_ml_c = std::atoi(next());
    }
    else if (arg == "--grow-ml-ub") {
      g_ml_ub = std::strcmp(next(), "on") == 0;
    }
    else if (arg == "--grow-ml-diag") {
      g_ml_diag = true;
    }
    else if (arg == "--grow-rec-check") {
      g_rec_check = true;
    }
    else if (arg == "--pr-U-nobml") {
      g_U_nobml = true;
    }
    else if (arg == "--grow-trunc") {
      g_tr_K = std::atoi(next());
    }
    else if (arg == "--grow-trunc-eps") {
      g_tr_eps = std::strtod(next(), nullptr);
    }
    else if (arg == "--grow-trunc-mode") {
      const std::string s = next();
      if (s != "side" && s != "pool") {
        std::fprintf(stderr, "--grow-trunc-mode must be side or pool\n");
        return 2;
      }
      g_tr_pool = s == "pool";
    }
    else if (arg == "--pk-l2r") {
      const std::string s = next();
      if (s != "on" && s != "off" && s != "check") {
        std::fprintf(stderr, "--pk-l2r must be on, off or check\n");
        return 2;
      }
      pk_l2r = s == "on" ? 1 : s == "check" ? 2 : 0;
    }
    else if (arg == "--pk-l2r-checkR") {
      pk_l2r_R = true;
    }
    else if (arg == "--l2r-beam") {
      const std::string b = next();
      l2r_beam = (b == "inf" || b == "none") ? 0 : std::atoi(b.c_str());
    }
    else if (arg == "--l2r-diag") {
      l2r_diag = true;
    }
    else if (arg == "--l2r-tbeam") {
      l2r_tbeam = next();
    }
    else if (arg == "--l2r-spawn-k") {
      const std::string b = next();
      l2r_spawn_k = b == "all" ? 0 : std::atoi(b.c_str());
    }
    else if (arg == "--l2r-spawn") {
      const std::string v = next();
      if (v != "all" && v != "cube" && v != "cube-static" && v != "inc") {
        std::fprintf(stderr, "--l2r-spawn must be all, cube, cube-static or inc\n");
        return 2;
      }
      l2r_spawn_cube = v == "cube" || v == "cube-static";
      l2r_spawn_static = v == "cube-static";
      l2r_spawn_inc = v == "inc";
    }
    else if (arg == "--l2r-p3-split") {
      if (std::sscanf(next(), "%d,%d", &l2r_p3_split[0], &l2r_p3_split[1]) != 2 || l2r_p3_split[0] < 1 ||
          l2r_p3_split[1] < 1) {
        std::fprintf(stderr, "--l2r-p3-split needs Ba,bh >= 1\n");
        return 2;
      }
    }
    else if (arg == "--l2r-pbeam") {
      std::stringstream ss(next());
      std::string item;
      int k = 0;
      while (std::getline(ss, item, ',') && k < 3)
        l2r_pbeam[k++] = (item == "inf" || item == "none") ? 0 : std::atoi(item.c_str());
      if (k != 3) {
        std::fprintf(stderr, "--l2r-pbeam needs b1,b2,b3\n");
        return 2;
      }
      pbeam_explicit_at = a;
    }
    else if (arg == "--l2r-p-hk") {
      l2r_p_hk = std::atoi(next());
      if (l2r_p_hk < 1) {
        std::fprintf(stderr, "--l2r-p-hk must be >= 1\n");
        return 2;
      }
    }
    else if (arg == "--l2r-p-aopen") {
      const std::string v = next();
      if (v != "q" && v != "raw") {
        std::fprintf(stderr, "--l2r-p-aopen must be q or raw\n");
        return 2;
      }
      l2r_p_raw = v == "raw";
    }
    else if (arg == "--l2r-p3-hex") {
      l2r_p3_hex = std::string(next()) == "on";
    }
    else if (arg == "--l2r-p3-lin") {
      l2r_p3_lin = std::string(next()) == "on";
    }
    else if (arg == "--l2r-inc-check") {
      l2r_inc_check = std::string(next()) == "on";
    }
    else if (arg == "--l2r-ninf-mass") {
      l2r_ninf_mass = std::string(next()) == "on";
    }
    else if (arg == "--l2r-keep-ninf") {
      l2r_keep_ninf = std::string(next()) == "on";
    }
    else if (arg == "--l2r-qz-audit") {
      l2r_qz_audit = std::string(next()) == "on";
    }
    else if (arg == "--grow-node-trunc-z") {
      const std::string v = next();
      if (v != "scan" && v != "bar") {
        std::fprintf(stderr, "--grow-node-trunc-z must be scan or bar\n");
        return 2;
      }
      ntr_z_set = true;
      g_ntr_zbar = v == "bar";
    }
    else if (arg == "--l2r-spawn-w") {
      l2r_spawn_wk = std::string(next()) != "suf";
    }
    else if (arg == "--l2r-spawn-pool") {
      l2r_spawn_pool = std::atoi(next());
    }
    else if (arg == "--l2r-spawn-maxpop") {
      l2r_spawn_maxpop = std::atoi(next());
    }
    else if (arg == "--l2r-trace") {
      l2r_trace = next();
    }
    else if (arg == "--l2r-perf") {
      l2r_perf = std::strcmp(next(), "on") == 0;
    }
    else if (arg == "--l2r-trace-key") {
      l2r_trace_key = next();
    }
    else if (arg == "--l2r-oracle-dump") {
      l2r_odump = next();
    }
    else if (arg == "--l2r-oracle") {
      l2r_oload = next();
    }
    else if (arg == "--l2r-sub") {
      l2r_sub = next();
    }
    else if (arg == "--l2r-ctx") {
      const std::string s = next();
      if (s != "off" && s != "on" && s != "p1") {
        std::fprintf(stderr, "--l2r-ctx must be off, on or p1\n");
        return 2;
      }
      if (s != "off") {
        std::fprintf(stderr, "[unsupported] --l2r-ctx %s is off; it builds Bml\n", s.c_str());
        return 2;
      }
      l2r_ctx = false;
      l2r_ctx_p1 = false;
    }
    else if (arg == "--l2r-ctx-tau") {
      l2r_ctx_tau = std::atof(next());
    }
    else if (arg == "--l2r-ctx-rounds") {
      l2r_ctx_rounds = std::max(1, std::atoi(next()));
    }
    else if (arg == "--l2r-kx") {
      l2r_kx = std::string(next()) == "on";
    }
    else if (arg == "--l2r-p1") {
      l2r_p1 = std::atoi(next());
    }
    else if (arg == "--l2r-p1-qb") {
      l2r_p1_own = std::string(next()) == "pass1";
    }
    else if (arg == "--l2r-p1-beam") {
      l2r_p1_beam = std::atoi(next());
    }
    else if (arg == "--l2r-split") {
      l2r_split = std::string(next()) == "on";
    }
    else if (arg == "--pk-dbeta1" || arg == "--pk-dbeta1m" || arg == "--pk-dbeta1p") {
      (arg == "--pk-dbeta1" ? pk_db[0] : arg == "--pk-dbeta1m" ? pk_db[1] : pk_db[2]) = std::atof(next());
    }
    else if (arg == "--pr-score-dump") {
      pr_set = true;
      pr_dump_path = next();
    }
    else if (arg == "--pr-armq-sep") {
      pr_set = true;
      pr_aq_sep = std::string(next()) == "on";
    }
    else if (arg == "--pr-v-ratio") {
      pr_set = true;
      pr_v_ratio = std::string(next()) == "on";
    }
    else if (arg == "--pr-v-diag") {
      pr_set = true;
      pr_v_diag = next();
    }
    else if (arg == "--pr-v-ratio2") {
      pr_set = true;
      g_v2ratio = std::string(next()) == "on";
      if (g_v2ratio) pr_v_ratio = true;
    }
    else if (arg == "--pr-v2-dump") {
      pr_set = true;
      g_v2dump = next();
    }
    else if (arg == "--pr-static-only") {
      pr_set = true;
      g_static_only = true;
    }
    else if (arg == "--pr-v-sample") {
      pr_set = true;
      pr_v_sample = std::strtod(next(), nullptr);
    }
    else if (arg == "--pr-aq-max") {
      pr_set = true;
      pr_aq_max = static_cast<size_t>(std::max(1, std::atoi(next())));
    }
    else if (arg == "--pr-arm-persist" || arg == "--pr-arm-gaponly") {
      pr_set = true;
      const std::string s = next();
      if (s != "on" && s != "off") {
        std::fprintf(stderr, "%s on|off\n", arg.c_str());
        return 2;
      }
      (arg == "--pr-arm-persist" ? pr_arm_persist : pr_arm_gaponly) = s == "on";
    }
    else if (arg == "--pr-partner-maxpairs") {
      pr_set = true;
      const std::string s = next();
      pr_pmax = s == "off" ? 0 : std::atoi(s.c_str());
      if (pr_pmax < 0 || pr_pmax > 255) {
        std::fprintf(stderr, "--pr-partner-maxpairs must be off or 1..255\n");
        return 2;
      }
    }
    else if (arg == "--pr-partner-check") {
      pr_set = true;
      pr_pcheck = true;
    }
    else if (arg == "--pr-partner-ml") {
      pr_set = true;
      const std::string s = next();
      if (s != "on" && s != "off") {
        std::fprintf(stderr, "--pr-partner-ml must be on or off\n");
        return 2;
      }
      pr_pml = s == "on";
    }
    else if (arg == "--pr-fast") {
      pr_set = true;
      const std::string s = next();
      if (s != "off" && s != "acc" && s != "ml") {
        std::fprintf(stderr, "--pr-fast must be off, acc or ml\n");
        return 2;
      }
      pr_fast = s == "off" ? 0 : s == "acc" ? 1 : 2;
    }
    else if (arg == "--pr-abl") {
      pr_set = true;
      std::stringstream ss(next());
      std::string t;
      while (std::getline(ss, t, ',')) {
        if (t == "bml") pr_abl |= DP03::PA_BML;
        else if (t == "q") pr_abl |= DP03::PA_Q;
        else if (t == "p") pr_abl |= DP03::PA_P;
        else if (t == "xs") pr_abl |= DP03::PA_XS;
        else if (t == "ks") pr_abl |= DP03::PA_KS;
        else if (t == "xb") pr_abl |= DP03::PA_XB;
        else if (t == "qr") pr_abl |= DP03::PA_QR;
        else if (t == "sep") pr_abl |= DP03::PA_SEP | DP03::PA_BML | DP03::PA_Q;
        else if (t == "sepv") pr_abl |= DP03::PA_SEP | DP03::PA_BML;
        else if (t != "none") {
          std::fprintf(stderr, "--pr-abl: comma list of bml, q, p, xs, xb, ks, qr, sep, sepv, none\n");
          return 2;
        }
      }
    }
    else if (arg == "--pr-abl-pk") {
      pr_set = true;
      pr_abl_pk = std::stoi(next());
    }
    else if (arg == "--pr-partner-ml-check") {
      pr_set = true;
      pr_pmlcheck = true;
    }
    else if (arg == "--pr-partner-span") {
      pr_set = true;
      pr_pspan = std::atoi(next());
    }
    else if (arg == "--pr-v-eps") {
      pr_set = true;
      pr_veps = std::atof(next());
    }
    else if (arg == "--pr-v-trim") {
      pr_set = true;
      const std::string s = next();
      if (s != "global" && s != "suffix") {
        std::fprintf(stderr, "--pr-v-trim must be global or suffix\n");
        return 2;
      }
      pr_vsuffix = s == "suffix";
    }
    else if (arg == "--pr-v-cap" || arg.rfind("--pr-v-cap=", 0) == 0) {
      std::fprintf(stderr, "--pr-v-cap is removed; use --bV\n");
      return 2;
    }
    else if (arg == "--bV") {
      pr_set = true;
      bV_set = true;
      const int n = std::max(0, std::atoi(next()));
      pr_vcap = n;
      p2_v_cap = n;  // one flag: xs-keyed cap and C_pre left-end cap
    }
    else if (arg == "--pr-v-mlsteps") {
      pr_vmlsteps = std::string(next()) == "on";
    }
    else if (arg == "--pr-v-mlk") {
      pr_vmlk = std::max(0, std::atoi(next()));
    }
    else if (arg == "--pr-v-key") {
      pr_set = true;
      const std::string s = next();
      if (s != "xsys" && s != "xs") {
        std::fprintf(stderr, "--pr-v-key must be xsys or xs\n");
        return 2;
      }
      pr_vxs = s == "xs";
    }
    else if (arg == "--pr-v-stream") {
      pr_set = true;
      pr_vstream = std::string(next()) == "on";
    }
    else if (arg == "--grow-node-trunc" || arg.rfind("--grow-node-trunc=", 0) == 0) {
      // Exact --grow-node-trunc only (not -r / -z).
      std::fprintf(stderr, "--grow-node-trunc is removed; use --bG\n");
      return 2;
    }
    else if (arg == "--bG") {
      g_ntr_K = std::atoi(next());
    }
    else if (arg == "--g-cut") {
      const std::string s = next();
      if (s == "topq") g_cut = G_CUT_TOPQ;
      else if (s == "collapse") g_cut = G_CUT_COLLAPSE;
      else if (s == "entries") g_cut = G_CUT_ENTRIES;
      else {
        std::fprintf(stderr, "--g-cut must be topq, collapse, or entries\n");
        return 2;
      }
    }
    else if (arg == "--g-sweep") {
      const std::string s = next();
      if (s == "direct") g_sweep = G_SWEEP_DIRECT;
      else if (s == "states") g_sweep = G_SWEEP_STATES;
      else {
        std::fprintf(stderr, "--g-sweep must be direct or states\n");
        return 2;
      }
    }
    else if (arg == "--out-id") g_out_id = true;
    else if (arg == "--g-cut-at") {
      const std::string s = next();
      if (s == "pair") g_cut_pair_only = true;
      else if (s == "all") g_cut_pair_only = false;
      else {
        std::fprintf(stderr, "--g-cut-at needs all or pair\n");
        return 2;
      }
    }
    else if (arg == "--g-trace") {
      const std::string s = next();
      if (std::sscanf(s.c_str(), "%d,%d,%d,%d", &g_trace_p, &g_trace_q, &g_trace_x, &g_trace_y) < 2) {
        std::fprintf(stderr, "--g-trace needs P,Q[,X,Y]\n");
        return 2;
      }
    }
    else if (arg == "--grow-node-trunc-r") {
      g_ntr_KR = std::max(-1, std::atoi(next()));
    }
    else if (arg == "--grow-r-from-v") {
      g_rfv = std::string(next()) == "on";
    }
    else if (arg == "--malloc-trim") {
      g_mtrim = std::string(next()) == "on";
    }
    else if (arg == "--l2r-arm-fb") {
      l2r_arm_fb = std::string(next()) == "on";
    }
    else if (arg == "--l2r-lean") {
      l2r_lean = std::string(next()) == "on";
    }
    else if (arg == "--mem-log") {
      g_mem_log = std::string(next()) == "on";
    }
    else if (arg == "--pr-arms-eps") {
      pr_set = true;
      pr_arms_eps = std::strtod(next(), nullptr);
    }
    else if (arg == "--pr-arms-max") {
      pr_set = true;
      pr_arms_max = std::atoi(next());
    }
    else if (arg == "--pr-oracle") pr_oracle_path = next();
    else if (arg == "--pr-oracle-set") pr_oracle_set = next();
    else if (arg == "--stem-est-ext-max") {
      stem_ext_set = true;
      stem_ext_max = std::atoi(next());
    }
    else if (arg == "--stem-est-qg-combine") {
      std::string s = next();
      stem_ext_set = true;
      stem_split = s == "split";
      if (s != "sum" && s != "split") {
        std::fprintf(stderr, "--stem-est-qg-combine must be 'sum' or 'split'\n");
        return 2;
      }
    }
    else if (arg == "--stem-ext-check") stem_ext_check = true;
    else if (arg == "--stem-bg-check") stem_bg_check = true;
    else if (arg == "--stem-diag") stem_diag_path = next();
    else if (arg == "--stem-diag-loops") {
      stem_diag_loops.clear();
      for (double x : parse_list(next())) stem_diag_loops.push_back(static_cast<int>(x));
    }
    else if (arg == "--gate-diag") gate_diag_path = next();
    else if (arg == "--gate-diag-taus") diag_taus = parse_list(next());
    else if (arg == "--gate-long-span") long_span = std::atoi(next());
    else if (arg == "--beam-table") {
      std::string s = next();
      size_t eq = s.find('=');
      int bt = eq == std::string::npos ? -1 : beam_table_id(s.substr(0, eq));
      std::string b = eq == std::string::npos ? "" : s.substr(eq + 1);
      if (bt < 0 || b.empty()) {
        std::fprintf(stderr, "--beam-table needs T=B with T in Qb,Qm,Qm2,Qmp,PKM,Qp,Qg,Qgls,"
                             "Qgl,Qgr,Qgrs and B a beam or inf\n");
        return 2;
      }
      table_beam[bt] = (b == "inf" || b == "infinity" || b == "none") ? 0
                                                                      : std::max(0, std::atoi(b.c_str()));
    }
    else if (arg == "--dump-states") states_path = next();
    else if (arg == "--dump-scores") scores_path = next();
    else if (arg == "--beam-stats") beam_stats_path = next();
    else if (arg == "--dump-ext") ext_path = next();
    else if (arg == "--dump-min-w") dump_min_w = std::strtod(next(), nullptr);
    else if (arg == "--dump-pruned") pruned_path = next();
    else if (arg == "--dump-kept") kept_path = next();
    else if (arg == "--dump-exact") exact_path = next();
    else if (arg == "--dump-pruned-min-w") pruned_min_w = std::strtod(next(), nullptr);
    else if (arg == "--remove-state") {
      std::string s = next();
      size_t c = s.find(',');
      int bt = c == std::string::npos ? -1 : beam_table_id(s.substr(0, c));
      int v[4] = {-1, -1, -1, -1};
      if (bt >= 0) {
        std::vector<double> xs = parse_list(s.substr(c + 1));
        for (size_t k = 0; k < xs.size() && k < 4; ++k) v[k] = static_cast<int>(xs[k]);
      }
      if (bt < 0) {
        std::fprintf(stderr, "--remove-state needs T,a,b,c,d (keys as in --dump-states)\n");
        return 2;
      }
      remove_keys.insert(state_key(bt, v[0], v[1], v[2], v[3]));
    }
    else if (arg == "--qz-beam" || arg.rfind("--qz-beam=", 0) == 0) {
      std::fprintf(stderr, "--qz-beam is removed; Q^z uses b2 at every length\n");
      return 2;
    }
    else if (arg == "--qz-key" || arg.rfind("--qz-key=", 0) == 0) {
      std::fprintf(stderr, "--qz-key is removed; Q^z is ranked by C(k-1)+Q^z\n");
      return 2;
    }
    else if (arg == "--v-pairs" || arg.rfind("--v-pairs=", 0) == 0) {
      const std::string s = arg == "--v-pairs" ? next() : arg.substr(10);
      if (s == "P1" || s == "p1") {
        p2_v_mode = 1;
        p2_v_bp = 0;
      } else if ((s.rfind("own:", 0) == 0 || s.rfind("inside:", 0) == 0) &&
                 std::atoi(s.c_str() + (s[0] == 'o' ? 4 : 7)) > 0) {
        p2_v_mode = 2;
        p2_v_bp = std::atoi(s.c_str() + (s[0] == 'o' ? 4 : 7));
      } else if (s.rfind("prefix:", 0) == 0 && std::atoi(s.c_str() + 7) > 0) {
        p2_v_mode = 3;
        p2_v_bp = std::atoi(s.c_str() + 7);
      } else {
        std::fprintf(stderr, "--v-pairs must be P1, inside:b_P, prefix:b_P, or own:b_P\n");
        return 2;
      }
    }
    else if (arg == "--v-cap" || arg.rfind("--v-cap=", 0) == 0) {
      std::fprintf(stderr, "--v-cap is removed; use --bV\n");
      return 2;
    }
    else if (arg == "--v-trim" || arg.rfind("--v-trim=", 0) == 0) {
      const std::string s = arg == "--v-trim" ? next() : arg.substr(9);
      if (s == "both") p2_v_trim = 0;
      else if (s == "rel") p2_v_trim = 1;
      else if (s == "cap") p2_v_trim = 2;
      else {
        std::fprintf(stderr, "--v-trim must be rel, cap, or both\n");
        return 2;
      }
    }
    else if (arg == "--g-pairs" || arg.rfind("--g-pairs=", 0) == 0) {
      const std::string s = arg == "--g-pairs" ? next() : arg.substr(10);
      if (s == "P1" || s == "p1") {
        p2_g_mode = 1;
        p2_g_bp = 0;
      } else if ((s.rfind("own:", 0) == 0 || s.rfind("inside:", 0) == 0) &&
                 std::atoi(s.c_str() + (s[0] == 'o' ? 4 : 7)) > 0) {
        p2_g_mode = 2;
        p2_g_bp = std::atoi(s.c_str() + (s[0] == 'o' ? 4 : 7));
      } else if (s.rfind("prefix:", 0) == 0 && std::atoi(s.c_str() + 7) > 0) {
        p2_g_mode = 3;
        p2_g_bp = std::atoi(s.c_str() + 7);
      } else {
        std::fprintf(stderr, "--g-pairs must be P1, inside:b_P, prefix:b_P, or own:b_P\n");
        return 2;
      }
    }
    else if (arg == "--g-mode" || arg.rfind("--g-mode=", 0) == 0) {
      const std::string s = arg == "--g-mode" ? next() : arg.substr(9);
      if (s == "eager") p2_g_eager = true;
      else if (s == "lazy") p2_g_eager = false;
      else {
        std::fprintf(stderr, "--g-mode must be eager or lazy\n");
        return 2;
      }
    }
    else if (arg == "--g-cap" || arg.rfind("--g-cap=", 0) == 0) {
      const std::string s = arg == "--g-cap" ? next() : arg.substr(8);
      p2_g_cap = std::max(0, std::atoi(s.c_str()));
      if (p2_g_cap > 0) {
        g_ntr_K = p2_g_cap;
        p2_g_suffix = true;
      }
    }
    else if (arg == "--g-fallback" || arg.rfind("--g-fallback=", 0) == 0) {
      const std::string s = arg == "--g-fallback" ? next() : arg.substr(13);
      if (s == "onestep" || s == "on") p2_g_fb = true;
      else if (s == "off") p2_g_fb = false;
      else {
        std::fprintf(stderr, "--g-fallback must be onestep or off\n");
        return 2;
      }
    }
    else if (arg == "--profile") {
      profile = true;
    }
    else if (arg == "--g-mass" || arg.rfind("--g-mass=", 0) == 0) {
      const std::string s = arg == "--g-mass" ? next() : arg.substr(9);
      g_mass_f = std::fopen(s.c_str(), "wb");
      if (!g_mass_f) {
        std::fprintf(stderr, "--g-mass cannot open %s\n", s.c_str());
        return 2;
      }
    }
    else if (arg == "--g-key" || arg.rfind("--g-key=", 0) == 0) {
      const std::string s = arg == "--g-key" ? next() : arg.substr(8);
      if (s != "suffix") {
        std::fprintf(stderr, "--g-key must be suffix\n");
        return 2;
      }
      p2_g_suffix = true;
    }
    else if (arg == "--sweep-check") sw_check = true;
    else if (arg == "--sweep") {
      sweep_set = true;
      if (a + 1 < argc && (std::string(argv[a + 1]) == "off" || std::string(argv[a + 1]) == "on"))
        p2_sweep = std::string(argv[++a]) == "on";
      else
        p2_sweep = true;
    }
    else if (arg == "--v-sweep") {
      vsweep_set = true;
      if (a + 1 < argc && (std::string(argv[a + 1]) == "off" || std::string(argv[a + 1]) == "on"))
        p2_vsweep = std::string(argv[++a]) == "on";
      else
        p2_vsweep = true;
    }
    else if (arg == "--k-memo") k_memo = true;
    else if (arg == "--k-exact") k_exact = true;
    else if (arg == "--k-bound") k_bound = true;
    else if (arg == "--k-exact-ub") k_exact_ub = true;
    else if (arg == "--k-exact-check") k_exact_check = true;
    else if (arg == "--rlt") rlt = true;
    else if (arg == "--pass1-only") pass1_only = true;
    else if (arg == "--c-trace") c_trace = true;
    else if (arg == "--grammar") {
      if (a + 1 >= argc) {
        std::fprintf(stderr, "--grammar needs full or simple\n");
        return 2;
      }
      const std::string g = argv[++a];
      if (g == "simple") grammar_simple = true;
      else if (g == "full") grammar_simple = false;
      else {
        std::fprintf(stderr, "--grammar must be full or simple\n");
        return 2;
      }
    }
    else if (arg == "--next-pair") next_pair = true;
    else if (arg == "--rlt-beam") {
      rlt = true;
      rlt_beam = std::max(0, std::atoi(next()));
      rlt_beam_set = true;
    }
    else if (arg == "--sample-prof") sample_prof = next();
    else if (arg == "--v-floor" || arg.rfind("--v-floor=", 0) == 0) {
      const std::string s = arg == "--v-floor" ? std::string("min") : arg.substr(10);
      if (s == "min") p2_v_floor = 1;
      else if (s == "0") p2_v_floor = 2;
      else {
        std::fprintf(stderr, "--v-floor must be min or 0\n");
        return 2;
      }
    }
    else if (arg == "--v-cut-fb") {
      if (a + 1 < argc && (std::string(argv[a + 1]) == "off" || std::string(argv[a + 1]) == "on"))
        p2_v_cut_fb = std::string(argv[++a]) == "on";
      else
        p2_v_cut_fb = true;
    }
    else if (arg == "--mu-dump") {
      g_mu_dump = next();
    }
    else if (arg == "--g-eq") {
      g_eq_path = next();
    }
    else if (arg == "--miss-sample") {
      g_miss_cap = std::max(0, std::atoi(next()));
      g_miss_path = next();
    }
    else if (arg == "--v-query") {
      const std::string in_path = next();
      g_vquery_out = next();
      if (a + 1 < argc && std::string(argv[a + 1]).rfind("--", 0) != 0) g_vtrace_path = next();
      std::ifstream vin(in_path);
      std::string line;
      while (std::getline(vin, line)) {
        if (line.empty() || line[0] == 's') continue;
        std::istringstream ls(line);
        std::string tag;
        int st, x, y, b;
        if (!(ls >> tag >> st >> x >> y >> b)) continue;
        int kind = tag == "K" ? 1 : tag == "H" ? 2 : 0;
        g_vquery.push_back({kind, st, x, y, b});
      }
    }
    else if (arg == "--lazy-ml-cap") {
      if (a + 1 < argc && (std::string(argv[a + 1]) == "off" || std::string(argv[a + 1]) == "on"))
        lazy_ml_cap = std::string(argv[++a]) == "on";
      else
        lazy_ml_cap = true;
    }
    else if (arg == "--nest-beam" || arg.rfind("--nest-beam=", 0) == 0) {
      const std::string s = arg == "--nest-beam" ? next() : arg.substr(12);
      nest_beam = std::max(0, std::atoi(s.c_str()));
      b1_set = true;
    }
    else if (arg == "--b1") {
      nest_beam = std::max(0, std::atoi(next()));
      b1_set = true;
    }
    else if (arg == "--b2") {
      const int nbeam = std::max(0, std::atoi(next()));
      table_beam[BT_QB] = table_beam[BT_QM] = table_beam[BT_QM2] = table_beam[BT_QMP] =
          table_beam[BT_QP] = nbeam;
      b2_at = a;
    }
    else if (arg == "--qzbar-old") {
      qzbar_old = true;
    }
    else if (arg == "--pbeam-ratio") {
      int ra, rb, rc;
      if (std::sscanf(next(), "%d:%d:%d", &ra, &rb, &rc) != 3 || ra < 0 || rb < 0 || rc < 0) {
        std::fprintf(stderr, "--pbeam-ratio needs a:b:c of non-negative integers\n");
        return 2;
      }
      pbr[0] = ra;
      pbr[1] = rb;
      pbr[2] = rc;
      pbeam_ratio_at = a;
    }
    else if (arg == "--pre-beam" || arg.rfind("--pre-beam=", 0) == 0) {
      const std::string s = arg == "--pre-beam" ? next() : arg.substr(11);
      pre_beam = std::max(0, std::atoi(s.c_str()));
    }
    else if (arg == "--pre-key" || arg.rfind("--pre-key=", 0) == 0) {
      const std::string s = arg == "--pre-key" ? next() : arg.substr(10);
      if (s == "prefix") pre_key = 0;
      else if (s == "inside") pre_key = 1;
      else {
        std::fprintf(stderr, "--pre-key must be prefix or inside\n");
        return 2;
      }
    }
    else if (arg == "--dump-qb") dump_qb_on = true;
    else if (arg == "-h" || arg == "--help") {
      usage();
      return 0;
    } else {
      usage();
      return 2;
    }
  }

