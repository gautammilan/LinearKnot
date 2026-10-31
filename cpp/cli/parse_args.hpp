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
    else if (arg == "--ml-lead-cap") {
      std::string v = next();
      if (v == "n" || v == "off") ml_lead_cap = 0;
      else ml_lead_cap = std::atoi(v.c_str());
    }
    else if (arg == "--max-loop") max_loop = std::atoi(next());
    else if (arg == "--wc-ends") wc_ends = true;
    else if (arg == "--rt") RT = std::atof(next());
    else if (arg == "--cutoff") cutoff = std::atof(next());
    else if (arg == "--no-bpp") do_bpp = false;
    else if (arg == "--gap-beam") {
      std::string m = next();
      if (m == "group") gap_step = false;
      else if (m == "step") gap_step = true;
      else {
        std::fprintf(stderr, "--gap-beam must be 'group' or 'step'\n");
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
    else if (arg == "--stem-est-helpers") {
      std::string s = next();
      helpers_set = true;
      stem_helpers = s == "grow" ? 1 : s == "nested" ? 2 : 0;
      if (s != "off" && s != "grow" && s != "nested") {
        std::fprintf(stderr, "--stem-est-helpers must be 'off', 'grow' or 'nested'\n");
        return 2;
      }
    }
    else if (arg == "--stem-est-principled") {
      std::string s = next();
      stem_pr = s == "on";
      if (s != "off" && s != "on") {
        std::fprintf(stderr, "--stem-est-principled must be 'off' or 'on'\n");
        return 2;
      }
    }
    else if (arg == "--pr-cut") {
      pr_set = true;
      pr_cut = std::strtod(next(), nullptr);
    }
    else if (arg == "--pr-ml-span") {
      pr_set = true;
      pr_ml_span = std::atoi(next());
    }
    else if (arg == "--pr-partner-floor") {
      pr_set = true;
      pr_partner_floor = std::strtod(next(), nullptr);
      if (!(pr_partner_floor >= 0.0 && pr_partner_floor < 1.0)) {
        std::fprintf(stderr, "--pr-partner-floor must be in [0, 1)\n");
        return 2;
      }
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
    else if (arg == "--grow-ml-c") {
      g_ml_c = std::atoi(next());
    }
    else if (arg == "--grow-ml-ub") {
      g_ml_ub = std::strcmp(next(), "on") == 0;
    }
    else if (arg == "--grow-trunc") {
      g_tr_K = std::atoi(next());
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
    else if (arg == "--l2r-spawn") {
      const std::string v = next();
      if (v != "all" && v != "cube" && v != "cube-static" && v != "inc") {
        std::fprintf(stderr, "--l2r-spawn must be all, cube, cube-static or inc\n");
        return 2;
      }
      l2r_spawn_cube = v == "cube" || v == "cube-static";
      l2r_spawn_inc = v == "inc";
    }
    else if (arg == "--l2r-p3-hex") {
      l2r_p3_hex = std::string(next()) == "on";
    }
    else if (arg == "--l2r-p3-lin") {
      l2r_p3_lin = std::string(next()) == "on";
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
    else if (arg == "--pr-armq-sep") {
      pr_set = true;
      pr_aq_sep = std::string(next()) == "on";
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
    else if (arg == "--bV") {
      pr_set = true;
      bV_set = true;
      const int n = std::max(0, std::atoi(next()));
      pr_vcap = n;
      p2_v_cap = n;  // one flag: xs-keyed cap and C_pre left-end cap
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
    else if (arg == "--bG") {
      g_ntr_K = std::atoi(next());
    }
    else if (arg == "--g-sweep") {
      const std::string s = next();
      if (s == "direct") g_sweep = G_SWEEP_DIRECT;
      else if (s == "states") g_sweep = G_SWEEP_STATES;
      else if (s == "lazy") g_sweep = G_SWEEP_LAZY;
      else {
        std::fprintf(stderr, "--g-sweep must be direct, states or lazy\n");
        return 2;
      }
    }
    else if (arg == "--bG-node") {
      g_lz_K = std::atoi(next());
    }
    else if (arg == "--bG-read") {
      g_lz_Kread = std::atoi(next());
    }
    else if (arg == "--bG-node-key") {
      const std::string s = next();
      if (s == "sbar") g_lz_key = 0;
      else if (s == "ul") g_lz_key = 1;
      else {
        std::fprintf(stderr, "--bG-node-key must be sbar or ul\n");
        return 2;
      }
    }
    else if (arg == "--out-id") g_out_id = true;
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
    else if (arg == "-h" || arg == "--help") {
      usage();
      return 0;
    } else {
      usage();
      return 2;
    }
  }

