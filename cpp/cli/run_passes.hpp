// main(), per sequence, part 2: pass-2 setup on the pre-pass and the scan (stamp_p2), then
// pass 3 (pfunc) with timing.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

    if (pk_l2r != 0) {
      if (pk_l2r == 2 && beam != 0) {
        std::fprintf(stderr, "--pk-l2r check needs --beam inf\n");
        return 2;
      }
      dp.l2r_mode = pk_l2r;
      dp.l2r_check_R = pk_l2r_R;
      dp.l2r_beam = pk_l2r == 1 ? (l2r_beam < 0 ? beam : l2r_beam) : 0;
      dp.l2r_spawn_k = dp.l2r_beam > 0 ? l2r_spawn_k : 0;
      dp.sp_cube = l2r_spawn_cube && dp.l2r_beam > 0;
      dp.l2r_arm_fb = l2r_arm_fb && dp.l2r_beam > 0;
      if (dp.sp_cube && l2r_spawn_k <= 0) {
        std::fprintf(stderr, "--l2r-spawn cube needs --l2r-spawn-k k > 0\n");
        return 2;
      }
      dp.sp_pool = l2r_spawn_pool;
      dp.sp_wk = l2r_spawn_wk;
      dp.sp_static = l2r_spawn_static;
      dp.sp_maxpop = l2r_spawn_maxpop;
      if (l2r_spawn_inc && pk_l2r != 1) {
        std::fprintf(stderr, "--l2r-spawn inc needs --pk-l2r on\n");
        return 2;
      }
      if (l2r_inc_check && (!l2r_spawn_inc || dp.l2r_beam > 0)) {
        std::fprintf(stderr, "--l2r-inc-check needs --l2r-spawn inc and no l2r beam\n");
        return 2;
      }
      if (l2r_spawn_inc && g_ntr_K > 0 && !g_ntr_zbar) {
        if (ntr_z_set) {
          std::fprintf(stderr, "--l2r-spawn inc ranks growth nodes from step a2 on: needs --grow-node-trunc-z bar\n");
          return 2;
        }
        g_ntr_zbar = true;
      }
      dp.l2r_inc = l2r_spawn_inc;
      dp.l2r_ninf_mass = l2r_ninf_mass;
      dp.l2r_keep_ninf = l2r_keep_ninf;
      dp.l2r_inc_check = l2r_inc_check;
      dp.qz_audit_req = l2r_qz_audit;
      dp.lp_hk = l2r_p_hk;
      dp.lp_aopen_raw = l2r_p_raw;
      dp.lp_hex = l2r_p3_hex;
      dp.lp_lin = l2r_p3_lin;
      for (int k = 0; k < 3; ++k) dp.l2r_tb[DP03::PT_P1 + k] = dp.l2r_beam > 0 ? l2r_pbeam[k] : 0;
      if (l2r_p3_split[0] > 0 && dp.l2r_beam > 0) {
        dp.lp_split_a = l2r_p3_split[0];
        dp.lp_split_h = l2r_p3_split[1];
        dp.l2r_tb[DP03::PT_P3] = l2r_p3_split[0] * l2r_p3_split[1];
      }
      {
        static const char* pn[DP03::PT_N] = {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us",
                                            "Uls", "Ulss", "Urs", "S", "P1", "P2", "P3"};
        std::stringstream ss(l2r_tbeam);
        std::string item;
        while (std::getline(ss, item, ',')) {
          const size_t eq = item.find('=');
          int t = 0;
          while (t < DP03::PT_N && (eq == std::string::npos || item.substr(0, eq) != pn[t])) ++t;
          if (t == DP03::PT_N) {
            std::fprintf(stderr, "--l2r-tbeam: bad item %s\n", item.c_str());
            return 2;
          }
          const std::string v = item.substr(eq + 1);
          dp.l2r_tb[t] = v == "inf" ? 0 : std::atoi(v.c_str());
        }
      }
      if (dp.l2r_beam > 0 && !(stem_pr && g_grow_rec)) {
        std::fprintf(stderr, "--pk-l2r on with a beam needs --stem-est-principled on and --grow-rec on\n");
        return 2;
      }
      dp.l2r_diag = l2r_diag && dp.l2r_beam > 0;
      dp.rlt = rlt;
      dp.next_pair = next_pair;
      dp.rlt_beam = rlt_beam;
      dp.grammar_simple = grammar_simple;
      dp.c_trace = c_trace;
      if (grammar_simple && !l2r_spawn_inc) {
        std::fprintf(stderr, "--grammar simple needs --l2r-spawn inc\n");
        return 2;
      }
      dp.l2r_init();
      if (g_mem_log) {
        std::vector<MemItem> it;
        dp.mem_items(it, "", true, true);
        stem_mem_items(stem, it);
        mem_log("l2r_init", it);
      }
      // --l2r-p1 R: R rounds of a beamed pseudoknot pass (a copy of this run, beams capped at
      // --l2r-p1-beam) scored with the current Sbar; its stored Qp items with the pre-pass
      // nested Qb (--l2r-p1-qb nested) or its own Qb (pass1) give Sbar_pk, which replaces Sbar
      // for the next round and the final pass.
      // --l2r-ctx p1: every pass-1 round also gives lambda(a') from its Qp items' inside and
      // outside (applied to the next round's scores and the final pass).
      if (l2r_ctx_p1 && dp.l2r_beam > 0) l2r_p1 = std::max(l2r_p1, l2r_ctx_rounds);
      if (l2r_p1 > 0 && dp.l2r_beam > 0) {
        if (!dp.l2r_stem_is(&stem)) {
          std::fprintf(stderr, "--l2r-p1: the run's stem estimator is not main's\n");
          return 2;
        }
        const std::vector<double> S0 = stem.prS;
        for (int r = 1; r <= l2r_p1; ++r) {
          const double c0 = cpu_sec();
          std::vector<double> sn, sp;
          double z1 = NEG_INF, c1 = c0;
          long long nqp = 0, nqb = 0;
          std::string prof;
          std::vector<double> lam;
          std::vector<int> lcnt;
          double c3 = 0.0;
          {
            DP03 p1(dp);
            p1.l2r_cap_beam(l2r_p1_beam);
            p1.l2r_diag = false;
            p1.drops = nullptr;
            p1.mem_scan = false;
            p1.l2r_yw_rebuild();
            p1.l2r_init_scores();
            z1 = p1.pfunc();
            c1 = cpu_sec();
            prof = p1.l2r_prof(c1 - c0);
            dp.l2r_sfx_pass(p1, sn, sp, nqp, nqb);
            if (l2r_ctx_p1) {
              const double c4 = cpu_sec();
              dp.l2r_ctx_lambda(p1, l2r_p1_own ? sp : sn, l2r_ctx_tau, lam, lcnt);
              c3 = cpu_sec() - c4;
            }
            if (g_mem_log) {
              std::vector<MemItem> it;
              p1.mem_items(it, "p1_", true, true);
              mem_log(("p1_round" + std::to_string(r)).c_str(), it);
            }
          }
          if (l2r_ctx_p1) {
            std::vector<double> lv;
            int hi = 0, lo = 0;
            std::string pos;
            char b2[64];
            for (size_t a = 0; a < lam.size(); ++a) {
              if (lcnt[a] == 0) continue;
              lv.push_back(lam[a]);
              if (lam[a] > std::log(2.0)) ++hi;
              if (lam[a] < std::log(0.5)) ++lo;
              std::snprintf(b2, sizeof b2, " %zu:%.3f", a, lam[a]);
              pos += b2;
            }
            std::sort(lv.begin(), lv.end());
            char b[400];
            std::snprintf(b, sizeof b,
                          "[l2r-ctx] %s round %d tau %.3g a_items %zu lambda_gt2 %d lambda_lt0.5 %d "
                          "loglam_min %.4f loglam_med %.4f loglam_max %.4f ctx_cpu %.3f | a:loglam",
                          sid.c_str(), r, l2r_ctx_tau, lv.size(), hi, lo, lv.empty() ? 0.0 : lv.front(),
                          lv.empty() ? 0.0 : lv[lv.size() / 2], lv.empty() ? 0.0 : lv.back(), c3);
            std::fprintf(stderr, "%s%s\n", b, pos.c_str());
            dp.lLam = std::move(lam);
          }
          const double c2 = cpu_sec();
          std::vector<double>& spk = l2r_p1_own ? sp : sn;
          std::string ln = "[l2r-p1] " + sid + " round " + std::to_string(r) + " beam " +
                           std::to_string(l2r_p1_beam) + " qb " + (l2r_p1_own ? "pass1" : "nested");
          char b[400];
          std::snprintf(b, sizeof b,
                        " logZ_pass1 %.9f logS_pk1_nestedqb %.9f logS_pk1_pass1qb %.9f "
                        "logS_nested1 %.9f qp_items %lld qb_items %lld p1_cpu %.3f sfx_cpu %.3f "
                        "hwm_mb %.0f | %s | y:log(Sbar_pk[y+1]/Sbar_nested[y+1])",
                        z1, sn[0], sp[0], S0[0], nqp, nqb, c1 - c0, c2 - c1, hwm_mb(), prof.c_str());
          ln += b;
          for (size_t y = 0; y + 1 < spk.size(); ++y) {
            std::snprintf(b, sizeof b, " %zu:%.2f", y, spk[y + 1] - S0[y + 1]);
            ln += b;
          }
          std::fprintf(stderr, "%s\n", ln.c_str());
          stem.prS = std::move(spk);
        }
        dp.l2r_yw_rebuild();
      }
      dp.lpf = l2r_perf;
      dp.k_memo = k_memo || k_exact_ub;
      dp.k_bound = k_bound;
      dp.k_exact_ub = k_exact_ub;
      dp.k_exact_check = k_exact_check;
      if (k_exact) dp.sx_k_exact = true;
      dp.rlt = rlt;
      dp.next_pair = next_pair;
      dp.rlt_beam = rlt_beam;
      dp.grammar_simple = grammar_simple;
      dp.c_trace = c_trace;
      dp.l2r_init_scores();
      if (p2_g_eager) dp.growth_eager();
      if (g_mem_log) {
        std::vector<MemItem> it;
        dp.mem_items(it, "", true, true);
        stem_mem_items(stem, it);
        mem_log("init_scores", it);
        dp.mem_scan = true;
      }
      dp.l2r_split = l2r_split && dp.l2r_diag;
      dp.lsp_rec = dp.l2r_split;
      dp.l2r_ctx = l2r_ctx && dp.l2r_beam > 0;
      dp.l2r_ctx_init();
      dp.lsub_a = l2r_kx || l2r_sub.find('a') != std::string::npos;
      dp.lsub_c = l2r_sub.find('c') != std::string::npos;
      dp.lsub_d = l2r_sub.find('d') != std::string::npos;
      dp.lsub_e = l2r_sub.find('e') != std::string::npos;
      if ((dp.lsub_c || dp.lsub_d || dp.lsub_e) && !dp.l2r_oracle_load(l2r_oload)) {
        std::fprintf(stderr, "--l2r-sub c/d/e needs a readable --l2r-oracle for this sequence\n");
        return 1;
      }
      if (!l2r_trace.empty()) {
        if (std::sscanf(l2r_trace_key.c_str(), "%d,%d,%d,%d,%d", &dp.ltr_ap, &dp.ltr_xs, &dp.ltr_xp, &dp.ltr_a2,
                        &dp.ltr_c1) < 3 ||
            !(dp.ltr_f = std::fopen(l2r_trace.c_str(), "w"))) {
          std::fprintf(stderr, "--l2r-trace needs a writable file and --l2r-trace-key a',xs,x'\n");
          return 1;
        }
      }
    }
    double logz = NEG_INF;
    bool partial = false;
    const double cpu_p2 = cpu_sec(), hwm_p2 = hwm_mb();
    if (!sample_prof.empty()) sample_start(sample_prof.c_str());
    const auto t_scan = std::chrono::steady_clock::now();
    g_beam_stats.set_pass("scan");
    apply_max_loop(loop_cap_scan);
    dp.set_loop_cap(loop_cap_scan);
    try {
      logz = dp.pfunc();
    } catch (const CpuStop&) {
      partial = true;
    }
    scan_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_scan).count();
    {
      const double cpu_scan = ru_cpu() - cpu0 - cpu_pass1 - cpu_pass2;
      const double hwm_scan = hwm_mb();
      const int nseq = static_cast<int>(seq.size());
      const double ntd = nseq > 0 ? static_cast<double>(nseq) : 1.0;
      const double idx_mb =
          (12.0 * static_cast<double>(stem.entries + stem.entries2 + stem.entries_av) +
           16.0 * static_cast<double>(stem.index.nodes() + stem.index2.nodes() +
                                      stem.index_av.nodes())) /
          1048576.0;
      std::fprintf(stderr,
                   "[stage] scan n %d cpu %.3f cpu_s_nt %.6f hwm_mb %.1f hwm_mb_nt %.4f\n",
                   nseq, cpu_scan, cpu_scan / ntd, hwm_scan, hwm_scan / ntd);
      std::fprintf(stderr,
                   "[tables] n %d stem_idx %zu stem_idx_mb %.1f pr_Mbar %zu pr_Mbar_mb %.1f "
                   "pr_QbIn %zu pr_QbIn_mb %.1f n2_mb %.1f\n",
                   nseq, stem.entries, idx_mb, stem.prMbar.size(),
                   static_cast<double>(stem.prMbar.size()) * 8.0 / 1048576.0, stem.prQbIn.size(),
                   static_cast<double>(stem.prQbIn.size()) * 8.0 / 1048576.0,
                   8.0 * ntd * ntd / 1048576.0);
    }
    if (!sample_prof.empty()) sample_stop();
    if (g_mem_log) {
      std::vector<MemItem> it;
      dp.mem_items(it, "", true, pk_l2r != 0);
      mem_log("scan_done", it);
    }
    if (pk_l2r == 1 && dp.l2r_beam > 0) {
      const double fin = cpu_sec() - cpu_p2;
      std::fprintf(stderr, "[l2r-time] %s before_final_cpu %.3f final_cpu %.3f hwm_before_mb %.0f "
                   "hwm_mb %.0f | final %s\n", sid.c_str(), cpu_p2, fin, hwm_p2, hwm_mb(),
                   dp.l2r_prof(fin + dp.l2r_isec).c_str());
      dp.profile_dump(sid.c_str(), fin, stem.prof_qz_pre);
    }
    if (g_mass_f) g_mass_flush();
    if (b2_width > 0 || p2_qz_beam > 0 || p2_v_mode > 0 || p2_g_mode > 0 || p2_v_cap > 0 ||
        p2_g_suffix || p2_g_eager || p2_g_fb || p2_sweep || p2_vsweep) {
      std::fprintf(stderr,
                   "[pass2] %s qz %lld qz_mb %.3f qz_s %.4f v_pairs %lld v_ent %lld v_s %.4f "
                   "g_L %lld prep_s %.4f\n",
                   sid.c_str(), stem.p2_qz_n,
                   static_cast<double>(stem.p2_qz_n * static_cast<long long>(sizeof(StemEst::QzEnt))) /
                       1048576.0,
                   stem.p2_qz_sec, dp.v_pairs, dp.pr_V_entries, dp.pr_static_sec, dp.ntr_kept[0],
                   stem.pr_prep_sec);
    }
