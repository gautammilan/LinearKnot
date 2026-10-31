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
      dp.l2r_beam = pk_l2r == 1 ? beam : 0;
      if (l2r_spawn_cube && dp.l2r_beam > 0) {
        std::fprintf(stderr, "--l2r-spawn cube is not available\n");
        return 2;
      }
      if (l2r_spawn_inc && pk_l2r != 1) {
        std::fprintf(stderr, "--l2r-spawn inc needs --pk-l2r on\n");
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
      dp.lp_hk = l2r_p_hk;
      dp.lp_hex = l2r_p3_hex;
      dp.lp_lin = l2r_p3_lin;
      for (int k = 0; k < 3; ++k) dp.l2r_tb[DP03::PT_P1 + k] = dp.l2r_beam > 0 ? l2r_pbeam[k] : 0;
      if (dp.l2r_beam > 0 && !stem_pr) {
        std::fprintf(stderr, "--pk-l2r on with a beam needs --stem-est-principled on\n");
        return 2;
      }
      dp.grammar_simple = grammar_simple;
      if (grammar_simple && !l2r_spawn_inc) {
        std::fprintf(stderr, "--grammar simple needs --l2r-spawn inc\n");
        return 2;
      }
      dp.l2r_init();
      dp.l2r_init_scores();
    }
    double logz = NEG_INF;
    bool partial = false;
    const double cpu_p2 = cpu_sec(), hwm_p2 = hwm_mb();
    const auto t_scan = std::chrono::steady_clock::now();
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
    if (pk_l2r == 1 && dp.l2r_beam > 0) {
      const double fin = cpu_sec() - cpu_p2;
      std::fprintf(stderr, "[l2r-time] %s before_final_cpu %.3f final_cpu %.3f hwm_before_mb %.0f "
                   "hwm_mb %.0f | final %s\n", sid.c_str(), cpu_p2, fin, hwm_p2, hwm_mb(),
                   dp.l2r_prof(fin + dp.l2r_isec).c_str());
    }
    if (b2_width > 0 || p2_v_cap > 0 || p2_sweep || p2_vsweep) {
      std::fprintf(stderr,
                   "[pass2] %s qz %lld qz_mb %.3f qz_s %.4f v_pairs %lld v_ent %lld v_s %.4f "
                   "g_L %lld prep_s %.4f\n",
                   sid.c_str(), stem.p2_qz_n,
                   static_cast<double>(stem.p2_qz_n * static_cast<long long>(sizeof(StemEst::QzEnt))) /
                       1048576.0,
                   stem.p2_qz_sec, dp.v_pairs, dp.pr_V_entries, dp.pr_static_sec, dp.ntr_kept[0],
                   stem.pr_prep_sec);
    }
