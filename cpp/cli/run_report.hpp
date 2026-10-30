// main(), per sequence, part 3: pass-2 / scan reports, checks, bpp, dumps, and the output
// record (>id, logZ, sec, pairs).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

    if (pk_l2r != 0) {
      const std::vector<long long> sc = dp.l2r_states();
      std::fprintf(stderr,
                   "[l2r] %s n %d mode %s states R %lld Rs %lld Rls %lld Rlss %lld Rrs %lld U %lld "
                   "Us %lld Uls %lld Ulss %lld Urs %lld S %lld tr R %lld Rs %lld Rls %lld Rrs %lld "
                   "U %lld Us %lld Uls %lld Urs %lld S %lld check_s %.3f",
                   sid.c_str(), static_cast<int>(seq.size()), pk_l2r == 1 ? "on" : "check", sc[0],
                   sc[1], sc[2], sc[3], sc[4], sc[5], sc[6], sc[7], sc[8], sc[9], sc[10],
                   dp.l2r_tr[DP03::LT_AG], dp.l2r_tr[DP03::LT_AGS], dp.l2r_tr[DP03::LT_ALS],
                   dp.l2r_tr[DP03::LT_ARS], dp.l2r_tr[DP03::LT_BG], dp.l2r_tr[DP03::LT_BGS],
                   dp.l2r_tr[DP03::LT_BLS], dp.l2r_tr[DP03::LT_BRS], dp.l2r_tr[DP03::LT_S],
                   dp.l2r_sec[5]);
      if (pk_l2r == 2)
        std::fprintf(stderr, " checkR %lld mis %lld max_rel %.3e checkQp %lld mis %lld max_rel %.3e",
                     dp.lc_R_n, dp.lc_R_mis, dp.lc_R_rel, dp.lc_Qp_n, dp.lc_Qp_mis, dp.lc_Qp_rel);
      std::fprintf(stderr, "\n");
      const double nn = static_cast<double>(seq.size());
      auto per = [&](long long x) { return nn > 0.0 ? static_cast<double>(x) / nn : 0.0; };
      const long long r_ml = dp.ext_n[1] + dp.ext_n[2] + dp.ext_n[3];
      const long long u_ml = dp.ext_n[5] + dp.ext_n[6] + dp.ext_n[7];
      std::fprintf(stderr,
                   "[ext] n %d R_int %.3f R_ml %.3f R_f %.3f R_g %.3f R_h %.3f "
                   "U_int %.3f U_ml %.3f U_f %.3f U_g %.3f U_h %.3f\n",
                   static_cast<int>(seq.size()), per(dp.ext_n[0]), per(r_ml), per(dp.ext_n[1]),
                   per(dp.ext_n[2]), per(dp.ext_n[3]), per(dp.ext_n[4]), per(u_ml), per(dp.ext_n[5]),
                   per(dp.ext_n[6]), per(dp.ext_n[7]));
      std::fprintf(stderr,
                   "[it] n %d R_int %.3f R_h %.3f R_ls %.3f R_rs %.3f "
                   "U_int %.3f U_h %.3f U_ls %.3f U_rs %.3f\n",
                   static_cast<int>(seq.size()), per(dp.it_int[0]), per(dp.it_h[0]),
                   per(dp.it_ls[0]), per(dp.it_rs[0]), per(dp.it_int[1]), per(dp.it_h[1]),
                   per(dp.it_ls[1]), per(dp.it_rs[1]));
      std::fprintf(stderr,
                   "[qp] n %d beam %d pre_per %.3f kept_per %.3f pre_max %d kept_max %d bind %lld\n",
                   static_cast<int>(seq.size()), dp.qp_beam, per(dp.qp_pre_n),
                   per(dp.qp_kept_n), dp.qp_pre_max, dp.qp_kept_max, dp.qp_bind_n);
      std::fprintf(stderr,
                   "[qm] n %d beam %d pre_per %.3f kept_per %.3f pre_max %d kept_max %d bind %lld\n",
                   static_cast<int>(seq.size()), dp.qm_beam, per(dp.qm_pre_n), per(dp.qm_kept_n),
                   dp.qm_pre_max, dp.qm_kept_max, dp.qm_bind_n);
      long long scored = 0;
      for (int t = DP03::PT_R; t < DP03::PT_N; ++t) scored += dp.l2r_pre[t];
      const double wall_s =
          std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      const double ntd = nn > 0.0 ? nn : 1.0;
      std::fprintf(stderr,
                   "[pos] n %d per_nt interior %.3f h %.3f ls %.3f rs %.3f tail %.3f "
                   "S %.3f P1 %.3f P2 %.3f P3 %.3f Qz %.3f Qz_pre %.3f K %.3f scored %.3f\n",
                   static_cast<int>(seq.size()),
                   per(dp.it_int[0] + dp.it_int[1]), per(dp.it_h[0] + dp.it_h[1]),
                   per(dp.it_ls[0] + dp.it_ls[1]), per(dp.it_rs[0] + dp.it_rs[1]), per(dp.rlt_ext),
                   per(dp.l2r_tr[DP03::LT_S]), per(dp.l2r_pre[DP03::PT_P1]),
                   per(dp.l2r_pre[DP03::PT_P2]), per(dp.l2r_pre[DP03::PT_P3]), per(dp.qz_cand),
                   per(stem.p2_qz_n), per(dp.sx_k_n), per(scored));
      std::fprintf(stderr,
                   "[vmiss] n %d floor %d cut_fb %d look %lld miss %lld look_U %lld miss_U %lld "
                   "look_R %lld miss_R %lld look_S %lld miss_S %lld key0_U %lld part_U %lld "
                   "key0_R %lld fb %lld\n",
                   static_cast<int>(seq.size()), dp.p2_v_floor, dp.p2_v_cut_fb ? 1 : 0, dp.rfv_look,
                   dp.rfv_miss, dp.rfv_look_u, dp.rfv_miss_u, dp.rfv_look_r, dp.rfv_miss_r,
                   dp.rfv_look_s, dp.rfv_miss_s, dp.rfv_key0_u, dp.rfv_part_u, dp.rfv_key0_r,
                   dp.vfb_n);
      {
        static const char* sn[DP03::ST_N] = {"P1", "P2", "P3", "R", "S", "U"};
        for (int s = 0; s < DP03::ST_N; ++s)
          std::fprintf(stderr,
                       "[vmiss3] %s V %lld %lld %lld K %lld %lld %lld H %lld %lld %lld\n",
                       sn[s], dp.vm_i[s], dp.vm_ii[s], dp.vm_iii[s], dp.km_i[s], dp.km_ii[s],
                       dp.km_iii[s], dp.hm_i[s], dp.hm_ii[s], dp.hm_iii[s]);
      }
      if (g_miss_cap > 0 && !g_miss_path.empty()) {
        std::ofstream mf(g_miss_path);
        mf << "seen V " << g_miss_seen[0] << " K " << g_miss_seen[1] << " H " << g_miss_seen[2]
           << "\n";
        static const char kn[3] = {'V', 'K', 'H'};
        for (int k = 0; k < 3; ++k)
          for (const MissRec& r : g_miss[k])
            mf << kn[k] << " " << r.st << " " << r.x << " " << r.y << " " << r.b << "\n";
      }
      std::fprintf(stderr,
                   "[times] n %d pass1 %.3f pass2 %.3f scan %.3f wall %.3f "
                   "pass1_s_nt %.6f pass2_s_nt %.6f scan_s_nt %.6f wall_s_nt %.6f\n",
                   static_cast<int>(seq.size()), pass1_s, pass2_s, scan_s - dp.rfv_sec, wall_s,
                   pass1_s / ntd, pass2_s / ntd, (scan_s - dp.rfv_sec) / ntd, wall_s / ntd);
      std::fprintf(stderr,
                   "[vfb] sec %.4f calls %lld  Qz cand_per %.3f kept_per %.3f reads_per %.3f\n",
                   dp.rfv_sec, dp.rfv_od, per(dp.qz_cand), per(dp.qz_kept), per(dp.qz_reads));
      if (!dp.col_cnt.empty()) {
        const int nc = static_cast<int>(dp.col_cnt.size());
        const int a = (nc * 4) / 5;
        auto tail = [&](auto get) {
          long long prev = a > 0 ? get(dp.col_cnt[static_cast<size_t>(a - 1)]) : 0;
          long long s = 0;
          for (int j = a; j < nc; ++j) {
            const long long v = get(dp.col_cnt[static_cast<size_t>(j)]);
            s += v - prev;
            prev = v;
          }
          return nc > a ? static_cast<double>(s) / (nc - a) : 0.0;
        };
        std::fprintf(stderr,
                     "[pos20] n %d from %d int %.3f h %.3f ls %.3f rs %.3f tail %.3f "
                     "edge %.3f qz_reads %.3f K %.3f qm_kept %.3f\n",
                     nc, a, tail([](const DP03::ColCnt& c) { return c.iti; }),
                     tail([](const DP03::ColCnt& c) { return c.ith; }),
                     tail([](const DP03::ColCnt& c) { return c.itls; }),
                     tail([](const DP03::ColCnt& c) { return c.itrs; }),
                     tail([](const DP03::ColCnt& c) { return c.tail; }),
                     tail([](const DP03::ColCnt& c) { return c.edge; }),
                     tail([](const DP03::ColCnt& c) { return c.qz; }),
                     tail([](const DP03::ColCnt& c) { return c.k; }),
                     tail([](const DP03::ColCnt& c) { return c.qm; }));
      }
      {
        std::vector<MemItem> mit;
        dp.mem_items(mit, "", true, pk_l2r != 0);
        stem_mem_items(stem, mit);
        double tot = 0;
        for (const MemItem& m : mit) tot += m.bytes;
        std::fprintf(stderr, "[mem-nt] n %d total_mb %.3f total_mb_nt %.6f",
                     static_cast<int>(seq.size()), tot / 1048576.0, tot / 1048576.0 / ntd);
        for (const MemItem& m : mit) {
          if (m.bytes < 1024.0 && m.ents == 0) continue;
          std::fprintf(stderr, " %s %lld %.3f %.6f", m.name.c_str(), m.ents, m.bytes / 1048576.0,
                       m.bytes / 1048576.0 / ntd);
        }
        std::fprintf(stderr, "\n");
      }
      dp.log_nest("scan");
      if (dp.rlt) {
        std::fprintf(stderr, "[rlt] n %d beam %d kept_per %.3f max %d\n", static_cast<int>(seq.size()),
                     dp.rlt_beam, per(dp.rlt_kept), dp.rlt_max);
        std::fprintf(stderr, "[rlt-fam] R %lld max %d U %lld max %d\n", dp.rlt_kept_f[0], dp.rlt_max_f[0],
                     dp.rlt_kept_f[1], dp.rlt_max_f[1]);
        const int rn = static_cast<int>(dp.rlt_pos[0].size());
        const int un = static_cast<int>(dp.rlt_pos[1].size());
        const int pn = std::max(rn, un);
        for (int j = 0; j < pn; ++j) {
          const int rv = j < rn ? dp.rlt_pos[0][j] : 0;
          const int uv = j < un ? dp.rlt_pos[1][j] : 0;
          if (rv || uv)
            std::fprintf(stderr, "[rlt-pos] %d %d %d\n", j, rv, uv);
        }
      }
      if (dp.l2r_inc) {
        long long st[3] = {};
        for (int k = 0; k < 3; ++k)
          for (const auto& m : k == 0 ? dp.LP1 : k == 1 ? dp.LP2 : dp.LP3)
            st[k] += static_cast<long long>(m.keys.size());
        std::string b;
        char buf[256];
        for (int k = 0; k < 3; ++k) {
          const int t = DP03::PT_P1 + k;
          std::snprintf(buf, sizeof buf,
                        " P%d beam %d pre %lld kept %lld max %lld stored %lld tr %lld build_s %.3f "
                        "score_s %.3f miss_arm %lld miss_H %lld miss_K %lld miss_open %lld",
                        k + 1, dp.ltb(t), dp.l2r_pre[t], dp.l2r_kept[t], dp.l2r_kmax[t], st[k],
                        dp.l2r_tr[DP03::LT_P1 + k], dp.l2r_bsec[t], dp.l2r_psec[t], dp.lp_miss[k][0],
                        dp.lp_miss[k][1], dp.lp_miss[k][2], dp.lp_miss[k][3]);
          b += buf;
        }
        std::fprintf(stderr,
                     "[l2r-inc] %s n %d hk %d aopen %s%s seed_checks %lld seeds %lld Rs_pre %lld "
                     "Rs_kept %lld seed_s %.3f H_lists %lld H_roots %lld H_ents %lld H_s %.3f "
                     "check %lld mis %lld max_rel %.3e p3_split %d,%d a2_groups %lld a2_kept %lld "
                     "p3_hex %d hex_builds %lld hex_terms %lld hex_ents %lld hex_s %.3f p3_lin %d\n",
                     sid.c_str(), static_cast<int>(seq.size()), dp.lp_hk,
                     dp.lp_aopen_raw ? "raw" : "q", b.c_str(), dp.lp_seed_checks, dp.lp_seeds,
                     dp.l2r_pre[DP03::PT_RS], dp.l2r_kept[DP03::PT_RS], dp.l2r_bsec[DP03::PT_RS],
                     dp.lp_h_lists, dp.lp_h_roots, dp.lp_h_ents, dp.lp_h_sec, dp.lc_I_n, dp.lc_I_mis,
                     dp.lc_I_rel, dp.lp_split_a, dp.lp_split_h, dp.lp_split_grp, dp.lp_split_grp_kept,
                     dp.lp_hex ? 1 : 0, dp.lp_hx_builds, dp.lp_hx_terms, dp.lp_hx_ents, dp.lp_hx_sec,
                     dp.lp_lin_ok() ? 1 : 0);
      }
      if (dp.lpf) {
        std::string b;
        char buf[160];
        auto I = [&](const char* k, long long v) {
          std::snprintf(buf, sizeof buf, " %s %lld", k, v);
          b += buf;
        };
        auto D = [&](const char* k, double v) {
          std::snprintf(buf, sizeof buf, " %s %.4f", k, v);
          b += buf;
        };
        static const char* pn[DP03::PT_N] = {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us",
                                            "Uls", "Ulss", "Urs", "S", "P1", "P2", "P3"};
        std::string k;
        for (int t = 0; t < DP03::PT_N; ++t) {
          if (dp.lpf_npos[t] == 0) continue;
          k = std::string("a_") + pn[t] + "_";
          I((k + "beam").c_str(), dp.ltb(t));
          I((k + "pos").c_str(), dp.lpf_npos[t]);
          I((k + "pre").c_str(), dp.l2r_pre[t]);
          I((k + "premax").c_str(), dp.lpf_premax[t]);
          I((k + "kept").c_str(), dp.l2r_kept[t]);
          I((k + "keptmax").c_str(), dp.l2r_kmax[t]);
          I((k + "full").c_str(), dp.lpf_full[t]);
          I((k + "scored").c_str(), dp.lpf_nsc[t]);
          D((k + "build_s").c_str(), dp.l2r_bsec[t]);
          D((k + "prune_s").c_str(), dp.lpf_prune_sec[t]);
        }
        static const char* fn[3] = {"S", "P1", "P3"};
        for (int f = 0; f < 3; ++f) {
          if (dp.lpf_fx_pos[f] == 0) continue;
          k = std::string("c_") + fn[f] + "_";
          I((k + "pos").c_str(), dp.lpf_fx_pos[f]);
          I((k + "items").c_str(), dp.lpf_fx_items[f]);
          I((k + "pair").c_str(), dp.lpf_fx_pair[f]);
          I((k + "pairmax").c_str(), dp.lpf_fx_pairmax[f]);
          I((k + "unp").c_str(), dp.lpf_fx_unp[f]);
          D((k + "s").c_str(), dp.lpf_fx_sec[f]);
        }
        for (int t = 0; t < 3; ++t) {
          k = std::string("d_P") + std::to_string(t + 1) + "_";
          I((k + "calls").c_str(), dp.lpf_sc_calls[t]);
          I((k + "hsteps").c_str(), dp.lpf_sc_hsteps[t]);
          I((k + "dot").c_str(), dp.lpf_sc_dot[t]);
          I((k + "kcalls").c_str(), dp.lpf_sc_kcalls[t]);
          I((k + "ksteps").c_str(), dp.lpf_sc_ksteps[t]);
          I((k + "osteps").c_str(), dp.lpf_sc_osteps[t]);
          D((k + "h_s").c_str(), dp.lpf_sc_hsec[t]);
          D((k + "k_s").c_str(), dp.lpf_sc_ksec[t]);
          D((k + "o_s").c_str(), dp.lpf_sc_osec[t]);
        }
        I("d_pspan", dp.pr_pspan);
        I("d_aq_sep", dp.pr_aq_sep ? 1 : 0);
        I("e_lists", dp.lp_h_lists);
        I("e_roots", dp.lp_h_roots);
        I("e_newroots", dp.lpf_h_newroots);
        I("e_lam", dp.lpf_h_lam);
        I("e_merge", dp.lpf_h_merge);
        I("e_trunc", dp.lpf_h_trunc);
        I("e_trunc_in", dp.lpf_h_trunc_in);
        I("e_ents", dp.lp_h_ents);
        I("e_helpers", dp.lpf_h_helpers);
        I("e_clears", dp.lpf_h_clears);
        D("e_s", dp.lp_h_sec);
        D("e_root_s", dp.lpf_h_root_sec);
        D("e_help_s", dp.lpf_h_help_sec);
        D("e_merge_s", dp.lpf_h_merge_sec);
        I("f_hg_clears", dp.lpf_hg_clears);
        D("f_hg_s", dp.lpf_hg_sec);
        D("f_mtrim_s", dp.lpf_mtrim_sec);
        I("f_help_built", dp.lpf_help_built);
        I("f_help_distinct", static_cast<long long>(dp.lpf_help_keys.size()));
        I("f_help_rebuilt", dp.lpf_help_rebuilt);
        D("f_help_s", dp.lpf_help_sec);
        I("f_ml_nodes", dp.lpf_ml_nodes);
        I("f_ml_cand", dp.lpf_ml_cand);
        I("f_ml_max", dp.lpf_ml_max);
        I("f_ub_nodes", dp.ml_ub_nodes);
        I("f_ub_cand", dp.ml_ub_cand);
        I("f_ub_built", dp.ml_ub_built);
        I("f_ub_keys", static_cast<long long>(dp.prRHub.size()));
        I("f_racc_dense", dp.lpf_racc_dense);
        I("f_racc_fin_n", dp.lpf_racc_n);
        I("f_racc_fin", dp.lpf_racc_fin);
        I("f_rn_roots", dp.rn_roots);
        I("f_rn_nodes", dp.rn_nodes);
        I("f_stc_n", dp.lpf_stc_n);
        I("f_stc_buckets", dp.lpf_stc_buckets);
        I("f_stc_scan", dp.lpf_stc_scan);
        D("f_stc_s", dp.lpf_stc_sec);
        I("v_pairs", dp.vi_pairs);
        D("v_pairs_s", dp.vi_pairs_s);
        I("v_nodes_out", dp.vi_nodes[0]);
        I("v_nodes_in", dp.vi_nodes[1]);
        I("v_nb_int", dp.vi_nb_int);
        I("v_nb_ml", dp.vi_nb_ml);
        I("v_ml_nodes", dp.vi_ml_nodes);
        I("v_mlc", dp.vi_mlc);
        I("v_mlc_max", dp.vi_mlc_max);
        I("v_help_built", dp.pr_hbuilt[1]);
        I("v_help_terms", dp.pr_hterms[1]);
        I("v_help_ents", dp.pr_hentries[1]);
        I("v_add_int", dp.vi_add_int);
        I("v_add_ml", dp.vi_add_ml);
        I("v_add_help", dp.vi_help_adds);
        I("v_ent_fin", dp.vi_ent_fin);
        I("v_ent_kept", dp.vi_ent_kept);
        I("v_ats", dp.vi_ats);
        D("v_build_s", dp.vi_vbuild_s);
        D("v_steps_s", dp.vi_steps_s);
        D("v_terms_s", dp.vi_terms_s);
        D("v_help_s", dp.vi_help_s);
        D("v_add_int_s", dp.vi_add_int_s);
        D("v_add_ml_s", dp.vi_add_ml_s);
        D("v_fin_s", dp.vi_fin_s);
        D("v_ats_s", dp.vi_ats_s);
        D("v_arm_s", dp.vi_arm_s);
        D("v_rest_s", dp.vi_rest_s);
        D("v_init_s", dp.l2r_isec);
        D("v_prep_s", stem.pr_prep_sec);
        std::fprintf(stderr, "[l2r-perf] %s n %d%s\n", sid.c_str(), static_cast<int>(seq.size()), b.c_str());
      }
      if (dp.qz_audit || g_ntr_zbar)
        std::fprintf(stderr,
                     "[l2r-qz] %s n %d ntr_z %s ntr_first_j %d ntr_arm_skip %lld audit %d qz_reads %lld "
                     "qz_unfilled %lld qz_unfilled_pr %lld\n",
                     sid.c_str(), static_cast<int>(seq.size()), g_ntr_zbar ? "bar" : "scan", dp.ntr_first_j,
                     dp.ntr_arm_skip, static_cast<int>(dp.qz_audit), dp.qz_reads, dp.qz_unf, dp.qz_unf_pr);
      // per table: states before / after the beam (summed over steps), max kept per step,
      // build and score+prune seconds
      static const char* pn[DP03::PT_N] = {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us",
                                          "Uls", "Ulss", "Urs", "S", "P1", "P2", "P3"};
      std::string b;
      char buf[256];
      double bs = 0.0, ps = 0.0;
      for (int t = 0; t < (dp.l2r_inc ? DP03::PT_N : DP03::PT_P1); ++t) {
        std::snprintf(buf, sizeof buf, " %s %lld/%lld/%lld/%.3f/%.3f", pn[t], dp.l2r_pre[t],
                      dp.l2r_kept[t], dp.l2r_kmax[t], dp.l2r_bsec[t], dp.l2r_psec[t]);
        b += buf;
        bs += dp.l2r_bsec[t];
        ps += dp.l2r_psec[t];
      }
      std::fprintf(stderr,
                   "[l2r-beam] %s n %d beam %d spawn_k %d steps %d tables(pre/kept/max/build_s/"
                   "score_s)%s close_s %.3f build_s %.3f score_s %.3f init_s %.3f spawn events "
                   "%lld cand %lld kept %lld arms %d\n",
                   sid.c_str(), static_cast<int>(seq.size()), dp.l2r_beam, dp.l2r_spawn_k,
                   static_cast<int>(seq.size()), b.c_str(), dp.l2r_csec, bs, ps, dp.l2r_isec,
                   dp.sp_events, dp.sp_cand, dp.sp_kept, dp.l2r_beam > 0 ? dp.l2r_arms() : 0);
      dp.scan_ct_dump();
      if (dp.l2r_beam > 0 && !partial)
        std::fprintf(stderr,
                     "[l2r-spawn] %s n %d mode %s w %s k %d pool %d maxpop %d events %lld cand %lld "
                     "pops %lld kept %lld row_ents %zu row_mb %.3f col_ents %zu col_max %zu "
                     "cols_s %.3f rows_s %.3f heap_s %.3f rank_s %.3f spawn_s %.3f rs_mass %.9g\n",
                     sid.c_str(), static_cast<int>(seq.size()),
                     dp.l2r_inc ? "inc" : !dp.sp_cube ? "all" : dp.sp_static ? "cube-static" : "cube",
                     dp.sp_wk ? "k" : "suf",
                     dp.l2r_spawn_k, dp.sp_pool, dp.sp_maxpop > 0 ? dp.sp_maxpop : 4 * dp.sp_pool,
                     dp.sp_events, dp.sp_cand, dp.sp_pops, dp.sp_kept, dp.sp_row_ents,
                     dp.sp_row_ents * sizeof(int) / 1048576.0, dp.sp_col_ents, dp.sp_col_max,
                     dp.sp_col_sec, dp.sp_row_sec,
                     dp.sp_heap_sec, dp.sp_rank_sec, dp.l2r_bsec[DP03::PT_RS] + dp.l2r_psec[DP03::PT_RS],
                     dp.l2r_rs_mass());
      if (dp.l2r_beam > 0 && !partial && dp.l2r_ninf_mass)
        std::fprintf(stderr, "[l2r-ninf] %s n %d keep_ninf %d%s\n",
                     sid.c_str(), static_cast<int>(seq.size()), static_cast<int>(dp.l2r_keep_ninf),
                     dp.l2r_ninf_report().c_str());
      std::fprintf(stderr,
                   "[l2r-v] %s n %d trim %s eps %.3g cap %d v_before %lld v_after %lld v_capped %lld "
                   "arms %d arm_entries %zu arm_groups %zu arm_fb %d fb_hit %lld fb_floor %lld\n",
                   sid.c_str(), static_cast<int>(seq.size()), dp.pr_vsuffix ? "suffix" : "global",
                   dp.pr_veps, dp.pr_vcap, dp.v_before, dp.v_after, dp.v_capped, dp.l2r_arms(),
                   dp.l2r_arm_ents(), dp.l2r_arm_groups(), static_cast<int>(dp.l2r_arm_fb), dp.l2r_fb_hit,
                   dp.l2r_fb_floor);
      if (dp.pr_vmlk > 0)
        std::fprintf(stderr,
                     "[l2r-vk] %s n %d k %d nodes %lld cand %lld drop %lld resid_max %.4g resid_mean %.4g "
                     "mass_s %.3f init_s %.3f\n",
                     sid.c_str(), static_cast<int>(seq.size()), dp.pr_vmlk, dp.vk_nodes, dp.vk_cand,
                     dp.vk_drop, dp.vk_rmax, dp.vk_nodes > 0 ? dp.vk_rsum / dp.vk_nodes : 0.0,
                     dp.vk_mass_s, dp.l2r_isec);
      std::fprintf(stderr,
                   "[l2r-mem] %s n %d v_key %s v_stream %d v_pairs %lld v_entries %lld v_per_pair %.1f "
                   "v_live_max %lld v_sec %.3f arm_entries %zu node_trunc %d ntr_nodes %lld "
                   "ntr_L_in %lld ntr_L_kept %lld ntr_R_in %lld ntr_R_kept %lld ntr_R_K %d ntr_R_miss %lld "
                   "malloc_trim %d\n",
                   sid.c_str(), static_cast<int>(seq.size()), dp.pr_vxs ? "xs" : "xsys",
                   static_cast<int>(dp.pr_vstream), dp.v_pairs, dp.pr_V_entries,
                   dp.v_pairs > 0 ? static_cast<double>(dp.pr_V_entries) / dp.v_pairs : 0.0, dp.v_live_max,
                   dp.pr_static_sec, dp.l2r_arm_ents(), g_ntr_K, dp.ntr_nodes, dp.ntr_in[0], dp.ntr_kept[0],
                   dp.ntr_in[1], dp.ntr_kept[1], g_ntr_KR > 0 ? g_ntr_KR : g_ntr_K, dp.ntr_rmiss,
                   static_cast<int>(g_mtrim));
      if (dp.sw_as_rn_calls > 0 || dp.sw_fill_sec > 0.0) {
        const double mean_in =
            dp.sw_as_rn_calls > 0
                ? static_cast<double>(dp.sw_as_rn_ent_in) / static_cast<double>(dp.sw_as_rn_calls)
                : 0.0;
        const double mean_out =
            dp.sw_as_rn_calls > 0
                ? static_cast<double>(dp.sw_as_rn_ent_out) / static_cast<double>(dp.sw_as_rn_calls)
                : 0.0;
        const double mean_ga =
            dp.sw_ga_nodes > 0 ? static_cast<double>(dp.sw_ga_ent_sum) / static_cast<double>(dp.sw_ga_nodes)
                              : 0.0;
        const double work_per_out =
            mean_out > 0.0 ? mean_in / mean_out : 0.0;
        const char* gcut =
            g_cut == G_CUT_COLLAPSE ? "collapse" : (g_cut == G_CUT_ENTRIES ? "entries" : "topq");
        std::fprintf(stderr,
                     "[g-read] n %d bG %d g_cut %s ntr_init_sec %.6f sw_fill_sec %.6f "
                     "sw_as_rn_calls %lld sw_as_rn_sec %.6f ms_per_read %.4f ent_in_mean %.1f "
                     "ent_out_mean %.1f work_per_out %.1f ga_nodes %lld ga_ent_mean %.1f "
                     "ga_ent_max %lld rn_keys %lld grow_sec %.6f\n",
                     static_cast<int>(seq.size()), g_ntr_K, gcut, dp.sw_ntr_init_sec, dp.sw_fill_sec,
                     dp.sw_as_rn_calls, dp.sw_as_rn_sec,
                     dp.sw_as_rn_calls > 0 ? 1e3 * dp.sw_as_rn_sec / static_cast<double>(dp.sw_as_rn_calls)
                                          : 0.0,
                     mean_in, mean_out, work_per_out, dp.sw_ga_nodes, mean_ga, dp.sw_ga_ent_max,
                     dp.rn_keys, dp.pr_grow_sec);
        auto meanll = [](long long s, long long n) {
          return n > 0 ? static_cast<double>(s) / static_cast<double>(n) : 0.0;
        };
        std::fprintf(stderr,
                     "[vg-fill] V pushes %lld pre_mean %.1f pre_max %lld post_mean %.1f post_max %lld "
                     "merge_sec %.6f cut_sec %.6f | G pushes %lld pre_mean %.1f pre_max %lld "
                     "post_mean %.1f post_max %lld merge_sec %.6f cut_sec %.6f\n",
                     dp.v_pushes, meanll(dp.v_pre_sum, dp.v_pre_n), dp.v_pre_max,
                     meanll(dp.v_post_sum, dp.v_post_n), dp.v_post_max, dp.v_merge_sec, dp.v_cut_sec,
                     dp.g_pushes, meanll(dp.g_pre_sum, dp.g_pre_n), dp.g_pre_max,
                     meanll(dp.g_post_sum, dp.g_post_n), dp.g_post_max, dp.g_merge_sec, dp.g_cut_sec);
        auto yo_mean = [](long long s, long long n) {
          return n > 0 ? static_cast<double>(s) / static_cast<double>(n) : 0.0;
        };
        std::fprintf(stderr,
                     "[ml-yo] qm_ix %d S %d path %s pairs %lld yo_mean %.2f yo_max %lld | "
                     "scan_path pairs %lld yo_mean %.2f yo_max %lld\n",
                     dp.ml_yo_qm_ix ? 1 : 0, dp.ml_yo_S, dp.ml_yo_qm_ix ? "qm_ix" : "scan",
                     dp.ml_yo_qm_n, yo_mean(dp.ml_yo_qm_sum, dp.ml_yo_qm_n), dp.ml_yo_qm_max,
                     dp.ml_yo_sc_n, yo_mean(dp.ml_yo_sc_sum, dp.ml_yo_sc_n), dp.ml_yo_sc_max);
        double vf_med = 0.0;
        if (!dp.arm_Vf_samples.empty()) {
          auto v = dp.arm_Vf_samples;
          std::nth_element(v.begin(), v.begin() + static_cast<long>(v.size() / 2), v.end());
          vf_med = v[v.size() / 2];
        }
        std::fprintf(stderr,
                     "[kl-m] read_sec %.6f kl_in_read_sec %.6f arms %lld Vf_mean %.1f Vf_med %.1f "
                     "Vf_max %lld Kq_mean %.1f Kq_max %lld\n",
                     dp.sw_as_rn_sec, dp.kl_in_read_sec, dp.arm_read_n,
                     meanll(dp.arm_Vf_sum, dp.arm_read_n), vf_med, dp.arm_Vf_max,
                     meanll(dp.arm_Kq_sum, dp.arm_read_n), dp.arm_Kq_max);
      }
      if (g_mem_log) {
        static const char* pn[DP03::PT_N] = {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us",
                                             "Uls", "Ulss", "Urs", "S", "P1", "P2", "P3"};
        std::string tb;
        char buf[160];
        for (int t = 0; t < DP03::PT_N; ++t) {
          std::snprintf(buf, sizeof buf, " %s %lld/%lld/%lld", pn[t], dp.l2r_nsc[t], dp.l2r_ninf[t],
                        dp.l2r_noarm[t]);
          tb += buf;
        }
        long long re = 0, ve = 0;
        const double rb = dp.rside_bytes(re), vb = dp.rfv_bytes(ve);
        auto mean = [](long long s, long long c) { return c > 0 ? static_cast<double>(s) / c : 0.0; };
        std::fprintf(stderr,
                     "[l2r-mem2] %s n %d v_cap_pairs %lld v_cap_entries %lld v_cap_frac_mean %.3g "
                     "v_cap_frac_max %.3g v_cap_vfrac_max %.3g v_cap_vfrac_tot %.3g R_node %.1f/%lld "
                     "R_root %.1f/%lld R_key %.1f/%lld R_entries %lld R_mb %.1f rfv %d V_keep %lld "
                     "V_keep_mb %.1f rfv_look %lld rfv_miss %lld rfv_seed %lld rfv_od %lld rfv_od_ent %lld "
                     "rfv_child_miss %lld rfv_keys %lld rfv_key_ent %lld | scored/ninf/noarm%s\n",
                     sid.c_str(), static_cast<int>(seq.size()), dp.v_cap_pairs, dp.v_capped,
                     dp.v_cap_pairs > 0 ? dp.v_cap_fsum / dp.v_cap_pairs : 0.0,
                     dp.v_cap_fmax, dp.v_cap_fwmax,
                     dp.v_cap_pairs > 0 ? std::exp(dp.v_cap_drop - dp.v_cap_all) : 0.0,
                     mean(dp.rz_sum[0], dp.rz_n[0]), dp.rz_max[0], mean(dp.rz_sum[1], dp.rz_n[1]),
                     dp.rz_max[1], mean(dp.rz_sum[2], dp.rz_n[2]), dp.rz_max[2], re, rb / 1048576.0,
                     static_cast<int>(g_rfv), ve, vb / 1048576.0, dp.rfv_look, dp.rfv_miss,
                     dp.rfv_seed, dp.rfv_od, dp.rfv_od_ent, dp.rfv_child_miss, dp.rfv_keys,
                     dp.rfv_key_ent, tb.c_str());
      }
    }
    dp.pr_v_finish(sid);
    if (pruned_f) {
      for (int t = 0; t < BT_N; ++t) {
        std::fprintf(pruned_f, "P %s %s events %lld dropped %lld matched %lld w %.9e pos",
                     sid.c_str(), BT_NAMES[t], drops.events[t], drops.dropped[t],
                     drops.matched[t], drops.w[t]);
        for (int k = 0; k < 10; ++k) std::fprintf(pruned_f, " %.6e", drops.wpos[t][k]);
        std::fprintf(pruned_f, "\n");
      }
      std::fflush(pruned_f);
      if (kept_f) std::fflush(kept_f);
    }
    if (gcfg.on || gate_stats_out) {
      const GateStats& gs = dp.gate_stats;
      std::string body;
      const char* names[2] = {"Qg", "Qgl"};
      for (int t = 0; t < 2; ++t)
        body += std::string(" ") + names[t] + " cand " + std::to_string(gs.cand[t]) + " elig " +
                std::to_string(gs.elig[t]) + " ranked " + std::to_string(gs.ranked[t]) +
                " kept " + std::to_string(gs.kept[t]) + " reserved " +
                std::to_string(gs.reserved[t]);
      std::fprintf(stderr, "[gate] %s n %d on %d tau %.17g span %d%s\n",
                   id.empty() ? "-" : id.c_str(), static_cast<int>(seq.size()), gcfg.on ? 1 : 0,
                   gcfg.tau, gcfg.span, body.c_str());
    }
    if (arm_score) {
      const ArmStats& st = dp.arm_stats;
      std::string hist;
      for (int k = 0; k < LK_HIST_BINS; ++k)
        hist += (k ? "," : "") + std::to_string(st.lk_hist[k]);
      std::fprintf(stderr,
                   "[arm] %s n %d scored %lld arm_wins %lld frac %.6f mean_log_kappa %.6f "
                   "max_log_kappa %.6f windows %lld floor_windows %lld chosen_floor %lld "
                   "lk_hist %s\n",
                   id.empty() ? "-" : id.c_str(), static_cast<int>(seq.size()), st.scored,
                   st.arm_wins, st.scored ? static_cast<double>(st.arm_wins) / st.scored : 0.0,
                   st.arm_wins ? st.sum_log_kappa / st.arm_wins : 0.0, st.max_log_kappa,
                   st.windows, st.floor_windows, st.chosen_floor, hist.c_str());
    }
    if (stem_rank)
      std::fprintf(stderr,
                   "[stem] %s n %d loop %d role %s entries %zu dropped %zu nodes %zu bytes %zu "
                   "removed_qg %lld removed_qgl %lld removed_qgr %lld sec %.6f\n",
                   sid.c_str(), static_cast<int>(seq.size()), stem.loop,
                   stem_qg2 ? "all" : stem_both ? "both" : "stem1", stem.entries, stem.dropped,
                   stem.index.nodes() + (stem.mirror ? stem.index2.nodes() : 0),
                   stem.bytes(), dp.stem_removed[0], dp.stem_removed[1], dp.stem_removed[2],
                   stem_sec);
    if (stem_rank && (stem_hr_av || stem_hr_phi))
      std::fprintf(stderr,
                   "[stem-hr] %s n %d mode %s gamma %.6g entries %zu dropped %zu nodes %zu "
                   "qg_states %lld av_only %lld orig_only %lld mean_delta %.9g "
                   "min_delta %.9g max_delta %.9g check_states %lld check_mismatch %lld "
                   "check_max_rel %.3e\n",
                   sid.c_str(), static_cast<int>(seq.size()), stem_hr_av ? "av" : "phi",
                   stem_hr_gamma, stem.entries_av,
                   stem.dropped_av, stem.index_av.nodes(), dp.hr_n, dp.hr_av_only,
                   dp.hr_orig_only, dp.hr_n ? dp.hr_dsum / dp.hr_n : 0.0, dp.hr_dmin,
                   dp.hr_dmax, dp.hrc_states, dp.hrc_mis, dp.hrc_rel);
    if (stem_rank && (stem_qg_r_set || stem_consts_set)) {
      // Pseudoknot factors each role adds on top of alpha and the outside estimate (H_R,
      // the mirrored H_L, or the partner tables); at_in = at(u, v), at_out = at(p, q).
      const char* qg_l = stem_consistent ? "beta1+2beta2+at_in+at_out" : "beta1+2beta2";
      const char* qg_r = stem_qg_partner
                             ? (stem_consistent ? "Bg[beta2+at_in]+beta1+at_out" : "Bg[beta2+at_in]")
                             : (stem_consistent ? "beta1+2beta2+at_in+at_out" : "beta1+2beta2");
      const char* gl_l = stem_consistent ? "beta1+beta2+at_out" : "beta1+2beta2";
      std::fprintf(stderr,
                   "[stem-consts] %s n %d qg_r %s consts %s Qg.L %s Qg.R %s Qgl.L %s Qgl.R %s "
                   "Qgr.R %s qg_score_sec %.6f\n",
                   sid.c_str(), static_cast<int>(seq.size()),
                   stem_qg_partner ? "partner" : "mirror", stem_consistent ? "consistent" : "as-is",
                   qg_l, stem_qg2 ? qg_r : "-", gl_l, stem_both ? "beta1+at_out" : "-",
                   stem_both ? "beta1+at_out" : "-", dp.stem_qg_sec);
    }
    if (stem_rank && stem_ext) {
      // loop / helpers only: growth and helper counts (Result 19)
      char grow[1024] = "";
      if (stem_ext == 3 || stem_helpers) {
        const long long hs = dp.hlp_states[0] + dp.hlp_states[1];
        const double lay = static_cast<double>(dp.grow_layer);
        std::snprintf(
            grow, sizeof grow,
            " lg %d cut %.3g cap %d helpers %s ml_span %d ranked_qg %lld ranked_qgls %lld "
            "ranked_qgrs %lld lists %lld list_nodes %lld nodes_per_list %.4f layer_nodes %lld "
            "cut_drop %lld cap_drop %lld cut_frac %.6f cap_frac %.6f helper_closings %lld "
            "closings_per_helper %.4f helper_nodes %lld helper_finite %lld helper_sec %.6f",
            stem_ext == 3 ? grow_lg : -1, grow_cut, grow_nodes,
            stem_helpers == 1 ? "grow" : stem_helpers == 2 ? "nested" : "off",
            ml_span, dp.ext_states, dp.hlp_states[0], dp.hlp_states[1], dp.grow_lists,
            dp.grow_list_nodes,
            dp.grow_lists ? static_cast<double>(dp.grow_list_nodes) / dp.grow_lists : 0.0,
            dp.grow_layer, dp.grow_cut_drop, dp.grow_cap_drop,
            lay > 0 ? dp.grow_cut_drop / lay : 0.0, lay > 0 ? dp.grow_cap_drop / lay : 0.0,
            dp.hlp_closings, hs ? static_cast<double>(dp.hlp_closings) / hs : 0.0,
            dp.hlp_nodes, dp.hlp_finite, dp.hlp_sec);
      }
      // --stem-est-profile (Result 21): keys per step, mean key-profile support, entries
      // evaluated per scored state, Bg lookups / raw computations / partner entries visited,
      // profile build and per-state scoring time
      char prof[1024] = "";
      if (stem_profile) {
        const long long keys = dp.pf_keys[0] + dp.pf_keys[1] + dp.pf_keys[2];
        std::snprintf(
            prof, sizeof prof,
            " profile on eps %.3g tail %s steps %lld keys_qg %lld keys_qgls %lld keys_qgrs %lld "
            "keys_per_step %.4f pair_profiles %lld mean_ly %.4f mean_rx %.4f trunc_drop %lld "
            "scored %lld evals_per_state %.4f tail_closings %lld bg_lookups %lld bg_misses %lld "
            "bg_visits %lld bg_visits_per_miss %.4f prof_build_sec %.6f prof_score_sec %.6f",
            prof_eps, ml_tail ? "nested" : "off", dp.pf_steps, dp.pf_keys[0], dp.pf_keys[1],
            dp.pf_keys[2], dp.pf_steps ? static_cast<double>(keys) / dp.pf_steps : 0.0,
            dp.pf_pairs, keys ? static_cast<double>(dp.pf_ly) / keys : 0.0,
            keys ? static_cast<double>(dp.pf_rx) / keys : 0.0, dp.pf_trunc, dp.pf_states,
            dp.pf_states ? static_cast<double>(dp.pf_evals) / dp.pf_states : 0.0, dp.pf_tail_cl,
            dp.bg_lookups, dp.bg_misses, dp.bg_visits,
            dp.bg_misses ? static_cast<double>(dp.bg_visits) / dp.bg_misses : 0.0,
            dp.pf_build_sec, dp.pf_total_sec - dp.pf_build_sec);
        if (bg_reorder) {
          const size_t len = std::strlen(prof);
          std::snprintf(prof + len, sizeof prof - len,
                        " bg_reorder on m_rows %lld m_build_visits %lld m_eval_visits %lld",
                        dp.ro_rows, dp.ro_build, dp.ro_eval);
        }
      }
      // --stem-est-chain (Result 24): kept pairs, chain nodes, steps by kind, root profiles and
      // their mean support, entries cut; chain time (part of prof_build_sec), nested pass time
      char chn[512] = "";
      if (stem_chain)
        std::snprintf(
            chn, sizeof chn,
            " chain on floor %.3g span %d loop %d cut %.3g ml %s kept %zu chain_nodes %lld "
            "edges_b %lld edges_f %lld edges_g %lld edges_h %lld roots %lld root_ly %.4f "
            "root_rx %.4f chain_cut_drop %lld chain_sec %.6f chain_prep_sec %.6f",
            chain_floor, chain_span, chain_loop, chain_cut, chain_ml ? "on" : "off", stem.n_kept,
            dp.ch_nodes, dp.ch_edges[0], dp.ch_edges[1], dp.ch_edges[2], dp.ch_edges[3],
            dp.ch_roots, dp.ch_roots ? static_cast<double>(dp.ch_ly) / dp.ch_roots : 0.0,
            dp.ch_roots ? static_cast<double>(dp.ch_rx) / dp.ch_roots : 0.0, dp.ch_cut_drop,
            dp.ch_sec, chain_prep_sec);
      // --stem-est-principled: ranked states and scoring time per table, key-profile entries
      // (G endpoints) per state, arms per T, T / K / Stop builds and hits, per-q' arm
      // structures, static sizes, pre-pass kept pairs and time, Bml check
      std::string prs;
      if (stem_pr) {
        long long st_all = 0;
        for (int t = 0; t < 6; ++t) st_all += dp.pr_states[t];
        const long long st4 = st_all - dp.pr_states[5];
        char b[2048];
        std::snprintf(
            b, sizeof b,
            " principled on floor %.3g span %d ml_span %d ml_close %d cut %.3g arms_k %d "
            "kept %zu states_qg %lld states_qgl %lld states_qgr %lld states_qgls %lld "
            "states_qgrs %lld states_qp %lld sec_qg %.6f sec_qgl %.6f sec_qgr %.6f "
            "sec_qgls %.6f sec_qgrs %.6f sec_qp %.6f endpoints_per_state %.4f "
            "arms_per_T %.4f arms_skipped %lld ys_per_T %.4f T_built %lld T_hits %lld "
            "K_built %lld K_hits %lld i_per_K %.4f stop_built %lld stop_hits %lld "
            "armq_built %lld V_entries %lld arms %lld arm_entries %lld G_nodes %lld "
            "G_entries %lld roots %lld keys %lld key_entries %lld cut_drop %lld "
            "static_sec %.6f pr_prep_sec %.6f bml_check_cells %lld bml_check_mis %lld "
            "bml_check_rel %.3e",
            pr_floor, 0, pr_ml_span, pr_ml_close, pr_cut, pr_arms_k, stem.pr_n_kept,
            dp.pr_states[0], dp.pr_states[1], dp.pr_states[2], dp.pr_states[3],
            dp.pr_states[4], dp.pr_states[5], dp.pr_sec[0], dp.pr_sec[1], dp.pr_sec[2],
            dp.pr_sec[3], dp.pr_sec[4], dp.pr_sec[5],
            st4 ? static_cast<double>(dp.pr_evals) / st4 : 0.0,
            dp.pr_T_built ? static_cast<double>(dp.pr_T_arms) / dp.pr_T_built : 0.0,
            dp.pr_T_skip,
            dp.pr_T_built ? static_cast<double>(dp.pr_T_ys) / dp.pr_T_built : 0.0, dp.pr_T_built,
            dp.pr_T_hits, dp.pr_K_built, dp.pr_K_hits,
            dp.pr_K_built ? static_cast<double>(dp.pr_K_is) / dp.pr_K_built : 0.0,
            dp.pr_slr_miss, dp.pr_slr_hits, dp.pr_aq_built, dp.pr_V_entries, dp.pr_arms,
            dp.pr_arm_entries, dp.pr_G_nodes, dp.pr_G_entries, dp.pr_roots, dp.pr_keys,
            dp.pr_key_entries, dp.pr_cut_drop, dp.pr_static_sec, stem.pr_prep_sec,
            stem.bml_check_cells, stem.bml_check_mis, stem.bml_check_rel);
        prs = b;
        if (stem.pr_partner) {
          std::snprintf(b, sizeof b, " partner_floor %.3g partner_pairs %zu", pr_partner_floor,
                        stem.pr_n_partner);
          prs += b;
        }
        if (pr_qz_unp) prs += " test_qz unpaired";
        if (pr_pmax > 0) {
          std::snprintf(b, sizeof b,
                        " partner_maxpairs %d V_pairs %lld V_layer_entries %lld V_map_per_pair %.4f"
                        " V_layers_per_pair %.4f pcheck_n %lld pcheck_key_mis %lld pcheck_maxd %.3e"
                        " penum_n %lld penum_key_mis %lld penum_maxd %.3e",
                        pr_pmax, dp.pr_V_pairs, dp.pr_V_layer_entries,
                        dp.pr_V_pairs ? static_cast<double>(dp.pr_V_entries) / dp.pr_V_pairs : 0.0,
                        dp.pr_V_pairs ? static_cast<double>(dp.pr_V_layer_entries) / dp.pr_V_pairs : 0.0,
                        dp.pc_n, dp.pc_key_mis, dp.pc_maxd, dp.pe_n, dp.pe_key_mis, dp.pe_maxd);
          prs += b;
        }
        if (!pr_pml) {
          std::snprintf(b, sizeof b, " partner_ml off pmlcheck_n %lld pmlcheck_key_mis %lld"
                        " pmlcheck_maxd %.3e", dp.pm_n, dp.pm_key_mis, dp.pm_maxd);
          prs += b;
        }
        {
          auto avg = [](long long a, long long c) { return c ? static_cast<double>(a) / c : 0.0; };
          auto shr = [](double a, long long c) { return c ? a / c : 0.0; };
          std::snprintf(
              b, sizeof b,
              " grow_select %s grow_eps %.3g grow_B %d grow_sec %.6f sel_sec %.6f sel_profiles "
              "%lld sel_before %.4f sel_after %.4f sel_iota_share %.6f sel_zero %lld "
              "node_sel_before %.4f node_sel_after %.4f node_iota_share %.6f node_zero %lld "
              "uL_cols %lld uR_entries %lld U_evals %lld",
              pr_gsel ? "bound" : "floor", pr_geps, pr_gB, dp.pr_grow_sec, dp.pr_sel_sec,
              dp.gs_prof[1], avg(dp.gs_before[1], dp.gs_prof[1]),
              avg(dp.gs_after[1], dp.gs_prof[1]), shr(dp.gs_share[1], dp.gs_prof[1]),
              dp.gs_zero[1], avg(dp.gs_before[0], dp.gs_prof[0]),
              avg(dp.gs_after[0], dp.gs_prof[0]), shr(dp.gs_share[0], dp.gs_prof[0]),
              dp.gs_zero[0], dp.gs_cols, dp.gs_rent, dp.gs_U);
          prs += b;
          std::snprintf(b, sizeof b,
                        " grow_built %lld grow_int %lld grow_ml %lld fut_built %lld fut_int %lld "
                        "fut_ml %lld partner_sec %.6f aq_sep %d aq_sec %.6f aq_max %zu "
                        "T_sec %.6f K_sec %.6f SLR_sec %.6f W_sec %.6f W_built %lld W_ent %lld",
                        dp.st_nodes[0], dp.st_kind[0][0], dp.st_kind[0][1], dp.st_nodes[1],
                        dp.st_kind[1][0], dp.st_kind[1][1], dp.pr_partner_sec,
                        static_cast<int>(dp.pr_aq_sep), dp.pr_aq_sec, dp.pr_aq_max, dp.pr_sec_T,
                        dp.pr_sec_K, dp.pr_sec_SLR, dp.pr_sec_W, dp.pr_W_built, dp.pr_W_ent);
          prs += b;
          const double pidx_mb = dp.pidx_mb;
          std::snprintf(b, sizeof b,
                        " arm_persist %d arm_gaponly %d partner_span %d v_eps %.3g pidx_sec %.6f "
                        "pidx_nodes %lld pidx_mb %.3f pidx_upd %lld aq_arms %lld "
                        "aq_arms_per_build %.4f gap_epochs %lld gap_miss %lld v_before %lld "
                        "v_after %lld",
                        static_cast<int>(pr_arm_persist), static_cast<int>(pr_arm_gaponly), pr_pspan,
                        pr_veps, dp.pidx_sec, dp.pidx_nodes, pidx_mb, dp.pidx_upd, dp.aq_arms,
                        avg(dp.aq_arms, dp.pr_aq_built), dp.gap_epochs, dp.gap_miss, dp.v_before,
                        dp.v_after);
          prs += b;
          std::snprintf(b, sizeof b,
                        " fast %d hg_built %lld hg_entries %lld hg_terms %lld hg_clears %lld "
                        "hv_built %lld hv_entries %lld hv_terms %lld node_terms_g %lld "
                        "node_terms_v %lld",
                        pr_fast, dp.pr_hbuilt[0], dp.pr_hentries[0], dp.pr_hterms[0],
                        dp.pr_hg_clears, dp.pr_hbuilt[1], dp.pr_hentries[1], dp.pr_hterms[1],
                        dp.pr_nterms[0], dp.pr_nterms[1]);
          prs += b;
        }
        if (pr_arms_mode == 1) {
          const long long steps = dp.pr_sep_steps_all();
          std::snprintf(
              b, sizeof b,
              " arms sep arm_z %s sep_queries %lld sep_cols %lld sep_groups %lld "
              "groups_per_step %.4f arms_per_group %.4f queries_per_group %.4f sep_sec %.6f",
              pr_arm_z == 1 ? "unpaired" : "qz", dp.pr_sep_q, dp.pr_sep_cols, dp.pr_sep_groups,
              steps ? static_cast<double>(dp.pr_sep_groups) / steps : 0.0,
              dp.pr_sep_groups ? static_cast<double>(dp.pr_sep_garms) / dp.pr_sep_groups : 0.0,
              dp.pr_sep_groups ? static_cast<double>(dp.pr_sep_q) / dp.pr_sep_groups : 0.0,
              dp.pr_sep_sec);
          prs += b;
        }
      }
      std::fprintf(stderr,
                   "[stem-ext] %s n %d ext %s kmax %d combine %s states %lld nodes %lld "
                   "nodes_per_state %.4f finite %lld ext_sec %.6f qg_score_sec %.6f%s%s%s%s\n",
                   sid.c_str(), static_cast<int>(seq.size()),
                   stem_ext == 1 ? "stack" : stem_ext == 2 ? "bulge1" : stem_ext == 3 ? "loop" : "off",
                   stem_ext_max, stem_split ? "split" : "sum", dp.ext_states, dp.ext_nodes,
                   dp.ext_states ? static_cast<double>(dp.ext_nodes) / dp.ext_states : 0.0,
                   dp.ext_finite, dp.ext_sec, dp.stem_qg_sec, grow, prof, chn, prs.c_str());
      if (stem_pr) {
        auto qs = [](std::vector<int> v) {
          char t[96];
          if (v.empty()) {
            std::snprintf(t, sizeof t, "0/0/0");
            return std::string(t);
          }
          std::sort(v.begin(), v.end());
          double s = 0;
          for (int x : v) s += x;
          std::snprintf(t, sizeof t, "%.1f/%d/%d", s / static_cast<double>(v.size()),
                        v[static_cast<size_t>(0.9 * static_cast<double>(v.size() - 1))], v.back());
          return std::string(t);
        };
        std::string g = g_grow_rec ? "rec" : "prof";
        char b[1024];
        std::snprintf(b, sizeof b,
                      " ml_L %d ml_R %d ml_c %d diag %d U_nobml %d bml_built %d grow_sec %.6f "
                      "sel_sec %.6f pr_prep_sec %.6f partner_sec %.6f mem_nodes_mb %.3f "
                      "mem_roots_mb %.3f mem_help_peak_mb %.3f mem_step_peak_mb %.3f "
                      "mem_total_peak_mb %.3f help_entries_peak %lld help_clears %lld",
                      g_ml_L, g_ml_R, g_ml_c, static_cast<int>(g_ml_diag),
                      static_cast<int>(g_U_nobml), static_cast<int>(!stem.prBml.empty()),
                      dp.pr_grow_sec, dp.pr_sel_sec, stem.pr_prep_sec, dp.pr_partner_sec,
                      dp.gm_nodes / 1048576.0, dp.gm_roots / 1048576.0, dp.gm_help_peak / 1048576.0,
                      dp.gm_step_peak / 1048576.0, dp.gm_total_peak / 1048576.0,
                      dp.gm_help_entries_peak, g_grow_rec ? dp.rh_clears : dp.pr_hg_clears);
        g += b;
        const char* nm[3] = {"node", "root", "key"};
        for (int w = 0; w < 3; ++w) {
          g += std::string(" ") + nm[w] + "_n " + std::to_string(dp.gz_after[w].size());
          if (!g_grow_rec) g += std::string(" ") + nm[w] + "_before " + qs(dp.gz_before[w]);
          g += std::string(" ") + nm[w] + (g_grow_rec ? "_elems " : "_after ") + qs(dp.gz_after[w]);
          g += std::string(" ") + nm[w] + "_dp " + qs(dp.gz_dp[w]);
          g += std::string(" ") + nm[w] + "_dq " + qs(dp.gz_dq[w]);
          if (!g_grow_rec) {
            g += std::string(" ") + nm[w] + "_dp_before " + qs(dp.gz_dpb[w]);
            g += std::string(" ") + nm[w] + "_dq_before " + qs(dp.gz_dqb[w]);
          }
        }
        if (g_grow_rec) {
          std::snprintf(b, sizeof b,
                        " rn_nodes %lld rn_roots %lld rn_keys %lld rn_terms %lld rh_built %lld "
                        "rh_terms %lld mlc_drop_left %lld mlc_drop_right %lld",
                        dp.rn_nodes, dp.rn_roots, dp.rn_keys, dp.rn_terms, dp.rh_built, dp.rh_terms,
                        dp.rn_mlc_drop_l, dp.rn_mlc_drop_r);
          g += b;
        }
        if (g_tr_K > 0) {
          std::snprintf(b, sizeof b, " trunc_K %d trunc_eps %g trunc_mode %s", g_tr_K, g_tr_eps,
                        g_tr_pool ? "pool" : "side");
          g += b;
          const char* tn[2] = {"root", "key"};
          for (int w = 0; w < 2; ++w) {
            const double nn = static_cast<double>(dp.tr_before[w].size());
            g += std::string(" ") + tn[w] + "_tr_n " + std::to_string(dp.tr_before[w].size());
            g += std::string(" ") + tn[w] + "_tr_before " + qs(dp.tr_before[w]);
            g += std::string(" ") + tn[w] + "_tr_after " + qs(dp.tr_after[w]);
            g += std::string(" ") + tn[w] + "_tr_L " + qs(dp.tr_kL[w]);
            g += std::string(" ") + tn[w] + "_tr_R " + qs(dp.tr_kR[w]);
            std::snprintf(b, sizeof b, " %s_tr_share %.8f %s_tr_zero %.0f", tn[w],
                          nn > 0 ? dp.tr_share[w] / nn : 1.0, tn[w], dp.tr_zero[w]);
            g += b;
          }
        }
        if (g_rec_check) {
          std::snprintf(b, sizeof b, " check_states %lld check_max_dscore %.3e check_entries %lld "
                        "check_max_dvec %.3e", dp.rc_n, dp.rc_max_score, dp.rc_vec_n, dp.rc_max_vec);
          g += b;
        }
        if (g_ml_diag && dp.rd_n > 0) {
          const char* cn[DP03::RD_CFG] = {"all_all", "all_none", "depth1_none", "depth1_depth1",
                                          "none_none"};
          std::snprintf(b, sizeof b, " diag_states %lld", dp.rd_n);
          g += b;
          for (int c = 0; c < DP03::RD_CFG; ++c) {
            // kept-fraction quantiles from the -log10 histogram (bin lower edges)
            auto q = [&](double frac) {
              const long long goal = static_cast<long long>(std::ceil(frac * static_cast<double>(dp.rd_n)));
              long long cum = 0;
              for (int k = DP03::RD_BINS; k >= 0; --k) {
                cum += dp.rd_hist[c][k];
                if (cum >= goal) return k == DP03::RD_BINS ? 0.0 : std::pow(10.0, -k / 100.0);
              }
              return 1.0;
            };
            std::snprintf(b, sizeof b,
                          " kept_%s_mean %.6f kept_%s_p50 %.4g kept_%s_p10 %.4g kept_%s_p1 %.4g "
                          "kept_%s_min %.3g",
                          cn[c], dp.rd_sum[c] / static_cast<double>(dp.rd_n), cn[c], q(0.5), cn[c],
                          q(0.1), cn[c], q(0.01), cn[c], dp.rd_min[c]);
            g += b;
          }
          for (int side = 0; side < 2; ++side)
            for (int m = 0; m < 3; ++m) {
              std::snprintf(b, sizeof b, " share_%s%d %.6f", side ? "R" : "L", m,
                            dp.rd_layer_sum[side][m] / static_cast<double>(dp.rd_n));
              g += b;
            }
        }
        std::fprintf(stderr, "[growth] %s n %d mode %s\n", sid.c_str(), static_cast<int>(seq.size()),
                     g.c_str());
      }
    }
    if (partial) {
      std::fprintf(stderr, "[partial] %s n %d stopped_at_step %d phase %s\n", sid.c_str(),
                   static_cast<int>(seq.size()), dp.stop_j, dp.stop_phase == 1 ? "vbuild" : "scan");
      std::fflush(stderr);
      return 3;
    }
    if (pr_gcheck)
      std::fprintf(stderr,
                   "[pr-grow-check] %s n %d endpoints %lld viol_L %lld viol_R %lld "
                   "max_log_ratio_L %.6e max_log_ratio_R %.6e sel_above_all %lld\n",
                   sid.c_str(), static_cast<int>(seq.size()), dp.gc_n, dp.gc_vL, dp.gc_vR,
                   dp.gc_maxL, dp.gc_maxR, dp.gc_sel_gt);
    if (pr_check)
      std::fprintf(stderr,
                   "[pr-check] %s n %d states %lld/%lld/%lld/%lld/%lld/%lld mismatch "
                   "%lld/%lld/%lld/%lld/%lld/%lld max_rel %.3e/%.3e/%.3e/%.3e/%.3e/%.3e "
                   "bml_cells %lld bml_mismatch %lld bml_max_rel %.3e\n",
                   sid.c_str(), static_cast<int>(seq.size()), dp.prc_states[0], dp.prc_states[1],
                   dp.prc_states[2], dp.prc_states[3], dp.prc_states[4], dp.prc_states[5],
                   dp.prc_mis[0], dp.prc_mis[1], dp.prc_mis[2], dp.prc_mis[3], dp.prc_mis[4],
                   dp.prc_mis[5], dp.prc_rel[0], dp.prc_rel[1], dp.prc_rel[2], dp.prc_rel[3],
                   dp.prc_rel[4], dp.prc_rel[5], stem.bml_check_cells, stem.bml_check_mis,
                   stem.bml_check_rel);
    if (ro_check)
      std::fprintf(stderr,
                   "[bg-reorder-check] %s n %d states %lld mismatch %lld max_rel %.3e events %lld "
                   "flip_events %lld flip_states %lld\n",
                   sid.c_str(), static_cast<int>(seq.size()), dp.roc_states, dp.roc_mis,
                   dp.roc_rel, dp.roc_events, dp.roc_flip_events, dp.roc_flip_states);
    if (prof_check)
      std::fprintf(stderr,
                   "[profile-check] %s n %d qg_states %lld qg_mismatch %lld qg_max_rel %.3e "
                   "helper_states %lld helper_mismatch %lld helper_max_rel %.3e\n",
                   sid.c_str(), static_cast<int>(seq.size()), dp.pc_states[0], dp.pc_mis[0],
                   dp.pc_rel[0], dp.pc_states[1], dp.pc_mis[1], dp.pc_rel[1]);
    if (stem_ext_check) {
      char hc[256] = "";
      if (stem_helpers)
        std::snprintf(hc, sizeof hc, " helper_states %lld helper_mismatch %lld helper_max_rel %.3e",
                      dp.hc_states, dp.hc_mis, dp.hc_rel);
      std::fprintf(stderr, "[ext-check] %s n %d states %lld mismatch %lld max_rel %.3e%s\n",
                   sid.c_str(), static_cast<int>(seq.size()), dp.xc_states, dp.xc_mis, dp.xc_rel,
                   hc);
    }
    std::vector<std::tuple<int, int, double>> P;
    if (do_bpp) P = dp.bpp(cutoff);
    if (g_mem_log) {
      std::vector<MemItem> it;
      dp.mem_items(it, "", true, pk_l2r != 0);
      mem_log("outside_done", it);
    }
    if (dp.l2r_diag) {
      dp.outside();
      if (dp.l2r_split) dp.l2r_split_run();
      std::fprintf(stderr, "[l2r-diag] %s n %d logZ %.12f%s\n", sid.c_str(), static_cast<int>(seq.size()),
                   logz, dp.l2r_diag_report({200, 500, 1000, 2000, 4000}).c_str());
    }
    if (!l2r_odump.empty() && pk_l2r == 1) {
      if (!dp.l2r_diag) dp.outside();
      dp.l2r_oracle_dump(l2r_odump);
    }
    if (dp.ltr_f) {
      dp.outside();
      dp.ltr_final();
    }
    if (stem_bg_check) std::fprintf(stderr, "%s\n", dp.stem_bg_check(sid, logz).c_str());
    double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (stem_diag) {
      std::vector<const StemEst*> ests;
      for (const StemEst& s : diag_ests) ests.push_back(&s);
      std::string js = dp.stem_diag(ests, sid);
      std::fprintf(stem_diag_f, "%s\n", js.c_str());
      std::fflush(stem_diag_f);
    }
    if (gate_diag) {
      std::string js = dp.gate_diag(*gate, diag_taus, id.empty() ? "-" : id, gcfg.span,
                                    gcfg.delta, long_span);
      std::fprintf(diag_f, "%s\n", js.c_str());
      std::fflush(diag_f);
    }
    if (states_f) dp.dump_states(states_f, sid, dump_min_w);
    if (oracle_f) dp.pr_oracle_diag(stem, oracle_f, sid, pr_oracle_set);
    if (pr_sep_check) dp.pr_sep_report(stderr, sid);
    if (scores_f) dp.dump_scores(scores_f, sid, stem, dump_min_w);
    if (ext_f) {
      auto it = exact_states.find(sid);
      if (it == exact_states.end()) std::fprintf(stderr, "--dump-exact: no states for %s\n", sid.c_str());
      dp.dump_ext(ext_f, sid, stem, it == exact_states.end() ? no_states : it->second,
                  stem_ext_max);
    }

    if (const char* p2dump = std::getenv("DUMP_P2_ALL")) {
      dp.outside();
      dp.dump_p2_all(p2dump);
    }
    if (!id.empty()) std::printf(">%s\n", id.c_str());
    std::printf("logZ %.17g\n", logz);
    std::printf("sec %.6f\n", sec);
    std::printf("pairs %zu\n", P.size());
    for (const auto& t : P)
      std::printf("%d %d %.17g\n", std::get<0>(t), std::get<1>(t), std::get<2>(t));
    if (dump_qb_on) dp.dump_qb(stdout);
    std::fflush(stdout);
    id.clear();
