// main(), per sequence, part 1: read the record, build DP03, pass 1 (nested pre-pass at b1)
// and the stem estimate.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

    if (!line.empty() && line[0] == '>') {
      id = line.substr(1);
      continue;
    }
    std::string seq = normalize(line);
    if (seq.empty()) continue;
    if (beam_audit(static_cast<int>(seq.size()))) return 2;

    auto t0 = std::chrono::steady_clock::now();
    const double cpu0 = ru_cpu();
    double pass1_s = 0, pass2_s = 0, scan_s = 0;
    double cpu_pass1 = 0, cpu_pass2 = 0;
    apply_max_loop(loop_cap_scan);
    DP03 dp(seq, beam, allow_pk, loop_cap_scan, score_qz, gap_step, table_beam);
    dp.qz_beam_fill = true;
    dp.p2_qz_beam = b2_width;
    const std::string sid = id.empty() ? "-" : id;
    // Stem estimate, principled pre-pass and chain all use b1. Q^z uses b2.
    StemEst stem;
    double stem_sec = 0.0, chain_prep_sec = 0.0;
    if (stem_rank) {
      auto ts = std::chrono::steady_clock::now();
      std::unique_ptr<DP03> nest;
      {
        const int scan_nb = nest_beam;
        apply_max_loop(loop_cap_p1);
        nest.reset(new DP03(seq, scan_nb, false, loop_cap_p1, score_qz, true));
        nest->chain_qmp = true;
        nest->qz_beam_fill = true;
        nest->p2_qz_beam = b2_width;
        std::fprintf(stderr, "[nest-beam] b1 %d\n", scan_nb);
        const auto t_stem = std::chrono::steady_clock::now();
        nest->pfunc();
        pass1_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_stem).count();
        nest->log_nest("stem");
      }
      if (stem_rank) {
        {
          const auto t_se = std::chrono::steady_clock::now();
          nest->build_stem_est(stem_loop, stem, stem_qg2);
          const double se_dn = static_cast<double>(seq.size());
          std::fprintf(stderr,
                       "[stem-est] sec %.6f entries %zu per_nt %.3f idx_mb %.1f t_mb %.1f hwm_mb %.1f\n",
                       std::chrono::duration<double>(std::chrono::steady_clock::now() - t_se).count(),
                       stem.entries, se_dn > 0 ? static_cast<double>(stem.entries) / se_dn : 0.0,
                       (12.0 * static_cast<double>(stem.entries + stem.entries2 + stem.entries_av) +
                        16.0 * static_cast<double>(stem.index.nodes() + stem.index2.nodes() +
                                                   stem.index_av.nodes())) /
                           1048576.0,
                       8.0 * se_dn * se_dn / 1048576.0, hwm_mb());
          if (stem_helpers == 2) {
            nest->build_nested_beta(stem);
          }
        }
        dp.set_stem(&stem, stem_both, stem_qg2, stem_qg_partner);
        dp.set_stem_ext(stem_ext);
        dp.set_stem_grow(grow_lg, stem_helpers);
        dp.set_stem_profile(true);
        // exact nested-only run, shared by the chain and the principled pre-pass when b1 is 0
        std::unique_ptr<DP03> nexact;
        auto exact_nested = [&]() -> DP03& {
          if (!nexact) {
            apply_max_loop(loop_cap_p1);
            nexact.reset(new DP03(seq, 0, false, loop_cap_p1, score_qz, true));
            nexact->chain_qmp = true;
            nexact->qz_beam_fill = true;
            nexact->p2_qz_beam = b2_width;
            const auto t_ex = std::chrono::steady_clock::now();
            nexact->pfunc();
            pass1_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_ex).count();
            nexact->log_nest("pass1");
          }
          return *nexact;
        };
        if (stem_pr) {
          // Nested pre-pass at b1.
          const auto tp = std::chrono::steady_clock::now();
          std::unique_ptr<DP03> nbeam;
          const int nb = nest_beam;
          const bool use_beam = nb > 0;
          std::fprintf(stderr, "[nest-beam] pre %d\n", use_beam ? nb : 0);
          if (use_beam) {
            apply_max_loop(loop_cap_p1);
            nbeam.reset(new DP03(seq, nb, false, loop_cap_p1, score_qz, true));
            nbeam->chain_qmp = true;
            nbeam->qz_beam_fill = true;
            nbeam->p2_qz_beam = b2_width;
            nbeam->p2_sweep = p2_sweep;
            const auto t_pre = std::chrono::steady_clock::now();
            nbeam->pfunc();
            pass1_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_pre).count();
            nbeam->log_nest("pass1");
          }
          DP03& np = nbeam ? *nbeam : exact_nested();
          apply_max_loop(loop_cap_p2);
          np.set_loop_cap(loop_cap_p2);
          auto stamp_p2 = [&](DP03& o) {
            o.qz_beam_fill = true;
            o.p2_qz_beam = b2_width;
            o.p2_qz_key = 1;
            o.p2_v_cap = p2_v_cap;
            o.p2_sweep = p2_sweep;
            o.p2_vsweep = p2_vsweep;
          };
          stamp_p2(np);
          stamp_p2(dp);
          cpu_pass1 = ru_cpu() - cpu0;
          const double hwm_pass1 = hwm_mb();
          const int nseq = static_cast<int>(seq.size());
          std::fprintf(stderr,
                       "[stage] pass1 n %d cpu %.3f cpu_s_nt %.6f hwm_mb %.1f hwm_mb_nt %.4f\n",
                       nseq, cpu_pass1, nseq > 0 ? cpu_pass1 / nseq : 0.0, hwm_pass1,
                       nseq > 0 ? hwm_pass1 / nseq : 0.0);
          const auto t_p2 = std::chrono::steady_clock::now();
          const double cpu_p2_0 = ru_cpu();
          np.build_principled(stem, pr_floor, pr_partner_floor);
          if (g_sweep == G_SWEEP_STATES && nbeam) {
            dp.g_pushes = nbeam->g_pushes;
            dp.g_pre_n = nbeam->g_pre_n;
            dp.g_pre_sum = nbeam->g_pre_sum;
            dp.g_pre_max = nbeam->g_pre_max;
            dp.g_post_n = nbeam->g_post_n;
            dp.g_post_sum = nbeam->g_post_sum;
            dp.g_post_max = nbeam->g_post_max;
            dp.g_merge_sec = nbeam->g_merge_sec;
            dp.g_cut_sec = nbeam->g_cut_sec;
          }
          pass2_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_p2).count();
          cpu_pass2 = ru_cpu() - cpu_p2_0;
          const double hwm_pass2 = hwm_mb();
          std::fprintf(stderr,
                       "[stage] pass2 n %d cpu %.3f cpu_s_nt %.6f hwm_mb %.1f hwm_mb_nt %.4f\n",
                       nseq, cpu_pass2, nseq > 0 ? cpu_pass2 / nseq : 0.0, hwm_pass2,
                       nseq > 0 ? hwm_pass2 / nseq : 0.0);
          stem.prof_qz_pre = np.prof_qz_sec;
          if (stem_pr) {
            dp.set_stem_principled(true, pr_ml_span, pr_cut);
            dp.set_stem_pr_select(pr_gapfree);
            dp.set_stem_pr_grow(pr_gsel, pr_geps);
            dp.set_stem_pr_armq_sep(pr_aq_sep);
            dp.set_stem_pr_arms(pr_veps);
            dp.pr_vsuffix = pr_vsuffix;
            dp.pr_vcap = pr_vcap;
            dp.pr_vmlk = pr_vmlk;
            dp.pr_vxs = pr_vxs;
            dp.pr_vstream = pr_vstream;
            dp.pr_fast = pr_fast;
            dp.pr_abl = pr_abl;
            dp.rnL = g_ml_L + 1;
            dp.rnR = g_ml_R + 1;
            dp.rec_noml = dp.rnL == 1 && dp.rnR == 1;
          }
          stem.pr_prep_sec =
              std::chrono::duration<double>(std::chrono::steady_clock::now() - tp).count();
        }
        {
          const auto tc = std::chrono::steady_clock::now();
          // Chain and hairpin tables use the same b1 nested DP.
          if (nest) {
            nest->build_chain(stem, chain_floor);
          }
          dp.set_stem_chain(true);
          chain_prep_sec =
              std::chrono::duration<double>(std::chrono::steady_clock::now() - tc).count();
        }
      }
      stem_sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - ts).count();
    }
