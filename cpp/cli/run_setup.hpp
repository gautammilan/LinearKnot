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
    if (restrict_pairs) {
      std::string pl;
      if (!std::getline(std::cin, pl)) {
        std::fprintf(stderr, "--restrict-pairs: missing pair line after the sequence\n");
        return 2;
      }
      restrict_list.clear();
      std::istringstream ps(pl);
      int a, b;
      while (ps >> a >> b) restrict_list.emplace_back(std::min(a, b) - 1, std::max(a, b) - 1);
      allowed_pairs = &restrict_list;
    }
    if (beam_audit(static_cast<int>(seq.size()))) return 2;

    auto t0 = std::chrono::steady_clock::now();
    const double cpu0 = ru_cpu();
    double pass1_s = 0, pass2_s = 0, scan_s = 0;
    double cpu_pass1 = 0, cpu_pass2 = 0;
    AccessRec access;
    if (need_access) access = load_access(access_idx, seq);
    std::unique_ptr<Gate> gate;
    if (gcfg.on || gate_diag)
      gate.reset(new Gate(static_cast<int>(seq.size()), access, gcfg.span, gcfg.delta, gcfg.both));
    apply_max_loop(loop_cap_scan);
    DP03 dp(seq, beam, allow_pk, loop_cap_scan, score_qz, arm_score, arm_score ? &access : nullptr,
            max_gap_scan, max_arm_len, gap_step, &gcfg, gcfg.on ? gate.get() : nullptr,
            table_beam, states_f != nullptr || dump_scores);
    dp.qz_beam_fill = !qzbar_old;
    dp.p2_qz_beam = b2_width;
    dp.qzbar_old = qzbar_old;
    if (qzbar_old)
      std::fprintf(stderr, "[qzbar-old] dense Qz→prQz (no beamed Qz_pre)\n");
    const std::string sid = id.empty() ? "-" : id;
    g_mem_sid = sid;
    if (g_mem_log) {
      std::vector<MemItem> it;
      dp.mem_items(it, "", true, false);
      mem_log("dp_alloc", it);
    }
    // Stem estimate, principled pre-pass and chain all use b1. Q^z uses b2.
    StemEst stem;
    std::vector<StemEst> diag_ests;
    double stem_sec = 0.0, chain_prep_sec = 0.0;
    // b1 is the nested beam at every length. The n<=1000 lean skip stays off once b1 is set.
    const bool lean = l2r_lean && pk_l2r == 1 && (l2r_beam < 0 ? beam : l2r_beam) > 0 &&
                      stem_rank && stem_pr && g_grow_rec && !stem_diag && !dump_scores &&
                      !dump_ext && !stem_hr_av && !stem_hr_phi && !stem_bg_check && nest_beam < 0;
    if (l2r_lean && !lean)
      std::fprintf(stderr, "[l2r-lean] %s not applicable to this run's flags: full pre-pass\n",
                   sid.c_str());
    stem.lean = lean;
    if (stem_rank || stem_diag || dump_scores || dump_ext) {
      auto ts = std::chrono::steady_clock::now();
      std::unique_ptr<DP03> nest;
      if (!lean) {
        const int scan_nb = nest_beam;
        apply_max_loop(loop_cap_p1);
        nest.reset(new DP03(seq, scan_nb, false, loop_cap_p1, score_qz, false, nullptr,
                            max_gap_scan, max_arm_len, true, nullptr, nullptr, nullptr, false));
        nest->chain_qmp = stem_chain;
        nest->qz_beam_fill = !qzbar_old;
        nest->p2_qz_beam = b2_width;
        nest->qzbar_old = qzbar_old;
        std::fprintf(stderr, "[nest-beam] b1 %d\n", scan_nb);
        if (pass1_only && !sample_prof.empty()) sample_start(sample_prof.c_str());
        const auto t_stem = std::chrono::steady_clock::now();
        g_beam_stats.set_pass("stem");
        nest->pfunc();
        pass1_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_stem).count();
        nest->log_nest("stem");
        if (g_mem_log) {
          std::vector<MemItem> it;
          nest->mem_items(it, "nest_", true, false);
          mem_log("nest_pfunc", it);
        }
      }
      if (stem_rank) {
        if (!lean) {
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
          if (g_mem_log) {
            std::vector<MemItem> it;
            stem_mem_items(stem, it);
            mem_log("stem_est", it);
          }
          if (stem_helpers == 2 || ml_tail) {
            nest->build_nested_beta(stem);
            if (g_mem_log) {
              std::vector<MemItem> it;
              nest->mem_items(it, "nest_", true, false);
              stem_mem_items(stem, it);
              mem_log("nested_beta", it);
            }
          }
        }
        dp.set_stem(&stem, stem_both, stem_qg2, stem_qg_partner, stem_consistent);
        dp.set_stem_ext(stem_ext, stem_ext_max, stem_split, stem_ext_check);
        dp.set_stem_grow(grow_lg, grow_cut, grow_nodes, stem_helpers, ml_span);
        if (stem_profile) dp.set_stem_profile(true, prof_eps, ml_tail, prof_check);
        if (bg_reorder) dp.set_stem_bg_reorder(true, ro_check);
        // exact nested-only run, shared by the chain and the principled pre-pass when b1 is 0
        std::unique_ptr<DP03> nexact;
        auto exact_nested = [&]() -> DP03& {
          if (!nexact) {
            apply_max_loop(loop_cap_p1);
            nexact.reset(new DP03(seq, 0, false, loop_cap_p1, score_qz, false, nullptr, max_gap_scan,
                                  max_arm_len, true, nullptr, nullptr, nullptr, false));
            nexact->chain_qmp = true;
            nexact->keep_closed_b = stem_pr_nested;
            nexact->prof_on = profile;
            nexact->qz_beam_fill = !qzbar_old;
            nexact->p2_qz_beam = b2_width;
            nexact->qzbar_old = qzbar_old;
            const auto t_ex = std::chrono::steady_clock::now();
            g_beam_stats.set_pass("pass1");
            nexact->pfunc();
            pass1_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_ex).count();
            nexact->log_nest("pass1");
            if (g_mem_log) {
              std::vector<MemItem> it;
              nexact->mem_items(it, "exact_", true, false);
              mem_log("exact_pfunc", it);
            }
          }
          return *nexact;
        };
        if (stem_pr || stem_pr_nested) {
          // Nested pre-pass at b1. --pre-beam replaces b1 when it is set.
          const auto tp = std::chrono::steady_clock::now();
          std::unique_ptr<DP03> nbeam;
          const int nb = pre_beam > 0 ? pre_beam : nest_beam;
          const bool use_beam = nb > 0;
          std::fprintf(stderr, "[nest-beam] pre %d\n", use_beam ? nb : 0);
          if (use_beam) {
            apply_max_loop(loop_cap_p1);
            nbeam.reset(new DP03(seq, nb, false, loop_cap_p1, score_qz, false, nullptr, max_gap_scan,
                                 max_arm_len, true, nullptr, nullptr, nullptr, false));
            nbeam->chain_qmp = true;
            nbeam->keep_closed_b = stem_pr_nested;
            nbeam->prof_on = profile;
            nbeam->qz_beam_fill = !qzbar_old;
            nbeam->p2_qz_beam = b2_width;
            nbeam->qzbar_old = qzbar_old;
            nbeam->pre_key_inside = pre_key == 1;
            nbeam->p2_sweep = p2_sweep;
            const auto t_pre = std::chrono::steady_clock::now();
            g_beam_stats.set_pass("pass1");
            nbeam->pfunc();
            pass1_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_pre).count();
            nbeam->log_nest("pass1");
            if (pass1_only) {
              if (!sample_prof.empty()) sample_stop();
              std::fprintf(stderr,
                           "[pass1-only] pass1_s %.6f stem_idx %zu pr_Mbar 0 pr_QbIn 0 grow_roots 0 "
                           "hwm_mb %.1f\n",
                           pass1_s, stem.entries, hwm_mb());
              nbeam->log_pr_side();
              return 0;
            }
          }
          DP03& np = nbeam ? *nbeam : exact_nested();
          apply_max_loop(loop_cap_p2);
          np.set_loop_cap(loop_cap_p2);
          auto stamp_p2 = [&](DP03& o) {
            o.qz_beam_fill = !qzbar_old;
            o.p2_qz_beam = b2_width;
            o.qzbar_old = qzbar_old;
            o.p2_qz_key = 1;
            o.p2_v_mode = p2_v_mode;
            o.p2_v_bp = p2_v_bp;
            o.p2_v_cap = p2_v_cap;
            o.p2_v_trim = p2_v_trim;
            o.p2_g_mode = p2_g_mode;
            o.p2_g_bp = p2_g_bp;
            o.p2_g_suffix = p2_g_suffix;
            o.p2_g_eager = p2_g_eager;
            o.p2_g_fb = p2_g_fb;
            o.p2_g_cap = p2_g_cap;
            o.sw_check = sw_check;
            o.p2_sweep = p2_sweep;
            o.p2_vsweep = p2_vsweep;
            o.p2_v_floor = p2_v_floor;
            o.p2_v_cut_fb = p2_v_cut_fb;
            o.lazy_ml_cap = lazy_ml_cap;
            o.prof_on = profile;
          };
          stamp_p2(np);
          stamp_p2(dp);
          stem.want_lp = oracle_f != nullptr || !pr_v_diag.empty() || !g_v2dump.empty() || g_v2ratio;
          stem.skip_kept = lean && oracle_f == nullptr && !pr_check;
          stem.skip_qbin = lean && l2r_p1 == 0 && !l2r_ctx_p1 && !pr_v_ratio &&
                           pr_v_diag.empty() && pr_v_sample <= 0.0 && g_v2dump.empty();
          cpu_pass1 = ru_cpu() - cpu0;
          const double hwm_pass1 = hwm_mb();
          const int nseq = static_cast<int>(seq.size());
          std::fprintf(stderr,
                       "[stage] pass1 n %d cpu %.3f cpu_s_nt %.6f hwm_mb %.1f hwm_mb_nt %.4f\n",
                       nseq, cpu_pass1, nseq > 0 ? cpu_pass1 / nseq : 0.0, hwm_pass1,
                       nseq > 0 ? hwm_pass1 / nseq : 0.0);
          const auto t_p2 = std::chrono::steady_clock::now();
          const double cpu_p2_0 = ru_cpu();
          g_beam_stats.set_pass("pass2");
          np.build_principled(stem, pr_floor, stem_pr_nested, pr_check, pr_partner_floor);
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
          if (g_mem_log) {
            std::vector<MemItem> it;
            np.mem_items(it, nbeam ? "beam_" : "exact_", true, false);
            stem_mem_items(stem, it);
            mem_log("principled", it);
          }
          if (stem_pr) {
            dp.set_stem_principled(true, chain_loop, chain_ml, pr_ml_span, pr_ml_close, pr_cut,
                                   pr_arms_k, pr_check);
            dp.set_stem_pr_select(pr_ratio, pr_arms_eps, pr_arms_max, pr_gapfree, pr_cut_roots);
            dp.set_stem_pr_sep(pr_arms_mode, pr_arm_z, pr_qz_unp, pr_sep_check);
            dp.set_stem_pr_grow(pr_gsel, pr_geps, pr_gB, pr_gcheck, pr_dump_f);
            dp.set_stem_pr_armq_sep(pr_aq_sep, pr_aq_max);
            dp.set_stem_pr_arms(pr_arm_persist, pr_arm_gaponly, pr_pspan, pr_veps);
            dp.pr_vsuffix = pr_vsuffix;
            dp.pr_vcap = pr_vcap;
            dp.pr_vmlk = pr_vmlk;
            dp.pr_vmlsteps = pr_vmlsteps;
            dp.pr_vxs = pr_vxs;
            dp.pr_vstream = pr_vstream;
            dp.set_stem_pr_partner_cap(pr_pmax, pr_pcheck);
            dp.set_stem_pr_partner_ml(pr_pml, pr_pmlcheck);
            dp.set_stem_pr_v(pr_v_ratio, pr_v_diag, pr_v_sample);
            dp.pr_fast = pr_fast;
            dp.pr_abl = pr_abl;
            dp.pr_abl_pk = pr_abl_pk;
            if (g_grow_rec) {
              dp.rnL = g_ml_diag ? 3 : g_ml_L + 1;
              dp.rnR = g_ml_diag ? 3 : g_ml_R + 1;
              dp.rec_noml = dp.rnL == 1 && dp.rnR == 1;
            }
          }
          if (stem_pr_nested) dp.set_stem_pr_nested(true);
          stem.pr_prep_sec =
              std::chrono::duration<double>(std::chrono::steady_clock::now() - tp).count();
        }
        if ((stem_chain && !lean) || stem_hr_av || stem_hr_phi) {
          const auto tc = std::chrono::steady_clock::now();
          // Chain and hairpin tables use the same b1 nested DP.
          if (!lean && nest) {
            if (stem_chain) nest->build_chain(stem, chain_floor);
            if (stem_hr_av) nest->build_hr_av(stem, stem_hr_gamma);
            if (stem_hr_phi) nest->build_hr_phi(stem, stem_hr_gamma);
          }
          if (stem_chain) dp.set_stem_chain(true, chain_span, chain_cut, chain_ml, chain_loop);
          if (stem_hr_av || stem_hr_phi) dp.set_stem_hr(true, stem_hr_av, stem_hr_check);
          chain_prep_sec =
              std::chrono::duration<double>(std::chrono::steady_clock::now() - tc).count();
          if (g_mem_log) {
            std::vector<MemItem> it;
            stem_mem_items(stem, it);
            mem_log("chain", it);
          }
        } else if (stem_chain) {
          dp.set_stem_chain(true, chain_span, chain_cut, chain_ml, chain_loop);
        }
      }
      if (dump_scores || dump_ext) nest->build_stem_est(stem_loop, stem, true);
      diag_ests.resize(stem_diag ? stem_diag_loops.size() : 0);
      for (size_t k = 0; k < diag_ests.size(); ++k)
        nest->build_stem_est(stem_diag_loops[k], diag_ests[k], true);
      stem_sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - ts).count();
    }
    if (g_mem_log) {
      std::vector<MemItem> it;
      stem_mem_items(stem, it);
      mem_log("stem_done", it);
    }
    DropLog drops;
    if (pruned_f) {
      drops.f = pruned_f;
      drops.id = sid;
      drops.n = std::max(1, static_cast<int>(seq.size()));
      drops.min_w = pruned_min_w;
      drops.kf = kept_f;
      if (!exact_path.empty()) {
        auto it = exact_states.find(sid);
        drops.exact = it == exact_states.end() ? &no_states : &it->second;
        if (it == exact_states.end())
          std::fprintf(stderr, "--dump-exact: no states for %s\n", sid.c_str());
      }
      dp.drops = &drops;
    }
    if (!remove_keys.empty()) dp.remove_keys = &remove_keys;
