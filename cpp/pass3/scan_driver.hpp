// class DP03, public part 1: constructor (beams b4 and per-table beams), pfunc() = the inside
// pass over j (nested Qb/Qm/Qz/C, then the pseudoknot part: the 4-index exact tables or the
// L2R scan), bpp(), and declarations of the out-of-class dumps.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

 public:
  DP03(const std::string& seq, int beam, bool allow_pk, int loop_cap, bool score_qz,
       bool gap_step = false, const int* table_beam = nullptr)
      : E(seq), n(E.n), beam(beam > 0 ? beam : 0), allow_pk(allow_pk),
        cap(loop_cap), score_qz(score_qz), gap_step(gap_step) {
    // table_beam[t] < 0: the global beam; 0: infinite
    for (int t = 0; t < BT_N; ++t)
      tbeam[t] = table_beam && table_beam[t] >= 0 ? table_beam[t] : this->beam;
    alloc();
  }


  GateStats gate_stats;

  double pfunc() {
    // The nested pass (no pseudoknots) walks mates and does not build the pair list.
    if (l2r_mode != 0 || allow_pk) ensure_pairs();
    for (int j = 0; j < n; ++j) {
      cur_j = j;
      cpu_stop_check(0);
      if (allow_pk && l2r_mode != 1) {
        fill_Qg(j);
        prune4(Qg[j], BT_QG, j, GATE_QG);
        fill_Qgls(j);
        prune4(Qgls[j], BT_QGLS, j);
        fill_Qgl(j);
        prune4(Qgl[j], BT_QGL, j, GATE_QGL);
        fill_Qgr(j);
        prune4(Qgr[j], BT_QGR, j);
        fill_Qp(j);
        prune2(Qp, Qp_at[j], j, BT_QP, nullptr, nullptr);
      }
      if (allow_pk && l2r_mode != 0) {
        l2r_step_pk(j);
        if (l2r_mode == 1) {
          const int pre = Qp_at[j].live;
          prune2(Qp, Qp_at[j], j, BT_QP, nullptr, nullptr);
          const int kept = Qp_at[j].live;
          qp_pre_n += pre;
          qp_kept_n += kept;
          qp_pre_max = std::max(qp_pre_max, pre);
          qp_kept_max = std::max(qp_kept_max, kept);
          qp_beam = tbeam[BT_QP];
          if (tbeam[BT_QP] > 0 && pre > tbeam[BT_QP]) ++qp_bind_n;
        }
      }
      if (l2r_mode != 0) l2r_done[LT_QP] = j;
      {
        const auto tq = std::chrono::steady_clock::now();
        fill_Qb(j);
        sec_qb += std::chrono::duration<double>(std::chrono::steady_clock::now() - tq).count();
      }
      prune2(Qb, Qb_at[j], j, BT_QB, nullptr, nullptr);
      if (l2r_mode != 0) l2r_done[LT_QB] = j;
      if (l2r_mode == 0 && !allow_pk) {
        const auto ti = std::chrono::steady_clock::now();
        lp_push_interior(j);
        sec_qb += std::chrono::duration<double>(std::chrono::steady_clock::now() - ti).count();
      } else if (l2r_mode != 0) {
        build_qlead(j);
      }
      fill_Qm(j);
      // Qm2 is dropped with Qm (its keys are a subset of Qm's), unless Qm2's own beam is
      // infinite; then it gets its own beam.  With equal beams the second call is a no-op.
      {
        const int pre = Qm_at[j].live;
        if (tbeam[BT_QM2] > 0) prune2(Qm, Qm_at[j], j, BT_QM, &Qm2, &Qm2_at[j]);
        else prune2(Qm, Qm_at[j], j, BT_QM, nullptr, nullptr);
        const int kept = Qm_at[j].live;
        qm_pre_n += pre;
        qm_kept_n += kept;
        qm_pre_max = std::max(qm_pre_max, pre);
        qm_kept_max = std::max(qm_kept_max, kept);
        qm_beam = tbeam[BT_QM];
        if (tbeam[BT_QM] > 0 && pre > tbeam[BT_QM]) ++qm_bind_n;
      }
      prune2(Qm2, Qm2_at[j], j, BT_QM2, nullptr, nullptr);
      if (l2r_mode == 0 && !allow_pk) lp_emit_multi(j);
      if ((allow_pk || chain_qmp) && E.stem_ml) {
        fill_Qmp(j);
        prune2(Qmp, Qmp_at[j], j, BT_QMP, nullptr, nullptr);
      }
      qm_reindex(j);
      if (l2r_mode != 0) l2r_done[LT_QM] = j;
      {
        const auto tq = std::chrono::steady_clock::now();
        if (qz_beam_fill) fill_Qz_kept(j);
        else fill_Qz(j);
        const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - tq).count();
        sec_qz += dt;
      }
      if (l2r_mode != 0) l2r_done[LT_QZ] = j;
      {
        const auto tq = std::chrono::steady_clock::now();
        fill_C(j);
        sec_c += std::chrono::duration<double>(std::chrono::steady_clock::now() - tq).count();
      }
      advance_pk_inside(j);
      if (allow_pk && l2r_mode != 1) {
        fill_Qgrs(j);
        prune4(Qgrs[j], BT_QGRS, j);
      }
      if (allow_pk && l2r_mode != 0) l2r_step_end(j);
      if (l2r_mode != 0) note_col();
    }
    inside_done = true;
    logZ = n > 0 ? C[n - 1] : 0.0;
    if (n > 0 && l2r_mode == 0) {
      long long qb_n = 0, qm_n = 0;
      for (int j = 0; j < n; ++j) {
        qb_n += Qb_at[j].live;
        qm_n += Qm_at[j].live;
      }
      const double dn = static_cast<double>(n);
      std::fprintf(stderr,
                   "[p1rule] n %d lead %d hair %.3f int %.3f m2xP %.3f m2multi %.3f multi_ext %.3f "
                   "m_m_u %.3f m_p %.3f extC %.3f qz_cand %.3f qb_live %.3f qm_live %.3f\n",
                   n, ml_lead_cap, r_hair / dn, r_int / dn, r_m2xp / dn, r_m2multi / dn,
                   r_muext / dn, r_mmu / dn, r_mp / dn, r_extc / dn,
                   static_cast<double>(qz_cand) / dn, qb_n / dn, qm_n / dn);
      std::fprintf(stderr, "[p1sec] qb %.6f qm %.6f branch %.6f qz %.6f c %.6f\n", sec_qb, sec_qm,
                   sec_branch, sec_qz, sec_c);
      const double ms = 1000.0 / dn;
      const double rule = sec_hair + sec_int + sec_m2xp + sec_emit + sec_mclose + sec_mp + sec_c + sec_qz;
      const double pfunc_s = sec_qb + sec_qm + sec_branch + sec_qz + sec_c;
      std::fprintf(stderr,
                   "[p1ms] hair %.4f int %.4f m2xp %.4f emit %.4f mclose %.4f m %.4f c %.4f "
                   "qz %.4f sum %.4f pfunc %.4f\n",
                   sec_hair * ms, sec_int * ms, sec_m2xp * ms, sec_emit * ms, sec_mclose * ms,
                   sec_mp * ms, sec_c * ms, sec_qz * ms, rule * ms, pfunc_s * ms);
      std::fprintf(stderr, "[p1z] %.17g\n", logZ);
    }
    if (grammar_simple && n > 0) {
      const double dn = static_cast<double>(n);
      std::fprintf(stderr,
                   "[grammar] simple n %d fill %lld int %lld M1 %lld M2 %lld M3 %lld M4 %lld "
                   "M5 %lld M6 %lld M7 %lld M8 %lld phase %lld lead %lld overcap %lld\n",
                   n, gs_fill, it_int[0] + it_int[1], gs_m[0], gs_m[1], gs_m[2], gs_m[3], gs_m[4],
                   gs_m[5], gs_m[6], gs_m[7], gs_phase, gs_lead, gs_over);
      std::fprintf(stderr,
                   "[grammar] per_nt fill %.3f int %.3f M1 %.3f M2 %.3f M3 %.3f M4 %.3f M5 %.3f "
                   "M6 %.3f M7 %.3f M8 %.3f phase %.3f lead %.3f\n",
                   gs_fill / dn, (it_int[0] + it_int[1]) / dn, gs_m[0] / dn, gs_m[1] / dn,
                   gs_m[2] / dn, gs_m[3] / dn, gs_m[4] / dn, gs_m[5] / dn, gs_m[6] / dn,
                   gs_m[7] / dn, gs_phase / dn, gs_lead / dn);
      static const char* tn[6] = {"R^M", "R^M2", "R^Multi", "U^M", "U^M2", "U^Multi"};
      std::fprintf(stderr, "[grammar-kept] n %d", n);
      for (int t = 0; t < 6; ++t)
        std::fprintf(stderr, " %s %.3f", tn[t], gs_tab_kept[t] / dn);
      std::fprintf(stderr, "\n");
    }
    return logZ;
  }

  std::vector<std::tuple<int, int, double>> bpp(double cutoff) {
    const auto t_out = std::chrono::steady_clock::now();
    outside();
    sec_outside += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_out).count();
    const auto t_pairs = std::chrono::steady_clock::now();
    std::vector<std::tuple<int, int, double>> p = collect_bpp(cutoff);
    sec_bpp += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_pairs).count();
    return p;
  }

