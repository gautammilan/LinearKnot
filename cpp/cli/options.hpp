// main(), part 1: option variables and their defaults.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  std::setvbuf(stderr, nullptr, _IOLBF, 0);
  std::signal(SIGXCPU, on_cpu_stop);
  std::signal(SIGUSR1, on_cpu_stop);
  int beam = 0;
  int loop_cap = MAX_LOOP;
  // Per-pass loop caps; -1 = use --loop-cap. NO_LOOP_CAP = uncapped. Also drives max_loop
  // for that pass.
  constexpr int NO_LOOP_CAP = 1000000;
  int loop_cap_p1 = -1, loop_cap_p2 = -1, loop_cap_scan = -1;
  bool score_qz = true;
  bool allow_pk = true;
  double cutoff = 1e-6;
  bool do_bpp = true;
  int pk_l2r = 0;
  bool l2r_spawn_cube = false;  // --l2r-spawn all|cube|cube-static|inc
  bool l2r_spawn_inc = false;
  int l2r_pbeam[3] = {1000, 1000, 1000};  // P1, P2, P3 beams, set to (b2, b4, b4)
  int l2r_p_hk = 64;                      // H entries per suffix
  bool l2r_p3_hex = false;                // --l2r-p3-hex on
  bool l2r_p3_lin = false;                // --l2r-p3-lin on
  bool ntr_z_set = false;
  bool grammar_simple = false;
  bool gap_step = false;
  bool loop_cap_set = false;
  int table_beam[BT_N];
  std::fill(table_beam, table_beam + BT_N, -1);
  bool stem_rank = false, stem_both = false, stem_qg2 = false;
  bool stem_qg_partner = false, stem_consistent = false, stem_qg_r_set = false;
  int stem_ext = 0, stem_ext_max = 30;
  bool stem_split = false, stem_ext_set = false;
  int grow_lg = 3, grow_nodes = 2000, ml_span = 30;
  double grow_cut = 1e-8;
  int stem_helpers = 0;
  bool grow_set = false, helpers_set = false;
  bool ml_tail = false;
  double prof_eps = 0.0;
  bool chain_ml = true;
  double chain_floor = 1e-7, chain_cut = 1e-8;
  int chain_span = 150, chain_loop = 30;
  bool stem_pr = false, pr_set = false;
  int pr_ml_span = 150, pr_ml_close = 30, pr_arms_k = 64;
  double pr_floor = 1e-7, pr_cut = 1e-8;
  bool pr_gapfree = false;
  // --pr-partner-floor (off: the kept pairs)
  double pr_partner_floor = -1.0;
  // --pr-grow-select floor|bound, --pr-grow-eps
  bool pr_gsel = false;
  bool pr_aq_sep = false;  // --pr-armq-sep on|off
  bool pr_arm_persist = false, pr_arm_gaponly = false;
  int pr_pspan = 0;
  int pr_fast = 0;                          // --pr-fast off|acc|ml
  int pr_abl = 0;                           // --pr-abl bml,q,p,xs,ks,sep
  double pr_veps = 1e-8;                                // --pr-v-eps (0: no trimming)
  bool pr_vsuffix = false;                              // --pr-v-trim global|suffix
  int pr_vcap = 0;                                      // --bV; default 256
  int p2_v_cap = 0;
  bool bV_set = false;
  int nest_beam = -1;
  bool b1_set = false, sweep_set = false, vsweep_set = false;
  int b2_width = 0;
  int b2_at = -1;
  bool p2_sweep = false, p2_vsweep = false;
  int pr_vmlk = 0;                                      // --pr-v-mlk K (0: all)
  bool pr_vxs = false;                                  // --pr-v-key xsys|xs
  bool pr_vstream = false;                              // --pr-v-stream off|on
  double pr_geps = 1e-6;
  int pr_gB = 512;
  int stem_loop = 2;

