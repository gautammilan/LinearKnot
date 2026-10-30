// main(), part 1: option variables and their defaults.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  std::setvbuf(stderr, nullptr, _IOLBF, 0);
  std::signal(SIGXCPU, on_cpu_stop);
  std::signal(SIGUSR1, on_cpu_stop);
  int beam = 0;
  int loop_cap = MAX_LOOP;
  // Temporary per-pass loop caps (--loop-cap-p1 / -p2 / -scan). -1 = use --loop-cap.
  // 0 / off / n = uncapped (NO_CAP). Also drives max_loop for that pass.
  constexpr int NO_LOOP_CAP = 1000000;
  int loop_cap_p1 = -1, loop_cap_p2 = -1, loop_cap_scan = -1;
  bool score_qz = true;
  bool allow_pk = true;
  double cutoff = 1e-6;
  bool do_bpp = true;
  int pk_l2r = 0, l2r_beam = -1, l2r_spawn_k = 0;
  bool l2r_spawn_cube = false, l2r_spawn_static = false;  // --l2r-spawn all|cube|cube-static|inc
  bool l2r_spawn_inc = false;
  int l2r_pbeam[3] = {1000, 1000, 1000};  // --l2r-pbeam b1,b2,b3 (P1, P2, P3; 0 or inf: none)
  int l2r_p_hk = 64;                      // --l2r-p-hk K (H entries per suffix)
  int l2r_p3_split[2] = {0, 0};           // --l2r-p3-split Ba,bh
  bool l2r_p_raw = false;                 // --l2r-p-aopen q|raw
  bool l2r_p3_hex = false;                // --l2r-p3-hex on
  bool l2r_p3_lin = false;                // --l2r-p3-lin on
  bool l2r_inc_check = false;             // --l2r-inc-check on
  bool l2r_ninf_mass = false, l2r_keep_ninf = false;  // --l2r-ninf-mass on, --l2r-keep-ninf on
  bool l2r_qz_audit = false, ntr_z_set = false;
  bool l2r_spawn_wk = false;                // --l2r-spawn-w suf|k
  int l2r_spawn_pool = 128, l2r_spawn_maxpop = 0;  // --l2r-spawn-pool, --l2r-spawn-maxpop (0: 4 pool)
  std::string l2r_tbeam;  // --l2r-tbeam NAME=B[,NAME=B..], B = 0 or inf: no beam
  std::string l2r_odump, l2r_oload, l2r_sub, l2r_trace, l2r_trace_key;
  bool l2r_perf = false;
  bool k_memo = false, k_bound = false, k_exact_ub = false, k_exact_check = false, k_exact = false;
  bool rlt = false, next_pair = false, grammar_simple = false, c_trace = false, pass1_only = false;
  int rlt_beam = 0;
  std::string sample_prof;
  bool l2r_ctx = false, l2r_kx = false;  // --l2r-ctx on: span context (needs Bml); --l2r-kx on: exact K
  bool l2r_ctx_p1 = false;               // --l2r-ctx p1: pass-1 context ratio lambda(a')
  double l2r_ctx_tau = 0.1;              // --l2r-ctx-tau
  int l2r_ctx_rounds = 1;                // --l2r-ctx-rounds: pass-1 rounds that recompute lambda
  bool l2r_split = false;                // --l2r-split on: consumer split of the outsides (diag)
  int l2r_p1 = 0;                        // --l2r-p1 R: rounds of the beamed pass-1 sweep
  int l2r_p1_beam = 250;                 // --l2r-p1-beam B
  bool l2r_p1_own = true;                // --l2r-p1-qb pass1|nested
  double pk_db[3] = {0.0, 0.0, 0.0};     // --pk-dbeta1 / -1m / -1p: kcal/mol added after --params
  bool pk_l2r_R = false, l2r_diag = false;
  bool arm_score = false;
  std::string access_file;
  int max_gap_scan = 200;
  int max_arm_len = 0;
  bool gap_step = false;
  bool loop_cap_set = false;
  bool restrict_pairs = false;
  std::vector<std::pair<int, int>> restrict_list;
  GateCfg gcfg;
  bool gate_stats_out = false;
  std::string gate_diag_path;
  std::vector<double> diag_taus = {0.0, 1e-4, 1e-3, 1e-2, 0.05, 0.1};
  int long_span = 100;
  int table_beam[BT_N];
  std::fill(table_beam, table_beam + BT_N, -1);
  std::string states_path, pruned_path, exact_path, scores_path, ext_path, kept_path, beam_stats_path;
  double dump_min_w = 0.0, pruned_min_w = 0.0;
  std::unordered_set<uint64_t> remove_keys;
  bool stem_rank = false, stem_both = false, stem_qg2 = false;
  bool stem_qg_partner = false, stem_consistent = false, stem_qg_r_set = false,
       stem_consts_set = false, stem_bg_check = false;
  int stem_ext = 0, stem_ext_max = 30;
  bool stem_split = false, stem_ext_set = false, stem_ext_check = false;
  int grow_lg = 3, grow_nodes = 2000, ml_span = 30;
  double grow_cut = 1e-8;
  int stem_helpers = 0;
  bool grow_set = false, helpers_set = false;
  bool stem_profile = false, profile_set = false, ml_tail = false, prof_check = false;
  bool bg_reorder = false, ro_check = false;
  double prof_eps = 0.0;
  bool stem_chain = false, chain_set = false, chain_ml = true;
  double chain_floor = 1e-7, chain_cut = 1e-8;
  int chain_span = 150, chain_loop = 30;
  bool stem_hr_av = false, stem_hr_phi = false, stem_hr_set = false, stem_hr_check = false;
  double stem_hr_gamma = 1.0;
  bool stem_pr = false, stem_pr_nested = false, pr_set = false, pr_check = false;
  int pr_span = 300, pr_ml_span = 150, pr_ml_close = 30, pr_arms_k = 64;
  double pr_floor = 1e-7, pr_cut = 1e-8;
  std::string pr_oracle_path, pr_oracle_set = "step1";
  bool pr_ratio = false, pr_gapfree = false, pr_cut_roots = false;
  // --pr-arms exact|sep, --pr-arm-z qz|unpaired, --pr-partner-floor (sep default 1e-20; off:
  // the kept pairs), --pr-test-qz unpaired, --pr-sep-check; sep defaults --pr-cut-by gapfree
  int pr_arms_mode = 0, pr_arm_z = 0;
  bool pr_qz_unp = false, pr_sep_check = false, pr_cut_by_set = false;
  double pr_partner_floor = -1.0;
  double pr_arms_eps = 1e-3;
  int pr_arms_max = 256;
  // --pr-grow-select floor|bound, --pr-grow-eps, --pr-grow-B, --pr-grow-check, --pr-score-dump
  bool pr_gsel = false, pr_gcheck = false;
  bool pr_aq_sep = false;  // --pr-armq-sep on|off
  size_t pr_aq_max = 128;  // --pr-aq-max: Step A entries kept
  bool pr_arm_persist = false, pr_arm_gaponly = false;  // --pr-arm-persist, --pr-arm-gaponly
  int pr_pspan = 0;                                     // --pr-partner-span W
  int pr_pmax = 0;         // --pr-partner-maxpairs K (0 = off)
  bool pr_pcheck = false;  // --pr-partner-check
  bool pr_pml = true, pr_pmlcheck = false;  // --pr-partner-ml, --pr-partner-ml-check
  int pr_fast = 0;                          // --pr-fast off|acc|ml
  int pr_abl = 0;                           // --pr-abl bml,q,p,xs,ks,sep
  int pr_abl_pk = 0;                        // --pr-abl-pk
  double pr_veps = 1e-8;                                // --pr-v-eps (0: no trimming)
  bool pr_vsuffix = false;                              // --pr-v-trim global|suffix
  int pr_vcap = 0;                                      // --bV (was --pr-v-cap); default 256
  int p2_qz_beam = 0, p2_qz_key = 1, p2_v_mode = 0, p2_v_bp = 0, p2_v_cap = 0, p2_v_trim = 0;
  bool bV_set = false;
  int p2_g_mode = 0, p2_g_bp = 0, p2_g_cap = 0, nest_beam = -1, pre_beam = 0, pre_key = 0;
  bool b1_set = false, rlt_beam_set = false, sweep_set = false, vsweep_set = false;
  int b2_width = 0;
  int b2_at = -1, pbeam_explicit_at = -1, pbeam_ratio_at = -1, pbr[3] = {1, 2, 2};
  bool user_pbeam = false, user_pratio = false;
  bool sw_check = false, p2_sweep = false, p2_vsweep = false;
  int p2_v_floor = 0;  // 0 off, 1 smallest kept mu, 2 log-weight 0
  bool p2_v_cut_fb = false;
  bool lazy_ml_cap = false;
  bool p2_g_suffix = false, p2_g_eager = false, p2_g_fb = false, profile = false, dump_qb_on = false;
  bool qzbar_old = false;  // --qzbar-old: dense Qz→prQz like Oct-8 incperf (no sparse Qz_pre)
  int pr_vmlk = 0;                                      // --pr-v-mlk K (0: all)
  bool pr_vmlsteps = false;                             // --pr-v-mlsteps on
  bool pr_vxs = false;                                  // --pr-v-key xsys|xs
  bool pr_vstream = false;                              // --pr-v-stream off|on
  bool l2r_arm_fb = false;                              // --l2r-arm-fb on|off
  // --l2r-lean on (--pk-l2r on with a beam): skip the pre-pass structures only the b16
  // scoring reads (beamed nested run with its stem index, S, C0, nested outside; chain;
  // prKept; prQbIn unless pass 1 or the V diagnostics need it)
  bool l2r_lean = false;
  bool pr_v_ratio = false;  // --pr-v-ratio on|off
  std::string pr_v_diag;    // --pr-v-diag PREFIX
  double pr_v_sample = 0.0; // --pr-v-sample P
  double pr_geps = 1e-6;
  int pr_gB = 512;
  std::string pr_dump_path;
  int stem_loop = 2;
  std::string stem_diag_path;
  std::vector<int> stem_diag_loops = {0, 2, 6};

