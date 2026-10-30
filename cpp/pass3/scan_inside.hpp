// Pass 3, the L2R scan = inside of LR-DP03: advance_pk_inside (PKM), l2r_init, tails, growth
// of outer pairs (l2r_grow), spawns, closings A / B, stem / filler extensions, the P1-P3 steps
// and the incremental seeds (--l2r-spawn inc).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  // PKM[j] = PKM[j-1] + one unpaired, plus Qp ending at j; then beam.
  void advance_pk_inside(int j) {
    if (!allow_pk) return;
    const double w3 = E.w_alpha3;
    OMap& at = PKM_at[j];
    if (dz) {
      RunSt& cu = rsPKM[j & 1];
      cu.clear();
      std::vector<int> surv;
      if (j > 0) {
        surv = live_keys(PKM_at[j - 1]);
        run_ext(rsPKM[(j - 1) & 1], cu, surv, j, w3);
      }
      const OMap& qp = Qp_at[j];
      for (size_t r = 0; r < qp.slots(); ++r)
        if (qp.keys[r] >= 0) cu.P.put(qp.keys[r], qp.vals[r]);
      run_fin(nullptr, at, cu, surv, j);
      prune1(at, j);
      return;
    }
    if (j > 0) {
      const OMap& prev = PKM_at[j - 1];
      for (size_t r = 0; r < prev.slots(); ++r) {
        int d = prev.keys[r];
        if (d >= 0) at.set(d, prev.vals[r] + w3);
      }
    }
    const OMap& qp = Qp_at[j];
    for (size_t r = 0; r < qp.slots(); ++r) {
      int d = qp.keys[r];
      if (d >= 0) at.set(d, log_add(at.get(d), qp.vals[r]));
    }
    prune1(at, j);
  }

  // ============================================ --pk-l2r: H-type pseudoknots in scan order
  // Qp(a', ys) = sum A F1 F2 F3 B (stemest/l2r/r_stage0.md) with the partner's pieces as states:
  //   A (LA): seeds Gs(a2, xs, x', c1) = pb + w_stem + at(a2, c1) + F1 + F2, grown G(a, xs, x', j);
  //     Tab4 group x', entry (p = a, u = xs).
  //   S(a', xs, x', j) = sum_q' [Afull(q') + at(a', q')] Qz(q'+1, j) through the Qz run automaton
  //     (P, U1, U2; c = the close-A term at j; F = closed + c).
  //   B (LB): seeds Gs(a', xs, x', y') = S_F(y'-1) + w_stem + at(x', y') + pb, grown G(a', xs, x, y),
  //     x >= xs; Tab4 group xs, entry (p = x, u = a').  Qp(a', y) += Bfull(xs, y) + at(xs, y).
  //   Gls / Glss (seed-derived) / Grs: Qgls / Qgrs with the carried keys.
  // l2r_mode 1: Qp from the recursion; 2 (check): both, the engine's Qp feeds the nested tables.
  int l2r_mode = 0;
  bool l2r_check_R = false;
  struct LFam {
    bool isU = false;
    std::vector<Tab4> G, Gs, Gls, Glss, Grs, Glt;
    // --grammar simple: M, M2, Multi. M2 groups with v >= n are M5 (stem then nested);
    // v < n are M4 (nested then stem). The logical group key is v % n.
    std::vector<Tab4> SM, SM2, SMulti;
  };
  LFam LA, LB;
  bool grammar_simple = false;
  bool c_trace = false;
  long long gs_fill = 0, gs_phase = 0, gs_lead = 0, gs_over = 0, gs_m[8] = {};
  // Kept after each type's own beam: R^M, R^M2, R^Multi, U^M, U^M2, U^Multi.
  long long gs_tab_pre[6] = {}, gs_tab_kept[6] = {}, gs_tab_max[6] = {};
  struct LeadCol {
    std::vector<double> w;
  };
  std::vector<LeadCol> gsLead;
  std::unordered_set<uint64_t> gs_seedA, gs_seedB;
  struct SSt {
    double P, U1, U2, c, F;
    double bP, bU1, bU2, bc, bF;
  };
  struct SMap {
    std::vector<uint64_t> keys;
    std::vector<SSt> v;
    std::unordered_map<uint64_t, int> idx;
    int at(uint64_t k) {
      auto it = idx.find(k);
      if (it != idx.end()) return it->second;
      idx.emplace(k, static_cast<int>(keys.size()));
      keys.push_back(k);
      v.push_back(SSt{NEG_INF, NEG_INF, NEG_INF, NEG_INF, NEG_INF, NEG_INF, NEG_INF, NEG_INF,
                      NEG_INF, NEG_INF});
      return static_cast<int>(keys.size()) - 1;
    }
    int find(uint64_t k) const {
      auto it = idx.find(k);
      return it == idx.end() ? -1 : it->second;
    }
  };
  std::vector<SMap> LS;
  std::vector<double> LQp;
  std::vector<int> lmaxR;  // largest y with pairs(x, y) and end_ok(x, y), -1 if none
  enum { LT_QZ, LT_QM, LT_QB, LT_QP, LT_AG, LT_AGS, LT_ALS, LT_ARS, LT_BG, LT_BGS, LT_BLS,
         LT_BRS, LT_S, LT_P1, LT_P2, LT_P3, LT_N };
  // --l2r-spawn inc: the A seeds built in scan order from three filler tables. The column
  // is the right end, so the index counts are P1(a, t) = 2, P2(a, b, t) = 3, P3(a, b, c, t) = 4.
  // key3 packs the indices other than t; unused slots stay 0.
  //   P1(a2, j) = Qz(a2+1, j), opened with c = 1 at j = a2, extended as S (fill_Qz automaton);
  //   P2(a2, xs, j) = P1(a2, xs-1) carried (value in F): B's arm xs.. is weighed later by U;
  //   P3(a2, xs, x', j) = P2(a2, xs, x') Qz(x'+1, j), opened with c = P2 at j = x', extended as S;
  //   Rs(a2, xs, x', c1) = P3(a2, xs, x', c1-1) + pb + w_stem + at(a2, c1) at a pair a2.c1.
  bool l2r_inc = false;
  std::vector<SMap> LP1, LP2, LP3;
  int l2r_done[LT_N];
  long long l2r_tr[LT_N] = {};
  // l2r_grow edges that land in R or U: 0..3 R interior/(f)/(g)/(h), 4..7 the same for U.
  int ext_slot = -1;
  long long ext_n[8] = {};
  // Loop iterations, including ones that fail the pair check. [0] R, [1] U.
  long long it_int[2] = {}, it_h[2] = {}, it_ls[2] = {}, it_rs[2] = {};
  // Cumulative counters at the end of each scan column, for the last-20% means.
  struct ColCnt {
    long long iti, ith, itls, itrs, tail, edge, qz, k, qm;
  };
  std::vector<ColCnt> col_cnt;
  bool next_pair = false;
  void note_col() {
    long long edge = 0;
    for (int s = 0; s < 8; ++s) edge += ext_n[s];
    col_cnt.push_back({it_int[0] + it_int[1], it_h[0] + it_h[1], it_ls[0] + it_ls[1],
                       it_rs[0] + it_rs[1], rlt_ext, edge, qz_reads, sx_k_n, qm_kept_n});
  }
  // Qp after the per-position beam (preset beam). pre is the column before prune2.
  long long qp_pre_n = 0, qp_kept_n = 0, qp_bind_n = 0;
  int qp_pre_max = 0, qp_kept_max = 0, qp_beam = 0;
  // Qm keys after the per-column beam.
  long long qm_pre_n = 0, qm_kept_n = 0, qm_bind_n = 0;
  int qm_pre_max = 0, qm_kept_max = 0, qm_beam = 0;
  // --rlt: case (f) right flank kept as a one-base tail. rlt_beam 0 does not prune.
  bool rlt = false;
  int rlt_beam = 0;
  long long rlt_kept = 0;
  int rlt_max = 0;
  long long rlt_ext = 0;
  long long rlt_kept_f[2] = {};
  int rlt_max_f[2] = {};
  std::vector<int> rlt_pos[2];
  // seconds: A grow+helpers, S, B grow+helpers, closes, spawn, checks
  double l2r_sec[6] = {};
  long long lc_R_n = 0, lc_R_mis = 0, lc_Qp_n = 0, lc_Qp_mis = 0;
  double lc_R_rel = 0.0, lc_Qp_rel = 0.0;

  static int k3a(uint64_t k) { return static_cast<int>(k >> 42); }
  static int k3b(uint64_t k) { return static_cast<int>((k >> 21) & 0x1FFFFF); }
  static int k3c(uint64_t k) { return static_cast<int>(k & 0x1FFFFF); }

  void l2r_init() {
    LA.isU = false;
    LB.isU = true;
    for (LFam* F : {&LA, &LB}) {
      F->G.assign(n, Tab4{});
      F->Gs.assign(n, Tab4{});
      F->Gls.assign(n, Tab4{});
      F->Glss.assign(n, Tab4{});
      F->Grs.assign(n, Tab4{});
      if (rlt) F->Glt.assign(n, Tab4{});
      if (grammar_simple) {
        F->SM.assign(n, Tab4{});
        F->SM2.assign(n, Tab4{});
        F->SMulti.assign(n, Tab4{});
      }
    }
    if (grammar_simple) gsLead.assign(n, LeadCol{});
    LS.assign(n, SMap{});
    if (l2r_inc) {
      LP1.assign(n, SMap{});
      LP2.assign(n, SMap{});
      LP3.assign(n, SMap{});
    }
    if (l2r_mode == 2) LQp.assign(static_cast<size_t>(n) * n, NEG_INF);
    ensure_pairs();
    lmaxR.assign(n, -1);
    for (int x = 0; x < n; ++x)
      for (int y : partners_right[x])
        if (E.end_ok(x, y)) lmaxR[x] = std::max(lmaxR[x], y);
    std::fill(l2r_done, l2r_done + LT_N, -1);
    qz_audit = qz_audit_req;
  }

  void lchk(int t, int pos) const {
    if (pos <= l2r_done[t]) return;
    static const char* nm[LT_N] = {"Qz", "Qm", "Qb", "Qp", "R", "Rs", "Rls", "Rrs",
                                   "U",  "Us", "Uls", "Urs", "S", "P1", "P2", "P3"};
    std::fprintf(stderr, "[l2r] abort: read of %s at %d, finished through %d (scan %d)\n", nm[t],
                 pos, l2r_done[t], cur_j);
    std::abort();
  }
  // B with left arm xs..x' can still close right of j
  bool l2r_arm(int xs, int xp, int j) const { return lmaxR[xs] > j && lmaxR[xp] > j; }
  void lacc(Tab4& t, int lt, int v, int p, int u, double val) {
    if (val < HALF_NEG) return;
    ++l2r_tr[lt];
    if (ext_slot >= 0) ++ext_n[ext_slot];
    t.acc(v, p, u, val);
  }
  // Every transition of the grow / helper tables is target (v, p, u) += s.a + w (+ Qm at m
  // when qm_i >= 0).  l2r_adj: the same enumeration pushes the target's outside value back.
  bool l2r_adj = false;
  void ledge(Tab4& T, int lt, int v, int p, int u, Ent& s, double w, int qm_i = -1,
             int qm_j = -1, double qm = 0.0) {
    if (!l2r_adj) {
      lacc(T, lt, v, p, u, s.a + (w + qm));
      return;
    }
    Ent* t = T.find(v, p, u);
    if (!t || t->b < HALF_NEG) return;
    add4(s, t->b + (w + qm));
    if (qm_i >= 0) add2b(qmp_b(), qm_i, qm_j, t->b + s.a + w);
  }

  // Tail key packs the run's left end s with the carried u (xs for R, a' for U).
  int tail_pack(int s, int u) const { return s * n + u; }
  void tail_unpack(int packed, int& s, int& u) const {
    s = packed / n;
    u = packed - s * n;
  }
  // Log weight added when the unpaired flank's right end becomes `right` (s == right: a new
  // one-base run). Dangles some: one base is s1, a longer run is d5(left) + d3(right).
  double tail_extend(int s, int right) const {
    if (!dz) return E.w_cp;
    if (s == right) return E.w_cp + E.s1w[right];
    if (right - s == 1) return E.w_cp + (E.d5w[s] + E.d3w[right]) - E.s1w[s];
    return E.w_cp + E.d3w[right] - E.d3w[right - 1];
  }
  void l2r_prune_by_inside(Tab4& t, int beam) {
    size_t tot = 0;
    for (const Grp& g : t.grps) tot += g.ents.size();
    if (beam <= 0 || tot <= static_cast<size_t>(beam)) return;
    std::vector<double> sc;
    sc.reserve(tot);
    for (const Grp& g : t.grps)
      for (const Ent& e : g.ents) sc.push_back(e.a);
    std::vector<char> keep;
    keep_mask(sc, keep, beam);
    filter_step(t, keep);
  }
  // Rank by C(p-1) + inside. Keys only; the stored inside is unchanged.
  void l2r_prune_prefix(Tab4& t, int b, int slot = -1) {
    size_t tot = 0;
    for (const Grp& g : t.grps) tot += g.ents.size();
    size_t kept = tot;
    if (b > 0 && tot > static_cast<size_t>(b)) {
      std::vector<double> sc;
      sc.reserve(tot);
      for (const Grp& g : t.grps)
        for (const Ent& e : g.ents) sc.push_back(prefix_log(e.p) + e.a);
      std::vector<char> keep;
      keep_mask(sc, keep, b);
      filter_step(t, keep);
      kept = 0;
      for (const Grp& g : t.grps) kept += g.ents.size();
    }
    if (slot >= 0) {
      gs_tab_pre[slot] += static_cast<long long>(tot);
      gs_tab_kept[slot] += static_cast<long long>(kept);
      if (static_cast<long long>(kept) > gs_tab_max[slot])
        gs_tab_max[slot] = static_cast<long long>(kept);
    }
  }
  // Same score as Rls / Uls (kind 1, not a seed). The packed run start is not part of the key.
  void l2r_prune_tail(LFam& F, Tab4& t, int j) {
    size_t tot = 0;
    for (const Grp& g : t.grps) tot += g.ents.size();
    if (rlt_beam <= 0 || tot <= static_cast<size_t>(rlt_beam)) return;
    std::vector<double> sc;
    sc.reserve(tot);
    for (const Grp& g : t.grps)
      for (const Ent& e : g.ents) {
        int run_s, uu;
        tail_unpack(e.u, run_s, uu);
        (void)run_s;
        if (!F.isU) {
          const RN& r = rec_key(true, e.p, false, j);
          const int arm = larm(uu, g.v);
          const double s = lsc_L(r, arm, 1);
          sc.push_back(s > HALF_NEG ? e.a + s + lkb_arm(e.p, arm) : NEG_INF);
        } else {
          double s;
          if (g_rfv) s = rfv_key(true, e.p, false, j, g.v);
          else s = lsc_R(rec_key(true, e.p, false, j), g.v);
          const double add = prefix_log(uu) + (lsub_c ? lrho(uu) : lkb_xs(uu, g.v));
          sc.push_back(s > HALF_NEG ? e.a + s + add + (lLam.empty() ? 0.0 : lLam[uu]) : NEG_INF);
        }
      }
    std::vector<char> keep;
    keep_mask(sc, keep, rlt_beam);
    filter_step(t, keep);
  }
  // Rlt/Ult at j: Rls/Uls at j-1 plus one unpaired base, and every previous tail plus one base.
  void l2r_tail(LFam& F, int j) {
    if (!rlt || j < 1) return;
    Tab4& dst = F.Glt[j];
    const double add_new = tail_extend(j, j);
    for (const Grp& g : F.Gls[j - 1].grps)
      for (const Ent& e : g.ents) {
        ++rlt_ext;
        dst.acc(g.v, e.p, tail_pack(j, e.u), e.a + add_new);
      }
    for (const Grp& g : F.Glt[j - 1].grps)
      for (const Ent& e : g.ents) {
        ++rlt_ext;
        int s, uu;
        tail_unpack(e.u, s, uu);
        dst.acc(g.v, e.p, e.u, e.a + tail_extend(s, j));
      }
    if (rlt_beam > 0) l2r_prune_tail(F, dst, j);
    size_t tot = 0;
    for (const Grp& g : dst.grps) tot += g.ents.size();
    const int fam = F.isU ? 1 : 0;
    if (rlt_pos[fam].size() != static_cast<size_t>(n)) rlt_pos[fam].assign(n, 0);
    rlt_pos[fam][j] = static_cast<int>(tot);
    rlt_kept_f[fam] += static_cast<long long>(tot);
    rlt_max_f[fam] = std::max(rlt_max_f[fam], static_cast<int>(tot));
    rlt_kept += static_cast<long long>(tot);
    rlt_max = std::max(rlt_max, static_cast<int>(tot));
  }
  // Adjoint of l2r_tail(j): Glt[j] back into Gls[j-1] (new one-base runs) and Glt[j-1].
  void l2r_tail_back(LFam& F, int j) {
    if (!rlt || j < 1) return;
    const double add_new = tail_extend(j, j);
    for (Grp& g : F.Glt[j].grps)
      for (Ent& e : g.ents) {
        if (e.b < HALF_NEG) continue;
        int s, uu;
        tail_unpack(e.u, s, uu);
        Tab4& src_tab = s == j ? F.Gls[j - 1] : F.Glt[j - 1];
        Grp* src = src_tab.get(g.v);
        if (!src) continue;
        const uint64_t key = s == j ? key2(e.p, uu) : key2(e.p, e.u);
        auto it = src->idx.find(key);
        if (it == src->idx.end()) continue;
        add4(src->ents[it->second], e.b + (s == j ? add_new : tail_extend(s, j)));
      }
  }
  // Case (f) from the tail at j-1, the empty run on Rls/Rlss at j-1, and Rlss for runs 1..cap.
  // Rlss is not copied into the tail, so a seed derivation is closed once.
  void l2r_close_tails(LFam& F, int j, int tg) {
    const int o = F.isU ? 4 : 0;
    ext_slot = (F.isU ? 4 : 0) + 1;
    auto close = [&](Grp& g, Ent& s, int uu, double w_pre) {
      const int i = s.p - 1;
      const int lo = F.isU ? g.v : 0;
      if (i < lo || !E.pairs(i, j) || !E.end_ok(i, j)) return;
      if (!F.isU && !l2r_arm(uu, g.v, j)) return;
      ledge(F.G[j], tg, g.v, i, uu, s, w_pre + E.at(i, j));
    };
    if (j < 1) return;
    lchk(LT_ALS + o, j - 1);
    const double w_empty = (E.w_ap + E.w_bp) + pair_bonus;
    for (Grp& g : F.Glt[j - 1].grps)
      for (Ent& s : g.ents) {
        int run_s, uu;
        tail_unpack(s.u, run_s, uu);
        (void)run_s;
        close(g, s, uu, w_empty);
      }
    for (int src = 0; src < 2; ++src)
      for (Grp& g : (src ? F.Glss[j - 1] : F.Gls[j - 1]).grps)
        for (Ent& s : g.ents) close(g, s, s.u, w_empty);
    const int flo = std::max(0, j - 1 - cap);
    for (int f = flo; f < j - 1; ++f) {
      lchk(LT_ALS + o, f);
      const int unpaired = j - f - 1;
      double w_ext = E.w_ap + E.w_bp + E.w_cp * unpaired + pair_bonus;
      if (dz) w_ext += E.run(f + 1, j - 1);
      for (Grp& g : F.Glss[f].grps)
        for (Ent& s : g.ents) close(g, s, s.u, w_ext);
    }
  }
  // New outer pairs i.j of family F (fill_Qg (b), (f), (g), (h); (a) = (b) on the seeds Gs).
  // B: every new left end >= the reserved xs (group key).
  void l2r_grow(LFam& F, int j) {
    const int o = F.isU ? 4 : 0;
    Tab4& G = F.G[j];
    const int tg = LT_AG + o;
    auto alive = [&](const Grp& g, const Ent& s) { return F.isU || l2r_arm(s.u, g.v, j); };
    const int base = F.isU ? 4 : 0;
    ext_slot = base;
    for (int f = std::max(0, j - 1 - cap); f < j; ++f) {
      lchk(LT_AG + o, f);
      lchk(LT_AGS + o, f);
      for (int src = 0; src < 2; ++src)
        for (Grp& g : (src ? F.Gs[f] : F.G[f]).grps) {
          const int lo = F.isU ? g.v : 0;
          for (Ent& s : g.ents) {
            if (!alive(g, s)) continue;
            const int c = s.p;
            const int iu = F.isU ? 1 : 0;
            const int i0 = std::max(lo, c - 1 - cap);
            auto take = [&](int i) {
              const double w = E.stem_interior(i, c, f, j);
              if (w < HALF_NEG) return;
              ledge(G, tg, g.v, i, s.u, s, w + pair_bonus);
            };
            if (next_pair) {
              const std::vector<int>& pl = partners_left[static_cast<size_t>(j)];
              for (auto it = std::lower_bound(pl.begin(), pl.end(), i0); it != pl.end() && *it < c;
                   ++it) {
                if (!l2r_adj) ++it_int[iu];
                take(*it);
              }
            } else {
              for (int i = i0; i < c; ++i) {
                if (!l2r_adj) ++it_int[iu];
                if (!E.pairs(i, j)) continue;
                take(i);
              }
            }
        }
      }
    }
    ext_slot = base + 1;
    // --grammar simple: the multiloop is M1..M8, not cases (f)(g)(h) and not the tail.
    if (grammar_simple) {
      ext_slot = -1;
      return;
    }
    if (rlt) {
      l2r_close_tails(F, j, tg);
    } else {
      for (int f = 0; f < j; ++f) {
        lchk(LT_ALS + o, f);
        const int unpaired = j - f - 1;
        double w_ext = E.w_ap + E.w_bp + E.w_cp * unpaired + pair_bonus;
        if (dz) w_ext += E.run(f + 1, j - 1);
        for (int src = 0; src < 2; ++src) {
          if (src == 1 && unpaired > cap) continue;
          for (Grp& g : (src ? F.Glss[f] : F.Gls[f]).grps) {
            const int lo = F.isU ? g.v : 0;
            for (Ent& s : g.ents) {
              const int i = s.p - 1;
              if (i < lo || !alive(g, s) || !E.pairs(i, j) || !E.end_ok(i, j)) continue;
              ledge(G, tg, g.v, i, s.u, s, w_ext + E.at(i, j));
            }
          }
        }
      }
    }
    if (j < 1) {
      ext_slot = -1;
      return;
    }
    lchk(LT_ARS + o, j - 1);
    for (Grp& g : F.Grs[j - 1].grps) {
      const int lo = F.isU ? g.v : 0;
      for (Ent& s : g.ents) {
        if (!alive(g, s)) continue;
        const int c = s.p;
        ext_slot = base + 2;
        for (int i = std::max(lo, lead_i_lo(c)); i < c; ++i) {
          if (!E.pairs(i, j) || !E.end_ok(i, j)) continue;
          double w = E.w_ap + E.w_bp + E.w_cp * (c - i - 1) + E.at(i, j) + pair_bonus;
          if (dz) w += E.run(i + 1, c - 1);
          ledge(G, tg, g.v, i, s.u, s, w);
        }
        if (c - 1 < 0) continue;
        lchk(LT_QM, c - 1);
        ext_slot = base + 3;
        const OMap& at = qmp_at(c - 1);
        const int iu = F.isU ? 1 : 0;
        qm_pair_slots(c - 1, j, [&](int r) {
          if (!l2r_adj) ++it_h[iu];
          const int ii = at.keys[r];
          if (ii < 0) return;
          const int outer = ii - 1;
          if (outer < lo || outer >= c || !E.pairs(outer, j) || !E.end_ok(outer, j)) return;
          ledge(G, tg, g.v, outer, s.u, s, E.w_ap + E.w_bp + E.at(outer, j) + pair_bonus, ii,
                c - 1, at.vals[r]);
        });
      }
    }
    ext_slot = -1;
  }

  // Gls / Glss at j: Qm(i, c-1) left of a branch c.j (fill_Qgls; seeds from Gs).
  void l2r_ls(LFam& F, int j) {
    const int o = F.isU ? 4 : 0;
    for (int src = 0; src < 2; ++src) {
      Tab4& out = src ? F.Glss[j] : F.Gls[j];
      for (Grp& g : (src ? F.Gs[j] : F.G[j]).grps) {
        const int lo = F.isU ? g.v + 1 : 1;
        for (Ent& s : g.ents) {
          const int c = s.p;
          if (c - 1 < 0 || !E.end_ok(c, j)) continue;
          lchk(LT_QM, c - 1);
          const double w = E.w_bp + E.at(c, j);
          const OMap& at = qmp_at(c - 1);
          const int iu = F.isU ? 1 : 0;
          for (size_t r = 0; r < at.slots(); ++r) {
            if (!l2r_adj) ++it_ls[iu];
            const int i = at.keys[r];
            if (i < lo) continue;
            ledge(out, LT_ALS + o, g.v, i, s.u, s, w, i, c - 1, at.vals[r]);
          }
        }
      }
    }
  }

  // Grs at the end of step j: a branch c.f, f < j, then Qm(f+1, j) (fill_Qgrs).
  void l2r_rs(LFam& F, int j) {
    const int o = F.isU ? 4 : 0;
    lchk(LT_QM, j);
    Tab4& out = F.Grs[j];
    const OMap& at = qmp_at(j);
    const int iu = F.isU ? 1 : 0;
    for (size_t r = 0; r < at.slots(); ++r) {
      if (!l2r_adj) ++it_rs[iu];
      const int p = at.keys[r];
      if (p < 1) continue;
      const int f = p - 1;
      lchk(LT_AG + o, f);
      lchk(LT_AGS + o, f);
      for (int src = 0; src < 2; ++src)
        for (Grp& g : (src ? F.Gs[f] : F.G[f]).grps)
          for (Ent& s : g.ents) {
            if (!E.end_ok(s.p, f) || (!F.isU && !l2r_arm(s.u, g.v, j))) continue;
            ledge(out, LT_ARS + o, g.v, s.p, s.u, s, E.w_bp + E.at(s.p, f), p, j, at.vals[r]);
          }
    }
  }

  // Outside of the A seeds: Rs(a2, xs, x', j) -> Qz(a2+1, xs-1), Qz(x'+1, j-1).
  void l2r_spawn_back(int j) {
    for (Grp& g : LA.Gs[j].grps) {
      const int xp = g.v;
      for (Ent& s : g.ents) {
        if (s.b < HALF_NEG) continue;
        const int a2 = s.p, xs = s.u;
        const double k0 = pair_bonus + E.w_stem + E.at(a2, j);
        const double f1 = qz_log(a2 + 1, xs - 1), f2 = qz_log(xp + 1, j - 1);
        if (xs - 1 >= a2 + 1) add2b(Qz_b, a2 + 1, xs - 1, s.b + k0 + f2);
        if (j - 1 >= xp + 1) add2b(Qz_b, xp + 1, j - 1, s.b + k0 + f1);
      }
    }
  }

  // --l2r-spawn cube (Huang & Chiang 2007 cube pruning) at the inner pair a2.c1 (c1 = j).
  // Rows: xs > a2 by Qzbar(a2+1, xs-1) + Wmax(xs) (pre-pass Qzbar: the scan's Qz right of j is
  // not filled yet), built on a2's first spawn and cached.  Columns: xs's V arms by the static
  // weight log sum_{y'} E_in V (lYwS at the arm's first y').  Priority of cell (i, c): exact
  // F1(a2, xs_i) + F2(x'_ic, c1) + W_c1(arm), W_c1 = lS(arm, c1) = log sum_{y' > c1} E_in V
  // (--l2r-spawn-w suf, default) or lK(arm, c1), the same sum with the gap Qz(c1+1, y'-1) (k):
  // R's K at q' = c1.  Pop the best; push (i, c+1),
  // and the next row's (i', 0) when c == 0; invalid cells (x' >= c1, no y' > c1, a B that
  // cannot close right of c1, F = 0) are counted as pops and expanded at once.  Stop after
  // sp_pool valid pops or sp_maxpop pops; the valid ones re-ranked by the Rs score (alpha +
  // lsc_L + context, as l2r_pruneA) and the top l2r_spawn_k spawned with value k0 F1 F2.
  // --l2r-spawn cube (default form): per step c1 = j every arm valid at c1 gets G = F2(x', c1) +
  // W_c1 (a2-free), columns are xs's arms sorted by G, and per pair every row head (exact
  // F1(a2, xs) + best G of xs) enters one heap: the axes are exactly monotone, so the pops are
  // the exact top sp_pool by F1 + F2 + W_c1.  cube-static: the cached rows / static columns
  // above (sp_static).
  bool sp_cube = false, sp_wk = false, sp_static = false;
  int sp_pool = 128, sp_maxpop = 0;
  std::vector<std::vector<std::pair<double, int>>> spG;
  std::vector<int> spGx;
  int spG_j = -1;
  size_t sp_col_ents = 0, sp_col_max = 0;
  void l2r_spawn_step_cols(int j) {
    SecTimer t(sp_col_sec);
    if (spG.empty()) spG.assign(n, {});
    for (int xs : spGx) spG[xs].clear();
    spGx.clear();
    size_t tot = 0;
    for (int r = 0; r < prXsFirst[j]; ++r) {
      const int xs = prArmXs[r], xp = prArmXp[r];
      if (xp >= j || lmaxR[xs] <= j || lmaxR[xp] <= j) continue;
      const double f2 = qz_log(xp + 1, j - 1);
      if (f2 < HALF_NEG) continue;
      const double w = sp_wk ? lK(r, j) : lS(r, j);
      if (w < HALF_NEG) continue;
      if (spG[xs].empty()) spGx.push_back(xs);
      spG[xs].emplace_back(f2 + w, r);
      ++tot;
    }
    for (int xs : spGx)
      std::sort(spG[xs].begin(), spG[xs].end(),
                [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                  return a.first > b.first || (a.first == b.first && a.second < b.second);
                });
    std::sort(spGx.begin(), spGx.end());
    sp_col_ents += tot;
    sp_col_max = std::max(sp_col_max, tot);
    spG_j = j;
  }
  void l2r_spawn_cube_step(int a2, int j, double k0) {
    if (spG_j != j) l2r_spawn_step_cols(j);
    struct Cell {
      double pr;
      int xs, c;
      bool operator<(const Cell& o) const { return pr < o.pr; }
    };
    std::vector<Cell> heap;
    std::vector<std::pair<double, int>> pool;  // (F1 + F2, arm)
    long long pops = 0;
    {
      SecTimer t(sp_heap_sec);
      for (auto it = std::upper_bound(spGx.begin(), spGx.end(), a2); it != spGx.end(); ++it) {
        const int xs = *it;
        const double f1 = qz_log(a2 + 1, xs - 1);
        if (f1 > HALF_NEG) heap.push_back({f1 + spG[xs][0].first, xs, 0});
      }
      std::make_heap(heap.begin(), heap.end());
      const int maxpop = sp_maxpop > 0 ? sp_maxpop : 4 * sp_pool;
      while (!heap.empty() && static_cast<int>(pool.size()) < sp_pool && pops < maxpop) {
        std::pop_heap(heap.begin(), heap.end());
        const Cell x = heap.back();
        heap.pop_back();
        ++pops;
        const auto& col = spG[x.xs];
        const double f1 = x.pr - col[x.c].first;
        const int r = col[x.c].second;
        pool.emplace_back(f1 + qz_log(prArmXp[r] + 1, j - 1), r);
        if (x.c + 1 < static_cast<int>(col.size())) {
          heap.push_back({f1 + col[x.c + 1].first, x.xs, x.c + 1});
          std::push_heap(heap.begin(), heap.end());
        }
      }
    }
    sp_pops += pops;
    sp_cand += static_cast<long long>(pool.size());
    if (ltr_f && a2 == ltr_a2 && j == ltr_c1) ltr_pool(a2, j, pool, pops, heap.size());
    l2r_spawn_rank(a2, j, k0, pool);
  }
  int ltr_a2 = -1, ltr_c1 = -1;
  void ltr_pool(int a2, int j, const std::vector<std::pair<double, int>>& pool, long long pops, size_t left) {
    const int rt = larm(ltr_xs, ltr_xp);
    const double f1 = qz_log(a2 + 1, ltr_xs - 1), f2 = qz_log(ltr_xp + 1, j - 1);
    const double w = rt >= 0 ? (sp_wk ? lK(rt, j) : lS(rt, j)) : NEG_INF;
    int col = -1;
    for (size_t c = 0; c < spG[ltr_xs].size(); ++c)
      if (spG[ltr_xs][c].second == rt) col = static_cast<int>(c);
    const double pt = f1 + f2 + w;
    int in = -1, rank = 1;
    double lo = INFINITY;
    for (size_t q = 0; q < pool.size(); ++q) {
      const int r = pool[q].second;
      const double p = pool[q].first + (sp_wk ? lK(r, j) : lS(r, j));
      lo = std::min(lo, p);
      rank += p > pt;
      if (r == rt) in = static_cast<int>(q);
    }
    std::fprintf(ltr_f, "spawn-pool a2 %d c1 %d arm %d f1 %.6f f2 %.6f w %.6f pr %.6f col %d/%zu in_pool %d "
                 "rank %d pool %zu pool_min %.6f pops %lld heap_left %zu sp_pool %d\n",
                 a2, j, rt, f1, f2, w, pt, col, spG[ltr_xs].size(), in, rank, pool.size(), lo, pops, left,
                 sp_pool);
  }
  void ltr_rank(int a2, int j, double k0, const RN& rn, const std::vector<std::pair<double, int>>& pool,
                const std::vector<std::pair<double, int>>& rk) {
    const int rt = larm(ltr_xs, ltr_xp);
    const double f1 = qz_log(a2 + 1, ltr_xs - 1), f2 = qz_log(ltr_xp + 1, j - 1);
    const double a = k0 + f1 + f2, s = lsc_L(rn, rt, 0), kb = lkb_arm(a2, rt);
    const double tot = a + s + kb;
    int rank = 1;
    std::vector<double> v;
    for (const auto& e : rk) {
      rank += e.first + k0 > tot;
      v.push_back(e.first + k0);
    }
    const size_t k = std::min(rk.size(), static_cast<size_t>(l2r_spawn_k));
    double cut = NEG_INF;
    if (k > 0) {
      std::nth_element(v.begin(), v.begin() + (k - 1), v.end(), std::greater<double>());
      cut = v[k - 1];
    }
    double full = NEG_INF;
    for (int m = 0; m <= g_ml_L && m < rnL; ++m) {
      const RV& x = rn.L[m];
      for (size_t i = 0; i < x.v.size(); ++i)
        if (x.v[i] > 0.0) {
          const double kk = lK(rt, x.lo + static_cast<int>(i));
          if (kk > HALF_NEG) full = log_add(full, x.mx + std::log(x.v[i]) + kk);
        }
    }
    std::fprintf(ltr_f, "spawn-rank a2 %d c1 %d arm %d a %.6f lsc_L %.6f ctx %.6f score %.6f rank %d of %zu k %zu "
                 "cut %.6f lsc_L_untrunc %.6f pool %zu\n",
                 a2, j, rt, a, s, kb, tot, rank, rk.size(), k, cut, full, pool.size());
    for (const auto& [qp, g] : rn.T[0])
      std::fprintf(ltr_f, "spawn-T0 qp %d lam %.6f K %.6f Kx %.6f prefix_a? zbar0q %.6f\n", qp, g, lK(rt, qp),
                   lKx(rt, qp), stem_->Qzbar(0, qp));
    if (rt >= 0)
      for (int k = prArmY0[rt]; k < prArmY0[rt + 1]; ++k)
        std::fprintf(ltr_f, "spawn-V y %d Yw %.6f zbar0y %.6f\n", prYy[k], prYw[k], stem_->Qzbar(0, prYy[k] - 1));
    for (int m = 0; m <= g_ml_L && m < rnL; ++m) {
      const RV& x = rn.L[m];
      for (size_t i = 0; i < x.v.size(); ++i)
        if (x.v[i] > 0.0)
          std::fprintf(ltr_f, "spawn-L m %d qp %d lam %.6f K %.6f\n", m, x.lo + static_cast<int>(i),
                       x.mx + std::log(x.v[i]), lK(rt, x.lo + static_cast<int>(i)));
    }
  }
  void l2r_spawn_rank(int a2, int j, double k0, const std::vector<std::pair<double, int>>& pool) {
    if (pool.empty()) return;
    ProfReg reg(*this, 2);
    std::vector<std::pair<double, int>> rk;
    {
      SecTimer t(sp_rank_sec);
      const RN& rn = rec_root(a2, j, true);
      for (size_t q = 0; q < pool.size(); ++q) {
        const int r = pool[q].second;
        if (lsub_e) {
          const double s = lors(prArmXs[r], prArmXp[r], a2, j);
          if (s > HALF_NEG) rk.emplace_back(pool[q].first + s, static_cast<int>(q));
          continue;
        }
        sx_ksrc = 3;
        const double s = lsc_L(rn, r, 0);
        sx_ksrc = 0;
        if (s > HALF_NEG) rk.emplace_back(pool[q].first + s + lkb_arm(a2, r), static_cast<int>(q));
      }
      if (ltr_f && a2 == ltr_a2 && j == ltr_c1) ltr_rank(a2, j, k0, rn, pool, rk);
    }
    const size_t k = std::min(rk.size(), static_cast<size_t>(l2r_spawn_k));
    std::partial_sort(rk.begin(), rk.begin() + static_cast<long>(k), rk.end(),
                      [](const std::pair<double, int>& x, const std::pair<double, int>& y) {
                        return x.first > y.first || (x.first == y.first && x.second < y.second);
                      });
    sp_kept += static_cast<long long>(k);
    for (size_t c = 0; c < k; ++c) {
      const auto& p = pool[rk[c].second];
      lacc(LA.Gs[j], LT_AGS, prArmXp[p.second], a2, prArmXs[p.second], k0 + p.first);
    }
  }
  long long sp_pops = 0;
  double sp_row_sec = 0.0, sp_col_sec = 0.0, sp_heap_sec = 0.0, sp_rank_sec = 0.0;
  std::vector<std::vector<int>> spCol, spRow;
  std::vector<char> spRowDone;
  std::vector<double> spWmax;
  size_t sp_row_ents = 0;
  void l2r_spawn_cols() {
    SecTimer t(sp_col_sec);
    spCol.assign(n, {});
    spWmax.assign(n, NEG_INF);
    spRow.assign(n, {});
    spRowDone.assign(n, 0);
    std::vector<double> w(prA, NEG_INF);
    for (int r = 0; r < prA; ++r) {
      if (prArmY0[r] < prArmY0[r + 1]) w[r] = lYwS[prArmY0[r]];
      if (w[r] < HALF_NEG) continue;
      spCol[prArmXs[r]].push_back(r);
      spWmax[prArmXs[r]] = std::max(spWmax[prArmXs[r]], w[r]);
    }
    for (auto& c : spCol)
      std::stable_sort(c.begin(), c.end(), [&](int a, int b) { return w[a] > w[b]; });
  }
  const std::vector<int>& l2r_spawn_row(int a2) {
    if (spRowDone[a2]) return spRow[a2];
    SecTimer t(sp_row_sec);
    spRowDone[a2] = 1;
    std::vector<std::pair<double, int>> v;
    for (int xs = a2 + 1; xs < n; ++xs) {
      if (spCol[xs].empty()) continue;
      const double z = stem_->Qzbar(a2 + 1, xs - 1);
      if (z > HALF_NEG) v.emplace_back(z + spWmax[xs], xs);
    }
    std::stable_sort(v.begin(), v.end(), [](const std::pair<double, int>& x,
                                           const std::pair<double, int>& y) { return x.first > y.first; });
    std::vector<int>& row = spRow[a2];
    row.reserve(v.size());
    for (const auto& e : v) row.push_back(e.second);
    sp_row_ents += row.size();
    return row;
  }
  void l2r_spawn_cube(int a2, int j, double k0) {
    const std::vector<int>& row = l2r_spawn_row(a2);
    const int R = static_cast<int>(row.size());
    const int maxpop = sp_maxpop > 0 ? sp_maxpop : 4 * sp_pool;
    std::vector<double> f1(R, 1.0);  // 1.0: not yet computed
    auto next_row = [&](int i) {
      for (; i < R; ++i) {
        const int xs = row[i];
        if (xs < j && lmaxR[xs] > j) return i;
      }
      return -1;
    };
    struct Cell {
      double pr;
      int i, c;
      bool operator<(const Cell& o) const { return pr < o.pr; }
    };
    std::priority_queue<Cell> heap;
    std::vector<std::pair<int, int>> todo;
    std::vector<std::pair<double, int>> pool;  // (F1 + F2, arm)
    long long pops = 0;
    auto offer = [&](int i0, int c0) {
      todo.emplace_back(i0, c0);
      while (!todo.empty()) {
        const auto [i, c] = todo.back();
        todo.pop_back();
        if (i < 0 || c >= static_cast<int>(spCol[row[i]].size())) continue;
        const int xs = row[i], r = spCol[xs][c], xp = prArmXp[r];
        if (f1[i] > 0.5) f1[i] = qz_log(a2 + 1, xs - 1);
        double pr = NEG_INF;
        if (xp < j && lmaxR[xp] > j && f1[i] > HALF_NEG) {
          const double f2 = qz_log(xp + 1, j - 1), w = sp_wk ? lK(r, j) : lS(r, j);
          if (f2 > HALF_NEG && w > HALF_NEG) pr = f1[i] + f2 + w;
        }
        if (pr > HALF_NEG) {
          heap.push({pr, i, c});
          continue;
        }
        if (++pops >= maxpop) return;
        todo.emplace_back(i, c + 1);
        if (c == 0) todo.emplace_back(next_row(i + 1), 0);
      }
    };
    {
      SecTimer t(sp_heap_sec);
      offer(next_row(0), 0);
      while (!heap.empty() && static_cast<int>(pool.size()) < sp_pool && pops < maxpop) {
        const Cell x = heap.top();
        heap.pop();
        ++pops;
        const int r = spCol[row[x.i]][x.c];
        pool.emplace_back(f1[x.i] + qz_log(prArmXp[r] + 1, j - 1), r);
        offer(x.i, x.c + 1);
        if (x.c == 0) offer(next_row(x.i + 1), 0);
      }
    }
    sp_pops += pops;
    sp_cand += static_cast<long long>(pool.size());
    l2r_spawn_rank(a2, j, k0, pool);
  }

  // A seeds at the inner pair a2.j: every arm a2 < xs <= x' < j whose B can close right of j.
  void l2r_spawn(int j) {
    std::vector<std::pair<double, int>> cand;
    for (int a2 : partners_left[j]) {
      if (!E.end_ok(a2, j)) continue;
      const double k0 = pair_bonus + E.w_stem + E.at(a2, j);
      ++sp_events;
      if (sp_cube) {
        lchk(LT_QZ, j - 1);
        if (sp_static) l2r_spawn_cube(a2, j, k0);
        else l2r_spawn_cube_step(a2, j, k0);
        continue;
      }
      if (l2r_spawn_k > 0) {
        cand.clear();
        lchk(LT_QZ, j - 1);
        for (int r = prXsFirst[a2 + 1]; r < prXsFirst[j]; ++r) {
          const int xs = prArmXs[r], xp = prArmXp[r];
          if (xp >= j || lmaxR[xs] <= j || lmaxR[xp] <= j) continue;
          const double f1 = qz_log(a2 + 1, xs - 1), f2 = qz_log(xp + 1, j - 1);
          if (f1 < HALF_NEG || f2 < HALF_NEG) continue;
          const double s = lS(r, j);
          if (s < HALF_NEG) continue;
          cand.emplace_back(f1 + s + f2, r);
        }
        sp_cand += static_cast<long long>(cand.size());
        const size_t k = std::min(cand.size(), static_cast<size_t>(l2r_spawn_k));
        std::partial_sort(cand.begin(), cand.begin() + static_cast<long>(k), cand.end(),
                          [](const std::pair<double, int>& x, const std::pair<double, int>& y) {
                            return x.first > y.first || (x.first == y.first && x.second < y.second);
                          });
        sp_kept += static_cast<long long>(k);
        for (size_t c = 0; c < k; ++c) {
          const int r = cand[c].second, xs = prArmXs[r], xp = prArmXp[r];
          lacc(LA.Gs[j], LT_AGS, xp, a2, xs, k0 + qz_log(a2 + 1, xs - 1) + qz_log(xp + 1, j - 1));
        }
        continue;
      }
      const long long before = l2r_tr[LT_AGS];
      for (int xs = a2 + 1; xs < j; ++xs) {
        if (lmaxR[xs] <= j) continue;
        lchk(LT_QZ, xs - 1);
        const double f1 = qz_log(a2 + 1, xs - 1);
        if (f1 < HALF_NEG) continue;
        for (int xp = xs; xp < j; ++xp) {
          if (lmaxR[xp] <= j) continue;
          lchk(LT_QZ, j - 1);
          const double f2 = qz_log(xp + 1, j - 1);
          if (f2 < HALF_NEG) continue;
          lacc(LA.Gs[j], LT_AGS, xp, a2, xs, k0 + f1 + f2);
        }
      }
      sp_cand += l2r_tr[LT_AGS] - before;
      sp_kept += l2r_tr[LT_AGS] - before;
    }
  }

  void l2r_closeA(int j) {
    SMap& S = LS[j];
    for (int src = 0; src < 2; ++src) {
      if (src == 1 && !single_stems) continue;
      for (const Grp& g : (src ? LA.Gs[j] : LA.G[j]).grps)
        for (const Ent& s : g.ents) {
          if (!E.end_ok(s.p, j)) continue;
          if (grammar_simple && !single_stems && gs_seedA.count(key3(s.p, s.u, g.v))) continue;
          SSt& st = S.v[S.at(key3(s.p, s.u, g.v))];
          st.c = log_add(st.c, s.a + E.at(s.p, j));
          ++l2r_tr[LT_S];
          if (grammar_simple) ++gs_phase;
        }
    }
  }

  // B seeds x'.j from S at j-1.
  void l2r_bseed(int j) {
    if (j < 1) return;
    lchk(LT_S, j - 1);
    const SMap& pv = LS[j - 1];
    for (size_t k = 0; k < pv.keys.size(); ++k) {
      const double f = pv.v[k].F;
      if (f < HALF_NEG) continue;
      const uint64_t key = pv.keys[k];
      const int xp = k3c(key);
      if (!E.pairs(xp, j) || !E.end_ok(xp, j)) continue;
      if (grammar_simple) {
        lacc(LB.Gs[j], LT_BGS, k3b(key), xp, k3a(key), f + E.w_stem + E.at(xp, j) + pair_bonus);
        gs_seedB.insert(key3(xp, k3a(key), k3b(key)));
        ++gs_phase;
      } else {
        lacc(LB.Gs[j], LT_BGS, k3b(key), xp, k3a(key), f + E.w_stem + E.at(xp, j) + pair_bonus);
      }
    }
  }

  void lc_cmp(double x, double y, long long& cnt, long long& mis, double& rel) {
    ++cnt;
    const bool fx = x > HALF_NEG, fy = y > HALF_NEG;
    if (fx != fy) {
      ++mis;
      rel = INFINITY;
      return;
    }
    if (!fx) return;
    const double r = std::expm1(std::fabs(x - y));
    rel = std::max(rel, r);
    if (r > 1e-10) ++mis;
  }

  void l2r_closeB(int j) {
    for (int src = 0; src < 2; ++src) {
      if (src == 1 && !single_stems) continue;
      for (const Grp& g : (src ? LB.Gs[j] : LB.G[j]).grps) {
        const int xs = g.v;
        if (!E.end_ok(xs, j)) continue;
        for (const Ent& s : g.ents) {
          if (s.p != xs) continue;
          if (grammar_simple && !single_stems && gs_seedB.count(key3(s.p, s.u, xs))) continue;
          const double val = s.a + E.at(xs, j);
          if (grammar_simple) ++gs_phase;
          if (l2r_mode == 1) set_pair(Qp, Qp_at, s.u, j, val);
          else LQp[ix(s.u, j)] = log_add(LQp[ix(s.u, j)], val);
        }
      }
    }
    if (l2r_mode == 2)
      for (int i = 0; i <= j; ++i) {
        const double x = Qp[ix(i, j)], y = LQp[ix(i, j)];
        if (x > HALF_NEG || y > HALF_NEG) lc_cmp(x, y, lc_Qp_n, lc_Qp_mis, lc_Qp_rel);
      }
  }

  // Check 1.1: Afull(a, xs, x', j) against sum_{a2, c1} Qg(a, a2, c1, j) w_stem at(a2, c1) F1 F2
  // from the engine's Qg[j] (same live-arm filter).
  void l2r_checkR(int j) {
    std::unordered_map<uint64_t, double> def, rec;
    for (const Grp& g : Qg[j].grps) {
      const int c1 = g.v;
      for (const Ent& s : g.ents) {
        const int a = s.p, a2 = s.u;
        if (!E.end_ok(a2, c1)) continue;
        const double base = s.a + E.w_stem + E.at(a2, c1);
        for (int xs = a2 + 1; xs < c1; ++xs) {
          if (lmaxR[xs] <= j) continue;
          const double f1 = qz_log(a2 + 1, xs - 1);
          if (f1 < HALF_NEG) continue;
          for (int xp = xs; xp < c1; ++xp) {
            if (lmaxR[xp] <= j) continue;
            const double f2 = qz_log(xp + 1, c1 - 1);
            if (f2 < HALF_NEG) continue;
            double& d = def.try_emplace(key3(a, xs, xp), NEG_INF).first->second;
            d = log_add(d, base + f1 + f2);
          }
        }
      }
    }
    for (int src = 0; src < 2; ++src) {
      if (src == 1 && !single_stems) continue;
      for (const Grp& g : (src ? LA.Gs[j] : LA.G[j]).grps)
        for (const Ent& s : g.ents) {
          double& d = rec.try_emplace(key3(s.p, s.u, g.v), NEG_INF).first->second;
          d = log_add(d, s.a);
        }
    }
    for (const auto& [k, x] : def) {
      auto it = rec.find(k);
      lc_cmp(x, it == rec.end() ? NEG_INF : it->second, lc_R_n, lc_R_mis, lc_R_rel);
    }
    for (const auto& [k, y] : rec)
      if (!def.count(k)) lc_cmp(NEG_INF, y, lc_R_n, lc_R_mis, lc_R_rel);
  }

  // S at the end of step j (after fill_Qz(j)): the Qz recursion of fill_Qz summed over q'.
  void l2r_sext(int j) {
    l2r_fext(LS, LT_S, j, [&](uint64_t key) { return l2r_arm(k3b(key), k3c(key), j); });
  }
  // One fill_Qz step of the filler table T at j (S, P1, P3): every state of T[j-1] / T[d-1]
  // with alive(key) at j takes one unpaired base / the item d.j; F = closed run + c (the
  // filler that starts after j, entered by the caller before this call).
  template <class Alive>
  void l2r_fext(std::vector<SMap>& T, int lt, int j, Alive alive) {
    const int fi = lt == LT_S ? 0 : lt == LT_P1 ? 1 : 2;
    std::unique_ptr<SecTimer> ptm(lpf ? new SecTimer(lpf_fx_sec[fi]) : nullptr);
    long long fpair = 0;
    if (lpf) ++lpf_fx_pos[fi];
    lchk(LT_QB, j);
    if (allow_pk) lchk(LT_QP, j);
    SMap& cur = T[j];
    auto branch = [&](int d, double bw) {
      if (d < 1) return;
      lchk(lt, d - 1);
      const SMap& pv = T[d - 1];
      if (lpf) {
        ++lpf_fx_items[fi];
        fpair += static_cast<long long>(pv.keys.size());
      }
      for (size_t k = 0; k < pv.keys.size(); ++k) {
        const double f = pv.v[k].F;
        const uint64_t key = pv.keys[k];
        if (f < HALF_NEG || !alive(key)) continue;
        SSt& s = cur.v[cur.at(key)];
        s.P = log_add(s.P, f + bw);
        ++l2r_tr[lt];
        if (grammar_simple) ++gs_fill;
      }
    };
    const OMap& atb = Qb_at[j];
    for (size_t r = 0; r < atb.slots(); ++r) {
      const int d = atb.keys[r];
      if (d < 0 || !E.end_ok(d, j)) continue;
      branch(d, atb.vals[r] + (E.w_nest + E.at(d, j)));
    }
    if (allow_pk) {
      const OMap& atp = Qp_at[j];
      for (size_t r = 0; r < atp.slots(); ++r) {
        const int d = atp.keys[r];
        if (d < 0) continue;
        branch(d, atp.vals[r] + E.w_beta1p + E.nest_pk_k * E.w_nest);
      }
    }
    if (lpf) {
      lpf_fx_pair[fi] += fpair;
      lpf_fx_pairmax[fi] = std::max(lpf_fx_pairmax[fi], fpair);
    }
    if (j >= 1) {
      lchk(lt, j - 1);
      const SMap& pv = T[j - 1];
      if (lpf) lpf_fx_unp[fi] += static_cast<long long>(pv.keys.size());
      for (size_t k = 0; k < pv.keys.size(); ++k) {
        const uint64_t key = pv.keys[k];
        if (!alive(key)) continue;
        const SSt& o = pv.v[k];
        if (dz) {
          const double u1 = log_add(o.P, o.c), u2 = log_add(o.U2, o.U1 + E.d5w[j - 1]);
          if (u1 < HALF_NEG && u2 < HALF_NEG) continue;
          SSt& s = cur.v[cur.at(key)];
          s.U1 = u1 > HALF_NEG ? u1 + E.wu : NEG_INF;
          s.U2 = u2 > HALF_NEG ? u2 + E.wu : NEG_INF;
        } else {
          if (o.F < HALF_NEG) continue;
          SSt& s = cur.v[cur.at(key)];
          s.P = log_add(s.P, o.F + E.wu);
        }
        ++l2r_tr[lt];
        if (grammar_simple) ++gs_fill;
      }
    }
    for (SSt& s : cur.v)
      s.F = log_add(dz ? run_close(s.P, s.U1, s.U2, E.s1w[j], E.d3w[j]) : s.P, s.c);
  }

  // ---- --l2r-spawn inc: P1 / P2 / P3 (forward).  alive: A (a2) and B (xs, x') can still
  // pair right of every c1 the state can reach (c1 >= j+1; the arm needs lmaxR > c1, as spawn).
  bool lpa1(int a2, int j) const { return lmaxR[a2] > j + 1; }
  bool lpa2(int a2, int xs, int j) const { return lmaxR[a2] > j && lmaxR[xs] > j + 1; }
  bool lpa3(int a2, int xs, int xp, int j) const {
    return lmaxR[a2] > j && lmaxR[xs] > j + 1 && lmaxR[xp] > j + 1;
  }
  bool l2r_pairable(int a2, int c1) const {
    const std::vector<int>& pl = partners_left[c1];
    return std::binary_search(pl.begin(), pl.end(), a2) && E.end_ok(a2, c1);
  }
  long long lp_seed_checks = 0, lp_seeds = 0;
  void l2r_p1_step(int j) {
    if (lpa1(j, j)) {
      SSt& s = LP1[j].v[LP1[j].at(key3(j, 0, 0))];
      s.c = 0.0;
      ++l2r_tr[LT_P1];
    }
    l2r_fext(LP1, LT_P1, j, [&](uint64_t key) { return lpa1(k3a(key), j); });
  }
  void l2r_p2_step(int j) {
    if (j < 1) return;
    SMap& cur = LP2[j];
    lchk(LT_P2, j - 1);
    const SMap& pv = LP2[j - 1];
    for (size_t k = 0; k < pv.keys.size(); ++k) {
      const uint64_t key = pv.keys[k];
      const double f = pv.v[k].F;
      if (f < HALF_NEG || !lpa2(k3a(key), k3b(key), j)) continue;
      cur.v[cur.at(key)].F = f;
      ++l2r_tr[LT_P2];
      if (grammar_simple) ++gs_fill;
    }
    lchk(LT_P1, j - 1);
    const SMap& p1 = LP1[j - 1];
    for (size_t k = 0; k < p1.keys.size(); ++k) {
      const int a2 = k3a(p1.keys[k]);
      const double f = p1.v[k].F;
      if (f < HALF_NEG || !lpa2(a2, j, j)) continue;
      cur.v[cur.at(key3(a2, j, 0))].F = f;
      ++l2r_tr[LT_P2];
      if (grammar_simple) ++gs_phase;
    }
  }
  void l2r_p3_step(int j) {
    SMap& cur = LP3[j];
    lchk(LT_P2, j);
    const SMap& p2 = LP2[j];
    for (size_t k = 0; k < p2.keys.size(); ++k) {
      const int a2 = k3a(p2.keys[k]), xs = k3b(p2.keys[k]);
      const double f = p2.v[k].F;
      if (f < HALF_NEG || !lpa3(a2, xs, j, j)) continue;
      SSt& s = cur.v[cur.at(key3(a2, xs, j))];
      s.c = log_add(s.c, f);
      ++l2r_tr[LT_P3];
      if (grammar_simple) ++gs_phase;
    }
    l2r_fext(LP3, LT_P3, j, [&](uint64_t key) { return lpa3(k3a(key), k3b(key), k3c(key), j); });
  }
  // A seeds at the inner pair a2.j from P3 at j-1 (spawn's arm rule: lmaxR[xs], lmaxR[x'] > j).
  void l2r_inc_seed(int j) {
    if (j < 1) return;
    lchk(LT_P3, j - 1);
    const SMap& pv = LP3[j - 1];
    for (size_t k = 0; k < pv.keys.size(); ++k) {
      const uint64_t key = pv.keys[k];
      const double f = pv.v[k].F;
      if (f < HALF_NEG) continue;
      const int a2 = k3a(key), xs = k3b(key), xp = k3c(key);
      if (lmaxR[xs] <= j || lmaxR[xp] <= j) continue;
      ++lp_seed_checks;
      if (!l2r_pairable(a2, j)) continue;
      ++lp_seeds;
      if (grammar_simple) {
        // Seeds are their own 4-index type (Rs), beamed at b4, not mixed into R.
        lacc(LA.Gs[j], LT_AGS, xp, a2, xs, (pair_bonus + E.w_stem + E.at(a2, j)) + f);
        gs_seedA.insert(key3(a2, xs, xp));
        ++gs_phase;
      } else {
        lacc(LA.Gs[j], LT_AGS, xp, a2, xs, (pair_bonus + E.w_stem + E.at(a2, j)) + f);
      }
    }
  }
  // --l2r-inc-check on (unbeamed): every A seed at j against spawn's k0 F1 F2 from the Qz table.
  bool l2r_inc_check = false;
  long long lc_I_n = 0, lc_I_mis = 0;
  double lc_I_rel = 0.0;
  void l2r_inc_check_step(int j) {
    std::unordered_map<uint64_t, double> def;
    for (int a2 : partners_left[j]) {
      if (!E.end_ok(a2, j)) continue;
      const double k0 = pair_bonus + E.w_stem + E.at(a2, j);
      for (int xs = a2 + 1; xs < j; ++xs) {
        if (lmaxR[xs] <= j) continue;
        const double f1 = qz_log(a2 + 1, xs - 1);
        if (f1 < HALF_NEG) continue;
        for (int xp = xs; xp < j; ++xp) {
          if (lmaxR[xp] <= j) continue;
          const double f2 = qz_log(xp + 1, j - 1);
          if (f2 < HALF_NEG) continue;
          def.emplace(key3(a2, xs, xp), k0 + f1 + f2);
        }
      }
    }
    std::unordered_map<uint64_t, double> rec;
    for (const Grp& g : LA.Gs[j].grps)
      for (const Ent& s : g.ents) rec.emplace(key3(s.p, s.u, g.v), s.a);
    for (const auto& [k, x] : def) {
      auto it = rec.find(k);
      lc_cmp(x, it == rec.end() ? NEG_INF : it->second, lc_I_n, lc_I_mis, lc_I_rel);
    }
    for (const auto& [k, y] : rec)
      if (!def.count(k)) lc_cmp(NEG_INF, y, lc_I_n, lc_I_mis, lc_I_rel);
  }

  // ---- beam (--pk-l2r on with a finite --l2r-beam; needs the principled pre-pass, --grow-rec)
  // Scores (log): R-type (R, Rs, Rls, Rlss, Rrs): alpha + log sum_q' g_L[q'] K_q'(arm), g_L the
  // (truncated) left growth vector of rec_root / rec_key, which already holds C[p'-1];
  // K_q'(arm) = sum_{y' > q'} Qz(q'+1, y'-1) E_in V (the engine's arm_phi, separable with
  // --pr-armq-sep).  U-type: alpha + g_R[xs] (V of the outer pair, with e^beta1 Sbar) +
  // C[a'-1].  S: F + K_j(arm) + C[a'-1].  Spawn: arms ranked by F1 + s_j(arm) + F2,
  // s_j = log sum_{y' > j} E_in V; top l2r_spawn_k (0: every arm, unranked).
  enum { PT_R, PT_RS, PT_RLS, PT_RLSS, PT_RRS, PT_U, PT_US, PT_ULS, PT_ULSS, PT_URS, PT_S,
         PT_P1, PT_P2, PT_P3, PT_N };
  int l2r_beam = 0, l2r_spawn_k = 0;
  // -1: l2r_beam; 0: none
  int l2r_tb[PT_N] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
  int ltb(int pt) const { return l2r_tb[pt] < 0 ? l2r_beam : l2r_tb[pt]; }
  std::unordered_map<uint64_t, int> lArm;
  std::vector<double> lYwS;
  long long l2r_pre[PT_N] = {}, l2r_kept[PT_N] = {}, l2r_kmax[PT_N] = {};
  long long l2r_pos[PT_N] = {}, l2r_bind[PT_N] = {};
  // scored states (pruning steps only), scores of -inf, R / S scores with no V arm
  long long l2r_nsc[PT_N] = {}, l2r_ninf[PT_N] = {}, l2r_noarm[PT_N] = {};
  // --l2r-ninf-mass on: per table, the states a beam cut scores -inf: log sum of their inside
  // values, kept / dropped counts, and the kept ones' keys per step (their posterior after the
  // outside pass).  --l2r-keep-ninf on: they are all kept, beyond the beam.
  bool l2r_ninf_mass = false, l2r_keep_ninf = false;
  std::vector<double> ninf_in = std::vector<double>(PT_N, NEG_INF);
  long long ninf_kept[PT_N] = {}, ninf_drop[PT_N] = {};
  std::vector<std::unordered_set<uint64_t>> ninf_keys[PT_N];
  static uint64_t ninf_key(int a, int b, int c) {
    return (static_cast<uint64_t>(a + 1) << 42) | (static_cast<uint64_t>(b + 1) << 21) |
           static_cast<uint64_t>(c + 1);
  }
  void ninf_note(int pt, int j, uint64_t key, double in, char& keep) {
    if (l2r_keep_ninf) keep = 1;
    if (!l2r_ninf_mass) return;
    if (in > HALF_NEG) ninf_in[pt] = log_add(ninf_in[pt], in);
    if (!keep) {
      ++ninf_drop[pt];
      return;
    }
    ++ninf_kept[pt];
    if (ninf_keys[pt].empty()) ninf_keys[pt].resize(n);
    ninf_keys[pt][j].insert(key);
  }
  // After the outside pass: per table "name ninf kept dropped log10(sum inside / Z) posterior"
  // of the -inf-scored states (posterior of the kept ones), then the posterior total.
  std::string l2r_ninf_report() {
    static const char* pn[PT_N] = {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us", "Uls", "Ulss",
                                   "Urs", "S", "P1", "P2", "P3"};
    outside();
    char buf[200];
    std::snprintf(buf, sizeof buf, " logZ %.9f | ninf kept dropped log10_in/Z post", logZ);
    std::string out = buf;
    double all = 0.0;
    for (int pt = 0; pt < PT_N; ++pt) {
      double post = 0.0;
      if (!ninf_keys[pt].empty())
        for (int j = 0; j < n; ++j) {
          const std::unordered_set<uint64_t>& ks = ninf_keys[pt][j];
          if (ks.empty()) continue;
          if (pt >= PT_S) {
            const SMap& m = l2r_smap(pt, j);
            for (size_t k = 0; k < m.keys.size(); ++k)
              if (ks.count(m.keys[k]) && m.v[k].F > HALF_NEG && m.v[k].bF > HALF_NEG)
                post += std::exp(m.v[k].F + m.v[k].bF - logZ);
          } else {
            for (const Grp& g : l2r_tab(pt, j).grps)
              for (const Ent& e : g.ents)
                if (ks.count(ninf_key(g.v, e.u, e.p)) && e.a > HALF_NEG && e.b > HALF_NEG)
                  post += std::exp(e.a + e.b - logZ);
          }
        }
      all += post;
      std::snprintf(buf, sizeof buf, " %s %lld %lld %lld %.4f %.6g", pn[pt], l2r_ninf[pt], ninf_kept[pt],
                    ninf_drop[pt], ninf_in[pt] > HALF_NEG ? (ninf_in[pt] - logZ) / std::log(10.0) : -999.0,
                    post);
      out += buf;
    }
    std::snprintf(buf, sizeof buf, " | post_total %.6g", all);
    return out + buf;
  }
  double l2r_bsec[PT_N] = {}, l2r_psec[PT_N] = {}, l2r_csec = 0.0, l2r_isec = 0.0;
  long long sp_events = 0, sp_cand = 0, sp_kept = 0;

