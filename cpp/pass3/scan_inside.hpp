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
      prune1(at);
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
    prune1(at);
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
  struct LFam {
    bool isU = false;
    std::vector<Tab4> G, Gs, Gls, Glss, Grs;
    // --grammar simple: M, M2, Multi. M2 groups with v >= n are M5 (stem then nested);
    // v < n are M4 (nested then stem). The logical group key is v % n.
    std::vector<Tab4> SM, SM2, SMulti;
  };
  LFam LA, LB;
  bool grammar_simple = false;
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
    FastMap<uint64_t, int> idx;
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
  // Right-flank tail extensions, recorded in the scan stats (always 0).
  long long rlt_ext = 0;
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
            for (int i = i0; i < c; ++i) {
              if (!l2r_adj) ++it_int[iu];
              if (!E.pairs(i, j)) continue;
              take(i);
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

  // A seeds at the inner pair a2.j: every arm a2 < xs <= x' < j whose B can close right of j.
  void l2r_spawn(int j) {
    std::vector<std::pair<double, int>> cand;
    for (int a2 : partners_left[j]) {
      if (!E.end_ok(a2, j)) continue;
      const double k0 = pair_bonus + E.w_stem + E.at(a2, j);
      ++sp_events;
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

  // S at the end of step j (after fill_Qz(j)): the Qz recursion of fill_Qz summed over q'.
  void l2r_sext(int j) {
    l2r_fext(LS, LT_S, j, [&](uint64_t key) { return l2r_arm(k3b(key), k3c(key), j); });
  }
  // One fill_Qz step of the filler table T at j (S, P1, P3): every state of T[j-1] / T[d-1]
  // with alive(key) at j takes one unpaired base / the item d.j; F = closed run + c (the
  // filler that starts after j, entered by the caller before this call).
  template <class Alive>
  void l2r_fext(std::vector<SMap>& T, int lt, int j, Alive alive) {
    lchk(LT_QB, j);
    if (allow_pk) lchk(LT_QP, j);
    SMap& cur = T[j];
    auto branch = [&](int d, double bw) {
      if (d < 1) return;
      lchk(lt, d - 1);
      const SMap& pv = T[d - 1];
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
    if (j >= 1) {
      lchk(lt, j - 1);
      const SMap& pv = T[j - 1];
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

  // ---- beam (--pk-l2r on with a finite --b4; needs the principled pre-pass)
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
  FastMap<uint64_t, int> lArm;
  std::vector<double> lYwS;
  long long l2r_pre[PT_N] = {}, l2r_kept[PT_N] = {}, l2r_kmax[PT_N] = {};
  long long l2r_pos[PT_N] = {}, l2r_bind[PT_N] = {};
  // scored states (pruning steps only), scores of -inf, R / S scores with no V arm
  long long l2r_nsc[PT_N] = {}, l2r_ninf[PT_N] = {}, l2r_noarm[PT_N] = {};
  double l2r_bsec[PT_N] = {}, l2r_psec[PT_N] = {}, l2r_csec = 0.0, l2r_isec = 0.0;
  long long sp_events = 0, sp_cand = 0, sp_kept = 0;

