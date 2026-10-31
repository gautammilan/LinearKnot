// Pass 3, the pseudoknot part of step j (Lead, simple grammar M1-M8, l2r_step_pk, step end)
// and the scan outside (--pk-l2r on): each step mirrored in reverse (l2r_*_back).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.


  // M1 and M4 from the stems kept at j. M3 records each M4 state as an M as well.
  // M4's group key is the carried v. Unpaired bases inside Qhat^M are already in that table.
  void l2r_simple_m4(LFam& F, int j) {
    auto from = [&](const Tab4& tab) {
      for (const Grp& g : tab.grps) {
        const int lo = F.isU ? g.v + 1 : 1;
        for (const Ent& s : g.ents) {
          if (E.end_ok(s.p, j)) {
            F.SM[j].acc(g.v, s.p, s.u, s.a + E.w_bp + E.at(s.p, j));
            ++gs_m[0];
          }
          const int c = s.p;
          if (c - 1 < 0 || !E.end_ok(c, j)) continue;
          lchk(LT_QM, c - 1);
          const double w = E.w_bp + E.at(c, j);
          const OMap& at = qmp_at(c - 1);
          for (size_t r = 0; r < at.slots(); ++r) {
            const int i = at.keys[r];
            if (i < lo) continue;
            F.SM2[j].acc(g.v, i, s.u, s.a + w + at.vals[r]);
            ++gs_m[3];
          }
        }
      }
    };
    from(F.G[j]);
    from(F.Gs[j]);
    for (const Grp& g : F.SM2[j].grps) {
      if (g.v >= n) continue;
      for (const Ent& s : g.ents) {
        F.SM[j].acc(g.v, s.p, s.u, s.a);
        ++gs_m[2];
      }
    }
  }

  // M5: stem ended at f = p-1, nested segment p..j on its right. Group key is v + n.
  void l2r_simple_m5(LFam& F, int j) {
    lchk(LT_QM, j);
    const OMap& at = qmp_at(j);
    const int stem = F.isU ? LT_BG : LT_AG;
    for (size_t r = 0; r < at.slots(); ++r) {
      const int p = at.keys[r];
      if (p < 1) continue;
      const int f = p - 1;
      const double qm = at.vals[r];
      lchk(stem, f);
      auto take = [&](const Tab4& tab) {
        for (const Grp& g : tab.grps)
          for (const Ent& s : g.ents) {
            if (!E.end_ok(s.p, f)) continue;
            if (!F.isU && !l2r_arm(s.u, g.v, j)) continue;
            F.SM2[j].acc(g.v + n, s.p, s.u, s.a + E.w_bp + E.at(s.p, f) + qm);
            ++gs_m[4];
          }
      };
      take(F.G[f]);
      take(F.Gs[f]);
    }
  }

  // Closes into X at j. (f) is M6+M7+M8 from M4. (g) is M6+M8 from M5 through Lead.
  // (h) is M5 plus a nested segment on the left. Survivors of the Multi beam are written.
  void l2r_simple_close(LFam& F, int j) {
    Tab4& Mu = F.SMulti[j];
    const int tg = F.isU ? LT_BG : LT_AG;
    for (int f = 0; f < j; ++f) {
      const int unpaired = j - f - 1;
      double w_ext = E.w_ap + E.w_bp + E.w_cp * unpaired + pair_bonus;
      if (dz) w_ext += E.run(f + 1, j - 1);
      for (const Grp& g : F.SM2[f].grps) {
        if (g.v >= n) continue;
        const int lo = F.isU ? g.v : 0;
        for (const Ent& s : g.ents) {
          const int i = s.p - 1;
          if (i < lo || !E.pairs(i, j) || !E.end_ok(i, j)) continue;
          if (!F.isU && !l2r_arm(s.u, g.v, j)) continue;
          Mu.acc(g.v, i, s.u, s.a + w_ext + E.at(i, j));
          ++gs_m[5];
          gs_m[6] += unpaired;
          ++gs_m[7];
          if (unpaired > cap) ++gs_over;
        }
      }
    }
    if (j >= 1) {
      for (const Grp& g : F.SM2[j - 1].grps) {
        if (g.v < n) continue;
        const int real_v = g.v - n;
        const int lo = F.isU ? real_v : 0;
        for (const Ent& s : g.ents) {
          const int c = s.p;
          if (!F.isU && !l2r_arm(s.u, real_v, j)) continue;
          // Leading run is the same window as pass 1: p in [c - ml_lead_cap, c).
          // ml_lead_cap 0 keeps every p. No Qhat^Lead table and no beam on the run.
          {
            const int i_lo = std::max(lo, lead_i_lo(c));
            for (int i = i_lo; i < c; ++i) {
              if (!E.pairs(i, j) || !E.end_ok(i, j)) continue;
              const int gap = c - i - 1;
              double lw = 0.0;
              if (gap > 0) {
                lw = E.w_cp * gap;
                if (dz) lw += E.run(i + 1, c - 1);
              }
              Mu.acc(real_v, i, s.u, s.a + lw + E.w_ap + E.w_bp + E.at(i, j) + pair_bonus);
              ++gs_m[5];
              ++gs_m[7];
              ++gs_lead;
              if (gap > cap) ++gs_over;
            }
          }
          if (c - 1 < 0) continue;
          lchk(LT_QM, c - 1);
          const OMap& at = qmp_at(c - 1);
          qm_pair_slots(c - 1, j, [&](int r) {
            const int ii = at.keys[r];
            if (ii < 0) return;
            const int outer = ii - 1;
            if (outer < lo || outer >= c || !E.pairs(outer, j) || !E.end_ok(outer, j)) return;
            Mu.acc(real_v, outer, s.u,
                   at.vals[r] + s.a + E.w_ap + E.w_bp + E.at(outer, j) + pair_bonus);
            ++gs_m[7];
          });
        }
      }
    }
    l2r_prune_prefix(Mu, l2r_beam, F.isU ? 5 : 2);
    for (const Grp& g : Mu.grps)
      for (const Ent& s : g.ents) lacc(F.G[j], tg, g.v, s.p, s.u, s.a);
  }

  void l2r_simple_pk(int j) {
    if (!l2r_inc) {
      std::fprintf(stderr, "[grammar] simple needs --l2r-spawn inc\n");
      std::abort();
    }
    if (l2r_beam > 0) stem_cache_step(j);
    gs_seedA.clear();
    gs_seedB.clear();
    {
      SecTimer t(l2r_bsec[PT_R]);
      l2r_grow(LA, j);
    }
    {
      SecTimer t(l2r_bsec[PT_RS]);
      l2r_inc_seed(j);
    }
    l2r_pruneA(LA.Gs[j], PT_RS, 0, true, j);
    l2r_simple_close(LA, j);
    l2r_pruneA(LA.G[j], PT_R, 0, false, j);
    l2r_done[LT_AG] = l2r_done[LT_AGS] = j;
    l2r_simple_m4(LA, j);
    {
      SecTimer t(l2r_csec);
      l2r_closeA(j);
    }
    l2r_done[LT_ALS] = j;
    {
      SecTimer t(l2r_bsec[PT_US]);
      l2r_bseed(j);
    }
    l2r_pruneB(LB.Gs[j], PT_US, 0, true, j);
    {
      SecTimer t(l2r_bsec[PT_U]);
      l2r_grow(LB, j);
    }
    l2r_simple_close(LB, j);
    l2r_pruneB(LB.G[j], PT_U, 0, false, j);
    l2r_done[LT_BG] = l2r_done[LT_BGS] = j;
    l2r_simple_m4(LB, j);
    l2r_done[LT_BLS] = j;
    {
      SecTimer t(l2r_csec);
      l2r_closeB(j);
    }
  }

  void l2r_simple_end(int j) {
    {
      SecTimer t(l2r_bsec[PT_S]);
      l2r_sext(j);
    }
    l2r_pruneS(j);
    l2r_done[LT_S] = j;
    {
      SecTimer t(l2r_bsec[PT_P1]);
      l2r_p1_step(j);
    }
    l2r_pruneP(PT_P1, j);
    l2r_done[LT_P1] = j;
    {
      SecTimer t(l2r_bsec[PT_P2]);
      l2r_p2_step(j);
    }
    l2r_pruneP(PT_P2, j);
    l2r_done[LT_P2] = j;
    {
      SecTimer t(l2r_bsec[PT_P3]);
      l2r_p3_step(j);
    }
    l2r_pruneP(PT_P3, j);
    l2r_done[LT_P3] = j;
    l2r_simple_m5(LA, j);
    l2r_simple_m5(LB, j);
    l2r_prune_prefix(LA.SM[j], l2r_beam, 0);
    l2r_prune_prefix(LB.SM[j], l2r_beam, 3);
    l2r_prune_prefix(LA.SM2[j], l2r_beam, 1);
    l2r_prune_prefix(LB.SM2[j], l2r_beam, 4);
    l2r_done[LT_ARS] = l2r_done[LT_BRS] = j;
  }

  void l2r_step_pk(int j) {
    if (grammar_simple) {
      l2r_simple_pk(j);
      return;
    }
    if (l2r_beam > 0) stem_cache_step(j);
    {
      SecTimer t(l2r_bsec[PT_R]);
      l2r_grow(LA, j);
    }
    {
      SecTimer t(l2r_bsec[PT_RS]);
      if (l2r_inc) l2r_inc_seed(j);
      else l2r_spawn(j);
    }
    l2r_pruneA(LA.G[j], PT_R, 0, false, j);
    l2r_pruneA(LA.Gs[j], PT_RS, 0, true, j);
    l2r_done[LT_AG] = l2r_done[LT_AGS] = j;
    {
      SecTimer t(l2r_bsec[PT_RLS]);
      l2r_ls(LA, j);
    }
    l2r_pruneA(LA.Gls[j], PT_RLS, 1, false, j);
    l2r_pruneA(LA.Glss[j], PT_RLSS, 1, true, j);
    l2r_done[LT_ALS] = j;
    {
      SecTimer t(l2r_csec);
      l2r_closeA(j);
    }
    {
      SecTimer t(l2r_bsec[PT_US]);
      l2r_bseed(j);
    }
    {
      SecTimer t(l2r_bsec[PT_U]);
      l2r_grow(LB, j);
    }
    l2r_pruneB(LB.G[j], PT_U, 0, false, j);
    l2r_pruneB(LB.Gs[j], PT_US, 0, true, j);
    l2r_done[LT_BG] = l2r_done[LT_BGS] = j;
    {
      SecTimer t(l2r_bsec[PT_ULS]);
      l2r_ls(LB, j);
    }
    l2r_pruneB(LB.Gls[j], PT_ULS, 1, false, j);
    l2r_pruneB(LB.Glss[j], PT_ULSS, 1, true, j);
    l2r_done[LT_BLS] = j;
    {
      SecTimer t(l2r_csec);
      l2r_closeB(j);
    }
  }

  // End of step j: S (needs Qb, Qp, Qz at j), then Grs (needs Qm(., j)).
  void l2r_step_end(int j) {
    if (grammar_simple) {
      l2r_simple_end(j);
      return;
    }
    {
      SecTimer t(l2r_bsec[PT_S]);
      l2r_sext(j);
    }
    l2r_pruneS(j);
    l2r_done[LT_S] = j;
    if (l2r_inc) {
      {
        SecTimer t(l2r_bsec[PT_P1]);
        l2r_p1_step(j);
      }
      l2r_pruneP(PT_P1, j);
      l2r_done[LT_P1] = j;
      {
        SecTimer t(l2r_bsec[PT_P2]);
        l2r_p2_step(j);
      }
      l2r_pruneP(PT_P2, j);
      l2r_done[LT_P2] = j;
      {
        SecTimer t(l2r_bsec[PT_P3]);
        l2r_p3_step(j);
      }
      l2r_pruneP(PT_P3, j);
      l2r_done[LT_P3] = j;
    }
    {
      SecTimer t(l2r_bsec[PT_RRS]);
      l2r_rs(LA, j);
    }
    l2r_pruneA(LA.Grs[j], PT_RRS, 2, false, j);
    {
      SecTimer t(l2r_bsec[PT_URS]);
      l2r_rs(LB, j);
    }
    l2r_pruneB(LB.Grs[j], PT_URS, 2, false, j);
    l2r_done[LT_ARS] = l2r_done[LT_BRS] = j;
  }

  // ---- outside (--pk-l2r on): each step mirrors the forward one in reverse order.
  void l2r_sext_back(int j) {
    l2r_fext_back(LS, j, [&](uint64_t key) { return l2r_arm(k3b(key), k3c(key), j); });
  }
  // Adjoint of l2r_fext at j: bF of T[j] (complete) into bP / bU1 / bU2 / bc, then into the
  // items d.j and T[j-1], T[d-1].
  template <class Alive>
  void l2r_fext_back(std::vector<SMap>& T, int j, Alive alive) {
    SMap& cur = T[j];
    for (SSt& s : cur.v) {
      if (s.bF < HALF_NEG) continue;
      s.bP = log_add(s.bP, s.bF);
      s.bc = log_add(s.bc, s.bF);
      if (dz) {
        s.bU1 = log_add(s.bU1, s.bF + E.s1w[j]);
        s.bU2 = log_add(s.bU2, s.bF + E.d3w[j]);
      }
    }
    auto branch = [&](int d, double qv, double bw, SpMat& qb) {
      if (d < 1) return;
      SMap& pv = T[d - 1];
      for (size_t k = 0; k < pv.keys.size(); ++k) {
        SSt& o = pv.v[k];
        const uint64_t key = pv.keys[k];
        if (o.F < HALF_NEG || !alive(key)) continue;
        const int c = cur.find(key);
        if (c < 0 || cur.v[c].bP < HALF_NEG) continue;
        const double bp = cur.v[c].bP;
        o.bF = log_add(o.bF, bp + qv + bw);
        add2b(qb, d, j, bp + o.F + bw);
      }
    };
    const OMap& atb = Qb_at[j];
    for (size_t r = 0; r < atb.slots(); ++r) {
      const int d = atb.keys[r];
      if (d < 0 || !E.end_ok(d, j)) continue;
      branch(d, atb.vals[r], E.w_nest + E.at(d, j), Qb_b);
    }
    if (allow_pk) {
      const OMap& atp = Qp_at[j];
      for (size_t r = 0; r < atp.slots(); ++r) {
        const int d = atp.keys[r];
        if (d < 0) continue;
        branch(d, atp.vals[r], E.w_beta1p + E.nest_pk_k * E.w_nest, Qp_b);
      }
    }
    if (j < 1) return;
    SMap& pv = T[j - 1];
    for (size_t k = 0; k < pv.keys.size(); ++k) {
      const uint64_t key = pv.keys[k];
      if (!alive(key)) continue;
      const int c = cur.find(key);
      if (c < 0) continue;
      const SSt& s = cur.v[c];
      SSt& o = pv.v[k];
      if (dz) {
        if (s.bU1 > HALF_NEG) {
          o.bP = log_add(o.bP, s.bU1 + E.wu);
          o.bc = log_add(o.bc, s.bU1 + E.wu);
        }
        if (s.bU2 > HALF_NEG) {
          o.bU2 = log_add(o.bU2, s.bU2 + E.wu);
          o.bU1 = log_add(o.bU1, s.bU2 + E.wu + E.d5w[j - 1]);
        }
      } else if (s.bP > HALF_NEG && o.F > HALF_NEG) {
        o.bF = log_add(o.bF, s.bP + E.wu);
      }
    }
  }

  // --l2r-spawn inc, outside at the start of step j: P3 (filler, then the hole closed at
  // x' = j into P2), P2 (carried into P2[j-1], or opened at xs = j from P1[j-1]), P1.
  void l2r_inc_back(int j) {
    l2r_fext_back(LP3, j, [&](uint64_t key) { return lpa3(k3a(key), k3b(key), k3c(key), j); });
    SMap& p2 = LP2[j];
    for (size_t k = 0; k < LP3[j].keys.size(); ++k) {
      const uint64_t key = LP3[j].keys[k];
      const SSt& s = LP3[j].v[k];
      if (k3c(key) != j || s.c < HALF_NEG || s.bc < HALF_NEG) continue;
      const int t = p2.find(key3(k3a(key), k3b(key), 0));
      if (t >= 0) p2.v[t].bF = log_add(p2.v[t].bF, s.bc);
    }
    if (j >= 1) {
      SMap& pv = LP2[j - 1];
      SMap& p1 = LP1[j - 1];
      for (size_t k = 0; k < p2.keys.size(); ++k) {
        const double b = p2.v[k].bF;
        if (b < HALF_NEG) continue;
        const uint64_t key = p2.keys[k];
        SMap& src = k3b(key) < j ? pv : p1;
        const int t = src.find(k3b(key) < j ? key : key3(k3a(key), 0, 0));
        if (t >= 0) src.v[t].bF = log_add(src.v[t].bF, b);
      }
    }
    l2r_fext_back(LP1, j, [&](uint64_t key) { return lpa1(k3a(key), j); });
  }

  // Outside of the A seeds from P3: Rs(a2, xs, x', j) -> P3(a2, xs, x', j-1).
  void l2r_inc_seed_back(int j) {
    if (j < 1) return;
    SMap& pv = LP3[j - 1];
    for (Grp& g : LA.Gs[j].grps)
      for (Ent& s : g.ents) {
        if (s.b < HALF_NEG) continue;
        const int k = pv.find(key3(s.p, s.u, g.v));
        if (k < 0) continue;
        SSt& o = pv.v[k];
        o.bF = log_add(o.bF, s.b + (pair_bonus + E.w_stem + E.at(s.p, j)));
      }
  }

  void l2r_closeA_back(int j) {
    const SMap& S = LS[j];
    for (int src = 0; src < 2; ++src) {
      if (src == 1 && !single_stems) continue;
      for (Grp& g : (src ? LA.Gs[j] : LA.G[j]).grps)
        for (Ent& s : g.ents) {
          if (!E.end_ok(s.p, j)) continue;
          const int k = S.find(key3(s.p, s.u, g.v));
          if (k >= 0) add4(s, S.v[k].bc + E.at(s.p, j));
        }
    }
  }

  void l2r_bseed_back(int j) {
    if (j < 1) return;
    SMap& pv = LS[j - 1];
    for (Grp& g : LB.Gs[j].grps)
      for (Ent& s : g.ents) {
        if (s.b < HALF_NEG) continue;
        const int k = pv.find(key3(s.u, g.v, s.p));
        if (k < 0) continue;
        SSt& o = pv.v[k];
        o.bF = log_add(o.bF, s.b + E.w_stem + E.at(s.p, j) + pair_bonus);
      }
  }

  void l2r_closeB_back(int j) {
    for (int src = 0; src < 2; ++src) {
      if (src == 1 && !single_stems) continue;
      for (Grp& g : (src ? LB.Gs[j] : LB.G[j]).grps) {
        const int xs = g.v;
        if (!E.end_ok(xs, j)) continue;
        for (Ent& s : g.ents)
          if (s.p == xs) {
            add4(s, static_cast<double>(Qp_b[ix(s.u, j)]) + E.at(xs, j));
          }
      }
    }
  }

  // Start of the outside step j (before the engine's): Urs, Rrs, S.
  void l2r_back_end(int j) {
    l2r_adj = true;
    l2r_rs(LB, j);
    l2r_rs(LA, j);
    l2r_sext_back(j);
    if (l2r_inc) l2r_inc_back(j);
    l2r_adj = false;
  }

  // At the engine's Qp adjoint (Qp outside at j final): close B, U side, B seeds, close A, R side.
  void l2r_back_pk(int j) {
    l2r_adj = true;
    l2r_closeB_back(j);
    l2r_ls(LB, j);
    l2r_grow(LB, j);
    l2r_bseed_back(j);
    l2r_closeA_back(j);
    l2r_ls(LA, j);
    if (l2r_inc) l2r_inc_seed_back(j);
    else l2r_spawn_back(j);
    l2r_grow(LA, j);
    l2r_adj = false;
  }

  // Stem pairs: the outer pair of every R / U state and the inner pair of every seed.
  template <class AddPair>
  void l2r_pairs(AddPair& add_pair) {
    for (LFam* F : {&LA, &LB})
      for (int j = 0; j < n; ++j)
        for (const std::vector<Tab4>* T : {&F->G, &F->Gs})
          for (const Grp& g : (*T)[j].grps)
            for (const Ent& s : g.ents) add_pair(s.p, j, s.a + s.b - logZ);
  }

