// Outside and pair probabilities: outside (nested and DP03), outside_from_future,
// outside_same_j, collect_bpp (P(i, j) = Qb(i, j) Qb_b(i, j) / Z).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  void outside() {
    if (!inside_done) pfunc();
    if (n == 0 || outside_done) return;
    outside_done = true;
    alloc_beta();
    if (lsp_rec)
      for (auto& v : QpS) v.assign(static_cast<size_t>(n) * n, NEG_INF);
    if (keep_closed_b && dz) {
      QmC_b.init(n);
      Qm2C_b.init(n);
    }
    for (int j = n - 1; j >= 0; --j) {
      if (allow_pk && l2r_mode == 1) l2r_back_end(j);
      outside_from_future(j);
      outside_same_j(j);
    }
  }

  void outside_from_future(int j) {
    const double w3 = E.w_alpha3;
    const double wu = E.wu;

    // PKM[j] = PKM[j-1] + α3 and + Qp[d][j]; its β is complete (Qb at j+1, PKM[j+1])
    if (allow_pk && dz) {
      const RunSt& nx = rbPKM[(j + 1) & 1];
      RunSt& cu = rbPKM[j & 1];
      cu.clear();
      const OMap& pkm = PKM_at[j];
      for (size_t r = 0; r < pkm.slots(); ++r) {
        int d = pkm.keys[r];
        if (d < 0) continue;
        double bP, b1, b2;
        run_back(PKM_b[ix(j, d)], nx.U1.get(d), nx.U2.get(d), j, w3, bP, b1, b2);
        cu.U1.put(d, b1);
        cu.U2.put(d, b2);
        if (keep_runs) run_weight(RS_PKM, d, j, bP, b1, b2);
        if (Qp_at[j].has(d)) {
          add2b(Qp_b, d, j, bP);
          qps(1, d, j, bP);
        }
      }
    } else if (allow_pk) {
      const OMap& pkm = PKM_at[j];
      for (size_t r = 0; r < pkm.slots(); ++r) {
        int d = pkm.keys[r];
        if (d < 0) continue;
        double pb = PKM_b[ix(j, d)];
        if (pb < HALF_NEG) continue;
        if (Qp_at[j].has(d)) {
          add2b(Qp_b, d, j, pb);
          qps(1, d, j, pb);
        }
        if (j > 0 && PKM_at[j - 1].has(d))
          PKM_b[ix(j - 1, d)] = log_add(static_cast<double>(PKM_b[ix(j - 1, d)]), pb + w3);
      }
    }

    if (dz) {
      double b1, b2;
      run_back(C_b[j], bc1, bc2, j, 0.0, C_b[j], b1, b2);
      bc1 = b1;
      bc2 = b2;
      const int nx = (j + 1) & 1, cu = j & 1;
      for (int i = 0; i <= j; ++i) {
        bz1[cu][i] = bz2[cu][i] = NEG_INF;
        if (Qz[ix(i, j)] > HALF_NEG) {
          double bP;
          run_back(Qz_b[ix(i, j)], bz1[nx][i], bz2[nx][i], j, wu, bP, bz1[cu][i], bz2[cu][i]);
          Qz_b[ix(i, j)] = bP;
        }
      }
    } else if (j + 1 < n) {
      C_b[j] = log_add(C_b[j], C_b[j + 1]);
      const OMap& am = Qm_at[j];
      for (size_t r = 0; r < am.slots(); ++r) {
        int i = am.keys[r];
        if (i < 0) continue;
        add2b(Qm_b, i, j, Qm_b[ix(i, j + 1)] + w3);
      }
      const OMap& am2 = Qm2_at[j];
      for (size_t r = 0; r < am2.slots(); ++r) {
        int i = am2.keys[r];
        if (i < 0) continue;
        add2b(Qm2_b, i, j, Qm2_b[ix(i, j + 1)] + w3);
      }
      if (allow_pk && E.stem_ml) {
        const OMap& amp = Qmp_at[j];
        for (size_t r = 0; r < amp.slots(); ++r) {
          int i = amp.keys[r];
          if (i < 0) continue;
          add2b(Qmp_b, i, j, Qmp_b[ix(i, j + 1)] + E.w_cp);
        }
      }
      for (int i = 0; i <= j; ++i)
        if (Qz[ix(i, j)] > HALF_NEG) add2b(Qz_b, i, j, Qz_b[ix(i, j + 1)] + wu);
    }

    // Qb interior wrap
    int jmax = std::min(n, j + 2 + cap);  // outer q with q - j - 1 <= cap
    {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        int d = at.keys[r];
        if (d < 0) continue;
        if (l2r_mode == 0) {
          const int p_lo = std::max(0, d - 1 - cap);
          for (int p = d - 1; p >= p_lo; --p) {
            const int budget = cap - (d - p - 1);
            if (budget < 0) break;
            for (int q = next_mate(p, j); q != -1 && q - j - 1 <= budget; q = next_mate(p, q)) {
              if (!Qb_at[q].has(p)) continue;
              const double w = E.interior(p, d, j, q);
              if (w < HALF_NEG) continue;
              add2b(Qb_b, d, j, Qb_b[ix(p, q)] + w);
            }
          }
          continue;
        }
        for (int q = j + 1; q < jmax; ++q) {
          const int budget = cap - (q - j - 1);
          const int i_hi = d - 1;
          const int i_lo = std::max(0, d - 1 - budget);
          if (i_lo > i_hi) continue;
          const std::vector<int>& pl = partners_left[q];
          auto it = std::lower_bound(pl.begin(), pl.end(), i_lo);
          const auto iend = std::upper_bound(it, pl.end(), i_hi);
          for (; it != iend; ++it) {
            const int i = *it;
            double w = E.interior(i, d, j, q);
            if (w < HALF_NEG) continue;
            add2b(Qb_b, d, j, Qb_b[ix(i, q)] + w);
          }
        }
      }
    }

    if (!allow_pk) return;

    // Qg (b)
    for (Grp& g : Qg[j].grps) {
      int e = g.v;
      for (Ent& s : g.ents) {
        int c = s.p, d = s.u;
        if (c == d) continue;
        for (int q = j + 1; q < jmax; ++q) {
          int i_lo = std::max(0, c - 1 - cap);
          for (int i = i_lo; i < c; ++i) {
            double w = E.stem_interior(i, c, j, q);
            if (w < HALF_NEG) continue;
            double pb = g4(Qg[q], e, i, d);
            add4(s, pb + w);
          }
        }
      }
    }

    // Qg (f)
    for (int q = j + 1; q < n; ++q) {
      int unpaired = q - j - 1;
      double w_ext = E.w_ap + E.w_bp + E.w_cp * unpaired;
      if (dz) w_ext += E.run(j + 1, q - 1);
      for (Grp& g : Qgls[j].grps) {
        int e = g.v;
        if (e == j && unpaired > cap) continue;
        for (Ent& s : g.ents) {
          int i = s.p - 1;
          if (!E.pairs(i, q) || !E.end_ok(i, q)) continue;
          double pb = g4(Qg[q], e, i, s.u);
          add4(s, pb + (w_ext + E.at(i, q)));
        }
      }
    }

    // Qg (g)(h)
    int q = j + 1;
    if (q < n) {
      for (Grp& g : Qgrs[j].grps) {
        int e = g.v;
        for (Ent& s : g.ents) {
          int c = s.p, d = s.u;
          int i_lo = lead_i_lo(c);
          for (int i = i_lo; i < c; ++i) {
            if (!E.pairs(i, q) || !E.end_ok(i, q)) continue;
            double w = E.w_ap + E.w_bp + E.w_cp * (c - i - 1) + E.at(i, q);
            if (dz) w += E.run(i + 1, c - 1);
            double pb = g4(Qg[q], e, i, d);
            add4(s, pb + w);
          }
          if (c - 1 >= 0) {
            const OMap& at = qmp_at(c - 1);
            qm_pair_slots(c - 1, q, [&](int r) {
              int ii = at.keys[r];
              if (ii < 0) return;
              double qm = at.vals[r];
              int outer = ii - 1;
              if (!E.pairs(outer, q) || outer >= c) return;
              if (!E.end_ok(outer, q)) return;
              double at_o = E.at(outer, q);
              double w = qm + E.w_ap + E.w_bp + at_o;
              double pb = g4(Qg[q], e, outer, d);
              add4(s, pb + w);
              add2b(qmp_b(), ii, c - 1, pb + s.a + E.w_ap + E.w_bp + at_o);
            });
          }
        }
      }
    }
  }

  void outside_same_j(int j) {
    const double w2 = E.w_alpha2;
    const double w3 = E.w_alpha3;
    const double a1 = E.w_alpha1;

    // C = C_left × Qb / Qp
    double cb = C_b[j];
    {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        int i = at.keys[r];
        if (i < 0 || !E.end_ok(i, j)) continue;
        double qb = at.vals[r];
        double at_ij = E.at(i, j);
        double left_a = (i == 0) ? 0.0 : C[i - 1];
        add2b(Qb_b, i, j, cb + left_a + at_ij);
        if (i > 0) C_b[i - 1] = log_add(C_b[i - 1], cb + qb + at_ij);
      }
    }
    if (allow_pk) {
      const OMap& at = Qp_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        int i = at.keys[r];
        if (i < 0) continue;
        double qp = at.vals[r];
        double left_a = (i == 0) ? 0.0 : C[i - 1];
        add2b(Qp_b, i, j, cb + left_a + E.w_beta1);
        qps(0, i, j, cb + left_a + E.w_beta1);
        if (i > 0) C_b[i - 1] = log_add(C_b[i - 1], cb + qp + E.w_beta1);
      }
    }

    // Qgrs = Qg[f] × Qm[f+1, j]
    if (allow_pk) {
      const OMap& at = qmp_at(j);
      SpMat& mpb = qmp_b();
      for (size_t r = 0; r < at.slots(); ++r) {
        int p = at.keys[r];
        if (p < 0) continue;
        double qm = at.vals[r];
        int f = p - 1;
        if (f < 0) continue;
        for (Grp& g : Qg[f].grps) {
          int e = g.v;
          for (Ent& s : g.ents) {
            if (s.p == s.u || !E.end_ok(s.p, f)) continue;
            double pb = g4(Qgrs[j], e, s.p, s.u);
            double w = E.w_bp + E.at(s.p, f);
            add4(s, pb + qm + w);
            add2b(mpb, p, j, pb + s.a + w);
          }
        }
        for (int c : partners_left[f]) {
          if (!E.end_ok(c, f)) continue;
          double pb = g4(Qgrs[j], f, c, c);
          add2b(mpb, p, j, pb + (E.w_bp + E.at(c, f)));
        }
      }
    }

    if (dz) {
      if (keep_closed_b) {
        for (size_t r = 0; r < Qm_at[j].slots(); ++r)
          if (Qm_at[j].keys[r] >= 0) QmC_b[ix(Qm_at[j].keys[r], j)] = Qm_b[ix(Qm_at[j].keys[r], j)];
        for (size_t r = 0; r < Qm2_at[j].slots(); ++r)
          if (Qm2_at[j].keys[r] >= 0)
            Qm2C_b[ix(Qm2_at[j].keys[r], j)] = Qm2_b[ix(Qm2_at[j].keys[r], j)];
      }
      run_back_at(Qm_b, Qm_at[j], j, w3, rbQm, RS_QM);
      run_back_at(Qm2_b, Qm2_at[j], j, w3, rbQm2, RS_QM2);
      if (allow_pk && E.stem_ml) run_back_at(Qmp_b, Qmp_at[j], j, E.w_cp, rbQmp, RS_QMP);
    }

    // Pass 1 E3 lands on the branch-ending outside, after the run split.
    lp_out_e3(j);

    // Qm / Qm2 from Qb, Qp at this j
    auto rev_branch = [&](int d, double obj_a, SpMat& obj_beta, bool pk,
                          double w_at) {
      double w_single_base = (pk ? E.w_beta1m + E.ml_pk_k * w2 : w2) + w_at;
      const OMap& amj = Qm_at[j];
      for (size_t r = 0; r < amj.slots(); ++r) {
        int i = amj.keys[r];
        if (i < 0 || i > d || d - i > E.ml_left || !lead_gap_ok(d - i)) continue;
        double w = w_single_base;
        if (l2r_mode == 0 && !allow_pk) {
          if (i != d) continue;
        } else {
          w += w3 * (d - i);
          if (dz) w += E.run(i, d - 1);
        }
        double mb = Qm_b[ix(i, j)];
        add2b(obj_beta, d, j, mb + w);
        if (pk) qps(1, d, j, mb + w);
      }
      if (d == 0) return;
      double w_more = w_single_base;
      const OMap& amd = Qm_at[d - 1];
      for (size_t r = 0; r < amd.slots(); ++r) {
        int i = amd.keys[r];
        if (i < 0) continue;
        double qm = amd.vals[r];
        double mb = Qm_b[ix(i, j)];
        add2b(obj_beta, d, j, mb + qm + w_more);
        if (pk) qps(1, d, j, mb + qm + w_more);
        add2b(Qm_b, i, d - 1, mb + obj_a + w_more);
        if (!pk) {
          double m2b = Qm2_b[ix(i, j)];
          add2b(obj_beta, d, j, m2b + qm + w_more);
          add2b(Qm_b, i, d - 1, m2b + obj_a + w_more);
        }
      }
    };
    {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        int d = at.keys[r];
        if (d < 0 || !E.end_ok(d, j)) continue;
        rev_branch(d, at.vals[r], Qb_b, false, E.at(d, j));
      }
    }
    if (allow_pk) {
      const OMap& at = Qp_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        int d = at.keys[r];
        if (d < 0) continue;
        rev_branch(d, at.vals[r], Qp_b, true, 0.0);
      }
    }

    // Qmp from Qb, Qp at this j
    if (allow_pk && E.stem_ml) {
      const double bp = E.w_bp, cp = E.w_cp;
      auto rev_branch_p = [&](int d, double obj_a, SpMat& obj_beta, bool pk,
                              double w_at) {
        double w_single = (pk ? E.w_beta1m + E.ml_pk_k * bp : bp) + w_at;
        const OMap& amj = Qmp_at[j];
        for (size_t r = 0; r < amj.slots(); ++r) {
          int i = amj.keys[r];
          if (i < 0 || i > d || !lead_gap_ok(d - i)) continue;
          double mb = Qmp_b[ix(i, j)];
          if (dz) {
            add2b(obj_beta, d, j, mb + w_single + cp * (d - i) + E.run(i, d - 1));
            if (pk) qps(1, d, j, mb + w_single + cp * (d - i) + E.run(i, d - 1));
            continue;
          }
          add2b(obj_beta, d, j, mb + w_single + cp * (d - i));
          if (pk) qps(1, d, j, mb + w_single + cp * (d - i));
        }
        if (d == 0) return;
        const OMap& amd = Qmp_at[d - 1];
        for (size_t r = 0; r < amd.slots(); ++r) {
          int i = amd.keys[r];
          if (i < 0) continue;
          double qm = amd.vals[r];
          double mb = Qmp_b[ix(i, j)];
          add2b(obj_beta, d, j, mb + qm + w_single);
          if (pk) qps(1, d, j, mb + qm + w_single);
          add2b(Qmp_b, i, d - 1, mb + obj_a + w_single);
        }
      };
      const OMap& atb = Qb_at[j];
      for (size_t r = 0; r < atb.slots(); ++r) {
        int d = atb.keys[r];
        if (d < 0 || !E.end_ok(d, j)) continue;
        rev_branch_p(d, atb.vals[r], Qb_b, false, E.at(d, j));
      }
      const OMap& atp = Qp_at[j];
      for (size_t r = 0; r < atp.slots(); ++r) {
        int d = atp.keys[r];
        if (d < 0) continue;
        rev_branch_p(d, atp.vals[r], Qp_b, true, 0.0);
      }
    }

    // Qz[i,j] += Qz_left × Qb/Qp
    for (int i = 0; i <= j; ++i) {
      double zb = Qz_b[ix(i, j)];
      if (zb < HALF_NEG) continue;
      const OMap& atb = Qb_at[j];
      for (size_t r = 0; r < atb.slots(); ++r) {
        int d = atb.keys[r];
        if (d < 0 || d < i || !E.end_ok(d, j)) continue;
        double qb = atb.vals[r];
        double w = E.w_nest + E.at(d, j);
        double left_a = (d == i) ? 0.0 : Qz[ix(i, d - 1)];
        add2b(Qb_b, d, j, zb + left_a + w);
        if (d > i) add2b(Qz_b, i, d - 1, zb + qb + w);
      }
      if (allow_pk) {
        const OMap& atp = Qp_at[j];
        for (size_t r = 0; r < atp.slots(); ++r) {
          int d = atp.keys[r];
          if (d < 0 || d < i) continue;
          double qp = atp.vals[r];
          double left_a = (d == i) ? 0.0 : Qz[ix(i, d - 1)];
          double w = E.w_beta1p + E.nest_pk_k * E.w_nest;
          add2b(Qp_b, d, j, zb + left_a + w);
          qps(2, d, j, zb + left_a + w);
          if (d > i) add2b(Qz_b, i, d - 1, zb + qp + w);
        }
      }
    }

    // Qb[i,j] += Qm2[i+1, j-1] × (α1+α2), plus the PK-inside terms from PKM[j-1]
    const double w_ml = a1 + E.w_alpha2;
    const double w_left = a1 + E.w_beta1m + (1.0 + E.ml_pk_k) * E.w_alpha2;
    {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        int i = at.keys[r];
        if (i < 0 || !E.end_ok(i, j)) continue;
        double at_ij = E.at_close(i, j);
        double bb = Qb_b[ix(i, j)];
        if (l2r_mode != 0 && i + 1 <= j - 1) add2b(Qm2_b, i + 1, j - 1, bb + (w_ml + at_ij));
        if (bb < HALF_NEG || !allow_pk || j < 1) continue;
        double wl = w_left + at_ij;

        const OMap& pkm = PKM_at[j - 1];
        for (size_t s = 0; s < pkm.slots(); ++s) {
          int d = pkm.keys[s];
          if (d <= i) continue;
          double pk = pkm.vals[s];
          const size_t pkix = ix(j - 1, d);
          double v = bb + wl + w3 * (d - i - 1);
          if (dz) v += E.run(i + 1, d - 1);
          PKM_b[pkix] = log_add(static_cast<double>(PKM_b[pkix]), v);
          if (i + 1 <= d - 1) {
            double left = Qm[ix(i + 1, d - 1)];
            if (left > HALF_NEG) {
              add2b(Qm_b, i + 1, d - 1, bb + pk + wl);
              PKM_b[pkix] = log_add(static_cast<double>(PKM_b[pkix]), bb + left + wl);
            }
          }
        }
      }
    }

    lp_out_multi(j);

    if (!allow_pk) return;
    if (l2r_mode == 1) l2r_back_pk(j);

    // Qp = Qgl × Qgr
    for (Grp& g : Qgr[j].grps) {
      int e_gr = g.v;
      int f = e_gr - 1;
      if (f < 0) continue;
      Tab4& gl_right = Qgl[f];
      for (Ent& s : g.ents) {
        int d = s.p, em1 = s.u;
        double gr = s.a;
        if (!E.end_ok(d, j)) continue;
        int e = em1 + 1;
        Grp* gl_group = gl_right.get(e);
        if (!gl_group) continue;
        for (Ent& t : gl_group->ents) {
          int i = t.p, uu = t.u;
          if (uu != d - 1 || i > d - 1 || !E.end_ok(i, f)) continue;
          double pb = Qp_b[ix(i, j)] + (E.at(i, f) + E.at(d, j));
          add4(s, pb + t.a);
          add4(t, pb + gr);
        }
      }
    }

    // Qgr = Qgl × Qz
    for (Grp& g : Qgl[j].grps) {
      int f = g.v;
      for (Ent& s : g.ents) {
        int i = s.p, d = s.u;
        double gl = s.a;
        for (int e = d + 1; e <= f; ++e) {
          double qz = qz_log(e, f - 1);
          if (qz < HALF_NEG) continue;
          double pb = g4(Qgr[j], e, i, d);
          add4(s, pb + qz);
          if (e <= f - 1) add2b(Qz_b, e, f - 1, pb + gl);
        }
      }
    }

    // Qgl = Qg × Qz × β2
    for (Grp& g : Qg[j].grps) {
      int f = g.v;
      for (Ent& s : g.ents) {
        int i = s.p, d = s.u;
        double qg = s.a;
        if (!E.end_ok(d, f)) continue;
        double w = E.w_stem + E.at(d, f);
        for (int e = d; e < f; ++e) {
          double qz = qz_log(d + 1, e);
          if (qz < HALF_NEG) continue;
          double pb = g4(Qgl[j], f, i, e);
          add4(s, pb + qz + w);
          if (d + 1 <= e) add2b(Qz_b, d + 1, e, pb + qg + w);
        }
      }
    }

    // Qgls = Qm × Qg
    SpMat& mpb = qmp_b();
    for (Grp& g : Qg[j].grps) {
      int e = g.v;
      for (Ent& s : g.ents) {
        int c = s.p, d = s.u;
        double qg = s.a;
        if (c - 1 < 0 || c == d || !E.end_ok(c, j)) continue;
        double w = E.w_bp + E.at(c, j);
        const OMap& at = qmp_at(c - 1);
        for (size_t r = 0; r < at.slots(); ++r) {
          int i = at.keys[r];
          if (i < 0) continue;
          double qm = at.vals[r];
          double pb = g4(Qgls[j], e, i, d);
          add4(s, pb + qm + w);
          add2b(mpb, i, c - 1, pb + qg + w);
        }
      }
    }
    for (int c : partners_left[j]) {
      if (c - 1 < 0 || !E.end_ok(c, j)) continue;
      double w = E.w_bp + E.at(c, j);
      const OMap& at = qmp_at(c - 1);
      for (size_t r = 0; r < at.slots(); ++r) {
        int i = at.keys[r];
        if (i < 0) continue;
        double pb = g4(Qgls[j], j, i, c);
        add2b(mpb, i, c - 1, pb + w);
      }
    }
  }

  // --------------------------------------------------------------- bpp
  std::vector<std::tuple<int, int, double>> collect_bpp(double cutoff) {
    std::vector<std::tuple<int, int, double>> out;
    if (n == 0 || logZ < HALF_NEG) return out;
    std::vector<double> logP(static_cast<size_t>(n) * n, NEG_INF);
    std::vector<char> seen(static_cast<size_t>(n) * n, 0);
    auto add_pair = [&](int i, int j, double lp) {
      if (lp < HALF_NEG || i < 0 || j <= i) return;
      size_t k = ix(i, j);
      logP[k] = log_add(seen[k] ? logP[k] : NEG_INF, lp);
      seen[k] = 1;
    };

    for (int j = 0; j < n; ++j) {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        int i = at.keys[r];
        if (i < 0) continue;
        add_pair(i, j, at.vals[r] + Qb_b[ix(i, j)] - logZ);
      }
    }

    if (allow_pk) {
      for (int j = 0; j < n; ++j)
        for (Grp& g : Qg[j].grps)
          for (Ent& s : g.ents) add_pair(s.p, j, s.a + s.b - logZ);
      collect_qg_inner_pairs(add_pair);
      if (l2r_mode == 1) l2r_pairs(add_pair);
    }

    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j) {
        size_t k = ix(i, j);
        if (!seen[k]) continue;
        double p = std::exp(std::min(0.0, logP[k]));
        if (p >= cutoff) out.emplace_back(i + 1, j + 1, std::min(1.0, p));
      }
    return out;
  }

  template <class AddPair>
  void collect_qg_inner_pairs(AddPair& add_pair) {
    for (int j = 0; j < n; ++j) {
      const auto& pair_i = partners_left[j];
      // (a)
      for (int i : pair_i) {
        int d_hi = std::min(j - 1, i + 1 + cap);
        for (int d = i + 1; d <= d_hi; ++d) {
          int e_lo = std::max(d + MIN_PAIR, j - 1 - cap);
          for (int e = e_lo; e < j; ++e) {
            double w = E.stem_interior(i, d, e, j);
            if (w < HALF_NEG) continue;
            add_pair(d, e, w + g4(Qg[j], e, i, d) - logZ);
          }
        }
      }
      // seed c·j with Qm on its left (Qgls group j)
      if (Grp* g = Qgls[j].get(j))
        for (const Ent& s : g->ents) add_pair(s.u, j, s.a + s.b - logZ);
      // seed c·e with Qm on its right (Qgrs entries (c, c))
      for (const Grp& g : Qgrs[j].grps)
        for (const Ent& s : g.ents)
          if (s.p == s.u) add_pair(s.p, g.v, s.a + s.b - logZ);
    }
  }
