// Exact DP03 pseudoknot tables (the reference, --pk-l2r off at --beam inf): fill_Qg, Qgls,
// Qgl, Qgr, Qp, Qgrs (4-index recursions of Dirks-Pierce 2003).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  // ------------------------------------------------------------------- Qg
  void fill_Qg(int j) {
    Tab4& G = Qg[j];
    const auto& pair_i = partners_left[j];

    // (a) interior loop between i·j and d·e
    for (int i : pair_i) {
      int d_hi = std::min(j - 1, i + 1 + cap);
      for (int d = i + 1; d <= d_hi; ++d) {
        int e_lo = std::max(d + MIN_PAIR, j - 1 - cap);
        for (int e = e_lo; e < j; ++e) G.acc(e, i, d, E.stem_interior(i, d, e, j) + 2.0 * pair_bonus);
      }
    }

    // (b) interior wrapping a stored Qg_{c,d,e,f}
    int f_lo = std::max(0, j - 1 - cap);
    for (int f = f_lo; f < j; ++f) {
      for (const Grp& g : Qg[f].grps) {
        int e = g.v;
        for (const Ent& s : g.ents) {
          int c = s.p, d = s.u;
          if (c == d) continue;
          double x = s.a;
          int i_lo = std::max(0, c - 1 - cap);
          for (int i = i_lo; i < c; ++i) {
            double w = E.stem_interior(i, c, f, j);
            G.acc(e, i, d, w + x + pair_bonus);
          }
        }
      }
    }

    // (f) Qgls plus unpaired run to j; seed-derived Qgls (e == f) are DP03 (c), run capped
    for (int f = 0; f < j; ++f) {
      int unpaired = j - f - 1;
      double w_ext = E.w_ap + E.w_bp + E.w_cp * unpaired + pair_bonus;
      if (dz) w_ext += E.run(f + 1, j - 1);
      for (const Grp& g : Qgls[f].grps) {
        int e = g.v;
        if (e == f && unpaired > cap) continue;
        for (const Ent& s : g.ents) {
          int i = s.p - 1;
          if (E.pairs(i, j) && E.end_ok(i, j)) G.acc(e, i, s.u, s.a + (w_ext + E.at(i, j)));
        }
      }
    }

    // (g) Qgrs at j−1, unpaired on the left
    if (j >= 1) {
      for (const Grp& g : Qgrs[j - 1].grps) {
        int e = g.v;
        for (const Ent& s : g.ents) {
          int c = s.p, d = s.u;
          double x = s.a;
          int i_lo = lead_i_lo(c);
          for (int i = i_lo; i < c; ++i) {
            if (!E.pairs(i, j) || !E.end_ok(i, j)) continue;
            double w = x + E.w_ap + E.w_bp + E.w_cp * (c - i - 1) + E.at(i, j) + pair_bonus;
            if (dz) w += E.run(i + 1, c - 1);
            G.acc(e, i, d, w);
          }
        }
      }
    }

    // (h) Qgrs + Qm on the left
    if (j >= 1) {
      for (const Grp& g : Qgrs[j - 1].grps) {
        int e = g.v;
        for (const Ent& s : g.ents) {
          int c = s.p, d = s.u;
          double x = s.a;
          if (c - 1 < 0) continue;
          const OMap& at = qmp_at(c - 1);
          qm_pair_slots(c - 1, j, [&](int r) {
            int i = at.keys[r];
            if (i < 0) return;
            double qm = at.vals[r];
            int outer = i - 1;
            if (E.pairs(outer, j) && outer < c && E.end_ok(outer, j))
              G.acc(e, outer, d, qm + x + E.w_ap + E.w_bp + E.at(outer, j) + pair_bonus);
          });
        }
      }
    }

    // (0) the pair i·j alone as a stem, Qg_{i,i,j,j} = 1
    if (single_stems) {
      for (int i : pair_i)
        if (E.end_ok(i, j)) G.acc(j, i, i, pair_bonus);
    }
  }

  void fill_Qgls(int j) {
    Tab4& out = Qgls[j];
    for (const Grp& g : Qg[j].grps) {
      int e = g.v;
      for (const Ent& s : g.ents) {
        int c = s.p, d = s.u;
        double qg = s.a;
        if (c - 1 < 0 || c == d || !E.end_ok(c, j)) continue;
        double w = E.w_bp + E.at(c, j);
        const OMap& at = qmp_at(c - 1);
        for (size_t r = 0; r < at.slots(); ++r) {
          int i = at.keys[r];
          if (i < 0) continue;
          out.acc(e, i, d, w + at.vals[r] + qg);
        }
      }
    }
    // seeds c·j (DP03 (c) once closed by (f))
    for (int c : partners_left[j]) {
      if (c - 1 < 0 || !E.end_ok(c, j)) continue;
      double w = E.w_bp + E.at(c, j) + pair_bonus;
      const OMap& at = qmp_at(c - 1);
      for (size_t r = 0; r < at.slots(); ++r) {
        int i = at.keys[r];
        if (i < 0) continue;
        out.acc(j, i, c, w + at.vals[r]);
      }
    }
  }

  void fill_Qgl(int j) {
    Tab4& out = Qgl[j];
    for (const Grp& g : Qg[j].grps) {
      int f = g.v;
      for (const Ent& s : g.ents) {
        int i = s.p, d = s.u;
        double qg = s.a;
        if (!E.end_ok(d, f)) continue;
        double w = E.w_stem + E.at(d, f);
        for (int e = d; e < f; ++e) {
          double qz = qz_log(d + 1, e);
          if (qz < HALF_NEG) continue;
          out.acc(f, i, e, qg + qz + w);
        }
      }
    }
  }

  void fill_Qgr(int j) {
    Tab4& out = Qgr[j];
    for (const Grp& g : Qgl[j].grps) {
      int f = g.v;
      for (const Ent& s : g.ents) {
        int i = s.p, d = s.u;
        double gl = s.a;
        for (int e = d + 1; e <= f; ++e) {
          double qz = qz_log(e, f - 1);
          if (qz < HALF_NEG) continue;
          out.acc(e, i, d, gl + qz);
        }
      }
    }
  }

  void fill_Qp(int j) {
    for (const Grp& g : Qgr[j].grps) {
      int e_gr = g.v;
      int f = e_gr - 1;
      if (f < 0) continue;
      Tab4& gl_right = Qgl[f];
      for (const Ent& s : g.ents) {
        int d = s.p, em1 = s.u;
        double gr = s.a;
        if (!E.end_ok(d, j)) continue;
        int e = em1 + 1;
        Grp* gl_group = gl_right.get(e);
        if (!gl_group || gl_group->ents.empty()) continue;
        for (const Ent& t : gl_group->ents) {
          int i = t.p, uu = t.u;
          if (uu != d - 1 || i > d - 1 || !E.end_ok(i, f)) continue;
          set_pair(Qp, Qp_at, i, j, t.a + gr + (E.at(i, f) + E.at(d, j)));
        }
      }
    }
  }

  void fill_Qgrs(int j) {
    Tab4& out = Qgrs[j];
    const OMap& at = qmp_at(j);
    for (size_t r = 0; r < at.slots(); ++r) {
      int p = at.keys[r];
      if (p < 0) continue;
      double qm = at.vals[r];
      int f = p - 1;
      if (f < 0) continue;
      for (const Grp& g : Qg[f].grps) {
        int e = g.v;
        for (const Ent& s : g.ents) {
          if (s.p == s.u || !E.end_ok(s.p, f)) continue;
          out.acc(e, s.p, s.u, s.a + qm + (E.w_bp + E.at(s.p, f)));
        }
      }
      // seeds c·f (DP03 (d), (e) once closed by (g), (h))
      for (int c : partners_left[f])
        if (E.end_ok(c, f)) out.acc(f, c, c, qm + (E.w_bp + E.at(c, f)) + pair_bonus);
    }
  }

