// Stem-estimate diagnostics (out of class): tagged_outside, rank statistics, H_bar and
// hybrid outside (pr_build_hbar, pr_hybrid_outside), --pr-oracle, --stem-diag.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.


// Stem-role outside of the 4-index tables: the 4-index part of outside() with Qp_b fixed, where
// the Qp = Qgl x Qgr edge feeds only Qgl[f] (tag 1: the state ends up as the left stem, stem 1)
// or only Qgr[j] (tag 2: right stem).  Overwrites Ent::b of all 4-index tables.
void DP03::tagged_outside(int tag) {
  for (std::vector<Tab4>* tabs : {&Qg, &Qgls, &Qgl, &Qgr, &Qgrs})
    for (Tab4& t : *tabs)
      for (Grp& g : t.grps)
        for (Ent& e : g.ents) e.b = NEG_INF;
  for (int j = n - 1; j >= 0; --j) {
    const int jmax = std::min(n, j + 2 + cap);
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
            add4(s, g4(Qg[q], e, i, d) + w);
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
          add4(s, g4(Qg[q], e, i, s.u) + (w_ext + E.at(i, q)));
        }
      }
    }
    // Qg (g)(h)
    if (j + 1 < n) {
      const int q = j + 1;
      for (Grp& g : Qgrs[j].grps) {
        int e = g.v;
        for (Ent& s : g.ents) {
          int c = s.p, d = s.u;
          for (int i = lead_i_lo(c); i < c; ++i) {
            if (!E.pairs(i, q) || !E.end_ok(i, q)) continue;
            double w = E.w_ap + E.w_bp + E.w_cp * (c - i - 1) + E.at(i, q);
            if (dz) w += E.run(i + 1, c - 1);
            add4(s, g4(Qg[q], e, i, d) + w);
          }
          if (c - 1 >= 0) {
            const OMap& at = qmp_at(c - 1);
            qm_pair_slots(c - 1, q, [&](int r) {
              int ii = at.keys[r];
              if (ii < 0) return;
              int outer = ii - 1;
              if (!E.pairs(outer, q) || outer >= c || !E.end_ok(outer, q)) return;
              double w = at.vals[r] + E.w_ap + E.w_bp + E.at(outer, q);
              add4(s, g4(Qg[q], e, outer, d) + w);
            });
          }
        }
      }
    }
    // Qgrs = Qg[f] x Qm[f+1, j]
    {
      const OMap& at = qmp_at(j);
      for (size_t r = 0; r < at.slots(); ++r) {
        int p = at.keys[r];
        if (p < 0 || p - 1 < 0) continue;
        double qm = at.vals[r];
        for (Grp& g : Qg[p - 1].grps) {
          int e = g.v;
          for (Ent& s : g.ents) {
            if (s.p == s.u || !E.end_ok(s.p, p - 1)) continue;
            add4(s, g4(Qgrs[j], e, s.p, s.u) + qm + (E.w_bp + E.at(s.p, p - 1)));
          }
        }
      }
    }
    // Qp = Qgl x Qgr
    for (Grp& g : Qgr[j].grps) {
      int f = g.v - 1;
      if (f < 0) continue;
      for (Ent& s : g.ents) {
        int d = s.p;
        if (!E.end_ok(d, j)) continue;
        Grp* gl_group = Qgl[f].get(s.u + 1);
        if (!gl_group) continue;
        for (Ent& t : gl_group->ents) {
          int i = t.p;
          if (t.u != d - 1 || i > d - 1 || !E.end_ok(i, f)) continue;
          double pb = Qp_b[ix(i, j)] + (E.at(i, f) + E.at(d, j));
          if (tag == 2) add4(s, pb + t.a);
          else add4(t, pb + s.a);
        }
      }
    }
    // Qgr = Qgl x Qz
    for (Grp& g : Qgl[j].grps) {
      int f = g.v;
      for (Ent& s : g.ents)
        for (int e = s.u + 1; e <= f; ++e) {
          double qz = qz_log(e, f - 1);
          if (qz < HALF_NEG) continue;
          add4(s, g4(Qgr[j], e, s.p, s.u) + qz);
        }
    }
    // Qgl = Qg x Qz x beta2
    for (Grp& g : Qg[j].grps) {
      int f = g.v;
      for (Ent& s : g.ents) {
        int d = s.u;
        if (!E.end_ok(d, f)) continue;
        double w = E.w_stem + E.at(d, f);
        for (int e = d; e < f; ++e) {
          double qz = qz_log(d + 1, e);
          if (qz < HALF_NEG) continue;
          add4(s, g4(Qgl[j], f, s.p, e) + qz + w);
        }
      }
    }
    // Qgls = Qm x Qg
    for (Grp& g : Qg[j].grps) {
      int e = g.v;
      for (Ent& s : g.ents) {
        int c = s.p, d = s.u;
        if (c - 1 < 0 || c == d || !E.end_ok(c, j)) continue;
        double w = E.w_bp + E.at(c, j);
        const OMap& at = qmp_at(c - 1);
        for (size_t r = 0; r < at.slots(); ++r) {
          int i = at.keys[r];
          if (i < 0) continue;
          add4(s, g4(Qgls[j], e, i, d) + at.vals[r] + w);
        }
      }
    }
  }
}

namespace stemdiag {

// Ranks (0-based, ties averaged) of x.
inline void avg_ranks(const std::vector<double>& x, std::vector<double>& r) {
  const size_t m = x.size();
  std::vector<size_t> o(m);
  for (size_t k = 0; k < m; ++k) o[k] = k;
  std::stable_sort(o.begin(), o.end(), [&](size_t a, size_t b) { return x[a] < x[b]; });
  r.assign(m, 0.0);
  for (size_t i = 0; i < m;) {
    size_t j = i;
    while (j + 1 < m && x[o[j + 1]] == x[o[i]]) ++j;
    double rk = 0.5 * static_cast<double>(i + j);
    for (size_t k = i; k <= j; ++k) r[o[k]] = rk;
    i = j + 1;
  }
}

// Pearson correlation; NaN if either side is constant.
inline double pearson(const std::vector<double>& a, const std::vector<double>& b) {
  const size_t m = a.size();
  double ma = 0, mb = 0;
  for (size_t k = 0; k < m; ++k) ma += a[k], mb += b[k];
  ma /= m;
  mb /= m;
  double sab = 0, saa = 0, sbb = 0;
  for (size_t k = 0; k < m; ++k) {
    sab += (a[k] - ma) * (b[k] - mb);
    saa += (a[k] - ma) * (a[k] - ma);
    sbb += (b[k] - mb) * (b[k] - mb);
  }
  if (saa <= 0 || sbb <= 0) return std::nan("");
  return sab / std::sqrt(saa * sbb);
}

// Top b of s (threshold rule, ties to the earliest), as DP03::keep_mask.
inline void top_b(const std::vector<double>& s, size_t b, std::vector<char>& keep) {
  keep.assign(s.size(), 1);
  if (s.size() <= b) return;
  std::vector<double> tmp(s);
  std::nth_element(tmp.begin(), tmp.begin() + (b - 1), tmp.end(), std::greater<double>());
  const double thr = tmp[b - 1];
  long room = static_cast<long>(b);
  for (double x : s) room -= x > thr;
  for (size_t k = 0; k < s.size(); ++k) {
    if (s[k] > thr) keep[k] = 1;
    else if (s[k] == thr && room > 0) keep[k] = 1, --room;
    else keep[k] = 0;
  }
}

inline void num(std::string& o, double x) {
  char buf[40];
  if (std::isnan(x)) std::snprintf(buf, sizeof buf, "null");
  else std::snprintf(buf, sizeof buf, "%.10g", x);
  o += buf;
}

}  // namespace stemdiag

// --pr-oracle, nested partner stems: fill_Qg on the kept pairs km (kl[j]: kept i < j) with
// interior loops <= lc per side, multiloop steps (f, g, h) if ml with the multiloop interiors
// mb (dense, both sides), seeds and run caps as the engine.
void DP03::pr_build_hbar(const std::vector<char>& km, const std::vector<std::vector<int>>& kl,
                         const std::vector<double>& mb, bool ml, int lc, std::vector<Tab4>& hb) {
  hb.assign(n, Tab4());
  std::vector<Tab4> ls(ml ? n : 0);
  auto M = [&](int i, int j) { return mb[ix(i, j)]; };
  for (int j = 0; j < n; ++j) {
    Tab4& G = hb[j];
    for (int i : kl[j]) {
      const int d_hi = std::min(j - 1, i + 1 + lc);
      for (int d = i + 1; d <= d_hi; ++d)
        for (int e = std::max(d + MIN_PAIR, j - 1 - lc); e < j; ++e)
          if (km[ix(d, e)]) G.acc(e, i, d, E.stem_interior(i, d, e, j) + 2.0 * pair_bonus);
    }
    for (int f = std::max(0, j - 1 - lc); f < j; ++f)
      for (const Grp& g : hb[f].grps)
        for (const Ent& s : g.ents) {
          if (s.p == s.u) continue;
          for (int i = std::max(0, s.p - 1 - lc); i < s.p; ++i)
            if (km[ix(i, j)]) G.acc(g.v, i, s.u, E.stem_interior(i, s.p, f, j) + s.a + pair_bonus);
        }
    if (ml) {
      for (int f = 0; f < j; ++f) {
        const int unpaired = j - f - 1;
        double w_ext = E.w_ap + E.w_bp + E.w_cp * unpaired + pair_bonus;
        if (dz) w_ext += E.run(f + 1, j - 1);
        for (const Grp& g : ls[f].grps) {
          if (g.v == f && unpaired > cap) continue;
          for (const Ent& s : g.ents) {
            const int i = s.p - 1;
            if (i >= 0 && km[ix(i, j)] && E.end_ok(i, j))
              G.acc(g.v, i, s.u, s.a + (w_ext + E.at(i, j)));
          }
        }
      }
      if (j >= 1) {
        Tab4 rs;
        const int jj = j - 1;
        for (int p = 1; p <= jj; ++p) {
          const double qm = M(p, jj);
          if (qm < HALF_NEG) continue;
          const int f = p - 1;
          for (const Grp& g : hb[f].grps)
            for (const Ent& s : g.ents)
              if (s.p != s.u && E.end_ok(s.p, f))
                rs.acc(g.v, s.p, s.u, s.a + qm + (E.w_bp + E.at(s.p, f)));
          for (int c : kl[f])
            if (E.end_ok(c, f)) rs.acc(f, c, c, qm + (E.w_bp + E.at(c, f)) + pair_bonus);
        }
        for (const Grp& g : rs.grps)
          for (const Ent& s : g.ents) {
            const int c = s.p, d = s.u;
            for (int i = lead_i_lo(c); i < c; ++i) {
              if (!km[ix(i, j)] || !E.end_ok(i, j)) continue;
              double w = s.a + E.w_ap + E.w_bp + E.w_cp * (c - i - 1) + E.at(i, j) + pair_bonus;
              if (dz) w += E.run(i + 1, c - 1);
              G.acc(g.v, i, d, w);
            }
            for (int i = 1; i < c; ++i) {
              const double qm = M(i, c - 1);
              if (qm < HALF_NEG) continue;
              const int o = i - 1;
              if (km[ix(o, j)] && E.end_ok(o, j))
                G.acc(g.v, o, d, qm + s.a + E.w_ap + E.w_bp + E.at(o, j) + pair_bonus);
            }
          }
      }
    }
    if (single_stems)
      for (int i : kl[j])
        if (E.end_ok(i, j)) G.acc(j, i, i, pair_bonus);
    if (ml) {
      Tab4& out = ls[j];
      for (const Grp& g : G.grps)
        for (const Ent& s : g.ents) {
          const int c = s.p;
          if (c - 1 < 1 || c == s.u || !E.end_ok(c, j)) continue;
          const double w = E.w_bp + E.at(c, j) + s.a;
          for (int i = 1; i <= c - 1; ++i) {
            const double qm = M(i, c - 1);
            if (qm > HALF_NEG) out.acc(g.v, i, s.u, w + qm);
          }
        }
      for (int c : kl[j]) {
        if (c - 1 < 1 || !E.end_ok(c, j)) continue;
        const double w = E.w_bp + E.at(c, j) + pair_bonus;
        for (int i = 1; i <= c - 1; ++i) {
          const double qm = M(i, c - 1);
          if (qm > HALF_NEG) out.acc(j, i, c, w + qm);
        }
      }
    }
  }
}

// --pr-oracle: the principled Out_hat as a 4-index outside pass (tagged_outside's order).
// Qg: growth onto kept pairs (interior <= lc per side; single-pair states too, as the
// inner=true roots), Qgls / Qgrs closings onto kept pairs (if ml), the stop through Qgl.
// Qgl's L role: the partner's Qgr from exact Qgl (hx) or from hb (nested H), times the future
// filler Qz (zx) or Qzbar; R role: the exact Qgl alpha of the partner.  Ctx: Qp_b (cx) or
// pr_ctx.  The right multiloop branch of a growth step: exact Qm (mx) or Mbar.  Ent.b holds
// the result; own[0] / own[1]: Qgls / Qgrs over every closing (their Out_hat), flat in
// off[1] / off[4] order.
void DP03::pr_hybrid_outside(const PrVar& v, const StemEst& st, const std::vector<char>& km,
                             const std::vector<std::vector<int>>& kr, std::vector<Tab4>* hb,
                             std::vector<double> own[2], const std::vector<size_t> off[5]) {
  for (std::vector<Tab4>* tabs : {&Qg, &Qgls, &Qgl, &Qgr, &Qgrs})
    for (Tab4& t : *tabs)
      for (Grp& g : t.grps)
        for (Ent& e : g.ents) e.b = NEG_INF;
  std::fill(own[0].begin(), own[0].end(), NEG_INF);
  std::fill(own[1].begin(), own[1].end(), NEG_INF);
  const int lc = std::min(v.lc, cap);
  auto ctx = [&](int i, int j) { return v.cx ? Qp_b[ix(i, j)] : pr_ctx(i, j); };
  Tab4 gl, gr;
  for (int j = n - 1; j >= 0; --j) {
    // Qg growth
    for (Grp& g : Qg[j].grps) {
      const int e = g.v;
      for (Ent& s : g.ents) {
        const int c = s.p, d = s.u;
        for (int i = std::max(0, c - 1 - lc); i < c; ++i) {
          const std::vector<int>& r = kr[i];
          for (auto it = std::upper_bound(r.begin(), r.end(), j);
               it != r.end() && *it <= j + 1 + lc; ++it) {
            const double w = E.stem_interior(i, c, j, *it);
            if (w < HALF_NEG) continue;
            add4(s, g4(Qg[*it], e, i, d) + w + pair_bonus);
          }
        }
      }
    }
    // Qgls closings (f)
    {
      size_t k = off[1][j];
      for (Grp& g : Qgls[j].grps) {
        const int e = g.v;
        for (Ent& s : g.ents) {
          double& ow = own[0][k++];
          const int i = s.p - 1;
          if (i < 0) continue;
          const std::vector<int>& pr = partners_right[i];
          for (auto it = std::upper_bound(pr.begin(), pr.end(), j); it != pr.end(); ++it) {
            const int q = *it, unpaired = q - j - 1;
            if (e == j && unpaired > cap) break;
            if (!E.end_ok(i, q)) continue;
            double w = E.w_ap + E.w_bp + E.w_cp * unpaired + pair_bonus + E.at(i, q);
            if (dz) w += E.run(j + 1, q - 1);
            const double t = g4(Qg[q], e, i, s.u);
            if (t < HALF_NEG) continue;
            ow = log_add(ow, t + w);
            if (v.ml && km[ix(i, q)]) add4(s, t + w);
          }
        }
      }
    }
    // Qgrs closings (g)(h) at j + 1
    if (j + 1 < n) {
      const int q = j + 1;
      size_t k = off[4][j];
      for (Grp& g : Qgrs[j].grps) {
        const int e = g.v;
        for (Ent& s : g.ents) {
          double& ow = own[1][k++];
          const int c = s.p, d = s.u;
          auto close = [&](int i, double w) {
            const double t = g4(Qg[q], e, i, d);
            if (t < HALF_NEG) return;
            ow = log_add(ow, t + w);
            if (v.ml && km[ix(i, q)]) add4(s, t + w);
          };
          for (int i = lead_i_lo(c); i < c; ++i) {
            if (!E.pairs(i, q) || !E.end_ok(i, q)) continue;
            double w = E.w_ap + E.w_bp + E.w_cp * (c - i - 1) + E.at(i, q) + pair_bonus;
            if (dz) w += E.run(i + 1, c - 1);
            close(i, w);
          }
          if (c - 1 >= 0) {
            const OMap& at = qmp_at(c - 1);
            for (size_t r = 0; r < at.slots(); ++r) {
              const int o = at.keys[r] - 1;
              if (at.keys[r] < 0 || o < 0 || o >= c || !E.pairs(o, q) || !E.end_ok(o, q)) continue;
              close(o, at.vals[r] + E.w_ap + E.w_bp + E.at(o, q) + pair_bonus);
            }
          }
        }
      }
    }
    // Qg[p-1] -> Qgrs[j] through the right branch [p, j]
    if (v.ml) {
      auto push = [&](int p, double qm) {
        const int f = p - 1;
        for (Grp& g : Qg[f].grps) {
          Grp* rg = Qgrs[j].get(g.v);
          if (!rg) continue;
          for (Ent& s : g.ents) {
            if (!E.end_ok(s.p, f)) continue;
            auto it = rg->idx.find(key2(s.p, s.u));
            if (it == rg->idx.end()) continue;
            add4(s, rg->ents[it->second].b + qm + (E.w_bp + E.at(s.p, f)));
          }
        }
      };
      if (v.mx) {
        const OMap& at = qmp_at(j);
        for (size_t r = 0; r < at.slots(); ++r)
          if (at.keys[r] >= 1) push(at.keys[r], at.vals[r]);
      } else {
        for (int p = 1; p <= j; ++p) {
          const double qm = st.PrMbar(p, j);
          if (qm > HALF_NEG) push(p, qm);
        }
      }
    }
    // Qp, R role: the partner L-stem's exact Qgl alpha
    for (Grp& g : Qgr[j].grps) {
      const int f = g.v - 1;
      if (f < 0) continue;
      for (Ent& s : g.ents) {
        const int d = s.p;
        if (!E.end_ok(d, j)) continue;
        Grp* yg = Qgl[f].get(s.u + 1);
        if (!yg) continue;
        for (const Ent& t : yg->ents) {
          if (t.u != d - 1 || t.p > d - 1 || !E.end_ok(t.p, f)) continue;
          add4(s, t.a + ctx(t.p, j) + (E.at(t.p, f) + E.at(d, j)));
        }
      }
    }
    // Qp, L role: the partner R-stem's (hybrid) Qgr alpha
    Tab4* GR = &Qgr[j];
    if (!(v.hx && v.zx)) {
      gl = Tab4();
      gr = Tab4();
      const Tab4* L = &Qgl[j];
      if (!v.hx) {
        for (const Grp& g : (*hb)[j].grps) {
          const int f = g.v;
          for (const Ent& s : g.ents) {
            const int d = s.u;
            if (!E.end_ok(d, f)) continue;
            const double w = E.w_stem + E.at(d, f) + s.a;
            for (int e = d; e < f; ++e) {
              const double qz = qz_log(d + 1, e);
              if (qz > HALF_NEG) gl.acc(f, s.p, e, w + qz);
            }
          }
        }
        L = &gl;
      }
      for (const Grp& g : L->grps) {
        const int f = g.v;
        for (const Ent& s : g.ents)
          for (int e = s.u + 1; e <= f; ++e) {
            const double z = v.zx ? qz_log(e, f - 1) : st.Qzbar(e, f - 1);
            if (z > HALF_NEG) gr.acc(e, s.p, s.u, s.a + z);
          }
      }
      GR = &gr;
    }
    for (const Grp& g : GR->grps) {
      const int f = g.v - 1;
      if (f < 0) continue;
      for (const Ent& s : g.ents) {
        const int d = s.p;
        if (!E.end_ok(d, j)) continue;
        Grp* yg = Qgl[f].get(s.u + 1);
        if (!yg) continue;
        for (Ent& t : yg->ents) {
          if (t.u != d - 1 || t.p > d - 1 || !E.end_ok(t.p, f)) continue;
          add4(t, s.a + ctx(t.p, j) + (E.at(t.p, f) + E.at(d, j)));
        }
      }
    }
    // Qgr = Qgl x Qz
    for (Grp& g : Qgl[j].grps) {
      const int f = g.v;
      for (Ent& s : g.ents)
        for (int e = s.u + 1; e <= f; ++e) {
          const double qz = qz_log(e, f - 1);
          if (qz < HALF_NEG) continue;
          add4(s, g4(Qgr[j], e, s.p, s.u) + qz);
        }
    }
    // Qgl = Qg x Qz x beta2
    for (Grp& g : Qg[j].grps) {
      const int f = g.v;
      for (Ent& s : g.ents) {
        const int d = s.u;
        if (!E.end_ok(d, f)) continue;
        const double w = E.w_stem + E.at(d, f);
        for (int e = d; e < f; ++e) {
          const double qz = qz_log(d + 1, e);
          if (qz < HALF_NEG) continue;
          add4(s, g4(Qgl[j], f, s.p, e) + qz + w);
        }
      }
    }
    // Qgls = Qm x Qg (single-pair states: their seeds)
    if (v.ml)
      for (Grp& g : Qg[j].grps) {
        Grp* lg = Qgls[j].get(g.v);
        if (!lg) continue;
        for (Ent& s : g.ents) {
          const int c = s.p;
          if (c - 1 < 0 || !E.end_ok(c, j)) continue;
          const double w = E.w_bp + E.at(c, j);
          const OMap& at = qmp_at(c - 1);
          for (size_t r = 0; r < at.slots(); ++r) {
            if (at.keys[r] < 0) continue;
            auto it = lg->idx.find(key2(at.keys[r], s.u));
            if (it != lg->idx.end()) add4(s, lg->ents[it->second].b + at.vals[r] + w);
          }
        }
      }
  }
}

// --pr-oracle PATH (--beam inf, --stem-est-principled on, after the outside pass).  Variant
// sets (comma list): step1 (every piece exact, kept = nested P > 0: ex_B ml on lc 30, ex_noml
// ml off lc 30, ex_A ml off lc 2), step2 (one piece exact, the others nested, floors 1e-7 and
// 0, B growth), step2A (the same, A growth), check (all nested on the run's kept pairs,
// pr_loop / pr_ml, against pr_score4; run with --pr-arms-k 0 --pr-cut 0 and large
// --pr-ml-span / --pr-ml-close).  For every state, score - alpha = Ent.b (Qg, Qgl, Qgr) or
// own (Qgls, Qgrs) against the exact beta; w = exp(alpha + beta - logZ).  Lines:
//   V id var hb_sec pass_sec stats_sec n_kept
//   R id var tab q N W kept1000 oracle1000 rho nboth sumd maxd wsumd nmiss wmiss nspur
//   D id var Qg2 nboth sumd maxd nmiss    (Qg states with two or more pairs)
//   C id tab N nmis maxd    (check)
void DP03::pr_oracle_diag(StemEst& st, FILE* f, const std::string& id, const std::string& sets) {
  outside();
  stem_ = &st;
  const char* sid = id.c_str();
  const char* TN[5] = {"Qg", "Qgls", "Qgl", "Qgr", "Qgrs"};
  std::vector<Tab4>* T5[5] = {&Qg, &Qgls, &Qgl, &Qgr, &Qgrs};
  std::vector<size_t> off[5];
  std::vector<double> bx[5];
  for (int t = 0; t < 5; ++t) {
    off[t].assign(n + 1, 0);
    for (int q = 0; q < n; ++q) {
      size_t c = 0;
      for (const Grp& g : (*T5[t])[q].grps) c += g.ents.size();
      off[t][q + 1] = off[t][q] + c;
    }
    bx[t].resize(off[t][n]);
    size_t k = 0;
    for (int q = 0; q < n; ++q)
      for (const Grp& g : (*T5[t])[q].grps)
        for (const Ent& e : g.ents) bx[t][k++] = e.b;
  }
  std::vector<double> own[2];
  own[0].assign(off[1][n], NEG_INF);
  own[1].assign(off[4][n], NEG_INF);
  const size_t nn = static_cast<size_t>(n) * n;
  std::vector<double> mbx(nn, NEG_INF);
  for (int j = 0; j < n; ++j) {
    const OMap& at = qmp_at(j);
    for (size_t r = 0; r < at.slots(); ++r)
      if (at.keys[r] >= 0) mbx[ix(at.keys[r], j)] = at.vals[r];
  }

  std::vector<PrVar> vars;
  bool check = false;
  auto has = [&](const char* s) {
    return ("," + sets + ",").find(std::string(",") + s + ",") != std::string::npos;
  };
  if (has("step1")) {
    vars.push_back({"ex_B", true, true, true, true, 1, true, 30});
    vars.push_back({"ex_noml", true, true, true, true, 1, false, 30});
    vars.push_back({"ex_A", true, true, true, true, 1, false, 2});
  }
  for (int a = 0; a < 2; ++a) {
    if (!has(a == 0 ? "step2" : "step2A")) continue;
    const bool ml = a == 0;
    const int lc = a == 0 ? 30 : 2;
    const std::string g = a == 0 ? "" : "A";
    for (int fi = 0; fi < 2; ++fi) {
      const std::string s = g + (fi == 0 ? "_f7" : "_f0");
      vars.push_back({"nest" + s, false, false, false, false, fi, ml, lc});
      vars.push_back({"Z" + s, false, true, false, false, fi, ml, lc});
      vars.push_back({"C" + s, false, false, false, true, fi, ml, lc});
      vars.push_back({"M" + s, false, false, true, false, fi, ml, lc});
      vars.push_back({"H" + s, true, false, false, false, fi, ml, lc});
    }
  }
  if (has("check")) {
    check = true;
    vars.push_back({"check", false, false, false, false, -1, pr_ml, pr_loop});
  }

  std::vector<char> km;
  std::vector<std::vector<int>> kr, kl;
  int km_fi = -2;
  auto set_mask = [&](int fi) {
    if (fi == km_fi) return;
    km_fi = fi;
    if (fi < 0) {
      km = st.prKept;
    } else {
      const double lf = fi == 0 ? std::log(1e-7) : NEG_INF;
      km.assign(nn, 0);
      for (size_t k = 0; k < nn; ++k)
        if (st.prLP[k] > HALF_NEG && !(st.prLP[k] < lf)) km[k] = 1;
    }
    kr.assign(n, {});
    kl.assign(n, {});
    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j)
        if (km[ix(i, j)]) kr[i].push_back(j), kl[j].push_back(i);
  };
  std::vector<Tab4> hb;
  std::tuple<int, bool, bool, int> hb_key{-9, false, false, 0};

  std::vector<double> sc, w, rs, rw;
  std::vector<char> keep;
  for (const PrVar& v : vars) {
    const auto t0 = std::chrono::steady_clock::now();
    set_mask(v.fi);
    long long nk = 0;
    for (char c : km) nk += c;
    if (!v.hx) {
      const std::tuple<int, bool, bool, int> key{v.fi, v.mx, v.ml, std::min(v.lc, cap)};
      if (key != hb_key) {
        hb_key = key;
        pr_build_hbar(km, kl, v.mx ? mbx : st.prMbar, v.ml, std::min(v.lc, cap), hb);
      }
    }
    const auto t1 = std::chrono::steady_clock::now();
    pr_hybrid_outside(v, st, km, kr, v.hx ? nullptr : &hb, own, off);
    const auto t2 = std::chrono::steady_clock::now();
    {
      long long nb = 0, nm = 0;
      double sd = 0, md = 0;
      size_t k = 0;
      for (int q = 0; q < n; ++q)
        for (const Grp& g : Qg[q].grps)
          for (const Ent& e : g.ents) {
            const double b0 = bx[0][k++];
            if (e.p == e.u || b0 < HALF_NEG) continue;
            if (e.b < HALF_NEG) {
              ++nm;
              continue;
            }
            const double d = std::fabs(e.b - b0);
            ++nb;
            sd += d;
            md = std::max(md, d);
          }
      std::fprintf(f, "D %s %s Qg2 %lld %.9e %.9e %lld\n", sid, v.name.c_str(), nb, sd, md, nm);
    }
    for (int t = 0; t < 5; ++t)
      for (int q = 0; q < n; ++q) {
        const size_t lo = off[t][q], N = off[t][q + 1] - lo;
        if (N == 0) continue;
        sc.resize(N);
        w.resize(N);
        double W = 0, sumd = 0, maxd = 0, wsumd = 0, wmiss = 0;
        long long nboth = 0, nmiss = 0, nspur = 0;
        size_t k = 0;
        for (const Grp& g : (*T5[t])[q].grps)
          for (const Ent& e : g.ents) {
            const double b0 = bx[t][lo + k];
            const double b1 = t == 1 ? own[0][lo + k] : t == 4 ? own[1][lo + k] : e.b;
            const double wk = (b0 > HALF_NEG && logZ > HALF_NEG) ? std::exp(e.a + b0 - logZ) : 0.0;
            w[k] = wk;
            W += wk;
            sc[k] = b1 > HALF_NEG ? e.a + b1 : NEG_INF;
            if (b0 > HALF_NEG) {
              if (b1 > HALF_NEG) {
                const double d = std::fabs(b1 - b0);
                ++nboth;
                sumd += d;
                wsumd += wk * d;
                maxd = std::max(maxd, d);
              } else {
                ++nmiss;
                wmiss += wk;
              }
            } else if (b1 > HALF_NEG) {
              ++nspur;
            }
            ++k;
          }
        double kept = W, orc = W, rho = std::nan("");
        if (N > 1000) {
          stemdiag::top_b(sc, 1000, keep);
          kept = 0;
          for (size_t m = 0; m < N; ++m)
            if (keep[m]) kept += w[m];
          stemdiag::top_b(w, 1000, keep);
          orc = 0;
          for (size_t m = 0; m < N; ++m)
            if (keep[m]) orc += w[m];
        }
        if (N >= 3) {
          stemdiag::avg_ranks(sc, rs);
          stemdiag::avg_ranks(w, rw);
          rho = stemdiag::pearson(rs, rw);
        }
        std::fprintf(f, "R %s %s %s %d %zu %.9e %.9e %.9e %.6f %lld %.9e %.9e %.9e %lld %.9e %lld\n",
                     sid, v.name.c_str(), TN[t], q, N, W, kept, orc, rho, nboth, sumd, maxd, wsumd,
                     nmiss, wmiss, nspur);
      }
    const auto t3 = std::chrono::steady_clock::now();
    auto sec = [](auto a, auto b) { return std::chrono::duration<double>(b - a).count(); };
    std::fprintf(f, "V %s %s %.2f %.2f %.2f %lld\n", sid, v.name.c_str(), sec(t0, t1), sec(t1, t2),
                 sec(t2, t3), nk);
    std::fflush(f);
    std::fprintf(stderr, "[pr-oracle] %s %s kept %lld hb %.1fs pass %.1fs stats %.1fs\n", sid,
                 v.name.c_str(), nk, sec(t0, t1), sec(t1, t2), sec(t2, t3));

    if (check && v.name == "check") {
      stA.clear();
      stE.clear();
      prKc.clear();
      prSc.clear();
      pr_sc_n = 0;
      stc_q = -1;
      const int BT[5] = {BT_QG, BT_QGLS, BT_QGL, BT_QGR, BT_QGRS};
      long long cn[5] = {}, cm[5] = {};
      double cd[5] = {};
      for (int q = 0; q < n; ++q) {
        stem_cache_step(q);
        for (int t = 0; t < 5; ++t) {
          size_t k = off[t][q];
          for (const Grp& g : (*T5[t])[q].grps)
            for (const Ent& e : g.ents) {
              const double b1 = t == 1 ? own[0][k] : t == 4 ? own[1][k] : e.b;
              ++k;
              const double s = pr_score4(BT[t], e, g.v, q);
              const double x = b1 > HALF_NEG ? e.a + b1 : NEG_INF;
              ++cn[t];
              if (s < HALF_NEG && x < HALF_NEG) continue;
              const double d = (s < HALF_NEG || x < HALF_NEG) ? INFINITY : std::fabs(s - x);
              if (d > 1e-9) ++cm[t];
              cd[t] = std::max(cd[t], d);
            }
        }
        stem_add_A(q);
      }
      for (int t = 0; t < 5; ++t) {
        std::fprintf(f, "C %s %s %lld %lld %.3e\n", sid, TN[t], cn[t], cm[t], cd[t]);
        std::fprintf(stderr, "[pr-oracle check] %s %s states %lld mismatch %lld max_absdiff %.3e\n",
                     sid, TN[t], cn[t], cm[t], cd[t]);
      }
    }
  }
  // score: pr_score4 as configured (--pr-abl etc.) on every state against its exact w, per
  // table and step (beams 100, 300, 1000; see the header)
  if (has("score")) {
    const auto t0 = std::chrono::steady_clock::now();
    stA.clear();
    stE.clear();
    prKc.clear();
    prSc.clear();
    pr_sc_n = 0;
    stc_q = -1;
    const int BT[5] = {BT_QG, BT_QGLS, BT_QGL, BT_QGR, BT_QGRS};
    for (int q = 0; q < n; ++q) {
      stem_cache_step(q);
      for (int t = 0; t < 5; ++t) {
        const size_t lo = off[t][q], N = off[t][q + 1] - lo;
        if (N == 0) continue;
        sc.resize(N);
        w.resize(N);
        double W = 0, wmiss = 0;
        long long nmiss = 0;
        size_t k = 0;
        for (const Grp& g : (*T5[t])[q].grps)
          for (const Ent& e : g.ents) {
            const double b0 = bx[t][lo + k];
            const double wk = (b0 > HALF_NEG && logZ > HALF_NEG) ? std::exp(e.a + b0 - logZ) : 0.0;
            const double s = pr_score4(BT[t], e, g.v, q);
            w[k] = wk;
            sc[k] = s;
            W += wk;
            if (s < HALF_NEG && wk > 0) ++nmiss, wmiss += wk;
            ++k;
          }
        double kb[3] = {W, W, W}, ob[3] = {W, W, W};
        const size_t BB[3] = {100, 300, 1000};
        for (int bi = 0; bi < 3; ++bi) {
          if (N <= BB[bi]) continue;
          stemdiag::top_b(sc, BB[bi], keep);
          kb[bi] = 0;
          for (size_t m = 0; m < N; ++m)
            if (keep[m]) kb[bi] += w[m];
          stemdiag::top_b(w, BB[bi], keep);
          ob[bi] = 0;
          for (size_t m = 0; m < N; ++m)
            if (keep[m]) ob[bi] += w[m];
        }
        double rho = std::nan("");
        if (N >= 3) {
          stemdiag::avg_ranks(sc, rs);
          stemdiag::avg_ranks(w, rw);
          rho = stemdiag::pearson(rs, rw);
        }
        std::fprintf(f, "A %s %s %d %zu %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.6f %lld %.9e\n", sid,
                     TN[t], q, N, W, kb[0], ob[0], kb[1], ob[1], kb[2], ob[2], rho, nmiss, wmiss);
      }
      stem_add_A(q);
    }
    std::fprintf(stderr, "[pr-oracle score] %s abl %d pk %d sep_evals %lld pp_mean %.3f ppk_mean %.3f %.1fs\n",
                 sid, pr_abl, pr_abl_pk, pr_sep_evals,
                 pr_sep_evals ? static_cast<double>(pr_sep_pp) / pr_sep_evals : 0.0,
                 pr_sep_evals ? static_cast<double>(pr_sep_ppk) / pr_sep_evals : 0.0,
                 std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
  }
  for (int t = 0; t < 5; ++t) {
    size_t k = 0;
    for (int q = 0; q < n; ++q)
      for (Grp& g : (*T5[t])[q].grps)
        for (Ent& e : g.ents) e.b = bx[t][k++];
  }
}

// --stem-diag (run at --beam inf after the outside pass): for every Qg / Qgl state, the exact
// posterior weight w = inside x outside / Z, its role (stem1: >= 99.9% of w through the left-stem
// Qp edge, stem2: <= 0.1%, mixed), and the pruning scores
//   s1 = log C[p-1] + log X + log Qz[u+1][v-1]                         (current)
//   s1H_L = s1 + log H_L(u, v, q) - log C0(q) + beta1 + 2 beta2      (stem 1, loops <= L)
//   s2 = log B(p, u+1, v) + log X + log S(q+1) + at(p, q) + beta1     (stem 2, Qgl only)
//   s_L = logsumexp(s1H_L, s2)                                        (Qgl only)
//   s2H_L = log S(q+1) + log X + log Qz[u+1][v-1] + log H2_L(p, u, v) - log S(p) + beta1
//           + 2 beta2  (Qg only: mirrored stem-2 estimate of --stem-est-role all)
//   sH_L = logsumexp(s1H_L, s2H_L)                                    (Qg only)
// with A, B from the exact Qgl table.  Per table, score, role class (all, stem1, stem2, mixed):
// per-step Spearman with log w (sum of rho x states, sum of states; live states only), set
// recall@b of the per-step top b by w (steps with > b states) and mass recall@b (all steps),
// b = 100, 300, 1000; per L the H = 0 states and the states the engine would remove.
std::string DP03::stem_diag(const std::vector<const StemEst*>& ests, const std::string& id) {
  using namespace stemdiag;
  auto t0 = std::chrono::steady_clock::now();
  outside();
  const int NL = static_cast<int>(ests.size());
  const size_t BEAMS[3] = {100, 300, 1000};
  std::vector<Tab4>* tabs[2] = {&Qg, &Qgl};
  auto flat_b = [&](int t, std::vector<double>& out) {
    out.clear();
    for (Tab4& tb : *tabs[t])
      for (Grp& g : tb.grps)
        for (Ent& e : g.ents) out.push_back(e.b);
  };
  std::vector<double> b0[2], b1[2], b2[2];
  for (int t = 0; t < 2; ++t) flat_b(t, b0[t]);
  tagged_outside(1);
  for (int t = 0; t < 2; ++t) flat_b(t, b1[t]);
  tagged_outside(2);
  for (int t = 0; t < 2; ++t) flat_b(t, b2[t]);
  for (int t = 0; t < 2; ++t) {
    size_t k = 0;
    for (Tab4& tb : *tabs[t])
      for (Grp& g : tb.grps)
        for (Ent& e : g.ents) e.b = b0[t][k++];
  }
  double tag_err = 0.0;
  long long tag_bad = 0;
  for (int t = 0; t < 2; ++t)
    for (size_t k = 0; k < b0[t].size(); ++k) {
      double s = log_add(b1[t][k], b2[t][k]);
      if (b0[t][k] > HALF_NEG && s > HALF_NEG) tag_err = std::max(tag_err, std::fabs(s - b0[t][k]));
      else if ((b0[t][k] > HALF_NEG) != (s > HALF_NEG)) ++tag_bad;
    }
  stA.clear();
  stE.clear();
  stc_q = -1;
  first_f = stemidx::MinTree(n + 1);
  for (int q = 0; q < n; ++q) stem_add_A(q);

  const double w_pk = E.w_beta1 + 2.0 * E.w_stem;
  const std::vector<double>& S = ests[0]->S;
  std::string out = "{\"id\":\"" + id + "\",\"n\":" + std::to_string(n) + ",\"logZ\":";
  num(out, logZ);
  out += ",\"loops\":[";
  for (int l = 0; l < NL; ++l) out += (l ? "," : "") + std::to_string(ests[l]->loop);
  out += "],\"stem\":[";
  for (int l = 0; l < NL; ++l)
    out += std::string(l ? "," : "") + "{\"entries\":" + std::to_string(ests[l]->entries) +
           ",\"dropped\":" + std::to_string(ests[l]->dropped) +
           ",\"nodes\":" + std::to_string(ests[l]->index.nodes()) +
           ",\"bytes\":" + std::to_string(ests[l]->bytes()) + "}";
  out += "],\"tag_err\":";
  num(out, tag_err);
  out += ",\"tag_bad\":" + std::to_string(tag_bad) + ",\"tables\":{";

  for (int t = 0; t < 2; ++t) {
    const bool gl = t == 1;
    std::vector<std::string> names = {"s1"};
    for (int l = 0; l < NL; ++l) names.push_back("s1H_L" + std::to_string(ests[l]->loop));
    if (gl) {
      names.push_back("s2");
      for (int l = 0; l < NL; ++l) names.push_back("s_L" + std::to_string(ests[l]->loop));
    } else {
      for (int l = 0; l < NL; ++l) names.push_back("s2H_L" + std::to_string(ests[l]->loop));
      for (int l = 0; l < NL; ++l) names.push_back("sH_L" + std::to_string(ests[l]->loop));
    }
    const int NS = static_cast<int>(names.size());
    // accumulators [score][class]
    std::vector<double> sp_sum(NS * 4, 0.0), sp_w(NS * 4, 0.0);
    std::vector<double> hit(NS * 4 * 3, 0.0), tot(NS * 4 * 3, 0.0), mk(NS * 4 * 3, 0.0),
        mt(NS * 4 * 3, 0.0);
    double cnt[4] = {}, mass[4] = {};
    long long states = 0, live_states = 0, dead_states = 0;
    // per L: H = 0 states, their live count / mass, with stem-1 outside, removable (+ live)
    std::vector<double> h0n(NL, 0), h0live(NL, 0), h0mass(NL, 0), h0b1(NL, 0), rmn(NL, 0),
        rmlive(NL, 0), rmmass(NL, 0);
    size_t off = 0;
    std::vector<double> w, sc_flat, rw, rs, sub_w, sub_s;
    std::vector<int> cls;
    std::vector<char> ks, ke;
    for (int q = 0; q < n; ++q) {
      Tab4& tb = (*tabs[t])[q];
      size_t N = 0;
      for (Grp& g : tb.grps) N += g.ents.size();
      if (N == 0) continue;
      w.assign(N, NEG_INF);
      cls.assign(N, -1);
      sc_flat.assign(static_cast<size_t>(NS) * N, NEG_INF);
      size_t k = 0;
      for (Grp& g : tb.grps)
        for (Ent& e : g.ents) {
          const size_t fk = off + k;
          const int p = e.p, u = e.u, v = g.v;
          if (b0[t][fk] > HALF_NEG) {
            w[k] = e.a + b0[t][fk] - logZ;
            double r1 = b1[t][fk] - log_add(b1[t][fk], b2[t][fk]);
            double share = b1[t][fk] > HALF_NEG ? std::exp(r1) : 0.0;
            cls[k] = share >= 0.999 ? 1 : share <= 0.001 ? 2 : 3;
          }
          const double s1 = four_score(p, u, v, e.a);
          sc_flat[k] = s1;
          double bval = NEG_INF, s2 = NEG_INF;
          if (gl) {
            if (E.end_ok(p, q)) bval = stem_B(p, u + 1, v);
            if (bval > HALF_NEG) s2 = bval + e.a + S[q + 1] + E.at(p, q) + E.w_beta1;
            sc_flat[static_cast<size_t>(1 + NL) * N + k] = s2;
          }
          for (int l = 0; l < NL; ++l) {
            double h = ests[l]->H(u, v, q);
            double s1h = h < HALF_NEG ? NEG_INF : s1 + h - ests[l]->C0[q] + w_pk;
            sc_flat[static_cast<size_t>(1 + l) * N + k] = s1h;
            if (gl) {
              sc_flat[static_cast<size_t>(2 + NL + l) * N + k] = log_add(s1h, s2);
            } else {
              double h2 = ests[l]->H2(p, u, v), s2h = NEG_INF;
              if (h2 > HALF_NEG) {
                s2h = S[q + 1] + e.a;
                if (score_qz) s2h += qz_log(u + 1, v - 1);
                s2h = s2h + h2 - S[p] + w_pk;
              }
              sc_flat[static_cast<size_t>(1 + NL + l) * N + k] = s2h;
              sc_flat[static_cast<size_t>(1 + 2 * NL + l) * N + k] = log_add(s1h, s2h);
            }
            const bool live = w[k] > HALF_NEG;
            if (h < HALF_NEG) {
              h0n[l] += 1;
              if (live) h0live[l] += 1, h0mass[l] += std::exp(w[k]);
              if (b1[t][fk] > HALF_NEG) h0b1[l] += 1;
            }
            bool rm = gl ? (!E.end_ok(p, q) || (h < HALF_NEG && bval < HALF_NEG))
                         : (h < HALF_NEG && first_f.query(u + 1, v) > v - 1);
            if (rm) {
              rmn[l] += 1;
              if (live) rmlive[l] += 1, rmmass[l] += std::exp(w[k]);
            }
          }
          ++k;
        }
      off += N;
      states += N;
      // live counts / mass per class (0 = all)
      for (size_t i = 0; i < N; ++i) {
        if (cls[i] < 0) {
          ++dead_states;
          continue;
        }
        ++live_states;
        double ew = std::exp(w[i]);
        cnt[0] += 1, mass[0] += ew;
        cnt[cls[i]] += 1, mass[cls[i]] += ew;
      }
      // Spearman per score / class
      for (int c = 0; c < 4; ++c) {
        std::vector<size_t> idx;
        for (size_t i = 0; i < N; ++i)
          if (cls[i] >= 0 && (c == 0 || cls[i] == c)) idx.push_back(i);
        if (idx.size() < 3) continue;
        sub_w.resize(idx.size());
        for (size_t m = 0; m < idx.size(); ++m) sub_w[m] = w[idx[m]];
        avg_ranks(sub_w, rw);
        for (int s = 0; s < NS; ++s) {
          sub_s.resize(idx.size());
          for (size_t m = 0; m < idx.size(); ++m) sub_s[m] = sc_flat[static_cast<size_t>(s) * N + idx[m]];
          avg_ranks(sub_s, rs);
          double rho = pearson(rw, rs);
          if (std::isnan(rho)) continue;
          sp_sum[s * 4 + c] += rho * idx.size();
          sp_w[s * 4 + c] += idx.size();
        }
      }
      // recall
      std::vector<size_t> live;
      std::vector<double> lw;
      for (size_t i = 0; i < N; ++i)
        if (cls[i] >= 0) live.push_back(i), lw.push_back(w[i]);
      for (int bi = 0; bi < 3; ++bi) {
        const size_t b = BEAMS[bi];
        std::vector<char> inE(N, 0);
        top_b(lw, b, ke);
        for (size_t m = 0; m < live.size(); ++m) inE[live[m]] = ke[m];
        for (int s = 0; s < NS; ++s) {
          std::vector<double> col(sc_flat.begin() + static_cast<size_t>(s) * N,
                                  sc_flat.begin() + static_cast<size_t>(s + 1) * N);
          top_b(col, b, ks);
          for (size_t i : live) {
            const double ew = std::exp(w[i]);
            for (int c : {0, cls[i]}) {
              const size_t a = (static_cast<size_t>(s) * 4 + c) * 3 + bi;
              mt[a] += ew;
              if (ks[i]) mk[a] += ew;
              if (N > b && inE[i]) {
                tot[a] += 1;
                if (ks[i]) hit[a] += 1;
              }
            }
          }
        }
      }
    }
    auto arr = [&](const std::vector<double>& v) {
      std::string o = "[";
      for (size_t i = 0; i < v.size(); ++i) {
        if (i) o += ",";
        num(o, v[i]);
      }
      return o + "]";
    };
    out += std::string(t ? "," : "") + "\"" + (gl ? "Qgl" : "Qg") + "\":{\"scores\":[";
    for (int s = 0; s < NS; ++s) out += std::string(s ? "," : "") + "\"" + names[s] + "\"";
    out += "],\"states\":" + std::to_string(states) + ",\"live\":" + std::to_string(live_states) +
           ",\"dead\":" + std::to_string(dead_states) + ",\"count\":" +
           arr(std::vector<double>(cnt, cnt + 4)) + ",\"mass\":" +
           arr(std::vector<double>(mass, mass + 4)) + ",\"sp_sum\":" + arr(sp_sum) +
           ",\"sp_w\":" + arr(sp_w) + ",\"hit\":" + arr(hit) + ",\"tot\":" + arr(tot) +
           ",\"mkept\":" + arr(mk) + ",\"mtot\":" + arr(mt) + ",\"h0_n\":" + arr(h0n) +
           ",\"h0_live\":" + arr(h0live) + ",\"h0_mass\":" + arr(h0mass) + ",\"h0_b1\":" +
           arr(h0b1) + ",\"rm_n\":" + arr(rmn) + ",\"rm_live\":" + arr(rmlive) +
           ",\"rm_mass\":" + arr(rmmass) + "}";
  }
  out += "},\"diag_sec\":";
  num(out, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
  out += "}";
  return out;
}

