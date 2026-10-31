// Pass-3 scan pruning: l2r_prune4, pruneS, pruneP, prune_sm, pruneA (exact and bounded),
// pruneB.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  void l2r_prune4(Tab4& t, int pt, Score score) {
    size_t tot = 0;
    for (const Grp& g : t.grps) tot += g.ents.size();
    const size_t pre = tot;
    l2r_pre[pt] += static_cast<long long>(tot);
    ++l2r_pos[pt];
    const int b = ltb(pt);
    if (b > 0 && pre > static_cast<size_t>(b)) ++l2r_bind[pt];
    if (b > 0 && tot > static_cast<size_t>(b)) {
      SecTimer tm(l2r_psec[pt]);
      std::vector<double> sc;
      sc.reserve(tot);
      for (const Grp& g : t.grps)
        for (const Ent& e : g.ents) sc.push_back(score(g, e));
      l2r_nsc[pt] += static_cast<long long>(tot);
      for (double s : sc)
        if (s < HALF_NEG) ++l2r_ninf[pt];
      std::vector<char> keep;
      keep_mask(sc, keep, b);
      filter_step(t, keep);
      tot = static_cast<size_t>(b);
    }
    l2r_kept[pt] += static_cast<long long>(tot);
    l2r_kmax[pt] = std::max(l2r_kmax[pt], static_cast<long long>(tot));
  }
  void l2r_pruneS(int j) {
    l2r_prune_sm(LS[j], PT_S, [&](uint64_t key, double f) { return l2r_scoreS(key, f, j); });
  }
  // P3 split counters, printed in [l2r-inc]; the split itself is not run (always 0).
  int lp_split_a = 0, lp_split_h = 0;
  long long lp_split_grp = 0, lp_split_grp_kept = 0;
  void l2r_pruneP(int pt, int j) {
    std::vector<SMap>& T = pt == PT_P1 ? LP1 : pt == PT_P2 ? LP2 : LP3;
    l2r_prune_sm(T[j], pt, [&](uint64_t key, double f) { return lp_score(pt, key, f, j); });
  }
  template <class Score>
  void l2r_prune_sm(SMap& m, int pt, Score score) {
    size_t tot = m.keys.size();
    l2r_pre[pt] += static_cast<long long>(tot);
    ++l2r_pos[pt];
    const int b = ltb(pt);
    const size_t pre = tot;
    if (b > 0 && pre > static_cast<size_t>(b)) ++l2r_bind[pt];
    if (b > 0 && tot > static_cast<size_t>(b)) {
      SecTimer tm(l2r_psec[pt]);
      std::vector<double> sc(tot);
      for (size_t k = 0; k < tot; ++k) {
        const uint64_t key = m.keys[k];
        sc[k] = score(key, m.v[k].F);
      }
      l2r_nsc[pt] += static_cast<long long>(tot);
      for (double s : sc)
        if (s < HALF_NEG) ++l2r_ninf[pt];
      std::vector<char> keep;
      keep_mask(sc, keep, b);
      SMap o;
      for (size_t k = 0; k < tot; ++k)
        if (keep[k]) {
          o.idx.emplace(m.keys[k], static_cast<int>(o.keys.size()));
          o.keys.push_back(m.keys[k]);
          o.v.push_back(m.v[k]);
        }
      m = std::move(o);
      tot = m.keys.size();
    }
    l2r_kept[pt] += static_cast<long long>(tot);
    l2r_kmax[pt] = std::max(l2r_kmax[pt], static_cast<long long>(tot));
  }
  // R-type tables at j: kind 0 root (seed: inner), 1 ls key, 2 rs key
  void l2r_pruneA(Tab4& t, int pt, int kind, bool seed, int j) {
    l2r_prune4(t, pt, [&](const Grp& g, const Ent& e) {
      if (kind != 0) ++sx_rkey_n;
      const RN& r = kind == 0 ? rec_root(e.p, j, seed) : rec_key(kind == 1, e.p, kind == 1 && seed, j);
      sx_ksrc = kind == 0 ? 1 : 2;
      const int arm = larm(e.u, g.v);
      if (arm < 0) ++l2r_noarm[pt];
      const double s = lsc_L(r, arm, kind == 0 ? 0 : 1);
      sx_ksrc = 0;
      return s > HALF_NEG ? e.a + s : NEG_INF;
    });
  }
  void l2r_pruneB(Tab4& t, int pt, int kind, bool seed, int j) {
    l2r_prune4(t, pt, [&](const Grp& g, const Ent& e) {
      double s;
      if (g_rfv) {
        s = kind == 0 ? rfv_root(e.p, j, seed, g.v) : rfv_key(kind == 1, e.p, kind == 1 && seed, j, g.v);
      } else {
        const RN& r = kind == 0 ? rec_root(e.p, j, seed) : rec_key(kind == 1, e.p, kind == 1 && seed, j);
        s = lsc_R(r, g.v);
      }
      return s > HALF_NEG ? e.a + s + prefix_log(e.u) : NEG_INF;
    });
  }


