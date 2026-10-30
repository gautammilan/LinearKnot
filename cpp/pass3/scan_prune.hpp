// Pass-3 scan pruning: l2r_prune4, pruneS, pruneP, prune_sm, pruneA (exact and bounded),
// pruneB; l2r_diag_report (--l2r-diag).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  void l2r_prune4(Tab4& t, int pt, Score score, int j) {
    std::unique_ptr<SecTimer> ptm(lpf ? new SecTimer(lpf_prune_sec[pt]) : nullptr);
    size_t tot = 0;
    for (const Grp& g : t.grps) tot += g.ents.size();
    const size_t pre = tot;
    l2r_pre[pt] += static_cast<long long>(tot);
    ++l2r_pos[pt];
    if (l2r_diag) {
      SecTimer tm(l2r_psec[pt]);
      std::vector<double> sc;
      sc.reserve(tot);
      for (const Grp& g : t.grps)
        for (const Ent& e : g.ents) sc.push_back(score(g, e));
      ldiag_put(pt, j, std::move(sc));
      l2r_kept[pt] += static_cast<long long>(tot);
      l2r_kmax[pt] = std::max(l2r_kmax[pt], static_cast<long long>(tot));
      return;
    }
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
      if (drops) {
        static const char* pn[PT_N] = {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us", "Uls", "Ulss",
                                       "Urs", "S", "P1", "P2", "P3"};
        std::vector<size_t> rk = score_ranks(sc);
        size_t k = 0;
        for (const Grp& g : t.grps)
          for (const Ent& e : g.ents) {
            if (!keep[k])
              note_l2r_drop(pn[pt], j, g.v, e.p, e.u, g.v, j, rk[k], tot, sc[k]);
            ++k;
          }
      }
      if (l2r_ninf_mass || l2r_keep_ninf) {
        size_t k = 0;
        for (const Grp& g : t.grps)
          for (const Ent& e : g.ents) {
            if (sc[k] < HALF_NEG) ninf_note(pt, j, ninf_key(g.v, e.u, e.p), e.a, keep[k]);
            ++k;
          }
      }
      if (g_beam_stats.f) {
        static const char* pn[PT_N] = {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us", "Uls", "Ulss",
                                       "Urs", "S", "P1", "P2", "P3"};
        g_beam_stats.log_prune_mask(pn[pt], j, sc, keep, b);
      }
      if (ltr_f) ltr_step4(t, pt, j, &sc, &keep, tot, b);
      filter_step(t, keep);
      tot = l2r_keep_ninf ? static_cast<size_t>(std::count(keep.begin(), keep.end(), 1)) : static_cast<size_t>(b);
    } else if (ltr_f) {
      ltr_step4(t, pt, j, nullptr, nullptr, tot, b);
    }
    if (lpf) lpf_prune_stat(pt, pre, tot, b, pre != tot || (b > 0 && pre > static_cast<size_t>(b)));
    l2r_kept[pt] += static_cast<long long>(tot);
    l2r_kmax[pt] = std::max(l2r_kmax[pt], static_cast<long long>(tot));
  }
  void l2r_pruneS(int j) {
    l2r_prune_sm(LS[j], PT_S, j, [&](uint64_t key, double f) { return l2r_scoreS(key, f, j); });
  }
  // --l2r-p3-split Ba,bh: P3 kept as the top Ba open a2 values by log sum of their states'
  // scores, then the top bh holes (xs, x') of each kept a2 by score (ties: earliest entry);
  // applied when P3 holds more than its beam (Ba bh).
  int lp_split_a = 0, lp_split_h = 0;
  long long lp_split_grp = 0, lp_split_grp_kept = 0;
  void keep_split(const std::vector<uint64_t>& keys, const std::vector<double>& sc, std::vector<char>& keep) {
    std::unordered_map<int, int> gi;
    std::vector<std::vector<int>> mem;
    std::vector<double> gs;
    for (size_t k = 0; k < keys.size(); ++k) {
      auto [it, ins] = gi.emplace(k3a(keys[k]), static_cast<int>(mem.size()));
      if (ins) {
        mem.emplace_back();
        gs.push_back(NEG_INF);
      }
      mem[it->second].push_back(static_cast<int>(k));
      if (sc[k] > HALF_NEG) gs[it->second] = log_add(gs[it->second], sc[k]);
    }
    keep.assign(keys.size(), 0);
    const int ng = static_cast<int>(mem.size());
    lp_split_grp += ng;
    std::vector<char> gk;
    if (ng > lp_split_a) keep_mask(gs, gk, lp_split_a);
    else gk.assign(ng, 1);
    std::vector<double> ms;
    std::vector<char> mk;
    for (int g = 0; g < ng; ++g) {
      if (!gk[g]) continue;
      ++lp_split_grp_kept;
      const std::vector<int>& v = mem[g];
      if (static_cast<int>(v.size()) <= lp_split_h) {
        for (int k : v) keep[k] = 1;
        continue;
      }
      ms.clear();
      for (int k : v) ms.push_back(sc[k]);
      keep_mask(ms, mk, lp_split_h);
      for (size_t r = 0; r < v.size(); ++r)
        if (mk[r]) keep[v[r]] = 1;
    }
  }
  void l2r_pruneP(int pt, int j) {
    std::vector<SMap>& T = pt == PT_P1 ? LP1 : pt == PT_P2 ? LP2 : LP3;
    l2r_prune_sm(T[j], pt, j, [&](uint64_t key, double f) { return lp_score(pt, key, f, j); });
  }
  template <class Score>
  void l2r_prune_sm(SMap& m, int pt, int j, Score score) {
    std::unique_ptr<SecTimer> ptm(lpf ? new SecTimer(lpf_prune_sec[pt]) : nullptr);
    size_t tot = m.keys.size();
    l2r_pre[pt] += static_cast<long long>(tot);
    ++l2r_pos[pt];
    const int b = l2r_diag ? 0 : ltb(pt);
    const size_t pre = tot;
    if (b > 0 && pre > static_cast<size_t>(b)) ++l2r_bind[pt];
    if (l2r_diag || (b > 0 && tot > static_cast<size_t>(b))) {
      SecTimer tm(l2r_psec[pt]);
      std::vector<double> sc(tot);
      for (size_t k = 0; k < tot; ++k) {
        const uint64_t key = m.keys[k];
        sc[k] = score(key, m.v[k].F);
      }
      if (l2r_diag) {
        ldiag_put(pt, j, std::move(sc));
        l2r_kept[pt] += static_cast<long long>(tot);
        l2r_kmax[pt] = std::max(l2r_kmax[pt], static_cast<long long>(tot));
        return;
      }
      l2r_nsc[pt] += static_cast<long long>(tot);
      for (double s : sc)
        if (s < HALF_NEG) ++l2r_ninf[pt];
      std::vector<char> keep;
      if (pt == PT_P3 && lp_split_a > 0) keep_split(m.keys, sc, keep);
      else keep_mask(sc, keep, b);
      if (drops) {
        static const char* pn[PT_N] = {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us", "Uls", "Ulss",
                                       "Urs", "S", "P1", "P2", "P3"};
        std::vector<size_t> rk = score_ranks(sc);
        for (size_t k = 0; k < tot; ++k) {
          if (keep[k]) continue;
          const uint64_t key = m.keys[k];
          note_l2r_drop(pn[pt], j, -1, k3a(key), k3b(key), k3c(key), j, rk[k], tot, sc[k]);
        }
      }
      if (l2r_ninf_mass || l2r_keep_ninf)
        for (size_t k = 0; k < tot; ++k)
          if (sc[k] < HALF_NEG) ninf_note(pt, j, m.keys[k], m.v[k].F, keep[k]);
      if (g_beam_stats.f) {
        static const char* pn[PT_N] = {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us", "Uls", "Ulss",
                                       "Urs", "S", "P1", "P2", "P3"};
        g_beam_stats.log_prune_mask(pn[pt], j, sc, keep, b);
      }
      if (ltr_f && pt == PT_S)
        for (size_t k = 0; k < tot; ++k)
          if (ltr_match3(m.keys[k]))
            ltr_line(pt, j, k3a(m.keys[k]), k3b(m.keys[k]), k3c(m.keys[k]), m.v[k].F, &sc, &keep, k,
                     tot, b);
      SMap o;
      for (size_t k = 0; k < tot; ++k)
        if (keep[k]) {
          o.idx.emplace(m.keys[k], static_cast<int>(o.keys.size()));
          o.keys.push_back(m.keys[k]);
          o.v.push_back(m.v[k]);
        }
      m = std::move(o);
      tot = m.keys.size();
    } else if (ltr_f && pt == PT_S) {
      for (size_t k = 0; k < tot; ++k)
        if (ltr_match3(m.keys[k]))
          ltr_line(pt, j, k3a(m.keys[k]), k3b(m.keys[k]), k3c(m.keys[k]), m.v[k].F, nullptr, nullptr,
                   k, tot, b);
    }
    if (lpf) lpf_prune_stat(pt, pre, tot, b, b > 0 && pre > static_cast<size_t>(b));
    l2r_kept[pt] += static_cast<long long>(tot);
    l2r_kmax[pt] = std::max(l2r_kmax[pt], static_cast<long long>(tot));
  }
  // R-type tables at j: kind 0 root (seed: inner), 1 ls key, 2 rs key
  struct ProfReg {
    DP03& d;
    int prev;
    bool on;
    ProfReg(DP03& d, int r) : d(d), prev(d.prof_region), on(d.prof_on) {
      if (on) d.prof_region = r;
    }
    ~ProfReg() {
      if (on) d.prof_region = prev;
    }
  };
  // Rank one R column by the exact gap sum. Every candidate gets sigma_ub from lK.
  // Exact sigma (lKx, memoized) is computed in that order and stops once the next
  // upper bound is below the current b-th exact score. exact <= ub carries the
  // per-state additives, so the comparison is on the full score.
  void l2r_pruneA_exact(Tab4& t, int pt, int kind, bool seed, int j) {
    ProfReg reg(*this, 1);
    std::unique_ptr<SecTimer> ptm(lpf ? new SecTimer(lpf_prune_sec[pt]) : nullptr);
    size_t tot = 0;
    for (const Grp& g : t.grps) tot += g.ents.size();
    const size_t pre = tot;
    l2r_pre[pt] += static_cast<long long>(tot);
    ++l2r_pos[pt];
    const int b = ltb(pt);
    if (b > 0 && pre > static_cast<size_t>(b)) ++l2r_bind[pt];
    if (!(b > 0 && tot > static_cast<size_t>(b))) {
      if (ltr_f) ltr_step4(t, pt, j, nullptr, nullptr, tot, b);
      if (lpf) lpf_prune_stat(pt, pre, tot, b, false);
      l2r_kept[pt] += static_cast<long long>(tot);
      l2r_kmax[pt] = std::max(l2r_kmax[pt], static_cast<long long>(tot));
      return;
    }
    SecTimer tm(l2r_psec[pt]);
    struct Ref {
      const Grp* g;
      const Ent* e;
    };
    std::vector<Ref> ref;
    ref.reserve(tot);
    for (const Grp& g : t.grps)
      for (const Ent& e : g.ents) ref.push_back(Ref{&g, &e});
    auto score = [&](const Ref& c, bool exact) {
      if (kind != 0 && !exact) ++sx_rkey_n;
      const RN& r =
          kind == 0 ? rec_root(c.e->p, j, seed) : rec_key(kind == 1, c.e->p, kind == 1 && seed, j);
      sx_ksrc = kind == 0 ? 1 : 2;
      sx_k_exact = exact;
      const int arm = larm(c.e->u, c.g->v);
      if (!exact && arm < 0) ++l2r_noarm[pt];
      const double s = lsc_L(r, arm, kind == 0 ? 0 : 1);
      sx_k_exact = false;
      sx_ksrc = 0;
      if (!exact && g_mass_f && s > HALF_NEG)
        g_mass_note(kind, seed ? 1 : 0, c.e->p, j, s);
      return s > HALF_NEG ? c.e->a + s + lkb_arm(c.e->p, arm) : NEG_INF;
    };
    std::vector<double> ub(tot), ex(tot, NEG_INF);
    std::vector<char> got(tot, 0);
    for (size_t i = 0; i < tot; ++i) ub[i] = score(ref[i], false);
    l2r_nsc[pt] += static_cast<long long>(tot);
    for (double s : ub)
      if (s < HALF_NEG) ++l2r_ninf[pt];
    std::vector<int> ord(tot);
    for (size_t i = 0; i < tot; ++i) ord[i] = static_cast<int>(i);
    std::sort(ord.begin(), ord.end(), [&](int a, int c) { return ub[a] > ub[c]; });
    std::priority_queue<double, std::vector<double>, std::greater<double>> heap;
    for (int idx : ord) {
      if (static_cast<int>(heap.size()) >= b && ub[idx] < heap.top()) break;
      ex[idx] = score(ref[idx], true);
      got[idx] = 1;
      ++sx_kx_sigma;
      if (static_cast<int>(heap.size()) < b) heap.push(ex[idx]);
      else if (ex[idx] > heap.top()) {
        heap.pop();
        heap.push(ex[idx]);
      }
    }
    std::vector<char> keep;
    keep_mask(ex, keep, b);
    if (k_exact_check) {
      ++k_exact_cols;
      std::vector<double> full(tot);
      for (size_t i = 0; i < tot; ++i) full[i] = got[i] ? ex[i] : score(ref[i], true);
      std::vector<char> keep2;
      keep_mask(full, keep2, b);
      if (keep != keep2) ++k_exact_mis;
    }
    if (l2r_ninf_mass || l2r_keep_ninf) {
      for (size_t k = 0; k < tot; ++k)
        if (ex[k] < HALF_NEG)
          ninf_note(pt, j, ninf_key(ref[k].g->v, ref[k].e->u, ref[k].e->p), ref[k].e->a, keep[k]);
    }
    if (ltr_f) ltr_step4(t, pt, j, &ex, &keep, tot, b);
    filter_step(t, keep);
    tot = l2r_keep_ninf ? static_cast<size_t>(std::count(keep.begin(), keep.end(), 1))
                        : static_cast<size_t>(b);
    if (lpf) lpf_prune_stat(pt, pre, tot, b, true);
    l2r_kept[pt] += static_cast<long long>(tot);
    l2r_kmax[pt] = std::max(l2r_kmax[pt], static_cast<long long>(tot));
  }
  void l2r_pruneA(Tab4& t, int pt, int kind, bool seed, int j) {
    if (k_exact_ub && !l2r_diag && !l2r_arm_fb && !(lsub_e && kind == 0 && seed)) {
      l2r_pruneA_exact(t, pt, kind, seed, j);
      return;
    }
    ProfReg reg(*this, 1);
    l2r_prune4(t, pt, [&](const Grp& g, const Ent& e) {
      if (lsub_e && kind == 0 && seed) {
        const double b = lors(e.u, g.v, e.p, j);
        return b > HALF_NEG ? e.a + b : NEG_INF;
      }
      if (kind != 0) ++sx_rkey_n;
      const RN& r = kind == 0 ? rec_root(e.p, j, seed) : rec_key(kind == 1, e.p, kind == 1 && seed, j);
      sx_ksrc = kind == 0 ? 1 : 2;
      if (l2r_arm_fb) {
        const auto [arm, fx] = larm_fb(e.u, g.v);
        if (arm >= 0 && arm != larm(e.u, g.v)) ++l2r_fb_hit;
        const double s = lsc_L(r, arm, kind == 0 ? 0 : 1);
        sx_ksrc = 0;
        if (s < HALF_NEG) {
          ++l2r_fb_floor;
          return e.a + kFbFloor;
        }
        return e.a + (s + fx) + lkb_arm(e.p, arm);
      }
      const int arm = larm(e.u, g.v);
      if (arm < 0) ++l2r_noarm[pt];
      const double s = lsc_L(r, arm, kind == 0 ? 0 : 1);
      sx_ksrc = 0;
      if (g_mass_f && s > HALF_NEG) g_mass_note(kind, seed ? 1 : 0, e.p, j, s);
      return s > HALF_NEG ? e.a + s + lkb_arm(e.p, arm) : NEG_INF;
    }, j);
  }
  void l2r_pruneB(Tab4& t, int pt, int kind, bool seed, int j) {
    l2r_prune4(t, pt, [&](const Grp& g, const Ent& e) {
      if (lsub_d && kind == 0) {
        const auto& m = seed ? lOUs : lOU;
        auto it = m.find(key4(e.u, g.v, e.p, j));
        return it == m.end() ? NEG_INF : e.a + it->second - (lsub_c ? 0.0 : lrho(e.u));
      }
      double s;
      if (g_rfv) {
        s = kind == 0 ? rfv_root(e.p, j, seed, g.v) : rfv_key(kind == 1, e.p, kind == 1 && seed, j, g.v);
      } else {
        const RN& r = kind == 0 ? rec_root(e.p, j, seed) : rec_key(kind == 1, e.p, kind == 1 && seed, j);
        s = lsc_R(r, g.v);
      }
      if (!lLam.empty())
        return s > HALF_NEG ? e.a + s + prefix_log(e.u) + (lsub_c ? lrho(e.u) : lkb_xs(e.u, g.v)) +
                                  lLam[e.u]
                            : NEG_INF;
      return s > HALF_NEG ? e.a + s + prefix_log(e.u) + (lsub_c ? lrho(e.u) : lkb_xs(e.u, g.v))
                          : NEG_INF;
    }, j);
  }

  // lost[t][0..1]: posterior mass outside the top-b by score / by posterior, b = bs[0], bs[1]
  std::string l2r_diag_report(const std::vector<int>& bs) {
    static const char* pn[PT_N] = {"R", "Rs", "Rls", "Rlss", "Rrs", "U", "Us", "Uls", "Ulss",
                                   "Urs", "S", "P1", "P2", "P3"};
    std::string out;
    char buf[512];
    for (int pt = 0; pt < PT_N; ++pt) {
      if (ldS[pt].empty()) continue;
      std::vector<double> lostS(bs.size(), 0.0), lostO(bs.size(), 0.0);
      double tot = 0.0, minf = 0.0, dsum = 0.0, dsq = 0.0, dw = 0.0, hdrop[8] = {};
      long long steps_over = 0, ninf = 0, hin[8] = {};
      std::vector<std::pair<double, double>> dq;  // (d, posterior weight) for the quantiles
      std::vector<std::array<double, 5>> ds;      // (d, weight, exterior / ml / gap share of b)
      for (int j = 0; j < n; ++j) {
        const std::vector<double>& sc = ldS[pt][j];
        std::vector<double> po, pbv;
        if (pt >= PT_S) {
          const SMap& m = l2r_smap(pt, j);
          for (size_t k = 0; k < m.keys.size(); ++k) {
            po.push_back(m.v[k].F + m.v[k].bF);
            pbv.push_back(m.v[k].bF);
          }
        } else {
          for (const Grp& g : l2r_tab(pt, j).grps)
            for (const Ent& e : g.ents) {
              po.push_back(e.a + e.b);
              pbv.push_back(e.b);
            }
        }
        if (po.size() != sc.size()) {
          std::snprintf(buf, sizeof buf, " [%s size mismatch at %d]", pn[pt], j);
          out += buf;
          continue;
        }
        std::vector<double> pe(po.size());
        for (size_t k = 0; k < po.size(); ++k) {
          pe[k] = po[k] > HALF_NEG ? std::exp(std::min(0.0, po[k] - logZ)) : 0.0;
          tot += pe[k];
          if (sc[k] < HALF_NEG) {
            ++ninf;
            minf += pe[k];
          } else if (po[k] > HALF_NEG) {
            // score - log posterior weight (the score's outside estimate / the true one)
            const double d = sc[k] - po[k];
            dsum += pe[k] * d;
            dsq += pe[k] * d * d;
            dw += pe[k];
            if (pe[k] > 0.0) {
              dq.emplace_back(d, pe[k]);
              std::array<double, 5> t{d, pe[k], 0.0, 0.0, 0.0};
              const double b = pbv[k];
              for (int c = 0; c < 3 && !lspB[c][pt].empty(); ++c) {
                const double bc = lspB[c][pt][j][k];
                t[2 + c] = bc > HALF_NEG ? std::exp(bc - b) : 0.0;
              }
              ds.push_back(t);
            }
          }
        }
        for (size_t c = 0; c < bs.size(); ++c) {
          const size_t b = static_cast<size_t>(bs[c]);
          if (po.size() <= b) continue;
          if (c == 0) ++steps_over;
          std::vector<char> keep, ko;
          keep_mask(sc, keep, static_cast<int>(b));
          for (size_t k = 0; k < po.size(); ++k)
            if (!keep[k]) lostS[c] += pe[k];
          keep_mask(po, ko, static_cast<int>(b));
          for (size_t k = 0; k < po.size(); ++k)
            if (!ko[k]) lostO[c] += pe[k];
          if (bs[c] != 1000) continue;
          // dropped by score, kept by the oracle: mass by d = score - log post; intruders
          // (kept by score, not by the oracle): count by d
          for (size_t k = 0; k < po.size(); ++k) {
            const double d = sc[k] < HALF_NEG ? -1e9 : po[k] < HALF_NEG ? 1e9 : sc[k] - po[k];
            const int bin = d < -4 ? 0 : d < -2 ? 1 : d < -1 ? 2 : d < 0 ? 3 : d < 1 ? 4 : d < 2 ? 5 : d < 4 ? 6 : 7;
            if (!keep[k] && ko[k]) hdrop[bin] += pe[k];
            if (keep[k] && !ko[k]) ++hin[bin];
          }
        }
      }
      const double mu = dw > 0 ? dsum / dw : 0.0;
      std::snprintf(buf, sizeof buf,
                    " %s mass %.4g steps_over %lld ninf %lld inf_mass %.4g logratio_mean %.3f sd %.3f",
                    pn[pt], tot, steps_over, ninf, minf, mu,
                    dw > 0 ? std::sqrt(std::max(0.0, dsq / dw - mu * mu)) : 0.0);
      out += buf;
      // low mode: peak of the weighted d histogram (bins 0.05) below -3; mean, sd, share of the
      // mass within 0.3 of it; consumer split of b in the mode and at |d| < 0.3
      {
        std::vector<double> h(340, 0.0);
        for (const auto& t : ds)
          if (t[0] >= -20.0 && t[0] < -3.0) h[static_cast<size_t>((t[0] + 20.0) / 0.05)] += t[1];
        const size_t pk = static_cast<size_t>(std::max_element(h.begin(), h.end()) - h.begin());
        const double c0 = -20.0 + 0.05 * (static_cast<double>(pk) + 0.5);
        double w = 0.0, m1 = 0.0, m2 = 0.0, sp[3] = {}, w0 = 0.0, s0[3] = {};
        for (const auto& t : ds) {
          if (std::fabs(t[0] - c0) < 0.3) {
            w += t[1];
            m1 += t[1] * t[0];
            m2 += t[1] * t[0] * t[0];
            for (int c = 0; c < 3; ++c) sp[c] += t[1] * t[2 + c];
          }
          if (std::fabs(t[0]) < 0.3) {
            w0 += t[1];
            for (int c = 0; c < 3; ++c) s0[c] += t[1] * t[2 + c];
          }
        }
        const double mu1 = w > 0 ? m1 / w : 0.0;
        std::snprintf(buf, sizeof buf, " mode %.3f sd %.3f share %.3f", mu1,
                      w > 0 ? std::sqrt(std::max(0.0, m2 / w - mu1 * mu1)) : 0.0, dw > 0 ? w / dw : 0.0);
        out += buf;
        if (!lspB[0][pt].empty()) {
          auto f = [](double a, double b) { return b > 0 ? a / b : 0.0; };
          std::snprintf(buf, sizeof buf, " split_mode(ext/ml/gap) %.4f/%.4f/%.4f split_d0 %.4f/%.4f/%.4f",
                        f(sp[0], w), f(sp[1], w), f(sp[2], w), f(s0[0], w0), f(s0[1], w0), f(s0[2], w0));
          out += buf;
        }
      }
      // posterior-weighted quantiles of d = log(score / (inside x exact outside))
      std::sort(dq.begin(), dq.end());
      out += " q01/05/25/50/75/95/99";
      {
        const double qs[7] = {0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99};
        double c = 0.0;
        size_t k = 0;
        for (double q : qs) {
          while (k < dq.size() && c + dq[k].second < q * dw) c += dq[k++].second;
          std::snprintf(buf, sizeof buf, " %.2f", k < dq.size() ? dq[k].first : 0.0);
          out += buf;
        }
      }
      for (size_t c = 0; c < bs.size(); ++c) {
        std::snprintf(buf, sizeof buf, " b%d score %.4g oracle %.4g", bs[c], lostS[c], lostO[c]);
        out += buf;
      }
      out += " b1000 drop_mass_by_d(<-4,-4,-2,-1,0,1,2,>=4)";
      for (int k = 0; k < 8; ++k) {
        std::snprintf(buf, sizeof buf, " %.2g", hdrop[k]);
        out += buf;
      }
      out += " intruders_by_d";
      for (int k = 0; k < 8; ++k) {
        std::snprintf(buf, sizeof buf, " %lld", hin[k]);
        out += buf;
      }
      out += " |";
    }
    return out;
  }

