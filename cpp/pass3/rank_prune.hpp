// Pass-3 pruning by ranking score: stem profiles (prof_*), prune4_helper / prune4_stem /
// prune4_step, score ranks, prune1 / prune2 (top-b of a column).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.


  void prof_finish_key(Prof& p, int kind) {
    pf_ly += static_cast<long long>(p.ly.size());
    pf_rx += static_cast<long long>(p.rx.size());
    ++pf_keys[kind];
  }

  // Key profile of the Qg states with outer pair p.stc_q: the pair profile, truncated.
  const Prof& qg_prof(int p) {
    pf_note_step();
    auto it = stPg.find(p);
    if (it != stPg.end()) return it->second;
    const auto t0 = std::chrono::steady_clock::now();
    Prof k = pair_prof(p, stc_q);
    prof_finish_key(k, 0);
    pf_build_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return stPg.emplace(p, std::move(k)).first->second;
  }

  // Sum over the closings (closing order) of w + the pair profile of the closed pair.
  void prof_closings(const std::vector<Closing>& cl, Prof& out) {
    ProfAcc acc;
    for (const Closing& c : cl) {
      const Prof& pp = pair_prof(c.x, c.y);
      for (size_t k = 0; k < pp.ly.size(); ++k) acc.addL(pp.ly[k], c.w + pp.lw[k]);
      for (size_t k = 0; k < pp.rx.size(); ++k) acc.addR(pp.rx[k], c.w + pp.rw[k]);
    }
    acc.finish(out);
    out.ncl = static_cast<int>(cl.size());
  }


  // Key profile of the Qgls[q] states (i, ., v) with seed = (v == q).
  const Prof& ls_prof(int i, bool seed, int q) {
    pf_note_step();
    const int key = 2 * i + (seed ? 1 : 0);
    auto it = stPls.find(key);
    if (it != stPls.end()) return it->second;
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<Closing> cl;
    closings_ls(i, seed ? q : -1, q, cl);
    Prof k;
    prof_closings(cl, k);
    prof_finish_key(k, 1);
    pf_build_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return stPls.emplace(key, std::move(k)).first->second;
  }

  // Key profile of the Qgrs[q] states (c, ., .).
  const Prof& rs_prof(int c, int q) {
    pf_note_step();
    auto it = stPrs.find(c);
    if (it != stPrs.end()) return it->second;
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<Closing> cl;
    closings_rs(c, q, cl);
    Prof k;
    prof_closings(cl, k);
    prof_finish_key(k, 2);
    pf_build_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return stPrs.emplace(c, std::move(k)).first->second;
  }

  // Gap part of a profile score, alpha excluded (DP03._prof_eval): logsumexp of
  // Qz[u+1][v-1] sum_y' H_R(u, v, y') P_L(y') (consistent: times at(u, v)),
  // sum_x' Bg(x', u, v) P_R(x') and Qz[u+1][v-1] tail.
  double prof_eval(const Prof& p, int u, int v) {
    if (u != row_u || v != row_v || stc_q != row_q) {
      if (rowH.size() != static_cast<size_t>(n + 1)) {
        rowH.assign(n + 1, NEG_INF);
        rowB.assign(n + 1, NEG_INF);
        rowHs.assign(n + 1, 0);
        rowBs.assign(n + 1, 0);
      }
      row_u = u;
      row_v = v;
      row_q = stc_q;
      ++row_stamp;
    }
    double L = NEG_INF;
    for (size_t k = 0; k < p.ly.size(); ++k) {
      const int y = p.ly[k];
      if (rowHs[y] != row_stamp) {
        rowH[y] = stem_H_cached(u, v, y);
        rowHs[y] = row_stamp;
      }
      const double h = rowH[y];
      if (h > HALF_NEG) L = log_add(L, h + p.lw[k]);
    }
    const double qz = score_qz ? qz_log(u + 1, v - 1) : 0.0;
    if (L > HALF_NEG) {
      L += qz;
    }
    double R = NEG_INF, Ro = NEG_INF;
    if (E.end_ok(u, v))
      for (size_t k = 0; k < p.rx.size(); ++k) {
        const int x = p.rx[k];
        if (rowBs[x] != row_stamp) {
          rowB[x] = stem_Bg(x, u, v);
          rowBs[x] = row_stamp;
        } else {
          ++bg_lookups;
        }
        const double bg = rowB[x];
        if (bg > HALF_NEG) Ro = log_add(Ro, bg + p.rw[k]);
      }
    R = Ro;
    double s = log_add(L, R);
    if (p.tail > HALF_NEG) s = log_add(s, p.tail + qz);
    pf_evals += static_cast<long long>(p.ly.size() + p.rx.size());
    return s;
  }





  // Scoring order of m profile states: by gap (v, u), stable, so prof_eval's rows are reused.
  // Scores do not depend on the order.
  template <class F>
  static std::vector<size_t> gap_order(size_t m, F key) {
    std::vector<std::pair<std::pair<int, int>, size_t>> ks(m);
    for (size_t r = 0; r < m; ++r) ks[r] = {key(r), r};
    std::sort(ks.begin(), ks.end());
    std::vector<size_t> out(m);
    for (size_t r = 0; r < m; ++r) out[r] = ks[r].second;
    return out;
  }


  // --stem-est-helpers grow / nested: rank Qgls / Qgrs by helper_sigma / helper_sigma_nested
  // (step: whole table, group: per group v; threshold rule of keep_mask).  No state is removed.
  void prune4_helper(Tab4& t, int bt, int q, int beam) {
    if (beam <= 0) return;
    const auto t0 = std::chrono::steady_clock::now();
    stem_cache_step(q);
    std::vector<const Ent*> ent;
    std::vector<int> gv;
    std::vector<std::pair<size_t, size_t>> spans;
    for (const Grp& g : t.grps) {
      const size_t lo = ent.size();
      for (const Ent& e : g.ents) {
        ent.push_back(&e);
        gv.push_back(g.v);
      }
      if (!gap_step) spans.emplace_back(lo, ent.size());
    }
    if (gap_step) spans.emplace_back(0, ent.size());
    const bool ls = bt == BT_QGLS;
    std::vector<char> keep(ent.size(), 1), km;
    std::vector<double> sc, so;
    std::vector<Closing> cl;
    bool any = false;
    for (const auto& sp : spans) {
      const size_t m = sp.second - sp.first;
      if (m <= static_cast<size_t>(beam)) continue;
      sc.clear();
      if (pr_on) {
        const auto tp = std::chrono::steady_clock::now();
        const int ti = pr_tab(bt);
        sc.assign(m, NEG_INF);
        for (size_t r : gap_order(m, [&](size_t r) {
               return std::make_pair(gv[sp.first + r], ent[sp.first + r]->u);
             })) {
          const size_t k = sp.first + r;
          sc[r] = pr_score4(bt, *ent[k], gv[k], q);
          hlp_finite += sc[r] > HALF_NEG;
        }
        pr_states[ti] += static_cast<long long>(m);
        pr_sec[ti] += std::chrono::duration<double>(std::chrono::steady_clock::now() - tp).count();
      } else if (stem_profile) {
        sc.assign(m, NEG_INF);
        for (size_t r : gap_order(m, [&](size_t r) {
               return std::make_pair(gv[sp.first + r], ent[sp.first + r]->u);
             })) {
          const size_t k = sp.first + r;
          const Ent& e = *ent[k];
          const Prof& kp = ls ? ls_prof(e.p, gv[k] == q, q) : rs_prof(e.p, q);
          const double s = e.a + prof_eval(kp, e.u, gv[k]);
          sc[r] = s;
          hlp_closings += kp.ncl;
          hlp_finite += s > HALF_NEG;
          ++pf_states;
        }
      }
      for (size_t k = sp.first; k < sp.second && !stem_profile && !pr_on; ++k) {
        const Ent& e = *ent[k];
        if (stem_helpers == 2) {
          int ncl = 0;
          const double s = helper_sigma_nested(ls, e, gv[k], q, ncl);
          sc.push_back(s);
          hlp_closings += ncl;
          hlp_finite += s > HALF_NEG;
          continue;
        }
        if (ls) closings_ls(e.p, gv[k], q, cl);
        else closings_rs(e.p, q, cl);
        const double s = helper_sigma(cl, e.u, gv[k], e.a);
        sc.push_back(s);
        hlp_closings += static_cast<long long>(cl.size());
        hlp_finite += s > HALF_NEG;
      }
      hlp_states[ls ? 0 : 1] += static_cast<long long>(m);
      keep_mask(sc, km, beam);
      for (size_t r = 0; r < m; ++r) {
        keep[sp.first + r] = km[r];
        any |= !km[r];
      }
    }
    if (any) filter_step(t, keep);
    const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    hlp_sec += dt;
    if (stem_profile) pf_total_sec += dt;
  }



  double stem_A(int d, int e, int f) const {
    auto it = stA.find(key2(d, e));
    if (it == stA.end()) return NEG_INF;
    const StemA& a = it->second;
    auto k = std::upper_bound(a.fs.begin(), a.fs.end(), f) - a.fs.begin() - 1;
    return k >= 0 && a.fs[k] == f ? a.vals[k] : NEG_INF;
  }

  // A(u+1, v, q) += C[p-1] Qgl at(p, q) for every surviving Qgl state at step q (with ok:
  // only the flattened states k with ok[k]).
  void stem_add_A(int q, const std::vector<char>* ok = nullptr) {
    size_t k = 0;
    for (const Grp& g : Qgl[q].grps)
      for (const Ent& s : g.ents) {
        if (ok && !(*ok)[k++]) continue;
        if (!E.end_ok(s.p, q)) continue;
        double val = prefix_log(s.p) + s.a + E.at(s.p, q);
        uint64_t key = key2(s.u + 1, g.v);
        auto it = stA.find(key);
        if (it == stA.end()) {
          it = stA.emplace(key, StemA{{q}, {val}, {}}).first;
          if (pr_on) it->second.il.emplace_back();
          first_f.set_min(g.v, q);
          std::vector<int>& es = stE[s.u + 1];
          es.insert(std::lower_bound(es.begin(), es.end(), g.v), g.v);
        } else if (it->second.fs.back() == q) {
          it->second.vals.back() = log_add(it->second.vals.back(), val);
        } else {
          it->second.fs.push_back(q);
          it->second.vals.push_back(val);
          if (pr_on) it->second.il.emplace_back();
        }
        if (pr_on) it->second.il.back().emplace_back(s.p, s.a + E.at(s.p, q));
        if (pr_on && pr_gsel)
          prPend.push_back({s.u + 1, g.v, q, s.p, s.a + E.at(s.p, q)});
      }
  }

  // Feasibility and stem score of the flattened states of Qg / Qgl / Qgr at step q
  // (DP03._stem_eval in engine.py).
  struct StemEval {
    DP03& dp;
    int bt, q;
    std::vector<const Grp*> grp;
    std::vector<const Ent*> ent;
    std::vector<char> ok;
    std::vector<double> bval, aval;
    std::unordered_map<uint64_t, double> hcache;
    double w_pk, c0q, sq, cq;

    StemEval(DP03& d, const Tab4& t, int bt_, int q_) : dp(d), bt(bt_), q(q_) {
      const Energy& E = dp.E;
      dp.stem_cache_step(q);
      w_pk = E.w_beta1 + 2.0 * E.w_stem;
      c0q = dp.stem_->C0[q];
      sq = dp.stem_->S[q + 1];
      cq = dp.hr_c(q);
      for (const Grp& g : t.grps)
        for (const Ent& e : g.ents) {
          grp.push_back(&g);
          ent.push_back(&e);
        }
      const size_t N = ent.size();
      ok.assign(N, 0);
      bval.assign(N, std::nan(""));
      if (bt == BT_QG) {
        for (size_t k = 0; k < N; ++k) {
          int u = ent[k]->u, v = grp[k]->v;
          ok[k] = H(u, v) > HALF_NEG || dp.first_f.query(u + 1, v) <= v - 1;
        }
      } else if (bt == BT_QGL) {
        for (size_t k = 0; k < N; ++k) {
          if (!E.end_ok(ent[k]->p, q)) continue;
          ok[k] = H(ent[k]->u, grp[k]->v) > HALF_NEG || B(k) > HALF_NEG;
        }
      } else {
        aval.assign(N, NEG_INF);
        for (size_t k = 0; k < N; ++k) {
          int p = ent[k]->p, u = ent[k]->u, v = grp[k]->v;
          if (E.end_ok(p, q)) aval[k] = dp.stem_A(p, u + 1, v - 1);
          ok[k] = aval[k] > HALF_NEG;
        }
      }
    }

    double H(int u, int v) {
      uint64_t key = key2(u, v);
      auto it = hcache.find(key);
      if (it != hcache.end()) return it->second;
      double h = dp.stem_->H(u, v, q);
      hcache.emplace(key, h);
      return h;
    }

    double B(size_t k) {
      if (std::isnan(bval[k])) {
        const int d = ent[k]->p, e = ent[k]->u + 1, v = grp[k]->v;
        bval[k] = dp.stem_qg_partner ? dp.stem_B_cached(d, e, v) : dp.stem_B(d, e, v);
      }
      return bval[k];
    }

    double s1h(size_t k) {
      const Ent& e = *ent[k];
      double h = H(e.u, grp[k]->v);
      if (h < HALF_NEG) return NEG_INF;
      return dp.four_score(e.p, e.u, grp[k]->v, e.a) + h + cq + w_pk;
    }

    // Stem-1 score (state as the left stem, H_R); NEG_INF if none.
    double sigL(size_t k) {
      if (bt == BT_QGR) return NEG_INF;
      return s1h(k);
    }

    // Stem-2 score (state as the right stem); NEG_INF if none or the role is off.
    double sigR(size_t k) {
      const Energy& E = dp.E;
      const Ent& e = *ent[k];
      const int v = grp[k]->v;
      if (bt == BT_QG) {
        if (!dp.stem_qg2) return NEG_INF;
        if (dp.stem_qg_partner) {
          double bg = dp.stem_Bg(e.p, e.u, v);
          if (bg < HALF_NEG) return NEG_INF;
          double s2 = bg + e.a + sq;
          return s2;
        }
        double h2 = dp.stem_->H2(e.p, e.u, v);
        if (h2 < HALF_NEG) return NEG_INF;
        double s2 = sq + e.a;
        if (dp.score_qz) s2 += dp.qz_log(e.u + 1, v - 1);
        s2 = s2 + h2 - dp.stem_->S[e.p] + w_pk;
        return s2;
      }
      if (!dp.stem_both) return NEG_INF;
      if (bt == BT_QGL) {
        double b = B(k);
        return b < HALF_NEG ? NEG_INF : b + e.a + sq + E.at(e.p, q) + E.w_beta1;
      }
      return aval[k] + e.a + sq + E.at(e.p, q) + E.w_beta1;
    }

    double score(size_t k) {
      const Ent& e = *ent[k];
      if (bt == BT_QG) {
        double s = sigL(k);
        if (!dp.stem_qg2) return s;
        double s2 = sigR(k);
        if (s2 < HALF_NEG) return s;
        return log_add(s, s2);
      }
      if (bt == BT_QGL) {
        double s = sigL(k);
        if (!dp.stem_both) return s;
        return log_add(s, sigR(k));
      }
      if (!dp.stem_both) return dp.four_score(e.p, e.u, grp[k]->v, e.a);
      return sigR(k);
    }
  };

  // --stem-est rank: drop infeasible states, rank the rest by the stem score (step: whole
  // table, group: per group v; threshold rule of keep_mask) (DP03._prune_stem).
  void prune4_stem(Tab4& t, int bt, int q, int gtab, int beam) {
    const auto t_qg = std::chrono::steady_clock::now();
    StemEval ev(*this, t, bt, q);
    const size_t N = ev.ent.size();
    std::vector<char> keep(ev.ok);
    size_t ranked = 0;
    if (beam > 0) {
      std::vector<std::pair<size_t, size_t>> spans;
      if (gap_step) {
        spans.emplace_back(0, N);
      } else {
        size_t lo = 0;
        for (const Grp& g : t.grps) {
          spans.emplace_back(lo, lo + g.ents.size());
          lo += g.ents.size();
        }
      }
      std::vector<size_t> idx;
      std::vector<double> sc, xs, sl, sr, so;
      std::vector<char> km;
      const bool ext = bt == BT_QG && stem_ext;
      for (const auto& sp : spans) {
        idx.clear();
        for (size_t k = sp.first; k < sp.second; ++k)
          if (ev.ok[k]) idx.push_back(k);
        if (idx.size() <= static_cast<size_t>(beam)) continue;
        ranked += idx.size();
        sc.clear();
        if (pr_on) {
          const auto tp = std::chrono::steady_clock::now();
          const int ti = pr_tab(bt);
          sc.assign(idx.size(), NEG_INF);
          for (size_t m : gap_order(idx.size(), [&](size_t m) {
                 return std::make_pair(ev.grp[idx[m]]->v, ev.ent[idx[m]]->u);
               })) {
            const size_t k = idx[m];
            sc[m] = pr_score4(bt, *ev.ent[k], ev.grp[k]->v, q);
          }
          pr_states[ti] += static_cast<long long>(idx.size());
          pr_sec[ti] += std::chrono::duration<double>(std::chrono::steady_clock::now() - tp).count();
        } else if (ext && stem_profile) {
          const auto t_x = std::chrono::steady_clock::now();
          sc.assign(idx.size(), NEG_INF);
          for (size_t m : gap_order(idx.size(), [&](size_t m) {
                 return std::make_pair(ev.grp[idx[m]]->v, ev.ent[idx[m]]->u);
               })) {
            const size_t k = idx[m];
            const Ent& e = *ev.ent[k];
            const int v = ev.grp[k]->v;
            const double s = e.a + prof_eval(qg_prof(e.p), e.u, v);
            sc[m] = s;
            ext_finite += s > HALF_NEG;
          }
          const double dt =
              std::chrono::duration<double>(std::chrono::steady_clock::now() - t_x).count();
          ext_sec += dt;
          pf_total_sec += dt;
          ext_states += static_cast<long long>(idx.size());
          pf_states += static_cast<long long>(idx.size());
        } else {
          for (size_t k : idx) sc.push_back(ev.score(k));
        }
        if (ext && !stem_profile && !pr_on) {
          const auto t_x = std::chrono::steady_clock::now();
          xs.clear();
          for (size_t k : idx) {
            const Ent& e = *ev.ent[k];
            xs.push_back(stem_ext_score(e.p, e.u, ev.grp[k]->v, e.a));
          }
          ext_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_x).count();
          ext_states += static_cast<long long>(idx.size());
          for (size_t m = 0; m < idx.size(); ++m) {
            ext_finite += xs[m] > HALF_NEG;
            sc[m] = log_add(sc[m], xs[m]);
          }
        }
        keep_mask(sc, km, beam);
        for (size_t m = 0; m < idx.size(); ++m) keep[idx[m]] = km[m];
      }
    }
    if (bt == BT_QG)
      stem_qg_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_qg).count();
    size_t n_ok = 0, n_keep = 0;
    for (size_t k = 0; k < N; ++k) {
      n_ok += ev.ok[k];
      n_keep += keep[k];
    }
    stem_removed[bt == BT_QG ? 0 : bt == BT_QGL ? 1 : 2] += N - n_ok;
    if (gtab >= 0) {
      gate_stats.cand[gtab] += N;
      gate_stats.elig[gtab] += n_ok;
      gate_stats.ranked[gtab] += ranked;
      gate_stats.kept[gtab] += n_keep;
    }
    if (n_keep < N) filter_step(t, keep);
    if (bt == BT_QGL) stem_add_A(q);
  }



  // Top `beam` of `scores`; ties at the threshold go to the earliest (DP03._keep_mask).
  static void keep_mask(const std::vector<double>& scores, std::vector<char>& keep, int beam) {
    std::vector<double> tmp(scores);
    std::nth_element(tmp.begin(), tmp.begin() + (beam - 1), tmp.end(), std::greater<double>());
    const double thr = tmp[beam - 1];
    int room = beam;
    for (double s : scores) room -= s > thr;
    keep.assign(scores.size(), 0);
    for (size_t k = 0; k < scores.size(); ++k) {
      if (scores[k] > thr) {
        keep[k] = 1;
      } else if (scores[k] == thr && room > 0) {
        keep[k] = 1;
        --room;
      }
    }
  }

  // Keeps g.ents[r] with keep[off + r], in insertion order; rebuilds the index.
  // Releases the pre-prune capacity: a group can hold many times the beam
  // before pruning and would otherwise keep that memory for the whole run.
  static void filter_group(Grp& g, const std::vector<char>& keep, size_t off) {
    std::vector<Ent> ents;
    for (size_t r = 0; r < g.ents.size(); ++r)
      if (keep[off + r]) ents.push_back(g.ents[r]);
    ents.shrink_to_fit();
    g.ents.swap(ents);
    FastMap<uint64_t, int> idx;
    idx.reserve(g.ents.size());
    for (int r = 0; r < static_cast<int>(g.ents.size()); ++r)
      idx.emplace(key2(g.ents[r].p, g.ents[r].u), r);
    g.idx.swap(idx);
  }

  // Top `beam` states of the whole table at this step (DP03._prune, gap_beam="step").
  // Kept states stay in insertion order; emptied groups are dropped.
  void prune4_step(Tab4& t, int gtab, int beam) {
    size_t total = 0;
    for (const Grp& g : t.grps) total += g.ents.size();
    if (gtab >= 0) {
      gate_stats.cand[gtab] += total;
      gate_stats.elig[gtab] += total;
    }
    if (total <= static_cast<size_t>(beam)) {
      if (gtab >= 0) gate_stats.kept[gtab] += total;
      return;
    }
    std::vector<double> scores;
    scores.reserve(total);
    for (const Grp& g : t.grps)
      for (const Ent& e : g.ents)
        scores.push_back(four_score(e.p, e.u, g.v, e.a));
    std::vector<char> keep;
    keep_mask(scores, keep, beam);
    if (gtab >= 0) {
      gate_stats.ranked[gtab] += total;
      gate_stats.kept[gtab] += std::min(total, static_cast<size_t>(beam));
    }
    filter_step(t, keep);
  }





  // Keeps the flattened states with keep[k] (groups in order, entries in insertion order);
  // emptied groups are dropped.
  void filter_step(Tab4& t, const std::vector<char>& keep) {
    std::vector<Grp> grps;
    t.vidx.clear();
    size_t off = 0;
    for (Grp& g : t.grps) {
      size_t sz = g.ents.size();
      filter_group(g, keep, off);
      off += sz;
      if (g.ents.empty()) continue;
      t.vidx.emplace(g.v, static_cast<int>(grps.size()));
      grps.push_back(std::move(g));
    }
    t.grps.swap(grps);
  }

  // Beam on a right-end map {left: value} by prefix score (DP03._prune_1d).
  void prune1(OMap& at) {
    const int beam = tbeam[BT_PKM];
    if (beam <= 0 || at.live <= beam) return;
    std::vector<std::pair<double, int>> ranked;
    ranked.reserve(at.live);
    for (size_t r = 0; r < at.slots(); ++r) {
      int i = at.keys[r];
      if (i < 0) continue;
      ranked.emplace_back(prefix_log(i) + at.vals[r], i);
    }
    std::sort(ranked.begin(), ranked.end(),
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                if (a.first != b.first) return a.first > b.first;
                return a.second > b.second;
              });
    std::vector<char> keep(n, 0);
    for (int k = 0; k < beam; ++k) keep[ranked[k].second] = 1;
    for (size_t r = 0; r < at.slots(); ++r) {
      int i = at.keys[r];
      if (i >= 0 && !keep[i]) at.erase(i);
    }
    at.compact();
  }

  // extra_row / extra_at: Qm2, dropped with Qm (bt = BT_QM).
  void prune2(AtView& table, OMap& at, int j, int bt, AtView* extra_row, OMap* extra_at) {
    const int beam = tbeam[bt];
    if (beam <= 0 || at.live <= beam) {
      return;
    }
    std::vector<std::pair<double, int>> ranked;
    ranked.reserve(at.live);
    // --stem-est-principled: Qp by alpha + log Ctx(i, j)
    const bool pq = pr_on && bt == BT_QP;
    const auto tp = std::chrono::steady_clock::now();
    for (size_t r = 0; r < at.slots(); ++r) {
      int i = at.keys[r];
      if (i < 0) continue;
      if (pq) {
        const double s = at.vals[r] + pr_ctx(i, j);
        ranked.emplace_back(s < HALF_NEG ? NEG_INF : s, i);
      } else {
        ranked.emplace_back(prefix_log(i) + at.vals[r], i);
      }
    }
    if (pq) {
      pr_states[5] += at.live;
      pr_sec[5] += std::chrono::duration<double>(std::chrono::steady_clock::now() - tp).count();
    }
    std::sort(ranked.begin(), ranked.end(),
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                if (a.first != b.first) return a.first > b.first;
                return a.second > b.second;
              });
    std::vector<char> keep(n, 0);
    for (int k = 0; k < beam; ++k) keep[ranked[k].second] = 1;
    for (size_t r = 0; r < at.slots(); ++r) {
      int i = at.keys[r];
      if (i < 0 || keep[i]) continue;
      at.erase(i);
      table[ix(i, j)] = NEG_INF;
      if (extra_row) (*extra_row)[ix(i, j)] = NEG_INF;
      if (extra_at) extra_at->erase(i);
    }
    at.compact();
    if (extra_at) extra_at->compact();
  }

  void set_pair(AtView& table, std::vector<OMap>& at, int i, int j, double val) {
    if (val < HALF_NEG) return;
    double t = log_add(static_cast<double>(table[ix(i, j)]), val);
    table[ix(i, j)] = t;
    at[j].set(i, t);
  }

  void set_Qm(int i, int j, double val, bool is_m2) {
    if (val < HALF_NEG) return;
    if (dz) {
      rsQm[j & 1].P.add(i, val);
      if (is_m2) rsQm2[j & 1].P.add(i, val);
      return;
    }
    double t = log_add(static_cast<double>(Qm[ix(i, j)]), val);
    Qm[ix(i, j)] = t;
    Qm_at[j].set(i, t);
    if (is_m2) {
      double t2 = log_add(static_cast<double>(Qm2[ix(i, j)]), val);
      Qm2[ix(i, j)] = t2;
      Qm2_at[j].set(i, t2);
    }
  }

