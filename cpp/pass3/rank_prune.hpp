// Pass-3 pruning by ranking score: stem profiles (prof_*), prune4_helper / prune4_stem /
// prune4_gated / prune4_step, score ranks, drop notes, prune1 / prune2 (top-b of a column).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  // --stem-est-profile-eps: keep the smallest set of entries (largest first, ties by index)
  // holding 1 - eps of the part's weight, in index order; returns the entries dropped.
  static long long prof_trunc(std::vector<int>& ix, std::vector<double>& w, double eps) {
    const size_t m = w.size();
    if (m <= 1) return 0;
    double tot = NEG_INF;
    for (double x : w) tot = log_add(tot, x);
    const double target = tot + std::log1p(-eps);
    std::vector<size_t> order(m);
    for (size_t k = 0; k < m; ++k) order[k] = k;
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return w[a] > w[b]; });
    std::vector<char> keep(m, 0);
    double acc = NEG_INF;
    for (size_t k : order) {
      keep[k] = 1;
      acc = log_add(acc, w[k]);
      if (acc >= target) break;
    }
    size_t o = 0;
    for (size_t k = 0; k < m; ++k)
      if (keep[k]) {
        ix[o] = ix[k];
        w[o] = w[k];
        ++o;
      }
    ix.resize(o);
    w.resize(o);
    return static_cast<long long>(m - o);
  }

  void prof_finish_key(Prof& p, int kind) {
    if (prof_eps > 0.0) pf_trunc += prof_trunc(p.ly, p.lw, prof_eps) + prof_trunc(p.rx, p.rw, prof_eps);
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

  // --stem-est-ml-tail nested: log sum over the case-(f) closings j > q+1+L_ml of the Qgls key
  // (i, seed) of w_f(i-1, j) + beta_N(i-1, j) (Result 20's term), j ascending.
  double ml_tail_sum(int i, bool seed, int q) {
    const int x = i - 1;
    if (x < 0) return NEG_INF;
    double acc = NEG_INF;
    for (int j = q + 2 + stem_ml_span; j <= n - 1; ++j) {
      const int unpaired = j - q - 1;
      if (seed && unpaired > cap) continue;
      if (!E.pairs(x, j) || !E.end_ok(x, j)) continue;
      const double b = stem_->B(x, j);
      if (b < HALF_NEG) continue;
      double w_ext = E.w_ap + E.w_bp + E.w_cp * unpaired + pair_bonus;
      if (dz) w_ext += E.run(q + 1, j - 1);
      acc = log_add(acc, (w_ext + E.at(x, j)) + b);
      ++pf_tail_cl;
    }
    return acc;
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
    if (ml_tail) k.tail = ml_tail_sum(i, seed, q);
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
      if (stem_consistent) L += E.at(u, v);
    }
    double R = NEG_INF, Ro = NEG_INF;
    if (E.end_ok(u, v) && bg_reorder) R = prof_R_reorder(p, u, v);
    if (E.end_ok(u, v) && (!bg_reorder || ro_check))
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
    if (!bg_reorder) R = Ro;
    double s = log_add(L, R);
    if (p.tail > HALF_NEG) s = log_add(s, p.tail + qz);
    if (ro_check) {
      ro_off = log_add(L, Ro);
      if (p.tail > HALF_NEG) ro_off = log_add(ro_off, p.tail + qz);
    }
    pf_evals += static_cast<long long>(p.ly.size() + p.rx.size());
    return s;
  }

  // --stem-est-bg-reorder: the M row of key profile p at v (cached per step) from e = lo up,
  // built over the R entries x' in index order and, per x', the e of stE[x'] in [lo, v - 1]
  // ascending.  States are scored by gap (v, u), so the first lo = u + 1 serves the rest; a
  // smaller u rebuilds the row (M(e, v) of each e does not depend on lo).
  const MRow& prof_mrow(const Prof& p, int v, int lo) {
    std::unordered_map<int, MRow>& byv = stM[&p];
    auto it = byv.find(v);
    if (it != byv.end() && it->second.lo <= lo) return it->second;
    ++ro_rows;
    if (mAcc.size() != static_cast<size_t>(n + 1)) {
      mAcc.assign(n + 1, NEG_INF);
      mStamp.assign(n + 1, 0);
    }
    ++m_stamp;
    mTouched.clear();
    for (size_t k = 0; k < p.rx.size(); ++k) {
      const int x = p.rx[k];
      auto ie = stE.find(x);
      if (ie == stE.end()) continue;
      const std::vector<int>& es = ie->second;
      for (auto ke = std::lower_bound(es.begin(), es.end(), lo); ke != es.end(); ++ke) {
        const int e = *ke;
        if (e > v - 1) break;
        ++ro_build;
        const double b = stem_B_cached(x, e, v);
        if (b < HALF_NEG) continue;
        if (mStamp[e] != m_stamp) {
          mStamp[e] = m_stamp;
          mAcc[e] = NEG_INF;
          mTouched.push_back(e);
        }
        mAcc[e] = log_add(mAcc[e], p.rw[k] + b);
      }
    }
    std::sort(mTouched.begin(), mTouched.end());
    MRow r;
    r.lo = lo;
    r.es = mTouched;
    r.ms.reserve(mTouched.size());
    for (int e : mTouched) r.ms.push_back(mAcc[e]);
    MRow& slot = byv[v];
    slot = std::move(r);
    return slot;
  }

  // log sum_x' P_R(x') Bg(x', u, v) = beta2 + at(u, v) + log sum_{e in [u+1, v-1]}
  // Qz[u+1][e-1] M(e, v) (end_ok(u, v) checked by the caller).
  double prof_R_reorder(const Prof& p, int u, int v) {
    if (p.rx.empty()) return NEG_INF;
    const MRow& r = prof_mrow(p, v, u + 1);
    double acc = NEG_INF;
    for (size_t k = std::lower_bound(r.es.begin(), r.es.end(), u + 1) - r.es.begin();
         k < r.es.size(); ++k) {
      ++ro_eval;
      const double qz = qz_log(u + 1, r.es[k] - 1);
      if (qz < HALF_NEG) continue;
      acc = log_add(acc, qz + r.ms[k]);
    }
    if (acc < HALF_NEG) return NEG_INF;
    return acc + (E.w_stem + E.at(u, v));
  }

  // --stem-bg-reorder-check: one scored state (reordered s, original ro_off, alpha a), stored
  // at off[m]
  void ro_cmp(double s, double a, std::vector<double>& off, size_t m) {
    const double so = a + ro_off;
    off[m] = so;
    ++roc_states;
    if ((s > HALF_NEG) != (so > HALF_NEG)) ++roc_mis;
    else if (s > HALF_NEG) roc_rel = std::max(roc_rel, std::fabs(std::expm1(s - so)));
  }

  // --stem-bg-reorder-check: keep decisions of one ranking event under the original scores
  // (order of the event's score vector) against km.
  void ro_flips(std::vector<double>& off, const std::vector<char>& km, int beam) {
    std::vector<char> ko;
    keep_mask(off, ko, beam);
    ++roc_events;
    long long f = 0;
    for (size_t r = 0; r < km.size(); ++r) f += ko[r] != km[r];
    roc_flip_states += f;
    roc_flip_events += f > 0;
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

  // --stem-profile-check (t = 0 Qg, 1 helpers)
  void prof_cmp(int t, double s, double ref) {
    ++pc_states[t];
    if ((ref > HALF_NEG) != (s > HALF_NEG)) ++pc_mis[t];
    else if (ref > HALF_NEG) pc_rel[t] = std::max(pc_rel[t], std::fabs(std::expm1(s - ref)));
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
    if (pr_on && pr_arm_gaponly) {
      std::vector<std::pair<int, int>> iv;
      for (const auto& sp : spans)
        if (sp.second - sp.first > static_cast<size_t>(beam))
          for (size_t k = sp.first; k < sp.second; ++k) iv.emplace_back(ent[k]->u + 1, gv[k] - 1);
      pr_gap_begin(iv);
    }
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
          if (pr_dump)
            std::fprintf(pr_dump, "%d %d %d %d %d %.17g\n", bt, q, ent[k]->p, ent[k]->u, gv[k], sc[r]);
          hlp_finite += sc[r] > HALF_NEG;
          if (pr_check) prc_cmp(ti, sc[r], pr_naive4(bt, *ent[k], gv[k], q));
        }
        pr_states[ti] += static_cast<long long>(m);
        pr_sec[ti] += std::chrono::duration<double>(std::chrono::steady_clock::now() - tp).count();
      } else if (stem_profile) {
        sc.assign(m, NEG_INF);
        if (ro_check) so.assign(m, NEG_INF);
        for (size_t r : gap_order(m, [&](size_t r) {
               return std::make_pair(gv[sp.first + r], ent[sp.first + r]->u);
             })) {
          const size_t k = sp.first + r;
          const Ent& e = *ent[k];
          const Prof& kp = ls ? ls_prof(e.p, gv[k] == q, q) : rs_prof(e.p, q);
          const double s = e.a + prof_eval(kp, e.u, gv[k]);
          sc[r] = s;
          if (ro_check) ro_cmp(s, e.a, so, r);
          hlp_closings += kp.ncl;
          hlp_finite += s > HALF_NEG;
          ++pf_states;
          if (prof_check) {
            if (ls) closings_ls(e.p, gv[k], q, cl);
            else closings_rs(e.p, q, cl);
            prof_cmp(1, s, stem_chain ? chain_naive(cl, e.u, gv[k], e.a)
                                      : helper_sigma(cl, e.u, gv[k], e.a));
          }
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
          if (ext_check) {
            const double nv = helper_nested_naive(ls, e, gv[k], q);
            ++hc_states;
            if ((nv > HALF_NEG) != (s > HALF_NEG)) ++hc_mis;
            else if (nv > HALF_NEG) hc_rel = std::max(hc_rel, std::fabs(std::expm1(s - nv)));
          }
          continue;
        }
        if (ls) closings_ls(e.p, gv[k], q, cl);
        else closings_rs(e.p, q, cl);
        const double s = helper_sigma(cl, e.u, gv[k], e.a);
        sc.push_back(s);
        hlp_closings += static_cast<long long>(cl.size());
        hlp_finite += s > HALF_NEG;
        if (ext_check) {
          const double nv = helper_naive(cl, e.u, gv[k], e.a);
          ++hc_states;
          if ((nv > HALF_NEG) != (s > HALF_NEG)) ++hc_mis;
          else if (nv > HALF_NEG) hc_rel = std::max(hc_rel, std::fabs(std::expm1(s - nv)));
        }
      }
      hlp_states[ls ? 0 : 1] += static_cast<long long>(m);
      keep_mask(sc, km, beam);
      if (ro_check && stem_profile) ro_flips(so, km, beam);
      for (size_t r = 0; r < m; ++r) {
        keep[sp.first + r] = km[r];
        any |= !km[r];
      }
      if (drops) {
        std::vector<size_t> rank = score_ranks(sc);
        ++drops->events[bt];
        for (size_t r = 0; r < m; ++r)
          if (!km[r]) {
            const Ent& e = *ent[sp.first + r];
            const int v = gv[sp.first + r];
            note_drop(bt, bt, q, gap_step ? -1 : v, e.p, e.u, v, q, rank[r], m, sc[r]);
          } else if (drops->kf) {
            const Ent& e = *ent[sp.first + r];
            const int v = gv[sp.first + r];
            note_kept(bt, q, gap_step ? -1 : v, e.p, e.u, v, q, rank[r], m, sc[r]);
          }
      }
    }
    if (any) filter_step(t, keep);
    const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    hlp_sec += dt;
    if (stem_profile) pf_total_sec += dt;
  }

  // Ranking order of keep_mask: descending, earliest first on ties.
  static std::vector<size_t> order_desc(const std::vector<double>& s) {
    std::vector<size_t> order(s.size());
    for (size_t k = 0; k < s.size(); ++k) order[k] = k;
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return s[a] > s[b]; });
    return order;
  }

  // --stem-est-qg-combine split (DP03._split_mask): for each of parts (sigma_L, sigma_R,
  // sigma_ext) in turn, up to ceil(beam/3) more states by that score in order_desc (finite
  // scores only, chosen states skipped), then up to beam states by the combined score sc.
  static void split_mask(const std::vector<double>* const parts[3], const std::vector<double>& sc,
                         std::vector<char>& keep, int beam) {
    keep.assign(sc.size(), 0);
    const int m = (beam + 2) / 3;
    int taken = 0;
    for (int r = 0; r < 3; ++r) {
      const std::vector<double>& s = *parts[r];
      int got = 0;
      for (size_t k : order_desc(s)) {
        if (got == m || taken == beam || s[k] < HALF_NEG) break;
        if (keep[k]) continue;
        keep[k] = 1;
        ++got;
        ++taken;
      }
    }
    for (size_t k : order_desc(sc)) {
      if (taken == beam) break;
      if (keep[k]) continue;
      keep[k] = 1;
      ++taken;
    }
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
        if (pr_on && (pr_gsel || pr_gcheck))
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
    std::unordered_map<uint64_t, double> hcache, hacache;
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

    // The L-stem score's H: H itself, or Hav with --stem-est-hr av (ok[] keeps H).
    double Hs(int u, int v) {
      if (!dp.hr_av) return H(u, v);
      uint64_t key = key2(u, v);
      auto it = hacache.find(key);
      if (it != hacache.end()) return it->second;
      double h = dp.stem_->Hav(u, v, q);
      hacache.emplace(key, h);
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
      double h = Hs(e.u, grp[k]->v);
      if (h < HALF_NEG) return NEG_INF;
      return dp.four_score(e.p, e.u, grp[k]->v, e.a) + h + cq + w_pk;
    }

    // Stem-1 score (state as the left stem, H_R); NEG_INF if none.
    double sigL(size_t k) {
      const Energy& E = dp.E;
      const Ent& e = *ent[k];
      if (bt == BT_QGR) return NEG_INF;
      double s = s1h(k);
      if (!dp.stem_consistent || s < HALF_NEG) return s;
      if (bt == BT_QG) return s + (E.at(e.u, grp[k]->v) + E.at(e.p, q));
      return s + (E.at(e.p, q) - E.w_stem);
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
          if (dp.stem_consistent) s2 += (E.at(e.p, q) + E.w_beta1);
          return s2;
        }
        double h2 = dp.stem_->H2(e.p, e.u, v);
        if (h2 < HALF_NEG) return NEG_INF;
        double s2 = sq + e.a;
        if (dp.score_qz) s2 += dp.qz_log(e.u + 1, v - 1);
        s2 = s2 + h2 - dp.stem_->S[e.p] + w_pk;
        if (dp.stem_consistent) s2 += (E.at(e.u, v) + E.at(e.p, q));
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
      if (pr_on && pr_arm_gaponly) {
        std::vector<std::pair<int, int>> iv;
        for (const auto& sp : spans) {
          size_t m = 0;
          for (size_t k = sp.first; k < sp.second; ++k) m += ev.ok[k] ? 1 : 0;
          if (m <= static_cast<size_t>(beam)) continue;
          for (size_t k = sp.first; k < sp.second; ++k)
            if (ev.ok[k]) iv.emplace_back(ev.ent[k]->u + 1, ev.grp[k]->v - 1);
        }
        pr_gap_begin(iv);
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
            if (pr_dump)
              std::fprintf(pr_dump, "%d %d %d %d %d %.17g\n", bt, q, ev.ent[k]->p, ev.ent[k]->u,
                           ev.grp[k]->v, sc[m]);
            if (pr_check) prc_cmp(ti, sc[m], pr_naive4(bt, *ev.ent[k], ev.grp[k]->v, q));
          }
          pr_states[ti] += static_cast<long long>(idx.size());
          pr_sec[ti] += std::chrono::duration<double>(std::chrono::steady_clock::now() - tp).count();
        } else if (ext && stem_profile) {
          const auto t_x = std::chrono::steady_clock::now();
          sc.assign(idx.size(), NEG_INF);
          if (ro_check) so.assign(idx.size(), NEG_INF);
          for (size_t m : gap_order(idx.size(), [&](size_t m) {
                 return std::make_pair(ev.grp[idx[m]]->v, ev.ent[idx[m]]->u);
               })) {
            const size_t k = idx[m];
            const Ent& e = *ev.ent[k];
            const int v = ev.grp[k]->v;
            const double s = e.a + prof_eval(qg_prof(e.p), e.u, v);
            sc[m] = s;
            if (ro_check) ro_cmp(s, e.a, so, m);
            ext_finite += s > HALF_NEG;
            if (prof_check && stem_chain)
              prof_cmp(0, s, chain_naive({{e.p, stc_q, 0.0}}, e.u, v, e.a));
            else if (prof_check)
              prof_cmp(0, s, log_add(ev.score(k), stem_ext_score(e.p, e.u, v, e.a)));
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
        if ((hr_av || hr_check) && bt != BT_QGR)
          for (size_t k : idx) {
            const int u = ev.ent[k]->u, v = ev.grp[k]->v;
            if (hr_check) {
              const double hs = ev.Hs(u, v), nv = hr_naive(u, v, q);
              ++hrc_states;
              if ((hs > HALF_NEG) != (nv > HALF_NEG)) ++hrc_mis;
              else if (nv > HALF_NEG) hrc_rel = std::max(hrc_rel, std::fabs(std::expm1(hs - nv)));
            }
            if (hr_av && bt == BT_QG) {
              const double lo = ev.H(u, v), la = ev.Hs(u, v);
              if (lo > HALF_NEG && la > HALF_NEG) {
                const double d = (la + ev.cq) - (lo - ev.c0q);
                hr_dmin = hr_n == 0 ? d : std::min(hr_dmin, d);
                hr_dmax = hr_n == 0 ? d : std::max(hr_dmax, d);
                hr_dsum += d;
                ++hr_n;
              } else if (la > HALF_NEG) {
                ++hr_av_only;
              } else if (lo > HALF_NEG) {
                ++hr_orig_only;
              }
            }
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
          if (ext_check)
            for (size_t m = 0; m < idx.size(); ++m) {
              const Ent& e = *ev.ent[idx[m]];
              const double nv = stem_ext_naive(e.p, e.u, ev.grp[idx[m]]->v, e.a);
              ++xc_states;
              if ((nv > HALF_NEG) != (xs[m] > HALF_NEG)) ++xc_mis;
              else if (nv > HALF_NEG) xc_rel = std::max(xc_rel, std::fabs(std::expm1(xs[m] - nv)));
            }
        }
        if (ext && stem_split && !pr_on) {
          sl.clear();
          sr.clear();
          for (size_t k : idx) {
            sl.push_back(ev.sigL(k));
            sr.push_back(ev.sigR(k));
          }
          const std::vector<double>* parts[3] = {&sl, &sr, &xs};
          split_mask(parts, sc, km, beam);
        } else {
          keep_mask(sc, km, beam);
          if (ro_check && ext && stem_profile) ro_flips(so, km, beam);
        }
        for (size_t m = 0; m < idx.size(); ++m) keep[idx[m]] = km[m];
        if (drops) {
          std::vector<size_t> rank = score_ranks(sc);
          ++drops->events[bt];
          for (size_t m = 0; m < idx.size(); ++m)
            if (!km[m]) {
              const Ent& e = *ev.ent[idx[m]];
              const int v = ev.grp[idx[m]]->v;
              note_drop(bt, bt, q, gap_step ? -1 : v, e.p, e.u, v, q, rank[m], idx.size(), sc[m]);
            } else if (drops->kf) {
              const Ent& e = *ev.ent[idx[m]];
              const int v = ev.grp[idx[m]]->v;
              note_kept(bt, q, gap_step ? -1 : v, e.p, e.u, v, q, rank[m], idx.size(), sc[m]);
            }
        }
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

  // Gate selection over `n_all` candidates (one group, or the whole table in step mode):
  // (a) drop g < tau; (b) if more than `beam` remain, keep the top `beam` of them by the
  // pruning score (threshold rule of keep_mask); (c) add up to `reserve` more eligible
  // states by g, descending, ties by insertion order.  Returns false if every candidate
  // is kept (the caller then leaves the states untouched).
  template <class ScoreAt>
  bool gate_select(const std::vector<double>& gv, ScoreAt score_at, std::vector<char>& keep,
                   int gtab, int beam) {
    const size_t n_all = gv.size();
    std::vector<char> elig(n_all, 0);
    size_t ne = 0;
    for (size_t k = 0; k < n_all; ++k) {
      elig[k] = gv[k] >= gcfg.tau;
      ne += elig[k];
    }
    gate_stats.cand[gtab] += n_all;
    gate_stats.elig[gtab] += ne;
    if (beam <= 0 || ne <= static_cast<size_t>(beam)) {
      gate_stats.kept[gtab] += ne;
      if (ne == n_all) return false;
      keep.swap(elig);
      return true;
    }
    std::vector<double> scores;
    scores.reserve(ne);
    for (size_t k = 0; k < n_all; ++k)
      if (elig[k]) scores.push_back(score_at(k));
    gate_stats.ranked[gtab] += ne;
    std::vector<char> km;
    keep_mask(scores, km, beam);
    keep.assign(n_all, 0);
    size_t m = 0, kept = 0;
    for (size_t k = 0; k < n_all; ++k)
      if (elig[k]) {
        keep[k] = km[m++];
        kept += keep[k];
      }
    if (gcfg.reserve > 0) {
      std::vector<size_t> rest;
      for (size_t k = 0; k < n_all; ++k)
        if (elig[k] && !keep[k]) rest.push_back(k);
      std::stable_sort(rest.begin(), rest.end(),
                       [&](size_t a, size_t b) { return gv[a] > gv[b]; });
      size_t add = std::min(rest.size(), static_cast<size_t>(gcfg.reserve));
      for (size_t r = 0; r < add; ++r) keep[rest[r]] = 1;
      kept += add;
      gate_stats.reserved[gtab] += add;
    }
    gate_stats.kept[gtab] += kept;
    return true;
  }

  void prune4_gated(Tab4& t, int q, bool use_arm, int gtab, int beam) {
    const Gate& G = *gate_;
    std::vector<double> gv;
    std::vector<char> keep;
    if (!gap_step) {
      for (auto& g : t.grps) {
        gv.clear();
        for (const Ent& e : g.ents) gv.push_back(G.g(e.p, e.u, g.v, q));
        auto score_at = [&](size_t k) {
          const Ent& e = g.ents[k];
          return use_arm ? arm_score_of(e.p, e.u, g.v, q, e.a) : four_score(e.p, e.u, g.v, e.a);
        };
        if (gate_select(gv, score_at, keep, gtab, beam)) filter_group(g, keep, 0);
      }
      return;
    }
    std::vector<std::pair<int, int>> where;  // (group, entry) of each flattened candidate
    for (int gi = 0; gi < static_cast<int>(t.grps.size()); ++gi) {
      const Grp& g = t.grps[gi];
      for (int r = 0; r < static_cast<int>(g.ents.size()); ++r) {
        const Ent& e = g.ents[r];
        gv.push_back(G.g(e.p, e.u, g.v, q));
        where.emplace_back(gi, r);
      }
    }
    auto score_at = [&](size_t k) {
      const Grp& g = t.grps[where[k].first];
      const Ent& e = g.ents[where[k].second];
      return use_arm ? arm_score_of(e.p, e.u, g.v, q, e.a) : four_score(e.p, e.u, g.v, e.a);
    };
    if (gate_select(gv, score_at, keep, gtab, beam)) filter_step(t, keep);
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
    std::unordered_map<uint64_t, int> idx;
    idx.reserve(g.ents.size());
    for (int r = 0; r < static_cast<int>(g.ents.size()); ++r)
      idx.emplace(key2(g.ents[r].p, g.ents[r].u), r);
    g.idx.swap(idx);
  }

  // Top `beam` states of the whole table at this step (DP03._prune, gap_beam="step").
  // Kept states stay in insertion order; emptied groups are dropped.
  void prune4_step(Tab4& t, int bt, int q, bool use_arm, int gtab, int beam) {
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
        scores.push_back(use_arm ? arm_score_of(e.p, e.u, g.v, q, e.a)
                                 : four_score(e.p, e.u, g.v, e.a));
    std::vector<char> keep;
    keep_mask(scores, keep, beam);
    if (gtab >= 0) {
      gate_stats.ranked[gtab] += total;
      gate_stats.kept[gtab] += std::min(total, static_cast<size_t>(beam));
    }
    if (drops) {
      std::vector<size_t> rank = score_ranks(scores);
      ++drops->events[bt];
      size_t k = 0;
      for (const Grp& g : t.grps)
        for (const Ent& e : g.ents) {
          if (!keep[k])
            note_drop(bt, bt, q, -1, e.p, e.u, g.v, q, rank[k], total, scores[k]);
          ++k;
        }
    }
    filter_step(t, keep);
  }

  // Rank of each score (0 = best) in keep_mask's order: descending, earliest first on ties.
  static std::vector<size_t> score_ranks(const std::vector<double>& s) {
    std::vector<size_t> order(s.size()), rank(s.size());
    for (size_t k = 0; k < s.size(); ++k) order[k] = k;
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return s[a] > s[b]; });
    for (size_t k = 0; k < order.size(); ++k) rank[order[k]] = k;
    return rank;
  }

  // One dropped candidate (table bt, dropped by the ranking of table src at step j; v = its
  // group in group mode, else -1).  With an exact state file only its states are written.
  void note_drop(int bt, int src, int j, int v, int a, int b, int c, int d, size_t rank,
                 size_t cands, double score) {
    DropLog& L = *drops;
    ++L.dropped[bt];
    double w = std::nan("");
    if (L.exact) {
      auto it = L.exact->find(state_key(bt, a, b, c, d));
      if (it == L.exact->end()) return;
      w = it->second;
      ++L.matched[bt];
      L.w[bt] += w;
      L.wpos[bt][std::min(9, 10 * j / L.n)] += w;
      if (w < L.min_w) return;
    }
    std::fprintf(L.f, "D %s %s %s %d %d %d %d %d %d %zu %zu %.9g %.6e\n", L.id.c_str(),
                 BT_NAMES[bt], BT_NAMES[src], j, v, a, b, c, d, rank, cands, score, w);
  }

  // Pass-2 / l2r table drop (R,U,S,P1-P3,…): same D line shape; table/src are l2r names.
  void note_l2r_drop(const char* tname, int j, int v, int a, int b, int c, int d, size_t rank,
                     size_t cands, double score) {
    if (!drops) return;
    std::fprintf(drops->f, "D %s %s %s %d %d %d %d %d %d %zu %zu %.9g nan\n", drops->id.c_str(),
                 tname, tname, j, v, a, b, c, d, rank, cands, score);
  }

  // --dump-kept: one kept candidate of a stem / helper ranking event (fields of note_drop;
  // w = exact w, nan if not in the exact state file or without one).
  void note_kept(int bt, int j, int v, int a, int b, int c, int d, size_t rank, size_t cands,
                 double score) {
    DropLog& L = *drops;
    double w = std::nan("");
    if (L.exact) {
      auto it = L.exact->find(state_key(bt, a, b, c, d));
      if (it != L.exact->end()) w = it->second;
    }
    std::fprintf(L.kf, "K %s %s %s %d %d %d %d %d %d %zu %zu %.9g %.6e\n", L.id.c_str(),
                 BT_NAMES[bt], BT_NAMES[bt], j, v, a, b, c, d, rank, cands, score, w);
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
  void prune1(OMap& at, int j) {
    const int beam = tbeam[BT_PKM];
    if (remove_keys) {
      for (size_t r = 0; r < at.slots(); ++r) {
        int d = at.keys[r];
        if (d >= 0 && remove_keys->count(state_key(BT_PKM, d, j, -1, -1))) at.erase(d);
      }
      at.compact();
    }
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
    if (drops) {
      ++drops->events[BT_PKM];
      for (size_t k = beam; k < ranked.size(); ++k)
        note_drop(BT_PKM, BT_PKM, j, -1, ranked[k].second, j, -1, -1, k, ranked.size(),
                  ranked[k].first);
    }
    for (size_t r = 0; r < at.slots(); ++r) {
      int i = at.keys[r];
      if (i >= 0 && !keep[i]) at.erase(i);
    }
    at.compact();
  }

  // extra_row / extra_at: Qm2, dropped with Qm (bt = BT_QM).
  void prune2(AtView& table, OMap& at, int j, int bt, AtView* extra_row, OMap* extra_at) {
    const int beam = tbeam[bt];
    if (remove_keys) {
      for (size_t r = 0; r < at.slots(); ++r) {
        int i = at.keys[r];
        if (i < 0 || !remove_keys->count(state_key(bt, i, j, -1, -1))) continue;
        at.erase(i);
        table[ix(i, j)] = NEG_INF;
      }
      at.compact();
    }
    if (beam <= 0 || at.live <= beam) {
      if (g_beam_stats.f && beam > 0)
        g_beam_stats.log_prune(BT_NAMES[bt], j, static_cast<size_t>(at.live),
                               static_cast<size_t>(at.live), 0.0, 0);
      return;
    }
    std::vector<std::pair<double, int>> ranked;
    ranked.reserve(at.live);
    // --stem-est-principled: Qp by alpha + log Ctx(i, j); -nested: Qb / Qm / Qm2 by alpha +
    // their nested outside
    const std::vector<double>* nb = nullptr;
    if (pr_nested && stem_)
      nb = bt == BT_QB ? &stem_->prQb_b : bt == BT_QM ? &stem_->prQm_b
                                         : bt == BT_QM2 ? &stem_->prQm2_b : nullptr;
    const bool pq = pr_on && bt == BT_QP;
    const auto tp = std::chrono::steady_clock::now();
    for (size_t r = 0; r < at.slots(); ++r) {
      int i = at.keys[r];
      if (i < 0) continue;
      if (pq) {
        const double s = at.vals[r] + pr_ctx(i, j);
        ranked.emplace_back(s < HALF_NEG ? NEG_INF : s, i);
        if (pr_check) {
          const long double c = prn_ctx(i, j);
          prc_cmp(5, ranked.back().first,
                  c < HALF_NEG ? NEG_INF : static_cast<double>(at.vals[r] + c));
        }
      } else if (nb) {
        const double s = at.vals[r] + (*nb)[ix(i, j)];
        ranked.emplace_back(s < HALF_NEG ? NEG_INF : s, i);
      } else if (pre_key_inside) {
        ranked.emplace_back(at.vals[r], i);
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
    if (g_beam_stats.f) {
      std::vector<double> keys;
      keys.reserve(ranked.size());
      for (const auto& e : ranked) keys.push_back(e.first);
      const size_t kept_n = static_cast<size_t>(beam);
      g_beam_stats.log_prune(BT_NAMES[bt], j, ranked.size(), kept_n,
                             BeamStats::mass_share(keys, kept_n), 1);
    }
    std::vector<char> keep(n, 0);
    for (int k = 0; k < beam; ++k) keep[ranked[k].second] = 1;
    if (drops) {
      ++drops->events[bt];
      for (size_t k = beam; k < ranked.size(); ++k) {
        const int i = ranked[k].second;
        note_drop(bt, bt, j, -1, i, j, -1, -1, k, ranked.size(), ranked[k].first);
        if (extra_at && extra_at->has(i))
          note_drop(BT_QM2, bt, j, -1, i, j, -1, -1, k, ranked.size(), ranked[k].first);
      }
    }
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

