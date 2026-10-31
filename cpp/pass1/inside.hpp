// Pass 1, LinearPartition-style nested inside at beam b1: mates, lead / multiloop helpers,
// lp_beam, lp_step, fill_Qb, Qm, Qmp (Q^b, Q^m, Q^m2).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  // ------------------------------------------------------------ two-index
  // Next index strictly after `after` that pairs with i, or -1. O(1) when the
  // base-type successor is a real pair; a restricted pair list may skip a few.
  void lp_build_mates() {
    if (mates_ready) return;
    mates_ready = true;
    for (int b = 0; b < 4; ++b) {
      mate_after[b].assign(static_cast<size_t>(n), -1);
      int nxt = -1;
      for (int j = n - 1; j >= 0; --j) {
        mate_after[b][static_cast<size_t>(j)] = nxt;
        const int c = E.code[j];
        const int t = c < 0 ? -1 : pair_type(b, c);
        if (t >= 0 && (allow_gu || (t != 3 && t != 5))) nxt = j;
      }
    }
  }
  int next_mate(int i, int after) const {
    if (i < 0 || after >= n - 1) return -1;
    const int b = E.code[i];
    if (b < 0) return -1;
    int j = mate_after[b][static_cast<size_t>(after)];
    while (j != -1 && !E.pairs(i, j)) j = j + 1 < n ? mate_after[b][static_cast<size_t>(j)] : -1;
    return j;
  }
  double unpaired_gap(int s, int t) const {
    if (t < s) return 0.0;
    return E.w_alpha3 * static_cast<double>(t - s + 1) + (dz ? E.run(s, t) : 0.0);
  }
  // Closing-pair left end for a first branch at `branch` (unpaired = branch - left - 1).
  int lead_i_lo(int branch) const {
    if (ml_lead_cap <= 0) return 0;
    return std::max(0, branch - ml_lead_cap);
  }
  // Qm-start form: unpaired bases from `start` up to but not including the branch.
  int lead_qm_lo(int branch) const {
    if (ml_lead_cap <= 0) return 0;
    return std::max(0, branch - (ml_lead_cap - 1));
  }
  bool lead_gap_ok(int gap) const {
    return ml_lead_cap <= 0 || gap <= ml_lead_cap - 1;
  }
  int lead_span() const { return ml_lead_cap <= 0 ? n + 1 : ml_lead_cap; }
  // Left unpaired run: the window of `cap` bases (x-1-cap).
  int ml_run_lo(int branch) const { return std::max(0, branch - 1 - cap); }
  void lp_beam(OMap& at) {
    const int B = tbeam[BT_QB];
    if (B <= 0 || at.live <= B) return;
    std::vector<std::pair<double, int>> ranked;
    ranked.reserve(static_cast<size_t>(at.live));
    for (size_t r = 0; r < at.slots(); ++r) {
      const int i = at.keys[r];
      if (i < 0) continue;
      const double s = prefix_log(i) + at.vals[r];
      ranked.emplace_back(s, i);
    }
    std::sort(ranked.begin(), ranked.end(), [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
      if (a.first != b.first) return a.first > b.first;
      return a.second > b.second;
    });
    std::vector<char> keep(static_cast<size_t>(n), 0);
    const int keep_n = std::min(B, static_cast<int>(ranked.size()));
    for (int k = 0; k < keep_n; ++k) keep[static_cast<size_t>(ranked[k].second)] = 1;
    for (size_t r = 0; r < at.slots(); ++r) {
      const int i = at.keys[r];
      if (i < 0 || keep[static_cast<size_t>(i)]) continue;
      at.erase(i);
    }
    at.compact();
  }
  // Right-run class. Length >= 2 shares one class: its extension delta is
  // d3[new] - d3[old], independent of the start. Empty and one base do not.
  int lp_kind(int gap) const {
    if (gap <= 0) return 0;
    if (!dz) return 1;
    return gap == 1 ? 1 : 2;
  }
  // One next-mate step of the right run. nk is the class at qn.
  bool lp_step(int p, int q, int kind, int& qn, int& nk, double& delta) const {
    qn = next_mate(p, q);
    if (qn < 0) return false;
    const int add = qn - q;
    if (add <= 0) return false;
    delta = E.w_alpha3 * static_cast<double>(add);
    if (!dz) {
      nk = 1;
      return true;
    }
    if (kind == 0) {
      nk = add == 1 ? 1 : 2;
      delta += E.run(q, qn - 1);
    } else if (kind == 1) {
      nk = 2;
      delta += E.run(q - 1, qn - 1) - E.s1w[q - 1];
    } else {
      nk = 2;
      delta += E.d3w[qn - 1] - E.d3w[q - 1];
    }
    return true;
  }
  void mu_add(int q, int p, int kind, double w) {
    if (q < 0 || q >= n || p < 0 || w < HALF_NEG) return;
    if (lpM.empty()) lpM.assign(static_cast<size_t>(n), {});
    OMap& col = lpM[static_cast<size_t>(q)];
    const int k = p * 3 + kind;
    const double cur = col.get(k);
    col.set(k, cur < HALF_NEG ? w : log_add(cur, w));
  }
  // Beam a Multi column. The key is p*3+kind; the prefix is C[p-1].
  void lp_beam_multi(OMap& at) {
    const int B = tbeam[BT_QB];
    if (B <= 0 || at.live <= B) return;
    struct R {
      double s;
      int k, p, kind;
    };
    std::vector<R> ranked;
    ranked.reserve(static_cast<size_t>(at.live));
    for (size_t r = 0; r < at.slots(); ++r) {
      const int k = at.keys[r];
      if (k < 0) continue;
      const int p = k / 3;
      const double s = prefix_log(p) + at.vals[r];
      ranked.push_back({s, k, p, k % 3});
    }
    std::sort(ranked.begin(), ranked.end(), [](const R& a, const R& b) {
      if (a.s != b.s) return a.s > b.s;
      if (a.p != b.p) return a.p > b.p;
      return a.kind > b.kind;
    });
    std::vector<int> drop;
    for (size_t i = static_cast<size_t>(std::min(B, static_cast<int>(ranked.size())));
         i < ranked.size(); ++i)
      drop.push_back(ranked[i].k);
    for (int k : drop) at.erase(k);
    at.compact();
  }
  // Prune H[j], seed H(j, first mate), extend only the kept hairpins, then free H[j].
  void lp_fill_qb(int j) {
    lp_build_mates();
    if (lpH.empty()) lpH.assign(static_cast<size_t>(n), {});
    if (lpM.empty()) lpM.assign(static_cast<size_t>(n), {});
    OMap& col = lpH[static_cast<size_t>(j)];
    lp_beam(col);
    const auto th = std::chrono::steady_clock::now();
    const int q0 = next_mate(j, j);
    if (q0 >= 0) {
      const double h0 = E.hairpin(j, q0);
      if (h0 > HALF_NEG) lpH[static_cast<size_t>(q0)].set(j, h0);
    }
    for (size_t r = 0; r < col.slots(); ++r) {
      const int i = col.keys[r];
      if (i < 0) continue;
      ++r_hair;
      set_pair(Qb, Qb_at, i, j, col.vals[r]);
      const int q = next_mate(i, j);
      if (q >= 0) {
        const double h = E.hairpin(i, q);
        if (h > HALF_NEG) lpH[static_cast<size_t>(q)].set(i, h);
      }
    }
    sec_hair += std::chrono::duration<double>(std::chrono::steady_clock::now() - th).count();
    col.init(0);
    const auto tc = std::chrono::steady_clock::now();
    lp_close_col(j);
    sec_mclose += std::chrono::duration<double>(std::chrono::steady_clock::now() - tc).count();
    // pair_bonus is one add on the finished pair sum (hairpin + interior + close),
    // the same place as set_pair(acc + pair_bonus) in the scan fill.
    if (pair_bonus != 0.0) {
      OMap& qb = Qb_at[j];
      for (size_t r = 0; r < qb.slots(); ++r) {
        const int i = qb.keys[r];
        if (i < 0) continue;
        const double v = qb.vals[r] + pair_bonus;
        qb.vals[r] = v;
        Qb[ix(i, j)] = v;
      }
    }
  }
  // Kept Multi(p, j) closes into Qb(p, j), then steps to the next mate.
  // The column that remains is the beam output; discarded arrivals are erased.
  void lp_close_col(int j) {
    if (lpM.empty()) return;
    OMap& col = lpM[static_cast<size_t>(j)];
    lp_beam_multi(col);
    const double base = E.w_alpha1 + E.w_alpha2;
    for (size_t r = 0; r < col.slots(); ++r) {
      const int k = col.keys[r];
      if (k < 0) continue;
      const int p = k / 3;
      const int kind = k % 3;
      const double w = col.vals[r];
      if (E.end_ok(p, j)) set_pair(Qb, Qb_at, p, j, w + base + E.at_close(p, j));
      int qn, nk;
      double delta;
      if (lp_step(p, j, kind, qn, nk, delta)) mu_add(qn, p, nk, w + delta);
    }
  }
  // From each kept pair, the outer pair sits inside the loop cap. Push it forward.
  // Same shape as LinearPartition handleP's internal/stack loop, with max_loop as the cap.
  void lp_push_interior(int j) {
    const auto ti = std::chrono::steady_clock::now();
    lp_build_mates();
    const OMap& at = Qb_at[j];
    for (size_t r = 0; r < at.slots(); ++r) {
      const int i = at.keys[r];
      if (i < 0) continue;
      const double qb = at.vals[r];
      if (qb < HALF_NEG) continue;
      const int p_lo = std::max(0, i - 1 - cap);
      for (int p = i - 1; p >= p_lo; --p) {
        const int l1 = i - p - 1;
        const int budget = cap - l1;
        if (budget < 0) break;
        for (int q = next_mate(p, j); q != -1 && q - j - 1 <= budget; q = next_mate(p, q)) {
          ++r_int;
          const double w = E.interior(p, i, j, q);
          if (w < HALF_NEG) continue;
          set_pair(Qb, Qb_at, p, q, w + qb);
        }
      }
    }
    sec_int += std::chrono::duration<double>(std::chrono::steady_clock::now() - ti).count();
  }
  // M2 that ends on a branch (lpP, no trailing run) emits Multi(p, next_mate(p, j)).
  // The right run j+1..q-1 is on that edge. The leading run is p+1..ib-1, capped
  // by ml_lead_cap. Arrivals at the same (p, q, kind) merge.
  void lp_emit_multi(int j) {
    if (lpP.empty() || j < 0 || j >= n) return;
    const auto te = std::chrono::steady_clock::now();
    lp_build_mates();
    if (lpM.empty()) lpM.assign(static_cast<size_t>(n), {});
    OMap& src = lpP[static_cast<size_t>(j)];
    std::vector<int> drop;
    for (size_t r = 0; r < src.slots(); ++r) {
      const int i = src.keys[r];
      if (i < 0) continue;
      if (!Qm2_at[j].has(i)) drop.push_back(i);
    }
    for (int i : drop) src.erase(i);
    src.compact();
    const int L = lead_span();
    for (size_t r = 0; r < src.slots(); ++r) {
      const int ib = src.keys[r];
      if (ib < 0) continue;
      const double mv = src.vals[r];
      if (mv < HALF_NEG) continue;
      const int p_lo = std::max(0, ib - L);
      for (int p = ib - 1; p >= p_lo; --p) {
        const int q = next_mate(p, j);
        if (q < 0) continue;
        ++r_m2multi;
        const int kind = lp_kind(q - j - 1);
        const double w = mv + unpaired_gap(p + 1, ib - 1) + unpaired_gap(j + 1, q - 1);
        mu_add(q, p, kind, w);
      }
    }
    sec_emit += std::chrono::duration<double>(std::chrono::steady_clock::now() - te).count();
  }

  void lp_mb_add(int q, int k, double w) {
    if (q < 0 || q >= n || w < HALF_NEG) return;
    if (lpMb.empty()) lpMb.assign(static_cast<size_t>(n), {});
    OMap& col = lpMb[static_cast<size_t>(q)];
    const double cur = col.get(k);
    col.set(k, cur < HALF_NEG ? w : log_add(cur, w));
  }
  // Outside of Multi(p, j): the closing pair, plus the jump to the next mate.
  void lp_out_multi(int j) {
    if (l2r_mode != 0 || lpM.empty() || j < 0 || j >= static_cast<int>(lpM.size())) return;
    if (lpMb.empty()) lpMb.assign(static_cast<size_t>(n), {});
    const OMap& col = lpM[static_cast<size_t>(j)];
    for (size_t r = 0; r < col.slots(); ++r) {
      const int k = col.keys[r];
      if (k < 0) continue;
      const int p = k / 3;
      const int kind = k % 3;
      if (E.end_ok(p, j)) {
        const double bb = Qb_b[ix(p, j)];
        if (bb > HALF_NEG)
          lp_mb_add(j, k, bb + E.w_alpha1 + E.w_alpha2 + E.at_close(p, j));
      }
      int qn, nk;
      double delta;
      if (!lp_step(p, j, kind, qn, nk, delta)) continue;
      const int kn = p * 3 + nk;
      if (!lpM[static_cast<size_t>(qn)].has(kn)) continue;
      const double cb = lpMb[static_cast<size_t>(qn)].get(kn);
      if (cb > HALF_NEG) lp_mb_add(j, k, cb + delta);
    }
  }
  // E3: Multi(p, q) -> M2(ib, j), q = next_mate(p, j). Qm2_b is already the
  // branch-ending (P) outside, so this does not flow into the trailing run.
  void lp_out_e3(int j) {
    if (l2r_mode != 0 || lpP.empty() || lpM.empty() || lpMb.empty()) return;
    if (j < 0 || j >= n) return;
    const int L = lead_span();
    const OMap& src = lpP[static_cast<size_t>(j)];
    for (size_t r = 0; r < src.slots(); ++r) {
      const int ib = src.keys[r];
      if (ib < 0 || !Qm2_at[j].has(ib)) continue;
      const int p_lo = std::max(0, ib - L);
      for (int p = ib - 1; p >= p_lo; --p) {
        const int q = next_mate(p, j);
        if (q < 0) continue;
        const int k = p * 3 + lp_kind(q - j - 1);
        if (q >= n || !lpM[static_cast<size_t>(q)].has(k)) continue;
        const double mb = lpMb[static_cast<size_t>(q)].get(k);
        if (mb < HALF_NEG) continue;
        const double factor = unpaired_gap(p + 1, ib - 1) + unpaired_gap(j + 1, q - 1);
        add2b(Qm2_b, ib, j, mb + factor);
      }
    }
  }

  void fill_Qb(int j) {
    // Nested pass only. --pk-l2r off is Beam-DP03: hairpin, interior, Qm2, and PKM.
    if (l2r_mode == 0 && !allow_pk) {
      lp_fill_qb(j);
      return;
    }
    // Interior spans are bounded by the loop cap, so the outer left end sits in a
    // window of width cap next to a kept inner pair. Walk those kept pairs and the
    // partners inside the window; do not walk every left partner of j.
    if (static_cast<int>(qb_acc.size()) < n) qb_acc.assign(static_cast<size_t>(n), NEG_INF);
    const std::vector<int>& pl = partners_left[j];
    for (int i : pl) {
      ++r_hair;
      qb_acc[static_cast<size_t>(i)] = E.hairpin(i, j);
    }
    const int e0 = std::max(0, j - 1 - cap);
    for (int e = e0; e < j; ++e) {
      const int budget = cap - (j - e - 1);
      const OMap& at = Qb_at[e];
      for (size_t r = 0; r < at.slots(); ++r) {
        const int d = at.keys[r];
        if (d < 0) continue;
        const double qb = at.vals[r];
        const int i_hi = std::min(d - 1, e - MIN_PAIR);
        const int i_lo = std::max(0, d - 1 - budget);
        if (i_lo > i_hi) continue;
        auto it = std::lower_bound(pl.begin(), pl.end(), i_lo);
        const auto iend = std::upper_bound(it, pl.end(), i_hi);
        for (; it != iend; ++it) {
          ++r_int;
          const int i = *it;
          const double w = E.interior(i, d, e, j);
          if (w < HALF_NEG) continue;
          double& acc = qb_acc[static_cast<size_t>(i)];
          acc = log_add(acc, w + qb);
        }
      }
    }
    for (int i : pl) {
      double acc = qb_acc[static_cast<size_t>(i)];
      qb_acc[static_cast<size_t>(i)] = NEG_INF;

      const bool closing_ok = E.end_ok(i, j);
      const double at_ij = E.at_close(i, j);

      if (closing_ok && i + 1 <= j - 1) {
        ++r_m2multi;
        double m2 = Qm2[ix(i + 1, j - 1)];
        if (m2 > HALF_NEG) acc = log_add(acc, m2 + E.w_alpha1 + E.w_alpha2 + at_ij);
      }

      if (allow_pk && closing_ok && j >= 1) {
        const OMap& pkm = PKM_at[j - 1];
        for (size_t r = 0; r < pkm.slots(); ++r) {
          int d = pkm.keys[r];
          if (d <= i) continue;
          double pk = pkm.vals[r];
          double w_pk = pk + E.w_alpha1 + E.w_beta1m + (1.0 + E.ml_pk_k) * E.w_alpha2 +
                        E.w_alpha3 * (d - i - 1) + at_ij;
          if (dz) w_pk += E.run(i + 1, d - 1);
          acc = log_add(acc, w_pk);
          if (i + 1 <= d - 1) {
            double left = Qm[ix(i + 1, d - 1)];
            if (left > HALF_NEG)
              acc = log_add(acc, left + pk + E.w_alpha1 + E.w_beta1m +
                                     (1.0 + E.ml_pk_k) * E.w_alpha2 + at_ij);
          }
        }
      }
      set_pair(Qb, Qb_at, i, j, acc + pair_bonus);
    }
  }

  void build_qlead(int t) {
    if (qlead.empty()) qlead.assign(static_cast<size_t>(n), {});
    const double w3 = E.w_alpha3;
    std::vector<LeadEnt> cand;
    if (t == 0) {
      cand.push_back({0, 0.0});
    } else {
      cand.reserve(qlead[static_cast<size_t>(t - 1)].size() + 1);
      for (const LeadEnt& e : qlead[static_cast<size_t>(t - 1)]) {
        if (t - e.x0 > E.ml_left || !lead_gap_ok(t - e.x0)) continue;
        double w = e.w + w3;
        if (dz) w += E.run(e.x0, t - 1) - E.run(e.x0, t - 2);
        cand.push_back({e.x0, w});
      }
      cand.push_back({t, 0.0});
    }
    const int B = p2_qz_beam;
    if (B > 0 && static_cast<int>(cand.size()) > B) {
      auto key = [&](const LeadEnt& e) {
        return e.w + (e.x0 > 0 ? C[static_cast<size_t>(e.x0 - 1)] : 0.0);
      };
      std::nth_element(cand.begin(), cand.begin() + B, cand.end(), [&](const LeadEnt& a, const LeadEnt& b) {
        const double ka = key(a), kb = key(b);
        if (ka != kb) return ka > kb;
        return a.x0 < b.x0;
      });
      cand.resize(static_cast<size_t>(B));
    }
    std::sort(cand.begin(), cand.end(), [](const LeadEnt& a, const LeadEnt& b) { return a.x0 < b.x0; });
    qlead[static_cast<size_t>(t)] = std::move(cand);
  }

  void add_branch(int j, int d, double obj, bool pk, double w_at) {
    const double w2 = E.w_alpha2;
    double w_single = obj + ((pk ? E.w_beta1m + E.ml_pk_k * w2 : w2) + w_at);
    int lead_lim = E.ml_left;
    if (ml_lead_cap > 0 && ml_lead_cap - 1 < lead_lim) lead_lim = ml_lead_cap - 1;
    const int lo = std::max(0, d - lead_lim);
    if (l2r_mode == 0 && !allow_pk) {
      // LinearPartition M = P: the fragment starts at the pair. Unpaired bases sit to its right.
      const auto tm = std::chrono::steady_clock::now();
      ++r_mp;
      set_Qm(d, j, w_single, false);
      sec_mp += std::chrono::duration<double>(std::chrono::steady_clock::now() - tm).count();
    } else if (static_cast<int>(qlead.size()) == n) {
      for (const LeadEnt& e : qlead[static_cast<size_t>(d)]) {
        if (e.x0 < lo || e.x0 > d) continue;
        ++r_mp;
        set_Qm(e.x0, j, w_single + e.w, false);
      }
    } else {
      const double w3 = E.w_alpha3;
      for (int i = lo; i <= d; ++i) {
        double w = w_single + w3 * (d - i);
        if (dz && i < d) w += E.run(i, d - 1);
        ++r_mp;
        set_Qm(i, j, w, false);
      }
    }
    if (d == 0) return;
    double w_more = obj + ((pk ? E.w_beta1m + E.ml_pk_k * w2 : w2) + w_at);
    const OMap& at = Qm_at[d - 1];
    const auto tx = l2r_mode == 0 ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
    for (size_t r = 0; r < at.slots(); ++r) {
      int i = at.keys[r];
      if (i < 0) continue;
      ++r_m2xp;
      const double nv = at.vals[r] + w_more;
      set_Qm(i, j, nv, !pk);
      if (!pk && nv > HALF_NEG) {
        if (lpP.empty()) lpP.assign(static_cast<size_t>(n), {});
        OMap& pc = lpP[static_cast<size_t>(j)];
        const double cur = pc.get(i);
        pc.set(i, cur < HALF_NEG ? nv : log_add(cur, nv));
      }
    }
    if (l2r_mode == 0)
      sec_m2xp += std::chrono::duration<double>(std::chrono::steady_clock::now() - tx).count();
  }

  void fill_Qm(int j) {
    const double w3 = E.w_alpha3;
    std::vector<int> surv, surv2;
    const auto t_ext = std::chrono::steady_clock::now();
    if (dz) {
      RunSt& cu = rsQm[j & 1];
      RunSt& cu2 = rsQm2[j & 1];
      cu.clear();
      cu2.clear();
      if (j > 0) {
        surv = live_keys(Qm_at[j - 1]);
        surv2 = live_keys(Qm2_at[j - 1]);
        r_mmu += static_cast<long long>(surv.size());
        r_muext += static_cast<long long>(surv2.size());
        run_ext(rsQm[(j - 1) & 1], cu, surv, j, w3);
        run_ext(rsQm2[(j - 1) & 1], cu2, surv2, j, w3);
      }
    } else if (j > 0) {
      const OMap& at = Qm_at[j - 1];
      for (size_t r = 0; r < at.slots(); ++r) {
        int i = at.keys[r];
        if (i < 0) continue;
        ++r_mmu;
        set_Qm(i, j, at.vals[r] + w3, false);
      }
      const OMap& at2 = Qm2_at[j - 1];
      for (size_t r = 0; r < at2.slots(); ++r) {
        int i = at2.keys[r];
        if (i < 0) continue;
        ++r_muext;
        double nv = at2.vals[r] + w3;
        double t = log_add(static_cast<double>(Qm2[ix(i, j)]), nv);
        Qm2[ix(i, j)] = t;
        Qm2_at[j].set(i, t);
      }
    }
    sec_qm += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_ext).count();
    const auto t_br = std::chrono::steady_clock::now();
    {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        int d = at.keys[r];
        if (d < 0 || !E.end_ok(d, j)) continue;
        add_branch(j, d, at.vals[r], false, E.at(d, j));
      }
    }
    if (allow_pk) {
      const OMap& at = Qp_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        int d = at.keys[r];
        if (d < 0) continue;
        add_branch(j, d, at.vals[r], true, 0.0);
      }
    }
    sec_branch += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_br).count();
    if (dz) {
      const auto t_fin = std::chrono::steady_clock::now();
      run_fin(&Qm, Qm_at[j], rsQm[j & 1], surv, j);
      run_fin(&Qm2, Qm2_at[j], rsQm2[j & 1], surv2, j);
      sec_qm += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_fin).count();
    }
  }

  void set_Qmp(int i, int j, double val) {
    if (val < HALF_NEG) return;
    if (dz) {
      rsQmp[j & 1].P.add(i, val);
      return;
    }
    double t = log_add(static_cast<double>(Qmp[ix(i, j)]), val);
    Qmp[ix(i, j)] = t;
    Qmp_at[j].set(i, t);
  }

  void add_branch_p(int j, int d, double obj, bool pk, double w_at) {
    const double w2 = E.w_bp;
    const double w3 = E.w_cp;
    double w_single = obj + ((pk ? E.w_beta1m + E.ml_pk_k * w2 : w2) + w_at);
    const int i_lo = lead_qm_lo(d);
    if (dz) {
      for (int i = i_lo; i <= d; ++i) set_Qmp(i, j, w_single + w3 * (d - i) + E.run(i, d - 1));
    } else {
      for (int i = i_lo; i <= d; ++i) set_Qmp(i, j, w_single + w3 * (d - i));
    }
    if (d == 0) return;
    const OMap& at = Qmp_at[d - 1];
    for (size_t r = 0; r < at.slots(); ++r) {
      int i = at.keys[r];
      if (i < 0) continue;
      set_Qmp(i, j, at.vals[r] + w_single);
    }
  }

  // fill_Qm with the stem-multiloop penalties b_p, c_p (no M2 split needed).
  void fill_Qmp(int j) {
    const double w3 = E.w_cp;
    std::vector<int> surv;
    if (dz) {
      RunSt& cu = rsQmp[j & 1];
      cu.clear();
      if (j > 0) {
        surv = live_keys(Qmp_at[j - 1]);
        run_ext(rsQmp[(j - 1) & 1], cu, surv, j, w3);
      }
    } else if (j > 0) {
      const OMap& at = Qmp_at[j - 1];
      for (size_t r = 0; r < at.slots(); ++r) {
        int i = at.keys[r];
        if (i < 0) continue;
        set_Qmp(i, j, at.vals[r] + w3);
      }
    }
    {
      const OMap& at = Qb_at[j];
      for (size_t r = 0; r < at.slots(); ++r) {
        int d = at.keys[r];
        if (d < 0 || !E.end_ok(d, j)) continue;
        add_branch_p(j, d, at.vals[r], false, E.at(d, j));
      }
    }
    const OMap& at = Qp_at[j];
    for (size_t r = 0; r < at.slots(); ++r) {
      int d = at.keys[r];
      if (d < 0) continue;
      add_branch_p(j, d, at.vals[r], true, 0.0);
    }
    if (dz) run_fin(&Qmp, Qmp_at[j], rsQmp[j & 1], surv, j);
  }

