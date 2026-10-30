// Pass 2, arms: gap coverage, k_bound_check, Sbar_pk (pr_sf_build), arm_raw / arm_phi
// (K^pre_L(q'; b, c) = sum_{f>q'} Q^z_pre(q'+1, f-1) V), persistent argmax index, arm
// queries (pr_armq), arm selection and the separable bound (pr_sep_*).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  void pr_gap_begin(std::vector<std::pair<int, int>>& iv) {
    if (!pr_arm_gaponly) return;
    std::sort(iv.begin(), iv.end());
    prGapIv.clear();
    for (const auto& [l, r] : iv) {
      if (l > r) continue;
      if (!prGapIv.empty() && l <= prGapIv.back().second + 1)
        prGapIv.back().second = std::max(prGapIv.back().second, r);
      else
        prGapIv.emplace_back(l, r);
    }
    ++gap_epoch;
    ++gap_epochs;
  }
  bool pr_gap_covered(int l, int r) const {
    auto it = std::upper_bound(prGapIv.begin(), prGapIv.end(), std::make_pair(l, INT32_MAX));
    return it != prGapIv.begin() && std::prev(it)->second >= r;
  }

  // phi(arm r) at q' as Step A computes it (separable: from the suffix sums prSF; with
  // --pr-partner-span W only the entries q' < y' <= q' + W, summed in y' order).
  // --v-sweep: lK reads the swept b_V entries. mu is the profile weight; the arm
  // base (pair energy, stem weight, pair_bonus) stays the one prYw was built with.
  // A partner the sweep did not keep is -inf. The floor is not applied here.
  void pr_yw_from_sweep() {
    if (!p2_vsweep || !stem_ || stem_->prMu.empty()) return;
    long long rep = 0, miss = 0, ent = 0;
    for (int r = 0; r < prA; ++r) {
      const int xs = prArmXs[r], xp = prArmXp[r];
      for (int k = prArmY0[r]; k < prArmY0[r + 1]; ++k) {
        ++ent;
        const int yp = prYy[k];
        double mu = NEG_INF;
        auto it = stem_->prMu.find(key2(xp, yp));
        if (it != stem_->prMu.end()) {
          auto jt = it->second.find(key2(xs, 0));
          if (jt != it->second.end()) mu = jt->second;
        }
        if (mu < HALF_NEG && p2_v_cut_fb &&
            p2_v_cap > 0 && it != stem_->prMu.end() &&
            static_cast<int>(it->second.size()) >= p2_v_cap) {
          mu = v_cut_fb_mu(xp, yp, xs);
          if (mu > HALF_NEG) ++vfb_n;
        }
        if (mu < HALF_NEG) {
          prYw[k] = NEG_INF;
          ++miss;
          continue;
        }
        prYw[k] = (mu + (E.at(xp, yp) + E.w_stem)) + pair_bonus;
        ++rep;
      }
    }
    prSF.clear();
    std::fprintf(stderr, "[v-sweep-lk] arms %d entries %lld replaced %lld missing %lld\n", prA, ent,
                 rep, miss);
  }
  // lK(arm, q') >= lKx(arm, q') on the weights lK is about to use. One pass, not on the hot path.
  void k_bound_check() {
    pr_sf_build();
    long long checks = 0, bad = 0, bad_inf = 0, bad_fin = 0;
    double worst = 0.0, wub = 0.0, wex = 0.0;
    int wr = -1, wq = -1;
    for (int r = 0; r < prA; ++r) {
      for (int q = 0; q + 1 <= n - 1; ++q) {
        const double ub = lK(r, q);
        const double ex = lKx_sum(r, q);
        ++checks;
        const bool ub_fin = ub > HALF_NEG, ex_fin = ex > HALF_NEG;
        double gap = 0.0;
        if (ex_fin && !ub_fin) gap = 1e300;
        else if (ex_fin && ub_fin) gap = ex - ub;
        if (gap > 1e-8) {
          ++bad;
          if (!ub_fin) ++bad_inf;
          else ++bad_fin;
          if (gap > worst) {
            worst = gap;
            wub = ub;
            wex = ex;
            wr = r;
            wq = q;
          }
        }
      }
    }
    std::fprintf(stderr,
                 "[k-bound] checks %lld violations %lld inf %lld finite %lld worst %.6e ub %.6e ex %.6e arm %d q %d\n",
                 checks, bad, bad_inf, bad_fin, worst, wub, wex, wr, wq);
  }
  void pr_sf_build() {
    if (!pr_aq_sep || !prSF.empty()) return;
    const StemEst& st = *stem_;
    prSF.assign(prYy.size(), NEG_INF);
    for (int r = 0; r < prA; ++r) {
      double acc = NEG_INF;
      for (int k = prArmY0[r + 1] - 1; k >= prArmY0[r]; --k) {
        const double z = st.Qzbar(0, prYy[k] - 1);
        if (z > HALF_NEG && prYw[k] > HALF_NEG) acc = log_add(acc, z + prYw[k]);
        prSF[k] = acc;
      }
    }
  }
  // raw: phi + log Qzbar[0][q'] (separable) or phi (not); ranking only depends on raw.
  double arm_raw(int r, int qp) const {
    const StemEst& st = *stem_;
    const int hi = prArmY0[r + 1];
    int k = static_cast<int>(std::upper_bound(prYy.begin() + prArmY0[r], prYy.begin() + hi, qp) -
                             prYy.begin());
    if (pr_pspan > 0 || !pr_aq_sep) {
      double acc = NEG_INF;
      for (; k < hi && (pr_pspan <= 0 || prYy[k] - qp <= pr_pspan); ++k) {
        const double z = pr_aq_sep ? st.Qzbar(0, prYy[k] - 1) : st.Qzbar(qp + 1, prYy[k] - 1);
        if (z > HALF_NEG && prYw[k] > HALF_NEG) acc = log_add(acc, z + prYw[k]);
      }
      return acc;
    }
    return k < hi ? prSF[k] : NEG_INF;
  }
  double arm_phi(int r, int qp, double z0) const {
    const double raw = arm_raw(r, qp);
    if (!pr_aq_sep) return raw;
    return raw > HALF_NEG && z0 > HALF_NEG ? raw - z0 : NEG_INF;
  }

  // Persistent argmax trees over the arm positions of order 1 (arm id) and order 2 (prByXp),
  // one version per q' (pRoot[t][q']): sweeping q' down from n-1, an arm's leaf is set to its
  // raw value whenever an entry y' = q' + 1 enters (or, with W, y' = q' + W + 1 leaves) its
  // window.  Ties: lower position, as aq_better.  Node 0: the empty tree.
  struct PNode {
    int l, r, best, ver;
    double val;
  };
  std::vector<PNode> pnd;
  std::vector<int> pRoot[2], prXpPos;
  int pver = 0;
  bool pidx_built = false;
  void pidx_pull(int nd) {
    const PNode& L = pnd[pnd[nd].l];
    const PNode& R = pnd[pnd[nd].r];
    int b;
    double v;
    if (L.best < 0 || (R.best >= 0 && R.val > L.val)) {
      b = R.best;
      v = R.val;
    } else {
      b = L.best;
      v = L.val;
    }
    pnd[nd].best = b;
    pnd[nd].val = v;
  }
  int pidx_set(int nd, int lo, int hi, int pos, double val) {
    if (pnd[nd].ver != pver) {
      PNode c = pnd[nd];
      c.ver = pver;
      pnd.push_back(c);
      nd = static_cast<int>(pnd.size()) - 1;
    }
    if (lo == hi) {
      pnd[nd].best = val > HALF_NEG ? pos : -1;
      pnd[nd].val = val > HALF_NEG ? val : NEG_INF;
      return nd;
    }
    const int mid = (lo + hi) / 2;
    if (pos <= mid) {
      const int c = pidx_set(pnd[nd].l, lo, mid, pos, val);
      pnd[nd].l = c;
    } else {
      const int c = pidx_set(pnd[nd].r, mid + 1, hi, pos, val);
      pnd[nd].r = c;
    }
    pidx_pull(nd);
    return nd;
  }
  void pidx_q(int nd, int lo, int hi, int l, int r, int& bp, double& bv) const {
    if (nd == 0 || r < lo || hi < l) return;
    const PNode& x = pnd[nd];
    if (x.best < 0) return;
    if (l <= lo && hi <= r) {
      if (bp < 0 || x.val > bv || (x.val == bv && x.best < bp)) {
        bp = x.best;
        bv = x.val;
      }
      return;
    }
    const int mid = (lo + hi) / 2;
    pidx_q(x.l, lo, mid, l, r, bp, bv);
    pidx_q(x.r, mid + 1, hi, l, r, bp, bv);
  }
  void pidx_build() {
    if (pidx_built) return;
    pidx_built = true;
    SecTimer ts(pidx_sec);
    pr_sf_build();
    pnd.assign(1, PNode{0, 0, -1, -1, NEG_INF});
    for (auto& v : pRoot) v.assign(n, 0);
    prXpPos.assign(prA, 0);
    for (int pos = 0; pos < prA; ++pos) prXpPos[prByXp[pos]] = pos;
    std::vector<std::vector<int>> ev(n);
    for (int r = 0; r < prA; ++r)
      for (int k = prArmY0[r]; k < prArmY0[r + 1]; ++k) {
        if (prYy[k] - 1 >= 0) ev[prYy[k] - 1].push_back(r);
        if (pr_pspan > 0 && prYy[k] - pr_pspan - 1 >= 0) ev[prYy[k] - pr_pspan - 1].push_back(r);
      }
    std::vector<int> stamp(prA, -1);
    int root[2] = {0, 0};
    for (int qp = n - 1; qp >= 0; --qp) {
      ++pver;
      for (int r : ev[qp]) {
        if (stamp[r] == qp) continue;
        stamp[r] = qp;
        const double raw = arm_raw(r, qp);
        root[0] = pidx_set(root[0], 0, prA - 1, r, raw);
        root[1] = pidx_set(root[1], 0, prA - 1, prXpPos[r], raw);
        ++pidx_upd;
      }
      std::vector<int>().swap(ev[qp]);
      pRoot[0][qp] = root[0];
      pRoot[1][qp] = root[1];
    }
    pidx_nodes = static_cast<long long>(pnd.size());
    pidx_mb = (static_cast<double>(pnd.capacity()) * sizeof(PNode) +
               static_cast<double>(2 * n + prA) * sizeof(int)) / 1048576.0;
  }

  double aq_val(const ArmQ& a, int ord, int pos) const {
    return (pr_ratio ? a.rs : a.phi)[ord == 1 ? pos : prByXp[pos]];
  }
  int aq_better(const ArmQ& a, int ord, int p1, int p2) const {
    if (p1 < 0) return p2;
    if (p2 < 0) return p1;
    const double v1 = aq_val(a, ord, p1), v2 = aq_val(a, ord, p2);
    if (v1 != v2) return v1 > v2 ? p1 : p2;
    return std::min(p1, p2);
  }
  void aq_tree(ArmQ& a, int ord, std::vector<int>& t) const {
    t.assign(2 * static_cast<size_t>(a.sz), -1);
    for (int p = 0; p < prA; ++p) t[a.sz + p] = p;
    for (int k = a.sz - 1; k >= 1; --k) t[k] = aq_better(a, ord, t[2 * k], t[2 * k + 1]);
  }
  int aq_query(const ArmQ& a, int ord, int l, int r) const {
    const std::vector<int>& t = ord == 1 ? a.t1 : a.t2;
    int best = -1;
    for (int lo = l + a.sz, hi = r + a.sz + 1; lo < hi; lo >>= 1, hi >>= 1) {
      if (lo & 1) best = aq_better(a, ord, best, t[lo++]);
      if (hi & 1) best = aq_better(a, ord, best, t[--hi]);
    }
    return best;
  }

  ArmQ& pr_armq(int qp) {
    auto it = prAQ.find(qp);
    if (it != prAQ.end()) {
      if (!pr_arm_gaponly || it->second.epoch == gap_epoch) {
        it->second.use = ++pr_use;
        return it->second;
      }
      prAQ.erase(it);
    }
    if (prAQ.size() >= pr_aq_max) {
      auto old = prAQ.begin();
      for (auto jt = prAQ.begin(); jt != prAQ.end(); ++jt)
        if (jt->second.use < old->second.use) old = jt;
      prAQ.erase(old);
    }
    const auto t0 = std::chrono::steady_clock::now();
    const StemEst& st = *stem_;
    ArmQ a;
    if (pr_arm_persist) {
      pidx_build();
      a.use = ++pr_use;
      ++pr_aq_built;
      pr_aq_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      return prAQ.emplace(qp, std::move(a)).first->second;
    }
    pr_sf_build();
    const double z0 = pr_aq_sep ? st.Qzbar(0, qp) : 0.0;
    if (pr_arm_gaponly) {
      a.epoch = gap_epoch;
      for (const auto& [l, r] : prGapIv)
        for (int pos = prXpFirst[l]; pos < prXpFirst[std::min(r, n - 1) + 1]; ++pos) {
          a.g2.pos.push_back(pos);
          a.g1.pos.push_back(prByXp[pos]);
        }
      std::sort(a.g1.pos.begin(), a.g1.pos.end());
      for (int p : a.g1.pos) a.g1.val.push_back(arm_phi(p, qp, z0));
      for (int p : a.g2.pos) a.g2.val.push_back(arm_phi(prByXp[p], qp, z0));
      a.g1.build();
      a.g2.build();
      aq_arms += static_cast<long long>(a.g2.pos.size());
      a.use = ++pr_use;
      ++pr_aq_built;
      pr_aq_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      return prAQ.emplace(qp, std::move(a)).first->second;
    }
    a.phi.assign(prA, NEG_INF);
    aq_arms += prA;
    for (int r = 0; r < prA; ++r) {
      if (pr_pspan > 0) {
        a.phi[r] = arm_phi(r, qp, z0);
        continue;
      }
      const int hi = prArmY0[r + 1];
      int k = static_cast<int>(std::upper_bound(prYy.begin() + prArmY0[r], prYy.begin() + hi, qp) -
                               prYy.begin());
      if (pr_aq_sep) {
        if (k < hi && prSF[k] > HALF_NEG && z0 > HALF_NEG) a.phi[r] = prSF[k] - z0;
        continue;
      }
      double acc = NEG_INF;
      for (; k < hi; ++k) {
        const double z = st.Qzbar(qp + 1, prYy[k] - 1);
        if (z > HALF_NEG) acc = log_add(acc, z + prYw[k]);
      }
      a.phi[r] = acc;
    }
    if (pr_arms_mode == 1)
      for (int r = 0; r < prA; ++r) a.nfin += a.phi[r] > HALF_NEG;
    if (pr_ratio) {
      a.rs.resize(prA);
      for (int r = 0; r < prA; ++r)
        a.rs[r] = a.phi[r] > HALF_NEG ? a.phi[r] - E.wu * (prArmXp[r] - prArmXs[r] + 1) : NEG_INF;
    }
    if (pr_arms_k > 0 || pr_ratio) {
      a.sz = 1;
      while (a.sz < std::max(1, prA)) a.sz <<= 1;
      aq_tree(a, 1, a.t1);
      aq_tree(a, 2, a.t2);
    }
    a.use = ++pr_use;
    ++pr_aq_built;
    pr_aq_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return prAQ.emplace(qp, std::move(a)).first->second;
  }

  // Arms of gap (u, v) at q' in selection order: xs >= u+1 (two: xs = u+1), x' <= v-1, finite
  // phi; the top pr_arms_k by phi (ties: lower position), all if 0; ratio: by rs until the
  // adaptive stop (set_stem_pr_select).
  void pr_select(const ArmQ& a, int qp, int u, int v, bool two, std::vector<int>& out) {
    out.clear();
    if (u + 1 > v - 1) return;
    if (pr_arm_gaponly && !pr_gap_covered(u + 1, v - 1)) ++gap_miss;
    int lo, hi, ord;
    if (two) {
      lo = prXsFirst[u + 1];
      const int end = prXsFirst[u + 2];
      hi = static_cast<int>(std::upper_bound(prArmXp.begin() + lo, prArmXp.begin() + end, v - 1) -
                            prArmXp.begin());
      ord = 1;
    } else {
      lo = prXpFirst[u + 1];
      hi = prXpFirst[v];
      ord = 2;
    }
    auto arm_of = [&](int pos) { return ord == 1 ? pos : prByXp[pos]; };
    auto ok = [&](int r) { return two || prArmXs[r] >= u + 1; };
    if (pr_arms_k == 0 && !pr_ratio) {
      for (int pos = lo; pos < hi; ++pos) {
        const int r = arm_of(pos);
        if (ok(r) && a.phi[r] > HALF_NEG) out.push_back(r);
      }
      return;
    }
    // ratio: the gap's total r over its arms (scaled by the largest), for the adaptive stop
    double smax = NEG_INF, tot = 0.0, sel = 0.0;
    if (pr_ratio) {
      for (int pos = lo; pos < hi; ++pos) {
        const int r = arm_of(pos);
        if (ok(r)) smax = std::max(smax, a.rs[r]);
      }
      if (smax < HALF_NEG) return;
      for (int pos = lo; pos < hi; ++pos) {
        const int r = arm_of(pos);
        if (ok(r) && a.rs[r] > HALF_NEG) tot += std::exp(a.rs[r] - smax);
      }
    }
    const int n_max = pr_ratio ? pr_arms_max : pr_arms_k;
    struct HE {
      double v;
      int pos, l, r;
      bool operator<(const HE& o) const { return v != o.v ? v < o.v : pos > o.pos; }
    };
    std::priority_queue<HE> pq;
    double z0 = 0.0;
    int root = 0;
    if (pr_arm_persist) {
      z0 = stem_->Qzbar(0, qp);
      if (z0 < HALF_NEG) return;
      root = pRoot[ord - 1][qp];
    }
    auto push = [&](int l, int r) {
      if (l > r) return;
      int p;
      double val;
      if (pr_arm_persist) {
        p = -1;
        double raw = NEG_INF;
        pidx_q(root, 0, prA - 1, l, r, p, raw);
        if (p < 0) return;
        val = raw - z0;
      } else if (pr_arm_gaponly) {
        const CTree& c = ord == 1 ? a.g1 : a.g2;
        const int i = c.query(l, r);
        if (i < 0) return;
        p = c.pos[i];
        val = c.val[i];
      } else {
        p = aq_query(a, ord, l, r);
        if (p < 0) return;
        val = aq_val(a, ord, p);
      }
      if (val < HALF_NEG) return;
      pq.push({val, p, l, r});
    };
    push(lo, hi - 1);
    while (!pq.empty() && static_cast<int>(out.size()) < n_max) {
      const HE t = pq.top();
      pq.pop();
      const int r = arm_of(t.pos);
      if (ok(r)) {
        out.push_back(r);
        if (pr_ratio) {
          sel += std::exp(t.v - smax);
          if (tot - sel < pr_arms_eps * sel) break;
        }
      } else {
        ++pr_T_skip;
      }
      push(t.l, t.pos - 1);
      push(t.pos + 1, t.r);
    }
  }

  // sep: columns known+1..b of a's tables (arms with x' = c added to every D[c][x], x <= xs).
  void pr_sep_extend(ArmQ& a, int b) {
    const auto t0 = std::chrono::steady_clock::now();
    if (a.ra.empty()) {
      a.ra.assign(prA, NEG_INF);
      a.P1.assign(prA, NEG_INF);
    }
    std::vector<double> col;
    for (int c = a.known + 1; c <= b; ++c) {
      col.assign(c + 1, NEG_INF);
      for (int pos = prXpFirst[c]; pos < prXpFirst[c + 1]; ++pos) {
        const int r = prByXp[pos], xs = prArmXs[r];
        if (a.phi[r] > HALF_NEG) {
          const double z = pr_arm_z == 1 ? E.wu * (c - xs + 1) : pr_qz(xs, c);
          if (z > HALF_NEG) a.ra[r] = a.phi[r] - z;
        }
        col[xs] = log_add(col[xs], a.ra[r]);
        const bool cont = r > 0 && prArmXs[r - 1] == xs;
        a.P1[r] = log_add(cont ? a.P1[r - 1] : NEG_INF, a.ra[r]);
      }
      std::vector<double> d(c + 1, NEG_INF);
      double suf = NEG_INF;
      for (int x = c; x >= 0; --x) {
        suf = log_add(suf, col[x]);
        d[x] = log_add(x <= c - 1 ? a.D[c - 1][x] : NEG_INF, suf);
      }
      a.D.push_back(std::move(d));
      ++pr_sep_cols;
    }
    a.known = std::max(a.known, b);
    pr_sep_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }

  std::unordered_set<int> prSepTouched;
  std::unordered_set<uint64_t> pscSeen;

  // sep: log of the exterior arm term without Ctx's prefix, Qz[u+1][v-1] S_q'(u, v).
  double pr_sep_ext(int qp, int u, int v, bool two) {
    if (u + 1 > v - 1 || qp + 1 > n - 1) return NEG_INF;
    ArmQ& a = pr_armq(qp);
    if (prSepTouched.insert(qp).second) {
      ++pr_sep_groups;
      pr_sep_garms += a.nfin;
    }
    ++pr_sep_q;
    if (a.known < v - 1) pr_sep_extend(a, v - 1);
    double S = NEG_INF;
    if (two) {
      const int lo = prXsFirst[u + 1], end = prXsFirst[u + 2];
      const int hi = static_cast<int>(
          std::upper_bound(prArmXp.begin() + lo, prArmXp.begin() + end, v - 1) - prArmXp.begin());
      if (hi > lo) S = a.P1[hi - 1];
    } else {
      S = a.D[v - 1][u + 1];
    }
    const double g = pr_qz(u + 1, v - 1);
    const double se = S > HALF_NEG && g > HALF_NEG ? g + S : NEG_INF;
    if (pr_sep_check) pr_sep_diag(a, qp, u, v, two, S, se);
    return se;
  }

  // sep_check: S against the brute-force sum of ra, se against the exact placement over all
  // arms; the per-term log ratio (Qz[u+1][v-1] / Zarm) / (Qz[u+1][xs-1] Qz[x'+1][v-1]).
  void pr_sep_diag(const ArmQ& a, int qp, int u, int v, bool two, double S, double se) {
    const uint64_t key = (static_cast<uint64_t>(qp) << 43) | (static_cast<uint64_t>(u) << 22) |
                         (static_cast<uint64_t>(v) << 1) | (two ? 1u : 0u);
    if (!pscSeen.insert(key).second) return;
    ++psc_q;
    int lo, hi;
    if (two) {
      lo = prXsFirst[u + 1];
      hi = prXsFirst[u + 2];
    } else {
      lo = prXpFirst[u + 1];
      hi = prXpFirst[v];
    }
    const double g = pr_qz(u + 1, v - 1);
    std::vector<std::pair<double, double>> tm;  // (exact term, log ratio)
    double bs = NEG_INF, ex = NEG_INF;
    for (int pos = lo; pos < hi; ++pos) {
      const int r = two ? pos : prByXp[pos];
      if (prArmXs[r] < u + 1 || prArmXp[r] > v - 1 || a.phi[r] < HALF_NEG) continue;
      bs = log_add(bs, a.ra[r]);
      const double gl = two ? 0.0 : pr_qz(u + 1, prArmXs[r] - 1);
      const double gr = pr_qz(prArmXp[r] + 1, v - 1);
      if (gl < HALF_NEG || gr < HALF_NEG) continue;
      const double t = (gl + gr) + a.phi[r];
      ex = log_add(ex, t);
      if (a.ra[r] > HALF_NEG && g > HALF_NEG) tm.emplace_back(t, (g + (a.ra[r] - a.phi[r])) - (gl + gr));
    }
    if ((S > HALF_NEG) != (bs > HALF_NEG)) ++psc_tab_mis;
    else if (S > HALF_NEG) psc_tab_rel = std::max(psc_tab_rel, std::fabs(std::expm1(S - bs)));
    if ((se > HALF_NEG) != (ex > HALF_NEG)) {
      ++psc_sep_mis;
    } else if (se > HALF_NEG) {
      const double d = se - ex;
      psc_sep_max = std::max(psc_sep_max, std::fabs(d));
      psc_d_sum += d;
      psc_d_abs += std::fabs(d);
      psc_hq[psc_bin(d)] += 1.0;
      for (const auto& [t, lr] : tm) {
        const double w = std::exp(t - ex);
        psc_ht[psc_bin(lr)] += w;
        psc_lr_w += w * lr;
        psc_lr_wabs += w * std::fabs(lr);
        ++psc_terms;
      }
    }
  }

