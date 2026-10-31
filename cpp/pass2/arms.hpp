// Pass 2, arms: gap coverage, Sbar_pk (pr_sf_build), arm_raw / arm_phi
// (K^pre_L(q'; b, c) = sum_{f>q'} Q^z_pre(q'+1, f-1) V), persistent argmax index, arm
// queries (ArmQ cache prAQ, pr_armq), arm selection (pr_select).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  // Per future end q': phi(arm) = log sum over its entries with y' > q' of
  // Qzbar[q'+1][y'-1] w e^beta1 Sbar[ys+1] (the gap-free future factor), argmax trees over the
  // (xs, x') and (x', xs) orders (k > 0), and W(arm)(ys) = sum_{y' > q'} Qzbar[q'+1][y'-1] w
  // (lazy).  At most PR_AQ_MAX kept, least recently used rebuilt on demand.
  // Argmax tree over a sorted list of positions (ties: lower position, as aq_better).
  struct CTree {
    std::vector<int> pos, t;
    std::vector<double> val;
    int sz = 0;
    int better(int i, int j) const {
      if (i < 0) return j;
      if (j < 0) return i;
      if (val[i] != val[j]) return val[i] > val[j] ? i : j;
      return std::min(i, j);
    }
    void build() {
      sz = 1;
      while (sz < std::max<int>(1, static_cast<int>(pos.size()))) sz <<= 1;
      t.assign(2 * static_cast<size_t>(sz), -1);
      for (int i = 0; i < static_cast<int>(pos.size()); ++i) t[sz + i] = i;
      for (int k = sz - 1; k >= 1; --k) t[k] = better(t[2 * k], t[2 * k + 1]);
    }
    // index (into pos) of the best position in [l, r], -1 if none
    int query(int l, int r) const {
      const int a = static_cast<int>(std::lower_bound(pos.begin(), pos.end(), l) - pos.begin());
      const int b = static_cast<int>(std::upper_bound(pos.begin(), pos.end(), r) - pos.begin()) - 1;
      int best = -1;
      for (int lo = a + sz, hi = b + sz + 1; lo < hi; lo >>= 1, hi >>= 1) {
        if (lo & 1) best = better(best, t[lo++]);
        if (hi & 1) best = better(best, t[--hi]);
      }
      return best;
    }
  };
  struct ArmQ {
    std::vector<double> phi, rs;  // rs: phi - wu (x' - xs + 1)
    std::vector<int> t1, t2;
    int sz = 0;
    // sep: ra = phi - log Zarm; D[b][x] = log sum ra over arms xs >= x, x' <= b; P1[r] = log
    // sum ra over the arms of r's xs with x' <= x'_r; columns 0..known built (Qz final there)
    std::vector<double> ra, P1;
    std::vector<std::vector<double>> D;
    int known = -1;
    std::unordered_map<int, std::vector<std::pair<int, double>>> W;
    long long use = 0;
    // gaponly: argmax trees over the indexed arms only (positions of order 1 / 2), built for
    // the gap union `epoch`
  };
  std::unordered_map<int, ArmQ> prAQ;
  long long pr_use = 0;


  // phi(arm r) at q' as Step A computes it (separable: from the suffix sums prSF).
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
    if (!pr_aq_sep) {
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
    int l, r, best;
    double val;
  };

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
      it->second.use = ++pr_use;
      return it->second;
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
    pr_sf_build();
    const double z0 = pr_aq_sep ? st.Qzbar(0, qp) : 0.0;
    a.phi.assign(prA, NEG_INF);
    aq_arms += prA;
    for (int r = 0; r < prA; ++r) {
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
  void pr_select(const ArmQ& a, int u, int v, bool two, std::vector<int>& out) {
    out.clear();
    if (u + 1 > v - 1) return;
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
    auto push = [&](int l, int r) {
      if (l > r) return;
      int p;
      double val;
      p = aq_query(a, ord, l, r);
      if (p < 0) return;
      val = aq_val(a, ord, p);
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


  std::unordered_set<int> prSepTouched;



