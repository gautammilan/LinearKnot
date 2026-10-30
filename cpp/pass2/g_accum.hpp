// Pass 2, G accumulators: racc_* (per state, then the b_G cut in racc_finish_trunc / ntr_init),
// and the rec_* growth recursion (--grow-rec on: rec_helper, rec_node, rec_gbuild).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  static double rbytes(const RN& r) {
    size_t c = 0;
    for (int m = 0; m < 3; ++m) c += r.L[m].v.capacity() + r.R[m].v.capacity();
    c += 2 * (r.T[0].capacity() + r.T[1].capacity());
    return 8.0 * static_cast<double>(c) + static_cast<double>(sizeof(RN) + 24);
  }
  int rnL = 3, rnR = 3;
  std::unordered_map<uint64_t, RN> prRG, prRR, prRH;
  std::unordered_map<int, RN> prRKls, prRKrs;
  // R side (layer vectors and kept R lists) of the stored growth nodes, roots, helpers, keys
  double rside_bytes(long long& e) const {
    double b = 0.0;
    e = 0;
    auto add = [&](const RN& r) {
      for (int m = 0; m < 3; ++m) {
        e += static_cast<long long>(r.R[m].v.size());
        b += 8.0 * static_cast<double>(r.R[m].v.capacity());
      }
      e += static_cast<long long>(r.T[1].size());
      b += 16.0 * static_cast<double>(r.T[1].capacity());
    };
    for (const auto* M : {&prRG, &prRR, &prRH})
      for (const auto& kv : *M) add(kv.second);
    for (const auto* M : {&prRKls, &prRKrs})
      for (const auto& kv : *M) add(kv.second);
    return b;
  }
  long long rh_elems = 0, rh_elems0 = 0, rh_built = 0, rh_terms = 0, rh_clears = 0, rn_nodes = 0,
            rn_roots = 0, rn_keys = 0, rn_terms = 0, rn_mlc_drop_l = 0, rn_mlc_drop_r = 0;
  static long long relems(const RN& r) {
    long long c = 0;
    for (int m = 0; m < 3; ++m) c += static_cast<long long>(r.L[m].v.size() + r.R[m].v.size());
    return c + static_cast<long long>(r.T[0].size() + r.T[1].size());
  }
  // dense accumulators of the six output vectors (rel. to rRef), touched ranges
  std::vector<double> rA[6];
  double rRef[6];
  int rLo[6], rHi[6];
  std::vector<int> rTouch[6];
  bool rSparse = false;
  struct RTerm {
    const RN* c;
    double w;
    int sh;
  };
  std::vector<RTerm> rTerms;
  double rSelfR = NEG_INF;  // --grow-r-from-v: the node's own R stop (log)
  void racc_begin() {
    rSparse = g_ntr_K > 0 && (p2_g_cap > 0 || p2_g_mode > 0 || p2_g_eager);
    for (int s = 0; s < 6; ++s) {
      if (rA[s].size() < static_cast<size_t>(n)) rA[s].assign(n, 0.0);
      rRef[s] = NEG_INF;
      rLo[s] = n;
      rHi[s] = -1;
      if (rSparse) rTouch[s].clear();
    }
    rSelfR = NEG_INF;
  }
  // slot of layer m on side (0 L, 1 R), -1 if not built
  int rslot(int side, int m) const { return m < (side ? rnR : rnL) ? side * 3 + m : -1; }
  static int rdst(int m, int sh) { return sh ? std::min(m + 1, 2) : m; }
  // Eager growth freezes the prefix at C_pre. Lazy (the default) uses the scan's C.
  double g_pref(int x) const {
    if (p2_g_eager && stem_ && !stem_->prC.empty())
      return x > 0 ? stem_->prC[static_cast<size_t>(x - 1)] : 0.0;
    return prefix_log(x);
  }
  void racc_self(int x, int y, double at) {
    const StemEst& st = *stem_;
    const double px = g_pref(x);
    rRef[0] = std::max(rRef[0], at + (lLam.empty() ? px : px + lLam[x]));
    if (g_rfv) return;
    rRef[3] = std::max(rRef[3], at + E.w_beta1 + st.prS[y + 1]);
  }
  void racc_refs() {
    if (g_ntr_K > 0) {
      for (const RTerm& t : rTerms) {
        for (const auto& [i, g] : t.c->T[0]) rRef[0] = std::max(rRef[0], t.w + g);
        if (g_rfv) continue;
        if (g_ntr_KR >= 0) {
          for (const auto& [i, g] : t.c->T[1]) rRef[3] = std::max(rRef[3], t.w + g);
        } else if (!t.c->R[0].v.empty()) {
          rRef[3] = std::max(rRef[3], t.w + t.c->R[0].mx);
        }
      }
      return;
    }
    for (const RTerm& t : rTerms)
      for (int side = 0; side < 2; ++side)
        for (int m = 0; m < 3; ++m) {
          const RV& a = side ? t.c->R[m] : t.c->L[m];
          if (a.v.empty()) continue;
          const int s = rslot(side, rdst(m, t.sh));
          if (s >= 0) rRef[s] = std::max(rRef[s], t.w + a.mx);
        }
  }
  void racc_add_self(int x, int y, double at) {
    const StemEst& st = *stem_;
    auto put = [&](int s, int i, double w) {
      if (w < HALF_NEG || rRef[s] < HALF_NEG) return;
      if (rSparse && rA[s][i] == 0.0) rTouch[s].push_back(i);
      rA[s][i] += std::exp(w - rRef[s]);
      rLo[s] = std::min(rLo[s], i);
      rHi[s] = std::max(rHi[s], i);
    };
    const double px = g_pref(x);
    put(0, y, at + (lLam.empty() ? px : px + lLam[x]));
    if (g_rfv) {
      rSelfR = at + E.w_beta1 + st.prS[y + 1];
      return;
    }
    put(3, x, at + E.w_beta1 + st.prS[y + 1]);
  }
  void racc_add_terms() {
    if (g_ntr_K > 0) {
      for (const RTerm& t : rTerms)
        for (int side = 0; side < 2; ++side) {
          const int s = side * 3;
          if (side == 1 && g_rfv) continue;
          if (side == 1 && g_ntr_KR < 0) {
            const RV& a = t.c->R[0];
            if (a.v.empty()) continue;
            const double sc = std::exp((t.w + a.mx) - rRef[s]);
            if (sc == 0.0) continue;
            double* d = rA[s].data() + a.lo;
            for (size_t i = 0; i < a.v.size(); ++i) d[i] += sc * a.v[i];
            if (rSparse)
              for (size_t i = 0; i < a.v.size(); ++i)
                if (a.v[i] != 0.0) rTouch[s].push_back(a.lo + static_cast<int>(i));
            if (lpf) lpf_racc_dense += static_cast<long long>(a.v.size());
            rLo[s] = std::min(rLo[s], a.lo);
            rHi[s] = std::max(rHi[s], a.lo + static_cast<int>(a.v.size()) - 1);
            continue;
          }
          for (const auto& [i, g] : t.c->T[side]) {
            const double v = std::exp((t.w + g) - rRef[s]);
            if (v == 0.0) continue;
            if (rSparse && rA[s][i] == 0.0) rTouch[s].push_back(i);
            rA[s][i] += v;
            rLo[s] = std::min(rLo[s], i);
            rHi[s] = std::max(rHi[s], i);
          }
        }
      return;
    }
    for (const RTerm& t : rTerms)
      for (int side = 0; side < 2; ++side)
        for (int m = 0; m < 3; ++m) {
          const RV& a = side ? t.c->R[m] : t.c->L[m];
          if (a.v.empty()) continue;
          const int s = rslot(side, rdst(m, t.sh));
          if (s < 0) continue;
          const double sc = std::exp((t.w + a.mx) - rRef[s]);
          if (sc == 0.0) continue;
          double* d = rA[s].data() + a.lo;
          const double* src = a.v.data();
          const size_t len = a.v.size();
          for (size_t i = 0; i < len; ++i) d[i] += sc * src[i];
          rLo[s] = std::min(rLo[s], a.lo);
          rHi[s] = std::max(rHi[s], a.lo + static_cast<int>(len) - 1);
        }
  }
  // --grow-node-trunc ranking (static): L entry q' by g uL(q'), uL(q') = sum over the V groups
  // (arm r, y' > q') of Qzbar(0, y'-1) prYw / Qz[xs][x'] / Qzbar(0, q') (tr_uL with the
  // separable fill); entries with uL = 0 dropped (no arm's lK is nonzero there).  R entry p'
  // by g, the p' that are some arm's xs first (lsc_R reads R at an arm's xs).
  std::vector<double> ntrUL;
  std::vector<char> ntrXs;
  long long ntr_in[2] = {}, ntr_kept[2] = {}, ntr_nodes = 0;
  mutable long long ntr_rmiss = 0;  // lsc_R lookups answered by the rest bound
  struct NE {
    int i;
    double g, s;
  };
  std::vector<NE> ntrE;
  int ntr_first_j = -2;
  long long ntr_arm_skip = 0;
  // --l2r-qz-audit on: reads of the scan's Qz past the filled column, total and via pr_qz
  bool qz_audit = false, qz_audit_req = false;
  mutable long long qz_unf = 0, qz_unf_pr = 0, qz_reads = 0;
  void ntr_init() {
    if (!ntrUL.empty()) return;
    const StemEst& st = *stem_;
    pr_ub_init();
    ntr_first_j = cur_j;
    std::vector<double> acc(n + 1, NEG_INF);
    for (int r = 0; r < prA; ++r) {
      const int xs = prArmXs[r], xp = prArmXp[r];
      double iz = NEG_INF;
      if (xs >= 1) {
        if (g_ntr_zbar && !pr_qz_unp) {
          const double z = st.Qzbar(xs, xp);
          if (z > HALF_NEG) iz = (prKap[xs - 1] + prKap[xp]) - z;
        } else {
          iz = pr_iqz(xs, xp);
        }
      }
      if (iz < HALF_NEG) {
        ++ntr_arm_skip;
        continue;
      }
      for (int k = prArmY0[r]; k < prArmY0[r + 1]; ++k) {
        const double z = st.Qzbar(0, prYy[k] - 1);
        if (z > HALF_NEG) acc[prYy[k]] = log_add(acc[prYy[k]], (z + prYw[k]) + iz);
      }
    }
    ntrUL.assign(n, NEG_INF);
    double suf = NEG_INF;
    for (int q = n - 1; q >= 0; --q) {
      suf = log_add(suf, acc[q + 1]);
      const double z = st.Qzbar(0, q);
      ntrUL[q] = suf > HALF_NEG ? suf - (z > HALF_NEG ? z : 0.0) : NEG_INF;
    }
    ntrXs.assign(n, 0);
    for (int r = 0; r < prA; ++r) ntrXs[prArmXs[r]] = 1;
  }
  void racc_finish_trunc(RN& r) {
    ntr_init();
    double tot = NEG_INF;
    r.rest = NEG_INF;
    for (const RTerm& t : rTerms)
      if (t.c->rest > HALF_NEG) r.rest = log_add(r.rest, t.w + t.c->rest);
    for (int side = 0; side < 2; ++side) {
      const size_t K = static_cast<size_t>(side && g_ntr_KR > 0 ? g_ntr_KR : g_ntr_K);
      const int s = side * 3;
      std::vector<std::pair<int, double>>& o = r.T[side];
      o.clear();
      if (side == 1 && g_rfv) {
        double rs = rSelfR;
        for (const RTerm& t : rTerms)
          if (t.c->rsum > HALF_NEG) rs = log_add(rs, t.w + t.c->rsum);
        r.R[0] = RV();
        r.rsum = rs;
        tot = log_add(tot, rs);
        continue;
      }
      if (side == 1 && g_ntr_KR < 0) {
        tot = log_add(tot, rv_finish(s, r.R[0]));
        continue;
      }
      if (rHi[s] < rLo[s]) continue;
      double* a = rA[s].data();
      ntrE.clear();
      auto one = [&](int i) {
        if (a[i] > 0.0) {
          const double g = rRef[s] + std::log(a[i]);
          tot = log_add(tot, g);
          if (side == 0) {
            if (p2_g_suffix) {
              const double sb = stem_->prS[static_cast<size_t>(i + 1)];
              if (sb > HALF_NEG) ntrE.push_back({i, g, g + sb});
            } else if (ntrUL[i] > HALF_NEG) {
              ntrE.push_back({i, g, g + ntrUL[i]});
            }
          } else {
            ntrE.push_back({i, g, ntrXs[i] ? 0.0 : NEG_INF});
          }
        }
        a[i] = 0.0;
      };
      if (rSparse) {
        std::vector<int>& tt = rTouch[s];
        std::sort(tt.begin(), tt.end());
        tt.erase(std::unique(tt.begin(), tt.end()), tt.end());
        for (int i : tt) one(i);
      } else {
        for (int i = rLo[s]; i <= rHi[s]; ++i) one(i);
      }
      ntr_in[side] += static_cast<long long>(ntrE.size());
      auto better = [side](const NE& x, const NE& y) {
        if (side == 1 && x.s != y.s) return x.s > y.s;
        const double vx = side ? x.g : x.s, vy = side ? y.g : y.s;
        if (vx != vy) return vx > vy;
        return x.i < y.i;
      };
      if (ntrE.size() > K) {
        std::nth_element(ntrE.begin(), ntrE.begin() + static_cast<long>(K), ntrE.end(), better);
        if (side == 1)
          for (size_t k = K; k < ntrE.size(); ++k) r.rest = log_add(r.rest, ntrE[k].g);
        ntrE.resize(K);
      }
      std::sort(ntrE.begin(), ntrE.end(), [](const NE& x, const NE& y) { return x.i < y.i; });
      o.reserve(ntrE.size());
      for (const NE& e : ntrE) o.emplace_back(e.i, e.g);
      ntr_kept[side] += static_cast<long long>(o.size());
    }
    ++ntr_nodes;
    r.tot = tot;
    rTerms.clear();
  }
  void racc_finish(RN& r) {
    if (g_ntr_K > 0) {
      racc_finish_trunc(r);
      return;
    }
    double tot = NEG_INF;
    for (int s = 0; s < 6; ++s) {
      RV& o = s < 3 ? r.L[s] : r.R[s - 3];
      const double t = rv_finish(s, o);
      if (t > HALF_NEG) tot = log_add(tot, t);
    }
    r.tot = tot;
    rTerms.clear();
  }
  // slot s of the accumulator into o (dense, relative to its max), reset; log of its sum
  double rv_finish(int s, RV& o) {
    o = RV();
    if (rHi[s] < rLo[s]) return NEG_INF;
    if (lpf) {
      ++lpf_racc_n;
      lpf_racc_fin += rHi[s] - rLo[s] + 1;
    }
    int lo = rLo[s], hi = rHi[s];
    double* a = rA[s].data();
    while (lo <= hi && a[lo] == 0.0) ++lo;
    while (hi >= lo && a[hi] == 0.0) --hi;
    double t = NEG_INF;
    if (lo <= hi) {
      double mxv = 0.0, sum = 0.0;
      for (int i = lo; i <= hi; ++i) {
        mxv = std::max(mxv, a[i]);
        sum += a[i];
      }
      o.lo = lo;
      o.mx = rRef[s] + std::log(mxv);
      o.v.resize(static_cast<size_t>(hi - lo + 1));
      for (int i = lo; i <= hi; ++i) o.v[i - lo] = a[i] / mxv;
      t = rRef[s] + std::log(sum);
    }
    for (int i = rLo[s]; i <= rHi[s]; ++i) a[i] = 0.0;
    return t;
  }
  void rec_hg_trim() {
    if (rh_elems - rh_elems0 <= PR_HG_MAX) return;
    std::unique_ptr<SecTimer> ptm(lpf ? new SecTimer(lpf_hg_sec) : nullptr);
    if (lpf) ++lpf_hg_clears;
    prRH.clear();
    {
      std::unique_ptr<SecTimer> mtm(lpf ? new SecTimer(lpf_mtrim_sec) : nullptr);
      mtrim();
    }
    rh_elems0 = rh_elems;
    ++rh_clears;
    gm_help = 0;
  }
  // --grow-ml-c: keep the C terms with the largest w + tot (ties: first listed)
  void rec_topc(std::vector<RTerm>& t, size_t from, long long& drop) {
    const size_t C = static_cast<size_t>(g_ml_c);
    if (C == 0 || t.size() - from <= C) return;
    std::stable_sort(t.begin() + static_cast<long>(from), t.end(), [](const RTerm& a, const RTerm& b) {
      return a.w + a.c->tot > b.w + b.c->tot;
    });
    drop += static_cast<long long>(t.size() - from - C);
    t.resize(from + C);
  }
  // terms of the helper (xo, y, kind) after its --grow-ml-c cut
  void rec_helper_terms(int xo, int y, int kind, std::vector<RTerm>& t, long long& drop) {
    const StemEst& st = *stem_;
    // Pair set P when --g-pairs is on. Flags off keeps prKeptRight, the current list.
    const std::vector<int>& kr =
        (p2_g_mode > 0 ? st.prKeptRightG[xo] : st.prKeptRight[xo]);
    const int S = pr_ml_span;
    auto close = [&](int a, int b) { return (E.w_ap + E.w_bp + E.at(a, b)) + pair_bonus; };
    t.clear();
    if (kind != 0) {
      const int rf = kind == 1 ? S : cap;
      for (auto jt = std::upper_bound(kr.begin(), kr.end(), y); jt != kr.end() && *jt - y - 1 <= rf; ++jt) {
        const int yo = *jt;
        if (!E.end_ok(xo, yo)) continue;
        double w = E.w_cp * (yo - y - 1);
        if (dz) w += E.run(y + 1, yo - 1);
        t.push_back({&prRG.find(key2(xo, yo))->second, w + close(xo, yo), 0});
      }
    }
    for (auto jt = std::upper_bound(kr.begin(), kr.end(), y + 1); jt != kr.end() && *jt - y - 1 <= S; ++jt) {
      const int yo = *jt;
      const double mb = st.PrMbar(y + 1, yo - 1);
      if (mb < HALF_NEG || !E.end_ok(xo, yo)) continue;
      t.push_back({&prRG.find(key2(xo, yo))->second, mb + close(xo, yo), 0});
    }
    if (prof_on) prof_help_steps += static_cast<long long>(t.size());
    rec_topc(t, 0, drop);
  }
  bool g_in_set(int x, int y) const {
    if (p2_g_mode <= 0 || !stem_ || x < 0 || x >= n) return true;
    const std::vector<int>& kr = stem_->prKeptRightG[static_cast<size_t>(x)];
    return std::binary_search(kr.begin(), kr.end(), y);
  }
  // One enclosing pair in P, the nearest left end that has a right end past y.
  void g_one_step(int x, int y, std::vector<ChainNbr>& nb) const {
    if (!stem_ || p2_g_mode <= 0) return;
    const StemEst& st = *stem_;
    const std::vector<std::vector<int>>& KR = st.prKeptRightG;
    for (int xo = x - 1; xo >= 0; --xo) {
      const std::vector<int>& kr = KR[static_cast<size_t>(xo)];
      auto it = std::upper_bound(kr.begin(), kr.end(), y);
      if (it == kr.end()) continue;
      const int yo = *it;
      const int gapL = x - xo - 1, gapR = yo - y - 1;
      if (gapL + gapR <= cap) {
        const double s = E.stem_interior(xo, x, y, yo);
        if (s > HALF_NEG) {
          nb.push_back({xo, yo, 0, s + pair_bonus});
          return;
        }
      }
      if (E.end_ok(x, y) && E.end_ok(xo, yo) && gapR <= pr_ml_span && y + 1 <= yo - 1) {
        const double mb = st.PrMbar(y + 1, yo - 1);
        if (mb > HALF_NEG) {
          double rl = E.w_cp * gapL;
          if (dz) rl += E.run(xo + 1, x - 1);
          const double br = E.w_bp + E.at(x, y);
          const double close = (E.w_ap + E.w_bp + E.at(xo, yo)) + pair_bonus;
          nb.push_back({xo, yo, 2, ((mb + br) + rl) + close});
          return;
        }
      }
    }
  }
  // --grow-ml-ub: log sum over the helper's terms of e^(w + child tot), cached for the run.
  // Bounds the helper's tot from above: the helper adds, per term, the child's T[0] (a subset of
  // its L entries) and its R side, while the child's tot sums all its L entries and all of R.
  static constexpr double kMlUbMargin = 1e-7;
  std::unordered_map<uint64_t, double> prRHub;
  long long ml_ub_cand = 0, ml_ub_built = 0, ml_ub_nodes = 0;
  double rec_helper_ub(int xo, int y, int kind) {
    const uint64_t key = key3(xo, y, kind);
    auto it = prRHub.find(key);
    if (it != prRHub.end()) return it->second;
    std::vector<RTerm> t;
    long long drop = 0;
    rec_helper_terms(xo, y, kind, t, drop);
    double u = NEG_INF;
    for (const RTerm& e : t)
      if (e.c->tot > HALF_NEG) u = log_add(u, e.w + e.c->tot);
    if (u > HALF_NEG) u += kMlUbMargin * std::max(1.0, std::fabs(u));
    return prRHub.emplace(key, u).first->second;
  }
  // helper of pr_helper's (xo, y, kind) over the RN children (layers unshifted)
  const RN& rec_helper(int xo, int y, int kind) {
    const uint64_t key = key3(xo, y, kind);
    auto it = prRH.find(key);
    if (it != prRH.end()) return it->second;
    std::unique_ptr<SecTimer> ptm(lpf ? new SecTimer(lpf_help_sec) : nullptr);
    if (lpf) {
      ++lpf_help_built;
      if (!lpf_help_keys.insert(key).second) ++lpf_help_rebuilt;
    }
    std::vector<RTerm> t;
    rec_helper_terms(xo, y, kind, t, rn_mlc_drop_r);
    std::vector<RTerm> keep;
    keep.swap(rTerms);
    rTerms = std::move(t);
    racc_begin();
    racc_refs();
    racc_add_terms();
    ++rh_built;
    rh_terms += static_cast<long long>(rTerms.size());
    RN r;
    racc_finish(r);
    rTerms.swap(keep);
    // --grow-r-from-v: helpers hold no R side, so their fixed size counts toward the clear limit
    rh_elems += relems(r) + (g_rfv ? static_cast<long long>(sizeof(RN) + 64) / 8 : 0);
    gm_help += rbytes(r);
    gm_help_entries_peak = std::max(gm_help_entries_peak, rh_elems - rh_elems0);
    gm_update();
    return prRH.emplace(key, std::move(r)).first->second;
  }
  // The multiloop terms rec_topc keeps from s (same terms, same order), building a helper only
  // while its bound (cached tot, else rec_helper_ub) is not below the C-th best tot so far.
  void rec_ml_ub(int y, bool inner, const MlSplit& s, std::vector<RTerm>& ml) {
    struct Cand {
      int xo, kind, idx;
      double w, ub;
    };
    std::vector<Cand> cd;
    cd.reserve(s.qm.size() + s.lr.size());
    for (const auto& [xo, qv] : s.qm) cd.push_back({xo, inner ? 2 : 1, static_cast<int>(cd.size()), qv + s.br, 0.0});
    for (const auto& [xo, rl] : s.lr) cd.push_back({xo, 0, static_cast<int>(cd.size()), rl + s.br, 0.0});
    for (Cand& c : cd) {
      auto it = prRH.find(key3(c.xo, y, c.kind));
      c.ub = c.w + (it != prRH.end() ? it->second.tot : rec_helper_ub(c.xo, y, c.kind));
    }
    std::vector<int> ord(cd.size());
    for (size_t k = 0; k < cd.size(); ++k) ord[k] = static_cast<int>(k);
    std::sort(ord.begin(), ord.end(), [&](int a, int b) { return cd[a].ub > cd[b].ub || (cd[a].ub == cd[b].ub && a < b); });
    const size_t C = static_cast<size_t>(g_ml_c);
    struct Got {
      RTerm t;
      int idx;
      double v;
    };
    std::vector<Got> got;
    std::vector<double> top;  // the C largest v so far, ascending
    for (int k : ord) {
      const Cand& c = cd[k];
      if (top.size() == C && c.ub < top.front()) break;
      const RN& h = rec_helper(c.xo, y, c.kind);
      const RTerm t{&h, c.w, 1};
      const double v = t.w + t.c->tot;
      got.push_back({t, c.idx, v});
      if (top.size() < C) top.insert(std::upper_bound(top.begin(), top.end(), v), v);
      else if (v > top.front()) {
        top.erase(top.begin());
        top.insert(std::upper_bound(top.begin(), top.end(), v), v);
      }
    }
    std::sort(got.begin(), got.end(), [](const Got& a, const Got& b) { return a.v > b.v || (a.v == b.v && a.idx < b.idx); });
    ++ml_ub_nodes;
    ml_ub_cand += static_cast<long long>(cd.size());
    ml_ub_built += static_cast<long long>(got.size());
    rn_mlc_drop_l += static_cast<long long>(cd.size() - C);
    for (size_t k = 0; k < C && k < got.size(); ++k) ml.push_back(got[k].t);
  }
  // RN of x.y with steps nb (pr_steps(x, y, inner, true)): the stop at x.y, interior steps
  // into the children, multiloop steps through the helpers as pr_node_acc with --pr-fast 2.
  void rec_node(int x, int y, bool inner, const std::vector<ChainNbr>& nb, RN& out) {
    rec_hg_trim();
    rTerms.clear();
    for (const ChainNbr& t : nb)
      if (t.t == 0) rTerms.push_back({&prRG.find(key2(t.x, t.y))->second, t.w, 0});
    MlSplit s;
    pr_split(x, y, !geq_pre, s);
    if (s.on) {
      if (lpf) {
        lpf_ml_cand += static_cast<long long>(s.qm.size() + s.lr.size());
        lpf_ml_max = std::max(lpf_ml_max, static_cast<long long>(s.qm.size() + s.lr.size()));
        ++lpf_ml_nodes;
      }
      std::vector<RTerm> ml;
      if (g_ml_ub && g_ml_c > 0 && s.qm.size() + s.lr.size() > static_cast<size_t>(g_ml_c)) rec_ml_ub(y, inner, s, ml);
      else {
        for (const auto& [xo, qv] : s.qm) ml.push_back({&rec_helper(xo, y, inner ? 2 : 1), qv + s.br, 1});
        for (const auto& [xo, rl] : s.lr) ml.push_back({&rec_helper(xo, y, 0), rl + s.br, 1});
        rec_topc(ml, 0, rn_mlc_drop_l);
      }
      rTerms.insert(rTerms.end(), ml.begin(), ml.end());
    }
    rn_terms += static_cast<long long>(rTerms.size());
    const double at = E.end_ok(x, y) ? E.at(x, y) : NEG_INF;
    racc_begin();
    if (at > HALF_NEG) racc_self(x, y, at);
    racc_refs();
    if (at > HALF_NEG) racc_add_self(x, y, at);
    racc_add_terms();
    racc_finish(out);
  }
  // R side entries per built node / root / key: count, sum, max
  long long rz_n[3] = {}, rz_sum[3] = {}, rz_max[3] = {};
  void rec_stat(int who, const RN& r) {
    gz_after[who].push_back(static_cast<int>(relems(r)));
    {
      long long c = static_cast<long long>(r.T[1].size());
      for (int m = 0; m < 3; ++m) c += static_cast<long long>(r.R[m].v.size());
      ++rz_n[who];
      rz_sum[who] += c;
      rz_max[who] = std::max(rz_max[who], c);
    }
    if (g_ntr_K > 0) {
      gz_dq[who].push_back(static_cast<int>(r.T[0].size()));
      gz_dp[who].push_back(static_cast<int>(r.T[1].size()));
      return;
    }
    auto nz = [&](const RV* a, int nl) {
      int lo = n, hi = -1;
      for (int m = 0; m < nl; ++m)
        if (!a[m].v.empty()) {
          lo = std::min(lo, a[m].lo);
          hi = std::max(hi, a[m].lo + static_cast<int>(a[m].v.size()) - 1);
        }
      int c = 0;
      for (int i = lo; i <= hi; ++i) {
        bool z = true;
        for (int m = 0; m < nl && z; ++m)
          if (!a[m].v.empty() && i >= a[m].lo && i < a[m].lo + static_cast<int>(a[m].v.size()) &&
              a[m].v[i - a[m].lo] > 0.0)
            z = false;
        if (!z) ++c;
      }
      return c;
    };
    gz_dq[who].push_back(nz(r.L, rnL));
    gz_dp[who].push_back(nz(r.R, rnR));
  }
  void rec_gbuild(int x, int y) {
    if (prRG.count(key2(x, y))) return;
    std::vector<std::pair<int, int>> stk{{x, y}};
    std::vector<ChainNbr> nb;
    while (!stk.empty()) {
      const auto [a, b] = stk.back();
      if (prRG.count(key2(a, b))) {
        stk.pop_back();
        continue;
      }
      pr_steps(a, b, false, !geq_pre, nb);
      if (nb.empty() && p2_g_fb && p2_g_mode > 0 && !g_in_set(a, b)) g_one_step(a, b, nb);
      bool ready = true;
      for (const ChainNbr& t : nb)
        if (!prRG.count(key2(t.x, t.y))) {
          stk.emplace_back(t.x, t.y);
          ready = false;
        }
      if (!ready) continue;
      count_steps(0, nb);
      RN r;
      rec_node(a, b, false, nb, r);
      ++rn_nodes;
      rec_stat(0, r);
      gm_nodes += rbytes(r);
      gm_update();
      prRG.emplace(key2(a, b), std::move(r));
      stk.pop_back();
    }
  }
  // One root from the swept gamma. g at q' already includes the scan prefix, as rec_root stores it.
  bool sw_g_ready = false;  // one-shot gate for G availability (direct fill or states prGamma)
