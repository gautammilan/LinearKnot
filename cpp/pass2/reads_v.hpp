// Pass-2 reads of V for pass 3 (--grow-r-from-v): kind maps, rfv_* (growth roots from V),
// mu dumps, V queries and traces.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  int v_map_kind(int x, int y, int xs) const {
    if (!stem_) return 3;
    auto it = stem_->prMu.find(key2(x, y));
    if (it == stem_->prMu.end()) return 3;
    auto jt = it->second.find(key2(xs, 0));
    if (jt != it->second.end() && jt->second > HALF_NEG) return 0;
    if (p2_v_cap > 0 && static_cast<int>(it->second.size()) >= p2_v_cap) return 1;
    return 2;
  }
  void note_kind(long long* a, long long* b, long long* c, int st, int k) {
    if (st < 0 || st >= ST_N) return;
    if (k == 1) ++a[st];
    else if (k == 2) ++b[st];
    else if (k == 3) ++c[st];
  }
  int k_map_kind(int r, int qp) const {
    if (r < 0 || qp + 1 > n - 1 || prArmY0.empty()) return 3;
    const int xs = prArmXs[r], xp = prArmXp[r], hi = prArmY0[r + 1];
    const int* yy = prYy.data();
    const int k = static_cast<int>(std::upper_bound(yy + prArmY0[r], yy + hi, qp) - yy);
    if (k >= hi) return 3;
    const int mk = v_map_kind(xp, yy[k], xs);
    return mk == 0 ? 3 : mk;
  }
  int h_map_kind(int a2, int j) const {
    if (a2 < 0 || a2 >= n || !stem_) return 3;
    bool full = false, any = false;
    for (int c1 : partners_right[a2]) {
      if (c1 <= j) continue;
      auto it = stem_->prMu.find(key2(a2, c1));
      if (it == stem_->prMu.end()) continue;
      any = true;
      if (p2_v_cap > 0 && static_cast<int>(it->second.size()) >= p2_v_cap) full = true;
    }
    if (!any) return 3;
    return full ? 1 : 2;
  }
  void rfv_put(uint64_t key, const PProf& p) {
    prVkI[key] = {static_cast<uint32_t>(prVkX.size()), static_cast<uint32_t>(p.k.size())};
    for (size_t i = 0; i < p.k.size(); ++i) {
      prVkX.push_back(k_hi(p.k[i]));
      prVkW.push_back(static_cast<float>(p.w[i]));
    }
  }
  void rfv_keep(std::unordered_map<uint64_t, PProf>& V) {
    for (const auto& kv : V) rfv_put(kv.first, kv.second);
    prVkX.shrink_to_fit();
    prVkW.shrink_to_fit();
  }
  bool rfv_get(int x, int y, VSpan& s) const {
    auto it = prVkI.find(key2(x, y));
    if (it == prVkI.end()) return false;
    s = {prVkX.data() + it->second.first, prVkW.data() + it->second.first, it->second.second};
    return true;
  }
  static double rfv_at(const VSpan& s, int xs) {
    const int32_t* e = s.xs + s.len;
    const int32_t* it = std::lower_bound(s.xs, e, xs);
    return it != e && *it == xs ? static_cast<double>(s.lw[it - s.xs]) : NEG_INF;
  }
  static void rfv_prof(const VSpan& s, PProf& p) {
    p = PProf();
    for (uint32_t i = 0; i < s.len; ++i) p.mx = std::max(p.mx, static_cast<double>(s.lw[i]));
    p.k.resize(s.len);
    p.lin.resize(s.len);
    for (uint32_t i = 0; i < s.len; ++i) {
      p.k[i] = key2(s.xs[i], 0);
      p.lin[i] = std::exp(static_cast<double>(s.lw[i]) - p.mx);
    }
  }
  VSpan rfv_node(int x, int y, bool inner) {
    const uint64_t key = key3(x, y, inner ? 1 : 0);
    auto it = prVod.find(key);
    if (it != prVod.end()) return it->second.span();
    const auto rfv_t0 = std::chrono::steady_clock::now();
    const StemEst& st = *stem_;
    std::vector<ChainNbr> nb;
    pr_steps(x, y, inner, false, nb, pr_fast < 2);
    const double self = E.end_ok(x, y) ? E.at(x, y) + E.w_beta1 + st.prS[y + 1] : NEG_INF;
    static const PProf empty;
    auto child = [&](int a, int b) -> const PProf* {
      const uint64_t k = key2(a, b);
      auto jt = prVch.find(k);
      if (jt != prVch.end()) return &jt->second;
      VSpan s;
      if (!rfv_get(a, b, s)) {
        ++rfv_child_miss;
        return &empty;
      }
      PProf& p = prVch[k];
      rfv_prof(s, p);
      return &p;
    };
    if (prVodH.size() > 100000) prVodH.clear();
    pr_node_acc(x, y, inner, false, nb, self, child, prVodH);
    prVch.clear();
    PProf p;
    pr_finish_fut_d(p, x);
    VOwn o;
    o.xs.reserve(p.k.size());
    o.lw.reserve(p.k.size());
    for (size_t i = 0; i < p.k.size(); ++i) {
      o.xs.push_back(k_hi(p.k[i]));
      o.lw.push_back(static_cast<float>(p.w[i]));
    }
    ++rfv_od;
    rfv_od_ent += static_cast<long long>(o.xs.size());
    rfv_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - rfv_t0).count();
    return prVod.emplace(key, std::move(o)).first->second.span();
  }
  double mu_at(int x, int y, int xs, int st = -1) {
    const int kind = v_map_kind(x, y, xs);
    if (kind == 0) {
      auto it = stem_->prMu.find(key2(x, y));
      return it->second.find(key2(xs, 0))->second;
    }
    note_kind(vm_i, vm_ii, vm_iii, st, kind);
    return NEG_INF;
  }
  // V of the root x.y (seed: x.y is the stem's inner pair) at the outer left end xs
  double rfv_root(int x, int y, bool seed, int xs) {
    if (p2_vsweep) {
      ++rfv_look;
      ++rfv_look_r;
      const double v = mu_at(x, y, xs, ST_R);
      if (v < HALF_NEG) {
        ++rfv_miss;
        ++rfv_miss_r;
        ++rfv_key0_r;
      }
      return v;
    }
    ++rfv_look;
    VSpan p;
    const bool kept = rfv_get(x, y, p);
    double v = NEG_INF;
    if (!kept || (seed && !E.end_ok(x, y))) {
      v = rfv_at(rfv_node(x, y, seed), xs);
    } else if (!seed) {
      v = rfv_at(p, xs);
    } else if (xs == x) {
      v = E.at(x, y) + E.w_beta1 + stem_->prS[y + 1];
    } else {
      ++rfv_seed;
      const int r = larm(xs, x);
      if (r >= 0) {
        const auto b = prYy.begin();
        auto it = std::lower_bound(b + prArmY0[r], b + prArmY0[r + 1], y);
        if (it != b + prArmY0[r + 1] && *it == y)
          v = prYw[static_cast<size_t>(it - b)] - ((E.at(x, y) + E.w_stem) + pair_bonus);
      }
    }
    if (v < HALF_NEG) ++rfv_miss;
    return v;
  }
  // U key (ls / rs of p at step q): sum over its closings c of w_c V(c.x, c.y) at xs
  double rfv_key(bool ls, int p, bool seed, int q, int xs) {
    if (p2_vsweep) {
      ++rfv_look;
      ++rfv_look_u;
      std::vector<Closing> cl;
      if (ls) pr_ls_closings(p, seed, q, cl);
      else closings_rs(p, q, cl);
      double s = NEG_INF;
      int missed = 0;
      for (const Closing& c : cl) {
        const double v = mu_at(c.x, c.y, xs, ST_U);
        if (v > HALF_NEG) s = log_add(s, c.w + v);
        else ++missed;
      }
      if (s < HALF_NEG) {
        ++rfv_miss;
        ++rfv_miss_u;
        ++rfv_key0_u;
      } else if (missed > 0) {
        ++rfv_part_u;
      }
      return s;
    }
    auto& m = prVkey[ls ? 0 : 1];
    const int key = ls ? 2 * p + (seed ? 1 : 0) : p;
    auto it = m.find(key);
    if (it == m.end()) {
      std::vector<Closing> cl;
      if (ls) pr_ls_closings(p, seed, q, cl);
      else closings_rs(p, q, cl);
      if (rfvAcc.size() < static_cast<size_t>(n)) rfvAcc.assign(n, NEG_INF);
      rfvT.clear();
      for (const Closing& c : cl) {
        VSpan s;
        if (!rfv_get(c.x, c.y, s)) s = rfv_node(c.x, c.y, false);
        for (uint32_t i = 0; i < s.len; ++i) {
          const int a = s.xs[i];
          if (rfvAcc[a] < HALF_NEG) rfvT.push_back(a);
          rfvAcc[a] = log_add(rfvAcc[a], c.w + static_cast<double>(s.lw[i]));
        }
      }
      std::sort(rfvT.begin(), rfvT.end());
      std::vector<std::pair<int, double>> o;
      o.reserve(rfvT.size());
      for (int a : rfvT) {
        o.emplace_back(a, rfvAcc[a]);
        rfvAcc[a] = NEG_INF;
      }
      ++rfv_keys;
      rfv_key_ent += static_cast<long long>(o.size());
      it = m.emplace(key, std::move(o)).first;
    }
    ++rfv_look;
    const auto& v = it->second;
    auto jt = std::lower_bound(v.begin(), v.end(), xs,
                               [](const std::pair<int, double>& e, int x) { return e.first < x; });
    if (jt != v.end() && jt->first == xs) return jt->second;
    ++rfv_miss;
    return NEG_INF;
  }
  double rfv_bytes(long long& e) const {
    e = static_cast<long long>(prVkX.size());
    double b = msz::umap(prVkI) + msz::vec(prVkX) + msz::vec(prVkW) + msz::umap(prVod);
    for (const auto& kv : prVod) {
      e += static_cast<long long>(kv.second.xs.size());
      b += msz::vec(kv.second.xs) + msz::vec(kv.second.lw);
    }
    return b;
  }

 private:
