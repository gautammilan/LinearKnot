// Diagnostics of V (--pr-v-diag, --pr-vdump): per-group log Vr, interior-only V, T tables.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  // Per V group g (prYy index): log Vr, the interior-only V and log sum w len (len = loops
  // from xs.ys down to x'.y'); arm ids by (xs, x'); the full-support Vr groups not in V.
  void pr_vdiag_static(const std::vector<std::pair<int, int>>& pairs) {
    const StemEst& st = *stem_;
    auto t0 = std::chrono::steady_clock::now();
    pr_vsuf_build();
    const size_t G = prYy.size();
    prVrg.assign(G, NEG_INF);
    prVr2g.assign(G, NEG_INF);
    prVsepg.assign(G, NEG_INF);
    prVog.assign(G, NEG_INF);
    prVlbg.assign(G, NEG_INF);
    std::unordered_map<uint64_t, int> gix;
    gix.reserve(G * 2 + 1);
    for (int r = 0; r < prA; ++r) {
      prArmOf[key2(prArmXs[r], prArmXp[r])] = r;
      for (int g = prArmY0[r]; g < prArmY0[r + 1]; ++g) {
        prVrg[g] = pr_vr_group(prArmXs[r], prArmXp[r], prYy[g]);
        prVsepg[g] = pr_vb_group(prArmXs[r], prArmXp[r], prYy[g], 0);
        prVr2g[g] = pr_vb_group(prArmXs[r], prArmXp[r], prYy[g], 1);
        prVlbg[g] = pr_vb_group(prArmXs[r], prArmXp[r], prYy[g], 2);
        prVog[g] = pr_vb_group(prArmXs[r], prArmXp[r], prYy[g], 3);
        gix[key3(prArmXs[r], prArmXp[r], prYy[g])] = g;
      }
    }
    for (const auto& [x, y] : pairs) {
      if (!E.end_ok(x, y)) continue;
      for (int xs = 0; xs <= x; ++xs) {
        if (pr_vr_group(xs, x, y) < HALF_NEG) continue;
        ++pr_vr_groups;
        if (!gix.count(key3(xs, x, y))) ++pr_vr_extra;
      }
    }
    pr_vr_sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    t0 = std::chrono::steady_clock::now();
    // interior-only V (pr_ml off): profiles of w and of log sum w len
    const bool ml0 = pr_ml;
    pr_ml = false;
    std::unordered_map<uint64_t, std::pair<PProf, PProf>> VI;
    VI.reserve(pairs.size() * 2 + 1);
    std::vector<ChainNbr> nb;
    auto merge2 = [&](PAcc& aw, PAcc& al) {
      for (const ChainNbr& t : nb) {
        const auto& q = VI.find(key2(t.x, t.y))->second;
        for (size_t k = 0; k < q.first.k.size(); ++k) aw.add(q.first.k[k], t.w + q.first.w[k]);
        for (size_t k = 0; k < q.first.k.size(); ++k) al.add(q.first.k[k], t.w + q.first.w[k]);
        for (size_t k = 0; k < q.second.k.size(); ++k) al.add(q.second.k[k], t.w + q.second.w[k]);
      }
    };
    for (const auto& [x, y] : pairs) {
      pr_steps(x, y, false, false, nb);
      PAcc aw, al;
      if (E.end_ok(x, y)) aw.add(key2(x, y), E.at(x, y));
      merge2(aw, al);
      std::pair<PProf, PProf> p;
      aw.finish(p.first, 0.0);
      al.finish(p.second, 0.0);
      VI.emplace(key2(x, y), std::move(p));
    }
    prVint.assign(G, NEG_INF);
    prVlen.assign(G, NEG_INF);
    for (const auto& [x, y] : pairs) {
      if (!E.end_ok(x, y)) continue;
      pr_steps(x, y, true, false, nb);
      PAcc aw, al;
      if (single_stems) aw.add(key2(x, y), E.at(x, y));
      merge2(aw, al);
      PProf pw, pl;
      aw.finish(pw, 0.0);
      al.finish(pl, 0.0);
      const double base = (E.at(x, y) + E.w_stem) + pair_bonus;
      for (size_t k = 0; k < pw.k.size(); ++k) {
        auto it = gix.find(key3(k_hi(pw.k[k]), x, y));
        if (it == gix.end()) continue;
        const double f = base + (E.w_beta1 + st.prS[k_lo(pw.k[k]) + 1]);
        prVint[it->second] = log_add(prVint[it->second], pw.w[k] + f);
      }
      for (size_t k = 0; k < pl.k.size(); ++k) {
        auto it = gix.find(key3(k_hi(pl.k[k]), x, y));
        if (it == gix.end()) continue;
        const double f = base + (E.w_beta1 + st.prS[k_lo(pl.k[k]) + 1]);
        prVlen[it->second] = log_add(prVlen[it->second], pl.w[k] + f);
      }
    }
    pr_ml = ml0;
    pr_vint_sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    prVuse.assign(G, 0);
    prVcon.assign(G, NEG_INF);
    prVshare.assign(G, 0.0);
  }

  // A chosen arm set of an L-lookup (q'; u, v): per group of the selected arms with y' > q',
  // the use, log sum of Qzbar[q'+1][y'-1] V (no gap fillers), and its share of the lookup's
  // exterior total (with the fillers).  Sampled lookups: every arm of the gap ranked by its
  // partner value with V and with Vr (full support).
  std::vector<std::pair<int, double>> prVtmp;
  void pr_vdiag_T(int qp, int u, int v, bool two, const std::vector<int>& sel) {
    const StemEst& st = *stem_;
    ++pr_vT_n;
    if (!prVuse.empty()) {
      prVtmp.clear();
      double tot = NEG_INF;
      for (int r : sel) {
        const double gl = two ? 0.0 : pr_qz(u + 1, prArmXs[r] - 1);
        const double gr = pr_qz(prArmXp[r] + 1, v - 1);
        if (gl < HALF_NEG || gr < HALF_NEG) continue;
        const int hi = prArmY0[r + 1];
        for (int g = static_cast<int>(std::upper_bound(prYy.begin() + prArmY0[r], prYy.begin() + hi, qp) -
                                      prYy.begin());
             g < hi; ++g) {
          const double z = st.Qzbar(qp + 1, prYy[g] - 1);
          if (z < HALF_NEG) continue;
          const double c = z + prYw[g];
          prVtmp.emplace_back(g, c + gl + gr);
          tot = log_add(tot, c + gl + gr);
          ++prVuse[g];
          prVcon[g] = log_add(prVcon[g], c);
        }
      }
      if (tot > HALF_NEG)
        for (const auto& [g, t] : prVtmp) prVshare[g] += std::exp(t - tot);
    }
    if (!pr_vsel || pr_vsel_n >= pr_vsample_max || pr_vrand() >= pr_vsample_p) return;
    ++pr_vsel_n;
    // exact partner value per arm of the gap
    std::vector<std::pair<double, int>> pv;
    std::unordered_map<int, double> pvOf;
    int lo, hi;
    if (two) {
      lo = prXsFirst[u + 1];
      hi = prXsFirst[u + 2];
    } else {
      lo = prXpFirst[u + 1];
      hi = prXpFirst[v];
    }
    double totV = NEG_INF;
    for (int pos = lo; pos < hi; ++pos) {
      const int r = two ? pos : prByXp[pos];
      if (prArmXp[r] > v - 1 || (!two && prArmXs[r] < u + 1)) continue;
      const double gl = two ? 0.0 : pr_qz(u + 1, prArmXs[r] - 1);
      const double gr = pr_qz(prArmXp[r] + 1, v - 1);
      if (gl < HALF_NEG || gr < HALF_NEG) continue;
      double s = NEG_INF;
      const int gh = prArmY0[r + 1];
      for (int g = static_cast<int>(std::upper_bound(prYy.begin() + prArmY0[r], prYy.begin() + gh, qp) -
                                    prYy.begin());
           g < gh; ++g) {
        const double z = st.Qzbar(qp + 1, prYy[g] - 1);
        if (z > HALF_NEG) s = log_add(s, z + prYw[g]);
      }
      if (s < HALF_NEG) continue;
      s += gl + gr;
      pv.emplace_back(s, r);
      pvOf[r] = s;
      totV = log_add(totV, s);
    }
    // partner value per arm (xs, x') of the gap, full support: kind -1 Vr, else pr_vb_group's
    const auto& KR = pr_KRp();
    auto rank = [&](int kind, std::vector<std::tuple<double, int, int>>& pr) {
      pr.clear();
      for (int xp = u + 1; xp <= v - 1; ++xp) {
        const int xs_lo = u + 1, xs_hi = two ? std::min(u + 1, xp) : xp;
        if (xs_lo > xs_hi) continue;
        const double gr = pr_qz(xp + 1, v - 1);
        if (gr < HALF_NEG) continue;
        std::vector<double> acc(xs_hi - xs_lo + 1, NEG_INF);
        for (auto it = std::upper_bound(KR[xp].begin(), KR[xp].end(), qp); it != KR[xp].end(); ++it) {
          const int yp = *it;
          if (!E.end_ok(xp, yp)) continue;
          const double z = st.Qzbar(qp + 1, yp - 1);
          if (z < HALF_NEG) continue;
          for (int xs = xs_lo; xs <= xs_hi; ++xs) {
            const double w = kind < 0 ? pr_vr_group(xs, xp, yp) : pr_vb_group(xs, xp, yp, kind);
            if (w > HALF_NEG) acc[xs - xs_lo] = log_add(acc[xs - xs_lo], z + w);
          }
        }
        for (int xs = xs_lo; xs <= xs_hi; ++xs) {
          if (acc[xs - xs_lo] < HALF_NEG) continue;
          const double gl = two ? 0.0 : pr_qz(u + 1, xs - 1);
          if (gl < HALF_NEG) continue;
          pr.emplace_back(acc[xs - xs_lo] + gl + gr, xs, xp);
        }
      }
      std::sort(pr.begin(), pr.end(), [](const auto& a, const auto& b) { return std::get<0>(a) > std::get<0>(b); });
    };
    const size_t K = static_cast<size_t>(pr_arms_k > 0 ? pr_arms_k : 64);
    std::sort(pv.begin(), pv.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
    std::unordered_set<int> topV;
    double keepV = NEG_INF, keepE = NEG_INF;
    for (size_t k = 0; k < std::min(K, pv.size()); ++k) {
      topV.insert(pv[k].second);
      keepV = log_add(keepV, pv[k].first);
    }
    // per ranking: arms ranked, overlap with V's top K, of those in V, log kept exact V
    std::vector<std::tuple<double, int, int>> pr;
    size_t nr[3], ov[3], inv[3];
    double keep[3];
    const int KINDS[3] = {-1, 1, 0};
    std::unordered_map<uint64_t, std::array<double, 4>> armv;  // exact, Vr, Vr2, Vsep
    if (pr_varms)
      for (const auto& [s, r] : pv) armv[key2(prArmXs[r], prArmXp[r])] = {s, NEG_INF, NEG_INF, NEG_INF};
    for (int c = 0; c < 3; ++c) {
      rank(KINDS[c], pr);
      if (pr_varms)
        for (const auto& [val, a, b] : pr) {
          auto it = armv.try_emplace(key2(a, b), std::array<double, 4>{NEG_INF, NEG_INF, NEG_INF, NEG_INF}).first;
          it->second[1 + c] = val;
        }
      nr[c] = pr.size();
      ov[c] = inv[c] = 0;
      keep[c] = NEG_INF;
      for (size_t k = 0; k < std::min(K, pr.size()); ++k) {
        auto ia = prArmOf.find(key2(std::get<1>(pr[k]), std::get<2>(pr[k])));
        if (ia == prArmOf.end()) continue;
        auto ip = pvOf.find(ia->second);
        if (ip == pvOf.end()) continue;
        ++inv[c];
        keep[c] = log_add(keep[c], ip->second);
        if (topV.count(ia->second)) ++ov[c];
      }
    }
    for (int r : sel) {
      auto ip = pvOf.find(r);
      if (ip != pvOf.end()) keepE = log_add(keepE, ip->second);
    }
    auto sh = [&](double x) { return totV > HALF_NEG && x > HALF_NEG ? std::exp(x - totV) : 0.0; };
    std::fprintf(pr_vsel, "%d\t%d\t%d\t%d\t%zu\t%zu\t%zu\t%zu\t%.9g\t%.9g\t%.9g\t%zu\t%zu\t%.9g\t%zu\t%zu\t%.9g\n",
                 qp, u, v, two ? 1 : 0, pv.size(), nr[0], ov[0], inv[0], sh(keep[0]), sh(keepV), sh(keepE),
                 ov[1], inv[1], sh(keep[1]), ov[2], inv[2], sh(keep[2]));
    if (pr_varms)
      for (const auto& [k, a] : armv)
        std::fprintf(pr_varms, "%lld\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\n", pr_vsel_n, k_hi(k), k_lo(k), a[0],
                     a[1], a[2], a[3]);
  }

  void pr_vdump(const std::string& sid) {
    if (prVrg.empty()) return;
    const StemEst& st = *stem_;
    if (!pr_vdiag.empty()) {
      const std::string cmd = "gzip -1 > '" + pr_vdiag + ".ventries.tsv.gz'";
      FILE* f = popen(cmd.c_str(), "w");
      std::fprintf(f, "xs\txp\typ\tlogV\tlogVr\tlogVint\tlogVlen\tlpn\tuse\tlogcon\tshare\t"
                      "logVr2\tlogVsep\tlogVo\tlogLB\n");
      for (int r = 0; r < prA; ++r)
        for (int g = prArmY0[r]; g < prArmY0[r + 1]; ++g) {
          const double lp = st.prLP.empty() ? NAN : st.prLP[static_cast<size_t>(prArmXp[r]) * n + prYy[g]];
          std::fprintf(f, "%d\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.6g\t%lld\t%.9g\t%.6g\t%.9g\t%.9g\t%.9g\t%.9g\n",
                       prArmXs[r], prArmXp[r], prYy[g], prYw[g], prVrg[g], prVint[g], prVlen[g], lp,
                       prVuse[g], prVcon[g], prVshare[g], prVr2g[g], prVsepg[g], prVog[g], prVlbg[g]);
        }
      pclose(f);
    }
    std::fprintf(stderr,
                 "[vratio] %s groups %zu vr_groups %lld vr_extra %lld used %lld T_seen %lld "
                 "sel_samples %lld vr_sec %.3f vint_sec %.3f\n",
                 sid.c_str(), prYy.size(), pr_vr_groups, pr_vr_extra,
                 static_cast<long long>(std::count_if(prVuse.begin(), prVuse.end(),
                                                      [](long long u) { return u > 0; })),
                 pr_vT_n, pr_vsel_n, pr_vr_sec, pr_vint_sec);
  }

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
    std::vector<double> phi, rs;  // rs: --pr-arms-rank ratio, phi - wu (x' - xs + 1)
    std::vector<int> t1, t2;
    int sz = 0;
    // sep: ra = phi - log Zarm; D[b][x] = log sum ra over arms xs >= x, x' <= b; P1[r] = log
    // sum ra over the arms of r's xs with x' <= x'_r; columns 0..known built (Qz final there)
    std::vector<double> ra, P1;
    std::vector<std::vector<double>> D;
    int known = -1, nfin = 0;
    std::unordered_map<int, std::vector<std::pair<int, double>>> W;
    long long use = 0;
    // gaponly: argmax trees over the indexed arms only (positions of order 1 / 2), built for
    // the gap union `epoch`
    CTree g1, g2;
    long long epoch = -1;
  };
  static constexpr size_t PR_AQ_MAX = 128;
  std::unordered_map<int, ArmQ> prAQ;
  long long pr_use = 0;

  // gaponly: merged x' intervals [l, r] of the gaps (u, v) scored since the last
  // pr_gap_begin (u+1..v-1), and its epoch.
  std::vector<std::pair<int, int>> prGapIv;
  long long gap_epoch = 0;
