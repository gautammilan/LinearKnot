// Pass 2, static V^pre build (pr_build_static): per arm the V entries from the nested
// pre-pass.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  void pr_build_static() {
    if (pr_built) return;
    pr_built = true;
    SecTimer ts(pr_static_sec);
    const StemEst& st = *stem_;
    std::vector<std::pair<int, int>> pairs;
    const std::vector<std::vector<int>>& vsrc = pr_pair_list(false);
    for (int x = 0; x < n; ++x)
      for (int y : vsrc[x]) pairs.emplace_back(x, y);
    if (!pr_v2dump.empty()) {
      pr_v2_dump(pairs);
      std::fflush(nullptr);
      std::_Exit(0);
    }
    std::stable_sort(pairs.begin(), pairs.end(),
                     [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
                       return a.second - a.first > b.second - b.first;
                     });
    std::vector<AT> ats;
    std::vector<std::vector<AT>> bk;
    if (pr_vstream) bk.resize(n);
    if (lpf) {
      vi_pairs = static_cast<long long>(pairs.size());
      vi_pairs_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - ts.t0).count();
    }
    {
      std::unique_ptr<SecTimer> vtm(lpf ? new SecTimer(vi_vbuild_s) : nullptr);
      pr_static_ats(pr_pmax, true, pairs, ats, pr_vstream ? &bk : nullptr);
    }
    if (sw_check) {
      sw_pairs = pairs;
      sw_verify();
      if (p2_vsweep && stem_ && !stem_->prMu.empty() && !sw_ref.empty()) {
        auto rel = [](double a, double b) {
          if (a < HALF_NEG && b < HALF_NEG) return 0.0;
          if (a < HALF_NEG || b < HALF_NEG) return 1e300;
          return std::fabs(std::expm1(a - b));
        };
        double worst = 0, wgot = 0, wref = 0;
        int wx = -1, wy = -1, wb = -1;
        int fx = -1, fy = -1, fb = -1;
        double fgot = 0, fref = 0;
        long long nent = 0, miss = 0;
        for (const auto& [pk, p] : sw_ref) {
          const int x = k_hi(pk), y = k_lo(pk);
          auto it = stem_->prMu.find(pk);
          for (size_t i = 0; i < p.k.size(); ++i) {
            const int b = k_hi(p.k[i]);
            const double ref = p.w[i];
            double got = NEG_INF;
            if (it != stem_->prMu.end()) {
              auto jt = it->second.find(key2(b, 0));
              if (jt != it->second.end()) got = jt->second;
            }
            if (got < HALF_NEG && ref < HALF_NEG) continue;
            ++nent;
            if (got < HALF_NEG || ref < HALF_NEG) ++miss;
            const double d = rel(got, ref);
            if (d > worst) {
              worst = d;
              wx = x;
              wy = y;
              wb = b;
              wgot = got;
              wref = ref;
            }
            if (d > 1e-12 && (fx < 0 || y < fy || (y == fy && x < fx) || (y == fy && x == fx && b < fb))) {
              fx = x;
              fy = y;
              fb = b;
              fgot = got;
              fref = ref;
            }
          }
        }
        std::fprintf(stderr,
                     "[v-sweep-check] entries %lld one_sided %lld max_rel %.3e at %d.%d b %d "
                     "sweep %.9g ref %.9g first %d.%d b %d sweep %.9g ref %.9g\n",
                     nent, miss, worst, wx, wy, wb, wgot, wref, fx, fy, fb, fgot, fref);
      }
    }
    std::unique_ptr<SecTimer> atm(lpf ? new SecTimer(vi_arm_s) : nullptr);
    mtrim();
    long long bk_n = 0;
    double bk_b = 0.0;
    for (const std::vector<AT>& v : bk) {
      bk_n += static_cast<long long>(v.size());
      bk_b += msz::vec(v);
    }
    mem_log("V_ats_done", {{"V_ats", static_cast<long long>(ats.size()) + bk_n, msz::vec(ats) + bk_b}});
    std::sort(ats.begin(), ats.end(), at_less);
    if (pr_pmax > 0 && pr_pcheck && n <= 40 && pr_cut <= 0.0) {
      std::map<std::tuple<int, int, int, int>, double> en;
      std::vector<ChainNbr> nb0;
      std::function<void(int, int, int, int, int, double)> dfs = [&](int x0, int y0, int x,
                                                                    int y, int c, double w) {
        const bool inner = c == 1;
        if (inner ? single_stems : E.end_ok(x, y)) {
          auto key = std::make_tuple(x, x0, y0, y);
          auto it = en.find(key);
          const double v = w + E.at(x, y);
          if (it == en.end()) en.emplace(key, v);
          else it->second = log_add(it->second, v);
        }
        if (c >= pr_pmax) return;
        std::vector<ChainNbr> nb;
        pr_steps(x, y, inner, false, nb);
        for (const ChainNbr& t : nb) dfs(x0, y0, t.x, t.y, c + 1, w + t.w);
      };
      for (const auto& [x, y] : pairs)
        if (E.end_ok(x, y)) dfs(x, y, x, y, 1, (E.at(x, y) + E.w_stem) + pair_bonus);
      std::vector<AT> ref;
      for (const auto& [k, w] : en)
        ref.push_back({std::get<0>(k), std::get<1>(k), std::get<2>(k), std::get<3>(k), w});
      std::sort(ref.begin(), ref.end(), at_less);
      size_t a = 0, b = 0;
      while (a < ats.size() || b < ref.size()) {
        if (b == ref.size() || (a < ats.size() && at_less(ats[a], ref[b]))) {
          ++pe_key_mis;
          ++a;
        } else if (a == ats.size() || at_less(ref[b], ats[a])) {
          ++pe_key_mis;
          ++b;
        } else {
          ++pe_n;
          pe_maxd = std::max(pe_maxd, std::fabs(ats[a].w - ref[b].w));
          ++a;
          ++b;
        }
      }
    }
    if (pr_pmax > 0 && pr_pcheck) {
      std::vector<AT> ref;
      pr_static_ats(0, false, pairs, ref);
      std::sort(ref.begin(), ref.end(), at_less);
      size_t a = 0, b = 0;
      while (a < ats.size() || b < ref.size()) {
        if (b == ref.size() || (a < ats.size() && at_less(ats[a], ref[b]))) {
          ++pc_key_mis;
          ++a;
        } else if (a == ats.size() || at_less(ref[b], ats[a])) {
          ++pc_key_mis;
          ++b;
        } else {
          ++pc_n;
          pc_maxd = std::max(pc_maxd, std::fabs(ats[a].w - ref[b].w));
          ++a;
          ++b;
        }
      }
    }
    if (!pr_partner_ml && pr_pmlcheck && pr_pmax == 0 && n <= 40 && pr_cut <= 0.0 &&
        pr_veps <= 0.0) {
      std::map<std::tuple<int, int, int, int>, double> en;
      std::function<void(int, int, int, int, bool, double)> dfs = [&](int x0, int y0, int x,
                                                                     int y, bool inner, double w) {
        if (inner ? single_stems : E.end_ok(x, y)) {
          auto key = std::make_tuple(x, x0, y0, y);
          auto it = en.find(key);
          const double v = w + E.at(x, y);
          if (it == en.end()) en.emplace(key, v);
          else it->second = log_add(it->second, v);
        }
        for (const auto& [xo, yo] : pairs) {
          if (!(xo < x && yo > y)) continue;
          const double s = E.stem_interior(xo, x, y, yo);
          if (s > HALF_NEG) dfs(x0, y0, xo, yo, false, w + (s + pair_bonus));
        }
      };
      for (const auto& [x, y] : pairs)
        if (E.end_ok(x, y)) dfs(x, y, x, y, true, (E.at(x, y) + E.w_stem) + pair_bonus);
      std::vector<AT> ref;
      for (const auto& [k, w] : en)
        ref.push_back({std::get<0>(k), std::get<1>(k), std::get<2>(k), std::get<3>(k), w});
      std::sort(ref.begin(), ref.end(), at_less);
      size_t a = 0, b = 0;
      while (a < ats.size() || b < ref.size()) {
        if (b == ref.size() || (a < ats.size() && at_less(ats[a], ref[b]))) {
          ++pm_key_mis;
          ++a;
        } else if (a == ats.size() || at_less(ref[b], ats[a])) {
          ++pm_key_mis;
          ++b;
        } else {
          ++pm_n;
          pm_maxd = std::max(pm_maxd, std::fabs(ats[a].w - ref[b].w));
          ++a;
          ++b;
        }
      }
    }
    // pv: the previous entry in at_less order (null at the first)
    auto add = [&](const AT& a, const AT* pv) {
      const bool new_arm = !pv || a.xs != pv->xs || a.xp != pv->xp;
      if (new_arm) {
        prArmXs.push_back(a.xs);
        prArmXp.push_back(a.xp);
        prArmE0.push_back(static_cast<int>(pr_vxs ? prYy.size() : prEy.size()));
        prArmY0.push_back(static_cast<int>(prYy.size()));
      }
      if (!pr_vxs) {
        prEy.push_back(a.yp);
        prEys.push_back(a.ys);
        prEw.push_back(a.w);
      }
      const double yw = pr_vxs ? a.w : a.w + (E.w_beta1 + st.prS[a.ys + 1]);
      if (new_arm || a.yp != pv->yp) {
        prYy.push_back(a.yp);
        prYw.push_back(yw);
      } else {
        prYw.back() = log_add(prYw.back(), yw);
      }
    };
    for (size_t k = 0; k < ats.size(); ++k) add(ats[k], k == 0 ? nullptr : &ats[k - 1]);
    if (pr_vstream) {
      if (pr_vxs) {
        prYy.reserve(static_cast<size_t>(bk_n));
        prYw.reserve(static_cast<size_t>(bk_n));
      } else {
        prEy.reserve(static_cast<size_t>(bk_n));
        prEys.reserve(static_cast<size_t>(bk_n));
        prEw.reserve(static_cast<size_t>(bk_n));
      }
      for (int xs = 0; xs < n; ++xs) {
        std::vector<AT>& v = bk[xs];
        std::sort(v.begin(), v.end(), at_less);
        for (size_t k = 0; k < v.size(); ++k) add(v[k], k == 0 ? nullptr : &v[k - 1]);
        std::vector<AT>().swap(v);
        if ((xs & 31) == 31) mtrim();
      }
      std::vector<std::vector<AT>>().swap(bk);
      mtrim();
    }
    prA = static_cast<int>(prArmXs.size());
    prArmE0.push_back(static_cast<int>(pr_vxs ? prYy.size() : prEy.size()));
    prArmY0.push_back(static_cast<int>(prYy.size()));
    prXsFirst.assign(n + 1, prA);
    for (int r = prA - 1; r >= 0; --r) prXsFirst[prArmXs[r]] = r;
    for (int x = n - 1; x >= 0; --x) prXsFirst[x] = std::min(prXsFirst[x], prXsFirst[x + 1]);
    prByXp.resize(prA);
    for (int r = 0; r < prA; ++r) prByXp[r] = r;
    std::stable_sort(prByXp.begin(), prByXp.end(),
                     [&](int a, int b) { return prArmXp[a] < prArmXp[b]; });
    prXpFirst.assign(n + 1, prA);
    for (int r = prA - 1; r >= 0; --r) prXpFirst[prArmXp[prByXp[r]]] = r;
    for (int x = n - 1; x >= 0; --x) prXpFirst[x] = std::min(prXpFirst[x], prXpFirst[x + 1]);
    pr_arms = prA;
    pr_arm_entries = static_cast<long long>(pr_vxs ? prYy.size() : prEy.size());
    pr_yw_from_sweep();
    if (k_bound) k_bound_check();
    if (g_mem_log) {
      std::vector<MemItem> it{{"V_ats", static_cast<long long>(ats.size()), msz::vec(ats)}};
      mem_items(it, "", false, true);
      it.resize(2);
      mem_log("V_arms_built", it);
    }
    if (pr_vstream || g_mtrim) {
      std::vector<AT>().swap(ats);
      mtrim();
    }
    if (!pr_vdiag.empty() || pr_vsample_p > 0.0) pr_vdiag_static(pairs);
    if (g_static_only) {
      const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - ts.t0).count();
      std::fprintf(stderr,
                   "[vstatic] n %d pairs %zu V_entries %lld arm_entries %lld arms %d adds %lld "
                   "node_terms %lld hv_built %lld hv_entries %lld hv_terms %lld fut_built %lld "
                   "fut_int %lld fut_ml %lld v_before %lld sec %.6f\n",
                   n, pairs.size(), pr_V_entries, pr_arm_entries, prA, pr_adds, pr_nterms[1],
                   pr_hbuilt[1], pr_hentries[1], pr_hterms[1], st_nodes[1], st_kind[1][0],
                   st_kind[1][1], v_before, sec);
      std::fflush(nullptr);
      std::_Exit(0);
    }
    if (p2_sweep || p2_vsweep) {
      long long ent = 0;
      if (stem_)
        for (const auto& kv : stem_->prGamma) ent += static_cast<long long>(kv.second.size());
      std::fprintf(stderr,
                   "[sweep] pairs %zu gamma %lld g_states %lld v_states %lld pushes %lld "
                   "cap %d qm_max %d sec %.3f\n",
                   pairs.size(), ent, stem_ ? stem_->pr_sw_states_g : 0,
                   stem_ ? stem_->pr_sw_states_v : 0, stem_ ? stem_->pr_sw_pushes : 0, g_ntr_K,
                   stem_ ? stem_->pr_mbar_max : 0, stem_ ? stem_->pr_sw_sec : 0.0);
    }
  }

