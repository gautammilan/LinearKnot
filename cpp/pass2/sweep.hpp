// Pass 2, the linear sweep (--sweep on): sw_as_rn runs the V / G sweep over the nested
// pre-pass with the b_V and b_G cuts, and the --out-id identity check
// sum_{p',q'} C_pre(p'-1) G(p', q') S_pre(q'+1) = Qb_b.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  RN sw_as_rn(int x, int y, bool inner) {
    if (p2_sweep && l2r_mode != 0 && stem_ && !sw_g_ready) {
      sw_g_ready = true;
      if (g_sweep == G_SWEEP_DIRECT) {
        sw_filled = true;
        // ntrUL must exist before sw_fill's cap_ga (same L trunc key as racc_finish_trunc).
        if (g_ntr_K > 0 && !p2_g_suffix) {
          const auto t0 = std::chrono::steady_clock::now();
          ntr_init();
          sw_ntr_init_sec +=
              std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        }
        {
          const auto t0 = std::chrono::steady_clock::now();
          sw_fill(0, n);
          sw_fill(1, std::min(pr_loop, cap));
          sw_fill_sec +=
              std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        }
        // Entries per G node after fill (slot 0 = full right gap).
        for (const auto& kv : sw_ga[0]) {
          const long long e = static_cast<long long>(kv.second.size());
          ++sw_ga_nodes;
          sw_ga_ent_sum += e;
          if (e > sw_ga_ent_max) sw_ga_ent_max = e;
        }
      } else {
        // --g-sweep states: G already in stem_->prGamma from sw_hyper_push.
        // Fill time is the G merge+cut copied onto this DP from the nested pass.
        sw_filled = false;
        sw_fill_sec = g_merge_sec + g_cut_sec;
        if (stem_) {
          for (const auto& kv : stem_->prGamma) {
            const long long e = static_cast<long long>(kv.second.size());
            ++sw_ga_nodes;
            sw_ga_ent_sum += e;
            if (e > sw_ga_ent_max) sw_ga_ent_max = e;
          }
        }
        if (g_out_id && stem_ && !stem_->prQb_b.empty()) {
          // Side fill of the direct path so states and direct can be compared uncapped.
          sw_fill(0, n);
          auto rel = [](double a, double b) {
            if (a < HALF_NEG && b < HALF_NEG) return 0.0;
            if (a < HALF_NEG || b < HALF_NEG) return 1e300;
            return std::fabs(std::expm1(a - b));
          };
          auto sum_g = [&](const std::unordered_map<uint64_t, double>* src) {
            double got = NEG_INF;
            if (!src) return got;
            for (const auto& [bk, w] : *src) {
              if (w < HALF_NEG) continue;
              const int p = k_hi(bk), q = k_lo(bk);
              double t = w;
              if (p > 0 && !stem_->prC.empty()) t += stem_->prC[static_cast<size_t>(p - 1)];
              if (q + 1 < 0 || q + 1 >= static_cast<int>(stem_->prS.size())) continue;
              const double sb = stem_->prS[static_cast<size_t>(q + 1)];
              if (sb < HALF_NEG) continue;
              t += sb;
              got = log_add(got, t);
            }
            return got;
          };
          std::vector<double> err;
          const int nn = n * n;
          for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
              const double ref = stem_->prQb_b[static_cast<size_t>(i) * n + j];
              const auto it = sw_ga[0].find(key2(i, j));
              const std::unordered_map<uint64_t, double>* src = it == sw_ga[0].end() ? nullptr : &it->second;
              if ((src == nullptr || src->empty()) && ref < HALF_NEG) continue;
              if (ref < HALF_NEG && src == nullptr) continue;
              err.push_back(rel(sum_g(src), ref));
            }
          }
          (void)nn;
          double mx = 0, med = 0;
          if (!err.empty()) {
            mx = *std::max_element(err.begin(), err.end());
            std::nth_element(err.begin(), err.begin() + static_cast<long>(err.size() / 2), err.end());
            med = err[err.size() / 2];
          }
          std::fprintf(stderr, "[out-id] G-direct pairs %zu max_rel %.3e med_rel %.3e\n", err.size(),
                       mx, med);
          long long only_s = 0, only_d = 0, both = 0;
          double max_e = 0;
          std::unordered_set<uint64_t> seen;
          for (const auto& [pk, sg] : stem_->prGamma) {
            seen.insert(pk);
            auto it = sw_ga[0].find(pk);
            if (it == sw_ga[0].end()) {
              only_s += static_cast<long long>(sg.size());
              continue;
            }
            for (const auto& [k, w] : sg) {
              auto jt = it->second.find(k);
              if (jt == it->second.end()) ++only_s;
              else {
                ++both;
                max_e = std::max(max_e, rel(w, jt->second));
              }
            }
            for (const auto& kv : it->second)
              if (!sg.count(kv.first)) ++only_d;
          }
          for (const auto& [pk, dg] : sw_ga[0])
            if (!seen.count(pk)) only_d += static_cast<long long>(dg.size());
          std::fprintf(stderr,
                       "[g-entry] both %lld only_states %lld only_direct %lld max_rel %.3e\n", both,
                       only_s, only_d, max_e);
          {
            // Mismatched entries closest to their outer pair: (x - p') + (q' - y) smallest.
            struct Mm {
              int d, x, y, p, q;
              double s, dd;
            };
            std::vector<Mm> mm;
            const double mm_tol = std::getenv("GMM_TOL") ? std::atof(std::getenv("GMM_TOL")) : 1e-9;
            for (const auto& [pk, sg] : stem_->prGamma) {
              auto it = sw_ga[0].find(pk);
              if (it == sw_ga[0].end()) continue;
              const int x = k_hi(pk), y = k_lo(pk);
              for (const auto& [k, w] : sg) {
                auto jt = it->second.find(k);
                const double dv = jt == it->second.end() ? NEG_INF : jt->second;
                if (rel(w, dv) <= mm_tol) continue;
                const int p = k_hi(k), q = k_lo(k);
                mm.push_back({(x - p) + (q - y), x, y, p, q, w, dv});
              }
            }
            std::sort(mm.begin(), mm.end(), [](const Mm& a, const Mm& b) {
              if (a.d != b.d) return a.d < b.d;
              return a.s - a.dd > b.s - b.dd;
            });
            std::vector<long long> by_d(8, 0);
            for (const Mm& m : mm) ++by_d[static_cast<size_t>(std::min(m.d, 7))];
            std::fprintf(stderr, "[g-mm] mismatched %zu by outer distance 0..6,7+:", mm.size());
            for (long long c : by_d) std::fprintf(stderr, " %lld", c);
            std::fprintf(stderr, "\n");
            for (size_t t = 0; t < mm.size() && t < 12; ++t)
              std::fprintf(stderr,
                           "[g-mm] pair %d.%d entry %d.%d dist %d states %.6f direct %.6f end_ok %d/%d\n",
                           mm[t].x, mm[t].y, mm[t].p, mm[t].q, mm[t].d, mm[t].s, mm[t].dd,
                           E.end_ok(mm[t].x, mm[t].y) ? 1 : 0, E.end_ok(mm[t].p, mm[t].q) ? 1 : 0);
          }
          // Pairs direct holds that states does not: why.
          long long dp_n = 0, dp_notend = 0, dp_noref = 0, dp_sumneg = 0, dp_instates = 0;
          for (const auto& [pk, dg] : sw_ga[0]) {
            if (dg.empty()) continue;
            const int i = k_hi(pk), j = k_lo(pk);
            const double ref = stem_->prQb_b[static_cast<size_t>(i) * n + j];
            ++dp_n;
            if (stem_->prGamma.count(pk)) {
              ++dp_instates;
              continue;
            }
            if (!E.end_ok(i, j)) ++dp_notend;
            else if (ref < HALF_NEG) ++dp_noref;
            else if (sum_g(&dg) < HALF_NEG) ++dp_sumneg;
          }
          std::fprintf(stderr,
                       "[d-pairs] direct %lld in_states %lld extra_not_end_ok %lld extra_no_Qb_b %lld "
                       "extra_other %lld\n",
                       dp_n, dp_instates, dp_notend, dp_noref, dp_n - dp_instates - dp_notend - dp_noref);
          (void)dp_sumneg;
          // Entries states holds that direct does not: by the seed pair (p',q').
          long long so_seed_not_direct = 0, so_seed_not_end = 0, so_outer_missing = 0, so_other = 0;
          for (const auto& [pk, sg] : stem_->prGamma) {
            auto it = sw_ga[0].find(pk);
            for (const auto& [k, w] : sg) {
              if (it != sw_ga[0].end() && it->second.count(k)) continue;
              (void)w;
              const int p = k_hi(k), q = k_lo(k);
              if (it == sw_ga[0].end()) ++so_outer_missing;
              else if (!sw_ga[0].count(key2(p, q))) ++so_seed_not_direct;
              else if (!E.end_ok(p, q)) ++so_seed_not_end;
              else ++so_other;
            }
          }
          std::fprintf(stderr,
                       "[s-only] outer_pair_not_in_direct %lld seed_pair_not_in_direct %lld "
                       "seed_not_end_ok %lld reachable_seed_but_missing %lld\n",
                       so_outer_missing, so_seed_not_direct, so_seed_not_end, so_other);
          {
            // Per pair: share of states' outside mass in entries direct lacks.
            std::vector<double> frac;
            long long over1 = 0;
            for (const auto& [pk, sg] : stem_->prGamma) {
              const double all = sum_g(&sg);
              if (all < HALF_NEG) continue;
              auto it = sw_ga[0].find(pk);
              std::unordered_map<uint64_t, double> miss;
              for (const auto& [k, w] : sg)
                if (it == sw_ga[0].end() || !it->second.count(k)) miss.emplace(k, w);
              const double m = sum_g(&miss);
              const double f = m < HALF_NEG ? 0.0 : std::exp(m - all);
              frac.push_back(f);
              if (f > 0.01) ++over1;
            }
            double fm = 0, fmed = 0;
            if (!frac.empty()) {
              fm = *std::max_element(frac.begin(), frac.end());
              std::nth_element(frac.begin(), frac.begin() + static_cast<long>(frac.size() / 2), frac.end());
              fmed = frac[frac.size() / 2];
            }
            std::fprintf(stderr, "[s-only] mass share per pair: max %.4f med %.4g pairs>1%% %lld of %zu\n",
                         fm, fmed, over1, frac.size());
            // Largest missing entries by share of their pair's states mass.
            struct Mi {
              double f;
              int x, y, p, q;
            };
            std::vector<Mi> top;
            for (const auto& [pk, sg] : stem_->prGamma) {
              const double all = sum_g(&sg);
              if (all < HALF_NEG) continue;
              auto it = sw_ga[0].find(pk);
              for (const auto& [k, w] : sg) {
                if (it != sw_ga[0].end() && it->second.count(k)) continue;
                std::unordered_map<uint64_t, double> one{{k, w}};
                const double s1 = sum_g(&one);
                if (s1 < HALF_NEG) continue;
                top.push_back({std::exp(s1 - all), k_hi(pk), k_lo(pk), k_hi(k), k_lo(k)});
              }
            }
            std::sort(top.begin(), top.end(), [](const Mi& a, const Mi& b) { return a.f > b.f; });
            for (size_t t = 0; t < top.size() && t < 8; ++t)
              std::fprintf(stderr, "[s-only] top pair %d.%d entry %d.%d share %.4f end_ok %d/%d\n",
                           top[t].x, top[t].y, top[t].p, top[t].q, top[t].f,
                           E.end_ok(top[t].x, top[t].y) ? 1 : 0, E.end_ok(top[t].p, top[t].q) ? 1 : 0);
          }
          // Identity on the common pair set (states' 2948 pairs).
          std::vector<double> cerr;
          // Signed: log(direct sum) - log(Qb_b), and states' sum vs Qb_b, split by end_ok.
          long long c_up[2] = {}, c_dn[2] = {}, c_eq[2] = {}, c_dneg[2] = {}, c_rneg[2] = {};
          long long s_bad[2] = {};
          std::vector<double> dlog[2];
          for (const auto& [pk, sg] : stem_->prGamma) {
            const int i = k_hi(pk), j = k_lo(pk);
            const double ref = stem_->prQb_b[static_cast<size_t>(i) * n + j];
            auto it = sw_ga[0].find(pk);
            const double got = sum_g(it == sw_ga[0].end() ? nullptr : &it->second);
            cerr.push_back(rel(got, ref));
            const int e = E.end_ok(i, j) ? 1 : 0;
            if (rel(sum_g(&sg), ref) > 1e-9) ++s_bad[e];
            if (ref < HALF_NEG) ++c_rneg[e];
            else if (got < HALF_NEG) ++c_dneg[e];
            else {
              const double d = got - ref;
              dlog[e].push_back(d);
              if (d > 1e-9) ++c_up[e];
              else if (d < -1e-9) ++c_dn[e];
              else ++c_eq[e];
            }
          }
          for (int e = 1; e >= 0; --e) {
            double mn = 0, md = 0, mxd = 0;
            if (!dlog[e].empty()) {
              mn = *std::min_element(dlog[e].begin(), dlog[e].end());
              mxd = *std::max_element(dlog[e].begin(), dlog[e].end());
              std::nth_element(dlog[e].begin(), dlog[e].begin() + static_cast<long>(dlog[e].size() / 2),
                               dlog[e].end());
              md = dlog[e][dlog[e].size() / 2];
            }
            std::fprintf(stderr,
                         "[out-id-signed] %s states_bad %lld Qb_b_neg %lld direct_neg %lld "
                         "direct>Qb_b %lld direct<Qb_b %lld equal %lld dlog min %.4f med %.4f max %.4f\n",
                         e ? "end_ok" : "GU", s_bad[e], c_rneg[e], c_dneg[e], c_up[e], c_dn[e], c_eq[e],
                         mn, md, mxd);
          }
          if (!cerr.empty()) {
            const double cmx = *std::max_element(cerr.begin(), cerr.end());
            std::nth_element(cerr.begin(), cerr.begin() + static_cast<long>(cerr.size() / 2), cerr.end());
            std::fprintf(stderr, "[out-id] G-direct-common pairs %zu max_rel %.3e med_rel %.3e\n",
                         cerr.size(), cmx, cerr[cerr.size() / 2]);
          }
        }
      }
    }
    const auto t_read0 = std::chrono::steady_clock::now();
    RN r;
    std::unordered_map<int, double> byq;
    const std::unordered_map<uint64_t, double>* src = nullptr;
    std::unordered_map<uint64_t, std::unordered_map<uint64_t, double>>::const_iterator git;
    if (sw_filled) {
      auto it = sw_ga[inner ? 1 : 0].find(key2(x, y));
      if (it != sw_ga[inner ? 1 : 0].end()) src = &it->second;
    } else if (stem_ && !stem_->prGamma.empty()) {
      git = stem_->prGamma.find(key2(x, y));
      if (git != stem_->prGamma.end()) src = &git->second;
    } else {
      auto it = sw_ga[inner ? 1 : 0].find(key2(x, y));
      if (it != sw_ga[inner ? 1 : 0].end()) src = &it->second;
    }
    if (src) {
      sw_as_rn_ent_in += static_cast<long long>(src->size());
      for (const auto& [k, w] : *src) {
        const int q = k_lo(k);
        const double g = w + prefix_log(k_hi(k));
        auto jt = byq.find(q);
        if (jt == byq.end()) byq.emplace(q, g);
        else jt->second = log_add(jt->second, g);
      }
    }
    double tot = NEG_INF;
    for (const auto& [q, g] : byq) tot = log_add(tot, g);
    r.tot = tot;
    if (g_ntr_K > 0) {
      ntr_init();
      struct E {
        int i;
        double g, s;
      };
      std::vector<E> v;
      v.reserve(byq.size());
      for (const auto& [q, g] : byq) {
        double s = NEG_INF;
        if (p2_g_suffix) {
          const double sb = stem_->prS[static_cast<size_t>(q + 1)];
          if (sb > HALF_NEG) s = g + sb;
        } else if (q < static_cast<int>(ntrUL.size()) && ntrUL[q] > HALF_NEG) {
          s = g + ntrUL[q];
        }
        if (s > HALF_NEG) v.push_back({q, g, s});
      }
      const size_t K = static_cast<size_t>(g_ntr_K);
      if (v.size() > K) {
        std::nth_element(v.begin(), v.begin() + static_cast<long>(K), v.end(),
                         [](const E& a, const E& b) {
                           if (a.s != b.s) return a.s > b.s;
                           return a.i < b.i;
                         });
        v.resize(K);
      }
      std::sort(v.begin(), v.end(), [](const E& a, const E& b) { return a.i < b.i; });
      r.T[0].reserve(v.size());
      for (const E& e : v) r.T[0].emplace_back(e.i, e.g);
    } else if (!byq.empty()) {
      int lo = n, hi = -1;
      double mx = NEG_INF;
      for (const auto& [q, g] : byq) {
        lo = std::min(lo, q);
        hi = std::max(hi, q);
        mx = std::max(mx, g);
      }
      RV& a = r.L[0];
      a.lo = lo;
      a.mx = mx;
      a.v.assign(static_cast<size_t>(hi - lo + 1), 0.0);
      for (const auto& [q, g] : byq)
        if (g > HALF_NEG) a.v[static_cast<size_t>(q - lo)] = std::exp(g - mx);
    }
    ++sw_as_rn_calls;
    sw_as_rn_ent_out += static_cast<long long>(r.T[0].size());
    sw_as_rn_sec +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t_read0).count();
    return r;
  }
