// Pass 2, G through the multiloop states (--g-sweep states, paper default): sw_hyper_push
// pushes G entries (p', q') through Qm / Qm2 / pair states (edges E1-E9) of the nested
// pre-pass, cutting each state to its top b_G entries by C_pre(p'-1) G S_pre(q'+1)
// (cap_ga). pass2_apply rewrites Q^z_pre / pair lists / growth ranking when pass-2 flags ask.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  void sw_hyper_push(StemEst& st) {
    using Ga = std::unordered_map<uint64_t, double>;
    struct M {
      Ga g, v, b;
      bool empty() const { return g.empty() && v.empty() && b.empty(); }
    };
    // Nested pre-pass has l2r_mode==0 (scan's --pk-l2r does not apply here). Multi edges E3/E4
    // multi fire on that path. --g-sweep states turns G on beside V.
    const bool do_g = (g_sweep == G_SWEEP_STATES);
    const bool do_v = p2_vsweep;
    const bool do_b = sw_check && p2_vsweep && p2_v_cap <= 0 && g_ntr_K <= 0;
    const int Kg = g_ntr_K, Kv = p2_v_cap;
    const double a1 = E.w_alpha1, a2 = E.w_alpha2, a3 = E.w_alpha3;
    long long pushes = 0, st_g = 0, st_v = 0;
    int qm_max = 0;
    struct Hist {
      std::vector<uint32_t> c;
      long long n = 0, sum = 0;
      void add(size_t s) {
        if (c.size() <= s) c.resize(s + 1);
        ++c[s];
        ++n;
        sum += static_cast<long long>(s);
      }
      double mean() const { return n ? static_cast<double>(sum) / static_cast<double>(n) : 0.0; }
      long long max() const { return n && !c.empty() ? static_cast<long long>(c.size() - 1) : 0; }
      long long p99() const {
        if (!n) return 0;
        long long target = (n * 99 + 99) / 100;
        if (target < 1) target = 1;
        long long acc = 0;
        for (size_t i = 0; i < c.size(); ++i) {
          acc += c[i];
          if (acc >= target) return static_cast<long long>(i);
        }
        return c.empty() ? 0 : static_cast<long long>(c.size() - 1);
      }
    } g_at, v_at, g_pre, v_pre, g_kept, v_kept;
    auto scale_into = [&](M& d, const M& s, double a) {
      if (a < HALF_NEG || s.empty()) return;
      auto one = [&](Ga& dd, const Ga& ss, double& sec, long long* cnt, Hist* hist) {
        if (ss.empty()) return;
        if (hist) hist->add(ss.size());
        const auto t0 = std::chrono::steady_clock::now();
        for (const auto& [k, v] : ss) {
          if (v < HALF_NEG) continue;
          const double w = a + v;
          auto it = dd.find(k);
          if (it == dd.end()) dd.emplace(k, w);
          else it->second = log_add(it->second, w);
        }
        sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (cnt) ++*cnt;
      };
      if (do_g) one(d.g, s.g, g_merge_sec, &g_pushes, &g_at);
      if (do_v) one(d.v, s.v, v_merge_sec, &v_pushes, &v_at);
      if (do_b) {
        double ignore = 0;
        one(d.b, s.b, ignore, nullptr, nullptr);
      }
    };
    auto push_to = [&](M& d, const M& s, double a) {
      if (a < HALF_NEG || s.empty()) return;
      scale_into(d, s, a);
      ++pushes;
    };
    // State being capped: kind 0 pair, 1 Qm, 2 Qm2; (cap_i, cap_j) its span.
    int cap_kind = 0, cap_i = -1, cap_j = -1;
    const uint64_t trace_k = g_trace_p >= 0 ? key2(g_trace_p, g_trace_q) : ~0ULL;
    long long tr_seen[3] = {}, tr_drop[3] = {}, tr_lines = 0;
    auto cap_ga = [&](Ga& m) {
      const long long pre = static_cast<long long>(m.size());
      g_pre.add(m.size());
      ++g_pre_n;
      g_pre_sum += pre;
      if (pre > g_pre_max) g_pre_max = pre;
      const auto tc0 = std::chrono::steady_clock::now();
      const bool cut_here = !(g_cut_pair_only && cap_kind != 0);
      // --g-sweep states: always entries cut by C_pre * G * S_pre (bG = Kg).
      if (cut_here && Kg > 0 && static_cast<int>(m.size()) > Kg) {
        struct E {
          uint64_t k;
          double w, s;
        };
        std::vector<E> v;
        v.reserve(m.size());
        for (const auto& [k, w] : m) {
          if (w < HALF_NEG) continue;
          const int p = k_hi(k), q = k_lo(k);
          double s = w;
          if (p > 0 && !st.prC.empty()) s += st.prC[static_cast<size_t>(p - 1)];
          if (q + 1 < 0 || q + 1 >= static_cast<int>(st.prS.size())) continue;
          const double sb = st.prS[static_cast<size_t>(q + 1)];
          if (sb < HALF_NEG) continue;
          s += sb;
          v.push_back({k, w, s});
        }
        const bool on_path = g_trace_x < 0 || (cap_i >= g_trace_p && cap_i <= g_trace_x &&
                                               cap_j >= g_trace_y && cap_j <= g_trace_q);
        if (trace_k != ~0ULL && on_path) {
          const E* te = nullptr;
          for (const E& e : v)
            if (e.k == trace_k) te = &e;
          if (te) {
            ++tr_seen[cap_kind];
            long long rank = 0;
            for (const E& e : v)
              if (e.s > te->s || (e.s == te->s && e.k < te->k)) ++rank;
            if (rank >= Kg) {
              ++tr_drop[cap_kind];
              if (tr_lines++ < 40) {
                const double kth = [&]() {
                  std::vector<double> ss;
                  ss.reserve(v.size());
                  for (const E& e : v) ss.push_back(e.s);
                  std::nth_element(ss.begin(), ss.begin() + (Kg - 1), ss.end(), std::greater<double>());
                  return ss[static_cast<size_t>(Kg - 1)];
                }();
                std::fprintf(stderr,
                             "[g-trace] drop %s %d.%d rank %lld of %zu (bG %d) score %.6f kth %.6f\n",
                             cap_kind == 0 ? "pair" : cap_kind == 1 ? "Qm" : "Qm2", cap_i, cap_j,
                             rank, v.size(), Kg, te->s, kth);
              }
            }
          }
        }
        if (static_cast<int>(v.size()) > Kg) {
          std::nth_element(v.begin(), v.begin() + Kg, v.end(), [](const E& a, const E& b) {
            if (a.s != b.s) return a.s > b.s;
            return a.k < b.k;
          });
          v.resize(static_cast<size_t>(Kg));
        }
        m.clear();
        for (const E& e : v) m.emplace(e.k, e.w);
      } else if (trace_k != ~0ULL && m.count(trace_k) &&
                 (g_trace_x < 0 || (cap_i >= g_trace_p && cap_i <= g_trace_x && cap_j >= g_trace_y &&
                                    cap_j <= g_trace_q))) {
        ++tr_seen[cap_kind];
      }
      g_cut_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - tc0).count();
      const long long post = static_cast<long long>(m.size());
      ++g_post_n;
      g_post_sum += post;
      if (post > g_post_max) g_post_max = post;
    };
    auto cap_mu = [&](Ga& m) {
      const long long pre = static_cast<long long>(m.size());
      v_pre.add(m.size());
      ++v_pre_n;
      v_pre_sum += pre;
      if (pre > v_pre_max) v_pre_max = pre;
      const auto tc0 = std::chrono::steady_clock::now();
      if (Kv > 0 && static_cast<int>(m.size()) > Kv) {
        struct E {
          uint64_t k;
          double w, s;
        };
        std::vector<E> v;
        v.reserve(m.size());
        for (const auto& [k, w] : m) {
          if (w < HALF_NEG) continue;
          const int b = k_hi(k);
          double s = w;
          if (b > 0 && !st.prC.empty()) s += st.prC[static_cast<size_t>(b - 1)];
          v.push_back({k, w, s});
        }
        if (static_cast<int>(v.size()) > Kv) {
          std::nth_element(v.begin(), v.begin() + Kv, v.end(), [](const E& a, const E& b) {
            if (a.s != b.s) return a.s > b.s;
            return a.k < b.k;
          });
          v.resize(static_cast<size_t>(Kv));
        }
        m.clear();
        for (const E& e : v) m.emplace(e.k, e.w);
      }
      v_cut_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - tc0).count();
      const long long post = static_cast<long long>(m.size());
      ++v_post_n;
      v_post_sum += post;
      if (post > v_post_max) v_post_max = post;
    };
    auto cap_m = [&](M& s) {
      if (do_g && !s.g.empty()) {
        cap_ga(s.g);
        if (!s.g.empty()) {
          ++st_g;
          g_kept.add(s.g.size());
        }
      }
      if (do_v && !s.v.empty()) {
        cap_mu(s.v);
        if (!s.v.empty()) {
          ++st_v;
          v_kept.add(s.v.size());
        }
      }
    };
    auto add1 = [&](Ga& d, uint64_t k, double w) {
      if (w < HALF_NEG) return;
      auto it = d.find(k);
      if (it == d.end()) d.emplace(k, w);
      else it->second = log_add(it->second, w);
    };
    // A pruned target is not zeroed by the edge coefficient: dQ_s/dQ_c does not contain Q_c.
    long long hk_map[3] = {}, hk_bad[3] = {}, hk_e[12] = {}, hk_eb[12] = {};
    auto kept_at = [&](const OMap& at, int i) { return i >= 0 && at.has(i); };
    auto note_push = [&](int e, bool ok) {
      ++hk_e[e];
      if (!ok) ++hk_eb[e];
    };
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<double> Cb(static_cast<size_t>(n), NEG_INF);
    if (n > 0) Cb[static_cast<size_t>(n - 1)] = 0.0;
    double bc1 = NEG_INF, bc2 = NEG_INF;
    std::vector<std::unordered_map<int, M>> qm(static_cast<size_t>(n)), qm2(static_cast<size_t>(n));
    std::vector<std::unordered_map<int, M>> u1(static_cast<size_t>(n)), u2(static_cast<size_t>(n));
    std::vector<std::unordered_map<int, M>> u12(static_cast<size_t>(n)), u22(static_cast<size_t>(n));
    std::vector<std::unordered_map<int, M>> lpm(static_cast<size_t>(n));
    std::unordered_map<uint64_t, M> pairs;
    auto run_mix = [&](M& closed, const M* u1n, const M* u2n, int j, double w, M& ou1, M& ou2) {
      M p = closed;
      if (u1n) scale_into(p, *u1n, w);
      M a, b;
      scale_into(a, closed, E.s1w[j]);
      if (u2n) scale_into(a, *u2n, E.d5w[j] + w);
      scale_into(b, closed, E.d3w[j]);
      if (u2n) scale_into(b, *u2n, w);
      if ((u1n && !u1n->empty()) || (u2n && !u2n->empty())) ++pushes;
      closed = std::move(p);
      ou1 = std::move(a);
      ou2 = std::move(b);
    };
    for (int j = n - 1; j >= 0; --j) {
      qm_max = std::max(qm_max, static_cast<int>(Qm_at[j].live));
      if (dz) {
        double bP, b1, b2;
        run_back(Cb[static_cast<size_t>(j)], bc1, bc2, j, 0.0, bP, b1, b2);
        Cb[static_cast<size_t>(j)] = bP;
        bc1 = b1;
        bc2 = b2;
      } else if (j + 1 < n) {
        Cb[static_cast<size_t>(j)] = log_add(Cb[static_cast<size_t>(j)], Cb[static_cast<size_t>(j + 1)]);
      }
      auto extend_col = [&](int edge, std::unordered_map<int, M>& col, const OMap& at,
                            std::unordered_map<int, M>& U1, std::unordered_map<int, M>& U2,
                            std::unordered_map<int, M>* U1n, std::unordered_map<int, M>* U2n) {
        std::unordered_map<int, M> nu1, nu2;
        for (size_t r = 0; r < at.slots(); ++r) {
          const int i = at.keys[r];
          if (i < 0) continue;
          note_push(edge, kept_at(at, i));
          const M* a = nullptr;
          const M* b = nullptr;
          if (U1n) {
            auto it = U1n->find(i);
            if (it != U1n->end()) a = &it->second;
          }
          if (U2n) {
            auto it = U2n->find(i);
            if (it != U2n->end()) b = &it->second;
          }
          M& closed = col[i];
          if (closed.empty() && !a && !b) continue;
          run_mix(closed, a, b, j, a3, nu1[i], nu2[i]);
        }
        U1.swap(nu1);
        U2.swap(nu2);
      };
      if (dz) {
        extend_col(0, qm[static_cast<size_t>(j)], Qm_at[j], u1[static_cast<size_t>(j)],
                   u2[static_cast<size_t>(j)], j + 1 < n ? &u1[static_cast<size_t>(j + 1)] : nullptr,
                   j + 1 < n ? &u2[static_cast<size_t>(j + 1)] : nullptr);
        extend_col(1, qm2[static_cast<size_t>(j)], Qm2_at[j], u12[static_cast<size_t>(j)],
                   u22[static_cast<size_t>(j)], j + 1 < n ? &u12[static_cast<size_t>(j + 1)] : nullptr,
                   j + 1 < n ? &u22[static_cast<size_t>(j + 1)] : nullptr);
      } else if (j + 1 < n) {
        for (size_t r = 0; r < Qm_at[j].slots(); ++r) {
          const int i = Qm_at[j].keys[r];
          if (i < 0) continue;
          auto it = qm[static_cast<size_t>(j + 1)].find(i);
          if (it != qm[static_cast<size_t>(j + 1)].end()) {
            note_push(0, kept_at(Qm_at[j], i));
            push_to(qm[static_cast<size_t>(j)][i], it->second, a3);
          }
        }
        for (size_t r = 0; r < Qm2_at[j].slots(); ++r) {
          const int i = Qm2_at[j].keys[r];
          if (i < 0) continue;
          auto it = qm2[static_cast<size_t>(j + 1)].find(i);
          if (it != qm2[static_cast<size_t>(j + 1)].end()) {
            note_push(1, kept_at(Qm2_at[j], i));
            push_to(qm2[static_cast<size_t>(j)][i], it->second, a3);
          }
        }
      }
      if (l2r_mode == 0 && !lpP.empty() && !lpM.empty()) {
        const int L = lead_span();
        const OMap& src = lpP[static_cast<size_t>(j)];
        for (size_t s = 0; s < src.slots(); ++s) {
          const int ib = src.keys[s];
          if (ib < 0 || !Qm2_at[j].has(ib)) continue;
          const int p_lo = std::max(0, ib - L);
          for (int p = ib - 1; p >= p_lo; --p) {
            const int q = next_mate(p, j);
            if (q < 0 || q >= n) continue;
            const int k = p * 3 + lp_kind(q - j - 1);
            if (!lpM[static_cast<size_t>(q)].has(k)) continue;
            auto it = lpm[static_cast<size_t>(q)].find(k);
            if (it == lpm[static_cast<size_t>(q)].end() || it->second.empty()) continue;
            const double factor = unpaired_gap(p + 1, ib - 1) + unpaired_gap(j + 1, q - 1);
            note_push(10, true);
            push_to(qm2[static_cast<size_t>(j)][ib], it->second, factor);
          }
        }
      }
      const int qhi = std::min(n, j + 2 + cap);
      for (size_t r = 0; r < Qb_at[j].slots(); ++r) {
        const int d = Qb_at[j].keys[r];
        if (d < 0) continue;
        M& inner = pairs[key2(d, j)];
        if (l2r_mode == 0) {
          const int p_lo = std::max(0, d - 1 - cap);
          for (int p = d - 1; p >= p_lo; --p) {
            const int budget = cap - (d - p - 1);
            if (budget < 0) break;
            for (int q = next_mate(p, j); q != -1 && q - j - 1 <= budget; q = next_mate(p, q)) {
              if (!Qb_at[q].has(p)) continue;
              const double w = E.interior(p, d, j, q);
              if (w < HALF_NEG) continue;
              auto it = pairs.find(key2(p, q));
              note_push(2, it != pairs.end());
              if (it != pairs.end()) push_to(inner, it->second, w);
            }
          }
          continue;
        }
        for (int q = j + 1; q < qhi; ++q) {
          for (size_t s = 0; s < Qb_at[q].slots(); ++s) {
            const int i = Qb_at[q].keys[s];
            if (i < 0 || d <= i) continue;
            if ((d - i - 1) + (q - j - 1) > cap) continue;
            const double w = E.interior(i, d, j, q);
            if (w < HALF_NEG) continue;
            auto it = pairs.find(key2(i, q));
            if (it != pairs.end()) {
              note_push(2, kept_at(Qb_at[j], d) && kept_at(Qb_at[q], i));
              push_to(inner, it->second, w);
            }
          }
        }
      }
      for (size_t r = 0; r < Qb_at[j].slots(); ++r) {
        const int i = Qb_at[j].keys[r];
        if (i < 0 || !E.end_ok(i, j)) continue;
        note_push(3, kept_at(Qb_at[j], i));
        const double at = E.at(i, j);
        const double cb = Cb[static_cast<size_t>(j)];
        if (do_b && cb > HALF_NEG && at > HALF_NEG) add1(pairs[key2(i, j)].b, key2(i, 0), cb + at);
        const double left = i == 0 ? 0.0 : C[static_cast<size_t>(i - 1)];
        const double qb = Qb[ix(i, j)];
        if (i > 0 && cb > HALF_NEG && qb > HALF_NEG)
          Cb[static_cast<size_t>(i - 1)] = log_add(Cb[static_cast<size_t>(i - 1)], (cb + qb) + at);
        (void)left;
      }
      cap_j = j;
      cap_kind = 1;
      for (auto& kv : qm[static_cast<size_t>(j)]) {
        cap_i = kv.first;
        cap_m(kv.second);
      }
      cap_kind = 2;
      for (auto& kv : qm2[static_cast<size_t>(j)]) {
        cap_i = kv.first;
        cap_m(kv.second);
      }
      cap_kind = 0;
      for (size_t r = 0; r < Qb_at[j].slots(); ++r) {
        const int d = Qb_at[j].keys[r];
        if (d < 0 || !E.end_ok(d, j)) continue;
        const double w_at = E.at(d, j);
        const double w2 = a2 + w_at;
        const double obj = Qb[ix(d, j)];
        M& dst = pairs[key2(d, j)];
        for (size_t s = 0; s < Qm_at[j].slots(); ++s) {
          const int i = Qm_at[j].keys[s];
          if (i < 0 || i > d || d - i > E.ml_left || !lead_gap_ok(d - i)) continue;
          const LeadEnt* lead = nullptr;
          if (l2r_mode != 0 && d < static_cast<int>(qlead.size())) {
            const auto& col = qlead[static_cast<size_t>(d)];
            auto it = std::lower_bound(col.begin(), col.end(), i,
                                       [](const LeadEnt& e, int x) { return e.x0 < x; });
            if (it != col.end() && it->x0 == i) lead = &*it;
          }
          // Pass 1's M = P has no left unpaired run; the fragment starts on the pair.
          // --pk-l2r off is the DP03 reference and keeps the leading-run window.
          if (l2r_mode == 0 && !allow_pk) {
            if (i != d) continue;
          } else if (!lead) {
            continue;
          }
          const double w = (l2r_mode == 0 && !allow_pk) ? w2 : w2 + lead->w;
          auto it = qm[static_cast<size_t>(j)].find(i);
          if (it != qm[static_cast<size_t>(j)].end()) {
            note_push(4, kept_at(Qb_at[j], d) && kept_at(Qm_at[j], i));
            push_to(dst, it->second, w);
          }
        }
        if (d == 0) continue;
        for (size_t s = 0; s < Qm_at[d - 1].slots(); ++s) {
          const int i = Qm_at[d - 1].keys[s];
          if (i < 0) continue;
          const double qm_in = Qm[ix(i, d - 1)];
          if (qm_in < HALF_NEG) continue;
          auto it = qm[static_cast<size_t>(j)].find(i);
          auto it2 = qm2[static_cast<size_t>(j)].find(i);
          if (it != qm[static_cast<size_t>(j)].end()) {
            note_push(5, kept_at(Qb_at[j], d) && kept_at(Qm_at[j], i) && kept_at(Qm_at[d - 1], i));
            push_to(dst, it->second, qm_in + w2);
            if (obj > HALF_NEG) {
              note_push(6, kept_at(Qm_at[d - 1], i) && kept_at(Qm_at[j], i));
              push_to(qm[static_cast<size_t>(d - 1)][i], it->second, obj + w2);
            }
          }
          if (it2 != qm2[static_cast<size_t>(j)].end()) {
            note_push(7, kept_at(Qb_at[j], d) && kept_at(Qm2_at[j], i) && kept_at(Qm_at[d - 1], i));
            push_to(dst, it2->second, qm_in + w2);
            if (obj > HALF_NEG) {
              note_push(8, kept_at(Qm_at[d - 1], i) && kept_at(Qm2_at[j], i));
              push_to(qm[static_cast<size_t>(d - 1)][i], it2->second, obj + w2);
            }
          }
        }
      }
      const double w_ml = a1 + a2;
      for (size_t r = 0; r < Qb_at[j].slots(); ++r) {
        const int i = Qb_at[j].keys[r];
        if (i < 0) continue;
        M& m = pairs[key2(i, j)];
        if (E.end_ok(i, j)) {
          // Exterior-loop pair factor. --out-id uses pass 1's E.at, which is also the production seed.
          const double e_out = E.at(i, j);
          if (do_g) add1(m.g, key2(i, j), e_out);
          if (do_v) {
            const double sb = st.prS[static_cast<size_t>(j + 1)];
            if (sb > HALF_NEG) add1(m.v, key2(i, 0), e_out + sb);
          }
        }
        cap_kind = 0;
        cap_i = i;
        cap_j = j;
        cap_m(m);
        if (l2r_mode != 0 && E.end_ok(i, j) && i + 1 <= j - 1 && Qm2_at[j - 1].has(i + 1)) {
          note_push(9, kept_at(Qb_at[j], i) && kept_at(Qm2_at[j - 1], i + 1));
          push_to(qm2[static_cast<size_t>(j - 1)][i + 1], m, w_ml + E.at_close(i, j));
        }
      }
      if (l2r_mode == 0 && !lpM.empty()) {
        const OMap& col = lpM[static_cast<size_t>(j)];
        for (size_t r = 0; r < col.slots(); ++r) {
          const int k = col.keys[r];
          if (k < 0) continue;
          const int p = k / 3;
          const int kind = k % 3;
          M& mu = lpm[static_cast<size_t>(j)][k];
          if (Qb_at[j].has(p) && E.end_ok(p, j)) {
            auto it = pairs.find(key2(p, j));
            if (it != pairs.end() && !it->second.empty()) {
              note_push(9, true);
              push_to(mu, it->second, E.w_alpha1 + E.w_alpha2 + E.at_close(p, j));
            }
          }
          int qn, nk;
          double delta;
          if (!lp_step(p, j, kind, qn, nk, delta)) continue;
          const int kn = p * 3 + nk;
          if (qn < 0 || qn >= n || !lpM[static_cast<size_t>(qn)].has(kn)) continue;
          auto it = lpm[static_cast<size_t>(qn)].find(kn);
          if (it == lpm[static_cast<size_t>(qn)].end() || it->second.empty()) continue;
          note_push(11, true);
          push_to(mu, it->second, delta);
        }
      }
      auto tally = [&](const std::unordered_map<int, M>& col, const OMap& at, int slot) {
        for (const auto& kv : col) {
          ++hk_map[slot];
          if (!at.has(kv.first)) ++hk_bad[slot];
        }
      };
      tally(qm[static_cast<size_t>(j)], Qm_at[j], 0);
      tally(qm2[static_cast<size_t>(j)], Qm2_at[j], 1);
      std::unordered_map<int, M>().swap(qm[static_cast<size_t>(j)]);
      std::unordered_map<int, M>().swap(qm2[static_cast<size_t>(j)]);
      if (j + 1 < n) {
        std::unordered_map<int, M>().swap(u1[static_cast<size_t>(j + 1)]);
        std::unordered_map<int, M>().swap(u2[static_cast<size_t>(j + 1)]);
        std::unordered_map<int, M>().swap(u12[static_cast<size_t>(j + 1)]);
        std::unordered_map<int, M>().swap(u22[static_cast<size_t>(j + 1)]);
      }
    }
    st.prGamma.clear();
    st.prMu.clear();
    st.prMuFloor.clear();
    for (auto& [k, m] : pairs) {
      const int px = k_hi(k), py = k_lo(k);
      ++hk_map[2];
      if (py < 0 || py >= n || !Qb_at[py].has(px)) ++hk_bad[2];
      if (do_g && !m.g.empty()) st.prGamma.emplace(k, std::move(m.g));
      if (do_v && !m.v.empty()) {
        double mn = NEG_INF;
        bool any = false;
        for (const auto& [bk, w] : m.v) {
          if (w < HALF_NEG) continue;
          if (!any || w < mn) {
            mn = w;
            any = true;
          }
        }
        if (any) st.prMuFloor.emplace(k, mn);
        st.prMu.emplace(k, std::move(m.v));
      }
    }
    if (!g_mu_dump.empty()) {
      std::ofstream md(g_mu_dump);
      for (const auto& [k, mv] : st.prMu) {
        const int x = k_hi(k), y = k_lo(k);
        for (const auto& [bk, w] : mv) {
          if (!(w > HALF_NEG)) continue;
          md << "V " << x << " " << y << " " << k_hi(bk) << " " << w << "\n";
        }
      }
      for (const auto& [k, mg] : st.prGamma) {
        const int x = k_hi(k), y = k_lo(k);
        std::unordered_map<int, double> byq;
        for (const auto& [bk, w] : mg) {
          if (!(w > HALF_NEG)) continue;
          const int p = k_hi(bk), q = k_lo(bk);
          md << "G " << x << " " << y << " " << p << " " << q << " " << w << "\n";
          const double g = w + prefix_log(p);
          auto it = byq.find(q);
          if (it == byq.end()) byq.emplace(q, g);
          else it->second = log_add(it->second, g);
        }
        for (const auto& [q, g] : byq) md << "Q " << x << " " << y << " " << q << " " << g << "\n";
      }
    }
    st.pr_sw_states_g = st_g;
    st.pr_sw_states_v = st_v;
    st.pr_sw_pushes = pushes;
    if (g_out_id) st.prQb_b = Qb_b.dense();
    auto out_rel = [](double a, double b) {
      if (a < HALF_NEG && b < HALF_NEG) return 0.0;
      if (a < HALF_NEG || b < HALF_NEG) return 1e300;
      return std::fabs(std::expm1(a - b));
    };
    auto out_id_one = [&](const char* tag, bool fold_s,
                          const std::unordered_map<uint64_t, std::unordered_map<uint64_t, double>>& maps) {
      if (!g_out_id) return;
      std::vector<double> err;
      err.reserve(static_cast<size_t>(n));
      int worst_i = -1, worst_j = -1;
      for (int j = 0; j < n; ++j) {
        for (size_t r = 0; r < Qb_at[j].slots(); ++r) {
          const int i = Qb_at[j].keys[r];
          if (i < 0 || !E.end_ok(i, j)) continue;
          const double ref = Qb_b[ix(i, j)];
          double got = NEG_INF;
          auto it = maps.find(key2(i, j));
          if (it != maps.end()) {
            for (const auto& [bk, w] : it->second) {
              if (w < HALF_NEG) continue;
              const int p = k_hi(bk);
              const int q = k_lo(bk);
              double t = w;
              if (p > 0) t += C[static_cast<size_t>(p - 1)];
              if (fold_s) {
                if (q + 1 < 0 || q + 1 >= static_cast<int>(st.prS.size())) continue;
                const double sb = st.prS[static_cast<size_t>(q + 1)];
                if (sb < HALF_NEG) continue;
                t += sb;
              }
              got = log_add(got, t);
            }
          }
          if (got < HALF_NEG && ref < HALF_NEG) continue;
          const double d = out_rel(got, ref);
          if (err.empty() || d > err.front()) {
            // track max separately
          }
          err.push_back(d);
          if (worst_i < 0 || d > out_rel(0, 0) && d >= err.back()) {
            worst_i = i;
            worst_j = j;
          }
        }
      }
      double mx = 0, med = 0;
      if (!err.empty()) {
        mx = *std::max_element(err.begin(), err.end());
        std::nth_element(err.begin(), err.begin() + static_cast<long>(err.size() / 2), err.end());
        med = err[err.size() / 2];
      }
      std::fprintf(stderr, "[out-id] %s pairs %zu max_rel %.3e med_rel %.3e\n", tag, err.size(), mx,
                   med);
      (void)worst_i;
      (void)worst_j;
    };
    if (g_out_id && do_v) out_id_one("V", false, st.prMu);
    if (g_out_id && do_g) out_id_one("G-states", true, st.prGamma);
    st.pr_mbar_max = std::max(st.pr_mbar_max, qm_max);
    st.pr_sw_sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (do_g && trace_k != ~0ULL) {
      int fin_n = 0;
      double fin_w = NEG_INF;
      for (const auto& [pk, g] : st.prGamma) {
        auto it = g.find(trace_k);
        if (it == g.end()) continue;
        ++fin_n;
        fin_w = log_add(fin_w, it->second);
      }
      std::fprintf(stderr,
                   "[g-trace] key %d.%d seen pair %lld Qm %lld Qm2 %lld dropped pair %lld Qm %lld "
                   "Qm2 %lld final_pairs_holding %d\n",
                   g_trace_p, g_trace_q, tr_seen[0], tr_seen[1], tr_seen[2], tr_drop[0], tr_drop[1],
                   tr_drop[2], fin_n);
      (void)fin_w;
      if (g_trace_x >= 0) {
        auto it = st.prGamma.find(key2(g_trace_x, g_trace_y));
        const bool has_root = it != st.prGamma.end();
        const bool has_key = has_root && it->second.count(trace_k);
        std::fprintf(stderr, "[g-trace] root %d.%d in_G %d entries %zu holds %d.%d %d w %.6f\n",
                     g_trace_x, g_trace_y, has_root ? 1 : 0, has_root ? it->second.size() : 0,
                     g_trace_p, g_trace_q, has_key ? 1 : 0,
                     has_key ? it->second.at(trace_k) : NEG_INF);
      }
    }
    std::fprintf(stderr,
                 "[sweep-kept] l2r_mode %d allow_pk %d maps qm %lld bad %lld qm2 %lld bad %lld "
                 "pair %lld bad %lld\n",
                 l2r_mode, allow_pk ? 1 : 0, hk_map[0], hk_bad[0], hk_map[1], hk_bad[1], hk_map[2],
                 hk_bad[2]);
    static const char* ename[12] = {"E1 qm-ext", "E2 qm2-ext", "E3 interior", "E4 exterior",
                                    "E5 branch1", "E6 branch-qm-pair", "E7 branch-qm-tail",
                                    "E8 branch-qm2-pair", "E8 branch-qm2-tail", "E9 close",
                                    "E3 multi", "E4 multi-ext"};
    for (int e = 0; e < 12; ++e)
      std::fprintf(stderr, "[sweep-kept] %s push %lld bad %lld\n", ename[e], hk_e[e], hk_eb[e]);
    // E3 multi / E4 multi-ext are the nested pre-pass Multi tables (l2r_mode==0). The scan's
    // --pk-l2r on does not change np.l2r_mode; those edges are expected here.
    if (do_v && p2_v_cap <= 0 && g_ntr_K <= 0 && !Qb_b.empty()) {
      auto rel = [](double a, double b) {
        if (a < HALF_NEG && b < HALF_NEG) return 0.0;
        if (a < HALF_NEG || b < HALF_NEG) return 1e300;
        return std::fabs(std::expm1(a - b));
      };
      double worst = 0, dmin = 1e300, dmax = -1e300, wgot = 0, wref = 0;
      int wx = -1, wy = -1, wn = 0;
      long long nent = 0, nbad = 0;
      for (int j = 0; j < n; ++j) {
        for (size_t r = 0; r < Qb_at[j].slots(); ++r) {
          const int i = Qb_at[j].keys[r];
          if (i < 0 || !E.end_ok(i, j)) continue;
          const double ref = Qb_b[ix(i, j)];
          double got = NEG_INF;
          int nk = 0;
          auto it = st.prMu.find(key2(i, j));
          if (it != st.prMu.end()) {
            for (const auto& [bk, w] : it->second) {
              if (w < HALF_NEG) continue;
              ++nk;
              const int b = k_hi(bk);
              const double c = b > 0 ? C[static_cast<size_t>(b - 1)] : 0.0;
              got = log_add(got, c + w);
            }
          }
          if (got < HALF_NEG && ref < HALF_NEG) continue;
          ++nent;
          const double d = rel(got, ref);
          const double dl = (got > HALF_NEG && ref > HALF_NEG) ? got - ref : 1e300;
          if (dl < dmin) dmin = dl;
          if (dl > dmax) dmax = dl;
          if (d > 1e-9) ++nbad;
          if (d > worst) {
            worst = d;
            wx = i;
            wy = j;
            wgot = got;
            wref = ref;
            wn = nk;
          }
        }
      }
      std::fprintf(stderr,
                   "[v-out] pairs %lld bad %lld max_rel %.3e at %d.%d got %.6f ref %.6f dlog %.6f "
                   "keys %d dlog_min %.6f dlog_max %.6f\n",
                   nent, nbad, worst, wx, wy, wgot, wref, wgot - wref, wn, dmin, dmax);
    }
    std::fprintf(stderr, "[sweep-h] g_states %lld v_states %lld pushes %lld qm_max %d sec %.3f\n",
                 st_g, st_v, pushes, qm_max, st.pr_sw_sec);
    std::fprintf(stderr,
                 "[map] v_maps %lld v_mean %.4f v_max %lld g_maps %lld g_mean %.4f g_max %lld "
                 "pushes %lld\n",
                 v_kept.n, v_kept.mean(), v_kept.max(), g_kept.n, g_kept.mean(), g_kept.max(),
                 pushes);
    const double usec = pushes > 0 ? st.pr_sw_sec * 1e6 / static_cast<double>(pushes) : 0.0;
    std::fprintf(stderr,
                 "[sweep-push] g_n %lld g_mean %.3f g_p99 %lld g_pre_mean %.3f "
                 "v_n %lld v_mean %.3f v_p99 %lld v_pre_mean %.3f usec %.4f\n",
                 g_at.n, g_at.mean(), g_at.p99(), g_pre.mean(), v_at.n, v_at.mean(), v_at.p99(),
                 v_pre.mean(), usec);
    if (do_b && !Qb_b.empty()) {
      auto rel = [](double a, double b) {
        if (a < HALF_NEG && b < HALF_NEG) return 0.0;
        if (a < HALF_NEG || b < HALF_NEG) return 1e300;
        return std::fabs(std::expm1(a - b));
      };
      double worst = 0;
      int wx = -1, wy = -1;
      long long nent = 0;
      for (int j = 0; j < n; ++j) {
        for (size_t r = 0; r < Qb_at[j].slots(); ++r) {
          const int i = Qb_at[j].keys[r];
          if (i < 0 || !E.end_ok(i, j)) continue;
          const double ref = Qb_b[ix(i, j)];
          double got = NEG_INF;
          auto it = pairs.find(key2(i, j));
          if (it != pairs.end()) {
            for (const auto& [bk, w] : it->second.b) {
              if (w < HALF_NEG) continue;
              const int b = k_hi(bk);
              const double c = b > 0 ? C[static_cast<size_t>(b - 1)] : 0.0;
              got = log_add(got, c + w);
            }
          }
          if (got < HALF_NEG && ref < HALF_NEG) continue;
          ++nent;
          const double d = rel(got, ref);
          if (d > worst) {
            worst = d;
            wx = i;
            wy = j;
          }
        }
      }
      std::fprintf(stderr, "[sweep-4c] pairs %lld max_rel %.3e at %d.%d\n", nent, worst, wx, wy);
    }
  }

  void pass2_apply(StemEst& st) {
    if (p2_v_cap > 0 || p2_g_eager || p2_qz_beam > 0 || p2_sweep || p2_vsweep) st.prC = C;
    if (sw_check && !Qb_b.empty()) st.prQb_b = Qb_b.dense();
    if (p2_sweep || p2_vsweep) sw_hyper_push(st);
    if (p2_qz_beam <= 0 && p2_v_mode == 0 && p2_g_mode == 0) return;
    if (st.prC.empty()) st.prC = C;
    auto take = [&](int mode, int bp, std::vector<std::vector<int>>& byL) {
      byL.assign(n, {});
      for (int j = 0; j < n; ++j) {
        const OMap& at = Qb_at[j];
        std::vector<std::pair<double, int>> v;
        v.reserve(at.slots());
        for (size_t r = 0; r < at.slots(); ++r) {
          const int i = at.keys[r];
          if (i < 0 || at.vals[r] < HALF_NEG) continue;
          // mode 2: inside Q^b. mode 3: C_pre(i-1) * Q^b (prefix * inside).
          double key = at.vals[r];
          if (mode == 3 && i > 0) key += C[static_cast<size_t>(i - 1)];
          v.emplace_back(key, i);
        }
        if ((mode == 2 || mode == 3) && bp > 0 && static_cast<int>(v.size()) > bp) {
          std::nth_element(v.begin(), v.begin() + bp, v.end(),
                           [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                             return a.first > b.first;
                           });
          v.resize(static_cast<size_t>(bp));
        }
        for (const auto& e : v) byL[e.second].push_back(j);
      }
      for (std::vector<int>& kr : byL) {
        std::sort(kr.begin(), kr.end());
        kr.erase(std::unique(kr.begin(), kr.end()), kr.end());
      }
    };
    if (p2_v_mode > 0) take(p2_v_mode, p2_v_bp, st.prKeptRightV);
    if (p2_g_mode > 0) {
      if (p2_v_mode == p2_g_mode && p2_v_bp == p2_g_bp)
        st.prKeptRightG = st.prKeptRightV;
      else
        take(p2_g_mode, p2_g_bp, st.prKeptRightG);
    }
    if (p2_qz_beam <= 0) return;
    // --qzbar-old: keep dense prQz from build_principled; do not replace with beamed qz_col.
    if (qzbar_old) {
      st.qz_on = 0;
      st.qz_col.clear();
      st.qz_col.shrink_to_fit();
      return;
    }
    if (qz_beam_fill && p2_qz_key != 2 && !qz_keep.empty()) {
      st.qz_col.assign(static_cast<size_t>(n), {});
      st.p2_qz_n = 0;
      for (int l = 0; l < n; ++l) {
        std::vector<StemEst::QzEnt>& col = st.qz_col[static_cast<size_t>(l)];
        col.reserve(qz_keep[static_cast<size_t>(l)].size());
        for (const QkEnt& e : qz_keep[static_cast<size_t>(l)])
          col.push_back({e.k, e.q, e.p, e.u1, e.u2});
        st.p2_qz_n += static_cast<long long>(col.size());
      }
      st.qz_on = 1;
      st.prQz.clear();
      st.prQz.shrink_to_fit();
      st.p2_qz_sec = 0.0;
      qz_keep.clear();
      qz_keep.shrink_to_fit();
      return;
    }
    if (p2_qz_key == 2) {
      st.prRight.assign(n, NEG_INF);
      for (int j = 0; j < n; ++j) {
        const OMap& at = Qb_at[j];
        double acc = NEG_INF;
        for (size_t r = 0; r < at.slots(); ++r) {
          const int i = at.keys[r];
          if (i < 0 || at.vals[r] < HALF_NEG) continue;
          const double b = Qb_b[ix(i, j)];
          if (b < HALF_NEG) continue;
          acc = log_add(acc, at.vals[r] + b - logZ);
        }
        st.prRight[j] = acc;
      }
    }
    const auto t0 = std::chrono::steady_clock::now();
    struct Acc {
      double acc, u1, u2;
    };
    std::vector<int> seen(static_cast<size_t>(n), 0), touch;
    std::vector<Acc> cell(static_cast<size_t>(n));
    int gen = 0;
    st.qz_col.assign(static_cast<size_t>(n), {});
    st.p2_qz_n = 0;
    auto cell_at = [&](int k) -> Acc& {
      if (seen[static_cast<size_t>(k)] != gen) {
        seen[static_cast<size_t>(k)] = gen;
        cell[static_cast<size_t>(k)] = {NEG_INF, NEG_INF, NEG_INF};
        touch.push_back(k);
      }
      return cell[static_cast<size_t>(k)];
    };
    for (int l = 0; l < n; ++l) {
      if (++gen == 0) {
        std::fill(seen.begin(), seen.end(), 0);
        gen = 1;
      }
      touch.clear();
      const double d5 = dz && l > 0 ? E.d5w[l - 1] : 0.0;
      const double s1 = dz ? E.s1w[l] : 0.0, d3 = dz ? E.d3w[l] : 0.0;
      if (dz) cell_at(l).u1 = E.wu;
      else cell_at(l).acc = E.wu;
      if (l > 0) {
        for (const StemEst::QzEnt& e : st.qz_col[static_cast<size_t>(l - 1)]) {
          Acc& a = cell_at(e.k);
          if (dz) {
            a.u1 = e.p > HALF_NEG ? e.p + E.wu : NEG_INF;
            const double t = log_add(e.u2, e.u1 > HALF_NEG ? e.u1 + d5 : NEG_INF);
            a.u2 = t > HALF_NEG ? t + E.wu : NEG_INF;
          } else {
            a.acc = e.q + E.wu;
          }
        }
      }
      const OMap& atb = Qb_at[l];
      for (size_t r = 0; r < atb.slots(); ++r) {
        const int d = atb.keys[r];
        if (d < 0 || !E.end_ok(d, l)) continue;
        const double w = atb.vals[r] + (E.w_nest + E.at(d, l));
        if (w < HALF_NEG) continue;
        Acc& self = cell_at(d);
        self.acc = log_add(self.acc, w);
        if (d > 0) {
          for (const StemEst::QzEnt& e : st.qz_col[static_cast<size_t>(d - 1)]) {
            if (e.q < HALF_NEG) continue;
            Acc& a = cell_at(e.k);
            a.acc = log_add(a.acc, e.q + w);
          }
        }
      }
      struct Cand {
        int k;
        double q, p, u1, u2, key;
      };
      std::vector<Cand> cs;
      cs.reserve(touch.size());
      for (int k : touch) {
        const Acc& a = cell[static_cast<size_t>(k)];
        const double q = dz ? run_close(a.acc, a.u1, a.u2, s1, d3) : a.acc;
        if (q < HALF_NEG) continue;
        double key = q + (k > 0 ? st.prC[static_cast<size_t>(k - 1)] : 0.0);
        // key2 = key1 * Pr(k-1 is a pair right end). Pr = 0 ranks last and is kept
        // only when the beam still has room (a missing Q^z_pre read stays 0).
        if (p2_qz_key == 2 && k > 0) {
          const double pr = st.prRight[static_cast<size_t>(k - 1)];
          key = pr > HALF_NEG ? key + pr : NEG_INF;
        }
        cs.push_back({k, q, dz ? a.acc : 0.0, a.u1, a.u2, key});
      }
      if (g_beam_stats.f) {
        std::vector<double> keys;
        keys.reserve(cs.size());
        for (const Cand& c : cs) keys.push_back(c.key);
        std::sort(keys.begin(), keys.end(), [](double a, double b) { return a > b; });
        const size_t kept_n =
            static_cast<int>(cs.size()) > p2_qz_beam ? static_cast<size_t>(p2_qz_beam) : keys.size();
        const int full = static_cast<int>(cs.size()) > p2_qz_beam ? 1 : 0;
        g_beam_stats.log_prune("Qz_pre", l, keys.size(), kept_n, BeamStats::mass_share(keys, kept_n),
                               full);
      }
      if (static_cast<int>(cs.size()) > p2_qz_beam) {
        std::nth_element(cs.begin(), cs.begin() + p2_qz_beam, cs.end(),
                         [](const Cand& a, const Cand& b) { return a.key > b.key; });
        cs.resize(static_cast<size_t>(p2_qz_beam));
      }
      std::sort(cs.begin(), cs.end(), [](const Cand& a, const Cand& b) { return a.k < b.k; });
      std::vector<StemEst::QzEnt>& col = st.qz_col[static_cast<size_t>(l)];
      col.reserve(cs.size());
      for (const Cand& c : cs) col.push_back({c.k, c.q, c.p, c.u1, c.u2});
      st.p2_qz_n += static_cast<long long>(col.size());
    }
    st.qz_on = 1;
    st.prQz.clear();
    st.prQz.shrink_to_fit();
    st.p2_qz_sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }

