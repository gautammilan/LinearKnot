// Dumps: JSON helpers, --gate-diag, --dump-states, --dump-scores, --dump-ext.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.


std::string json_num(double x) {
  char buf[40];
  std::snprintf(buf, sizeof buf, "%.17g", x);
  return buf;
}

template <class T>
std::string json_arr(const std::vector<T>& v) {
  std::string s = "[";
  for (size_t k = 0; k < v.size(); ++k) {
    if (k) s += ", ";
    if constexpr (std::is_floating_point<T>::value)
      s += json_num(v[k]);
    else
      s += std::to_string(v[k]);
  }
  return s + "]";
}

std::string json_str(const std::string& s) {
  std::string out = "\"";
  for (char c : s) {
    if (c == '"' || c == '\\') out += '\\';
    out += c;
  }
  return out + "\"";
}

// Per table (Qg, Qgl), over every stored state (run at infinite beam with the gate off for
// exact weights): w = exp(log inside + log outside - log Z), the state's posterior weight;
// "long" = outer pair p.q with q - p + 1 > long_span (a lower bound on the span of any
// pseudoknot the state is part of).  Per role (stem1: g = g1, both: g = max(g1, g2)) and
// tau: states passing g >= tau, their weight, their long weight, and how many of each
// (v, q) group's top-5 / top-20 states by w pass.
std::string DP03::gate_diag(const Gate& gd, const std::vector<double>& taus,
                            const std::string& id, int span, double delta, int long_span) {
  outside();
  const size_t T = taus.size();
  const int TOPS[2] = {5, 20};
  std::string js = "{\"id\": " + json_str(id) + ", \"n\": " + std::to_string(n) +
                   ", \"logZ\": " + json_num(logZ) + ", \"span\": " + std::to_string(span) +
                   ", \"delta\": " + json_num(delta) + ", \"long_span\": " +
                   std::to_string(long_span) + ", \"taus\": " + json_arr(taus) + ", \"tables\": {";
  const char* names[2] = {"Qg", "Qgl"};
  std::vector<Tab4>* tabs[2] = {&Qg, &Qgl};
  struct Item {
    double w, g[2];
    bool lg;
  };
  for (int ti = 0; ti < 2; ++ti) {
    long long states = 0, states_long = 0, groups = 0, top_total[2] = {0, 0};
    double w_all = 0.0, w_long = 0.0;
    std::vector<long long> pass_n[2], top_pass[2][2];
    std::vector<double> pass_w[2], pass_wl[2];
    for (int r = 0; r < 2; ++r) {
      pass_n[r].assign(T, 0);
      pass_w[r].assign(T, 0.0);
      pass_wl[r].assign(T, 0.0);
      top_pass[r][0].assign(T, 0);
      top_pass[r][1].assign(T, 0);
    }
    std::vector<Item> items;
    std::vector<size_t> order;
    for (int q = 0; q < n; ++q)
      for (const Grp& g : (*tabs[ti])[q].grps) {
        if (g.ents.empty()) continue;
        ++groups;
        items.clear();
        for (const Ent& e : g.ents) {
          Item it;
          it.w = (e.b < HALF_NEG || logZ < HALF_NEG) ? 0.0 : std::exp(e.a + e.b - logZ);
          double g1, g2;
          gd.roles(e.p, e.u, g.v, q, g1, g2);
          it.g[0] = g1;
          it.g[1] = std::max(g1, g2);
          it.lg = q - e.p + 1 > long_span;
          items.push_back(it);
          ++states;
          w_all += it.w;
          if (it.lg) {
            ++states_long;
            w_long += it.w;
          }
          for (int r = 0; r < 2; ++r)
            for (size_t k = 0; k < T; ++k)
              if (it.g[r] >= taus[k]) {
                ++pass_n[r][k];
                pass_w[r][k] += it.w;
                if (it.lg) pass_wl[r][k] += it.w;
              }
        }
        order.resize(items.size());
        for (size_t k = 0; k < order.size(); ++k) order[k] = k;
        std::stable_sort(order.begin(), order.end(),
                         [&](size_t a, size_t b) { return items[a].w > items[b].w; });
        for (int t = 0; t < 2; ++t) {
          size_t top = std::min(order.size(), static_cast<size_t>(TOPS[t]));
          top_total[t] += top;
          for (size_t m = 0; m < top; ++m)
            for (int r = 0; r < 2; ++r)
              for (size_t k = 0; k < T; ++k)
                if (items[order[m]].g[r] >= taus[k]) ++top_pass[r][t][k];
        }
      }
    if (ti) js += ", ";
    js += json_str(names[ti]) + ": {\"states\": " + std::to_string(states) +
          ", \"states_long\": " + std::to_string(states_long) + ", \"groups\": " +
          std::to_string(groups) + ", \"w\": " + json_num(w_all) + ", \"w_long\": " +
          json_num(w_long) + ", \"top5_total\": " + std::to_string(top_total[0]) +
          ", \"top20_total\": " + std::to_string(top_total[1]);
    const char* roles[2] = {"stem1", "both"};
    for (int r = 0; r < 2; ++r)
      js += ", " + json_str(roles[r]) + ": {\"pass_n\": " + json_arr(pass_n[r]) +
            ", \"pass_w\": " + json_arr(pass_w[r]) + ", \"pass_w_long\": " + json_arr(pass_wl[r]) +
            ", \"top5_pass\": " + json_arr(top_pass[r][0]) + ", \"top20_pass\": " +
            json_arr(top_pass[r][1]) + "}";
    js += "}";
  }
  return js + "}}";
}

// After outside: "Z id n logZ", then "S id table a b c d w" for every stored state of every
// beamed table with w = exp(log inside + log outside - log Z) >= min_w (keys as in
// state_key), then per table "T id table states N w W below_states K below_w V" (all
// states with w > 0; below = w < min_w, not written).  With dangles the Qm, Qm2, Qmp, PKM
// weight sums the cell's run states P, U1, U2 (needs keep_runs).
void DP03::dump_states(FILE* f, const std::string& id, double min_w) {
  outside();
  const char* sid = id.c_str();
  std::fprintf(f, "Z %s %d %.17g\n", sid, n, logZ);
  long long cnt[BT_N] = {}, below[BT_N] = {};
  double tot[BT_N] = {}, belw[BT_N] = {};
  auto emit = [&](int bt, int a, int b, int c, int d, double lw) {
    if (lw < HALF_NEG || logZ < HALF_NEG) return;
    const double w = std::exp(lw - logZ);
    ++cnt[bt];
    tot[bt] += w;
    if (w < min_w) {
      ++below[bt];
      belw[bt] += w;
      return;
    }
    std::fprintf(f, "S %s %s %d %d %d %d %.9e\n", sid, BT_NAMES[bt], a, b, c, d, w);
  };
  auto two = [&](int bt, const std::vector<OMap>& at, const SpMat& beta, int rs) {
    for (int j = 0; j < n; ++j)
      for (size_t r = 0; r < at[j].slots(); ++r) {
        int i = at[j].keys[r];
        if (i < 0) continue;
        double lw = (dz && rs >= 0) ? runW[rs][ix(i, j)] : at[j].vals[r] + beta[ix(i, j)];
        emit(bt, i, j, -1, -1, lw);
      }
  };
  two(BT_QB, Qb_at, Qb_b, -1);
  two(BT_QM, Qm_at, Qm_b, RS_QM);
  two(BT_QM2, Qm2_at, Qm2_b, RS_QM2);
  if (allow_pk && E.stem_ml) two(BT_QMP, Qmp_at, Qmp_b, RS_QMP);
  if (allow_pk) {
    two(BT_QP, Qp_at, Qp_b, -1);
    for (int j = 0; j < n; ++j)
      for (size_t r = 0; r < PKM_at[j].slots(); ++r) {
        int d = PKM_at[j].keys[r];
        if (d < 0) continue;
        double lw = dz ? runW[RS_PKM][ix(d, j)] : PKM_at[j].vals[r] + PKM_b[ix(j, d)];
        emit(BT_PKM, d, j, -1, -1, lw);
      }
    const std::pair<int, std::vector<Tab4>*> tabs[5] = {
        {BT_QG, &Qg}, {BT_QGLS, &Qgls}, {BT_QGL, &Qgl}, {BT_QGR, &Qgr}, {BT_QGRS, &Qgrs}};
    for (const auto& tb : tabs)
      for (int q = 0; q < n; ++q)
        for (const Grp& g : (*tb.second)[q].grps)
          for (const Ent& e : g.ents) emit(tb.first, e.p, e.u, g.v, q, e.a + e.b);
  }
  for (int t = 0; t < BT_N; ++t)
    std::fprintf(f, "T %s %s states %lld w %.9e below_states %lld below_w %.9e\n", sid,
                 BT_NAMES[t], cnt[t], tot[t], below[t], belw[t]);
  std::fflush(f);
}

// After outside of a step-mode run with every table at infinite beam: "Z id n logZ", then per
// table and step j with N > 0 stored states (the candidates of that pruning point) "E table
// j N K" and "C table j a b c d w r0 r1 cls" for each candidate with w >= min_w (w as in
// dump_states, keys as in state_key).  r0 = rank (0 = best) by the default pruning score
// (Qb, Qm, Qm2, Qmp, Qp, PKM: log C[i-1] + alpha, prune2's order; 4-index: four_score,
// keep_mask's order); r1 = rank by the --stem-est rank score of st with role all (Qg, Qgl,
// Qgr, partner table A built from the feasible Qgl states as in prune4_stem; infeasible
// states last; K = feasible count) and r1 = r0, K = N for the other tables.  cls (Qg, Qgl):
// stem-1 share of the outside from the tagged outside, 1 (>= 0.999) / 2 (<= 0.001) / 3
// (mixed) / 0 (w = 0); -1 for the other tables.  Qg / Qgl / Qgr lines then add "r2 r3 r4"
// (ranks as r1 of the combined score with (--stem-est-qg-r, --stem-est-consts) = (mirror,
// consistent), (partner, as-is), (partner, consistent)) and 10 tie-averaged ranks "a0 .. a9"
// of: the combined score (mirror as-is, mirror consistent, partner as-is, partner consistent),
// sigma_L (as-is, consistent), sigma_R (mirror as-is, mirror consistent, partner as-is,
// partner consistent); infeasible states and missing roles score -1e300.  Their E lines add
// the sum and the sum of squares of each a over all N candidates.
void DP03::dump_scores(FILE* f, const std::string& id, const StemEst& st, double min_w) {
  outside();
  std::fprintf(f, "Z %s %d %.17g\n", id.c_str(), n, logZ);
  std::vector<Tab4>* t4[5] = {&Qg, &Qgls, &Qgl, &Qgr, &Qgrs};
  const int b4[5] = {BT_QG, BT_QGLS, BT_QGL, BT_QGR, BT_QGRS};
  auto flat_b = [&](int t, std::vector<double>& out) {
    out.clear();
    for (Tab4& tb : *t4[t])
      for (Grp& g : tb.grps)
        for (Ent& e : g.ents) out.push_back(e.b);
  };
  std::vector<double> b0[5], b1;
  std::vector<signed char> cls[5];
  for (int t = 0; t < 5; ++t) flat_b(t, b0[t]);
  tagged_outside(1);
  for (int t : {0, 2}) {
    flat_b(t, b1);
    cls[t].resize(b1.size());
    for (size_t k = 0; k < b1.size(); ++k) {
      if (b0[t][k] < HALF_NEG) {
        cls[t][k] = 0;
        continue;
      }
      double share = b1[k] > HALF_NEG ? std::exp(b1[k] - b0[t][k]) : 0.0;
      cls[t][k] = share >= 0.999 ? 1 : share <= 0.001 ? 2 : 3;
    }
  }
  for (int t = 0; t < 5; ++t) {
    size_t k = 0;
    for (Tab4& tb : *t4[t])
      for (Grp& g : tb.grps)
        for (Ent& e : g.ents) e.b = b0[t][k++];
    b0[t].clear();
    b0[t].shrink_to_fit();
  }

  auto weight = [&](double lw) { return lw < HALF_NEG ? 0.0 : std::exp(lw - logZ); };
  auto two = [&](int bt, const std::vector<OMap>& at, const SpMat* beta, int rs,
                 bool pkm) {
    std::vector<std::pair<double, int>> ranked;
    std::vector<size_t> rank(n);
    for (int j = 0; j < n; ++j) {
      ranked.clear();
      for (size_t r = 0; r < at[j].slots(); ++r) {
        int i = at[j].keys[r];
        if (i >= 0) ranked.emplace_back(prefix_log(i) + at[j].vals[r], i);
      }
      if (ranked.empty()) continue;
      std::sort(ranked.begin(), ranked.end(),
                [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                  if (a.first != b.first) return a.first > b.first;
                  return a.second > b.second;
                });
      for (size_t k = 0; k < ranked.size(); ++k) rank[ranked[k].second] = k;
      std::fprintf(f, "E %s %d %zu %zu\n", BT_NAMES[bt], j, ranked.size(), ranked.size());
      for (size_t r = 0; r < at[j].slots(); ++r) {
        int i = at[j].keys[r];
        if (i < 0) continue;
        double lw = (dz && rs >= 0) ? runW[rs][ix(i, j)]
                                    : at[j].vals[r] + (*beta)[pkm ? ix(j, i) : ix(i, j)];
        double w = weight(lw);
        if (w < min_w || w == 0.0) continue;
        std::fprintf(f, "C %s %d %d %d -1 -1 %.9e %zu %zu -1\n", BT_NAMES[bt], j, i, j, w,
                     rank[i], rank[i]);
      }
    }
  };
  two(BT_QB, Qb_at, &Qb_b, -1, false);
  two(BT_QM, Qm_at, &Qm_b, RS_QM, false);
  two(BT_QM2, Qm2_at, &Qm2_b, RS_QM2, false);
  if (E.stem_ml) two(BT_QMP, Qmp_at, &Qmp_b, RS_QMP, false);
  two(BT_QP, Qp_at, &Qp_b, -1, false);
  two(BT_PKM, PKM_at, &PKM_b, RS_PKM, true);

  stA.clear();
  stE.clear();
  stc_q = -1;
  set_stem(&st, true, true);
  // Stem tables also get the ranks of the combined score under (qg_r, consts) = (mirror,
  // consistent), (partner, as-is), (partner, consistent) and the raw sigma_L (as-is,
  // consistent) and sigma_R (mirror as-is / consistent, partner as-is / consistent).
  static const bool VAR[3][2] = {{false, true}, {true, false}, {true, true}};
  // Tie-averaged ranks (0 = best) of the descending values.
  auto avg_ranks = [](const std::vector<double>& x) {
    std::vector<size_t> o(x.size());
    for (size_t k = 0; k < o.size(); ++k) o[k] = k;
    std::stable_sort(o.begin(), o.end(), [&](size_t a, size_t b) { return x[a] > x[b]; });
    std::vector<double> r(x.size());
    for (size_t a = 0; a < o.size();) {
      size_t b = a;
      while (b + 1 < o.size() && x[o[b + 1]] == x[o[a]]) ++b;
      for (size_t k = a; k <= b; ++k) r[o[k]] = 0.5 * static_cast<double>(a + b);
      a = b + 1;
    }
    return r;
  };
  size_t off[5] = {};
  std::vector<double> s0, s1, sv[3], raw[6];
  for (int q = 0; q < n; ++q)
    for (int t = 0; t < 5; ++t) {
      Tab4& tb = (*t4[t])[q];
      const int bt = b4[t];
      const bool stem = bt == BT_QG || bt == BT_QGL || bt == BT_QGR;
      std::unique_ptr<StemEval> ev;
      if (stem) ev.reset(new StemEval(*this, tb, bt, q));
      s0.clear();
      s1.clear();
      for (int x = 0; x < 3; ++x) sv[x].clear();
      for (int x = 0; x < 6; ++x) raw[x].clear();
      size_t k = 0, K = 0;
      for (const Grp& g : tb.grps)
        for (const Ent& e : g.ents) {
          s0.push_back(four_score(e.p, e.u, g.v, e.a));
          if (stem) {
            const bool okk = ev->ok[k];
            s1.push_back(okk ? ev->score(k) : NEG_INF);
            for (int x = 0; x < 3; ++x) {
              stem_qg_partner = VAR[x][0];
              stem_consistent = VAR[x][1];
              sv[x].push_back(okk ? ev->score(k) : NEG_INF);
            }
            for (int c = 0; c < 2; ++c) {
              stem_consistent = c == 1;
              stem_qg_partner = false;
              raw[c].push_back(okk ? ev->sigL(k) : NEG_INF);
              raw[2 + c].push_back(okk ? ev->sigR(k) : NEG_INF);
              stem_qg_partner = true;
              raw[4 + c].push_back(okk ? ev->sigR(k) : NEG_INF);
            }
            stem_qg_partner = stem_consistent = false;
            K += okk;
          }
          ++k;
        }
      const size_t N = k;
      if (bt == BT_QGL) stem_add_A(q, &ev->ok);
      if (N == 0) continue;
      std::vector<size_t> r0 = score_ranks(s0);
      std::vector<size_t> r1 = stem ? score_ranks(s1) : r0;
      std::vector<size_t> rv[3];
      std::vector<double> av[10];
      if (stem) {
        for (int x = 0; x < 3; ++x) rv[x] = score_ranks(sv[x]);
        av[0] = avg_ranks(s1);
        for (int x = 0; x < 3; ++x) av[1 + x] = avg_ranks(sv[x]);
        for (int x = 0; x < 6; ++x) av[4 + x] = avg_ranks(raw[x]);
      }
      std::fprintf(f, "E %s %d %zu %zu", BT_NAMES[bt], q, N, stem ? K : N);
      if (stem)
        for (int x = 0; x < 10; ++x) {
          double s = 0.0, s2 = 0.0;
          for (double r : av[x]) {
            s += r;
            s2 += r * r;
          }
          std::fprintf(f, " %.17g %.17g", s, s2);
        }
      std::fputc('\n', f);
      k = 0;
      for (const Grp& g : tb.grps)
        for (const Ent& e : g.ents) {
          double w = weight(e.a + e.b);
          if (w >= min_w && w > 0.0) {
            std::fprintf(f, "C %s %d %d %d %d %d %.9e %zu %zu %d", BT_NAMES[bt], q, e.p, e.u,
                         g.v, q, w, r0[k], r1[k], t == 0 || t == 2 ? cls[t][off[t] + k] : -1);
            if (stem) {
              std::fprintf(f, " %zu %zu %zu", rv[0][k], rv[1][k], rv[2][k]);
              for (int x = 0; x < 10; ++x) std::fprintf(f, " %.1f", av[x][k]);
            }
            std::fputc('\n', f);
          }
          ++k;
        }
      off[t] += N;
    }
  stem_ = nullptr;
  std::fflush(f);
}

// --dump-ext: after the inside pass of a step-mode run with every table at infinite beam (no
// outside), the Qg candidates of every step under the Result-18 scores, with the partner
// table A built from the feasible Qgl states as in prune4_stem (role all, qg-r partner,
// as-is).  "Z id n", then per step with N > 0 Qg candidates "E Qg j N K xs xb" (K feasible,
// xs / xb extension nodes visited, stack / bulge1) followed by the sum and the sum of squares
// over all N of 7 tie-averaged ranks (0 = best) of: partner (logsumexp(sigma_L, sigma_R)),
// partner + stack, partner + bulge1 (logsumexp with sigma_ext), sigma_L, sigma_R, sigma_ext
// stack, sigma_ext bulge1; infeasible states and missing terms score -1e300.  Then
// "C Qg j a b c d w rp rs rb ks kb a0 .. a6" for each candidate in w_of (the --dump-exact
// states, w as there): rp / rs / rb = keep_mask-order ranks of the three combined scores, ks /
// kb = bit m set if split selection (feasible states only) keeps it at b = 300, 1000, 3000
// (m = 0, 1, 2), stack / bulge1; a0 .. a6 as on the E line.
void DP03::dump_ext(FILE* f, const std::string& id, const StemEst& st,
                    const std::unordered_map<uint64_t, double>& w_of, int kmax) {
  stA.clear();
  stE.clear();
  stc_q = -1;
  set_stem(&st, true, true, true, false);
  stem_ext_max = kmax;
  std::fprintf(f, "Z %s %d\n", id.c_str(), n);
  auto avg_ranks = [](const std::vector<double>& x) {
    std::vector<size_t> o(x.size());
    for (size_t k = 0; k < o.size(); ++k) o[k] = k;
    std::stable_sort(o.begin(), o.end(), [&](size_t a, size_t b) { return x[a] > x[b]; });
    std::vector<double> r(x.size());
    for (size_t a = 0; a < o.size();) {
      size_t b = a;
      while (b + 1 < o.size() && x[o[b + 1]] == x[o[a]]) ++b;
      for (size_t k = a; k <= b; ++k) r[o[k]] = 0.5 * static_cast<double>(a + b);
      a = b + 1;
    }
    return r;
  };
  static const int BEAMS[3] = {300, 1000, 3000};
  std::vector<double> v[7];  // cp, cs, cb, sL, sR, xs, xb
  for (int q = 0; q < n; ++q) {
    {
      StemEval ev(*this, Qg[q], BT_QG, q);
      const size_t N = ev.ent.size();
      if (N > 0) {
        for (auto& x : v) x.assign(N, NEG_INF);
        size_t K = 0;
        for (size_t k = 0; k < N; ++k) {
          if (!ev.ok[k]) continue;
          ++K;
          v[3][k] = ev.sigL(k);
          v[4][k] = ev.sigR(k);
          v[0][k] = ev.score(k);
        }
        long long nodes[2];
        for (int m = 0; m < 2; ++m) {
          stem_ext = m + 1;
          stX.clear();
          const long long n0 = ext_nodes;
          for (size_t k = 0; k < N; ++k)
            if (ev.ok[k]) {
              const Ent& e = *ev.ent[k];
              v[5 + m][k] = stem_ext_score(e.p, e.u, ev.grp[k]->v, e.a);
            }
          nodes[m] = ext_nodes - n0;
        }
        stX.clear();
        stem_ext = 0;
        for (size_t k = 0; k < N; ++k)
          if (ev.ok[k]) {
            v[1][k] = log_add(v[0][k], v[5][k]);
            v[2][k] = log_add(v[0][k], v[6][k]);
          }
        std::vector<size_t> rk[3];
        for (int x = 0; x < 3; ++x) rk[x] = score_ranks(v[x]);
        std::vector<double> av[7];
        for (int x = 0; x < 7; ++x) av[x] = avg_ranks(v[x]);
        // split selection over the feasible states, as prune4_stem
        std::vector<unsigned char> kept[2];
        std::vector<size_t> idx;
        for (size_t k = 0; k < N; ++k)
          if (ev.ok[k]) idx.push_back(k);
        std::vector<double> sl, sr, sx, sc;
        std::vector<char> km;
        for (int m = 0; m < 2; ++m) {
          kept[m].assign(N, 0);
          sl.clear();
          sr.clear();
          sx.clear();
          sc.clear();
          for (size_t k : idx) {
            sl.push_back(v[3][k]);
            sr.push_back(v[4][k]);
            sx.push_back(v[5 + m][k]);
            sc.push_back(v[1 + m][k]);
          }
          const std::vector<double>* parts[3] = {&sl, &sr, &sx};
          for (int bi = 0; bi < 3; ++bi) {
            if (idx.size() <= static_cast<size_t>(BEAMS[bi])) {
              for (size_t k : idx) kept[m][k] |= 1 << bi;
              continue;
            }
            split_mask(parts, sc, km, BEAMS[bi]);
            for (size_t r = 0; r < idx.size(); ++r)
              if (km[r]) kept[m][idx[r]] |= 1 << bi;
          }
        }
        std::fprintf(f, "E Qg %d %zu %zu %lld %lld", q, N, K, nodes[0], nodes[1]);
        for (int x = 0; x < 7; ++x) {
          double s = 0.0, s2 = 0.0;
          for (double r : av[x]) {
            s += r;
            s2 += r * r;
          }
          std::fprintf(f, " %.17g %.17g", s, s2);
        }
        std::fputc('\n', f);
        for (size_t k = 0; k < N; ++k) {
          const Ent& e = *ev.ent[k];
          const int vv = ev.grp[k]->v;
          auto it = w_of.find(state_key(BT_QG, e.p, e.u, vv, q));
          if (it == w_of.end()) continue;
          std::fprintf(f, "C Qg %d %d %d %d %d %.9e %zu %zu %zu %d %d", q, e.p, e.u, vv, q,
                       it->second, rk[0][k], rk[1][k], rk[2][k], kept[0][k], kept[1][k]);
          for (int x = 0; x < 7; ++x) std::fprintf(f, " %.1f", av[x][k]);
          std::fputc('\n', f);
        }
      }
    }
    StemEval evl(*this, Qgl[q], BT_QGL, q);
    stem_add_A(q, &evl.ok);
  }
  stem_ = nullptr;
  std::fflush(f);
}

