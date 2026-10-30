// Q^z (gap partition function) and the exterior C: fill_Qz, fill_Qz_kept (top b2 starts by
// C(k-1) Q^z(k, j)), fill_C; outside allocation.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  void fill_Qz(int j) {
    const int pv = (j + 1) & 1, cu = j & 1;  // (j - 1) & 1 without a negative j - 1
    const double d5 = dz && j > 0 ? E.d5w[j - 1] : 0.0;
    const double s1 = dz ? E.s1w[j] : 0.0, d3 = dz ? E.d3w[j] : 0.0;
    for (int i = 0; i <= j; ++i) {
      double acc, u1 = NEG_INF, u2 = NEG_INF;
      if (dz) {
        if (i == j) {
          u1 = E.wu;  // empty region before j
        } else {
          double p0 = zP[pv][i];
          u1 = p0 > HALF_NEG ? p0 + E.wu : NEG_INF;
          double t = log_add(zU2[pv][i], zU1[pv][i] + d5);
          u2 = t > HALF_NEG ? t + E.wu : NEG_INF;
        }
        acc = NEG_INF;
      } else {
        acc = (i == j) ? E.wu : Qz[ix(i, j - 1)] + E.wu;
      }
      const OMap& atb = Qb_at[j];
      for (size_t r = 0; r < atb.slots(); ++r) {
        int d = atb.keys[r];
        if (d < 0 || d < i || !E.end_ok(d, j)) continue;
        double left = (d == i) ? 0.0 : Qz[ix(i, d - 1)];
        if (left < HALF_NEG) continue;
        acc = log_add(acc, left + atb.vals[r] + (E.w_nest + E.at(d, j)));
      }
      if (allow_pk) {
        const OMap& atp = Qp_at[j];
        for (size_t r = 0; r < atp.slots(); ++r) {
          int d = atp.keys[r];
          if (d < 0 || d < i) continue;
          double left = (d == i) ? 0.0 : Qz[ix(i, d - 1)];
          if (left < HALF_NEG) continue;
          acc = log_add(acc, left + atp.vals[r] + E.w_beta1p + E.nest_pk_k * E.w_nest);
        }
      }
      if (dz) {
        zP[cu][i] = acc;
        zU1[cu][i] = u1;
        zU2[cu][i] = u2;
        acc = run_close(acc, u1, u2, s1, d3);
      }
      Qz[ix(i, j)] = acc;
    }
  }

  // Every length. Top b2 starts at this right end, by C(k-1) + Q^z(k, j)
  // (C_hat on the scan, C_pre on the pre-pass). Dropped starts stay absent.
  void fill_Qz_kept(int j) {
    if (qz_keep.empty()) qz_keep.assign(static_cast<size_t>(n), {});
    const int B = p2_qz_beam;
    const double d5 = dz && j > 0 ? E.d5w[j - 1] : 0.0;
    const double s1 = dz ? E.s1w[j] : 0.0, d3 = dz ? E.d3w[j] : 0.0;
    struct Acc {
      double acc, u1, u2;
    };
    std::vector<int> seen(static_cast<size_t>(j + 1), 0), touch;
    std::vector<Acc> cell(static_cast<size_t>(j + 1));
    auto cell_at = [&](int k) -> Acc& {
      if (!seen[static_cast<size_t>(k)]) {
        seen[static_cast<size_t>(k)] = 1;
        cell[static_cast<size_t>(k)] = {NEG_INF, NEG_INF, NEG_INF};
        touch.push_back(k);
      }
      return cell[static_cast<size_t>(k)];
    };
    if (dz) cell_at(j).u1 = E.wu;
    else cell_at(j).acc = E.wu;
    if (j > 0) {
      for (const QkEnt& e : qz_keep[static_cast<size_t>(j - 1)]) {
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
    auto add_split = [&](int d, double w) {
      if (d < 0 || w < HALF_NEG) return;
      Acc& self = cell_at(d);
      self.acc = log_add(self.acc, w);
      if (d <= 0) return;
      for (const QkEnt& e : qz_keep[static_cast<size_t>(d - 1)]) {
        if (e.q < HALF_NEG || e.k > j) continue;
        Acc& a = cell_at(e.k);
        a.acc = log_add(a.acc, e.q + w);
      }
    };
    const OMap& atb = Qb_at[j];
    for (size_t r = 0; r < atb.slots(); ++r) {
      const int d = atb.keys[r];
      if (d < 0 || !E.end_ok(d, j)) continue;
      add_split(d, atb.vals[r] + (E.w_nest + E.at(d, j)));
    }
    if (allow_pk) {
      const OMap& atp = Qp_at[j];
      for (size_t r = 0; r < atp.slots(); ++r) {
        const int d = atp.keys[r];
        if (d < 0) continue;
        add_split(d, atp.vals[r] + E.w_beta1p + E.nest_pk_k * E.w_nest);
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
      const double key = q + (k > 0 ? C[static_cast<size_t>(k - 1)] : 0.0);
      cs.push_back({k, q, dz ? a.acc : 0.0, a.u1, a.u2, key});
    }
    qz_cand += static_cast<long long>(cs.size());
    if (g_beam_stats.f) {
      std::vector<double> keys;
      keys.reserve(cs.size());
      for (const Cand& c : cs) keys.push_back(c.key);
      std::sort(keys.begin(), keys.end(), [](double a, double b) { return a > b; });
      const size_t kept_n =
          B > 0 && static_cast<int>(cs.size()) > B ? static_cast<size_t>(B) : keys.size();
      const int full = B > 0 && static_cast<int>(cs.size()) > B ? 1 : 0;
      g_beam_stats.log_prune("Qz", j, keys.size(), kept_n, BeamStats::mass_share(keys, kept_n),
                             full);
    }
    if (B > 0 && static_cast<int>(cs.size()) > B) {
      std::nth_element(cs.begin(), cs.begin() + B, cs.end(),
                       [](const Cand& a, const Cand& b) { return a.key > b.key; });
      cs.resize(static_cast<size_t>(B));
    }
    qz_kept += static_cast<long long>(cs.size());
    std::sort(cs.begin(), cs.end(), [](const Cand& a, const Cand& b) { return a.k < b.k; });
    std::vector<QkEnt>& col = qz_keep[static_cast<size_t>(j)];
    col.clear();
    col.reserve(cs.size());
    for (const Cand& c : cs) {
      Qz[ix(c.k, j)] = c.q;
      col.push_back({c.k, c.q, c.p, c.u1, c.u2});
    }
  }

  void fill_C(int j) {
    double acc, u1 = NEG_INF, u2 = NEG_INF;
    if (dz) {
      u1 = cP;
      u2 = j > 0 ? log_add(cU2, cU1 + E.d5w[j - 1]) : NEG_INF;
      acc = NEG_INF;
    } else {
      acc = (j == 0) ? 0.0 : C[j - 1];
    }
    const OMap& atb = Qb_at[j];
    for (size_t r = 0; r < atb.slots(); ++r) {
      int i = atb.keys[r];
      if (i < 0 || !E.end_ok(i, j)) continue;
      ++r_extc;
      double left = (i == 0) ? 0.0 : C[i - 1];
      if (left > HALF_NEG) acc = log_add(acc, left + atb.vals[r] + E.at(i, j));
    }
    if (allow_pk) {
      const OMap& atp = Qp_at[j];
      for (size_t r = 0; r < atp.slots(); ++r) {
        int i = atp.keys[r];
        if (i < 0) continue;
        double left = (i == 0) ? 0.0 : C[i - 1];
        if (left > HALF_NEG) acc = log_add(acc, left + atp.vals[r] + E.w_beta1);
      }
    }
    if (dz) {
      cP = acc;
      cU1 = u1;
      cU2 = u2;
      acc = run_close(acc, u1, u2, E.s1w[j], E.d3w[j]);
    }
    C[j] = acc;
    if (c_trace) {
      double qb = NEG_INF, qp = NEG_INF;
      const OMap& atb = Qb_at[j];
      for (size_t r = 0; r < atb.slots(); ++r) {
        int i = atb.keys[r];
        if (i < 0 || !E.end_ok(i, j)) continue;
        double left = (i == 0) ? 0.0 : C[i - 1];
        if (left > HALF_NEG) qb = log_add(qb, left + atb.vals[r] + E.at(i, j));
      }
      if (allow_pk) {
        const OMap& atp = Qp_at[j];
        for (size_t r = 0; r < atp.slots(); ++r) {
          int i = atp.keys[r];
          if (i < 0) continue;
          double left = (i == 0) ? 0.0 : C[i - 1];
          if (left > HALF_NEG) qp = log_add(qp, left + atp.vals[r] + E.w_beta1);
        }
      }
      std::fprintf(stderr, "[c] %d %.17g %.17g %.17g\n", j, C[j], qb, qp);
    }
  }

  // ============================================================= outside
  template <class T>
  void add2b(T& t, int i, int j, double val) {
    if (val < HALF_NEG || i < 0 || j < 0) return;
    const size_t k = ix(i, j);
    t[k] = log_add(static_cast<double>(t[k]), val);
  }

  void alloc_beta() {
    Qb_b.init(n);
    Qm_b.init(n);
    Qm2_b.init(n);
    if (E.stem_ml) Qmp_b.init(n);
    Qz_b.init(n);
    Qp_b.init(n);
    if (allow_pk) PKM_b.init(n);
    C_b.assign(n, NEG_INF);
    C_b[n - 1] = 0.0;
    if (dz) {
      for (int k = 0; k < 2; ++k) {
        rbQm[k].init(n);
        rbQm2[k].init(n);
        rbQmp[k].init(n);
        rbPKM[k].init(n);
        bz1[k].assign(n, NEG_INF);
        bz2[k].assign(n, NEG_INF);
      }
      bc1 = bc2 = NEG_INF;
    }
  }

