// class DP03, private data: the tables (Qb, Qm, Qm2, Qmp, Qz, C, Qg, Qgl, Qgr, Qp, PKM, ...),
// beams, run states (U1, U2, P) and their helpers, alloc, pair lists, the pass-3 Q^z filler
// qz_log / four_score, arm scores, and prune4 (top-b cut of a 4-index column).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

 private:
  Energy E;
  int n;
  int beam;
  bool allow_pk;
  int cap;
  bool score_qz;
  bool gap_step;
  // --stem-est: partner table A(d, e, f) of the surviving stem-1 Qgl states, per key2(d, e)
  // with f ascending; first_f[e] = first f of any entry (., e, f)
  struct StemA {
    std::vector<int> fs;
    std::vector<double> vals;
    std::vector<std::vector<std::pair<int, double>>> il;  // pr_on: (p, alpha + at(p, f)) per f
  };
  const StemEst* stem_ = nullptr;
  bool stem_both = false, stem_qg2 = false, stem_qg_partner = false, stem_consistent = false;
  std::unordered_map<uint64_t, StemA> stA;
  stemidx::MinTree first_f;
  // e values with an entry (d, e, .) per d, ascending; B(d, e, v) and Bg(d, u, v) caches of
  // the step stc_q (all their terms end before the step, so they are final)
  std::unordered_map<int, std::vector<int>> stE;
  std::unordered_map<uint64_t, double> stBc, stBgc;
  int stc_q = -1;
  // sigma_ext: growth nodes per outer pair key2(p, q) and H(u, v, q') per key3, of step
  // stc_q; helper outside O(x, y; u, v) per key4 (mode 3 = loop)
  struct ExtNode {
    int d, q;
    double w;
  };
  int stem_ext = 0, stem_ext_max = 30;
  bool stem_split = false;
  int grow_lg = 3, grow_nodes = 2000, stem_ml_span = 30;
  double grow_cut = 1e-8;
  int stem_helpers = 0;
  // nested helper sums for step hn_q: Qgls by 2 * i + (v == q), Qgrs by c; count < 0 = unset
  int hn_q = -1;
  std::vector<std::pair<double, int>> hnL, hnR;
  std::unordered_map<uint64_t, std::vector<ExtNode>> stX;
  std::unordered_map<uint64_t, double> stHc, stOc;
  // Stem profiles: sparse log profiles of a key over the final right ends y' (L part,
  // coefficient of Qz[u+1][v-1] H_R(u, v, y')) and left ends x' (R part, coefficient of
  // Bg(x', u, v)), index ascending; tail: the Qgls ml tail (coefficient of Qz[u+1][v-1]);
  // ncl: closings of a helper key.  Per step: untruncated pair profiles by key2(x, y), and
  // the scored key profiles of Qg (by p), Qgls (by 2 i + seed) and Qgrs (by c).
  struct Prof {
    std::vector<int> ly, rx;
    std::vector<double> lw, rw;
    double tail = NEG_INF;
    int ncl = 0;
  };
  bool stem_profile = false, ml_tail = false;
  double prof_eps = 0.0;
  int pf_q = -1;
  std::unordered_map<uint64_t, Prof> stP;
  std::unordered_map<int, Prof> stPg, stPls, stPrs;
  // prof_eval's dense H / Bg rows of the current gap (u, v) (stamped).
  std::vector<double> rowH, rowB;
  std::vector<int> rowHs, rowBs;
  int row_u = -1, row_v = -1, row_q = -1, row_stamp = 0;
  int tbeam[BT_N];
  bool inside_done = false;
  bool outside_done = false;
  double logZ = NEG_INF;

  AtView Qb, Qm, Qm2, Qp;
  SpMat Qz;
  std::vector<double> C;
  std::vector<OMap> Qb_at, Qm_at, Qm2_at, Qp_at;
  // Qmp: Qm with the stem-multiloop penalties b_p, c_p; qmp_at() / qmp_b() fall back to
  // Qm unless E.stem_ml
  AtView Qmp;
  SpMat Qmp_b;
  std::vector<OMap> Qmp_at;
  const OMap& qmp_at(int j) const { return E.stem_ml ? Qmp_at[j] : Qm_at[j]; }
  SpMat& qmp_b() { return E.stem_ml ? Qmp_b : Qm_b; }
  // Per column, Qm slot indices bucketed by the base at key-1 (the left end of a case-(h) pair).
  mutable std::vector<std::array<std::vector<int>, 4>> qm_bkt;
  mutable std::vector<char> qm_bkt_ok;
  void qm_reindex(int j) {
    if (static_cast<int>(qm_bkt.size()) != n) {
      qm_bkt.assign(n, {});
      qm_bkt_ok.assign(n, 0);
    }
    auto& b = qm_bkt[j];
    for (auto& v : b) v.clear();
    const OMap& at = qmp_at(j);
    for (size_t r = 0; r < at.slots(); ++r) {
      const int i = at.keys[r];
      if (i < 1) continue;
      const int base = E.code[i - 1];
      if (base < 0 || base > 3) continue;
      b[base].push_back(static_cast<int>(r));
    }
    qm_bkt_ok[j] = 1;
  }
  // Slots of Qm at `col` whose base at key-1 can pair with the base at `right`.
  // Falls back to every slot when the column has not been bucketed.
  template <class F>
  void qm_pair_slots(int col, int right, F&& fn) const {
    const OMap& at = qmp_at(col);
    if (static_cast<int>(qm_bkt_ok.size()) == n && qm_bkt_ok[col]) {
      const int bj = E.code[right];
      for (int bi = 0; bi < 4; ++bi) {
        const int t = pair_type(bi, bj);
        if (t < 0 || (!allow_gu && (t == 3 || t == 5))) continue;
        for (int r : qm_bkt[col][bi]) fn(r);
      }
      return;
    }
    for (size_t r = 0; r < at.slots(); ++r) fn(static_cast<int>(r));
  }
  std::vector<OMap> PKM_at;  // PKM_at[j][d] = sum_{e<=j} Qp[d][e] e^{-α3 (j-e)}, beamed
  std::vector<Tab4> Qg, Qgls, Qgrs, Qgl, Qgr;
  std::vector<std::vector<int>> partners_left, partners_right;
  std::vector<double> qb_acc;  // fill_Qb scratch, one row, reset on the partners just written
  std::vector<std::vector<LeadEnt>> qlead;
  // Pass 1 (l2r_mode == 0) enumerates like LinearPartition, via mate_next only.
  // Interior loops use the energy cap. The multiloop leading run uses ml_lead_cap.
  // Multi key is (p, q, kind): 0 empty right run, 1 one base, 2 length >= 2.
  // A length >= 2 extension does not depend on where the run started, so those
  // arrivals merge. Empty and one-base runs do: the dangle is s1w of that base.
  std::vector<int> mate_after[4];
  bool mates_ready = false;
  std::vector<OMap> lpH;   // hairpin, column = right end, key = left end
  std::vector<OMap> lpM;   // Multi inside, key = p * 3 + kind
  std::vector<OMap> lpMb;  // Multi outside
  std::vector<OMap> lpP;   // M2 weight that ends on a branch, no trailing run

  SpMat Qb_b, Qm_b, Qm2_b, Qp_b, PKM_b;
  SpMat Qz_b;
  std::vector<double> C_b;

  // Dangles: trailing-run states (P, U1, U2) of the tables that grow by unpaired
  // bases, indexed by column parity; the stored table value closes the run for a
  // paired j+1: X = P + U1 s1w[j] + U2 d3w[j].  Outside keeps the U1, U2 outsides.
  bool dz = false;
  RunSt rsQm[2], rsQm2[2], rsQmp[2], rsPKM[2];
  std::vector<double> zP[2], zU1[2], zU2[2];
  double cP = 0.0, cU1 = NEG_INF, cU2 = NEG_INF;
  RunSt rbQm[2], rbQm2[2], rbQmp[2], rbPKM[2];
  std::vector<double> bz1[2], bz2[2];
  double bc1 = NEG_INF, bc2 = NEG_INF;

  static std::vector<int> live_keys(const OMap& m) {
    std::vector<int> out;
    for (size_t r = 0; r < m.slots(); ++r)
      if (m.keys[r] >= 0) out.push_back(m.keys[r]);
    return out;
  }

  // U1, U2 at j from the run states pv at j-1 of the survivors surv.
  void run_ext(const RunSt& pv, RunSt& cu, const std::vector<int>& surv, int j, double w) {
    const double d5 = E.d5w[j - 1];
    for (int i : surv) {
      double p = pv.P.get(i);
      if (p > HALF_NEG) cu.U1.put(i, p + w);
      double t = log_add(pv.U2.get(i), pv.U1.get(i) + d5);
      if (t > HALF_NEG) cu.U2.put(i, t + w);
    }
  }

  static double run_close(double p, double u1, double u2, double s1, double d3) {
    double v = p;
    if (u1 > HALF_NEG) v = log_add(v, u1 + s1);
    if (u2 > HALF_NEG) v = log_add(v, u2 + d3);
    return v;
  }

  // Closed values at j into at (and tbl unless null); keys first, then new keys of P.
  void run_fin(AtView* tbl, OMap& at, const RunSt& c, const std::vector<int>& first,
               int j) {
    const double s1 = E.s1w[j], d3 = E.d3w[j];
    auto one = [&](int i) {
      if (at.has(i)) return;
      double v = run_close(c.P.get(i), c.U1.get(i), c.U2.get(i), s1, d3);
      if (v > HALF_NEG) {
        if (tbl) (*tbl)[ix(i, j)] = v;
        at.set(i, v);
      }
    };
    for (int i : first) one(i);
    for (int i : c.P.keys) one(i);
  }

  void run_back(double bF, double bU1n, double bU2n, int j, double w, double& bP, double& bU1,
                double& bU2) const {
    bP = log_add(bF, bU1n + w);
    bU1 = log_add(bF + E.s1w[j], bU2n + (E.d5w[j] + w));
    bU2 = log_add(bF + E.d3w[j], bU2n + w);
  }

  // Replaces the outside of the closed values at j by that of P; keeps U1's, U2's.
  template <class Beta>
  void run_back_at(Beta& beta, const OMap& at, int j, double w, RunSt* rb) {
    const RunSt& nx = rb[(j + 1) & 1];
    RunSt& cu = rb[j & 1];
    cu.clear();
    for (size_t r = 0; r < at.slots(); ++r) {
      int i = at.keys[r];
      if (i < 0) continue;
      double bP, b1, b2;
      run_back(static_cast<double>(beta[ix(i, j)]), nx.U1.get(i), nx.U2.get(i), j, w, bP, b1, b2);
      beta[ix(i, j)] = bP;
      cu.U1.put(i, b1);
      cu.U2.put(i, b2);
    }
  }



  size_t ix(int i, int j) const { return static_cast<size_t>(i) * n + j; }

  void alloc() {
    Qz.init(n);
    C.assign(n, NEG_INF);
    auto init_maps = [&](std::vector<OMap>& v) {
      v.resize(n);
      for (auto& m : v) m.init(n);
    };
    init_maps(Qb_at);
    init_maps(Qm_at);
    init_maps(Qm2_at);
    init_maps(Qp_at);
    init_maps(PKM_at);
    Qb.bind(&Qb_at, n);
    Qm.bind(&Qm_at, n);
    Qm2.bind(&Qm2_at, n);
    Qp.bind(&Qp_at, n);
    if (E.stem_ml) {
      init_maps(Qmp_at);
      Qmp.bind(&Qmp_at, n);
    }
    dz = E.dangle;
    if (dz) {
      for (int k = 0; k < 2; ++k) {
        rsQm[k].init(n);
        rsQm2[k].init(n);
        rsQmp[k].init(n);
        rsPKM[k].init(n);
        zP[k].assign(n, NEG_INF);
        zU1[k].assign(n, NEG_INF);
        zU2[k].assign(n, NEG_INF);
      }
    }
    Qg.resize(n);
    Qgls.resize(n);
    Qgrs.resize(n);
    Qgl.resize(n);
    Qgr.resize(n);
  }

  void build_pairs() {
    partners_left.assign(n, {});
    partners_right.assign(n, {});
    for (int j = 0; j < n; ++j)
      for (int i = 0; i <= j - MIN_PAIR; ++i)
        if (E.pairs(i, j)) {
          partners_left[j].push_back(i);
          partners_right[i].push_back(j);
        }
  }
  // Pass 1 does not build this list. The scan and the stem index do, on first use.
  void ensure_pairs() const {
    if (n <= 0 || !partners_left.empty()) return;
    const_cast<DP03*>(this)->build_pairs();
  }

  // ---------------------------------------------------------------- scores
  double qz_log(int left, int right) const {
    if (left > right) return 0.0;
    if (left < 0 || right >= n) return NEG_INF;
    ++qz_reads;
    return Qz[ix(left, right)];
  }

  double prefix_log(int p) const {
    if (p <= 0) return 0.0;
    return C[p - 1];
  }

  double four_score(int p, int u, int v, double log_x) const {
    double s = prefix_log(p) + log_x;
    if (score_qz) s += qz_log(u + 1, v - 1);
    return s;
  }

  // gtab: GATE_QG / GATE_QGL for the tables the gate may act on (stats are kept for them
  // either way), -1 for the others.
  void prune4(Tab4& t, int bt, int q, int gtab = -1) {
    const int beam = tbeam[bt];
    if (stem_ && (bt == BT_QG || bt == BT_QGL || bt == BT_QGR)) {
      prune4_stem(t, bt, q, gtab, beam);
      return;
    }
    if (stem_ && (stem_helpers || pr_on) && (bt == BT_QGLS || bt == BT_QGRS)) {
      prune4_helper(t, bt, q, beam);
      return;
    }
    if (beam <= 0) {
      if (gtab >= 0)
        for (const Grp& g : t.grps) {
          gate_stats.cand[gtab] += g.ents.size();
          gate_stats.elig[gtab] += g.ents.size();
          gate_stats.kept[gtab] += g.ents.size();
        }
      return;
    }
    if (gap_step) {
      prune4_step(t, gtab, beam);
      return;
    }
    std::vector<double> scores;
    std::vector<char> keep;
    for (auto& g : t.grps) {
      if (gtab >= 0) {
        gate_stats.cand[gtab] += g.ents.size();
        gate_stats.elig[gtab] += g.ents.size();
      }
      if (static_cast<int>(g.ents.size()) <= beam) {
        if (gtab >= 0) gate_stats.kept[gtab] += g.ents.size();
        continue;
      }
      scores.clear();
      for (const Ent& e : g.ents)
        scores.push_back(four_score(e.p, e.u, g.v, e.a));
      keep_mask(scores, keep, beam);
      filter_group(g, keep, 0);
      if (gtab >= 0) {
        gate_stats.ranked[gtab] += scores.size();
        gate_stats.kept[gtab] += g.ents.size();
      }
    }
  }

