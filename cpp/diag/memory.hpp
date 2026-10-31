// Memory accounting: log_nest ([nest] lines), state counts, *_sz byte estimates and mem_items
// ([mem-nt] line).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  // states per table: R, Rs, Rls, Rlss, Rrs, U, Us, Uls, Ulss, Urs, S
  void log_nest(const char* tag) const {
    auto oml = [](const std::vector<OMap>& v) {
      long long e = 0;
      for (const OMap& m : v) e += m.live;
      return e;
    };
    const long long nn = static_cast<long long>(n) * n;
    std::fprintf(stderr,
                 "[nest] %s n %d dense %lld Qb %lld Qm %lld Qm2 %lld Qp %lld Qmp %lld Qz %lld "
                 "qz_mb %.3f\n",
                 tag, n, nn, oml(Qb_at), oml(Qm_at), oml(Qm2_at), oml(Qp_at), oml(Qmp_at),
                 Qz.finite(), static_cast<double>(Qz.bytes()) / 1048576.0);
  }

  std::vector<long long> l2r_states() const {
    std::vector<long long> c(11, 0);
    auto cnt = [](const std::vector<Tab4>& v) {
      long long s = 0;
      for (const Tab4& t : v)
        for (const Grp& g : t.grps) s += static_cast<long long>(g.ents.size());
      return s;
    };
    const LFam* fs[2] = {&LA, &LB};
    for (int k = 0; k < 2; ++k) {
      c[5 * k + 0] = cnt(fs[k]->G);
      c[5 * k + 1] = cnt(fs[k]->Gs);
      c[5 * k + 2] = cnt(fs[k]->Gls);
      c[5 * k + 3] = cnt(fs[k]->Glss);
      c[5 * k + 4] = cnt(fs[k]->Grs);
    }
    for (const SMap& m : LS) c[10] += static_cast<long long>(m.keys.size());
    return c;
  }

  // Per structure entries and bytes (MemItem).  tabs: the engine tables (inside
  // dense + kept maps, outside, b16 4-index); pk: V / arms, K, growth, l2r tables, caches.
  static void tab4_sz(const std::vector<Tab4>& v, long long& e, double& b) {
    b += msz::vec(v);
    for (const Tab4& t : v) {
      b += msz::vec(t.grps) + msz::umap(t.vidx);
      for (const Grp& g : t.grps) {
        e += static_cast<long long>(g.ents.size());
        b += msz::vec(g.ents) + msz::umap(g.idx);
      }
    }
  }
  static void omap_sz(const std::vector<OMap>& v, long long& e, double& b) {
    b += msz::vec(v);
    for (const OMap& m : v) {
      e += m.live;
      b += msz::vec(m.keys) + msz::vec(m.vals) + m.pos.bytes();
    }
  }
  template <class M>
  static double pmap_sz(const M& m, long long& e) {
    double b = msz::umap(m);
    for (const auto& kv : m) {
      e += static_cast<long long>(kv.second.k.size());
      b += msz::vec(kv.second.k) + msz::vec(kv.second.w) + msz::vec(kv.second.lin);
    }
    return b;
  }
  template <class M>
  static double rmap_sz(const M& m, long long& e) {
    double b = msz::umap(m);
    for (const auto& kv : m) {
      e += relems(kv.second);
      for (int s = 0; s < 3; ++s) b += msz::vec(kv.second.L[s].v) + msz::vec(kv.second.R[s].v);
      b += msz::vec(kv.second.T[0]) + msz::vec(kv.second.T[1]);
    }
    return b;
  }
  template <class M>
  static double profmap_sz(const M& m, long long& e) {
    double b = msz::umap(m);
    for (const auto& kv : m) {
      const Prof& p = kv.second;
      e += static_cast<long long>(p.ly.size() + p.rx.size());
      b += msz::vec(p.ly) + msz::vec(p.rx) + msz::vec(p.lw) + msz::vec(p.rw);
    }
    return b;
  }
  template <class M>
  static double vmap_sz(const M& m, long long& e) {
    double b = msz::umap(m);
    for (const auto& kv : m) {
      e += static_cast<long long>(kv.second.size());
      b += msz::vec(kv.second);
    }
    return b;
  }
  void mem_items(std::vector<MemItem>& out, const std::string& pre, bool tabs, bool pk) const {
    if (tabs) {
      long long e = Qz.finite();
      double b = Qz.bytes() + msz::vec(C);
      for (int k = 0; k < 2; ++k) b += msz::vec(zP[k]) + msz::vec(zU1[k]) + msz::vec(zU2[k]);
      for (const std::vector<OMap>* v : {&Qb_at, &Qm_at, &Qm2_at, &Qp_at, &Qmp_at, &PKM_at}) omap_sz(*v, e, b);
      out.push_back({pre + "in", e, b});
      e = 0;
      b = Qz_b.bytes() + msz::vec(C_b);
      e += Qz_b.finite();
      for (const SpMat* v : {&Qb_b, &Qm_b, &Qm2_b, &Qp_b, &PKM_b, &Qmp_b, &QmC_b, &Qm2C_b}) {
        e += v->finite();
        b += v->bytes();
      }
      for (int k = 0; k < 2; ++k) b += msz::vec(bz1[k]) + msz::vec(bz2[k]);
      out.push_back({pre + "out", e, b});
      e = 0;
      b = 0.0;
      for (const std::vector<Tab4>* v : {&Qg, &Qgls, &Qgrs, &Qgl, &Qgr}) tab4_sz(*v, e, b);
      out.push_back({pre + "tab4", e, b});
    }
    if (!pk) return;
    {
      double b = msz::umap(lArm) + msz::vec(lYwS) + msz::vec(prYw) + msz::vec(prEw);
      for (const std::vector<int>* v : {&prArmXs, &prArmXp, &prArmY0, &prYy, &prArmE0, &prEy, &prEys,
                                        &prXsFirst, &prByXp, &prXpFirst})
        b += msz::vec(*v);
      out.push_back({"V_arms", static_cast<long long>(l2r_arm_ents()), b});
    }
    {
      double b = msz::vec(prSF) + msz::umap(prAQ);
      long long e = static_cast<long long>(prSF.size());
      for (const auto& kv : prAQ) {
        const ArmQ& a = kv.second;
        e += static_cast<long long>(a.phi.size());
        b += msz::vec(a.phi) + msz::vec(a.rs) + msz::vec(a.t1) + msz::vec(a.t2) + msz::vec(a.ra) +
             msz::vec(a.P1) + msz::vec2(a.D) + vmap_sz(a.W, e);
      }
      out.push_back({"K_sfx", e, b});
    }
    long long e = 0;
    double b = rmap_sz(prRG, e);
    out.push_back({"grow_nodes", e, b});
    e = 0;
    b = rmap_sz(prRR, e);
    out.push_back({"grow_roots", e, b});
    e = 0;
    b = rmap_sz(prRH, e);
    out.push_back({"grow_help", e, b});
    e = 0;
    b = rmap_sz(prRKls, e) + rmap_sz(prRKrs, e);
    out.push_back({"grow_keys", e, b});
    e = 0;
    b = rside_bytes(e);
    out.push_back({"grow_R_side", e, b});
    e = 0;
    b = rfv_bytes(e);
    out.push_back({"V_keep", e, b});
    e = 0;
    b = msz::umap(prRT);
    for (const auto& kv : prRT) {
      e += static_cast<long long>(kv.second.L.size() + kv.second.R.size());
      b += msz::vec(kv.second.L) + msz::vec(kv.second.R);
    }
    out.push_back({"grow_trunc", e, b});
    e = 0;
    b = 0.0;
    for (const LFam* F : {&LA, &LB})
      for (const std::vector<Tab4>* v : {&F->G, &F->Gs, &F->Gls, &F->Glss, &F->Grs}) tab4_sz(*v, e, b);
    b += msz::vec(LS);
    for (const SMap& m : LS) {
      e += static_cast<long long>(m.keys.size());
      b += msz::vec(m.keys) + msz::vec(m.v) + msz::umap(m.idx);
    }
    out.push_back({"l2r_tab", e, b});
    if (l2r_inc) {
      e = 0;
      b = 0.0;
      for (const std::vector<SMap>* T : {&LP1, &LP2, &LP3}) {
        b += msz::vec(*T);
        for (const SMap& m : *T) {
          e += static_cast<long long>(m.keys.size());
          b += msz::vec(m.keys) + msz::vec(m.v) + msz::umap(m.idx);
        }
      }
      out.push_back({"l2r_P", e, b});
      e = 0;
      b = msz::vec(lH) + msz::vec(lHdone) + msz::vec2(lArmOpen) + msz::vec(lBopen);
      for (const LH& h : lH) {
        b += msz::vec(h.c1) + msz::vec(h.aopen) + msz::vec2(h.v);
        for (const auto& v : h.v) e += static_cast<long long>(v.size());
      }
      out.push_back({"l2r_H", e, b});
    }
    e = 0;
    b = pmap_sz(prG, e) + pmap_sz(prGR, e) + pmap_sz(prHg, e) + pmap_sz(prGRs, e) + pmap_sz(prKls, e) +
        pmap_sz(prKrs, e);
    out.push_back({"b16_grow_prof", e, b});
    e = 0;
    b = msz::umap(prTc);
    for (const auto& kv : prTc) {
      b += msz::umap(kv.second);
      for (const auto& t : kv.second) {
        e += static_cast<long long>(t.second.ys.size());
        b += msz::vec(t.second.ys) + msz::vec(t.second.w);
      }
    }
    out.push_back({"b16_T", e, b});
    e = 0;
    b = msz::umap(prKc);
    for (const auto& kv : prKc) {
      e += static_cast<long long>(kv.second.is.size());
      b += msz::vec(kv.second.is) + msz::vec(kv.second.w);
    }
    out.push_back({"b16_kappa_K", e, b});
    e = 0;
    b = msz::umap(prUL) + msz::vec(prUR) + msz::umap(prUc) + msz::vec(prEBr) + msz::vec(prEBc) +
        msz::vec(prMBr) + msz::vec(prMBc) + msz::vec(prKap);
    for (const auto& kv : prUL) {
      e += static_cast<long long>(kv.second.w.size());
      b += msz::vec(kv.second.w) + msz::vec(kv.second.lin);
    }
    for (const PrUR& r : prUR) {
      e += static_cast<long long>(r.w.size());
      b += msz::vec(r.w) + msz::vec(r.lin);
    }
    out.push_back({"grow_sel_U", e, b});
    e = 0;
    b = msz::umap(stA) + msz::umap(stE) + msz::umap(stBc) + msz::umap(stBgc) + vmap_sz(stX, e) +
        msz::umap(stHc) + msz::umap(stOc) + profmap_sz(stP, e) + profmap_sz(stPg, e) + profmap_sz(stPls, e) +
        profmap_sz(stPrs, e);
    out.push_back({"b16_stem", e, b});
    e = 0;
    b = profmap_sz(chM, e) + profmap_sz(chR, e) + msz::umap(prSc) + msz::umap(prAc);
    out.push_back({"b16_chain_naive", e, b});
  }

