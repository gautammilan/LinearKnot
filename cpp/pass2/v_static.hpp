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
    std::stable_sort(pairs.begin(), pairs.end(),
                     [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
                       return a.second - a.first > b.second - b.first;
                     });
    std::vector<AT> ats;
    std::vector<std::vector<AT>> bk;
    if (pr_vstream) bk.resize(n);
    {
      pr_static_ats(true, pairs, ats, pr_vstream ? &bk : nullptr);
    }
    mtrim();
    long long bk_n = 0;
    for (const std::vector<AT>& v : bk) bk_n += static_cast<long long>(v.size());
    std::sort(ats.begin(), ats.end(), at_less);
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
    if (pr_vstream || g_mtrim) {
      std::vector<AT>().swap(ats);
      mtrim();
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

