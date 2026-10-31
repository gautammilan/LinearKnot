// Scan helpers: pseudoknot suffix Sbar_pk, items, per-pass profile and counters, beam caps,
// table views, b-collection; scan_ct_dump.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  // One pass's time split (s): init (l2r_init_scores: arm entries + prYw + prSF), growth and
  // selection (counted where they run), Rs = spawn (build + score), all tables' score and build
  // (Rs included), close; other = pass CPU - init - build - score - close (nested tables).
  std::string l2r_prof(double pass_cpu) const {
    double bs = 0.0, ps = 0.0;
    for (int t = 0; t < PT_N; ++t) {
      bs += l2r_bsec[t];
      ps += l2r_psec[t];
    }
    char b[512];
    std::snprintf(b, sizeof b,
                  "cpu %.2f init %.2f static %.2f grow %.2f sel %.2f spawn %.2f (cand %lld) "
                  "score %.2f build %.2f close %.2f other %.2f",
                  pass_cpu, l2r_isec, pr_static_sec, pr_grow_sec, pr_sel_sec,
                  l2r_bsec[PT_RS] + l2r_psec[PT_RS], sp_cand, ps, bs, l2r_csec,
                  pass_cpu - l2r_isec - bs - ps - l2r_csec);
    return b;
  }
  static int pct_p99(const long long* h, int nbin) {
    long long tot = 0;
    for (int i = 0; i < nbin; ++i) tot += h[i];
    if (tot <= 0) return -1;
    const long long need = (tot * 99 + 99) / 100;
    long long c = 0;
    for (int i = 0; i < nbin; ++i) {
      c += h[i];
      if (c >= need) return i;
    }
    return nbin - 1;
  }
  void scan_ct_dump() const {
    auto grp = [&](int a, int b, long long& pre, long long& kept, long long& mx, long long& bind,
                   long long& pos, int& beam) {
      pre = kept = mx = bind = pos = 0;
      beam = ltb(a);
      for (int t = a; t <= b; ++t) {
        pre += l2r_pre[t];
        kept += l2r_kept[t];
        mx = std::max(mx, l2r_kmax[t]);
        bind += l2r_bind[t];
        pos += l2r_pos[t];
      }
    };
    auto one = [&](const char* name, int a, int b) {
      long long pre, kept, mx, bind, pos;
      int beam;
      grp(a, b, pre, kept, mx, bind, pos, beam);
      const double per = pos > 0 ? static_cast<double>(pre) / static_cast<double>(pos) : 0.0;
      std::fprintf(stderr,
                   " %s beam %d pre %lld kept %lld max %lld bind %lld pos %lld cand_per_pos %.3f",
                   name, beam, pre, kept, mx, bind, pos, per);
    };
    long long m_n = 0, m_sum = 0, list_n = 0, list_sum = 0;
    for (int i = 0; i < 8; ++i) {
      m_n += sx_m_h[i];
      m_sum += static_cast<long long>(i) * sx_m_h[i];
    }
    for (int i = 0; i < 65; ++i) {
      list_n += sx_list_h[i];
      list_sum += static_cast<long long>(i) * sx_list_h[i];
    }
    const double k_per = sx_rkey_n > 0 ? static_cast<double>(sx_k_src[2]) / static_cast<double>(sx_rkey_n)
                                       : 0.0;
    long long rpre = 0;
    for (int t = PT_R; t <= PT_RRS; ++t) rpre += l2r_pre[t];
    const double rkey_per = rpre > 0 ? static_cast<double>(sx_rkey_n) / static_cast<double>(rpre) : 0.0;
    std::fprintf(stderr, "[scan-ct] n %d", n);
    one("P1", PT_P1, PT_P1);
    one("P2", PT_P2, PT_P2);
    one("P3", PT_P3, PT_P3);
    one("R", PT_R, PT_RRS);
    one("S", PT_S, PT_S);
    one("U", PT_U, PT_URS);
    std::fprintf(stderr,
                 " rkey_calls %lld r_cand %lld rkey_per_cand %.6f k_rkey %lld k_per_rkey %.3f "
                 "k_root %lld k_spawn %lld k_n %lld m_n %lld m_mean %.3f m_p99 %d "
                 "qlist_n %lld qlist_mean %.3f qlist_p99 %d\n",
                 sx_rkey_n, rpre, rkey_per, sx_k_src[2], k_per, sx_k_src[1], sx_k_src[3], sx_k_n,
                 m_n, m_n > 0 ? static_cast<double>(m_sum) / static_cast<double>(m_n) : 0.0,
                 pct_p99(sx_m_h, 8), list_n,
                 list_n > 0 ? static_cast<double>(list_sum) / static_cast<double>(list_n) : 0.0,
                 pct_p99(sx_list_h, 65));
  }


