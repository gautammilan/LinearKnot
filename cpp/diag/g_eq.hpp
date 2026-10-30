// Diagnostic: g_eq_write (--g-eq PATH) writes sweep prGamma and sw_fill against a fresh lazy
// growth vector on the final scan prefixes.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

  // Sweep prGamma (and the unused sw_fill) against a fresh lazy growth vector.
  // Both sides use the final scan prefix. Writes --g-eq.
  void g_eq_write() {
    std::ofstream out(g_eq_path);
    if (!out) {
      std::fprintf(stderr, "[g-eq] cannot write %s\n", g_eq_path.c_str());
      return;
    }
    auto rel = [](double a, double b) {
      if (a < HALF_NEG && b < HALF_NEG) return 0.0;
      if (a < HALF_NEG || b < HALF_NEG) return 1e300;
      return std::fabs(std::expm1(a - b));
    };
    auto agg = [&](const std::unordered_map<uint64_t, double>* src, std::vector<double>& v) {
      v.assign(static_cast<size_t>(n), NEG_INF);
      if (!src) return;
      for (const auto& [k, w] : *src) {
        if (w < HALF_NEG) continue;
        const int q = k_lo(k);
        if (q < 0 || q >= n) continue;
        v[static_cast<size_t>(q)] = log_add(v[static_cast<size_t>(q)], w + prefix_log(k_hi(k)));
      }
    };
    auto from_rn = [&](const RN& r, std::vector<double>& v) {
      v.assign(static_cast<size_t>(n), NEG_INF);
      if (g_ntr_K > 0) {
        for (const auto& [q, g] : r.T[0])
          if (q >= 0 && q < n && g > HALF_NEG) v[static_cast<size_t>(q)] = log_add(v[static_cast<size_t>(q)], g);
        return;
      }
      for (int m = 0; m < rnL; ++m) {
        const RV& a = r.L[m];
        for (size_t i = 0; i < a.v.size(); ++i) {
          if (!(a.v[i] > 0.0)) continue;
          const int q = a.lo + static_cast<int>(i);
          if (q < 0 || q >= n) continue;
          v[static_cast<size_t>(q)] = log_add(v[static_cast<size_t>(q)], a.mx + std::log(a.v[i]));
        }
      }
    };
    try {
      sw_fill(0, n);
    } catch (const std::exception& e) {
      out << "fill_error " << e.what() << "\n";
      std::fprintf(stderr, "[g-eq] sw_fill %s\n", e.what());
    }
    const bool saved = p2_sweep;
    p2_sweep = false;
    prRT.clear();
    prRR.clear();
    prRG.clear();
    prRH.clear();
    std::vector<std::pair<int, int>> pairs;
    {
      std::unordered_set<uint64_t> seen;
      const std::vector<std::vector<int>>& src = pr_pair_list(true);
      for (int x = 0; x < n; ++x)
        for (int y : src[x]) {
          const uint64_t k = key2(x, y);
          if (seen.insert(k).second) pairs.emplace_back(x, y);
        }
      if (stem_)
        for (const auto& kv : stem_->prGamma) {
          const uint64_t k = kv.first;
          if (seen.insert(k).second) pairs.emplace_back(k_hi(k), k_lo(k));
        }
    }
    std::fprintf(stderr, "[g-eq] pairs %zu\n", pairs.size());
    struct Bad {
      int x, y, q;
      double hy, gr, ad, r;
    };
    std::vector<Bad> top;
    auto consider = [&](int x, int y, int q, double hy, double gr) {
      const double r = rel(hy, gr);
      if (!(r > 1e-9)) return;
      const double ad = (hy < HALF_NEG || gr < HALF_NEG) ? 1e300 : std::fabs(hy - gr);
      top.push_back({x, y, q, hy, gr, ad, r});
    };
    long long ncmp = 0, nbad = 0, ncmp_f = 0, nbad_f = 0;
    std::unordered_set<uint64_t> bad_pair;
    std::vector<double> hy, gr, fl, spec_h, spec_g;
    bool have_spec = false;
    int built = 0;
    for (const auto& [x, y] : pairs) {
      const RN& rn = rec_root(x, y, false);
      from_rn(rn, gr);
      const std::unordered_map<uint64_t, double>* hs = nullptr;
      {
        auto it = sw_ga[0].find(key2(x, y));
        if (it != sw_ga[0].end()) hs = &it->second;
      }
      agg(hs, hy);
      const std::unordered_map<uint64_t, double>* fs = nullptr;
      auto fit = sw_ga[0].find(key2(x, y));
      if (fit != sw_ga[0].end()) fs = &fit->second;
      agg(fs, fl);
      bool any = false;
      for (int q = 0; q < n; ++q) {
        const double a = hy[static_cast<size_t>(q)], b = gr[static_cast<size_t>(q)];
        if (a < HALF_NEG && b < HALF_NEG) continue;
        ++ncmp;
        if (rel(a, b) > 1e-9) {
          ++nbad;
          any = true;
          consider(x, y, q, a, b);
        }
        const double c = fl[static_cast<size_t>(q)];
        if (c < HALF_NEG && b < HALF_NEG) continue;
        ++ncmp_f;
        if (rel(c, b) > 1e-9) ++nbad_f;
      }
      if (any) bad_pair.insert(key2(x, y));
      if (x == 4 && y == 45) {
        spec_h = hy;
        spec_g = gr;
        have_spec = true;
      }
      if (++built % 400 == 0) std::fprintf(stderr, "[g-eq] built %d\n", built);
    }
    p2_sweep = saved;
    long long ncmp_p = 0, nbad_p = 0, nbad_pf = 0;
    {
      geq_pre = true;
      p2_sweep = false;
      prRT.clear();
      prRR.clear();
      prRG.clear();
      prRH.clear();
      int built_p = 0;
      for (const auto& [x, y] : pairs) {
        const RN& rn = rec_root(x, y, false);
        from_rn(rn, gr);
        const std::unordered_map<uint64_t, double>* hs = nullptr;
        {
          auto it = sw_ga[0].find(key2(x, y));
          if (it != sw_ga[0].end()) hs = &it->second;
        }
        agg(hs, hy);
        const std::unordered_map<uint64_t, double>* fs = nullptr;
        auto fit = sw_ga[0].find(key2(x, y));
        if (fit != sw_ga[0].end()) fs = &fit->second;
        agg(fs, fl);
        for (int q = 0; q < n; ++q) {
          const double a = hy[static_cast<size_t>(q)], b = gr[static_cast<size_t>(q)],
                       c = fl[static_cast<size_t>(q)];
          if (a < HALF_NEG && b < HALF_NEG) continue;
          ++ncmp_p;
          if (rel(a, b) > 1e-9) ++nbad_p;
          if (!(c < HALF_NEG && b < HALF_NEG) && rel(c, b) > 1e-9) ++nbad_pf;
        }
        if (x == 1 && y == 5) {
          out << "pre15\n";
          for (int q = 0; q < n; ++q) {
            const double a = hy[static_cast<size_t>(q)], b = gr[static_cast<size_t>(q)],
                         c = fl[static_cast<size_t>(q)];
            if (a < HALF_NEG && b < HALF_NEG && c < HALF_NEG) continue;
            out << "pq " << q << " " << std::setprecision(17) << a << " " << c << " " << b << "\n";
          }
        }
        if (++built_p % 400 == 0) std::fprintf(stderr, "[g-eq] pre %d\n", built_p);
      }
      geq_pre = false;
      prRT.clear();
      prRR.clear();
      prRG.clear();
      prRH.clear();
      p2_sweep = saved;
    }
    out << "pre " << ncmp_p << " " << nbad_p << " " << nbad_pf << "\n";
    std::sort(top.begin(), top.end(), [](const Bad& a, const Bad& b) {
      if (a.ad != b.ad) return a.ad > b.ad;
      if (a.y != b.y) return a.y < b.y;
      if (a.x != b.x) return a.x < b.x;
      return a.q < b.q;
    });
    if (top.size() > 10) top.resize(10);
    out << "cmp " << ncmp << " " << nbad << "\n";
    out << "fill " << ncmp_f << " " << nbad_f << "\n";
    out << "pairs " << pairs.size() << " bad_pairs " << bad_pair.size() << "\n";
    for (size_t i = 0; i < top.size(); ++i)
      out << "top " << i << " " << top[i].x << " " << top[i].y << " " << top[i].q << " " << std::setprecision(17)
          << top[i].hy << " " << top[i].gr << " " << top[i].ad << " " << top[i].r << "\n";
    out << "vec " << (have_spec ? 1 : 0) << "\n";
    if (have_spec) {
      for (int q = 0; q < n; ++q) {
        const double a = spec_h[static_cast<size_t>(q)], b = spec_g[static_cast<size_t>(q)];
        if (a < HALF_NEG && b < HALF_NEG) continue;
        out << "vq " << q << " " << std::setprecision(17) << a << " " << b << "\n";
      }
    }
    std::unordered_map<uint64_t, int> depth;
    std::function<int(int, int)> dep = [&](int x, int y) -> int {
      const uint64_t k = key2(x, y);
      auto it = depth.find(k);
      if (it != depth.end()) return it->second;
      depth.emplace(k, 0);
      std::vector<ChainNbr> nb;
      pr_steps(x, y, false, true, nb);
      int d = 0;
      if (!nb.empty()) {
        d = n + 1;
        for (const ChainNbr& t : nb) d = std::min(d, 1 + dep(t.x, t.y));
      }
      depth[k] = d;
      return d;
    };
    int bx = -1, by = -1, bd = n + 2;
    for (uint64_t k : bad_pair) {
      const int x = k_hi(k), y = k_lo(k);
      const int d = dep(x, y);
      if (d < bd || (d == bd && (y < by || (y == by && x < bx)))) {
        bd = d;
        bx = x;
        by = y;
      }
    }
    out << "short " << bx << " " << by << " " << (bx < 0 ? -1 : bd) << "\n";
    p2_sweep = false;
    long long ibad = 0;
    int ix = -1, iy = -1, idp = n + 2;
    for (uint64_t k : bad_pair) {
      const int x = k_hi(k), y = k_lo(k);
      std::vector<ChainNbr> nb;
      pr_steps(x, y, false, true, nb);
      bool ml = false;
      for (const ChainNbr& t : nb)
        if (t.t != 0) ml = true;
      if (ml) continue;
      ++ibad;
      const int d = dep(x, y);
      if (d < idp || (d == idp && (y < iy || (y == iy && x < ix)))) {
        idp = d;
        ix = x;
        iy = y;
      }
    }
    out << "ibad " << ibad << " ishort " << ix << " " << iy << " " << (ix < 0 ? -1 : idp) << "\n";
    if (bx >= 0) {
      int x = bx, y = by;
      for (int hop = 0; hop < n && x >= 0; ++hop) {
        const RN& rn = rec_root(x, y, false);
        from_rn(rn, gr);
        const std::unordered_map<uint64_t, double>* hs = nullptr;
        {
          auto it = sw_ga[0].find(key2(x, y));
          if (it != sw_ga[0].end()) hs = &it->second;
        }
        agg(hs, hy);
        int nd = 0;
        for (int q = 0; q < n; ++q)
          if (rel(hy[static_cast<size_t>(q)], gr[static_cast<size_t>(q)]) > 1e-9) ++nd;
        std::vector<ChainNbr> nb;
        pr_steps(x, y, false, true, nb);
        out << "node " << hop << " " << x << " " << y << " depth " << dep(x, y) << " ndiff " << nd << " nstep "
            << nb.size() << " at " << std::setprecision(17) << (E.end_ok(x, y) ? E.at(x, y) : NEG_INF) << " pref "
            << prefix_log(x) << "\n";
        if (hop == 0) {
          for (int q = 0; q < n; ++q) {
            const double a = hy[static_cast<size_t>(q)], b = gr[static_cast<size_t>(q)];
            if (a < HALF_NEG && b < HALF_NEG) continue;
            out << "sq " << q << " " << std::setprecision(17) << a << " " << b << "\n";
            if (!(rel(a, b) > 1e-9)) continue;
            out << "diffq " << q << " " << std::setprecision(17) << a << " " << b << "\n";
            if (hs) {
              for (const auto& [k, w] : *hs)
                if (k_lo(k) == q)
                  out << "raw " << k_hi(k) << " " << q << " " << std::setprecision(17) << w << "\n";
            }
            for (size_t s = 0; s < nb.size(); ++s) {
              const RN& ch = rec_root(nb[s].x, nb[s].y, false);
              std::vector<double> cg;
              from_rn(ch, cg);
              const double cv = cg[static_cast<size_t>(q)];
              if (cv < HALF_NEG) continue;
              out << "via " << s << " " << q << " " << std::setprecision(17) << (cv + nb[s].w) << "\n";
            }
          }
        }
        int nxt = -1, ny = -1, bestd = n + 3;
        for (size_t s = 0; s < nb.size(); ++s) {
          const int d = dep(nb[s].x, nb[s].y);
          out << "step " << s << " t " << nb[s].t << " " << nb[s].x << " " << nb[s].y << " " << std::setprecision(17)
              << nb[s].w << " depth " << d << " bad " << (bad_pair.count(key2(nb[s].x, nb[s].y)) ? 1 : 0) << "\n";
          if (d < bestd) {
            bestd = d;
            nxt = nb[s].x;
            ny = nb[s].y;
          }
        }
        if (nxt < 0 || dep(x, y) == 0) break;
        x = nxt;
        y = ny;
      }
    }
    out.flush();
    std::fprintf(stderr, "[g-eq] cmp %lld bad %lld fill %lld bad %lld\n", ncmp, nbad, ncmp_f, nbad_f);
  }

