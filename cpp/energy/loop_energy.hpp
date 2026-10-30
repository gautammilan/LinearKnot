// struct Energy: Boltzmann weights of loops for one sequence, e^{-dG/RT}: hairpin, interior /
// bulge / stack, multiloop (a, b, c) and exterior terms, dangles (some | all | d2), terminal
// AU/GU, and the pseudoknot weights used by the DP03 recursions.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

struct Energy {
  std::string seq;
  int n = 0;
  std::vector<int> code;
  std::vector<char> cp;  // n*n, can_pair
  double w_alpha1, w_alpha2, w_alpha3, w_beta1, w_beta1m, w_beta1p, w_stem, w_nest,
      w_beta3, wu, nest_pk_k, ml_pk_k, st_mult, int_mult, w_ap, w_bp, w_cp;
  bool stem_mult, stem_ml;

  explicit Energy(const std::string& s) : seq(s), n(static_cast<int>(s.size())) {
    code.resize(n);
    for (int i = 0; i < n; ++i) code[i] = base_code(seq[i]);
    cp.assign(static_cast<size_t>(n) * n, 0);
    for (int i = 0; i < n; ++i)
      for (int j = i + MIN_PAIR; j < n; ++j) {
        int t = pair_type(code[i], code[j]);
        if (t < 0 || (!allow_gu && (t == 3 || t == 5))) continue;
        cp[static_cast<size_t>(i) * n + j] = 1;
      }
    if (allowed_pairs) {
      std::vector<char> keep(cp.size(), 0);
      for (const auto& p : *allowed_pairs) {
        size_t k = static_cast<size_t>(p.first) * n + p.second;
        if (0 <= p.first && p.first < p.second && p.second < n) keep[k] = cp[k];
      }
      cp.swap(keep);
    }
    w_alpha1 = log_weight(PK.alpha1);
    w_alpha2 = log_weight(PK.alpha2);
    w_alpha3 = log_weight(PK.alpha3);
    w_beta1 = log_weight(PK.beta1);
    w_beta1m = log_weight(PK.beta1m);
    w_beta1p = log_weight(PK.beta1p);
    w_stem = log_weight(PK.stem);
    w_nest = log_weight(PK.nest);
    w_beta3 = log_weight(PK.beta3);
    wu = w_beta3;
    nest_pk_k = PK.nest_pk_k;
    ml_pk_k = PK.ml_pk_k;
    st_mult = PK.st_mult;
    int_mult = PK.int_mult;
    stem_mult = st_mult != 1.0 || int_mult != 1.0;
    w_ap = log_weight(PK.a_p);
    w_bp = log_weight(PK.b_p);
    w_cp = log_weight(PK.c_p);
    stem_ml = PK.stem_ml;
    non_c.resize(n);
    crun.assign(n + 1, 0);
    for (int i = 0; i < n; ++i) {
      non_c[i] = code[i] != 1;
      crun[i + 1] = crun[i] + (code[i] == 1);
    }
    if (mat_rna1999) w_at_pen = log_weight(T99.at_penalty);
    if (mat_v04) {
      w_alpha1 = log_weight(V04.ml_closing / 100.0);
      w_alpha2 = log_weight(V04.ml_intern / 100.0);
      w_alpha3 = log_weight(V04.ml_base / 100.0);
      max_left = V04_MAX_LEFT;
      ml_left = V04_MAX_LEFT;
    }
    dangle = mat_rna1999 && dangles_mode != 0;
    if (dangle) {
      d5w.resize(n);
      d3w.resize(n);
      s1w.resize(n);
      for (int k = 0; k < n; ++k) {
        // Non-ACGU (N): dangle / mismatch terms that read it are 0 (LinearPartition-style).
        double d5 = 0.0, d3 = 0.0;
        if (k > 0 && code[k - 1] >= 0 && code[k] >= 0)
          d5 = T99.dangle[code[k - 1] * 4 + code[k]];
        if (k < n - 1 && code[k + 1] >= 0 && code[k] >= 0)
          d3 = T99.dangle[24 + (3 - code[k + 1]) * 4 + code[k]];
        d5w[k] = log_weight(d5);
        d3w[k] = log_weight(d3);
        s1w[k] = log_weight(dangles_mode == 1 && 0 < k && k < n - 1 ? std::min(d3, d5)
                                                                    : d3 + d5);
      }
    }
  }

  std::vector<char> non_c;
  std::vector<int> crun;  // crun[k] = number of C in seq[0, k)
  double w_at_pen = 0.0;

  int max_left = NO_LIMIT;  // interior loop: most unpaired on the left
  int ml_left = NO_LIMIT;   // multiloop: most unpaired before the first branch

  // Log weight of the terminal AU/GU penalty of pair i·j (0.0 for dp03); for vienna2004 also
  // its -d2 mismatch / dangle.
  double at(int i, int j) const {
    if (mat_v04) return log_weight(branch04(i, j) / 100.0);
    return non_c[i] && non_c[j] ? w_at_pen : 0.0;
  }

  // at() for pair i·j closing a multiloop.
  double at_close(int i, int j) const {
    if (mat_v04) return log_weight(closing04(i, j) / 100.0);
    return at(i, j);
  }

  // Dangles (ene.c ExplDangleRaw): an unpaired run s..t of an exterior loop, multiloop
  // or pseudoknot gap gets d5(s) + d3(t); a single base k gets s1w[k].
  bool dangle = false;
  std::vector<double> d5w, d3w, s1w;
  double run(int s, int t) const {
    if (t < s) return 0.0;
    if (s == t) return s1w[s];
    return d5w[s] + d3w[t];
  }

  bool pairs(int i, int j) const {
    return 0 <= i && i < j && j < n && cp[static_cast<size_t>(i) * n + j];
  }

  bool end_ok(int i, int j) const {
    if (!wc_ends) return true;
    int t = pair_type(code[i], code[j]);
    return t >= 0 && t != 3 && t != 5;
  }

  double hairpin(int i, int j) const {
    if (!pairs(i, j)) return NEG_INF;
    int u = j - i - 1;
    if (u < MIN_HAIRPIN) return NEG_INF;
    if (mat_rna1999) return log_weight(hairpin99(i, j));
    if (mat_v04) return log_weight(hairpin04(i, j) / 100.0);
    return log_weight(size_continuation(u, HAIRPIN_INIT, HAIRPIN_K));
  }

  double interior(int i, int d, int e, int j) const {
    if (!(pairs(i, j) && pairs(d, e))) return NEG_INF;
    if (!(i < d && d < e && e < j)) return NEG_INF;
    int l1 = d - i - 1;
    int l2 = j - e - 1;
    if (l1 + l2 > (mat_v04 ? V04_MAX_LOOP : max_loop) || l1 > max_left) return NEG_INF;
    if (mat_rna1999) return log_weight(interior99(i, d, e, j));
    if (mat_v04) return log_weight(interior04(i, j, d, e) / 100.0);
    if (l1 == 0 && l2 == 0) return stack(i, j);
    if (l1 == 0 || l2 == 0) return bulge(i, d, e, j, l1 + l2);
    return internal(l1, l2);
  }

  // interior() for consecutive pairs of a pseudoknot stem, scaled by st_mult / int_mult.
  double stem_interior(int i, int d, int e, int j) const {
    double w = interior(i, d, e, j);
    if (stem_mult && w > HALF_NEG) w *= (d == i + 1) ? st_mult : int_mult;
    return w;
  }

 private:
  // ---- rna1999 (kcal/mol), ene.c HairpinEnergy / InteriorEnergyFull / asymmetryEfn
  // code[k] is -1 for any non-ACGU base (N). Those bases cannot pair; loop lengths still
  // count them; every dangle / mismatch / special-hairpin term that would read an N is 0.
  bool okb(int k) const { return code[k] >= 0; }
  bool ok_range(int lo, int hi) const {
    for (int k = lo; k <= hi; ++k)
      if (code[k] < 0) return false;
    return true;
  }
  int b(int k) const { return code[k] + 1; }
  static double size_log(int size) { return 1.75 * RT * std::log(size / 30.0); }
  static double helix99(int a, int bb, int c, int d) {
    return T99.stack[mshift(a, bb) * 6 + mshift(c, d)];
  }
  static double interior_mm(int a, int bb, int x, int y) {
    return T99.mm_il[(4 * (x - 1) + (y - 1)) * 6 + mshift(a, bb)];
  }

  double hairpin99(int i, int j) const {
    const Rna1999Tables& T = T99;
    int size = j - i - 1;
    double g;
    if (size <= 30) {
      g = T.loop37[60 + size - 1];
    } else {
      g = T.loop37[60 + 30 - 1];
      g += size_log(size);
    }
    bool poly_c = crun[j] - crun[i + 1] == size;
    if (size == 3) {
      if (b(i) != 2 && b(j) != 2) g += T.at_penalty;
      if (ok_range(i, j)) {
        int idx = 0;
        for (int k = i; k <= j; ++k) idx = idx * 4 + code[k];
        g += T.tri[idx];
      }
      if (poly_c) g += T.polyc3;
    } else if (size == 4) {
      if (ok_range(i, j)) {
        int idx = 0;
        for (int k = i; k <= j; ++k) idx = idx * 4 + code[k];
        g += T.tetra[idx];
      }
      if (okb(i + 1) && okb(j - 1))
        g += T.mm_hp[(4 * (b(i + 1) - 1) + (b(j - 1) - 1)) * 6 + mshift(b(i), b(j))];
      if (poly_c) g += T.polyc_slope * size + T.polyc_int;
    } else {
      if (okb(i + 1) && okb(j - 1))
        g += T.mm_hp[(4 * (b(i + 1) - 1) + (b(j - 1) - 1)) * 6 + mshift(b(i), b(j))];
      if (poly_c) g += T.polyc_slope * size + T.polyc_int;
    }
    return g;
  }

  static double asymmetry99(int l1, int l2, int size) {
    const Rna1999Tables& T = T99;
    int asym = std::abs(l1 - l2);
    double g;
    if (size <= 30) {
      g = T.loop37[size - 1];
    } else {
      g = T.loop37[30 - 1];
      g += size_log(size);
    }
    int k = std::min(4, std::min(l1, l2));
    if (asym * T.asym[k - 1] < T.max_asym)
      g += asym * T.asym[k - 1];
    else
      g += T.max_asym;
    return g;
  }

  double interior99(int i, int d, int e, int j) const {
    const Rna1999Tables& T = T99;
    int l1 = d - i - 1;
    int l2 = j - e - 1;
    int size = l1 + l2;
    if (size == 0) return helix99(b(i), b(j), b(d), b(e));
    if (l1 == 0 || l2 == 0) {
      double g;
      if (size <= 30) {
        g = T.loop37[30 + size - 1];
      } else {
        g = T.loop37[30 + 30 - 1];
        g += size_log(size);
      }
      if (size == 1) {
        g += helix99(b(i), b(j), b(d), b(e));
      } else {
        if (b(i) != 2 && b(j) != 2) g += T.at_penalty;
        if (b(d) != 2 && b(e) != 2) g += T.at_penalty;
      }
      return g;
    }
    int asym = std::abs(l1 - l2);
    if (asym > 1 || size > 4) {
      double g = asymmetry99(l1, l2, size);
      if (l1 > 1 && l2 > 1) {
        if (okb(e + 1) && okb(d - 1)) g += interior_mm(b(e), b(d), b(e + 1), b(d - 1));
        if (okb(i + 1) && okb(j - 1)) g += interior_mm(b(i), b(j), b(i + 1), b(j - 1));
      } else {
        g += interior_mm(b(e), b(d), 1, 1);
        g += interior_mm(b(i), b(j), 1, 1);
      }
      return g;
    }
    if (asym == 0 && size == 2) {
      if (!(okb(i + 1) && okb(j - 1))) return 0.0;
      int cpi = mshift(b(i), b(j)), ipi = mshift(b(d), b(e));
      return 0.0 + T.sint2[96 * cpi + 16 * ipi + 4 * (b(i + 1) - 1) + (b(j - 1) - 1)];
    }
    if (asym == 0 && size == 4) {
      if (!(okb(i + 1) && okb(i + 2) && okb(j - 1) && okb(j - 2))) return 0.0;
      int cpi = mshift(b(i), b(j)), ipi = mshift(b(d), b(e));
      return 0.0 + T.sint4[cpi * 1536 + ipi * 256 + (4 * (b(i + 1) - 1) + b(j - 1) - 1) * 16 +
                           (4 * (b(i + 2) - 1) + b(j - 2) - 1)];
    }
    if (l1 == 1) {
      if (!(okb(i + 1) && okb(j - 1) && okb(j - 2))) return 0.0;
      int cpi = mshift(b(i), b(j)), ipi = mshift(b(d), b(e));
      return 0.0 + T.asint1x2[cpi * 384 + (b(j - 2) - 1) * 96 + (b(i + 1) - 1) * 24 + 4 * ipi +
                              ((b(j - 1) - 1) % 4)];
    }
    if (!(okb(i + 1) && okb(i + 2) && okb(j - 1))) return 0.0;
    int cpi = mshift(b(j), b(i)), ipi = mshift(b(e), b(d));
    return 0.0 + T.asint1x2[ipi * 384 + (b(i + 1) - 1) * 96 + (b(j - 1) - 1) * 24 + 4 * cpi +
                            ((b(i + 2) - 1) % 4)];
  }

  double stack(int i, int j) const {
    int a = pair_type(code[i], code[j]);
    int b = pair_type(code[i + 1], code[j - 1]);
    if (a < 0 || b < 0) return NEG_INF;
    return log_weight(STACK[a][b]);
  }

  double bulge(int i, int d, int e, int j, int nn) const {
    if (nn < 1) return NEG_INF;
    double g = size_continuation(nn, BULGE_INIT, BULGE_K);
    if (nn == 1) {
      int a = pair_type(code[i], code[j]);
      int b = pair_type(code[d], code[e]);
      g += (a >= 0 && b >= 0) ? STACK[a][b] : 0.0;
    }
    return log_weight(g);
  }

  double internal(int l1, int l2) const {
    int nn = l1 + l2;
    if (nn < 4) return NEG_INF;
    double g = size_continuation(nn, INTERIOR_INIT, INTERIOR_K);
    int asym = std::abs(l1 - l2);
    if (asym) g += std::min(3.0, 0.5 * asym);
    return log_weight(g);
  }

  // ---- vienna2004 (dcal/mol), LinearPartition-V's utility_v.h; bases A=1..U=4,
  // pair types CG=1 GC=2 GU=3 UG=4 AU=5 UA=6
  int vn(int k) const { return code[k] + 1; }

  int vtype(int i, int j) const {
    int a = code[i], c = code[j];
    if (a == 1 && c == 2) return 1;
    if (a == 2 && c == 1) return 2;
    if (a == 2 && c == 3) return 3;
    if (a == 3 && c == 2) return 4;
    if (a == 0 && c == 3) return 5;
    if (a == 3 && c == 0) return 6;
    return 0;
  }

  int vext(int size) const { return static_cast<int>(V04.lxc * std::log(size / 30.0)); }

  int hairpin04(int i, int j) const {
    const V04Tables& T = V04;
    int size = j - i - 1;
    int t = vtype(i, j);
    int e = size <= 30 ? T.hairpin[size] : T.hairpin[30] + vext(size);
    if (size == 4) {
      auto it = T.tetra.find(seq.substr(i, 6));
      if (it != T.tetra.end()) return it->second;
    } else if (size == 6) {
      auto it = T.hexa.find(seq.substr(i, 8));
      if (it != T.hexa.end()) return it->second;
    } else if (size == 3) {
      auto it = T.tri.find(seq.substr(i, 5));
      if (it != T.tri.end()) return it->second;
      return e + (t > 2 ? T.au : 0);
    }
    return e + T.mmH[t * 25 + vn(i + 1) * 5 + vn(j - 1)];
  }

  // Outer pair i·j, inner pair p·q (v_score_single).
  int interior04(int i, int j, int p, int q) const {
    const V04Tables& T = V04;
    int si1 = vn(i + 1), sj1 = vn(j - 1), sp1 = vn(p - 1), sq1 = vn(q + 1);
    int t = vtype(i, j);
    int t2 = vtype(q, p);
    int n1 = p - i - 1, n2 = j - q - 1;
    int nl = n1 > n2 ? n1 : n2;
    int ns = n1 > n2 ? n2 : n1;
    if (nl == 0) return T.stack[t * 8 + t2];
    if (ns == 0) {
      int e = nl <= 30 ? T.bulge[nl] : T.bulge[30] + vext(nl);
      if (nl == 1) {
        e += T.stack[t * 8 + t2];
      } else {
        if (t > 2) e += T.au;
        if (t2 > 2) e += T.au;
      }
      return e;
    }
    if (ns == 1) {
      if (nl == 1) return T.int11[((t * 8 + t2) * 5 + si1) * 5 + sj1];
      if (nl == 2) {
        if (n1 == 1) return T.int21[(((t * 8 + t2) * 5 + si1) * 5 + sq1) * 5 + sj1];
        return T.int21[(((t2 * 8 + t) * 5 + sq1) * 5 + si1) * 5 + sp1];
      }
      int e = nl + 1 <= 30 ? T.interior[nl + 1] : T.interior[30] + vext(nl + 1);
      e += std::min(T.max_ninio, (nl - ns) * T.ninio);
      return e + T.mm1n[t * 25 + si1 * 5 + sj1] + T.mm1n[t2 * 25 + sq1 * 5 + sp1];
    }
    if (ns == 2) {
      if (nl == 2) return T.int22[((((t * 8 + t2) * 5 + si1) * 5 + sp1) * 5 + sq1) * 5 + sj1];
      if (nl == 3)
        return T.interior[5] + T.ninio + T.mm23[t * 25 + si1 * 5 + sj1] +
               T.mm23[t2 * 25 + sq1 * 5 + sp1];
    }
    int u = nl + ns;
    int e = u <= 30 ? T.interior[u] : T.interior[30] + vext(u);
    e += std::min(T.max_ninio, (nl - ns) * T.ninio);
    return e + T.mmI[t * 25 + si1 * 5 + sj1] + T.mmI[t2 * 25 + sq1 * 5 + sp1];
  }

  // Pair i·j as an exterior-loop or multiloop branch, without ML_intern
  // (v_score_external_paired; E_MLstem uses the same values, mismatchM = mismatchExt).
  int branch04(int i, int j) const {
    const V04Tables& T = V04;
    int t = vtype(i, j);
    int e = t > 2 ? T.au : 0;
    if (dangles_mode == 3) {
      int si1 = i > 0 ? vn(i - 1) : -1;
      int sj1 = j < n - 1 ? vn(j + 1) : -1;
      if (si1 >= 0 && sj1 >= 0) e += T.mmExt[t * 25 + si1 * 5 + sj1];
      else if (si1 >= 0) e += T.d5[t * 5 + si1];
      else if (sj1 >= 0) e += T.d3[t * 5 + sj1];
    }
    return e;
  }

  // Pair i·j closing a multiloop, without ML_closing and ML_intern (v_score_multi).
  int closing04(int i, int j) const {
    int t = vtype(j, i);
    int e = t > 2 ? V04.au : 0;
    if (dangles_mode == 3) e += V04.mmM[t * 25 + vn(j - 1) * 5 + vn(i + 1)];
    return e;
  }
};

