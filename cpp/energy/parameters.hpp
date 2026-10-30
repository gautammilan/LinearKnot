// Energy parameters. PkParams: the DP03 pseudoknot penalties (beta1, beta2, beta3, ...) for
// --params dp03 | dp09*; nested-loop tables (stack, hairpin, bulge, interior); the NUPACK
// rna1999.dG loader (--material rna1999) and the ViennaRNA 2004 loader (--material vienna2004).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

// Multiloop and pseudoknot penalties (kcal/mol), selected with --params.
//   stem:      per pseudoknot stem (Qgl); NUPACK's beta2, HotKnots' Pb (per band)
//   nest:      per nested pair inside a gap (Qz); beta2 / HotKnots' Pps
//   nest_pk_k: a pseudoknot nested in a gap pays nest_pk_k * nest (+ beta1p)
//   ml_pk_k:   a pseudoknot as a multiloop branch pays ml_pk_k * alpha2 (+ beta1m)
//   st_mult, int_mult: scale the energy of a stack / interior loop inside a pseudoknot
//              stem (Qg (a), (b)); HotKnots' stP when the left side has no unpaired base,
//              else intP
//   a_p, b_p, c_p: multiloop inside a pseudoknot stem (Qg (f)-(h), Qgls, Qgrs);
//              stem_ml gives its branches their own Qm table (Qmp) so b_p, c_p apply there
struct PkParams {
  double alpha1, alpha2, alpha3, beta1, beta1m, beta1p, stem, nest, beta3;
  double nest_pk_k, ml_pk_k;
  double st_mult, int_mult, a_p, b_p, c_p;
  bool stem_ml;
};
// Dirks & Pierce 2003, as in NUPACK 3.
constexpr PkParams DP03_PARAMS{3.4, 0.4, 0.0, 9.6, 15.0, 15.0, 0.1, 0.1, 0.1, 2.0, 2.0,
                               1.0, 1.0, 3.4, 0.4, 0.0, false};
// Andronescu et al. 2010 (HotKnots 2.0 params/parameters_DP09.txt, pkmodelDP a, b, c,
// Ps, Psm, Psp, Pb, Pps, Pup, stP, intP, a_p, b_p, c_p); HotKnots counts a nested
// pseudoknot as one closed region and a pseudoknot branch as one multiloop branch.
// dp09 = penalties only; dp09st adds stP/intP; dp09ml adds a_p/b_p/c_p; dp09full both.
constexpr PkParams DP09_PARAMS{3.39, 0.03, 0.02, -1.38, 10.07, 15.00, 2.46, 0.96, 0.06, 1.0, 1.0,
                               1.0, 1.0, 3.39, 0.03, 0.02, false};
constexpr PkParams DP09ST_PARAMS{3.39, 0.03, 0.02, -1.38, 10.07, 15.00, 2.46, 0.96, 0.06, 1.0, 1.0,
                                 0.89, 0.74, 3.39, 0.03, 0.02, false};
constexpr PkParams DP09ML_PARAMS{3.39, 0.03, 0.02, -1.38, 10.07, 15.00, 2.46, 0.96, 0.06, 1.0, 1.0,
                                 1.0, 1.0, 3.41, 0.56, 0.12, true};
constexpr PkParams DP09FULL_PARAMS{3.39, 0.03, 0.02, -1.38, 10.07, 15.00, 2.46, 0.96, 0.06, 1.0,
                                   1.0, 0.89, 0.74, 3.41, 0.56, 0.12, true};
PkParams PK = DP03_PARAMS;

constexpr int MIN_HAIRPIN = 3;
constexpr int MIN_PAIR = MIN_HAIRPIN + 1;
constexpr int MAX_LOOP = 30;
// Interior-loop size cap (l1 + l2). Python's energy.py hard-codes MAX_LOOP;
// raise it only to compare against uncapped NUPACK.
int max_loop = MAX_LOOP;
bool allow_gu = true;
// One-pair pseudoknot stems (Qg_{i,i,j,j} = 1), stored as Qg[j][j][(i,i)].
bool single_stems = true;
// NUPACK 3: pairs bordering exterior / multi / Qz regions must be Watson-Crick.
bool wc_ends = false;
double pair_bonus = 0.0;
// Multiloop leading unpaired run, every pass: closing left p and first branch i
// satisfy i - p - 1 <= ml_lead_cap - 1 (29 when the flag is 30, p >= i - 30).
// 0 turns the cap off. Interior loops, the trailing run, between-branch runs,
// the exterior, and Q^z stay on their own rules.
int ml_lead_cap = 30;
const std::vector<std::pair<int, int>>* allowed_pairs = nullptr;  // 0-based, nullptr = all

// pair types: AU CG GC GU UA UG
inline int base_code(char c) {
  switch (c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'U': return 3;
    default: return -1;
  }
}

inline int pair_type(int a, int b) {
  static const int T[4][4] = {
      {-1, -1, -1, 0},
      {-1, -1, 1, -1},
      {-1, 2, -1, 3},
      {4, -1, 5, -1},
  };
  if (a < 0 || b < 0) return -1;
  return T[a][b];
}

const double STACK[6][6] = {
    {-0.9, -2.2, -2.1, -0.6, -1.1, -1.4},
    {-2.1, -3.3, -2.4, -1.4, -2.1, -2.1},
    {-2.4, -3.4, -3.3, -1.5, -2.2, -2.5},
    {-1.3, -2.5, -2.1, -0.5, -1.4, 1.3},
    {-1.3, -2.4, -2.1, -1.0, -0.9, -1.3},
    {-1.0, -1.5, -1.4, -0.3, -0.6, -0.5},
};

const double HAIRPIN_INIT[21] = {0,   0,   0,   5.4, 5.6, 5.7, 5.4,
                                 6.0, 5.5, 6.4, 6.5, 6.6, 6.7, 6.8,
                                 6.9, 6.9, 7.0, 7.1, 7.1, 7.2, 7.2};
constexpr int HAIRPIN_K = 20;
const double BULGE_INIT[7] = {0, 3.8, 2.8, 3.2, 3.6, 4.0, 4.4};
constexpr int BULGE_K = 6;
const double INTERIOR_INIT[21] = {0,   0,   0,   0,   1.1, 2.0, 2.0,
                                  2.1, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8,
                                  2.9, 2.9, 3.0, 3.1, 3.1, 3.2, 3.3};
constexpr int INTERIOR_K = 20;

inline double log_weight(double g) { return -g / RT; }

inline double size_continuation(int n, const double* table, int k) {
  if (n <= k) return table[n];
  return table[k] + 1.75 * RT * std::log(static_cast<double>(n) / k);
}

// --material rna1999: NUPACK 3.0.6's Turner 1999 loop energies (dp03/rna1999.py) read from
// --param-file, in NUPACK's table layout (init.c LoadEnergies), plus its terminal AU/GU
// penalty on pairs bordering an exterior loop, multiloop or pseudoknot gap.  Bases in the
// NUPACK formulas are coded A=1, C=2, G=3, U=4.
bool mat_rna1999 = false;
std::string param_file;
// --dangles: 0 none; rna1999: 1 some, 2 all; vienna2004: 3 d2
int dangles_mode = 0;

struct Rna1999Tables {
  double stack[36], loop37[90], asym[4], max_asym;
  double tri[1024] = {}, tetra[4096] = {};  // base-4 index of the loop, A=0
  double mm_hp[96], mm_il[96], dangle[48];
  double alpha[3], at_penalty;
  double sint2[576], sint4[9216], asint1x2[2304];
  double polyc3, polyc_slope, polyc_int, beta[5];
};
Rna1999Tables T99;

[[noreturn]] void param_error(const std::string& msg) {
  std::fprintf(stderr, "%s: %s\n", param_file.c_str(), msg.c_str());
  std::exit(2);
}

std::vector<int> line_ints(const std::string& line) {
  std::vector<int> out;
  size_t k = 0;
  while (k < line.size()) {
    while (k < line.size() && std::isspace(static_cast<unsigned char>(line[k]))) ++k;
    size_t b = k;
    while (k < line.size() && !std::isspace(static_cast<unsigned char>(line[k]))) ++k;
    if (b == k) break;
    std::string t = line.substr(b, k - b);
    size_t s = (t[0] == '-' || t[0] == '+') ? 1 : 0;
    if (s < t.size() && t.find_first_not_of("0123456789", s) == std::string::npos)
      out.push_back(std::atoi(t.c_str()));
  }
  return out;
}

int loop_index(const std::string& loop) {
  int idx = 0;
  for (char c : loop) {
    int b = base_code(c == 'T' ? 'U' : c);
    if (b < 0) param_error("bad loop sequence " + loop);
    idx = idx * 4 + b;
  }
  return idx;
}

void load_rna1999(const std::string& path) {
  FILE* f = std::fopen(path.c_str(), "r");
  if (!f) param_error("cannot open (use --param-file)");
  std::vector<std::vector<std::string>> secs;
  std::string line;
  bool open = false;
  int c;
  while (true) {
    line.clear();
    while ((c = std::fgetc(f)) != EOF && c != '\n') line.push_back(static_cast<char>(c));
    if (c == EOF && line.empty()) break;
    if (!line.empty() && line[0] == '>') {
      open = false;
    } else {
      if (!open) secs.emplace_back();
      open = true;
      secs.back().push_back(line);
    }
    if (c == EOF) break;
  }
  std::fclose(f);
  size_t si = 0;
  auto next = [&]() -> const std::vector<std::string>& {
    if (si >= secs.size()) param_error("too few sections");
    return secs[si++];
  };
  auto rows = [&](const std::vector<std::string>& sec, size_t k, size_t width, double* out) {
    if (sec.size() != k) param_error("bad section");
    for (size_t r = 0; r < k; ++r) {
      std::vector<int> v = line_ints(sec[r]);
      if (v.size() != width) param_error("bad row " + sec[r]);
      for (size_t q = 0; q < width; ++q) out[r * width + q] = v[q] / 100.0;
    }
  };
  Rna1999Tables& T = T99;
  rows(next(), 6, 6, T.stack);
  for (int i = 0; i < 3; ++i) {
    std::vector<int> v = line_ints(next()[0]);
    if (v.size() != 30) param_error("loop rows must list sizes 1..30");
    for (int j = 0; j < 30; ++j) T.loop37[30 * (2 - i) + j] = v[j] / 100.0;
  }
  {
    std::vector<int> v = line_ints(next()[0]);
    for (int k = 0; k < 4; ++k) T.asym[k] = v[k] / 100.0;
    T.max_asym = v[4] / 100.0;
  }
  for (int which = 0; which < 2; ++which) {
    double* bonus = which == 0 ? T.tri : T.tetra;
    size_t len = which == 0 ? 5 : 6;
    for (const std::string& s : next()) {
      size_t b = s.find_first_not_of(" \t");
      size_t e = b == std::string::npos ? b : s.find_first_of(" \t", b);
      std::vector<int> v = e == std::string::npos ? std::vector<int>() : line_ints(s.substr(e));
      if (e == std::string::npos || e - b != len || v.size() != 1)
        param_error("bad loop bonus " + s);
      bonus[loop_index(s.substr(b, e - b))] = v[0] / 100.0;
    }
  }
  rows(next(), 16, 6, T.mm_hp);
  rows(next(), 16, 6, T.mm_il);
  rows(next(), 6, 4, T.dangle);
  rows(next(), 6, 4, T.dangle + 24);
  {
    std::vector<int> v = line_ints(next()[0]);
    for (int k = 0; k < 3; ++k) T.alpha[k] = v[k] / 100.0;
    T.at_penalty = line_ints(next()[0])[0] / 100.0;
  }
  auto blocks = [&](size_t k, auto&& put) {
    const std::vector<std::string>& sec = next();
    if (sec.size() != 5 * k) param_error("bad loop table");
    for (size_t b = 0; b < k; ++b)
      for (int j = 0; j < 4; ++j) {
        std::vector<int> v = line_ints(sec[5 * b + 1 + j]);
        if (v.size() != 4) param_error("bad table row " + sec[5 * b + 1 + j]);
        for (int q = 0; q < 4; ++q) put(static_cast<int>(b), j, q, v[q] / 100.0);
      }
  };
  blocks(36, [&](int i, int j, int k, double x) { T.sint2[i * 16 + j * 4 + k] = x; });
  blocks(36 * 16, [&](int i, int j, int k, double x) {
    T.sint4[1536 * (i / 96) + 256 * ((i % 96) / 16) + 64 * ((i % 16) / 4) + 4 * (i % 4) +
            k * 16 + j] = x;
  });
  blocks(144, [&](int i, int j, int k, double x) {
    T.asint1x2[384 * (i / 24) + 4 * ((i % 24) / 4) + 24 * (i % 4) + 96 * j + k] = x;
  });
  {
    std::vector<int> v = line_ints(next()[0]);
    T.polyc3 = v[0] / 100.0;
    T.polyc_slope = v[1] / 100.0;
    T.polyc_int = v[2] / 100.0;
    v = line_ints(next()[0]);
    for (int k = 0; k < 5; ++k) T.beta[k] = v[k] / 100.0;
  }
}

// --material vienna2004: LinearPartition-V's Turner 2004 model (dp03/vienna2004.py), nested only,
// read from --param-file (dp03/params/vienna2004.par, written by dp03/make_vienna2004.py);
// energies in dcal/mol.  LinearPartition-V's loop limits: interior loops of at most 70 unpaired,
// 69 on the left; at most 69 unpaired before a multiloop's first branch.
bool mat_v04 = false;
constexpr int V04_MAX_LOOP = 70;
constexpr int V04_MAX_LEFT = 69;
constexpr int NO_LIMIT = 1 << 30;

struct V04Tables {
  double lxc = 0.0;
  int ml_intern = 0, ml_closing = 0, ml_base = 0, max_ninio = 0, ninio = 0, au = 0;
  std::vector<int> stack, hairpin, bulge, interior, mmI, mmH, mmM, mm1n, mm23, mmExt, d5, d3,
      int11, int21, int22;
  std::unordered_map<std::string, int> tri, tetra, hexa;
};
V04Tables V04;

void load_v04(const std::string& path) {
  FILE* f = std::fopen(path.c_str(), "r");
  if (!f) param_error("cannot open (use --param-file)");
  int c;
  while ((c = std::fgetc(f)) != EOF && c != '\n') {}
  auto tok = [&]() -> std::string {
    std::string t;
    while ((c = std::fgetc(f)) != EOF && std::isspace(c)) {}
    while (c != EOF && !std::isspace(c)) {
      t.push_back(static_cast<char>(c));
      c = std::fgetc(f);
    }
    return t;
  };
  V04Tables& T = V04;
  std::unordered_map<std::string, int*> scalars = {
      {"ML_intern37", &T.ml_intern}, {"ML_closing37", &T.ml_closing}, {"ML_BASE37", &T.ml_base},
      {"MAX_NINIO", &T.max_ninio},   {"ninio37", &T.ninio},           {"TerminalAU37", &T.au}};
  std::unordered_map<std::string, std::vector<int>*> arrays = {
      {"stack37", &T.stack},       {"hairpin37", &T.hairpin},    {"bulge37", &T.bulge},
      {"internal_loop37", &T.interior}, {"mismatchI37", &T.mmI}, {"mismatchH37", &T.mmH},
      {"mismatchM37", &T.mmM},     {"mismatch1nI37", &T.mm1n},   {"mismatch23I37", &T.mm23},
      {"mismatchExt37", &T.mmExt}, {"dangle5_37", &T.d5},        {"dangle3_37", &T.d3},
      {"int11_37", &T.int11},      {"int21_37", &T.int21},       {"int22_37", &T.int22}};
  bool has_lxc = false;
  size_t n_scalars = 0;
  while (true) {
    std::string kind = tok();
    if (kind.empty()) break;
    std::string name = tok();
    if (kind == "scalar") {
      std::string v = tok();
      if (name == "lxc37") {
        T.lxc = std::strtod(v.c_str(), nullptr);
        has_lxc = true;
      } else {
        auto it = scalars.find(name);
        if (it == scalars.end()) param_error("unknown scalar " + name);
        *it->second = std::atoi(v.c_str());
        ++n_scalars;
      }
    } else if (kind == "array") {
      auto it = arrays.find(name);
      if (it == arrays.end()) param_error("unknown array " + name);
      int nd = std::atoi(tok().c_str());
      size_t size = 1;
      for (int k = 0; k < nd; ++k) size *= static_cast<size_t>(std::atoi(tok().c_str()));
      it->second->resize(size);
      for (size_t q = 0; q < size; ++q) (*it->second)[q] = std::atoi(tok().c_str());
    } else if (kind == "loops") {
      std::unordered_map<std::string, int>& m =
          name == "Triloops" ? T.tri : name == "Tetraloops" ? T.tetra : T.hexa;
      if (name != "Triloops" && name != "Tetraloops" && name != "Hexaloops")
        param_error("unknown loops " + name);
      int cnt = std::atoi(tok().c_str());
      for (int k = 0; k < cnt; ++k) {
        std::string s = tok();
        m[s] = std::atoi(tok().c_str());
      }
    } else {
      param_error("bad token " + kind);
    }
  }
  std::fclose(f);
  if (!has_lxc || n_scalars != scalars.size()) param_error("missing scalars");
  for (const auto& a : arrays)
    if (a.second->empty()) param_error("missing array " + a.first);
}

inline int mshift(int b1, int b2) {  // GetMismatchShift: AU CG GC UA GU UG
  if (b1 == 1) return 0;
  if (b1 == 2) return 1;
  if (b1 == 3) return b2 == 2 ? 2 : 4;
  return b2 == 1 ? 3 : 5;
}

