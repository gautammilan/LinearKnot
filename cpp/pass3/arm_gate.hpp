// Legacy pass-3 ranking inputs, not used by the paper preset: the accessibility file
// (--arm-score, --access-file: AccessRec, AccessIndex) and the arm gate (--gate: GateCfg,
// Gate, GateStats).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

// ------------------------------------------ LinearAccess pre-pass (arm score, gate)
// Mirrors dp03/access.py.  Format 3 (natural log, 1-based): header line
// "#dp03-access-format 3 ...", then per sequence "S seq", "A a b logacc",
// "D a b y logQ", "P i j p" (every pair with p >= 1e-5).  Stored 0-based; a missing
// A line means logacc = ACC_LOG_FLOOR, missing D lines mean log D = -inf.
const char* const ACCESS_FORMAT_HEADER = "#dp03-access-format 3";
constexpr double ACC_LOG_FLOOR = -200.0;
constexpr int LK_HIST_BINS = 9;
constexpr double LK_HIST_EDGES[LK_HIST_BINS - 1] = {0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0};

// Same operations as access.log_add (true -inf, not NEG_INF).
inline double acc_log_add(double a, double b) {
  if (a == -INFINITY) return b;
  if (b == -INFINITY) return a;
  if (a > b) return a + std::log1p(std::exp(b - a));
  return b + std::log1p(std::exp(a - b));
}

struct AccWin {
  int t;
  std::vector<int> ys;       // innermost partners, ascending
  std::vector<double> lsuf;  // lsuf[k] = log sum_{m >= k} Q(ys[m])
};

struct AccessRec {
  int n = 0;
  std::unordered_map<uint64_t, double> logacc_tab, bpp;
  std::vector<std::vector<AccWin>> by_start;  // sorted by t

  double log_acc(int a, int b) const {
    auto it = logacc_tab.find(key2(a, b));
    return it == logacc_tab.end() ? ACC_LOG_FLOOR : it->second;
  }
  bool has_acc(int a, int b) const { return logacc_tab.count(key2(a, b)) != 0; }
  double pair_prob(int i, int j) const {
    auto it = bpp.find(key2(i, j));
    return it == bpp.end() ? 0.0 : it->second;
  }
  // The P lines as (i, j, p), 0-based, sorted by (i, j), one entry per pair.
  std::vector<gate::Pair> pairs() const {
    std::vector<gate::Pair> out;
    out.reserve(bpp.size());
    for (const auto& kv : bpp)
      out.push_back({static_cast<int>(kv.first >> 32), static_cast<int>(kv.first & 0xffffffffu),
                     kv.second});
    std::sort(out.begin(), out.end(), [](const gate::Pair& a, const gate::Pair& b) {
      return a.i != b.i ? a.i < b.i : a.j < b.j;
    });
    return out;
  }
};

// sequence -> byte span [start, end) of its record; the file stays open for seeks.
struct AccessIndex {
  FILE* f = nullptr;
  std::string path;
  std::unordered_map<std::string, std::pair<long, long>> span;
};

bool read_line(FILE* f, std::string& line) {
  line.clear();
  int c;
  while ((c = std::fgetc(f)) != EOF) {
    line.push_back(static_cast<char>(c));
    if (c == '\n') break;
  }
  return !line.empty();
}

AccessIndex index_access_file(const std::string& path) {
  AccessIndex idx;
  idx.path = path;
  idx.f = std::fopen(path.c_str(), "rb");
  if (!idx.f) {
    std::fprintf(stderr, "cannot open access file %s\n", path.c_str());
    std::exit(2);
  }
  std::string line;
  read_line(idx.f, line);
  if (line.compare(0, std::strlen(ACCESS_FORMAT_HEADER), ACCESS_FORMAT_HEADER) != 0) {
    std::fprintf(stderr,
                 "%s: not an access file of format 3 (older format?); "
                 "rebuild it with access/access_prepass.py\n",
                 path.c_str());
    std::exit(2);
  }
  std::string cur;
  bool have = false;
  long start = 0, pos = std::ftell(idx.f);
  while (read_line(idx.f, line)) {
    if (line.size() >= 2 && line[0] == 'S' && line[1] == ' ') {
      if (have) idx.span.emplace(cur, std::make_pair(start, pos));
      size_t b = 2, e = line.size();
      while (b < e && std::isspace(static_cast<unsigned char>(line[b]))) ++b;
      while (e > b && std::isspace(static_cast<unsigned char>(line[e - 1]))) --e;
      cur = line.substr(b, e - b);
      have = true;
      start = pos;
    }
    pos += static_cast<long>(line.size());
  }
  if (have) idx.span.emplace(cur, std::make_pair(start, pos));
  return idx;
}

AccessRec load_access(const AccessIndex& idx, const std::string& seq) {
  AccessRec rec;
  const int n = static_cast<int>(seq.size());
  rec.n = n;
  rec.by_start.assign(n, {});
  auto sp = idx.span.find(seq);
  if (sp == idx.span.end()) {
    std::fprintf(stderr, "[arm] warning: no access record for sequence of length %d in %s; D = 0\n", n,
                 idx.path.c_str());
    return rec;
  }
  std::fseek(idx.f, sp->second.first, SEEK_SET);
  std::string block(static_cast<size_t>(sp->second.second - sp->second.first), '\0');
  if (std::fread(&block[0], 1, block.size(), idx.f) != block.size()) {
    std::fprintf(stderr, "short read in access file %s\n", idx.path.c_str());
    std::exit(2);
  }
  std::unordered_map<uint64_t, std::vector<std::pair<int, double>>> qs;
  std::vector<uint64_t> order;
  size_t p0 = 0;
  while (p0 < block.size()) {
    size_t p1 = block.find('\n', p0);
    if (p1 == std::string::npos) p1 = block.size();
    const std::string line = block.substr(p0, p1 - p0);
    const char* ln = line.c_str();
    char tag[8], x1[64], x2[64], x3[64], x4[64], x5[8];
    int k = std::sscanf(ln, "%7s %63s %63s %63s %63s %7s", tag, x1, x2, x3, x4, x5);
    if (k == 4 && std::strcmp(tag, "A") == 0) {
      int a = std::atoi(x1) - 1, b = std::atoi(x2) - 1;
      if (0 <= a && a <= b && b < n) rec.logacc_tab[key2(a, b)] = std::strtod(x3, nullptr);
    } else if (k == 4 && std::strcmp(tag, "P") == 0) {
      int i = std::atoi(x1) - 1, j = std::atoi(x2) - 1;
      if (0 <= i && i < j && j < n) rec.bpp[key2(i, j)] = std::strtod(x3, nullptr);
    } else if (k == 5 && std::strcmp(tag, "D") == 0) {
      int a0 = std::atoi(x1) - 1, b0 = std::atoi(x2) - 1, y0 = std::atoi(x3) - 1;
      double lq = std::strtod(x4, nullptr);
      if (0 <= a0 && a0 <= b0 && b0 < y0 && y0 < n) {
        uint64_t key = key2(a0, b0);
        auto it = qs.find(key);
        if (it == qs.end()) {
          order.push_back(key);
          it = qs.emplace(key, std::vector<std::pair<int, double>>{}).first;
        }
        bool found = false;
        for (auto& e : it->second)
          if (e.first == y0) {
            e.second = acc_log_add(e.second, lq);
            found = true;
            break;
          }
        if (!found) it->second.emplace_back(y0, acc_log_add(-INFINITY, lq));
      }
    }
    p0 = p1 + 1;
  }
  for (uint64_t key : order) {
    auto slot = qs[key];
    std::sort(slot.begin(), slot.end());
    AccWin w;
    w.t = static_cast<int>(key & 0xffffffffu);
    w.ys.resize(slot.size());
    w.lsuf.assign(slot.size(), -INFINITY);
    double run = -INFINITY;
    for (int m = static_cast<int>(slot.size()) - 1; m >= 0; --m) {
      run = acc_log_add(run, slot[m].second);
      w.ys[m] = slot[m].first;
      w.lsuf[m] = run;
    }
    rec.by_start[static_cast<int>(key >> 32)].push_back(std::move(w));
  }
  for (auto& row : rec.by_start)
    std::stable_sort(row.begin(), row.end(),
                     [](const AccWin& x, const AccWin& y) { return x.t < y.t; });
  return rec;
}

struct ArmStats {
  long long scored = 0, arm_wins = 0, windows = 0, floor_windows = 0, chosen_floor = 0;
  double sum_log_kappa = 0.0, max_log_kappa = 0.0;
  long long lk_hist[LK_HIST_BINS] = {};
};

struct ArmCand {
  double a, z;
  bool floor;
};

// ------------------------------------------------------- eligibility gate (Qg / Qgl)
// A four-index state with outer pair p.q and gap [u+1, v-1] gets
//   g1 = min(1, max_{x in gap} r(x, q) / max(delta, 1 - P(p, q)))       (stem-1 role)
//   g2 = the same on the mirrored pairs, partners left of p              (stem-2 role)
//   g  = g1 (role stem1) or max(g1, g2) (role both); an empty gap gives g = 0,
// with r(x, q) = sum of P(x, y) over q < y <= q + span (gate_index.hpp).  At every
// pruning event of a gated table, states with g < tau are dropped before ranking; the
// rest are ranked by the unchanged pruning score, and `reserve` more eligible states
// with the highest g may be added.  With tau > 0 states are removed even at infinite
// beam, so infinite-beam exactness holds only for tau = 0.
enum GateTable { GATE_QG = 0, GATE_QGL = 1 };

struct GateCfg {
  bool on = false;
  double tau = 0.01;
  int span = 0;
  double delta = 1e-3;
  bool both = true;
  bool tables[2] = {true, true};
  int reserve = 0;
};

class Gate {
 public:
  Gate(int n, const AccessRec& acc, int span, double delta, bool both)
      : n_(n), acc_(acc), delta_(delta), both_(both) {
    std::vector<gate::Pair> pairs = acc.pairs();
    fwd_ = gate::PersistentGate(n, pairs, span);
    mir_ = gate::PersistentGate(n, gate::mirror_pairs(n, pairs), span);
  }

  void roles(int p, int u, int v, int q, double& g1, double& g2) const {
    const int lo = u + 1, hi = v - 1;
    if (lo > hi) {
      g1 = g2 = 0.0;
      return;
    }
    const double denom = std::max(delta_, 1.0 - acc_.pair_prob(p, q));
    g1 = std::min(1.0, fwd_.range_max(q, lo, hi) / denom);
    g2 = std::min(1.0, mir_.range_max(n_ - 1 - p, n_ - 1 - hi, n_ - 1 - lo) / denom);
  }

  double g(int p, int u, int v, int q) const {
    double g1, g2;
    roles(p, u, v, q, g1, g2);
    return both_ ? std::max(g1, g2) : g1;
  }

 private:
  int n_;
  const AccessRec& acc_;
  double delta_;
  bool both_;
  gate::PersistentGate fwd_, mir_;
};

// Per gated table, summed over pruning events: candidates seen, eligible (g >= tau),
// ranked (scored by the pruning score), kept, and kept through the reserve.
struct GateStats {
  long long cand[2] = {}, elig[2] = {}, ranked[2] = {}, kept[2] = {}, reserved[2] = {};
};

