// Log-space arithmetic: NEG_INF, log_add(a, b) = log(e^a + e^b); RT (kT in kcal/mol);
// the optional per-seed mass accumulator (--mass-log).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.


constexpr double NEG_INF = -1e300;
constexpr double HALF_NEG = NEG_INF / 2;

inline double log_add(double a, double b) {
  if (a < b) std::swap(a, b);
  if (a < NEG_INF / 2) return NEG_INF;
  if (a - b > 50.0) return a;
  return a + std::log1p(std::exp(b - a));
}
// --g-mass: log-sum of each R-key growth score, per (kind, seed, p, j). Arm is summed in.
static std::unordered_map<uint64_t, double> g_mass_acc;
static void g_mass_note(int kind, int seed, int p, int j, double s) {
  if (!g_mass_f || s < NEG_INF / 2) return;
  const uint64_t k = (static_cast<uint64_t>(kind) << 56) | (static_cast<uint64_t>(seed) << 48) |
                      (static_cast<uint64_t>(p) << 24) | static_cast<uint64_t>(j & 0xffffff);
  auto it = g_mass_acc.find(k);
  if (it == g_mass_acc.end()) g_mass_acc.emplace(k, s);
  else it->second = log_add(it->second, s);
}
static void g_mass_flush() {
  if (!g_mass_f) return;
  for (const auto& kv : g_mass_acc) {
    const int rec[5] = {static_cast<int>(kv.first >> 56), static_cast<int>((kv.first >> 48) & 0xff),
                         static_cast<int>((kv.first >> 24) & 0xffffff),
                         static_cast<int>(kv.first & 0xffffff), -1};
    const double s = kv.second;
    std::fwrite(rec, sizeof rec, 1, g_mass_f);
    std::fwrite(&s, sizeof s, 1, g_mass_f);
  }
  std::fflush(g_mass_f);
  g_mass_acc.clear();
}

// ------------------------------------------------------------------ energy
double RT = 0.616330318;  // --rt overrides (e.g. NUPACK's kB*T at 37 C)

