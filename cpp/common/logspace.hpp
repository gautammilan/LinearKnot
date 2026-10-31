// Log-space arithmetic: NEG_INF, log_add(a, b) = log(e^a + e^b); RT (kT in kcal/mol).
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

// ------------------------------------------------------------------ energy
double RT = 0.616330318;  // --rt overrides (e.g. NUPACK's kB*T at 37 C)

