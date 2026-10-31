// Log-space arithmetic: NEG_INF, log_add(a, b) = log(e^a + e^b); RT (kT in kcal/mol).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.


constexpr double NEG_INF = -1e300;
constexpr double HALF_NEG = NEG_INF / 2;

// Binary exponent of a (|a| in [2^e, 2^(e+1))): field bits, no branch, no libm. Zero and
// subnormals give e = -1023, inf/NaN give 1024; both only make the bound below looser.
inline int log_add_exp2(double a) {
  std::uint64_t u;
  std::memcpy(&u, &a, sizeof u);
  return static_cast<int>((u >> 52) & 0x7ff) - 1023;
}

// log_add(a, b) = log(e^a + e^b) for log-space weights.
//
// Two early outs before the transcendental pair, both of which return a bit for bit:
//   b - a <= -50 is the cheap one (exp underflows into the rounding of a anyway);
//   b - a < (e_a - 55) ln 2 is the ulp one: a quarter ulp of a is 2^(e_a - 55), and adding
//   anything below a quarter ulp to a rounds back to a for every rounding mode we build with
//   (log1p(x) < 1.01x for x < 1e-16, so the term stays under a quarter ulp).
// The comparison is done in the d/ln 2 form to keep it to one multiply; the factor 4 margin
// between the quarter-ulp bound and the half-ulp rounding threshold absorbs its rounding.
// -inf, +inf, 0, subnormals and NaN all fall through to the slow path unchanged.
inline double log_add(double a, double b) {
  if (a < b) std::swap(a, b);
  if (a < NEG_INF / 2) return NEG_INF;
  const double d = b - a;  // <= 0
  if (d < -50.0) return a;
  if (d * 1.4426950408889634074 + 55.0 < static_cast<double>(log_add_exp2(a))) return a;
  return a + std::log1p(std::exp(d));
}

// ------------------------------------------------------------------ energy
double RT = 0.616330318;  // --rt overrides (e.g. NUPACK's kB*T at 37 C)

