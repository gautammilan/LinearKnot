// Timers: cpu_sec (clock), ru_cpu (user + sys from getrusage), hwm_mb (peak RSS).
//
// Build with -DLK_FINE_TIMERS to enable the sub-stage clocks (SecTimer, PrGrowTimer,
// PrPartnerTimer, the K-read and rec_key charges). They read steady_clock twice per call and
// sit in per-column loops, so without the macro they are compiled out and their *_sec values
// stay 0. Only the per-pass `[stage]` lines (pass1 / pass2 / scan / outside, from getrusage)
// are always on.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

// Process CPU seconds and peak RSS so far.
static double cpu_sec() { return static_cast<double>(std::clock()) / CLOCKS_PER_SEC; }
static double ru_cpu() {
  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru) != 0) return 0.0;
  return static_cast<double>(ru.ru_utime.tv_sec) + ru.ru_utime.tv_usec * 1e-6 +
         static_cast<double>(ru.ru_stime.tv_sec) + ru.ru_stime.tv_usec * 1e-6;
}
static double hwm_mb() {
  std::ifstream f("/proc/self/status");
  std::string k;
  while (f >> k)
    if (k == "VmHWM:") {
      double kb = 0.0;
      f >> kb;
      return kb / 1024.0;
    }
  return 0.0;
}

