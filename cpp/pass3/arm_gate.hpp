// Per-table candidate counters of the Qg / Qgl rankings (GateStats; counted at every pruning
// event, not reported).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

enum GateTable { GATE_QG = 0, GATE_QGL = 1 };

// Per table, summed over pruning events: candidates seen, eligible, ranked (scored by the
// pruning score), kept.
struct GateStats {
  long long cand[2] = {}, elig[2] = {}, ranked[2] = {}, kept[2] = {};
};
