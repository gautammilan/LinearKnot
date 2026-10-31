// Beam bookkeeping shared by all passes: table ids BT_* (Qb, Qm, Qm2, Qmp, PKM, Qp, Qg, Qgls,
// Qgl, Qgr, ...).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.


// Beamed tables (per-table beams).
enum BeamTable { BT_QB, BT_QM, BT_QM2, BT_QMP, BT_PKM, BT_QP, BT_QG, BT_QGLS, BT_QGL, BT_QGR,
                 BT_QGRS, BT_N };
const char* const BT_NAMES[BT_N] = {"Qb", "Qm",  "Qm2", "Qmp", "PKM", "Qp",
                                    "Qg", "Qgls", "Qgl", "Qgr", "Qgrs"};


