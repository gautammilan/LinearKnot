// FastMap: the hash map type of the index maps on the hot paths.
//
// `std::unordered_map` chains: every lookup hashes, indexes a bucket array, then follows a
// pointer to a node that may be anywhere in memory. ankerl::unordered_dense is open addressing
// over flat arrays with the 64-bit hash stored next to the key, so a hit is one random access.
// Same interface (find / emplace / operator[] / clear / size / iteration over .first/.second),
// so a swap is a one-line change of declaration.
//
// Only use it for maps that are lookup-only or whose iteration order never reaches a
// floating-point sum: iteration order is the vector order here and was unspecified before, so
// anything that accumulates over a map in iteration order would change results bit for bit.
//
// Include at global scope (cpp/main.cpp does, before the anonymous namespace): the included
// third-party header opens ankerl::unordered_dense, which must not land in a local namespace.

#include "../third_party/ankerl_unordered_dense.h"

template <class K, class V>
using FastMap = ankerl::unordered_dense::map<K, V>;
