// Memory accounting: MemItem / msz (bytes per structure, for the [mem-nt] line).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

// Per structure entries and bytes.  Bytes: container capacities plus 16 per heap block and the
// hash nodes / buckets of the unordered maps (libstdc++ layout, no cached hash).
struct MemItem {
  std::string name;
  long long ents;
  double bytes;
};
namespace msz {
inline double node(size_t payload) {
  return static_cast<double>(std::max<size_t>(32, (payload + 16 + 15) & ~static_cast<size_t>(15)));
}
template <class T>
double vec(const std::vector<T>& v) {
  return v.capacity() ? static_cast<double>(v.capacity() * sizeof(T) + 16) : 0.0;
}
template <class T>
double vec2(const std::vector<std::vector<T>>& v) {
  double b = vec(v);
  for (const auto& x : v) b += vec(x);
  return b;
}
template <class M>
double umap(const M& m) {
  return static_cast<double>(m.bucket_count() * sizeof(void*)) +
         static_cast<double>(m.size()) * node(sizeof(typename M::value_type));
}
inline long long fin(const std::vector<double>& v) {
  long long c = 0;
  for (double x : v) c += x > HALF_NEG;
  return c;
}
}  // namespace msz
inline void stem_mem_items(const StemEst& s, std::vector<MemItem>& out) {
  const size_t se = s.entries + s.entries2 + s.entries_av;
  out.push_back({"stem_idx", static_cast<long long>(se),
                 12.0 * static_cast<double>(se) +
                     16.0 * static_cast<double>(s.index.nodes() + s.index2.nodes() + s.index_av.nodes())});
  out.push_back({"stem_S_C0", static_cast<long long>(s.S.size() + s.C0.size()), msz::vec(s.S) + msz::vec(s.C0)});
  out.push_back({"nested_beta", msz::fin(s.beta), msz::vec(s.beta)});
  out.push_back({"chain", static_cast<long long>(s.n_kept),
                 msz::vec(s.kept) + msz::vec2(s.kept_right) + msz::vec(s.mbar)});
  out.push_back({"pr_Qz", msz::fin(s.prQz), msz::vec(s.prQz)});
  out.push_back({"pr_Mbar", msz::fin(s.prMbar), msz::vec(s.prMbar)});
  out.push_back({"pr_QbIn", msz::fin(s.prQbIn), msz::vec(s.prQbIn)});
  out.push_back({"Bml", msz::fin(s.prBml), msz::vec(s.prBml)});
  out.push_back({"pr_nested_b", msz::fin(s.prQb_b), msz::vec(s.prQb_b)});
  out.push_back({"pr_kept", static_cast<long long>(s.pr_n_kept + s.pr_n_partner),
                 msz::vec(s.prKept) + msz::vec2(s.prKeptRight) + msz::vec2(s.prKeptRightP)});
  out.push_back({"pr_S", static_cast<long long>(s.prS.size()), msz::vec(s.prS)});
}


