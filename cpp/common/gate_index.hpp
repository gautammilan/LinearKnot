// Partner-beyond-q index for the eligibility gate (mirrors dp03/gate_index.py).
//
// For pseudoknot-free pair probabilities P(x, y) (0-based, x < y) and a span M (M <= 0 means n):
//
//     r(x, q) = sum of P(x, y) over partners y with q < y <= q + M
//
// PersistentGate answers max over x in [lo, hi] of r(x, q) for any q in O(log n): a max segment
// tree over x for every q (a version), built by path copying as q advances (partners y = q leave,
// y = q + M enter).  Leaf values are always recomputed as pre[x][b] - pre[x][a], with the prefix
// sums of x's partners accumulated in increasing y from 0.0, so the values are bit-identical to
// the Python index.
#pragma once

#include <algorithm>
#include <utility>
#include <vector>

namespace gate {

struct Pair {
  int i, j;
  double p;
};

class PartnerSums {
 public:
  PartnerSums() = default;
  PartnerSums(int n, const std::vector<Pair>& pairs, int span)
      : n_(n), span_(span), ys_(n), pre_(n), lefts_(n) {
    std::vector<std::vector<std::pair<int, double>>> rows(n);
    for (const Pair& e : pairs)
      if (0 <= e.i && e.i < e.j && e.j < n) {
        rows[e.i].emplace_back(e.j, e.p);
        lefts_[e.j].push_back(e.i);
      }
    for (int x = 0; x < n; ++x) {
      auto& row = rows[x];
      std::stable_sort(row.begin(), row.end(),
                       [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                         return a.first < b.first;
                       });
      double acc = 0.0;
      pre_[x].reserve(row.size() + 1);
      pre_[x].push_back(0.0);
      ys_[x].reserve(row.size());
      for (const auto& t : row) {
        acc += t.second;
        pre_[x].push_back(acc);
        ys_[x].push_back(t.first);
      }
    }
    for (auto& row : lefts_) std::sort(row.begin(), row.end());
  }

  int n() const { return n_; }
  int span() const { return span_; }
  const std::vector<int>& lefts(int y) const { return lefts_[y]; }

  double r(int x, int q) const {
    if (x < 0 || x >= n_) return 0.0;
    const std::vector<int>& ys = ys_[x];
    const std::vector<double>& pre = pre_[x];
    size_t a = std::upper_bound(ys.begin(), ys.end(), q) - ys.begin();
    size_t b = span_ <= 0 ? ys.size()
                          : static_cast<size_t>(std::upper_bound(ys.begin(), ys.end(), q + span_) -
                                                ys.begin());
    return pre[b] - pre[a];
  }

 private:
  int n_ = 0;
  int span_ = 0;
  std::vector<std::vector<int>> ys_;
  std::vector<std::vector<double>> pre_;
  std::vector<std::vector<int>> lefts_;
};

class PersistentGate {
 public:
  PersistentGate() = default;
  PersistentGate(int n, const std::vector<Pair>& pairs, int span) : n_(n), sums_(n, pairs, span) {
    if (n == 0) return;
    roots_.reserve(n);
    roots_.push_back(build(0, n - 1));
    std::vector<int> changed;
    for (int q = 1; q < n; ++q) {
      const std::vector<int>& out = sums_.lefts(q);
      changed.assign(out.begin(), out.end());
      if (span > 0 && q + span <= n - 1) {
        const std::vector<int>& in = sums_.lefts(q + span);
        changed.insert(changed.end(), in.begin(), in.end());
      }
      std::sort(changed.begin(), changed.end());
      changed.erase(std::unique(changed.begin(), changed.end()), changed.end());
      int root = roots_[q - 1];
      for (int x : changed) root = update(root, 0, n - 1, x, sums_.r(x, q));
      roots_.push_back(root);
    }
  }

  // max over x in [lo, hi] of r(x, q); 0.0 for an empty range or q outside [0, n).
  double range_max(int q, int lo, int hi) const {
    lo = std::max(lo, 0);
    hi = std::min(hi, n_ - 1);
    if (lo > hi || q < 0 || q >= n_) return 0.0;
    return query(roots_[q], 0, n_ - 1, lo, hi);
  }

  double value(int x, int q) const { return sums_.r(x, q); }
  size_t nodes() const { return mx_.size(); }

 private:
  int n_ = 0;
  PartnerSums sums_;
  std::vector<int> l_, r_;
  std::vector<double> mx_;
  std::vector<int> roots_;

  int node(int left, int right, double mx) {
    l_.push_back(left);
    r_.push_back(right);
    mx_.push_back(mx);
    return static_cast<int>(mx_.size()) - 1;
  }

  int build(int lo, int hi) {
    if (lo == hi) return node(-1, -1, sums_.r(lo, 0));
    int mid = (lo + hi) / 2;
    int a = build(lo, mid);
    int b = build(mid + 1, hi);
    return node(a, b, std::max(mx_[a], mx_[b]));
  }

  int update(int nd, int lo, int hi, int x, double val) {
    if (lo == hi) return node(-1, -1, val);
    int mid = (lo + hi) / 2;
    int a = l_[nd], b = r_[nd];
    if (x <= mid)
      a = update(a, lo, mid, x, val);
    else
      b = update(b, mid + 1, hi, x, val);
    return node(a, b, std::max(mx_[a], mx_[b]));
  }

  double query(int nd, int lo, int hi, int a, int b) const {
    if (a <= lo && hi <= b) return mx_[nd];
    int mid = (lo + hi) / 2;
    if (b <= mid) return query(l_[nd], lo, mid, a, b);
    if (a > mid) return query(r_[nd], mid + 1, hi, a, b);
    return std::max(query(l_[nd], lo, mid, a, b), query(r_[nd], mid + 1, hi, a, b));
  }
};

// (i, j) -> (n-1-j, n-1-i): partners left of p become partners right of n-1-p.
inline std::vector<Pair> mirror_pairs(int n, const std::vector<Pair>& pairs) {
  std::vector<Pair> out;
  out.reserve(pairs.size());
  for (const Pair& e : pairs) out.push_back({n - 1 - e.j, n - 1 - e.i, e.p});
  std::sort(out.begin(), out.end(), [](const Pair& a, const Pair& b) {
    if (a.i != b.i) return a.i < b.i;
    if (a.j != b.j) return a.j < b.j;
    return a.p < b.p;
  });
  return out;
}

}  // namespace gate
