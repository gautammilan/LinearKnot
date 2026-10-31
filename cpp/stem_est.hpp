// Persistent log-sum index of the --stem-est estimate (mirrors PartnerLogSums / PersistentSum in
// dp03/gate_index.py; bit-identical).
//
// Entries w(x, y) in log space (0-based, x < y); leaf(x, q) = log sum of w(x, y) over y > q,
// computed as the suffix log-sum of x's entries accumulated in decreasing y from NEG_INF.
// range_logsum(q, lo, hi) = log sum over x in [lo, hi] of leaf(x, q), O(log n): a log-sum segment
// tree over x per q, built by path copying as q advances (entries with y = q leave).
#pragma once

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace stemidx {

constexpr double NEG_INF = -1e300;

inline double lse(double a, double b) {
  if (a < b) std::swap(a, b);
  if (a < NEG_INF / 2) return NEG_INF;
  if (a - b > 50.0) return a;
  return a + std::log1p(std::exp(b - a));
}

struct Entry {
  int x, y;
  double lw;
};

class PartnerLogSums {
 public:
  PartnerLogSums() = default;
  PartnerLogSums(int n, const std::vector<Entry>& ents) : n_(n), ys_(n), suf_(n), lefts_(n) {
    std::vector<std::vector<std::pair<int, double>>> rows(n);
    for (const Entry& e : ents)
      if (0 <= e.x && e.x < e.y && e.y < n) {
        rows[e.x].emplace_back(e.y, e.lw);
        lefts_[e.y].push_back(e.x);
      }
    for (int x = 0; x < n; ++x) {
      auto& row = rows[x];
      std::stable_sort(row.begin(), row.end(),
                       [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                         return a.first < b.first;
                       });
      std::vector<double>& suf = suf_[x];
      suf.assign(row.size() + 1, NEG_INF);
      for (size_t k = row.size(); k-- > 0;) suf[k] = lse(suf[k + 1], row[k].second);
      ys_[x].reserve(row.size());
      for (const auto& t : row) ys_[x].push_back(t.first);
    }
    for (auto& row : lefts_) std::sort(row.begin(), row.end());
  }

  const std::vector<int>& lefts(int y) const { return lefts_[y]; }

  double value(int x, int q) const {
    if (x < 0 || x >= n_) return NEG_INF;
    const std::vector<int>& ys = ys_[x];
    size_t a = std::upper_bound(ys.begin(), ys.end(), q) - ys.begin();
    return suf_[x][a];
  }

  size_t bytes() const {
    size_t b = 0;
    for (int x = 0; x < n_; ++x) b += ys_[x].size() * 4 + suf_[x].size() * 8;
    return b;
  }

 private:
  int n_ = 0;
  std::vector<std::vector<int>> ys_;
  std::vector<std::vector<double>> suf_;
  std::vector<std::vector<int>> lefts_;
};

class PersistentSum {
 public:
  PersistentSum() = default;
  PersistentSum(int n, const std::vector<Entry>& ents) : n_(n), sums_(n, ents) {
    if (n == 0) return;
    roots_.reserve(n);
    roots_.push_back(build(0, n - 1));
    std::vector<int> changed;
    for (int q = 1; q < n; ++q) {
      const std::vector<int>& out = sums_.lefts(q);
      changed.assign(out.begin(), out.end());
      std::sort(changed.begin(), changed.end());
      changed.erase(std::unique(changed.begin(), changed.end()), changed.end());
      int root = roots_[q - 1];
      for (int x : changed) root = update(root, 0, n - 1, x, sums_.value(x, q));
      roots_.push_back(root);
    }
  }

  // log sum over x in [lo, hi] of leaf(x, q); NEG_INF for an empty range or q outside [0, n).
  double range_logsum(int q, int lo, int hi) const {
    lo = std::max(lo, 0);
    hi = std::min(hi, n_ - 1);
    if (lo > hi || q < 0 || q >= n_) return NEG_INF;
    return query(roots_[q], 0, n_ - 1, lo, hi);
  }

  size_t nodes() const { return s_.size(); }
  size_t bytes() const { return s_.size() * 16 + roots_.size() * 4 + sums_.bytes(); }

 private:
  int n_ = 0;
  PartnerLogSums sums_;
  std::vector<int> l_, r_;
  std::vector<double> s_;
  std::vector<int> roots_;

  int node(int left, int right, double s) {
    l_.push_back(left);
    r_.push_back(right);
    s_.push_back(s);
    return static_cast<int>(s_.size()) - 1;
  }

  int build(int lo, int hi) {
    if (lo == hi) return node(-1, -1, sums_.value(lo, 0));
    int mid = (lo + hi) / 2;
    int a = build(lo, mid);
    int b = build(mid + 1, hi);
    return node(a, b, lse(s_[a], s_[b]));
  }

  int update(int nd, int lo, int hi, int x, double val) {
    if (lo == hi) return node(-1, -1, val);
    int mid = (lo + hi) / 2;
    int a = l_[nd], b = r_[nd];
    if (x <= mid)
      a = update(a, lo, mid, x, val);
    else
      b = update(b, mid + 1, hi, x, val);
    return node(a, b, lse(s_[a], s_[b]));
  }

  double query(int nd, int lo, int hi, int a, int b) const {
    if (a <= lo && hi <= b) return s_[nd];
    int mid = (lo + hi) / 2;
    if (b <= mid) return query(l_[nd], lo, mid, a, b);
    if (a > mid) return query(r_[nd], mid + 1, hi, a, b);
    return lse(query(l_[nd], lo, mid, a, b), query(r_[nd], mid + 1, hi, a, b));
  }
};

// Point min-update, range-min query over [0, n) (integers; empty = INF).
class MinTree {
 public:
  static constexpr int INF = 1 << 30;
  explicit MinTree(int n = 0) {
    while (size_ < std::max(1, n)) size_ *= 2;
    t_.assign(2 * static_cast<size_t>(size_), INF);
  }
  void set_min(int k, int val) {
    k += size_;
    if (val >= t_[k]) return;
    t_[k] = val;
    for (k /= 2; k; k /= 2) t_[k] = std::min(t_[2 * k], t_[2 * k + 1]);
  }
  int query(int lo, int hi) const {
    int res = INF;
    lo += size_;
    hi += size_ + 1;
    while (lo < hi) {
      if (lo & 1) res = std::min(res, t_[lo++]);
      if (hi & 1) res = std::min(res, t_[--hi]);
      lo /= 2;
      hi /= 2;
    }
    return res;
  }

 private:
  int size_ = 1;
  std::vector<int> t_;
};

}  // namespace stemidx
