// Sparse DP storage: Ix (index), OMap (insertion-ordered hash map), AtView, SpMat, TriMat,
// Col, RunSt, and the 4-index tables Tab4 / Ent / Grp with key2, g4, add4. Insertion order
// is kept so log-space sums run in the same order as the Python engine.
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

// ---------------------------------------------------- insertion-ordered maps
// int key in [0, n) -> double; deleted keys leave a tombstone (key = -1) so
// re-insertion goes to the end, like a Python dict.
// Open-addressed key -> slot. Keys are >= 0. Empty is -1, tombstone is -2.
struct Ix {
  std::vector<int> k, v;
  uint32_t mask = 0;
  int live = 0, tombs = 0;

  static uint32_t hash(int x) { return static_cast<uint32_t>(x) * 0x9e3779b1u; }
  void clear() { k.clear(); v.clear(); mask = 0; live = tombs = 0; }
  double bytes() const {
    return (k.capacity() ? static_cast<double>(k.capacity() * sizeof(int) + 16) : 0.0) +
           (v.capacity() ? static_cast<double>(v.capacity() * sizeof(int) + 16) : 0.0);
  }
  void grow() {
    size_t cap = k.empty() ? 8 : k.size() * 2;
    std::vector<int> nk(cap, -1), nv(cap);
    uint32_t m = static_cast<uint32_t>(cap - 1);
    for (size_t i = 0; i < k.size(); ++i) {
      if (k[i] < 0) continue;
      uint32_t s = hash(k[i]) & m;
      while (nk[s] >= 0) s = (s + 1) & m;
      nk[s] = k[i];
      nv[s] = v[i];
    }
    k.swap(nk);
    v.swap(nv);
    mask = m;
    tombs = 0;
  }
  int slot_of(int key) const {
    if (k.empty() || key < 0) return -1;
    uint32_t s = hash(key) & mask;
    for (;;) {
      int c = k[s];
      if (c == key) return v[s];
      if (c == -1) return -1;
      s = (s + 1) & mask;
    }
  }
  void put(int key, int slot) {
    if (k.empty() || (live + tombs + 1) * 2 > static_cast<int>(k.size())) grow();
    uint32_t s = hash(key) & mask;
    int tomb = -1;
    for (;;) {
      int c = k[s];
      if (c == key) {
        v[s] = slot;
        return;
      }
      if (c == -2 && tomb < 0) tomb = static_cast<int>(s);
      if (c == -1) {
        if (tomb >= 0) {
          k[tomb] = key;
          v[tomb] = slot;
          --tombs;
        } else {
          k[s] = key;
          v[s] = slot;
        }
        ++live;
        return;
      }
      s = (s + 1) & mask;
    }
  }
  void del(int key) {
    if (k.empty() || key < 0) return;
    uint32_t s = hash(key) & mask;
    for (;;) {
      int c = k[s];
      if (c == key) {
        k[s] = -2;
        --live;
        ++tombs;
        return;
      }
      if (c == -1) return;
      s = (s + 1) & mask;
    }
  }
};

struct OMap {
  std::vector<int> keys;
  std::vector<double> vals;
  Ix pos;
  int live = 0;

  void init(int) {
    keys.clear();
    vals.clear();
    pos.clear();
    live = 0;
  }
  bool has(int k) const { return pos.slot_of(k) >= 0; }
  double get(int k) const {
    int p = pos.slot_of(k);
    return p >= 0 ? vals[p] : NEG_INF;
  }
  void set(int k, double v) {
    int p = pos.slot_of(k);
    if (p >= 0) {
      vals[p] = v;
      return;
    }
    pos.put(k, static_cast<int>(keys.size()));
    keys.push_back(k);
    vals.push_back(v);
    ++live;
  }
  void erase(int k) {
    int p = pos.slot_of(k);
    if (p < 0) return;
    keys[p] = -1;
    pos.del(k);
    --live;
  }
  // Drops tombstones, keeping live keys in insertion order.  Must not be
  // called while iterating over slots().
  void compact() {
    size_t w = 0;
    for (size_t r = 0; r < keys.size(); ++r) {
      int k = keys[r];
      if (k < 0) continue;
      keys[w] = k;
      vals[w] = vals[r];
      ++w;
    }
    keys.resize(w);
    vals.resize(w);
    pos.clear();
    for (size_t r = 0; r < w; ++r) pos.put(keys[r], static_cast<int>(r));
    if (keys.capacity() > 2 * w + 16) {
      keys.shrink_to_fit();
      vals.shrink_to_fit();
    }
  }
  bool sparse() const { return keys.size() > 2 * static_cast<size_t>(live) + 16; }
  size_t slots() const { return keys.size(); }
};

// Nested table stored as the kept keys of one OMap column. ix(i, j) = i * n + j.
// A read of a missing key is NEG_INF. Assigning NEG_INF drops the key.
struct AtView {
  std::vector<OMap>* cols = nullptr;
  int n = 0;
  void bind(std::vector<OMap>* c, int n_) { cols = c; n = n_; }
  bool empty() const { return n == 0 || cols == nullptr; }
  std::vector<double> dense() const {
    std::vector<double> o(static_cast<size_t>(n) * n, NEG_INF);
    if (!cols) return o;
    for (int j = 0; j < n; ++j)
      for (size_t r = 0; r < (*cols)[j].slots(); ++r) {
        int i = (*cols)[j].keys[r];
        if (i >= 0) o[static_cast<size_t>(i) * n + j] = (*cols)[j].vals[r];
      }
    return o;
  }
  struct Proxy {
    OMap* col;
    int i;
    operator double() const { return col->get(i); }
    Proxy& operator=(double v) {
      if (v < HALF_NEG) {
        if (col->has(i)) col->erase(i);
      } else {
        col->set(i, v);
      }
      return *this;
    }
    Proxy& operator=(const Proxy& o) { return *this = static_cast<double>(o); }
  };
  Proxy operator[](size_t k) { return Proxy{&(*cols)[k % static_cast<size_t>(n)], static_cast<int>(k / n)}; }
  double operator[](size_t k) const { return (*cols)[k % static_cast<size_t>(n)].get(static_cast<int>(k / n)); }
};

// Outside (and other) values that are finite only on kept keys. Same ix layout.
struct SpMat {
  int n = 0;
  std::vector<OMap> col;
  void init(int n_) {
    n = n_;
    col.assign(n_, {});
  }
  bool empty() const { return n == 0; }
  void clear() { col.clear(); n = 0; }
  long long finite() const {
    long long e = 0;
    for (const OMap& m : col) e += m.live;
    return e;
  }
  double bytes() const {
    double b = col.capacity() ? static_cast<double>(col.capacity() * sizeof(OMap) + 16) : 0.0;
    for (const OMap& m : col)
      b += (m.keys.capacity() ? static_cast<double>(m.keys.capacity() * sizeof(int) + 16) : 0.0) +
           (m.vals.capacity() ? static_cast<double>(m.vals.capacity() * sizeof(double) + 16) : 0.0) +
           m.pos.bytes();
    return b;
  }
  std::vector<double> dense() const {
    std::vector<double> o(static_cast<size_t>(n) * n, NEG_INF);
    for (int j = 0; j < n; ++j)
      for (size_t r = 0; r < col[j].slots(); ++r) {
        int i = col[j].keys[r];
        if (i >= 0) o[static_cast<size_t>(i) * n + j] = col[j].vals[r];
      }
    return o;
  }
  struct Proxy {
    OMap* col;
    int i;
    operator double() const { return col->get(i); }
    Proxy& operator=(double v) {
      if (v < HALF_NEG) {
        if (col->has(i)) col->erase(i);
      } else {
        col->set(i, v);
      }
      return *this;
    }
    Proxy& operator=(const Proxy& o) { return *this = static_cast<double>(o); }
  };
  Proxy operator[](size_t k) { return Proxy{&col[k % static_cast<size_t>(n)], static_cast<int>(k / n)}; }
  double operator[](size_t k) const { return col[k % static_cast<size_t>(n)].get(static_cast<int>(k / n)); }
};

// Qz is finite on every i <= j, so the column is a packed triangle rather than a hash.
struct TriMat {
  int n = 0;
  std::vector<std::vector<double>> col;
  void init(int n_) {
    n = n_;
    col.resize(n_);
    for (int j = 0; j < n_; ++j) col[j].assign(j + 1, NEG_INF);
  }
  bool empty() const { return n == 0; }
  long long finite() const {
    long long e = 0;
    for (const auto& c : col)
      for (double x : c) e += x > HALF_NEG;
    return e;
  }
  double bytes() const {
    double b = col.capacity() ? static_cast<double>(col.capacity() * sizeof(std::vector<double>) + 16) : 0.0;
    for (const auto& c : col)
      if (c.capacity()) b += static_cast<double>(c.capacity() * sizeof(double) + 16);
    return b;
  }
  struct Proxy {
    TriMat* m;
    int i, j;
    operator double() const {
      if (j < 0 || j >= m->n || i < 0 || i >= static_cast<int>(m->col[j].size())) return NEG_INF;
      return m->col[j][i];
    }
    Proxy& operator=(double v) {
      if (j < 0 || j >= m->n || i < 0) return *this;
      if (i >= static_cast<int>(m->col[j].size())) m->col[j].resize(i + 1, NEG_INF);
      m->col[j][i] = v;
      return *this;
    }
    Proxy& operator=(const Proxy& o) { return *this = static_cast<double>(o); }
  };
  Proxy operator[](size_t k) {
    return Proxy{this, static_cast<int>(k / n), static_cast<int>(k % static_cast<size_t>(n))};
  }
  double operator[](size_t k) const {
    int i = static_cast<int>(k / n), j = static_cast<int>(k % static_cast<size_t>(n));
    if (j < 0 || j >= n || i < 0 || i >= static_cast<int>(col[j].size())) return NEG_INF;
    return col[j][i];
  }
};

// Values over [0, n) with NEG_INF for absent keys; keys() in first-insertion order.
struct Col {
  std::vector<double> v;
  std::vector<char> in;
  std::vector<int> keys;
  void init(int n) {
    v.assign(n, NEG_INF);
    in.assign(n, 0);
    keys.clear();
  }
  void clear() {
    for (int k : keys) {
      v[k] = NEG_INF;
      in[k] = 0;
    }
    keys.clear();
  }
  double get(int k) const { return v[k]; }
  void put(int k, double x) {
    if (!in[k]) {
      in[k] = 1;
      keys.push_back(k);
    }
    v[k] = x;
  }
  void add(int k, double x) { put(k, log_add(in[k] ? v[k] : NEG_INF, x)); }
};

// Trailing-run states of one column (see DP03::run_ext); outside keeps U1, U2 only.
struct RunSt {
  Col P, U1, U2;
  void init(int n) {
    P.init(n);
    U1.init(n);
    U2.init(n);
  }
  void clear() {
    P.clear();
    U1.clear();
    U2.clear();
  }
};

struct Ent {
  int p, u;
  double a, b;  // inside, outside
};

struct Grp {
  int v;
  std::vector<Ent> ents;
  std::unordered_map<uint64_t, int> idx;
};

inline uint64_t key2(int p, int u) {
  return (static_cast<uint64_t>(static_cast<uint32_t>(p)) << 32) |
         static_cast<uint32_t>(u);
}

// groups[v][(p, u)] with insertion order at both levels.
struct Tab4 {
  std::vector<Grp> grps;
  std::unordered_map<int, int> vidx;

  Grp* get(int v) {
    auto it = vidx.find(v);
    return it == vidx.end() ? nullptr : &grps[it->second];
  }
  Ent* find(int v, int p, int u) {
    Grp* g = get(v);
    if (!g) return nullptr;
    auto it = g->idx.find(key2(p, u));
    return it == g->idx.end() ? nullptr : &g->ents[it->second];
  }
  void acc(int v, int p, int u, double val) {
    if (val < HALF_NEG) return;
    Grp* g;
    auto it = vidx.find(v);
    if (it == vidx.end()) {
      vidx.emplace(v, static_cast<int>(grps.size()));
      grps.push_back(Grp{v, {}, {}});
      g = &grps.back();
    } else {
      g = &grps[it->second];
    }
    uint64_t k = key2(p, u);
    auto jt = g->idx.find(k);
    if (jt == g->idx.end()) {
      g->idx.emplace(k, static_cast<int>(g->ents.size()));
      g->ents.push_back(Ent{p, u, log_add(NEG_INF, val), NEG_INF});
    } else {
      Ent& e = g->ents[jt->second];
      e.a = log_add(e.a, val);
    }
  }
};

inline double g4(Tab4& t, int v, int p, int u) {
  Ent* e = t.find(v, p, u);
  return e ? e->b : NEG_INF;
}

inline void add4(Ent& e, double val) {
  if (val < HALF_NEG) return;
  e.b = log_add(e.b, val);
}

