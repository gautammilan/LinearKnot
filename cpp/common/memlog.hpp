// Memory and CPU logging: MemItem / msz (bytes per structure), mem_rss_mb, mem_log
// (--mem-log), and the sampling profiler (--sample-prof).
//
// Fragment of the unity build: #included by cpp/main.cpp at this point of the original
// file (inside namespace / class DP03 / main() as main.cpp shows); not a standalone header.

// --mem-log on: stage checkpoints on stderr (RSS, live heap = malloc in use, CPU) with per
// structure entries and bytes.  Bytes: container capacities plus 16 per heap block and the
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
static double mem_rss_mb() {
  long pages = 0, res = 0;
  if (FILE* f = std::fopen("/proc/self/statm", "r")) {
    if (std::fscanf(f, "%ld %ld", &pages, &res) != 2) res = 0;
    std::fclose(f);
  }
  return static_cast<double>(res) * static_cast<double>(sysconf(_SC_PAGESIZE)) / 1048576.0;
}
static double mem_heap_mb() {
  const struct mallinfo2 mi = mallinfo2();
  return static_cast<double>(mi.uordblks + mi.hblkhd) / 1048576.0;
}
inline void stem_mem_items(const StemEst& s, std::vector<MemItem>& out) {
  const size_t se = s.entries + s.entries2 + s.entries_av;
  out.push_back({"stem_idx", static_cast<long long>(se),
                 12.0 * static_cast<double>(se) +
                     16.0 * static_cast<double>(s.index.nodes() + s.index2.nodes() + s.index_av.nodes())});
  out.push_back({"stem_S_C0", static_cast<long long>(s.S.size() + s.C0.size()), msz::vec(s.S) + msz::vec(s.C0)});
  out.push_back({"nested_beta", msz::fin(s.beta), msz::vec(s.beta)});
  out.push_back({"chain", static_cast<long long>(s.n_kept),
                 msz::vec(s.kept) + msz::vec2(s.kept_right) + msz::vec(s.mbar)});
  out.push_back({"hr_lphi", msz::fin(s.lphi), msz::vec(s.lphi)});
  out.push_back({"pr_Qz", msz::fin(s.prQz), msz::vec(s.prQz)});
  out.push_back({"pr_Mbar", msz::fin(s.prMbar), msz::vec(s.prMbar)});
  out.push_back({"pr_QbIn", msz::fin(s.prQbIn), msz::vec(s.prQbIn)});
  out.push_back({"Bml", msz::fin(s.prBml), msz::vec(s.prBml)});
  out.push_back({"pr_nested_b", msz::fin(s.prQb_b) + msz::fin(s.prQm_b) + msz::fin(s.prQm2_b),
                 msz::vec(s.prQb_b) + msz::vec(s.prQm_b) + msz::vec(s.prQm2_b)});
  out.push_back({"pr_kept", static_cast<long long>(s.pr_n_kept + s.pr_n_partner),
                 msz::vec(s.prKept) + msz::vec2(s.prKeptRight) + msz::vec2(s.prKeptRightP)});
  out.push_back({"pr_LP", msz::fin(s.prLP), msz::vec(s.prLP)});
  out.push_back({"pr_S", static_cast<long long>(s.prS.size()), msz::vec(s.prS)});
}
// One [mem] line: stage, CPU (total, since the previous checkpoint), RSS, peak RSS, live heap,
// then name:entries:MB:bytes_per_entry per item.
static bool g_mem_log = false;
static std::string g_mem_sid = "-";
static double g_mem_cpu = 0.0;
static void mem_log(const char* stage, const std::vector<MemItem>& items = {}) {
  if (!g_mem_log) return;
  const double c = static_cast<double>(std::clock()) / CLOCKS_PER_SEC;
  double hwm = 0.0;
  if (FILE* f = std::fopen("/proc/self/status", "r")) {
    char line[256];
    while (std::fgets(line, sizeof line, f))
      if (std::strncmp(line, "VmHWM:", 6) == 0) hwm = std::atof(line + 6) / 1024.0;
    std::fclose(f);
  }
  std::string s;
  char b[256];
  for (const MemItem& it : items) {
    if (it.ents == 0 && it.bytes == 0.0) continue;
    std::snprintf(b, sizeof b, " %s:%lld:%.1f:%.1f", it.name.c_str(), it.ents, it.bytes / 1048576.0,
                  it.ents > 0 ? it.bytes / static_cast<double>(it.ents) : 0.0);
    s += b;
  }
  std::fprintf(stderr, "[mem] %s stage %s cpu %.2f dcpu %.2f rss_mb %.0f hwm_mb %.0f heap_mb %.0f |%s\n",
               g_mem_sid.c_str(), stage, c, c - g_mem_cpu, mem_rss_mb(), hwm, mem_heap_mb(), s.c_str());
  std::fflush(stderr);
  g_mem_cpu = c;
}

// --sample-prof PATH: SIGPROF samples of the scan (dp.pfunc), 100 Hz. perf and
// gperftools are not installed here. Frames need a binary built with
// -fno-omit-frame-pointer. Off unless the flag is passed.
static int g_sample_fd = -1;
static void sample_on_sig(int) {
  void* stack[48];
  int n = backtrace(stack, 48);
  if (n < 1) n = 1;
  if (n > 48) n = 48;
  unsigned char nb = static_cast<unsigned char>(n);
  if (g_sample_fd >= 0) {
    ssize_t w = write(g_sample_fd, &nb, 1);
    if (w == 1) {
      const size_t bytes = static_cast<size_t>(n) * sizeof(void*);
      ssize_t left = static_cast<ssize_t>(bytes);
      const char* p = reinterpret_cast<const char*>(stack);
      while (left > 0) {
        ssize_t k = write(g_sample_fd, p, static_cast<size_t>(left));
        if (k <= 0) break;
        p += k;
        left -= k;
      }
    }
  }
}
static void sample_start(const char* path) {
  g_sample_fd = open(path, O_CREAT | O_TRUNC | O_WRONLY, 0644);
  if (g_sample_fd < 0) {
    std::perror(path);
    return;
  }
  const char magic[4] = {'S', 'A', 'M', 'P'};
  const uint64_t anchor = static_cast<uint64_t>(reinterpret_cast<uintptr_t>(&sample_on_sig));
  if (write(g_sample_fd, magic, 4) != 4 || write(g_sample_fd, &anchor, sizeof anchor) != static_cast<ssize_t>(sizeof anchor)) {
    std::perror(path);
    return;
  }
  struct sigaction sa;
  std::memset(&sa, 0, sizeof sa);
  sa.sa_handler = sample_on_sig;
  sigemptyset(&sa.sa_mask);
  sa.sa_flags = SA_RESTART;
  sigaction(SIGPROF, &sa, nullptr);
  struct itimerval it;
  std::memset(&it, 0, sizeof it);
  it.it_interval.tv_usec = 10000;
  it.it_value.tv_usec = 10000;
  setitimer(ITIMER_PROF, &it, nullptr);
  std::fprintf(stderr, "[sample] on %s 100 Hz\n", path);
}
static void sample_stop() {
  struct itimerval it;
  std::memset(&it, 0, sizeof it);
  setitimer(ITIMER_PROF, &it, nullptr);
  if (g_sample_fd >= 0) {
    close(g_sample_fd);
    g_sample_fd = -1;
  }
}

