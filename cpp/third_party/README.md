# Vendored third-party headers

Nothing here is modified; the files are byte for byte what upstream publishes, so they can be
refreshed by re-downloading the URL below.

## ankerl_unordered_dense.h

- Upstream: https://github.com/martinus/unordered_dense (single header distribution)
- URL: https://raw.githubusercontent.com/martinus/unordered_dense/v4.5.0/include/ankerl/unordered_dense.h
- Version: 4.5.0 (`ANKERL_UNORDERED_DENSE_VERSION_*` at the top of the file); 4.5.0 is the last
  tag whose header is one self-contained file -- `main` splits it into `stl.h` + `detail/`.
- License: MIT (the full text is in the file's header)
- sha256: b5f67c6895a58e059273cf24bae16cb07e6c72f0e3f88bbebecaddd987b102b3
- Why: the index maps on the hot paths (`Grp::idx`, `Tab4::vidx`, `SMap::idx`, `lArm`) are
  chained `std::unordered_map`s over integer keys. This is the same interface in an open
  addressing table with the hash stored next to the key, which is the cheapest drop-in that
  does not change any floating-point result (see `cpp/common/fast_map.hpp`).

Used through `cpp/common/fast_map.hpp`, which is the only file that includes it.
