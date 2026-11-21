# third_party

C++ 后端的系统依赖(不进 git,构建期专用)。

**注意:纯 Rust 构建(`cargo build --no-default-features`)完全不需要本目录** ——
`build.rs` 仅在 `cpp` feature 启用(默认)时才编译 C++。

- `boost/include/boost/` — Boost 1.74 头文件(仅头文件部分)。
  ```bash
  cd /tmp && mkdir boost_dl && cd boost_dl
  apt-get download libboost1.74-dev
  dpkg -x libboost1.74-dev_*.deb extracted
  mkdir -p <repo>/third_party/boost/include
  cp -r extracted/usr/include/boost <repo>/third_party/boost/include/boost
  ```

- `seqan/include/seqan/` — SeqAn 2.4.0 头文件(header-only)。
  ```bash
  cd /tmp
  curl -sL -o seqan240.tar.gz https://github.com/seqan/seqan/archive/refs/tags/seqan-v2.4.0.tar.gz
  tar xzf seqan240.tar.gz
  mkdir -p <repo>/third_party/seqan/include
  cp -r seqan-seqan-v2.4.0/include/seqan <repo>/third_party/seqan/include/seqan
  ```

`build.rs` 会把这两个目录作为 `-I` 传给 CMake;目录缺失时 build 仍会走系统 include 路径
(适用于已通过 apt 安装 `libboost-dev` 的环境)。
