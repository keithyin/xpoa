# Changelog

格式遵循 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/);
版本号遵循 [Semantic Versioning](https://semver.org/lang/zh-CN/)。

## [0.5.0] - 2026-11-21

C++ 核心算法完全以 Rust 重写,双后端共存、运行时可切换,结果逐字符一致,
性能约 2.1×。本次为 breaking release。

### Added

- **纯 Rust 后端**(`src/native/`,零 C++/Boost/SeqAn 依赖),模块与 C++ 文件
  一一对应,覆盖全部 FFI 可达路径:
  - `poa_graph.rs` — POA 图、稀疏 DP 对齐列、traceback 穿线、共识路径
  - `range_finder.rs` — `SdpRangeFinder` trait + V1/V2(新增 V3 无需改图代码)
  - `sparse_poa.rs` / `mod.rs` — read 编排与入口
  - `qgram.rs` / `find_seeds.rs` / `chain_seeds.rs` / `sparse_align.rs` /
    `seed.rs` — q-gram 索引与 SDP 种子链接(精确复刻 `std::set`/`multiset`
    等价拒绝语义)
  - `edlib.rs` — Myers 位向量编辑距离(NW + TASK_DISTANCE 子集)
  - `sequence.rs` — IUPAC 全互补表
- **运行时后端选择**:`PoaSetting` 新增 `backend` 字段
  (`BACKEND_CPP = 0` / `BACKEND_RUST = 1`),同一二进制内切换;
  `--no-default-features` 可裁掉 C++ 构建链
- **错误处理**:`src/error.rs` 定义 `XpoaError`(thiserror),公共 API 改返回
  `XpoaResult<Option<(String, usize)>>`;空 read、NUL 字节、非法 version 在
  进入任何后端前统一拦截(此前 C++ 会 segfault / 异常穿越 FFI)
- **差分测试体系**:
  - 161 个 golden 快照(`tests/golden/`,由 C++ 后端生成期望值,见
    `cargo run --release --example golden_gen`)
  - 属性式双后端实时对拍(`tests/diff_tests.rs`,错误路径同样差分)
  - 组件级单测:edlib 对暴力 DP、种子查找对穷举搜索等独立参照
- **基准**(`benches/poa_bench.rs`,criterion):3 种典型负载 × v1/v2 × 双后端
- **调试工具**:`examples/native_case.rs`、`examples/native_dump.rs`、
  `tools/dump_case.cpp`、`tools/dump_cols.cpp`(双端图/列/锚点差分定位)
- `third_party/README.md`:Boost 1.74 / SeqAn 2.4.0 头文件的免 root 获取步骤
- **GitHub Actions CI/CD**(`.github/workflows/`):
  - `ci.yml`:push/PR 触发 —— fmt + clippy(零警告强制)、纯 Rust 测试、
    双后端测试 + golden 快照跨环境复现校验
  - `release.yml`:`v*` tag 触发 —— tag/版本一致性校验、双后端 + 纯 Rust
    全量测试、crate 打包校验、发布 crates.io(需仓库 secret
    `CARGO_REGISTRY_TOKEN`)、创建 GitHub Release

### Changed

- **breaking**:公共 API `poa_consensus` / `poa_consensus_fwd` / `poa_consensus_v2`
  返回类型 `Option<(String, usize)>` → `XpoaResult<Option<(String, usize)>>`;
  `PoaSetting` 新增 `backend` 字段(手动构造字面量需补一行)
- `get_default_poa_setting()` 返回当前构建下可用的默认后端
- `build.rs` 仅在 `cpp` feature(默认启用)时构建 C++;`build.rs` 注入
  `third_party/` 头文件路径
- `cpp_code/Settings.h` 末尾新增 `backend` 字段(仅同步 ABI 布局,C++ 端忽略)

### Fixed

- 空 read:两后端统一返回 `Err(XpoaError::EmptyRead)`(C++ 原版会 segfault)
- 非法 version / 非法碱基:两后端统一返回 `Err`(C++ 原版的异常会穿越
  extern "C" 边界,属 UB)
- 现有 4 个 smoke 测试此前无断言,已补齐

### Known differences vs C++

- `RangeFinderV2` 在锚点位于共识路径末尾 KMER-1 个位置内时,C++ 存在越界读
  (UB);Rust 按预期语义跳过越界下标。C++ 触发 UB 的输入行为不可靠,不做对齐。

### Performance

- criterion 实测,全部场景(30×1.5kb / 100×500bp / 50×60bp reads,v1/v2)
  Rust 后端约 **2.1×** 快于 C++ 后端 —— 收益来自容器与内存布局
  (`Vec` 邻接表替代 BGL setS/listS 堆节点、id 直接索引替代 map 间接寻址、
  `BTreeMap`+`Vec` 替代 multiset 红黑树节点),而非算法差异。

## [0.4.0] - 2026-04-16

- Rust wrapper(bindgen FFI + CMake 构建 `libxpoa.a`)的初始公开形态
- `poa_consensus` / `poa_consensus_fwd` / `poa_consensus_v2` 公共 API
- RangeFinderV1/V2 分版本支持(`PoaSetting.version`)

## [0.3.2] - 2026-03-25

- 早期版本(Rust wrapper + C++ 核心)
