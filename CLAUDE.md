# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

xpoa 生成 DNA 测序的 partial-order alignment (POA) 共识序列,源自 PacBio pbcopper/unanimity。
**v0.5.0 起为双后端架构**:

- **Rust 后端**(`src/native/`,默认推荐):纯 Rust 重写,零 C++/Boost/SeqAn 依赖
- **C++ 后端**(`cpp_code/`):历史实现,经 bindgen FFI 链接 `libxpoa.a`

两后端**运行时可切换**(见下),结果保证逐字符一致(由差分测试体系保证)。

## Key Build Commands

```bash
cargo build --release          # 默认构建:Rust + C++ 双后端(需要 cmake + third_party 头文件)
cargo test --release           # 全部测试:单测 + 差分(golden + 属性式)
cargo test --release --no-default-features   # 纯 Rust 构建(无 cmake/Boost/SeqAn)
cargo bench                    # criterion 双后端性能对比
cargo clippy --release --all-targets && cargo fmt --check
```

### Backend Selection(运行时)

`PoaSetting` 末尾的 `backend` 字段:`xpoa::BACKEND_CPP`(0,历史默认)或 `xpoa::BACKEND_RUST`(1)。
`get_default_poa_setting()` 返回当前构建下可用的默认后端。
C++ 后端被裁掉(`--no-default-features`)时选 CPP 返回 `Err(XpoaError::CxxBackendUnavailable)`。

### Public API(返回 Result)

三个入口均返回 `XpoaResult<Option<(String, usize)>>`(`src/error.rs`):
`Ok(Some((consensus, n_passes)))` 为正常结果;`Ok(None)` 表示空 read 集合(历史语义);
`Err(XpoaError)` 为输入/配置错误 —— 空读、NUL 字节、非法 version、非法碱基(RC 路径)、
纯 Rust 构建下选 cpp 后端。校验在进入任何后端之前完成,两后端行为一致。

## Differential Testing(一致性保证体系)

迁移正确性由三层测试保证,任何算法改动后必须全部通过:

1. **Golden 快照**(`tests/golden/*.json`,161 case):`cargo run --release --example golden_gen`
   以 C++ 后端生成期望值,`tests/golden_tests.rs` 断言 Rust 后端重现。
2. **属性式对拍**(`tests/diff_tests.rs`):随机生成的用例在同一进程内实时对拍两后端。
3. **单元/穷举对拍**(`src/native/*` 各模块 `#[cfg(test)]`):edlib 对暴力 DP、种子查找对
   穷举搜索等独立参照。

### Debugging Tools(差分定位)

- `examples/native_case.rs` / `tools/dump_case.cpp`:参数化重放(读 TSV),输出每条 read
  后的图状态与中间共识,用于双端二分定位
- `tools/dump_cols.cpp` / `examples/native_dump.rs`:dump DP 列原始内容(score/move/prev)

## Architecture

### Rust Side

```
src/
├── lib.rs           # 公共 API(poa_consensus*,Result),按 setting.backend 分发
├── error.rs         # XpoaError 唯一定义点(thiserror)
├── setting.rs       # PoaSetting 唯一定义点(repr(C),含 backend 字段)
├── xpoa.rs          # TSubread trait、get_default_poa_setting、PoaResult(cpp)
├── xpoa_sys.rs      # bindgen FFI(仅 cpp feature)
├── dna.rs           # poa_consensus_v2 用的 ACGTN 互补表(既有行为,RC 返回 Result)
├── casegen.rs       # 差分测试用例生成器(含确定性 DetRng)
└── native/          # 纯 Rust 后端
    ├── mod.rs           # 入口层(≙ PoaDraft.cpp)
    ├── sparse_poa.rs    # read 编排(≙ SparsePoa.cpp)
    ├── poa_graph.rs     # 图 + 稀疏 DP 列 + traceback(≙ PoaGraphImpl/Traversals)
    ├── range_finder.rs  # SdpRangeFinder trait + V1/V2(≙ RangeFinder*)
    ├── sparse_align.rs  # 锚点入口(≙ SparseAlignment.h)
    ├── chain_seeds.rs   # SDP 种子链接(≙ ChainSeeds.h)
    ├── find_seeds.rs    # multiset<Seed> 语义(≙ FindSeeds.h/Seeds.h)
    ├── qgram.rs         # q-gram 索引(≙ Index/Hashing/IndexHits)
    ├── seed.rs          # Seed 类型(≙ Seed.h)
    ├── edlib.rs         # Myers 位向量 NW 距离(≙ edlib 子集)
    └── sequence.rs      # IUPAC 全互补表(≙ Sequence.cpp)
```

**模块边界原则**:每模块单一职责、与 C++ 文件一一对应(便于差分定位);
`poa_graph` 只依赖 `SdpRangeFinder` trait(V1/V2 经 `Box<dyn>` 注入);
后端分发用 `match`(两后端签名不同,不做 trait 抽象)。

### C++ Side(保留未改)

`cpp_code/` 原样保留(仅 `Settings.h` 末尾加了 `backend` 字段用于布局同步,C++ 端忽略其值)。
`build.rs` 在 `CARGO_FEATURE_CPP` 存在时经 CMake 构建 `libxpoa.a`。

## Fidelity Notes(移植必须知道的行为细节)

1. **FwdStrand 怪癖**:`PoaDraftGenWithAllFwdStrand` 用默认构造的 `SparsePoa`
   (**固定 V1 RangeFinder**,`setting.version` 被忽略)。
2. **V2 RangeFinder UB**:C++ `RangeFinderV2.cpp` 的 "add more bound" 循环在
   `cssPos+i` 超出共识路径时越界读(UB)。Rust 按**预期语义**跳过越界下标 ——
   这是算法行为上与 C++ 的唯一有意差异(C++ 触发 UB 的输入不可复现,属不可靠
   行为)。此外错误处理层面的统一拦截见第 7 条。
3. **BGL setS 去重**:`add_edge` 等价边自动去重 → `HashSet<(u,v)>` 查重。
4. **顶点 id 恒等**:external == internal == vertex_index == 插入序(PruneGraph 不可达)。
5. **std::set 等价拒绝**:种子层 sweepSet/colSet/Seeds(multiset)的排序与等价拒绝
   语义在 `chain_seeds.rs`/`find_seeds.rs` 有逐条对应实现。
6. **RC 表差异**:native 内部用 C++ 的 IUPAC 全表(`native/sequence.rs`);
   `poa_consensus_v2` 的 Rust 侧预处理保持 `dna.rs` 的 ACGTN 表(既有行为)。
7. **错误处理(两后端统一)**:公共 API 返回 `Result<Option<(String, usize)>, XpoaError>`
   (`src/error.rs`,thiserror)。空 read(`EmptyRead`)、NUL 字节(`NulByte`)、非法
   version(`InvalidVersion`)在进入任何后端前统一拦截 —— C++ 对这些输入会
   segfault / 异常穿越 FFI(UB);非法碱基(`InvalidBase`)在 RC 路径返回 Err
   (对应 C++ throw)。`Ok(None)` 保留空 read 集合的历史语义。
8. `consensusPath` 的 `-0.0001f` tie-break、`makeAlignmentColumnForExit` 的
   插入序 first-wins、`inEdges` 按 (src,dst) 排序 —— 均已精确复刻,改动前先看差分。

## Default Parameters

```rust
PoaSetting {
    min_identity: 0.82, match_score: 3, mismatch_score: -5,
    insertion_score: -2, deletion_score: -2,
    ed_unify_strand: 1, version: 1,   // 2 = RangeFinderV2
    backend: /* 构建决定 */,
}
```

## Third-party(C++ 构建期依赖,不进 git)

`third_party/README.md` 有完整获取步骤(boost 1.74 头文件、seqan 2.4.0 头文件);
`build.rs` 自动把存在的目录注入 `-I`。纯 Rust 构建完全不需要它们。

## Benchmark

`cargo bench`(criterion):Rust 后端在典型负载(30×1.5kb、100×500bp、50×60bp,
v1/v2)下稳定 **≈2.1× 快于 C++ 后端**(消除 BGL setS 节点分配与指针比较)。
