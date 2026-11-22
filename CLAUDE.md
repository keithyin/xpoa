# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

xpoa 生成 DNA 测序的 partial-order alignment (POA) 共识序列,算法源自 PacBio
pbcopper/unanimity。**v0.6.0 起为纯 Rust 实现**(`xpoa/src/native/`,零
C++/Boost/SeqAn 依赖);历史 C++ 后端(`cpp_code/`,经 FFI 调用)已于 0.6.0
移除 —— 其输出经 161 golden + 属性式差分验证与 Rust 实现逐字符一致,见
`xpoa/docs/CHANGELOG.md`。

**仓库布局(Cargo workspace)**:根为 virtual workspace;`xpoa/` 为 lib package,
`xpoa-cli/` 为 CLI package(按分组生成 consensus 的命令行工具,依赖
clap/rayon/indicatif,与 lib 依赖完全隔离)。virtual workspace 下带 feature
参数的 cargo 命令需显式 `-p` 指定包。lib 与 cli **版本锁步**:统一由根
`[workspace.package]` 继承;发版时同步更新 `[workspace.package]` 的 `version`
与 `[workspace.dependencies]` 中 xpoa 的 `version` 需求(跨小版本如 0.7.0 时
两者都要改,0.6.x 内升级 xpoa 无需动依赖需求)。

## Key Build Commands

```bash
cargo build --release          # workspace 全量:lib + xpoa-cli(纯 Rust,无外部依赖)
cargo test --release           # 全部测试:单测 + golden 回归基线(161 case)
cargo bench                    # criterion 性能基准(3 负载 × v1/v2)
cargo clippy --release --workspace --all-targets && cargo fmt --check
cargo run -p xpoa-cli -- input.fa -o out.fa   # CLI:按分组生成 consensus
```

## Public API(返回 Result)

三个入口均返回 `XpoaResult<Option<(String, usize)>>`(`src/error.rs`):
`Ok(Some((consensus, n_passes)))` 为正常结果;`Ok(None)` 表示空 read 集合(历史语义);
`Err(XpoaError)` 为输入/配置错误 —— 空读、NUL 字节、非法 version、
非法碱基(RC 路径)。校验在进入计算之前完成。

## Testing(回归保证体系)

1. **Golden 快照**(`xpoa/tests/golden/*.json`,161 case):回归基线。v0.5.0 时由
   C++ 后端生成、Rust 实现经差分验证一致;`cargo run --release -p xpoa --example golden_gen`
   用当前实现重生成,应与基线逐字节相同(CI 有此复现校验)。任何算法改动后
   `xpoa/tests/golden_tests.rs` 必须通过;**有意改变输出**时需重新生成基线并
   在 CHANGELOG 说明理由。
2. **单元/穷举对拍**(`src/native/*` 各模块 `#[cfg(test)]`):edlib 对暴力 DP、
   种子查找对穷举搜索等独立参照。

### Debugging Tools(流程定位)

- `examples/native_case.rs`(TSV 参数化重放)/ `examples/native_dump.rs`(case_004
  图状态 + DP 列 dump):输出每条 read 后的图状态与中间共识,供人工/脚本定位
  (历史上曾与已移除的 C++ 侧 dump 工具对拍)

## Architecture

```
xpoa/src/
├── lib.rs           # 公共 API(poa_consensus*,入口统一校验 validate)
├── error.rs         # XpoaError 唯一定义点(thiserror)
├── setting.rs       # PoaSetting 唯一定义点
├── xpoa.rs          # TSubread trait、get_default_poa_setting
├── dna.rs           # poa_consensus_v2 用的 ACGTN 互补表(既有行为,RC 返回 Result)
├── casegen.rs       # golden 用例生成器(含确定性 DetRng)
└── native/          # 纯 Rust 实现
    ├── mod.rs           # 入口层(≙ 历史 PoaDraft.cpp)
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

**模块边界原则**:每模块单一职责、与历史 C++ 文件一一对应(便于对照溯源);
`poa_graph` 只依赖 `SdpRangeFinder` trait(V1/V2 经 `Box<dyn>` 注入)。

## Behavior Notes(相对历史 C++ 基线的行为细节,改动前先看 golden)

1. **FwdStrand 怪癖**:FwdStrand 路径(`poa_consensus_fwd`/`poa_consensus_v2` 的
   底层)固定 V1 RangeFinder,`setting.version` 被忽略 —— 精确复刻历史 C++ 行为。
2. **V2 RangeFinder**:历史 C++ `RangeFinderV2.cpp` 的 "add more bound" 循环有
   越界读(UB);Rust 按**预期语义**跳过越界下标 —— 这是唯一有意算法差异
   (C++ 触发 UB 的输入不可复现,属不可靠行为)。
3. **BGL setS 去重**:历史 `add_edge` 等价边自动去重 → `HashSet<(u,v)>` 查重。
4. **顶点 id 恒等**:external == internal == vertex_index == 插入序
   (历史 PruneGraph 不可达,故可简化)。
5. **std::set 等价拒绝**:种子层 sweepSet/colSet/Seeds(multiset)的排序与等价拒绝
   语义在 `chain_seeds.rs`/`find_seeds.rs` 有逐条对应实现。
6. **RC 表差异**:native 内部用历史 C++ 的 IUPAC 全表(`native/sequence.rs`);
   `poa_consensus_v2` 的预处理保持 `dna.rs` 的 ACGTN 表(既有行为)。
7. **错误处理**:公共 API 返回 `Result<Option<(String, usize)>, XpoaError>`
   (`src/error.rs`,thiserror)。空 read(`EmptyRead`)、NUL 字节(`NulByte`)、
   非法 version(`InvalidVersion`)在入口统一拦截(历史 C++ 对这些输入会
   segfault / 异常穿越 FFI);非法碱基(`InvalidBase`)在 RC 路径返回 Err
   (对应历史 C++ throw)。`Ok(None)` 保留空 read 集合的历史语义。
8. `consensusPath` 的 `-0.0001f` tie-break、`makeAlignmentColumnForExit` 的
   插入序 first-wins、`inEdges` 按 (src,dst) 排序 —— 均已精确复刻,改动前先看 golden。

## Default Parameters

```rust
PoaSetting {
    min_identity: 0.82, match_score: 3, mismatch_score: -5,
    insertion_score: -2, deletion_score: -2,
    ed_unify_strand: 1, version: 1,   // 2 = RangeFinderV2
}
```

## Benchmark

`cargo bench`(criterion):3 种典型负载(30×1.5kb、100×500bp、50×60bp,
v1/v2)。历史参考:v0.5.0 对比实测本实现较 C++ 后端稳定快约 2.1×
(收益来自容器与内存布局,见 CHANGELOG 0.5.0)。
