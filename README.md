# xpoa

DNA 测序的 partial-order alignment (POA) 共识序列生成库(源自 PacBio
pbcopper/unanimity 算法线)。

**v0.5.0:纯 Rust 实现 + C++ 实现双后端共存,运行时可切换,结果逐字符一致。**

```rust
use xpoa::{poa_consensus, xpoa::get_default_poa_setting};

let reads = vec!["AACGGATCGGA", "AACGGATCGGA", "AACGGATCGGA"];
let (consensus, n_passes) = poa_consensus(&reads, &get_default_poa_setting()).unwrap();
```

## 后端选择

```rust
let mut setting = xpoa::get_default_poa_setting();
setting.backend = xpoa::BACKEND_RUST; // 纯 Rust(推荐,零 C++ 依赖)
setting.backend = xpoa::BACKEND_CPP;  // 历史 C++ 实现(需 cmake + Boost/SeqAn)
setting.version = 2;                  // RangeFinderV2
```

- 默认构建(`--features cpp`)两个后端都编译,Rust 后端经完整差分测试保证与
  C++ 输出**逐字符一致**。
- `--no-default-features` 纯 Rust 构建:不需要 cmake/Boost/SeqAn,选 CPP 后端返回
  `Err(XpoaError::CxxBackendUnavailable)`。

## CLI:xpoa-cli(按分组生成 consensus)

仓库为 Cargo workspace:`xpoa/`(lib,上述 API)与 `xpoa-cli/`(命令行工具)。
clap/rayon/indicatif 仅是 CLI 依赖,lib 依赖树不受影响。

```bash
cargo install --path xpoa-cli                  # 本地安装(默认含 cpp 后端)
cargo install --path xpoa-cli --no-default-features   # 纯 Rust 安装,零 C++ 依赖
xpoa-cli input.fa -o consensus.fa
```

输入 FASTA 的每个 qname 以 `|group=GROUP_ID` 结尾,**相同 GROUP_ID 的序列必须连续
出现**(否则 panic)。每组经 POA 生成 consensus,以 GROUP_ID 为 qname 输出,顺序与
输入组顺序一致:

```
>read_1|group=sampleA      →      >sampleA
ACGT...                           <consensus>
>read_2|group=sampleA
ACGT...
>read_3|group=sampleB      →      >sampleB
TTGA...                           <consensus>
```

| 参数 | 说明 |
|---|---|
| `-o, --output <PATH>` | 输出 FASTA(缺省 stdout) |
| `-t, --threads <N>` | 并行 consensus 线程数(0 = 全核) |
| `--backend <cpp\|rust>` | POA 后端(缺省 = 编译期默认) |
| `--poa-version <1\|2>` | POA 算法版本(默认 1) |
| `-w, --line-width <N>` | 输出折行宽度,默认 60(0 = 单行) |

读/算/写为流水线处理(reader 线程解析分组 → rayon 并行计算 → 按序写出),
进度条输出到 stderr。

## API

三个入口均返回 `XpoaResult<Option<(String, usize)>>`:

| 函数 | 说明 |
|---|---|
| `poa_consensus(reads, setting)` | 基本共识(edlib 链方向统一,`ed_unify_strand`) |
| `poa_consensus_fwd(reads, setting)` | 假定 read 已正向,不做定向 |
| `poa_consensus_v2(reads, setting)` | 混合方向:Rust 侧按 `is_fwd()` 预处理后再走 FwdStrand |

- `Ok(Some((consensus, n_passes)))` — 正常结果
- `Ok(None)` — 空 read 集合
- `Err(XpoaError)` — 输入/配置错误(见下),两后端行为一致

`TSubread` trait(`get_seq` / `get_cx` / `is_fwd`)抽象输入;`&str` / `String` 已实现。

### 错误处理(`xpoa::error::XpoaError`)

| 变体 | 触发条件 |
|---|---|
| `EmptyRead { index }` | 空 read 序列(C++ 旧版对此会 segfault,现两后端统一拒绝) |
| `NulByte { index }` | read 含 NUL 字节,无法穿越 C FFI |
| `InvalidBase { base }` | 反向互补遇到无法互补的碱基(对应 C++ throw) |
| `InvalidVersion { version }` | `version` 非 1/2(对应 C++ throw) |
| `CxxBackendUnavailable` | 纯 Rust 构建下选了 cpp 后端 |

校验在进入任何后端之前完成。

## 测试与验证

```bash
cargo test --release          # 单测 + golden(161 case)+ 属性式双后端对拍
cargo bench                   # criterion:rust vs cpp 性能对比
```

## 为什么 Rust 版本更快

两个后端的**算法与浮点评分完全相同**(edlib 位向量、稀疏 DP 列、traceback、
tie-break 逐行对应;edlib 部分两版性能持平)。差距全部来自数据结构与内存布局 ——
POA 的热点(图构建、traceback/tagSpan、consensusPath 重算)恰好是 C++ 容器
节点分配开销最集中的地方。criterion 实测(30×1.5kb、100×500bp、50×60bp reads,
v1/v2 全场景)Rust 稳定约 **2.1×** 快于 C++ 后端:

| 来源 | C++(BGL + SeqAn 容器) | Rust(本实现) |
|---|---|---|
| 图存储 | `adjacency_list<setS, listS, bidirectionalS>`:边与顶点各是 `std::set` / `std::list` 的独立堆节点;每次 `add_edge` 一次节点分配 + O(log n) 有序插入;遍历依赖指针追逐,cache 不友好 | `Vec<PoaNode>` + `Vec<Vec<usize>>` 邻接表;`add_edge` 为摊销 O(1) 的 Vec push(以 `HashSet` 复刻 BGL setS 的去重语义);遍历为连续内存线性扫描 |
| 顶点寻址 | external id / internal descriptor / vertex_index 三套标识,经 `std::map` 与 property map 间接转换 | 三套 id 恒等(均为插入序),顶点 id 直接作 `Vec` 下标,零转换层 |
| 对齐列 | 每列 `new AlignmentColumn`,列集合为 `unordered_map<VD, col*>` 散列节点 | 列集合为 `Vec<Option<AlignmentColumn>>`,按顶点 id 连续存放 |
| 种子层 | `std::multiset<Seed>` / `std::set<SDPHit>` / `std::set<SDPColumn>` 红黑树节点,每个种子一次堆分配 | `BTreeMap<diagonal, Vec<Seed>>`:同对角线种子连续存放,树节点数降至对角线数级别 |
| 拓扑序 | BGL `topological_sort` 经 graph_traits 迭代器抽象做 DFS | 直接数组上的显式栈 DFS,序与 BGL 一致 |

除上表外,两侧热路径均无锁;唯一的多态调用(RangeFinder,每顶点一次)与
edlib/DP 的位运算也逐行对应 —— 因此 2.1× 是容器与内存布局的净收益,而非算法差异。

## C++ 构建期依赖

默认构建需要 cmake 与 Boost/SeqAn 头文件(获取步骤见 `xpoa/third_party/README.md`);
纯 Rust 构建完全不需要。详见 [CLAUDE.md](CLAUDE.md) 与
[xpoa/docs/CHANGELOG.md](xpoa/docs/CHANGELOG.md)。

## bindgen 再生成(改 C++ 接口头文件后)

在 `xpoa/` 目录下执行:

```
bindgen cpp_code/PoaDraft.h \
  -o src/xpoa_sys.rs \
  --allowlist-function "PoaDraft.*" \
  -- -x c++ -std=c++17 -Icpp_code -Icpp_code/edlib \
     -I/usr/include/c++/13 -I/usr/include/x86_64-linux-gnu/c++/13
```
