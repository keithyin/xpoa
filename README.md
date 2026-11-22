# xpoa

DNA 测序的 partial-order alignment (POA) 共识序列生成库。算法源自 PacBio
pbcopper/unanimity,自 v0.5.0 起为纯 Rust 实现(零 C++/Boost/SeqAn 依赖),
算法与浮点评分与历史 C++ 实现逐字符一致(161 个 golden 快照基线)。

```rust
use xpoa::{poa_consensus, xpoa::get_default_poa_setting};

let reads = vec!["AACGGATCGGA", "AACGGATCGGA", "AACGGATCGGA"];
let (consensus, n_passes) = poa_consensus(&reads, &get_default_poa_setting()).unwrap();
```

## POA 版本

```rust
let mut setting = xpoa::get_default_poa_setting();
setting.version = 2; // RangeFinderV2(默认 1 = RangeFinderV1)
```

## CLI:xpoa-cli(按分组生成 consensus)

仓库为 Cargo workspace:`xpoa/`(lib,上述 API)与 `xpoa-cli/`(命令行工具)。
clap/rayon/indicatif 仅是 CLI 依赖,lib 依赖树不受影响。

```bash
cargo install --path xpoa-cli
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
| `poa_consensus_v2(reads, setting)` | 混合方向:按 `is_fwd()` 预处理(过滤 + RC)后走 FwdStrand |

- `Ok(Some((consensus, n_passes)))` — 正常结果
- `Ok(None)` — 空 read 集合
- `Err(XpoaError)` — 输入/配置错误(见下)

`TSubread` trait(`get_seq` / `get_cx` / `is_fwd`)抽象输入;`&str` / `String` 已实现。

### 错误处理(`xpoa::error::XpoaError`)

| 变体 | 触发条件 |
|---|---|
| `EmptyRead { index }` | 空 read 序列(旧版 C++ 对此会 segfault) |
| `NulByte { index }` | read 含 NUL 字节 |
| `InvalidBase { base }` | 反向互补遇到无法互补的碱基(旧版 C++ 对应 throw) |
| `InvalidVersion { version }` | `version` 非 1/2(旧版 C++ 对应 throw) |

校验在进入计算之前完成。

## 测试与验证

```bash
cargo test --release          # 单测 + golden 回归基线(161 case)
cargo bench                   # criterion 性能基准(3 负载 × v1/v2)
```

golden 快照(`xpoa/tests/golden/`)为回归基线:最初在 v0.5.0 由 C++ 实现生成,
纯 Rust 重写经双后端差分(161 golden + 随机属性对拍)验证逐字符一致后,
于 v0.6.0 移除 C++ 实现。`cargo run --release -p xpoa --example golden_gen`
可重生成快照做跨环境复现校验(CI 集成)。历史细节见
[xpoa/docs/CHANGELOG.md](xpoa/docs/CHANGELOG.md)。
