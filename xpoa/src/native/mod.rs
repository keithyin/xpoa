//! 纯 Rust 后端入口(对应 cpp_code/PoaDraft.cpp)。
//!
//! 职责:读取编排(minCov 计算、coverage 上限、方向统一、n_passes 统计);
//! 下层算法由 sparse_poa / poa_graph / range_finder / 种子层 / edlib 提供。
//!
//! 两个入口与 C++ 逐一对应,均返回 `XpoaResult`(空 read / 非法碱基 / 非法
//! version 返回 `Err`,见 `src/error.rs`):
//! - `poa_draft_gen` ≙ `PoaDraftGen`(edlib 链方向统一)
//! - `poa_draft_gen_fwd` ≙ `PoaDraftGenWithAllFwdStrand`
//!   **怪癖**:C++ 该路径用默认构造的 `SparsePoa`(固定 V1 RangeFinder),
//!   `setting.version` 被忽略 —— 精确复制。

pub mod chain_seeds;
pub mod edlib;
pub mod find_seeds;
pub mod qgram;
pub mod seed;
pub mod sequence;
pub mod sparse_align;
pub mod sparse_poa;

pub mod poa_graph;
pub mod range_finder;

use crate::error::{XpoaError, XpoaResult};
use crate::native::edlib::edlib_global_align_distance;
use crate::native::sequence::reverse_complement;
use crate::native::sparse_poa::{PoaAlignmentSummary, ReadKey, SparsePoa};
use crate::setting::PoaSetting;

/// LocalContextFlags 位(对应 PoaDraft.cpp 的枚举)
pub const ADAPTER_BEFORE: i32 = 1;
pub const ADAPTER_AFTER: i32 = 2;
/// 对应 `maxPoaCov`
pub const MAX_POA_COV: usize = 1_000_000;

/// native 后端的读取输入(seq 借用 + LocalContextFlags 位组合)
pub struct Read<'a> {
    pub seq: &'a str,
    pub flags: i32,
}

/// 对应 `_ComputeNPassesForNoPolish`
fn compute_n_passes(
    reads: &[Read<'_>],
    read_keys: &[ReadKey],
    summaries: &[PoaAlignmentSummary],
    min_identity: f32,
) -> usize {
    let mut n_passes = 0usize;
    for (i, read) in reads.iter().enumerate() {
        let key = read_keys[i];
        if key < 0 {
            continue;
        }
        if summaries[key as usize].alignment_identity < min_identity {
            continue;
        }
        let flags = read.flags;
        if !(flags & ADAPTER_BEFORE != 0 && flags & ADAPTER_AFTER != 0) {
            continue;
        }
        n_passes += 1;
    }
    n_passes
}

/// `_PoaDraftGenCore` 与 `_PoaDraftGenCoreWithAllFwdStrandCore` 的
/// read 定向策略(与 C++ 三个分支一一对应)
#[derive(Clone, Copy, PartialEq)]
enum DraftMode {
    /// ed_unify_strand != 0:read0 直接加,其余按 edlib 距离选正反向
    EdUnify,
    /// ed_unify_strand == 0:OrientAndAddRead(POA 分数定向)
    Orient,
    /// FwdStrand 入口:全部直接 AddRead,不做定向
    FwdDirect,
}

/// 对应 `_PoaDraftGenCore` / `_PoaDraftGenCoreWithAllFwdStrandCore` 的公共主体。
fn draft_gen_core(
    reads: &[Read<'_>],
    setting: &PoaSetting,
    mode: DraftMode,
    force_v1_range_finder: bool,
) -> XpoaResult<(Vec<u8>, usize, Vec<ReadKey>)> {
    for (i, read) in reads.iter().enumerate() {
        if read.seq.is_empty() {
            // C++ 端对空 read 会 segfault(threadFirstRead 的 UB);
            // 两后端统一在入口拒绝
            return Err(XpoaError::EmptyRead { index: i });
        }
    }

    let mut poa = if force_v1_range_finder {
        SparsePoa::new_v1()
    } else {
        SparsePoa::new(setting.version)?
    };
    let mut cov = 0usize;
    let mut read_keys = Vec::with_capacity(reads.len());

    let reads_0 = reads[0].seq.as_bytes();

    for (i, read) in reads.iter().enumerate() {
        let seq = read.seq.as_bytes();
        let key = match mode {
            DraftMode::EdUnify => {
                if i == 0 {
                    poa.add_read(seq, setting, true)
                } else {
                    let rc_seq = reverse_complement(seq)?;
                    let distance1 = edlib_global_align_distance(seq, reads_0);
                    let distance2 = edlib_global_align_distance(&rc_seq, reads_0);
                    // 严格大于:相等时保留正向(C++ 原样)
                    if distance1 > distance2 {
                        poa.add_read(&rc_seq, setting, true)
                    } else {
                        poa.add_read(seq, setting, false)
                    }
                }
            }
            DraftMode::Orient => poa.orient_and_add_read(seq, setting)?,
            DraftMode::FwdDirect => poa.add_read(seq, setting, false),
        };

        read_keys.push(key);
        if key >= 0 {
            cov += 1;
            if cov >= MAX_POA_COV {
                break;
            }
        }
    }

    // 至少一半 read 覆盖(C++ 的 minCov 公式)
    let min_cov: i32 = if cov < 5 {
        1
    } else {
        (cov.div_ceil(2) - 1) as i32
    };
    let (seq, summaries) = poa.find_consensus(min_cov, setting);
    let n_passes = compute_n_passes(reads, &read_keys, &summaries, setting.min_identity);
    Ok((seq, n_passes, read_keys))
}

/// 对应 C++ `PoaDraftGen`
pub fn poa_draft_gen(reads: &[Read<'_>], setting: &PoaSetting) -> XpoaResult<(String, usize)> {
    let mode = if setting.ed_unify_strand != 0 {
        DraftMode::EdUnify
    } else {
        DraftMode::Orient
    };
    let (seq, n_passes, _keys) = draft_gen_core(reads, setting, mode, false)?;
    // 输入经 `sequence.rs` 的 IUPAC 校验,共识字符必为 UTF-8 安全的 IUPAC 子集
    Ok((
        String::from_utf8(seq).expect("consensus must be UTF-8"),
        n_passes,
    ))
}

/// 对应 C++ `PoaDraftGenWithAllFwdStrand`
pub fn poa_draft_gen_fwd(reads: &[Read<'_>], setting: &PoaSetting) -> XpoaResult<(String, usize)> {
    // 怪癖:C++ 该入口固定 V1 RangeFinder(默认构造 SparsePoa),version 被忽略
    let (seq, n_passes, _keys) = draft_gen_core(reads, setting, DraftMode::FwdDirect, true)?;
    Ok((
        String::from_utf8(seq).expect("consensus must be UTF-8"),
        n_passes,
    ))
}
