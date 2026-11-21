//! SDP 对齐入口(对应 cpp_code/SparseAlignment.h 的移植子集)。
//!
//! `SparseAlign(q, seq1, seq2)` → 种子查找(q-gram + 同聚物过滤)→ 链接 →
//! 最佳链的 (beginH, beginV) 锚点对。

use crate::native::chain_seeds::{ChainSeedsConfig, chain_seeds};
use crate::native::find_seeds::find_seeds;
use crate::native::seed::Seed;

/// 对应 `CCS::SparseAlignSeeds`:返回最佳链(可能为空)
pub fn sparse_align_seeds(q_gram_size: usize, seq1: &[u8], seq2: &[u8]) -> Vec<Seed> {
    let config = ChainSeedsConfig::sparse_align();
    let seeds = find_seeds(q_gram_size, seq1, seq2, true);
    let chains = chain_seeds(&seeds, &config);
    if chains.is_empty() {
        Vec::new()
    } else {
        chains.into_iter().next().unwrap()
    }
}

/// 对应 `CCS::SparseAlign`:锚点对 (beginPosH, beginPosV) 列表
pub fn sparse_align(q_gram_size: usize, seq1: &[u8], seq2: &[u8]) -> Vec<(usize, usize)> {
    sparse_align_seeds(q_gram_size, seq1, seq2)
        .iter()
        .map(|s| (s.begin_pos_h() as usize, s.begin_pos_v() as usize))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn anchors_on_exact_match() {
        let consensus = b"ACGTACGTACGTTTTT";
        let read = b"AAAACGTACGTACGTTT";
        let anchors = sparse_align(6, consensus, read);
        assert!(!anchors.is_empty());
        // 锚点必须单调递增(H、V 同向)
        for w in anchors.windows(2) {
            assert!(w[0].0 < w[1].0);
            assert!(w[0].1 < w[1].1);
        }
    }

    #[test]
    fn too_short_returns_empty() {
        assert!(sparse_align(6, b"ACGTA", b"ACG").is_empty());
    }
}
