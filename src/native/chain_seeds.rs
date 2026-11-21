//! SDP 种子链接(对应 cpp_code/ChainSeeds.h 的移植子集)。
//!
//! 只移植 FFI 路径可达的部分:`ChainSeeds(seedSet, config)` 单 seedSet、
//! `numCandidates = 1`(SparseAlignSeeds 的配置)。多 seedSet 重载
//! (`ChainedSeedSets`、map 版)未使用,不移植。
//!
//! 精确性关键点(与 C++ 逐一对应):
//! 1. `std::set<SDPHit>`(sweepSet)按 (BeginDiagonal, BeginPositionH) 比较,
//!    **等价键拒绝插入** —— 用 BTreeMap + `entry().or_insert()` 复刻;
//! 2. `std::set<SDPColumn>`(colSet)**仅按 EndPositionH 比较**,等价拒绝,
//!    且 `insert` 失败时 `std::next(已存在元素)` 参与后续 erase 循环;
//! 3. `lower_bound`/`upper_bound`/`prev` 游标语义用 BTreeMap range 等价表达;
//! 4. `priority_queue` 在 `numCandidates == 1` 时退化为"严格大于才替换"
//!    (first-wins);
//! 5. HV/VH 排序的 tie 仅发生在完全相同的 seed 之间(size 恒 6),任意
//!    合法排序(这里用 stable sort)与 `std::sort` 等价。

use std::collections::BTreeMap;
use std::ops::Bound::{Excluded, Unbounded};

use crate::native::find_seeds::Seeds;
use crate::native::seed::Seed;

/// 对应 `ChainSeedsConfig`
#[derive(Debug, Clone, Copy)]
pub struct ChainSeedsConfig {
    pub num_candidates: usize,
    pub min_score: i64,
    pub match_score: i64,
    pub non_match_penalty: i64,
    pub insertion_penalty: i64,
    pub deletion_penalty: i64,
    pub max_seed_gap: i64,
}

impl ChainSeedsConfig {
    /// SparseAlignSeeds 使用的配置:`{1, 1, 3, -1, -1, -1, INT_MAX}`
    pub fn sparse_align() -> Self {
        Self {
            num_candidates: 1,
            min_score: 1,
            match_score: 3,
            non_match_penalty: -1,
            insertion_penalty: -1,
            deletion_penalty: -1,
            max_seed_gap: i32::MAX as i64,
        }
    }
}

/// 对应匿名空间 `LinkScore`
fn link_score(lhs: &Seed, rhs: &Seed, config: &ChainSeedsConfig) -> i64 {
    let l_h = lhs.begin_pos_h() as i64;
    let l_v = lhs.begin_pos_v() as i64;
    let r_h = rhs.begin_pos_h() as i64;
    let r_v = rhs.begin_pos_v() as i64;
    let k = (lhs.size().min(rhs.size())) as i64;
    let fwd = (l_h - r_h).min(l_v - r_v);

    // matchReward = anchor bases * matchScore
    let matches = k - 0i64.max(k - fwd);
    let match_reward = matches * config.match_score;

    // nonMatchPenalty = on-diagonal bases between seeds * penalty
    let non_matches = fwd - matches;
    if non_matches > config.max_seed_gap {
        return -1;
    }
    let non_match_penalty = non_matches * config.non_match_penalty;

    // indelPenalty = diagonal drift * penalty
    let drift = lhs.diagonal() - rhs.diagonal();
    let mut indel_penalty = 0i64;
    if drift > 0 {
        indel_penalty = drift * config.insertion_penalty;
    } else if drift < 0 {
        indel_penalty = -drift * config.deletion_penalty;
    }

    match_reward + non_match_penalty + indel_penalty
}

#[inline]
fn hv_compare(a: &Seed, b: &Seed) -> std::cmp::Ordering {
    a.begin_pos_h()
        .cmp(&b.begin_pos_h())
        .then(a.end_pos_v().cmp(&b.end_pos_v()))
}

#[inline]
fn vh_compare(a: &Seed, b: &Seed) -> std::cmp::Ordering {
    a.begin_pos_v()
        .cmp(&b.begin_pos_v())
        .then(a.end_pos_h().cmp(&b.end_pos_h()))
}

/// sweepSet 键:(BeginDiagonal, BeginPositionH) —— DiagonalCompare
type SweepKey = (i64, u64);

#[inline]
fn sweep_key(s: &Seed) -> SweepKey {
    (s.begin_diagonal(), s.begin_pos_h())
}

#[inline]
fn z_score(scores: &[i64], idx: usize, s: &Seed) -> i64 {
    scores[idx] + s.begin_pos_h() as i64 + s.begin_pos_v() as i64
}

/// 对应 `ChainSeeds`(单 seedSet,任意 numCandidates;当前调用方固定为 1)。
/// 返回按分数降序的候选链。
pub fn chain_seeds(seed_set: &Seeds, config: &ChainSeedsConfig) -> Vec<Vec<Seed>> {
    let num_candidates = config.num_candidates;
    let seeds: Vec<Seed> = seed_set.iter().copied().collect();
    let n = seeds.len();
    if n == 0 {
        return Vec::new();
    }

    let mut scores: Vec<i64> = seeds.iter().map(|s| s.size() as i64).collect();
    let mut chain_pred: Vec<Option<usize>> = vec![None; n];

    // 候选链:priority_queue<ChainHit> 在 numCandidates == 1 时退化为
    // "严格大于才替换"(first-wins)。本模块按移植范围只支持该配置,
    // C++ 的多候选堆语义不在 FFI 可达路径上。
    assert!(
        num_candidates == 1,
        "chain_seeds 仅支持 numCandidates == 1(SparseAlignSeeds 配置)"
    );
    let mut best: Option<(i64, usize)> = None; // (score, endIdx)
    #[inline]
    fn heap_push(best: &mut Option<(i64, usize)>, score: i64, idx: usize) {
        if best.is_none() || score > best.unwrap().0 {
            *best = Some((score, idx));
        }
    }

    // ---------------- ComputeVisibilityLeft(HV 排序) ----------------
    let mut h_order: Vec<usize> = (0..n).collect();
    h_order.sort_by(|&a, &b| hv_compare(&seeds[a], &seeds[b]));
    let mut visible: Vec<Option<usize>> = vec![None; n];
    {
        let mut sweep: BTreeMap<SweepKey, usize> = BTreeMap::new();
        let mut to_remove = 0usize;
        let mut it = 0usize;
        while it < n {
            let col = seeds[h_order[it]].begin_pos_h();
            let start = it;
            while it < n && seeds[h_order[it]].begin_pos_h() == col {
                let idx = h_order[it];
                // upper_bound:第一个严格大于当前键的 sweepSet 元素
                let key = sweep_key(&seeds[idx]);
                if let Some((_, &vi)) = sweep.range((Excluded(key), Unbounded)).next() {
                    visible[idx] = Some(vi);
                }
                it += 1;
            }
            // std::set::insert(first, last):等价键拒绝
            for &idx in &h_order[start..it] {
                sweep.entry(sweep_key(&seeds[idx])).or_insert(idx);
            }
            // 移除结束位置早于当前列的种子
            while to_remove < n && seeds[h_order[to_remove]].end_pos_h() < col {
                sweep.remove(&sweep_key(&seeds[h_order[to_remove]]));
                to_remove += 1;
            }
        }
        // sweepSet.clear() 在 VH 阶段前重新开始
    }

    // ---------------- 主扫描(VH 排序) ----------------
    let mut v_order: Vec<usize> = (0..n).collect();
    v_order.sort_by(|&a, &b| vh_compare(&seeds[a], &seeds[b]));
    let mut sweep: BTreeMap<SweepKey, usize> = BTreeMap::new();
    let mut col_set: BTreeMap<u64, usize> = BTreeMap::new(); // key = EndPositionH
    let mut to_remove = 0usize;
    let mut it = 0usize;
    while it < n {
        let row = seeds[v_order[it]].begin_pos_v();
        let start = it;
        while it < n && seeds[v_order[it]].begin_pos_v() == row {
            let idx = v_order[it];
            let s = &seeds[idx];
            let mut best_score = i64::MIN;
            let mut best_seed: Option<usize> = None;

            // 1) column predecessor:colSet.lower_bound(BeginH) 的 prev
            {
                let begin_h = s.begin_pos_h();
                if let Some((&_, &pi)) = col_set.range(..begin_h).next_back() {
                    let sc = scores[pi] + link_score(s, &seeds[pi], config);
                    if sc > best_score {
                        best_score = sc;
                        best_seed = Some(pi);
                    }
                }
            }
            // 2) visible above:sweepSet.lower_bound(seed) 的 prev
            {
                let key = sweep_key(s);
                if let Some((&_, &vi)) = sweep.range(..key).next_back() {
                    let sc = scores[vi] + link_score(s, &seeds[vi], config);
                    if sc > best_score {
                        best_score = sc;
                        best_seed = Some(vi);
                    }
                }
            }
            // 3) visible left
            if let Some(li) = visible[idx] {
                let sc = scores[li] + link_score(s, &seeds[li], config);
                if sc > best_score {
                    best_score = sc;
                    best_seed = Some(li);
                }
            }

            if best_seed.is_some() && best_score >= config.min_score {
                scores[idx] = best_score;
                chain_pred[idx] = best_seed;
                heap_push(&mut best, best_score, idx);
            } else if scores[idx] >= config.min_score {
                // C++ 原样:push 的 score 是本轮 bestScore(可能为 i64::MIN)
                heap_push(&mut best, best_score, idx);
            }
            it += 1;
        }

        // 当前行种子进入 sweepSet(等价拒绝)
        for &idx in &v_order[start..it] {
            sweep.entry(sweep_key(&seeds[idx])).or_insert(idx);
        }

        // 移除结束位置早于当前行的种子,并同步维护 colSet
        while to_remove < n && seeds[v_order[to_remove]].end_pos_v() < row {
            let idx = v_order[to_remove];
            let r = &seeds[idx];
            let col_key = r.end_pos_h();
            let z_remove = z_score(&scores, idx, r);

            // find(col):等价键判定
            let should_update = match col_set.get(&col_key) {
                None => true,
                Some(&existing) => z_score(&scores, existing, &seeds[existing]) < z_remove,
            };
            if should_update {
                // insert(等价拒绝:存在则保留原值,仅参与后续 erase 循环)
                col_set.entry(col_key).or_insert(idx);
                // it = next(insert 的位置);持续擦除 zScore 更小的后继列
                loop {
                    let next = col_set
                        .range((Excluded(col_key), Unbounded))
                        .next()
                        .map(|(k, v)| (*k, *v));
                    match next {
                        Some((k, vi)) if z_score(&scores, vi, &seeds[vi]) < z_remove => {
                            col_set.remove(&k);
                        }
                        _ => break,
                    }
                }
            }
            sweep.remove(&sweep_key(r));
            to_remove += 1;
        }
    }

    // ---------------- 重建链(分数降序) ----------------
    let mut chains = Vec::new();
    if let Some((_, end_idx)) = best {
        let mut chain = Vec::new();
        let mut cur = Some(end_idx);
        while let Some(i) = cur {
            chain.push(seeds[i]);
            cur = chain_pred[i];
        }
        chain.reverse();
        chains.push(chain);
    }
    chains
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native::find_seeds::find_seeds;

    #[test]
    fn link_score_math() {
        let config = ChainSeedsConfig::sparse_align();
        let a = Seed::new(0, 0, 6);
        let b = Seed::new(6, 6, 6); // 完美共对角
        // lhs=b, rhs=a:fwd = 6,matches = 6,reward = 18;nonMatches = 0;drift = 0
        assert_eq!(link_score(&b, &a, &config), 18);
        // lhs=a, rhs=b:fwd = min(-6,-6) = -6,matches = 6 - max(0, 12) = -6,
        // reward = -18;nonMatches = -6 - (-6) = 0;drift = 0 → -18
        assert_eq!(link_score(&a, &b, &config), -18);
    }

    #[test]
    fn chain_single_diagonal_seeds() {
        // 同对角线上的 3 个完美种子应链成一条
        let ref_seq = b"ACGTACGTACGTACGT";
        let query = b"ACGTACGT";
        let seeds = find_seeds(6, query, ref_seq, true);
        let chains = chain_seeds(&seeds, &ChainSeedsConfig::sparse_align());
        assert_eq!(chains.len(), 1);
        let c = &chains[0];
        assert!(!c.is_empty());
        // 链必须自左向右单调
        for w in c.windows(2) {
            assert!(w[0].begin_pos_h() < w[1].begin_pos_h());
        }
    }

    #[test]
    fn empty_seeds_return_empty() {
        let chains = chain_seeds(&Seeds::new(), &ChainSeedsConfig::sparse_align());
        assert!(chains.is_empty());
    }
}
