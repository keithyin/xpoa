//! 种子查找(对应 cpp_code/FindSeeds.h、Seeds.h / Seeds-inl.h)。
//!
//! `Seeds` 精确模拟 C++ `std::multiset<Seed>`(按 `operator<` = BeginDiagonal
//! 升序;等值元素保持插入序,C++11 起标准保证)。用 `BTreeMap<diagonal,
//! Vec<Seed>>` 表达:迭代序 = diagonal 升序 + 同 diagonal 插入序。

use crate::native::qgram::{QGramIndex, seed_from_hit};

#[derive(Debug, Clone, Default)]
pub struct Seeds {
    groups: std::collections::BTreeMap<i64, Vec<Seed>>,
    len: usize,
}

use crate::native::seed::Seed;

impl Seeds {
    pub fn new() -> Self {
        Self::default()
    }

    /// 对应 `AddSeed`(multiset::insert,不去重)
    pub fn add_seed(&mut self, s: Seed) {
        self.groups.entry(s.begin_diagonal()).or_default().push(s);
        self.len += 1;
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// multiset 迭代序:diagonal 升序、同 diagonal 按插入序
    pub fn iter(&self) -> impl Iterator<Item = &Seed> + '_ {
        self.groups.values().flatten()
    }
}

impl<'s> IntoIterator for &'s Seeds {
    type Item = &'s Seed;
    type IntoIter = std::vec::IntoIter<&'s Seed>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter().collect::<Vec<_>>().into_iter()
    }
}

/// 双序列种子查找(对应 `FindSeeds(qGramSize, seq1, seq2, filterHomopolymers)`)。
/// seq2 建索引(reference),seq1 为 query。
pub fn find_seeds(
    q_gram_size: usize,
    seq1: &[u8],
    seq2: &[u8],
    filter_homopolymers: bool,
) -> Seeds {
    let mut seeds = Seeds::new();
    if seq2.len() < q_gram_size {
        return seeds;
    }
    let index = QGramIndex::new(q_gram_size, seq2);
    find_seeds_from_index(&index, seq1, filter_homopolymers, &mut seeds);
    seeds
}

/// 对应 `FindSeeds(index, seq, filterHomopolymers)`(单序列 index → 取第一个
/// reference 的 Seeds,即 `multiSeeds.cbegin()->second`)
pub fn find_seeds_from_index(
    index: &QGramIndex,
    seq: &[u8],
    filter_homopolymers: bool,
    out: &mut Seeds,
) {
    for hits in index.hits(seq, filter_homopolymers) {
        let query_pos = hits.query_position();
        for hit in hits.iter() {
            out.add_seed(seed_from_hit(query_pos, hit.position, index.size()));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn multiset_order_by_diagonal() {
        let mut seeds = Seeds::new();
        seeds.add_seed(Seed::new(10, 2, 3)); // diag +8
        seeds.add_seed(Seed::new(0, 0, 3)); // diag 0
        seeds.add_seed(Seed::new(5, 0, 3)); // diag +5
        seeds.add_seed(Seed::new(5, 0, 3)); // 同 diag 重复(插入序保留)
        let diags: Vec<i64> = seeds.iter().map(|s| s.begin_diagonal()).collect();
        assert_eq!(diags, vec![0, 5, 5, 8]);
        assert_eq!(seeds.len(), 4);
    }

    #[test]
    fn short_reference_returns_empty() {
        let seeds = find_seeds(6, b"ACGTAC", b"ACG", true);
        assert!(seeds.is_empty());
    }

    #[test]
    fn finds_seeds_between_sequences() {
        let ref_seq = b"TTACGTACGATT";
        let query = b"ACGTACG";
        let seeds = find_seeds(6, query, ref_seq, true);
        // query 的 6-mers:ACGTAC(0) CGTACG(1);ref 中 ACGTAC@2, CGTACG@3
        assert_eq!(seeds.len(), 2);
        let all: Vec<(u64, u64)> = seeds
            .iter()
            .map(|s| (s.begin_pos_h(), s.begin_pos_v()))
            .collect();
        assert!(all.contains(&(0, 2)));
        assert!(all.contains(&(1, 3)));
    }

    /// 穷举参照:所有 q-mer 精确匹配,按 multiset 迭代序排列
    /// (diagonal 升序,同 diagonal 内按插入序 = h 升序)
    fn brute_seeds(q: usize, seq1: &[u8], seq2: &[u8]) -> Vec<(i64, u64, u64)> {
        let mut out = Vec::new();
        for h in 0..seq1.len().saturating_sub(q - 1) {
            for v in 0..seq2.len().saturating_sub(q - 1) {
                if seq1[h..h + q] == seq2[v..v + q] {
                    out.push(((h as i64 - v as i64), h as u64, v as u64));
                }
            }
        }
        out.sort_by_key(|&(d, h, _)| (d, h));
        out
    }

    #[test]
    fn matches_bruteforce_on_random_sequences() {
        let mut rng = crate::casegen::DetRng::new(99);
        for case in 0..40 {
            let l1 = 10 + rng.gen_range(0..120);
            let l2 = 10 + rng.gen_range(0..120);
            let seq1: Vec<u8> = (0..l1)
                .map(|_| crate::casegen::BASES[rng.gen_range(0..4)])
                .collect();
            let seq2: Vec<u8> = (0..l2)
                .map(|_| crate::casegen::BASES[rng.gen_range(0..4)])
                .collect();
            // 一半 case 让 seq2 由 seq1 切片变异而来,保证高命中
            let (seq1, seq2) = if case % 2 == 0 {
                let start = rng.gen_range(0..l1 / 2);
                let end = (start + 20).min(l1);
                let mut s2 = seq1[start..end].to_vec();
                for b in s2.iter_mut() {
                    if rng.gen_bool(0.15) {
                        *b = crate::casegen::BASES[rng.gen_range(0..4)];
                    }
                }
                (seq1, s2)
            } else {
                (seq1, seq2)
            };

            let expect: Vec<(u64, u64)> = brute_seeds(6, &seq1, &seq2)
                .into_iter()
                .map(|(_, h, v)| (h, v))
                .collect();
            let got: Vec<(u64, u64)> = find_seeds(6, &seq1, &seq2, false)
                .iter()
                .map(|s| (s.begin_pos_h(), s.begin_pos_v()))
                .collect();
            assert_eq!(got, expect, "case {case}: q={} r={}", l1, l2);
        }
    }
}
