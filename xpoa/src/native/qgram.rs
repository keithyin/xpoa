//! q-gram 索引(对应 cpp_code/Index.h、Index-inl.h、Hashing-inl.h、
//! IndexHit.h、IndexHits.h 的移植子集,单序列索引)。

use crate::native::seed::Seed;

const ALPHABET_SIZE: u64 = 4; // {A, C, G, T}

/// 非 ACGT 字符(含 N、小写、任意其他字节)一律映射为 0(与 C++ 一致)
#[inline]
pub fn base_code(c: u8) -> u64 {
    match c {
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 0,
    }
}

/// 滚动 q-gram 哈希器(对应 `QGram::internal::Shape`)
pub struct Shape<'a> {
    q: usize,
    hash_factor: u64, // 4^(q-1)
    seq: &'a [u8],
    iter: usize,
    current_hash: u64,
    left_char: u8,
}

impl<'a> Shape<'a> {
    /// 构造时只计算前 q-1 个字符的"伪窗口"哈希(C++ `HashImpl(BaseCode(*iter_),
    /// iter_, q - 1)`);首个 `hash_next()` 补入第 q 个字符,得到第一个完整
    /// q-mer 的哈希。`seq.len() >= q` 且 `1 <= q <= 16`。
    pub fn new(q: usize, seq: &'a [u8]) -> Self {
        assert!(seq.len() >= q, "sequence must be >= q");
        assert!((1..=16).contains(&q), "q must be in [1,16]");
        let mut h = base_code(seq[0]);
        for c in seq.iter().take(q - 1).skip(1) {
            h = h * ALPHABET_SIZE + base_code(*c);
        }
        Self {
            q,
            hash_factor: {
                // C++ 用 std::pow(4, q-1);此处用整数幂(结果精确一致)
                let mut f = 1u64;
                for _ in 1..q {
                    f *= ALPHABET_SIZE;
                }
                f
            },
            seq,
            iter: 0,
            current_hash: h,
            left_char: 0, // '\0' → base_code = 0,首个 HashNext 不减项
        }
    }

    /// 推进并返回下一个 q-mer 的哈希;可安全调用 `seq.len() - q + 1` 次
    pub fn hash_next(&mut self) -> u64 {
        self.current_hash = (self.current_hash - base_code(self.left_char) * self.hash_factor)
            * ALPHABET_SIZE
            + base_code(self.seq[self.iter + self.q - 1]);
        self.left_char = self.seq[self.iter];
        self.iter += 1;
        self.current_hash
    }
}

/// 同聚物哈希过滤器(对应 `HpHasher`)
pub struct HpHasher {
    hashes: [u64; 4],
}

impl HpHasher {
    pub fn new(q: usize) -> Self {
        let dna: [u8; 4] = *b"ACGT";
        let mut hashes = [0u64; 4];
        for (i, &base) in dna.iter().enumerate() {
            let mut h = base_code(base);
            for _ in 1..q {
                h = h * ALPHABET_SIZE + base_code(base);
            }
            hashes[i] = h;
        }
        Self { hashes }
    }

    #[inline]
    pub fn is_homopolymer(&self, h: u64) -> bool {
        self.hashes.contains(&h)
    }
}

/// q-gram 命中(对应 `IndexHit`)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IndexHit {
    pub id: u32,
    pub position: u64,
}

/// 索引中一段连续命中的视图(对应 `IndexHits`)
pub struct IndexHits<'a> {
    source: &'a [IndexHit],
    begin: usize,
    end: usize,
    query_pos: usize,
}

impl<'a> IndexHits<'a> {
    #[inline]
    pub fn query_position(&self) -> usize {
        self.query_pos
    }

    pub fn iter(&self) -> impl Iterator<Item = IndexHit> + '_ {
        self.source[self.begin..self.end].iter().copied()
    }
}

/// q-gram 查找索引(对应 `QGram::Index`,单序列)
pub struct QGramIndex {
    q: usize,
    /// 按 (id, position) 排序的命中数组,经 hash_lookup 分段
    suffix_array: Vec<IndexHit>,
    /// hash → suffix_array 起点的前缀和表(长度 4^q + 1)
    hash_lookup: Vec<u64>,
}

impl QGramIndex {
    /// 对应 `Index(q, seq)`
    pub fn new(q: usize, seq: &[u8]) -> Self {
        Self::new_multi(q, std::iter::once(seq))
    }

    /// 对应 `Index(q, seqs)`(多序列;本工程实际只用单序列,保留多序列以
    /// 完整对应 Init 的计数逻辑)
    pub fn new_multi<'s, I: IntoIterator<Item = &'s [u8]>>(q: usize, seqs: I) -> Self {
        assert!((1..=16).contains(&q), "qgram size must be in [1,16]");
        let seqs: Vec<&[u8]> = seqs.into_iter().collect();
        for s in &seqs {
            assert!(s.len() >= q, "sequence size must be >= q");
        }

        // 计数 + 累计和(hash_lookup[i] = Σ_{j<=i-2} counts[j],滞后 2;
        // 填充时用 [hash+1] 作起点,即标准滞后 1 前缀和)
        let lookup_size = (ALPHABET_SIZE.pow(q as u32) + 1) as usize;
        let mut hash_lookup = vec![0u64; lookup_size];
        let mut total_num_qgrams: u32 = 0;
        for seq in &seqs {
            let num_qgrams = seq.len() - q + 1;
            let mut shape = Shape::new(q, seq);
            for _ in 0..num_qgrams {
                hash_lookup[shape.hash_next() as usize] += 1;
            }
            total_num_qgrams += num_qgrams as u32;
        }

        let mut prev_diff = 0u64;
        let mut prev_diff2 = 0u64;
        let mut sum = 0u64;
        for h in hash_lookup.iter_mut() {
            sum += prev_diff2;
            prev_diff2 = prev_diff;
            prev_diff = *h;
            *h = sum;
        }

        // 计数排序填充 suffix_array
        let mut suffix_array = vec![IndexHit { id: 0, position: 0 }; total_num_qgrams as usize];
        // 与 C++ Init 的 seqNo 递增结构对应
        #[allow(clippy::explicit_counter_loop)]
        let mut seq_no: u32 = 0;
        #[allow(clippy::explicit_counter_loop)]
        for seq in &seqs {
            let num_qgrams = seq.len() - q + 1;
            let mut shape = Shape::new(q, seq);
            for i in 0..num_qgrams {
                let h = shape.hash_next();
                let slot = hash_lookup[(h + 1) as usize] as usize;
                suffix_array[slot] = IndexHit {
                    id: seq_no,
                    position: i as u64,
                };
                hash_lookup[(h + 1) as usize] += 1;
            }
            seq_no += 1;
        }

        Self {
            q,
            suffix_array,
            hash_lookup,
        }
    }

    #[inline]
    pub fn size(&self) -> usize {
        self.q
    }

    /// 对应 `Index::Hits(seq, filterHomopolymers)`:每个 query 位置一段视图
    pub fn hits(&self, seq: &[u8], filter_homopolymers: bool) -> Vec<IndexHits<'_>> {
        let mut result = Vec::new();
        if seq.len() < self.q {
            return result;
        }
        let end = seq.len() + 1 - self.q;
        let mut shape = Shape::new(self.q, seq);
        let hp = HpHasher::new(self.q);
        for i in 0..end {
            let h = shape.hash_next();
            if !filter_homopolymers || !hp.is_homopolymer(h) {
                result.push(self.hits_for_hash(h, i));
            }
        }
        result
    }

    #[inline]
    fn hits_for_hash(&self, hash: u64, query_pos: usize) -> IndexHits<'_> {
        let begin = self.hash_lookup[hash as usize] as usize;
        let end = self.hash_lookup[(hash + 1) as usize] as usize;
        IndexHits {
            source: &self.suffix_array,
            begin,
            end,
            query_pos,
        }
    }
}

/// 为单个 (queryPos, hit) 对构造 seed(对应 FindSeeds 里的 `Seed{queryPos,
/// hit.Position(), index.Size()}`)
#[inline]
pub fn seed_from_hit(query_pos: usize, hit_pos: u64, q: usize) -> Seed {
    Seed::new(query_pos as u64, hit_pos, q as u64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shape_rolling_hash() {
        // q=6, seq="AAAAACA"(7bp → 2 个 q-mer):base_code A=0, C=1
        // q-mer0 "AAAAAC" = 1;首个 hash_next 补入 seq[5]='C'
        let mut s = Shape::new(6, b"AAAAACA");
        assert_eq!(s.hash_next(), 1, "q-mer0 AAAAAC");
        assert_eq!(s.hash_next(), 4, "q-mer1 AAAACA = 1*4 + 0");
        // 长度 7,q=6 只有 2 个 q-mer;end = len+1-q = 2,不可再推进
    }

    #[test]
    fn index_finds_all_hits() {
        let ref_seq = b"ACGTACGTAC";
        let index = QGramIndex::new(6, ref_seq);
        assert_eq!(index.size(), 6);
        let query = b"ACGTAC";
        let hits = index.hits(query, false);
        // query 有 1 个 q-mer(ACGTAC),在 ref 的 0、4 两处出现
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].query_position(), 0);
        let positions: Vec<u64> = hits[0].iter().map(|h| h.position).collect();
        assert_eq!(positions, vec![0, 4]);
    }

    #[test]
    fn homopolymer_filter() {
        let index = QGramIndex::new(6, b"AAAAAAACCC");
        let query = b"AAAAAA";
        let unfiltered = index.hits(query, false);
        let filtered = index.hits(query, true);
        assert_eq!(unfiltered.len(), 1);
        assert!(filtered.is_empty(), "纯 A 的 q-mer 应被同聚物过滤");
    }

    #[test]
    fn non_acgt_maps_to_a() {
        let index = QGramIndex::new(4, b"ACGTNNNNACGT");
        let query = b"ACGTNNNN"; // N 按 0(A)处理
        let hits = index.hits(query, false);
        assert!(!hits.is_empty());
    }
}
