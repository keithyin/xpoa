//! Myers 位向量编辑距离 —— `cpp_code/edlib/edlib.h` 的 Rust 移植子集。
//!
//! 仅移植 FFI 路径用到的部分:`EDLIB_MODE_NW` + `EDLIB_TASK_DISTANCE`
//! (即 `EdlibGlobalAlignDistance`),k 未指定时动态放大(dynamicK)。
//! 不移植:semi-global/HW 模式、TASK_PATH/对齐重建、CIGAR。
//!
//! 移植要点:`calculateBlock` 的位运算、Ukkonen 带状块调整(条件逐行对应)、
//! `getBlockCellValues` 的分数展开顺序 —— 这些共同决定结果与 C++ 完全一致。

type Word = u64;

const WORD_SIZE: usize = 64;
const WORD_1: Word = 1;
const HIGH_BIT_MASK: Word = WORD_1 << (WORD_SIZE - 1);
/// 每隔多少列做一次强收缩(与 C++ 同值)
const STRONG_REDUCE_NUM: i32 = 2048;

#[derive(Clone, Copy, Default)]
struct Block {
    score: i32,
    p: Word,
    m: Word,
}

#[inline]
fn ceil_div(x: i32, y: i32) -> i32 {
    if x % y != 0 { x / y + 1 } else { x / y }
}

/// 对应 `Advance_Block`:计算一个 64bp 块的一列。
/// 返回 (PvOut, MvOut, hout)。
#[inline]
fn calculate_block(pv: Word, mv: Word, mut eq: Word, hin: i32) -> (Word, Word, i32) {
    // hin ∈ {1, 0, -1};-1 的补码表示使后续位运算自然生效
    let hin_is_neg = ((hin >> 2) as Word) & WORD_1;

    let xv = eq | mv;
    eq |= hin_is_neg;
    // 位向量递推(加法故意溢出,用 wrapping)
    let xh = (((eq & pv).wrapping_add(pv)) ^ pv) | eq;

    let mut ph = mv | !(xh | pv);
    let mut mh = pv & xh;

    let mut hout = ((ph & HIGH_BIT_MASK) >> (WORD_SIZE - 1)) as i32;
    hout -= ((mh & HIGH_BIT_MASK) >> (WORD_SIZE - 1)) as i32;

    ph <<= 1;
    mh <<= 1;

    mh |= hin_is_neg;
    ph |= ((hin + 1) >> 1) as Word;

    let pv_out = mh | !(xv | ph);
    let mv_out = ph & xv;

    (pv_out, mv_out, hout)
}

/// 块内 64 个 cell 的分数(自最高位起为"最底" cell),对应 `getBlockCellValues`
#[allow(clippy::needless_range_loop)] // 与 C++ getBlockCellValues 逐行对应
fn get_block_cell_values(block: Block) -> [i32; WORD_SIZE] {
    let mut scores = [0i32; WORD_SIZE];
    let mut score = block.score;
    let mut mask = HIGH_BIT_MASK;
    for i in 0..WORD_SIZE - 1 {
        scores[i] = score;
        if block.p & mask != 0 {
            score -= 1;
        }
        if block.m & mask != 0 {
            score += 1;
        }
        mask >>= 1;
    }
    scores[WORD_SIZE - 1] = score;
    scores
}

/// 对应 `transformSequences`:按首次出现顺序构建 alphabet,序列转为索引
fn transform_sequences(query: &[u8], target: &[u8]) -> (Vec<u8>, Vec<u8>, usize) {
    let mut letter_idx = [0u8; MAX_UCHAR + 1];
    let mut in_alphabet = [false; MAX_UCHAR + 1];
    let mut alphabet_len = 0usize;

    let mut q = Vec::with_capacity(query.len());
    for &c in query {
        if !in_alphabet[c as usize] {
            in_alphabet[c as usize] = true;
            letter_idx[c as usize] = alphabet_len as u8;
            alphabet_len += 1;
        }
        q.push(letter_idx[c as usize]);
    }

    let mut t = Vec::with_capacity(target.len());
    for &c in target {
        if !in_alphabet[c as usize] {
            in_alphabet[c as usize] = true;
            letter_idx[c as usize] = alphabet_len as u8;
            alphabet_len += 1;
        }
        t.push(letter_idx[c as usize]);
    }
    (q, t, alphabet_len)
}

const MAX_UCHAR: usize = 255;

/// 对应 `buildPeq`:`peq[symbol * max_num_blocks + b]`,最后一行是通配(全 1)
fn build_peq(query: &[u8], alphabet_len: usize, max_num_blocks: usize) -> Vec<Word> {
    let mut peq = vec![0u64; (alphabet_len + 1) * max_num_blocks];
    for symbol in 0..alphabet_len {
        for b in 0..max_num_blocks {
            let mut v: Word = 0;
            // 自高位向低位填充;queryLen 之后的 padding 视为通配匹配
            let mut r = (b + 1) * WORD_SIZE;
            while r > b * WORD_SIZE {
                r -= 1;
                v <<= 1;
                if r >= query.len() || query[r] as usize == symbol {
                    v += 1;
                }
            }
            peq[symbol * max_num_blocks + b] = v;
        }
    }
    for b in 0..max_num_blocks {
        peq[alphabet_len * max_num_blocks + b] = u64::MAX;
    }
    peq
}

/// 对应 `myersCalcEditDistanceNW`(findAlignment=false, targetStopPosition=-1)。
/// 返回编辑距离;无解(带状消失且未达尾列)返回 -1。
fn myers_calc_nw(
    peq: &[Word],
    w: i32,
    max_num_blocks: i32,
    query_length: i32,
    target: &[u8],
    mut k: i32,
) -> i32 {
    let target_length = target.len() as i32;

    if k < (target_length - query_length).abs() {
        return -1;
    }
    k = k.min(query_length.max(target_length)); // k 上界

    // Ukkonen 带的初始块范围
    let mut first_block: i32 = 0;
    let mut last_block: i32 = max_num_blocks.min(ceil_div(
        k.min((k + query_length - target_length) / 2) + 1,
        WORD_SIZE as i32,
    )) - 1;

    let mut blocks = vec![Block::default(); max_num_blocks as usize];
    for b in 0..=last_block {
        blocks[b as usize] = Block {
            score: (b + 1) * WORD_SIZE as i32,
            p: u64::MAX,
            m: 0,
        };
    }

    for c in 0..target_length {
        let peq_c = (target[c as usize] as usize) * max_num_blocks as usize;

        //----------------------- 计算一列 -------------------------//
        let mut hout: i32 = 1;
        for b in first_block..=last_block {
            let bl = &mut blocks[b as usize];
            let (p, m, h) = calculate_block(bl.p, bl.m, peq[peq_c + b as usize], hout);
            bl.p = p;
            bl.m = m;
            bl.score += h;
            hout = h;
        }
        let bl_score = blocks[last_block as usize].score;
        //----------------------------------------------------------------//

        // 更新 k(尾列补偿 padding)
        k = k.min(
            bl_score
                + (target_length - c - 1)
                    .max(query_length - ((1 + last_block) * WORD_SIZE as i32 - 1) - 1)
                + if last_block == max_num_blocks - 1 {
                    w
                } else {
                    0
                },
        );

        //---------- Ukkonen 带调整 ----------//
        //--- last block ---//
        if last_block + 1 < max_num_blocks
            && ((last_block + 1) * WORD_SIZE as i32 - 1
                <= k - bl_score + 2 * WORD_SIZE as i32 - 2 - target_length + c + query_length)
        {
            last_block += 1;
            let b = last_block as usize;
            let (p, m, new_hout) = calculate_block(u64::MAX, 0, peq[peq_c + b], hout);
            blocks[b].p = p;
            blocks[b].m = m;
            blocks[b].score = blocks[b - 1].score - hout + WORD_SIZE as i32 + new_hout;
            // 与 C++ 保持一致:后续循环头会重置 hout
            #[allow(unused_assignments)]
            {
                hout = new_hout;
            }
        }

        // 自底向上移除越带块
        while last_block >= first_block
            && (blocks[last_block as usize].score >= k + WORD_SIZE as i32
                || ((last_block + 1) * WORD_SIZE as i32 - 1
                    > k - blocks[last_block as usize].score + 2 * WORD_SIZE as i32
                        - 2
                        - target_length
                        + c
                        + query_length
                        + 1))
        {
            last_block -= 1;
        }
        //--- first block ---//
        while first_block <= last_block
            && (blocks[first_block as usize].score >= k + WORD_SIZE as i32
                || ((first_block + 1) * WORD_SIZE as i32 - 1
                    < blocks[first_block as usize].score - k - target_length + query_length + c))
        {
            first_block += 1;
        }

        // 周期性强收缩(与 C++ myersCalcEditDistanceNW 的强收缩逐行对应)
        if c % STRONG_REDUCE_NUM == 0 {
            while last_block >= first_block {
                let scores = get_block_cell_values(blocks[last_block as usize]);
                let num_cells = if last_block == max_num_blocks - 1 {
                    WORD_SIZE as i32 - w
                } else {
                    WORD_SIZE as i32
                };
                let mut r = last_block * WORD_SIZE as i32 + num_cells - 1;
                let mut reduce = true;
                #[allow(clippy::needless_range_loop)]
                for i in (WORD_SIZE as i32 - num_cells) as usize..WORD_SIZE {
                    if scores[i] <= k && r <= k - scores[i] - target_length + c + query_length + 1 {
                        reduce = false;
                        break;
                    }
                    r -= 1;
                }
                if !reduce {
                    break;
                }
                last_block -= 1;
            }

            while first_block <= last_block {
                let scores = get_block_cell_values(blocks[first_block as usize]);
                let num_cells = if first_block == max_num_blocks - 1 {
                    WORD_SIZE as i32 - w
                } else {
                    WORD_SIZE as i32
                };
                let mut r = first_block * WORD_SIZE as i32 + num_cells - 1;
                let mut reduce = true;
                #[allow(clippy::needless_range_loop)]
                for i in (WORD_SIZE as i32 - num_cells) as usize..WORD_SIZE {
                    if scores[i] <= k && r >= scores[i] - k - target_length + query_length + c {
                        reduce = false;
                        break;
                    }
                    r -= 1;
                }
                if !reduce {
                    break;
                }
                first_block += 1;
            }
        }

        // 带消失 → 无解
        if last_block < first_block {
            return -1;
        }
    }

    // 尾列取分数(补偿 padding 位)
    if last_block == max_num_blocks - 1 {
        let best_score = get_block_cell_values(blocks[last_block as usize])[w as usize];
        if best_score <= k {
            return best_score;
        }
    }
    -1
}

/// 对应 `edlibAlign`(NW + TASK_DISTANCE + dynamicK):返回编辑距离
pub fn edlib_distance_nw(query: &[u8], target: &[u8]) -> i32 {
    // 空序列特判(NW 模式):distance = 较长者长
    if query.is_empty() || target.is_empty() {
        return query.len().max(target.len()) as i32;
    }

    let (q, t, alphabet_len) = transform_sequences(query, target);
    let max_num_blocks = ceil_div(q.len() as i32, WORD_SIZE as i32);
    let w = (max_num_blocks as usize * WORD_SIZE - q.len()) as i32;
    let peq = build_peq(&q, alphabet_len, max_num_blocks as usize);

    // C++ 侧 config.k = -1 → dynamicK:从 WORD_SIZE 起翻倍直到有解
    let mut k = WORD_SIZE as i32;
    loop {
        let dist = myers_calc_nw(&peq, w, max_num_blocks, q.len() as i32, &t, k);
        if dist >= 0 {
            return dist;
        }
        k *= 2;
    }
}

/// 对应 `Geneus::Align::EdlibGlobalAlignDistance`:失败时返回 INT_MAX
pub fn edlib_global_align_distance(query: &[u8], target: &[u8]) -> i32 {
    let dist = edlib_distance_nw(query, target);
    if dist >= 0 { dist } else { 2147483647 }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 经典 Levenshtein DP(小规模参照实现)
    fn reference_distance(a: &[u8], b: &[u8]) -> usize {
        let mut prev: Vec<usize> = (0..=b.len()).collect();
        let mut cur = vec![0usize; b.len() + 1];
        for (i, &ca) in a.iter().enumerate() {
            cur[0] = i + 1;
            for (j, &cb) in b.iter().enumerate() {
                cur[j + 1] = (prev[j + 1] + 1)
                    .min(cur[j] + 1)
                    .min(prev[j] + usize::from(ca != cb));
            }
            std::mem::swap(&mut prev, &mut cur);
        }
        prev[b.len()]
    }

    fn check(a: &str, b: &str) {
        let expect = reference_distance(a.as_bytes(), b.as_bytes());
        let got = edlib_distance_nw(a.as_bytes(), b.as_bytes());
        assert_eq!(got as usize, expect, "query={a:?} target={b:?}");
    }

    #[test]
    fn basic_distances() {
        check("", "");
        check("", "ACGT");
        check("ACGT", "");
        check("ACGT", "ACGT");
        check("ACGT", "ACGA");
        check("AAAA", "TTTT");
        check("AAAA", "A");
        check("A", "AAAA");
        check("AAAA", "AAAAA");
        check("ACGTACGT", "ACGACGT");
        check("ACGT", "TGCA");
    }

    #[test]
    fn cross_block_lengths() {
        // 触发多块路径(> 64bp)
        let base = "ACGT".repeat(20); // 80bp,跨 2 块
        check(&base, &base);
        let mut mutated = base.clone();
        mutated.replace_range(10..12, "");
        mutated.push('A');
        check(&base, &mutated);
        check(&mutated, &base);
        // 精确 64 与 65 边界
        let s64 = "A".repeat(64);
        let s65 = format!("{s64}G");
        check(&s64, &s64);
        check(&s64, &s65);
        check(&s65, &s64);
    }

    #[test]
    fn long_targets_many_chars() {
        let mut rng = crate::casegen::DetRng::new(7);
        for _ in 0..50 {
            let qlen = 1 + rng.gen_range(0..200);
            let tlen = 1 + rng.gen_range(0..260);
            let q: Vec<u8> = (0..qlen)
                .map(|_| crate::casegen::BASES[rng.gen_range(0..4)])
                .collect();
            let mut t: Vec<u8> = (0..tlen)
                .map(|_| crate::casegen::BASES[rng.gen_range(0..4)])
                .collect();
            // 以一定概率让 target 由 query 变异而来(保证多种距离值)
            if rng.gen_bool(0.5) {
                let start = rng.gen_range(0..q.len());
                let end = (start + 1 + rng.gen_range(0..80)).min(q.len());
                t = q[start..end].to_vec();
                for b in t.iter_mut() {
                    if rng.gen_bool(0.1) {
                        *b = crate::casegen::BASES[rng.gen_range(0..4)];
                    }
                }
            }
            let expect = reference_distance(&q, &t);
            let got = edlib_distance_nw(&q, &t);
            assert_eq!(got as usize, expect, "q={} t={}", q.len(), t.len());
        }
    }

    #[test]
    fn failure_path_matches_cpp_semantics() {
        // EdlibGlobalAlignDistance 的失败值是 INT_MAX(实际 NW 路径不会失败)
        assert_eq!(edlib_global_align_distance(b"A", b"T"), 1);
        assert_eq!(edlib_global_align_distance(b"", b"AAA"), 3);
    }
}
