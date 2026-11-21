//! DNA 碱基互补(对应 cpp_code/Sequence.cpp)。
//!
//! 精确复制 C++ 的 IUPAC 全表 —— 包括 V/G/H/D/M/K/N/Y/S/A/B/W/R 等模糊码。
//! 注意:这与 `src/dna.rs` 的 ACGTN 表不同;`dna.rs` 服务于 `poa_consensus_v2`
//! 的既有行为,本模块服务于 native 后端内部(ed_unify_strand / OrientAndAdd)。

/// 对应 `PacBio::Data::Complement`:非法字符返回 `None`(C++ 抛异常)
#[inline]
pub fn complement(base: u8) -> Option<u8> {
    const LOOKUP: [u8; 256] = {
        let mut t = [0u8; 256];
        // 45 == '-'
        t[45] = b'-';
        // 大写(65-89),与 Sequence.cpp 的 64-95 段一致
        t[65] = b'T'; // A -> T
        t[66] = b'V'; // B -> V
        t[67] = b'G'; // C -> G
        t[68] = b'H'; // D -> H
        t[71] = b'C'; // G -> C
        t[72] = b'D'; // H -> D
        t[75] = b'M'; // K -> M
        t[77] = b'K'; // M -> K
        t[78] = b'N'; // N -> N
        t[82] = b'Y'; // R -> Y
        t[83] = b'S'; // S -> S
        t[84] = b'A'; // T -> A
        t[86] = b'B'; // V -> B
        t[87] = b'W'; // W -> W
        t[89] = b'R'; // Y -> R
        // 小写(97-121),96-127 段同构
        t[97] = b't'; // a
        t[98] = b'v'; // b
        t[99] = b'g'; // c
        t[100] = b'h'; // d
        t[103] = b'c'; // g
        t[104] = b'd'; // h
        t[107] = b'm'; // k
        t[109] = b'k'; // m
        t[110] = b'n'; // n
        t[114] = b'y'; // r
        t[115] = b's'; // s
        t[116] = b'a'; // t
        t[118] = b'b'; // v
        t[119] = b'w'; // w
        t[121] = b'r'; // y
        t
    };
    // 大写 Y(89) 的补是 'R';上面块内注释已列,直接以 C++ 逐字节表为准
    match LOOKUP[base as usize] {
        0 => None,
        c => Some(c),
    }
}

/// 对应 `PacBio::Data::ReverseComplement`:反向 + 互补。
/// 非法碱基(C++ 抛 `invalid_argument`)返回 `Err`。
pub fn reverse_complement(seq: &[u8]) -> Result<Vec<u8>, crate::error::XpoaError> {
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        match complement(b) {
            Some(c) => out.push(c),
            None => return Err(crate::error::XpoaError::InvalidBase { base: b }),
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn acgt_roundtrip() {
        assert_eq!(reverse_complement(b"ACGTAA").unwrap(), b"TTACGT");
        assert_eq!(reverse_complement(b"ACGT").unwrap(), b"ACGT");
        assert_eq!(reverse_complement(b"A").unwrap(), b"T");
    }

    #[test]
    fn iupac_codes_match_cpp_table() {
        // 与 Sequence.cpp LOOKUP 逐项对拍
        assert_eq!(complement(b'A'), Some(b'T'));
        assert_eq!(complement(b'C'), Some(b'G'));
        assert_eq!(complement(b'G'), Some(b'C'));
        assert_eq!(complement(b'T'), Some(b'A'));
        assert_eq!(complement(b'-'), Some(b'-'));
        assert_eq!(complement(b'N'), Some(b'N'));
        assert_eq!(complement(b'B'), Some(b'V'));
        assert_eq!(complement(b'D'), Some(b'H'));
        assert_eq!(complement(b'H'), Some(b'D'));
        assert_eq!(complement(b'V'), Some(b'B'));
        assert_eq!(complement(b'M'), Some(b'K'));
        assert_eq!(complement(b'K'), Some(b'M'));
        assert_eq!(complement(b'S'), Some(b'S'));
        assert_eq!(complement(b'W'), Some(b'W'));
        assert_eq!(complement(b'Y'), Some(b'R'));
        assert_eq!(complement(b'R'), Some(b'Y'));
        assert_eq!(complement(b'a'), Some(b't'));
        assert_eq!(complement(b'y'), Some(b'r'));
        assert_eq!(complement(b'X'), None);
        assert_eq!(complement(b'1'), None);
    }
}
