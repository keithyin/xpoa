use std::borrow::Cow;
mod dna;

pub mod casegen;
pub mod error;
pub mod native;
pub mod setting;
pub mod xpoa;

use crate::{dna::reverse_complement, error::XpoaError, error::XpoaResult, xpoa::TSubread};
pub use setting::PoaSetting;

/// 入口统一校验(对应历史 C++ 后端的"失败路径"语义,0.6.0 起单后端保留):
/// 空 read、NUL 字节、非法 version。
///
/// 这些输入在已移除的历史 C++ 后端会触发 segfault / 异常穿越 FFI 等 UB,
/// 统一在此拦截为 `Err`。
fn validate<T: TSubread>(read_infos: &[T], setting: &PoaSetting) -> XpoaResult<()> {
    for (i, read) in read_infos.iter().enumerate() {
        let seq = read.get_seq();
        if seq.is_empty() {
            return Err(XpoaError::EmptyRead { index: i });
        }
        if seq.as_bytes().contains(&0) {
            return Err(XpoaError::NulByte { index: i });
        }
    }
    if setting.version != 1 && setting.version != 2 {
        return Err(XpoaError::InvalidVersion {
            version: setting.version,
        });
    }
    Ok(())
}

/// 将 `TSubread` 集合转成 native 后端的零拷贝输入
fn native_reads<T: TSubread>(read_infos: &[T]) -> Vec<native::Read<'_>> {
    read_infos
        .iter()
        .map(|r| native::Read {
            seq: r.get_seq(),
            flags: r.get_cx() as i32,
        })
        .collect()
}

pub fn poa_consensus<T: TSubread>(
    read_infos: &[T],
    setting: &PoaSetting,
) -> XpoaResult<Option<(String, usize)>> {
    if read_infos.is_empty() {
        return Ok(None);
    }
    validate(read_infos, setting)?;
    native::poa_draft_gen(&native_reads(read_infos), setting).map(Some)
}

pub fn poa_consensus_fwd<T: TSubread>(
    read_infos: &[T],
    setting: &PoaSetting,
) -> XpoaResult<Option<(String, usize)>> {
    if read_infos.is_empty() {
        return Ok(None);
    }
    validate(read_infos, setting)?;
    native::poa_draft_gen_fwd(&native_reads(read_infos), setting).map(Some)
}

/// reads info 的 is_fwd 需要被设置
/// rust端 调整 query 的方向，然后调用 PoaDraftGenWithAllFwdStrand
pub fn poa_consensus_v2<T: TSubread>(
    read_infos: &[T],
    setting: &PoaSetting,
) -> XpoaResult<Option<(String, usize)>> {
    if read_infos.is_empty() {
        return Ok(None);
    }
    validate(read_infos, setting)?;

    // 预处理:过滤 is_fwd == None 的 read;非 forward 的取反向互补
    let prepared: Vec<(String, i32)> = read_infos
        .iter()
        .filter(|read_info| read_info.is_fwd().is_some())
        .map(|read_info| {
            let seq = if read_info.is_fwd().unwrap() {
                Cow::Borrowed(read_info.get_seq())
            } else {
                let rc = reverse_complement(read_info.get_seq().as_bytes())?;
                Cow::Owned(String::from_utf8(rc).expect("ACGTN RC stays UTF-8"))
            };
            Ok::<_, XpoaError>((seq.into_owned(), read_info.get_cx() as i32))
        })
        .collect::<XpoaResult<Vec<_>>>()?;

    let cx3_cnt = prepared.iter().filter(|(_, flags)| *flags == 3).count();

    let result = {
        let reads: Vec<native::Read> = prepared
            .iter()
            .map(|(seq, flags)| native::Read { seq, flags: *flags })
            .collect();
        native::poa_draft_gen_fwd(&reads, setting)
    };

    // v2 的既有不变量:cx3 reads 数不低于 n_passes
    if let Ok((_, npasses)) = &result {
        assert!(cx3_cnt >= *npasses, "cx3={}, npasses={}", cx3_cnt, npasses);
    }
    result.map(Some)
}

#[cfg(test)]
mod tests {
    use crate::{
        casegen::TestRead, error::XpoaError, poa_consensus, poa_consensus_fwd, poa_consensus_v2,
        xpoa::get_default_poa_setting,
    };

    fn default_setting() -> crate::PoaSetting {
        get_default_poa_setting()
    }

    #[test]
    fn test_poa_consensus() {
        let sequences = vec!["AACGGATCGGA", "AACGGATCGGA", "AACGGATCGGA", "AACGGATCGGA"];

        let res = poa_consensus(&sequences, &default_setting()).unwrap();
        assert_eq!(
            res,
            Some(("AACGGATCGGA".to_string(), 4)),
            "v1 全同 reads 的共识与 n_passes"
        );
    }

    #[test]
    fn test_poa_consensus_range_finder2() {
        let sequences = vec!["AACGGATCGGA", "AACGGATCGGA", "AACGGATCGGA", "AACGGATCGGA"];
        let mut setting = default_setting();
        setting.version = 2;

        let res = poa_consensus(&sequences, &setting).unwrap();
        assert_eq!(
            res,
            Some(("AACGGATCGGA".to_string(), 4)),
            "v2(range finder 2)全同 reads 的共识与 n_passes"
        );
    }

    #[test]
    fn test_poa_consensus_very_short() {
        let sequences = vec!["A", "C", "A", "C", "A", "A"];

        let res = poa_consensus(&sequences, &default_setting()).unwrap();
        assert_eq!(res, Some(("A".to_string(), 1)));

        let sequences = vec!["AC", "CA", "AC", "AC", "G", "T"];

        let res = poa_consensus(&sequences, &default_setting()).unwrap();
        assert_eq!(res, Some(("A".to_string(), 0)));
    }

    #[test]
    fn test_poa_consensus_fwd() {
        let sequences = vec!["AACGGATCGGA", "AACGGATCGGA", "AACGGATCGGA", "AACGGATCGGA"];
        let mut setting = default_setting();
        setting.version = 2;

        let res = poa_consensus_fwd(&sequences, &setting).unwrap();
        // 注意:FwdStrand 路径在历史 C++ 端固定使用 RangeFinderV1(见
        // PoaDraft.cpp _PoaDraftGenCoreWithAllFwdStrandCore 的默认构造
        // SparsePoa),setting.version 被忽略 —— 行为与 v1 相同,已复刻
        assert_eq!(res, Some(("AACGGATCGGA".to_string(), 4)));
    }

    #[test]
    fn test_empty_reads_rejected() {
        let sequences = vec!["ACGT", "", "ACGT"];
        let err = poa_consensus(&sequences, &default_setting()).unwrap_err();
        assert_eq!(err, XpoaError::EmptyRead { index: 1 });
    }

    #[test]
    fn test_empty_input_set_is_none() {
        let sequences: Vec<&str> = vec![];
        let res = poa_consensus(&sequences, &default_setting()).unwrap();
        assert_eq!(res, None);
    }

    #[test]
    fn test_invalid_version_rejected() {
        let sequences = vec!["ACGTACGT", "ACGTACGT"];
        let mut setting = default_setting();
        setting.version = 3;
        let err = poa_consensus(&sequences, &setting).unwrap_err();
        assert_eq!(err, XpoaError::InvalidVersion { version: 3 });
    }

    #[test]
    fn test_invalid_base_rc_rejected() {
        // 'X' 无法反向互补(历史 C++ 侧对应 throw invalid_argument)
        let mut setting = default_setting();
        setting.ed_unify_strand = 0; // OrientAndAddRead 路径会 RC
        let sequences = vec!["ACGTACGTAC", "ACGTXCGTAC", "ACGTACGTAC"];
        let err = poa_consensus(&sequences, &setting).unwrap_err();
        assert_eq!(err, XpoaError::InvalidBase { base: b'X' });
    }

    #[test]
    fn test_v2_invalid_base_rejected() {
        let mut setting = default_setting();
        setting.ed_unify_strand = 1;
        // 'X' 不在 dna.rs 互补表中,预处理 RC 时返回 Err
        let reads = vec![TestRead {
            seq: "ACGTXCGT".into(),
            flags: 3,
            is_fwd: Some(false), // 预处理会 RC
        }];
        let err = poa_consensus_v2(&reads, &setting).unwrap_err();
        assert_eq!(err, XpoaError::InvalidBase { base: b'X' });
    }
}
