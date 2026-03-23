use std::{
    borrow::Cow,
    ffi::{CString, c_int},
};
mod dna;
mod xpoa_sys;

use crate::{
    dna::reverse_complement,
    xpoa::{PoaResult, TSubread},
    xpoa_sys::{PoaDraftGen, PoaDraftGenWithAllFwdStrand, PoaSetting, Subread},
};
pub mod xpoa;

pub fn poa_consensus<T: TSubread>(
    read_infos: &[T],
    setting: &PoaSetting,
) -> Option<(String, usize)> {
    if read_infos.is_empty() {
        return None;
    }
    let read_infos = read_infos
        .iter()
        .map(|read_info| {
            (
                CString::new(read_info.get_seq()).unwrap(),
                read_info.get_cx() as c_int,
            )
        })
        .collect::<Vec<_>>();
    let reads = read_infos
        .iter()
        .map(|(seq, flag)| Subread {
            seq: seq.as_ptr() as *mut i8,
            flags: *flag as c_int,
        })
        .collect::<Vec<_>>();

    unsafe {
        let poa_res: PoaResult = PoaDraftGen(reads.as_ptr(), reads.len(), setting).into();
        let seq = poa_res.seq();
        let npasses = poa_res.n_passes();
        // FreePoaResult(poa_res);
        Some((seq, npasses))
    }
}

pub fn poa_consensus_fwd<T: TSubread>(
    read_infos: &[T],
    setting: &PoaSetting,
) -> Option<(String, usize)> {
    if read_infos.is_empty() {
        return None;
    }
    let read_infos = read_infos
        .iter()
        .map(|read_info| {
            (
                CString::new(read_info.get_seq()).unwrap(),
                read_info.get_cx() as c_int,
            )
        })
        .collect::<Vec<_>>();
    let reads = read_infos
        .iter()
        .map(|(seq, flag)| Subread {
            seq: seq.as_ptr() as *mut i8,
            flags: *flag as c_int,
        })
        .collect::<Vec<_>>();

    unsafe {
        let poa_res: PoaResult =
            PoaDraftGenWithAllFwdStrand(reads.as_ptr(), reads.len(), setting).into();
        let seq = poa_res.seq();
        let npasses = poa_res.n_passes();
        // FreePoaResult(poa_res);
        Some((seq, npasses))
    }
}

pub fn poa_consensus_v2<T: TSubread>(
    read_infos: &[T],
    setting: &PoaSetting,
) -> Option<(String, usize)> {
    if read_infos.is_empty() {
        return None;
    }
    let read_infos = read_infos
        .iter()
        .filter(|read_info| read_info.is_fwd().is_some())
        .map(|read_info| {
            let seq = if read_info.is_fwd().unwrap() {
                Cow::Borrowed(read_info.get_seq())
            } else {
                Cow::Owned(
                    String::from_utf8(reverse_complement(read_info.get_seq().as_bytes())).unwrap(),
                )
            };

            (
                CString::new(seq.as_bytes()).unwrap(),
                read_info.get_cx() as c_int,
            )
        })
        .collect::<Vec<_>>();
    let reads = read_infos
        .iter()
        .map(|(seq, flag)| Subread {
            seq: seq.as_ptr() as *mut i8,
            flags: *flag as c_int,
        })
        .collect::<Vec<_>>();

    let cx3_cnt = reads.iter().filter(|info| info.flags == 3).count();

    unsafe {
        let poa_res: PoaResult =
            PoaDraftGenWithAllFwdStrand(reads.as_ptr(), reads.len(), setting).into();
        // let seq = CString::from_raw(poa_res.seq).to_str().unwrap().to_string();
        let seq = poa_res.seq();
        let npasses = poa_res.n_passes();
        // FreePoaResult(poa_res);
        assert!(cx3_cnt >= npasses, "cx3={}, npasses={}", cx3_cnt, npasses);
        Some((seq, npasses))
    }
}

#[cfg(test)]
mod tests {
    use crate::{poa_consensus, xpoa::get_default_poa_setting};

    #[test]
    fn test_poa_consensus() {
        let sequences = vec!["AACGGATCGGA", "AACGGATCGGA", "AACGGATCGGA", "AACGGATCGGA"];

        let res = poa_consensus(&sequences, &get_default_poa_setting());
        println!("res:{:?}", res);
    }
}
