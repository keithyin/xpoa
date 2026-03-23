use std::ffi::{CString, c_int};
mod xpoa_sys;

use crate::{
    xpoa::{PoaResult, TSubread},
    xpoa_sys::{PoaDraftGen, PoaSetting, Subread},
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
