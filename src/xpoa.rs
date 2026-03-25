use std::ffi::CStr;

use crate::xpoa_sys::{self, PoaDraftFreeResult, PoaSetting};

pub trait TSubread {
    fn get_seq(&self) -> &str;
    fn get_cx(&self) -> u8;
    fn is_fwd(&self) -> Option<bool> {
        None
    }
}

impl TSubread for &String {
    fn get_seq(&self) -> &str {
        self
    }

    fn get_cx(&self) -> u8 {
        3
    }
}

impl TSubread for String {
    fn get_seq(&self) -> &str {
        self
    }

    fn get_cx(&self) -> u8 {
        3
    }
}

impl TSubread for &str {
    fn get_seq(&self) -> &str {
        self
    }

    fn get_cx(&self) -> u8 {
        3
    }
}

pub fn get_default_poa_setting() -> PoaSetting {
    PoaSetting {
        min_identity: 0.82,
        match_score: 3,
        mismatch_score: -5,
        insertion_score: -2,
        deletion_score: -2,
        ed_unify_strand: 1,
        version: 1,
    }
}

pub struct PoaResult {
    result: xpoa_sys::Result,
}

impl PoaResult {
    pub fn new(result: xpoa_sys::Result) -> Self {
        Self { result }
    }

    pub fn seq(&self) -> String {
        unsafe {
            CStr::from_ptr(self.result.seq)
                .to_str()
                .unwrap()
                .to_string()
        }
    }

    pub fn n_passes(&self) -> usize {
        self.result.n_passes
    }
}

impl From<xpoa_sys::Result> for PoaResult {
    fn from(value: xpoa_sys::Result) -> Self {
        Self { result: value }
    }
}

impl Drop for PoaResult {
    fn drop(&mut self) {
        unsafe {
            PoaDraftFreeResult(self.result);
        }
    }
}
