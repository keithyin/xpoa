use crate::setting::{BACKEND_CPP, BACKEND_RUST, PoaSetting};

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
        version: 1, // version 2 意味着使用 range_finder2
        // 默认返回当前构建下可用的后端:双后端构建默认 C++(迁移期保持既有
        // 行为);纯 Rust 构建(default-features=false)只能是 Rust 后端
        backend: if cfg!(feature = "cpp") {
            BACKEND_CPP
        } else {
            BACKEND_RUST
        },
    }
}

#[cfg(feature = "cpp")]
pub struct PoaResult {
    result: crate::xpoa_sys::Result,
}

#[cfg(feature = "cpp")]
impl PoaResult {
    pub fn new(result: crate::xpoa_sys::Result) -> Self {
        Self { result }
    }

    pub fn seq(&self) -> String {
        unsafe {
            std::ffi::CStr::from_ptr(self.result.seq)
                .to_str()
                .unwrap()
                .to_string()
        }
    }

    pub fn n_passes(&self) -> usize {
        self.result.n_passes
    }
}

#[cfg(feature = "cpp")]
impl From<crate::xpoa_sys::Result> for PoaResult {
    fn from(value: crate::xpoa_sys::Result) -> Self {
        Self { result: value }
    }
}

#[cfg(feature = "cpp")]
impl Drop for PoaResult {
    fn drop(&mut self) {
        unsafe {
            crate::xpoa_sys::PoaDraftFreeResult(self.result);
        }
    }
}
