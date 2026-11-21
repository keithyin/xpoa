//! `PoaSetting` 的唯一定义点。
//!
//! 纯 Rust 后端与 C++ FFI 后端共用;`repr(C)` 保持与 `cpp_code/Settings.h`
//! 中 `PoaSetting` 的 ABI 布局一致(字段顺序、大小、对齐),使该结构可直接
//! 传给 C 接口。

/// 使用 C++ 后端(libxpoa.a,需 cmake/Boost/SeqAn 构建环境)
pub const BACKEND_CPP: i32 = 0;
/// 使用纯 Rust 后端(src/native,零 C++ 依赖)
pub const BACKEND_RUST: i32 = 1;

/// POA 对齐与共识参数。`backend` 为最后一项,C++ 侧忽略其值(仅用于
/// Rust 端分发)。
#[repr(C)]
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct PoaSetting {
    pub min_identity: f32,
    pub match_score: i32,
    pub mismatch_score: i32,
    pub insertion_score: i32,
    pub deletion_score: i32,
    pub ed_unify_strand: i32,
    pub version: i32,
    pub backend: i32,
}

#[allow(clippy::identity_op)] // bindgen 风格的偏移断言表
const _: () = {
    ["Size of PoaSetting"][::std::mem::size_of::<PoaSetting>() - 32usize];
    ["Alignment of PoaSetting"][::std::mem::align_of::<PoaSetting>() - 4usize];
    ["Offset of field: PoaSetting::min_identity"]
        [::std::mem::offset_of!(PoaSetting, min_identity) - 0usize];
    ["Offset of field: PoaSetting::match_score"]
        [::std::mem::offset_of!(PoaSetting, match_score) - 4usize];
    ["Offset of field: PoaSetting::mismatch_score"]
        [::std::mem::offset_of!(PoaSetting, mismatch_score) - 8usize];
    ["Offset of field: PoaSetting::insertion_score"]
        [::std::mem::offset_of!(PoaSetting, insertion_score) - 12usize];
    ["Offset of field: PoaSetting::deletion_score"]
        [::std::mem::offset_of!(PoaSetting, deletion_score) - 16usize];
    ["Offset of field: PoaSetting::ed_unify_strand"]
        [::std::mem::offset_of!(PoaSetting, ed_unify_strand) - 20usize];
    ["Offset of field: PoaSetting::version"][::std::mem::offset_of!(PoaSetting, version) - 24usize];
    ["Offset of field: PoaSetting::backend"][::std::mem::offset_of!(PoaSetting, backend) - 28usize];
};
