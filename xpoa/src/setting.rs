//! `PoaSetting` 的唯一定义点。
//!
//! 历史说明:曾保留 `#[repr(C)]` 与 `backend` 字段以同 C++ 后端(`cpp_code/`,
//! 已于 0.6.0 移除)的 `Settings.h` 保持 ABI 对齐并支持运行时后端分发;
//! C++ 后端移除后该动机消失,`repr(C)` 仅作为布局稳定性的保留。

/// POA 对齐与共识参数。
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
}

#[allow(clippy::identity_op)] // bindgen 风格的偏移断言表
const _: () = {
    ["Size of PoaSetting"][::std::mem::size_of::<PoaSetting>() - 28usize];
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
};
