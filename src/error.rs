//! 错误类型(单一定义点)。
//!
//! 对应 C++ 侧的"失败路径"语义,并消除其穿越 FFI 的 UB:
//!
//! - 空 read:C++ `PoaDraftGen` 会 segfault(`threadFirstRead` 的 UB)
//! - 非法 version:C++ `SparsePoa(int version)` 会 throw 穿越 extern "C"(UB)
//! - 非法碱基:C++ `ReverseComplement` 会 throw `invalid_argument`(UB)
//!
//! 以上在两后端统一为返回 `Err`,由 Rust 层在进入任何后端之前完成校验。

use thiserror::Error;

/// POA 计算错误。公共 API 返回 `Result<Option<(String, usize)>, XpoaError>`:
/// `Ok(None)` 表示空 read 集合(历史语义),`Err` 表示输入或配置错误。
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum XpoaError {
    /// 空 read 序列。C++ 后端对此会 segfault,两后端统一拒绝。
    #[error("empty read sequence at index {index}")]
    EmptyRead { index: usize },

    /// read 含 NUL 字节,无法穿越 C FFI 边界(仅影响 cpp 后端的传参,
    /// 统一在入口校验以保证两后端行为一致)。
    #[error("read sequence at index {index} contains a NUL byte")]
    NulByte { index: usize },

    /// 反向互补时遇到无法处理的碱基(C++ `Sequence.cpp` 对应 throw)。
    #[error("invalid base 0x{base:02X} cannot be reverse-complemented")]
    InvalidBase { base: u8 },

    /// POA version 仅支持 1(V1 RangeFinder)与 2(V2)。
    #[error("POA version must be 1 or 2, got {version}")]
    InvalidVersion { version: i32 },

    /// 在 `--no-default-features`(纯 Rust)构建下选择了 C++ 后端。
    #[error("cpp backend not compiled; rebuild with --features cpp")]
    CxxBackendUnavailable,
}

/// 便捷别名
pub type XpoaResult<T> = Result<T, XpoaError>;
