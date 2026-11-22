//! 错误类型(单一定义点)。
//!
//! 对应历史 C++ 实现(已于 0.6.0 移除)的"失败路径"语义,并消除其穿越 FFI
//! 的 UB:
//!
//! - 空 read:历史 C++ `PoaDraftGen` 会 segfault(`threadFirstRead` 的 UB)
//! - 非法 version:历史 C++ `SparsePoa(int version)` 会 throw 穿越 extern "C"(UB)
//! - 非法碱基:历史 C++ `ReverseComplement` 会 throw `invalid_argument`(UB)
//!
//! 以上统一为返回 `Err`,由入口层在进入计算前完成校验。

use thiserror::Error;

/// POA 计算错误。公共 API 返回 `Result<Option<(String, usize)>, XpoaError>`:
/// `Ok(None)` 表示空 read 集合(历史语义),`Err` 表示输入或配置错误。
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum XpoaError {
    /// 空 read 序列。历史 C++ 后端对此会 segfault,入口统一拒绝。
    #[error("empty read sequence at index {index}")]
    EmptyRead { index: usize },

    /// read 含 NUL 字节(历史 C++ FFI 无法传参;校验保留以保证输入语义不变)。
    #[error("read sequence at index {index} contains a NUL byte")]
    NulByte { index: usize },

    /// 反向互补时遇到无法处理的碱基(历史 C++ `Sequence.cpp` 对应 throw)。
    #[error("invalid base 0x{base:02X} cannot be reverse-complemented")]
    InvalidBase { base: u8 },

    /// POA version 仅支持 1(V1 RangeFinder)与 2(V2)。
    #[error("POA version must be 1 or 2, got {version}")]
    InvalidVersion { version: i32 },
}

/// 便捷别名
pub type XpoaResult<T> = Result<T, XpoaError>;
