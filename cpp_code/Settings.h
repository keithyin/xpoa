#pragma once

extern "C"
{
    // backend 字段由 Rust 端用于运行时选择实现(0=cpp, 1=rust);C++ 端忽略其值。
    // 放在末尾以保持 ABI 兼容(见 src/setting.rs)。
    struct PoaSetting
    {
        float min_identity = 0.82;
        int match_score = 3;
        int mismatch_score = -5;
        int insertion_score = -2;
        int deletion_score = -2;
        int ed_unify_strand = 1;
        int version = 1;
        int backend = 0;
    };
}