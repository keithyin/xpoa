//! Golden 一致性测试:读取 `tests/golden/*.json`(期望值由 C++ 后端生成,
//! 见 `examples/golden_gen.rs`),断言当前默认后端产出完全一致的
//! `(consensus, n_passes)`。
//!
//! 依赖 serde_json/serde,仅测试期使用。

use std::fs;
use std::path::Path;

use xpoa::casegen::TestCase;
use xpoa::{BACKEND_RUST, poa_consensus, poa_consensus_fwd, poa_consensus_v2};

fn run_case(case: &TestCase) -> Result<Option<(String, usize)>, xpoa::error::XpoaError> {
    // golden 期望值由 C++ 后端生成;此处用 **Rust 后端** 重现,即
    // "结果一致性"的验收断言
    let setting = golden_gen_setting(case);
    match case.mode.as_str() {
        "gen" => poa_consensus(&case.reads, &setting),
        "fwd" => poa_consensus_fwd(&case.reads, &setting),
        "v2" => poa_consensus_v2(&case.reads, &setting),
        other => panic!("unknown mode: {other}"),
    }
}

/// 与 examples/golden_gen.rs 中的 to_poa_setting 保持一致,但强制 Rust 后端
fn golden_gen_setting(case: &TestCase) -> xpoa::PoaSetting {
    let s = &case.setting;
    xpoa::PoaSetting {
        min_identity: s.min_identity,
        match_score: s.match_score,
        mismatch_score: s.mismatch_score,
        insertion_score: s.insertion_score,
        deletion_score: s.deletion_score,
        ed_unify_strand: s.ed_unify_strand,
        version: s.version,
        backend: BACKEND_RUST,
    }
}

fn collect_golden_files() -> Vec<std::path::PathBuf> {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
    let mut files: Vec<_> = fs::read_dir(&dir)
        .expect("tests/golden 目录不存在 —— 先运行 `cargo run --release --example golden_gen`")
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|e| e == "json"))
        .collect();
    files.sort();
    assert!(!files.is_empty(), "tests/golden 中没有任何 golden 文件");
    files
}

// golden 对拍为全量重负载(161 case,数千 bp×数十 reads),debug 模式慢一个
// 数量级以上 —— 请用 `cargo test --release` 运行
#[test]
#[cfg_attr(debug_assertions, ignore = "heavy: run with --release")]
fn rust_backend_matches_golden() {
    let files = collect_golden_files();
    let mut checked = 0usize;
    for path in &files {
        let raw = fs::read_to_string(path).expect("read golden");
        let value: serde_json::Value = serde_json::from_str(&raw).expect("parse golden");
        let case: TestCase =
            serde_json::from_value(value["case"].clone()).expect("parse case from golden");
        let expected = &value["result"];

        let actual = run_case(&case).expect("rust backend error");

        let actual_json = actual
            .as_ref()
            .map(|(seq, n)| serde_json::json!({"seq": seq, "n_passes": n}));

        if let Some(actual_json) = &actual_json {
            assert_eq!(
                actual_json["seq"],
                expected["seq"],
                "consensus mismatch in {} ({})",
                path.display(),
                case.name
            );
            assert_eq!(
                actual_json["n_passes"],
                expected["n_passes"],
                "n_passes mismatch in {} ({})",
                path.display(),
                case.name
            );
        } else {
            panic!(
                "backend returned None for {} ({}), expected {expected}",
                path.display(),
                case.name
            );
        }
        checked += 1;
    }
    println!("checked {checked} golden cases");
}
