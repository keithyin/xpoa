//! 生成 golden 回归基线:`cargo run --release --example golden_gen`
//!
//! 对 `src/casegen.rs` 产出的全部用例跑当前实现,把
//! `(consensus, n_passes)` 固化到 `tests/golden/case_NNN.json`。
//!
//! 历史说明:这些快照最初在 v0.5.0 由 C++ 后端生成,当时 Rust 后端经
//! golden + 属性式对拍验证与其逐字符一致;0.6.0 移除 C++ 后端后,本工具
//! 改用当前实现重生成,输出与既有基线逐字节相同(CI 有重生成校验)。

use std::fs;
use std::path::Path;

use xpoa::casegen::{TestCase, TestSetting, generate_cases};
use xpoa::{PoaSetting, poa_consensus, poa_consensus_fwd, poa_consensus_v2};

pub fn to_poa_setting(s: &TestSetting) -> PoaSetting {
    PoaSetting {
        min_identity: s.min_identity,
        match_score: s.match_score,
        mismatch_score: s.mismatch_score,
        insertion_score: s.insertion_score,
        deletion_score: s.deletion_score,
        ed_unify_strand: s.ed_unify_strand,
        version: s.version,
    }
}

fn run_case(case: &TestCase) -> Result<Option<(String, usize)>, xpoa::error::XpoaError> {
    let setting = to_poa_setting(&case.setting);
    match case.mode.as_str() {
        "gen" => poa_consensus(&case.reads, &setting),
        "fwd" => poa_consensus_fwd(&case.reads, &setting),
        "v2" => poa_consensus_v2(&case.reads, &setting),
        other => panic!("unknown mode: {other}"),
    }
}

fn main() {
    let cases = generate_cases();
    // 以 package 目录为基准,避免从 workspace 根执行 cargo run 时写错位置
    let out_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
    fs::create_dir_all(&out_dir).expect("create tests/golden");

    // 清理旧 golden,保证集合与生成器同步
    if let Ok(entries) = fs::read_dir(&out_dir) {
        for entry in entries.flatten() {
            if entry.path().extension().is_some_and(|e| e == "json") {
                fs::remove_file(entry.path()).expect("remove stale golden");
            }
        }
    }

    let mut failures = Vec::new();
    for (i, case) in cases.iter().enumerate() {
        let result = run_case(case).expect("backend error");
        let json = serde_json::json!({
            "case": case,
            "result": result.as_ref().map(|(seq, n)| {
                serde_json::json!({"seq": seq, "n_passes": n})
            }),
        });
        let path = out_dir.join(format!("case_{i:03}.json"));
        fs::write(&path, serde_json::to_string_pretty(&json).unwrap()).expect("write golden file");

        // 运行时 sanity:C++ 后端必须返回有效共识(除非全部 read 被拒)
        match &result {
            None => failures.push(format!("case_{i:03} ({}) returned None", case.name)),
            Some((seq, _)) if seq.is_empty() => failures.push(format!(
                "case_{i:03} ({}) produced empty consensus",
                case.name
            )),
            _ => {}
        }
    }

    println!(
        "generated {} golden files into {}",
        cases.len(),
        out_dir.display()
    );
    if !failures.is_empty() {
        for f in &failures {
            println!("WARN: {f}");
        }
    }
}
