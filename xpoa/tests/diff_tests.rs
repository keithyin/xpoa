//! 属性式双后端差分测试:随机生成用例,同一输入实时对拍
//! C++ 后端与 Rust 后端的 `(consensus, n_passes)`。
//!
//! 与 golden_tests(固化快照)互补:这里覆盖无限输入空间。
//! 需要 `--features cpp`(默认开启);纯 Rust 构建下跳过。

#![cfg(feature = "cpp")]

use xpoa::casegen::TestRead;
use xpoa::{BACKEND_CPP, BACKEND_RUST, PoaSetting};

type ModeResults = Vec<Result<Option<(String, usize)>, xpoa::error::XpoaError>>;

fn run_all_modes(reads: &[TestRead], setting: &PoaSetting) -> ModeResults {
    let mut out: ModeResults = Vec::new();
    for backend in [BACKEND_CPP, BACKEND_RUST] {
        let mut s = *setting;
        s.backend = backend;
        out.push(xpoa::poa_consensus(reads, &s));
        out.push(xpoa::poa_consensus_fwd(reads, &s));
        out.push(xpoa::poa_consensus_v2(reads, &s));
    }
    out
}

fn assert_backends_agree(reads: &[TestRead], setting: &PoaSetting, tag: &str) {
    let cpp = &run_all_modes(reads, setting)[..3];
    let native = &run_all_modes(reads, setting)[3..];
    for (mode, (c, n)) in ["gen", "fwd", "v2"].iter().zip(cpp.iter().zip(native)) {
        // Result 整体相等:Ok 值一致,或 Err 错误一致(错误路径同样差分)
        assert_eq!(
            c, n,
            "{tag} mode={mode}: cpp 与 rust 后端结果不一致\n  cpp={c:?}\n  rust={n:?}"
        );
    }
}

/// 确定性随机 read 集生成(与 src/casegen 的 DetRng 一致的思路)
fn make_reads(
    rng: &mut xpoa::casegen::DetRng,
    n_reads: usize,
    ref_len: usize,
    err_rate: f64,
    mixed_dir: bool,
) -> (Vec<TestRead>, Vec<String>) {
    use xpoa::casegen::BASES;

    let mut reference = Vec::with_capacity(ref_len);
    for _ in 0..ref_len {
        reference.push(BASES[rng.gen_range(0..4)]);
    }

    let mutate = |rng: &mut xpoa::casegen::DetRng, seq: &[u8]| -> Vec<u8> {
        let mut out = Vec::with_capacity(seq.len());
        for &base in seq {
            if rng.gen_bool(err_rate) {
                match rng.gen_range(0..3) {
                    0 => out.push(BASES[rng.gen_range(0..4)]),
                    1 => {
                        out.push(BASES[rng.gen_range(0..4)]);
                        out.push(base);
                    }
                    _ => {}
                }
            } else {
                out.push(base);
            }
        }
        out
    };

    let reads: Vec<TestRead> = (0..n_reads)
        .map(|i| {
            let mut seq = mutate(rng, &reference);
            if mixed_dir && rng.gen_bool(0.35) {
                // 模拟 v2 的反向 read:直接给 RC 后的序列 + is_fwd=false,
                // 两种后端拿到的输入完全一致(预处理都在 Rust 侧完成)
                seq = xpoa::casegen::revcomp(&seq);
            }
            let is_fwd = if mixed_dir {
                Some(i % 2 == 0)
            } else {
                Some(true)
            };
            TestRead {
                seq: String::from_utf8(seq).unwrap(),
                flags: xpoa::casegen::ADAPTER_BEFORE | xpoa::casegen::ADAPTER_AFTER,
                is_fwd,
            }
        })
        .collect();
    (reads, Vec::<String>::new())
}

#[test]
#[cfg_attr(debug_assertions, ignore = "heavy: run with --release")]
fn backends_agree_on_random_cases() {
    let mut rng = xpoa::casegen::DetRng::new(20240101);
    let mut cases = 0usize;

    for &(ref_len, n_reads) in [
        (40usize, 3usize),
        (50, 10),
        (120, 5),
        (200, 20),
        (500, 10),
        (1000, 4),
        (2000, 3),
    ]
    .iter()
    {
        for &err_rate in [0.03f64, 0.08, 0.15].iter() {
            for &(unify, version) in [(1i32, 1i32), (1, 2), (0, 1)].iter() {
                let (reads, _) = make_reads(&mut rng, n_reads, ref_len, err_rate, false);
                let setting = PoaSetting {
                    min_identity: 0.82,
                    match_score: 3,
                    mismatch_score: -5,
                    insertion_score: -2,
                    deletion_score: -2,
                    ed_unify_strand: unify,
                    version,
                    backend: BACKEND_CPP,
                };
                assert_backends_agree(
                    &reads,
                    &setting,
                    &format!("grid_l{ref_len}_n{n_reads}_e{err_rate}_u{unify}_v{version}"),
                );
                cases += 1;
            }
        }
    }

    // 混合方向(v2 模式)
    for &(ref_len, n_reads) in [(100usize, 8usize), (500, 15)].iter() {
        for &err_rate in [0.05f64, 0.12].iter() {
            for &version in [1i32, 2i32].iter() {
                let (reads, _) = make_reads(&mut rng, n_reads, ref_len, err_rate, true);
                let setting = PoaSetting {
                    min_identity: 0.82,
                    match_score: 3,
                    mismatch_score: -5,
                    insertion_score: -2,
                    deletion_score: -2,
                    ed_unify_strand: 1,
                    version,
                    backend: BACKEND_CPP,
                };
                assert_backends_agree(
                    &reads,
                    &setting,
                    &format!("mixed_l{ref_len}_n{n_reads}_e{err_rate}_v{version}"),
                );
                cases += 1;
            }
        }
    }

    // 边缘:超短 read、单条、同聚物、非 ACGT 字符
    {
        let mut setting = PoaSetting {
            min_identity: 0.82,
            match_score: 3,
            mismatch_score: -5,
            insertion_score: -2,
            deletion_score: -2,
            ed_unify_strand: 1,
            version: 2,
            backend: BACKEND_CPP,
        };

        let edge_sets: Vec<Vec<TestRead>> = vec![
            vec![
                TestRead {
                    seq: "A".into(),
                    flags: 3,
                    is_fwd: Some(true),
                },
                TestRead {
                    seq: "C".into(),
                    flags: 3,
                    is_fwd: Some(true),
                },
                TestRead {
                    seq: "A".into(),
                    flags: 3,
                    is_fwd: Some(true),
                },
            ],
            vec![
                TestRead {
                    seq: "AC".into(),
                    flags: 3,
                    is_fwd: Some(true),
                },
                TestRead {
                    seq: "CA".into(),
                    flags: 3,
                    is_fwd: Some(true),
                },
            ],
            vec![
                TestRead {
                    seq: "A".repeat(80),
                    flags: 3,
                    is_fwd: Some(true),
                },
                TestRead {
                    seq: "A".repeat(80),
                    flags: 3,
                    is_fwd: Some(true),
                },
                TestRead {
                    seq: "A".repeat(80),
                    flags: 3,
                    is_fwd: Some(true),
                },
            ],
            vec![
                TestRead {
                    seq: "ACGTNNNNACGTACGTNNACGT".into(),
                    flags: 3,
                    is_fwd: Some(true),
                },
                TestRead {
                    seq: "ACGTNNNNACGTACGTNNACGT".into(),
                    flags: 3,
                    is_fwd: Some(true),
                },
            ],
            vec![
                TestRead {
                    seq: "ACGTACGTAC".into(),
                    flags: 3,
                    is_fwd: Some(true),
                },
                TestRead {
                    seq: "ACGTACGTAC".into(),
                    flags: 0,
                    is_fwd: Some(false),
                },
                TestRead {
                    seq: "ACGTACGTAC".into(),
                    flags: 1,
                    is_fwd: Some(true),
                },
            ],
        ];
        for (i, reads) in edge_sets.iter().enumerate() {
            for &version in [1i32, 2i32].iter() {
                setting.version = version;
                assert_backends_agree(reads, &setting, &format!("edge{i}_v{version}"));
                cases += 1;
            }
        }
    }

    println!("backends agree on {cases} random cases × 3 modes × 2 backends");
}
