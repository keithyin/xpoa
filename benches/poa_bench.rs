//! cpp / rust 双后端性能对比基准:`cargo bench`
//!
//! 场景(近似 PacBio CCS reads):
//! - A:30 reads × 1.5kb × 12% 错误(单碱基流共识的典型负载)
//! - B:100 reads × 500bp × 10% 错误(高覆盖短读)
//! - C:50 reads × 60bp × 8% 错误(短读,RangeFinder 锚点稀疏)
//!
//! 每个场景分别在 version 1/2 下跑 poa_consensus(cpp 后端与 rust 后端,
//! 同一二进制运行时切换)。

use std::time::Duration;

use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};

use xpoa::casegen::TestRead;
use xpoa::{BACKEND_CPP, BACKEND_RUST, PoaSetting};

/// 与 src/casegen::DetRng 相同的确定性 RNG(基准输入可复现)
struct DetRng {
    state: u64,
}

impl DetRng {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn gen_range(&mut self, range: std::ops::Range<usize>) -> usize {
        range.start + (self.next_u64() as usize % (range.end - range.start))
    }
    fn gen_bool(&mut self, p: f64) -> bool {
        (self.next_u64() as f64 / u64::MAX as f64) < p
    }
}

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

fn make_reads(seed: u64, n_reads: usize, ref_len: usize, err_rate: f64) -> Vec<TestRead> {
    let mut rng = DetRng::new(seed);
    let reference: Vec<u8> = (0..ref_len).map(|_| BASES[rng.gen_range(0..4)]).collect();
    let mut reads = Vec::with_capacity(n_reads);
    for i in 0..n_reads {
        let mut seq = Vec::with_capacity(ref_len + 16);
        for &base in &reference {
            if rng.gen_bool(err_rate) {
                match rng.gen_range(0..3) {
                    0 => seq.push(BASES[rng.gen_range(0..4)]),
                    1 => {
                        seq.push(BASES[rng.gen_range(0..4)]);
                        seq.push(base);
                    }
                    _ => {}
                }
            } else {
                seq.push(base);
            }
        }
        reads.push(TestRead {
            seq: String::from_utf8(seq).unwrap(),
            flags: 3,
            is_fwd: Some(true),
        });
        let _ = i;
    }
    reads
}

fn bench_scenario(c: &mut Criterion, name: &str, reads: &[TestRead]) {
    let mut group = c.benchmark_group(name);
    group.measurement_time(Duration::from_secs(10));

    for &version in [1i32, 2i32].iter() {
        for (tag, backend) in [("cpp", BACKEND_CPP), ("rust", BACKEND_RUST)] {
            let setting = PoaSetting {
                min_identity: 0.82,
                match_score: 3,
                mismatch_score: -5,
                insertion_score: -2,
                deletion_score: -2,
                ed_unify_strand: 1,
                version,
                backend,
            };
            group.bench_with_input(
                BenchmarkId::new(format!("v{version}"), tag),
                reads,
                |b, reads| {
                    b.iter(|| {
                        let result = xpoa::poa_consensus(black_box(reads), black_box(&setting))
                            .expect("benchmark inputs are valid");
                        debug_assert!(result.is_some());
                    })
                },
            );
        }
    }
    group.finish();
}

fn bench_poa(c: &mut Criterion) {
    bench_scenario(c, "A_30x1500bp_e12", &make_reads(1, 30, 1500, 0.12));
    bench_scenario(c, "B_100x500bp_e10", &make_reads(2, 100, 500, 0.10));
    bench_scenario(c, "C_50x60bp_e8", &make_reads(3, 50, 60, 0.08));
}

criterion_group!(benches, bench_poa);
criterion_main!(benches);
