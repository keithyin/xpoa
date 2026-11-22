//! 调试工具:用 native 后端重建 case_004 的流程,输出每步图节点表,
//! 与 tools/cpp_dump.cpp 的输出对拍。

use std::fs;

use xpoa::native::poa_graph::PoaGraph;
use xpoa::native::range_finder::SdpRangeFinderV1;
use xpoa::native::sequence::reverse_complement;
use xpoa::native::sparse_poa::{ReadKey, SparsePoa};
use xpoa::{BACKEND_RUST, PoaSetting};

const ADAPTER_BEFORE: i32 = 1;
const ADAPTER_AFTER: i32 = 2;

fn dump_graph(g: &PoaGraph, path: &str) {
    let mut s = String::from("Id,Base,Reads,SpanningReads,Score,ReachingScore\n");
    for &v in &g.topological_order() {
        let n = g.node(v);
        s.push_str(&format!(
            "{},{},{},{},{},{}\n",
            n.id, n.base as char, n.reads, n.spanning_reads, n.score, n.reaching_score
        ));
    }
    fs::write(path, s).unwrap();
}

fn dump_edges(g: &xpoa::native::poa_graph::PoaGraph, header: &str) {
    println!("{header}");
    for (u, v) in g.edges() {
        println!("{u} -> {v}");
    }
}

fn move_code(m: xpoa::native::poa_graph::MoveType) -> i32 {
    use xpoa::native::poa_graph::MoveType::*;
    match m {
        Invalid => 0,
        Start => 1,
        End => 2,
        Match => 3,
        Mismatch => 4,
        Delete => 5,
        Extra => 6,
    }
}

fn dump_cols(
    g: &xpoa::native::poa_graph::PoaGraph,
    mat: &xpoa::native::poa_graph::PoaAlignmentMatrix,
) {
    for id in 0..g.num_vertices() {
        if let Some(col) = mat.columns.get_opt(id) {
            for row in col.begin_row()..col.end_row() {
                println!(
                    "{} {} {:.4} {} {}",
                    id,
                    row,
                    col.score.get(row),
                    move_code(col.reaching_move.get(row)),
                    col.previous_vertex.get(row)
                );
            }
        }
    }
}

fn main() {
    let reads = [
        "CCCAAGAGGGACTTCAGACGGACGCGGGGCGGGTAGTTGTATCTAACCAGG",
        "CCCAAGAGGGACTTCAGAACTGACGCGGGGCGGGTAGTTGTATCTAACCACGG",
        "CCCAAGAGGGACTTCAAACGACGCGGGGCGTGGTAGGTTGTATCTAACCAGG",
    ];
    let flags = [3i32, 3, 3];

    let setting = PoaSetting {
        min_identity: 0.82,
        match_score: 3,
        mismatch_score: -5,
        insertion_score: -2,
        deletion_score: -2,
        ed_unify_strand: 0,
        version: 1,
        backend: BACKEND_RUST,
    };

    let mut poa = SparsePoa::new(setting.version).unwrap();
    let mut keys: Vec<ReadKey> = Vec::new();
    let mut cov = 0usize;

    for (i, seq) in reads.iter().enumerate() {
        let seq = seq.as_bytes();
        // try_add_read 的中间共识(与 C++ TryAddRead 内部一致)
        if poa.num_reads() > 0 {
            let css_path = poa.graph().consensus_path(
                xpoa::native::poa_graph::AlignMode::Local,
                -xpoa::native::poa_graph::INT_MAX_HALF,
            );
            let css_seq = poa.graph().sequence_along_path(&css_path);
            println!(
                "intermediate_css for read{i}: {}",
                String::from_utf8_lossy(&css_seq)
            );
            let anchors = xpoa::native::sparse_align::sparse_align(6, &css_seq, seq);
            println!("anchors for read{i}: {:?}", anchors);
        }
        let key = if setting.ed_unify_strand != 0 {
            if i == 0 {
                poa.add_read(seq, &setting, true)
            } else {
                let rc = reverse_complement(seq).unwrap();
                let d1 = xpoa::native::edlib::edlib_global_align_distance(seq, reads[0].as_bytes());
                let d2 = xpoa::native::edlib::edlib_global_align_distance(&rc, reads[0].as_bytes());
                if d1 > d2 {
                    poa.add_read(&rc, &setting, true)
                } else {
                    poa.add_read(seq, &setting, false)
                }
            }
        } else {
            poa.orient_and_add_read(seq, &setting).unwrap()
        };
        keys.push(key);
        if key >= 0 {
            cov += 1;
        }
        println!("read{i} key={key}");
        dump_graph(poa.graph(), &format!("rust_graph_after_{i}.csv"));
        dump_edges(poa.graph(), &format!("=== after read{i} ==="));
    }

    let min_cov: i32 = if cov < 5 {
        1
    } else {
        (cov.div_ceil(2) - 1) as i32
    };
    let (css, summaries) = poa.find_consensus(min_cov, &setting);
    dump_graph(poa.graph(), "rust_graph_final.csv");

    // 与 tools 侧 dump_cols.cpp 对拍:手动重放并 dump read2 的 DP 列
    {
        let mut g = PoaGraph::new();
        let mut path = Vec::new();
        g.add_first_read(reads[0].as_bytes(), Some(&mut path));
        let config = xpoa::native::poa_graph::AlignConfig {
            match_score: 3,
            mismatch_score: -5,
            insert: -2,
            delete: -2,
            mode: xpoa::native::poa_graph::AlignMode::Local,
        };
        let mut rf = SdpRangeFinderV1::new();
        let m1 = g.try_add_read(reads[1].as_bytes(), &config, &mut rf);
        g.commit_add(m1, Some(&mut path));
        let m2 = g.try_add_read(reads[2].as_bytes(), &config, &mut rf);
        dump_cols(&g, &m2);
    }

    println!("consensus: {}", String::from_utf8(css.clone()).unwrap());
    for (i, s) in summaries.iter().enumerate() {
        println!(
            "read{i} identity={:.6} rc={} read=[{},{}] css=[{},{}]",
            s.alignment_identity,
            s.reverse_complemented_read as i32,
            s.extent_on_read.0,
            s.extent_on_read.1,
            s.extent_on_consensus.0,
            s.extent_on_consensus.1
        );
    }
    let mut n_passes = 0;
    for i in 0..reads.len() {
        if keys[i] < 0 {
            continue;
        }
        if summaries[keys[i] as usize].alignment_identity < setting.min_identity {
            continue;
        }
        if !(flags[i] & ADAPTER_BEFORE != 0 && flags[i] & ADAPTER_AFTER != 0) {
            continue;
        }
        n_passes += 1;
    }
    println!("n_passes: {n_passes}");
    let _ = fs::remove_file("unused.csv");
}
