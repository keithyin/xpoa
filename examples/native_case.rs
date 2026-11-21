//! 参数化 case 重放(native 侧):与 tools/dump_case.cpp 输出对齐。
//! 用法:native_case <tsv> <n_reads> <mode:0=orient|1=unify|2=fwd>

use std::fs;

use xpoa::native::edlib::edlib_global_align_distance;
use xpoa::native::poa_graph::{AlignConfig, AlignMode, PoaGraph};
use xpoa::native::range_finder::SdpRangeFinderV1;
use xpoa::native::sequence::reverse_complement;

fn dump_graph(g: &PoaGraph, i: usize) {
    let mut s = String::from("Id,Base,Reads,SpanningReads\n");
    for id in 0..g.num_vertices() {
        let n = g.node(id);
        s.push_str(&format!(
            "{},{},{},{}\n",
            n.id, n.base as char, n.reads, n.spanning_reads
        ));
    }
    fs::write(format!("rustN{i}.csv"), s).unwrap();
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let tsv = &args[1];
    let n: usize = args[2].parse().unwrap();
    let mode: i32 = args[3].parse().unwrap();

    let raw = fs::read_to_string(tsv).unwrap();
    let mut lines = raw.lines();
    let first: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let _version: i32 = first[0].parse().unwrap();
    let _unify: i32 = first[1].parse().unwrap();
    let mut reads: Vec<(String, i32)> = Vec::new();
    let mut fwd_flags: Vec<Option<bool>> = Vec::new();
    for line in lines {
        let mut parts = line.split('\t');
        let seq = parts.next().unwrap();
        let flags: i32 = parts.next().unwrap().parse().unwrap();
        let fwd = parts
            .next()
            .and_then(|s| if s.is_empty() { None } else { Some(s == "1") });
        reads.push((seq.to_string(), flags));
        fwd_flags.push(fwd);
    }
    let reads: Vec<(String, i32)> = reads.into_iter().take(n).collect();
    let fwd_flags: Vec<Option<bool>> = fwd_flags.into_iter().take(n).collect();

    // mode=3:直接走公共入口 poa_consensus_v2(含 v2 预处理与后端分发)
    if mode == 3 {
        use xpoa::casegen::TestRead;

        let test_reads: Vec<TestRead> = reads
            .iter()
            .zip(&fwd_flags)
            .map(|(r, f)| TestRead {
                seq: r.0.clone(),
                flags: r.1,
                is_fwd: *f,
            })
            .collect();
        let setting = xpoa::PoaSetting {
            min_identity: 0.82,
            match_score: 3,
            mismatch_score: -5,
            insertion_score: -2,
            deletion_score: -2,
            ed_unify_strand: 1,
            version: 1,
            backend: xpoa::BACKEND_RUST,
        };
        let result = xpoa::poa_consensus_v2(&test_reads, &setting).expect("native backend error");
        if let Some((seq, np)) = result {
            println!("consensus: {seq}");
            println!("n_passes: {np}");
        } else {
            println!("None");
        }
        return;
    }

    let config = AlignConfig {
        match_score: 3,
        mismatch_score: -5,
        insert: -2,
        delete: -2,
        mode: AlignMode::Local,
    };
    let mut g = PoaGraph::new();
    let mut rf = SdpRangeFinderV1::new();
    let mut path = Vec::new();

    for i in 0..reads.len() {
        let (seq, _) = (&reads[i].0, reads[i].1);
        if i > 0 {
            let css_path =
                g.consensus_path(AlignMode::Local, -xpoa::native::poa_graph::INT_MAX_HALF);
            let css_seq = g.sequence_along_path(&css_path);
            println!("midcss{i}: {}", String::from_utf8_lossy(&css_seq));
        }
        if i == 0 {
            g.add_first_read(seq.as_bytes(), Some(&mut path));
        } else {
            match mode {
                2 => {
                    let m = g.try_add_read(seq.as_bytes(), &config, &mut rf);
                    println!("try{i} score={}", m.score);
                    if m.score >= 0.0 {
                        g.commit_add(m, Some(&mut path));
                    }
                }
                1 => {
                    let rc = reverse_complement(seq.as_bytes()).unwrap();
                    let d1 = edlib_global_align_distance(seq.as_bytes(), reads[0].0.as_bytes());
                    let d2 = edlib_global_align_distance(&rc, reads[0].0.as_bytes());
                    let chosen: &[u8] = if d1 > d2 { &rc } else { seq.as_bytes() };
                    let m = g.try_add_read(chosen, &config, &mut rf);
                    g.commit_add(m, Some(&mut path));
                }
                _ => {
                    let rc = reverse_complement(seq.as_bytes()).unwrap();
                    let m1 = g.try_add_read(seq.as_bytes(), &config, &mut rf);
                    let m2 = g.try_add_read(&rc, &config, &mut rf);
                    if m1.score >= m2.score && m1.score >= 0.0 {
                        g.commit_add(m1, Some(&mut path));
                    } else if m2.score >= m1.score && m2.score >= 0.0 {
                        g.commit_add(m2, Some(&mut path));
                    }
                }
            }
        }
        dump_graph(&g, i);
    }

    if !reads.is_empty() {
        let min_cov: i32 = args
            .get(4)
            .map(|s| s.parse().unwrap())
            .unwrap_or(-xpoa::native::poa_graph::INT_MAX_HALF);
        let best = g.consensus_path(AlignMode::Local, min_cov);
        let css = g.sequence_along_path(&best);
        println!(
            "consensus(minCov={min_cov}): {}",
            String::from_utf8_lossy(&css)
        );
        println!(
            "path: {}",
            best.iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(" ")
        );
    }
}
