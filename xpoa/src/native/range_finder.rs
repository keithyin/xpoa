//! SDP 可对齐区间定位(对应 cpp_code/RangeFinder.h/.cpp、RangeFinderV1、
//! RangeFinderV2)。
//!
//! trait 是图与定位算法之间唯一的抽象缝(对应 C++ 的 SdpRangeFinder 抽象类):
//! `poa_graph` 只依赖该 trait,V1/V2 通过 `Box<dyn SdpRangeFinder>` 注入,
//! 新增 V3 无需改动图代码。
//!
//! V2 与 V1 的唯一差异:锚点处沿共识路径向前再扩展 KMER-1 个顶点的
//! "直接区间"。C++ 原版在 `cssPos + i` 超出共识路径时越界读 vector(UB),
//! 此处按预期语义**跳过越界下标** —— 这是与 C++ 的唯一有意差异。

use crate::native::poa_graph::{INT_MAX_HALF, PoaGraph, Vertex};
use crate::native::sparse_align::sparse_align;

/// 区间哨兵:对应 `emptyInterval = (INT_MAX/2, -INT_MAX/2)`
pub const EMPTY_INTERVAL: Interval = Interval(INT_MAX_HALF, -INT_MAX_HALF);
const WIDTH: i32 = 30;
const KMER: usize = 6;

/// 闭开风格区间(对应 `std::pair<int, int>`)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Interval(pub i32, pub i32);

#[inline]
fn range_union2(a: Interval, b: Interval) -> Interval {
    Interval(a.0.min(b.0), a.1.max(b.1))
}

#[inline]
fn range_union(ranges: &[Interval]) -> Interval {
    let mut result = EMPTY_INTERVAL;
    for &r in ranges {
        result = range_union2(result, r);
    }
    result
}

#[inline]
fn next_interval(v: Interval, upper_bound: i32) -> Interval {
    if v == EMPTY_INTERVAL {
        v
    } else {
        Interval(
            v.0.saturating_add(1).min(upper_bound),
            v.1.saturating_add(1).min(upper_bound),
        )
    }
}

#[inline]
fn prev_interval(v: Interval, lower_bound: i32) -> Interval {
    if v == EMPTY_INTERVAL {
        v
    } else {
        Interval(
            v.0.saturating_sub(1).max(lower_bound),
            v.1.saturating_sub(1).max(lower_bound),
        )
    }
}

/// 对应 `binarySearchAnchors`:按 css 位置精确匹配(lower_bound)
fn binary_search_anchors(
    anchors: &[(usize, usize)],
    css_position: usize,
) -> Option<(usize, usize)> {
    let idx = anchors.partition_point(|a| a.0 < css_position);
    if idx < anchors.len() && anchors[idx].0 == css_position {
        Some(anchors[idx])
    } else {
        None
    }
}

/// V1/V2 共享的区间计算骨架(对应基类 `InitRangeFinder` 的主体)。
/// `extend_kmer`:V2 的锚点向前扩展。
fn compute_intervals(
    graph: &PoaGraph,
    consensus_path: &[Vertex],
    read_sequence: &[u8],
    anchors: &[(usize, usize)],
    extend_kmer: bool,
) -> Vec<Interval> {
    let read_length = read_sequence.len() as i32;
    let n = graph.num_vertices();

    let mut direct_ranges: Vec<Option<Interval>> = vec![None; n];
    let mut fwd_marks: Vec<Interval> = vec![EMPTY_INTERVAL; n];
    let mut rev_marks: Vec<Interval> = vec![EMPTY_INTERVAL; n];

    // 锚点 → 直接区间
    for (css_pos, &v) in consensus_path.iter().enumerate() {
        match binary_search_anchors(anchors, css_pos) {
            Some((_, read_pos)) => {
                let rp = read_pos as i32;
                direct_ranges[v] =
                    Some(Interval((rp - WIDTH).max(0), (rp + WIDTH).min(read_length)));
                if extend_kmer {
                    for i in 1..KMER {
                        let idx = css_pos + i;
                        if idx >= consensus_path.len() {
                            // C++ 原版此处越界读 consensusPath[cssPos+i](UB);
                            // 按预期语义跳过
                            continue;
                        }
                        let v2 = consensus_path[idx];
                        let rp = (read_pos + i) as i32;
                        direct_ranges[v2] =
                            Some(Interval((rp - WIDTH).max(0), (rp + WIDTH).min(read_length)));
                    }
                }
            }
            None => {
                direct_ranges[v] = None;
            }
        }
    }

    let sorted = graph.topological_order();

    // 前向递推:无直接区间的顶点取前驱"步进"区间的并
    for &v in &sorted {
        if let Some(dr) = direct_ranges[v] {
            fwd_marks[v] = dr;
        } else {
            let stepped: Vec<Interval> = graph
                .sorted_in_edges(v)
                .iter()
                .map(|&pred| next_interval(fwd_marks[pred], read_length))
                .collect();
            fwd_marks[v] = range_union(&stepped);
        }
    }

    // 后向递推(顺序无关的并集,故 out_edges 用插入序)
    for &v in sorted.iter().rev() {
        if let Some(dr) = direct_ranges[v] {
            rev_marks[v] = dr;
        } else {
            let stepped: Vec<Interval> = graph
                .out_edges(v)
                .iter()
                .map(|&succ| prev_interval(rev_marks[succ], 0))
                .collect();
            rev_marks[v] = range_union(&stepped);
        }
    }

    // 前向/后向取包络(按顶点 id 索引,与 find_alignable_range 的
    // `intervals[v]` 查询一致;不能按拓扑序 collect —— 那会错位!)
    let mut out = vec![EMPTY_INTERVAL; n];
    for &v in &sorted {
        out[v] = range_union2(fwd_marks[v], rev_marks[v]);
    }
    out
}

/// 对应 C++ `SdpRangeFinder` 抽象类。`poa_graph::try_add_read` 仅依赖此
/// 接口的两个方法。
pub trait SdpRangeFinder {
    /// 每次对新 read 对齐前调用(状态随下一次 init 重置)
    fn init(
        &mut self,
        graph: &PoaGraph,
        consensus_path: &[Vertex],
        consensus_sequence: &[u8],
        read_sequence: &[u8],
    );

    /// 查询顶点 v 的可对齐 read 区间 [start, end);end == -INT_MAX/2 表示空
    fn find_alignable_range(&self, v: Vertex) -> (i32, i32);
}

/// 对应 `SdpRangeFinderV1`
pub struct SdpRangeFinderV1 {
    intervals: Vec<Interval>,
}

impl SdpRangeFinderV1 {
    pub fn new() -> Self {
        Self {
            intervals: Vec::new(),
        }
    }
}

impl Default for SdpRangeFinderV1 {
    fn default() -> Self {
        Self::new()
    }
}

impl SdpRangeFinder for SdpRangeFinderV1 {
    fn init(
        &mut self,
        graph: &PoaGraph,
        consensus_path: &[Vertex],
        consensus_sequence: &[u8],
        read_sequence: &[u8],
    ) {
        let anchors = sparse_align(KMER, consensus_sequence, read_sequence);
        self.intervals = compute_intervals(graph, consensus_path, read_sequence, &anchors, false);
    }

    fn find_alignable_range(&self, v: Vertex) -> (i32, i32) {
        let iv = self.intervals[v];
        (iv.0, iv.1)
    }
}

/// 对应 `SdpRangeFinderV2`(锚点向前扩展 KMER-1 个顶点)
pub struct SdpRangeFinderV2 {
    intervals: Vec<Interval>,
}

impl SdpRangeFinderV2 {
    pub fn new() -> Self {
        Self {
            intervals: Vec::new(),
        }
    }
}

impl Default for SdpRangeFinderV2 {
    fn default() -> Self {
        Self::new()
    }
}

impl SdpRangeFinder for SdpRangeFinderV2 {
    fn init(
        &mut self,
        graph: &PoaGraph,
        consensus_path: &[Vertex],
        consensus_sequence: &[u8],
        read_sequence: &[u8],
    ) {
        let anchors = sparse_align(KMER, consensus_sequence, read_sequence);
        self.intervals = compute_intervals(graph, consensus_path, read_sequence, &anchors, true);
    }

    fn find_alignable_range(&self, v: Vertex) -> (i32, i32) {
        let iv = self.intervals[v];
        (iv.0, iv.1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native::poa_graph::{AlignConfig, AlignMode};

    fn make_chain_graph(seq: &[u8]) -> PoaGraph {
        let mut g = PoaGraph::new();
        g.add_first_read(seq, None);
        g
    }

    #[test]
    fn v1_covers_chain_read() {
        let ref_seq = b"ACGTACGTACGT";
        let g = make_chain_graph(ref_seq);
        let read = b"ACGTACGTACGT";
        let mut rf = SdpRangeFinderV1::new();
        // TryAddRead 中 C++ 的 consensusPath 默认 minCoverage = -INT_MAX
        let css_path = g.consensus_path(AlignMode::Local, -INT_MAX_HALF);
        let css_seq = g.sequence_along_path(&css_path);
        rf.init(&g, &css_path, &css_seq, read);
        // 每个 read 碱基的锚点使中间顶点的区间覆盖整个 read
        for v in 2..g.num_vertices() - 1 {
            let (s, e) = rf.find_alignable_range(v);
            assert!(s <= e, "vertex {v}: interval ({s},{e})");
        }
    }

    #[test]
    fn enter_exit_have_empty_intervals() {
        // 短 read(< KMER)无锚点:所有顶点的区间为 EMPTY
        // enter 无前驱 → fwd EMPTY;exit 无后继 → rev EMPTY;中间全 EMPTY
        let g = make_chain_graph(b"ACGTACGT");
        let mut rf = SdpRangeFinderV1::new();
        let css_path = g.consensus_path(AlignMode::Local, -INT_MAX_HALF);
        let css_seq = g.sequence_along_path(&css_path);
        rf.init(&g, &css_path, &css_seq, b"AC"); // < KMER,无锚
        for v in 0..g.num_vertices() {
            let (s, e) = rf.find_alignable_range(v);
            assert_eq!((s, e), (INT_MAX_HALF, -INT_MAX_HALF), "vertex {v}");
        }
    }

    #[test]
    fn anchors_fill_enter_exit_hulls() {
        // 完全匹配 read:enter 的 hull 由后继 rev 递推填充(非 EMPTY),
        // exit 的 hull 由前驱 fwd 递推填充 —— 与 C++ 一致
        let g = make_chain_graph(b"ACGTACGT");
        let mut rf = SdpRangeFinderV1::new();
        let css_path = g.consensus_path(AlignMode::Local, -INT_MAX_HALF);
        let css_seq = g.sequence_along_path(&css_path);
        rf.init(&g, &css_path, &css_seq, b"ACGTACGT");
        for v in 0..g.num_vertices() {
            let (s, e) = rf.find_alignable_range(v);
            assert!(s <= e, "vertex {v}: ({s},{e})");
        }
    }

    #[test]
    fn config_local_mode_wiring() {
        // 保证 poa_graph 与 range_finder 的接口能编译协同
        let mut g = make_chain_graph(b"ACGTACGTACGTACGT");
        let config = AlignConfig {
            match_score: 3,
            mismatch_score: -5,
            insert: -2,
            delete: -2,
            mode: AlignMode::Local,
        };
        let mut rf = SdpRangeFinderV2::new();
        let mat = g.try_add_read(b"ACGTACGTACGTACGT", &config, &mut rf);
        assert!(mat.score > 0.0);
        g.commit_add(mat, None);
        assert_eq!(g.num_reads(), 2);
    }
}
