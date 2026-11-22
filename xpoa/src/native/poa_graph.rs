//! POA 图:结构 + 稀疏 DP 列 + 回溯穿线 + 共识路径。
//!
//! 对应 cpp_code 的 PoaGraphImpl.cpp、PoaGraphTraversals.cpp、
//! PoaAlignmentMatrix.h、VectorL.h、AlignConfig.h、BoostGraph.h。
//!
//! 与 C++ 的结构对应(一致性关键):
//! - 顶点句柄 `Vertex = usize`,与插入序一致。C++ 中 external id ==
//!   internal id == vertex_index(`PruneGraph` 在 FFI 路径上从不被调用),
//!   因此省去 internalize/externalize 层。
//! - BGL `adjacency_list<setS, listS, bidirectionalS>`:`setS` 的边存储
//!   **自动去重**(`add_edge` 可能对同一 (u, fork) 重复调用)→
//!   `HashSet<(u, v)>` 查重;`listS` 顶点迭代 = 插入序 → Vec 下标序。
//! - `inEdges()` 显式按 (source, target) 排序(EdgeComparator);
//!   out_edges 迭代序不影响结果,这里用插入序。
//! - 所有可能的 float 值均为整数值(match/mismatch/insert/delete 为整数,
//!   consensusPath 加 -0.0001f 偏移),不涉及舍入顺序敏感的运算。

use std::collections::HashSet;

use crate::native::range_finder::SdpRangeFinder;

pub const INT_MAX_HALF: i32 = i32::MAX / 2;

pub const NULL_VERTEX: usize = usize::MAX;
pub const ENTER_VERTEX_BASE: u8 = b'^';
pub const EXIT_VERTEX_BASE: u8 = b'$';

/// 对应 `AlignMode`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignMode {
    Global,
    SemiGlobal,
    Local,
}

/// 对应 `AlignParams` + `AlignConfig`
#[derive(Debug, Clone, Copy)]
pub struct AlignConfig {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub insert: i32,
    pub delete: i32,
    pub mode: AlignMode,
}

/// 对应 `MoveType`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MoveType {
    Invalid,
    Start,
    End,
    Match,
    Mismatch,
    Delete,
    Extra,
}

/// 对应 `PoaNode`
#[derive(Debug, Clone)]
pub struct PoaNode {
    pub id: Vertex,
    pub base: u8,
    pub reads: i32,
    pub spanning_reads: i32,
    pub score: f32,
    pub reaching_score: f32,
}

pub type Vertex = usize;

/// 行区间向量(对应 `VectorL<T>`):只存储 [begin_row, end_row)。
/// 所有使用场景的 T 均为 Copy(f32 / MoveType / Vertex)。
#[derive(Debug, Clone)]
pub struct VectorL<T: Copy> {
    storage: Vec<T>,
    begin_row: usize,
    end_row: usize,
}

impl<T: Copy> VectorL<T> {
    pub fn new(begin_row: i32, end_row: i32, default_val: T) -> Self {
        assert!(end_row >= begin_row, "VectorL: end < begin");
        Self {
            storage: vec![default_val; (end_row - begin_row) as usize],
            begin_row: begin_row as usize,
            end_row: end_row as usize,
        }
    }

    #[inline]
    pub fn get(&self, pos: usize) -> T {
        debug_assert!(self.begin_row <= pos && pos < self.end_row);
        self.storage[pos - self.begin_row]
    }

    #[inline]
    pub fn set(&mut self, pos: usize, value: T) {
        debug_assert!(self.begin_row <= pos && pos < self.end_row);
        self.storage[pos - self.begin_row] = value;
    }

    #[inline]
    pub fn begin_row(&self) -> usize {
        self.begin_row
    }

    #[inline]
    pub fn end_row(&self) -> usize {
        self.end_row
    }

    #[inline]
    pub fn has_row(&self, i: usize) -> bool {
        self.begin_row <= i && i < self.end_row
    }
}

/// 对应 `ArgMax(VectorL<T>)`:返回**第一个**最大值的行号
/// (C++ `std::max_element` 语义,注意 Rust `max_by` 返回最后一个)
impl VectorL<f32> {
    pub fn arg_max(&self) -> usize {
        let mut best_idx = 0usize;
        let mut best_val = f32::NEG_INFINITY;
        for (i, &v) in self.storage.iter().enumerate() {
            if v > best_val {
                best_val = v;
                best_idx = i;
            }
        }
        self.begin_row + best_idx
    }
}

/// 对应 `AlignmentColumn`
#[derive(Debug, Clone)]
pub struct AlignmentColumn {
    pub current_vertex: Vertex,
    pub score: VectorL<f32>,
    pub reaching_move: VectorL<MoveType>,
    pub previous_vertex: VectorL<Vertex>,
}

impl AlignmentColumn {
    pub fn new(vertex: Vertex, begin_row: i32, end_row: i32) -> Self {
        Self {
            current_vertex: vertex,
            score: VectorL::new(begin_row, end_row, f32::NEG_INFINITY),
            reaching_move: VectorL::new(begin_row, end_row, MoveType::Invalid),
            previous_vertex: VectorL::new(begin_row, end_row, NULL_VERTEX),
        }
    }

    #[inline]
    pub fn begin_row(&self) -> usize {
        self.score.begin_row()
    }

    #[inline]
    pub fn end_row(&self) -> usize {
        self.score.end_row()
    }

    #[inline]
    pub fn has_row(&self, i: usize) -> bool {
        self.begin_row() <= i && i < self.end_row()
    }
}

/// 对齐列集合(对应 `AlignmentColumnMap`;C++ 是 unordered_map<VD, col*>,
/// 仅按键访问,这里用顶点索引的 Vec)
pub struct AlignmentColumnMap {
    columns: Vec<Option<AlignmentColumn>>,
}

impl AlignmentColumnMap {
    pub fn new(num_vertices: usize) -> Self {
        Self {
            columns: (0..num_vertices).map(|_| None).collect(),
        }
    }

    pub fn insert(&mut self, v: Vertex, col: AlignmentColumn) {
        self.columns[v] = Some(col);
    }

    pub fn get(&self, v: Vertex) -> &AlignmentColumn {
        self.columns[v].as_ref().expect("alignment column missing")
    }

    /// 调试/验证用
    pub fn get_opt(&self, v: Vertex) -> Option<&AlignmentColumn> {
        self.columns[v].as_ref()
    }
}

/// 对应 `PoaAlignmentMatrixImpl`
pub struct PoaAlignmentMatrix {
    pub read_sequence: Vec<u8>,
    pub mode: AlignMode,
    pub columns: AlignmentColumnMap,
    pub score: f32,
}

/// POA 图(对应 `PoaGraphImpl`)
pub struct PoaGraph {
    nodes: Vec<PoaNode>,
    /// out-edges,插入序
    out_edges: Vec<Vec<Vertex>>,
    /// in-edges,插入序(排序视图按需生成)
    in_edges: Vec<Vec<Vertex>>,
    /// BGL setS 边存储自动去重的等价物
    edge_set: HashSet<(Vertex, Vertex)>,
    enter_vertex: Vertex,
    exit_vertex: Vertex,
    num_reads: usize,
}

impl PoaGraph {
    pub fn new() -> Self {
        let mut g = Self {
            nodes: Vec::new(),
            out_edges: Vec::new(),
            in_edges: Vec::new(),
            edge_set: HashSet::new(),
            enter_vertex: 0,
            exit_vertex: 0,
            num_reads: 0,
        };
        g.enter_vertex = g.add_vertex(ENTER_VERTEX_BASE, 0, 0);
        g.exit_vertex = g.add_vertex(EXIT_VERTEX_BASE, 0, 0);
        g
    }

    pub fn enter_vertex(&self) -> Vertex {
        self.enter_vertex
    }

    pub fn exit_vertex(&self) -> Vertex {
        self.exit_vertex
    }

    pub fn num_vertices(&self) -> usize {
        self.nodes.len()
    }

    pub fn num_reads(&self) -> usize {
        self.num_reads
    }

    pub fn node(&self, v: Vertex) -> &PoaNode {
        &self.nodes[v]
    }

    fn add_vertex(&mut self, base: u8, reads: i32, spanning: i32) -> Vertex {
        let v = self.nodes.len();
        self.nodes.push(PoaNode {
            id: v,
            base,
            reads,
            spanning_reads: spanning,
            score: 0.0,
            reaching_score: 0.0,
        });
        self.out_edges.push(Vec::new());
        self.in_edges.push(Vec::new());
        v
    }

    /// 对应 `add_edge`(BGL setS:等价边自动去重)
    fn add_edge(&mut self, u: Vertex, v: Vertex) {
        if self.edge_set.insert((u, v)) {
            self.out_edges[u].push(v);
            self.in_edges[v].push(u);
        }
    }

    /// 对应 `inEdges(v, g)`:按 (source, target) 排序(此处 target 恒同,
    /// 即按 source 升序)
    pub fn sorted_in_edges(&self, v: Vertex) -> Vec<Vertex> {
        let mut preds = self.in_edges[v].clone();
        preds.sort_unstable();
        preds
    }

    pub fn out_edges(&self, v: Vertex) -> &[Vertex] {
        &self.out_edges[v]
    }

    /// 全部边(调试/验证用):(src, dst) 升序
    pub fn edges(&self) -> Vec<(Vertex, Vertex)> {
        let mut out = Vec::new();
        for (u, dsts) in self.out_edges.iter().enumerate() {
            for &v in dsts {
                out.push((u, v));
            }
        }
        out.sort_unstable();
        out
    }

    /// 对应 `topological_sort`:DFS(顶点按插入序、out_edges 按插入序)的
    /// 逆后序,与 BGL 行为一致
    pub fn topological_order(&self) -> Vec<Vertex> {
        let n = self.nodes.len();
        let mut color = vec![0u8; n]; // 0=white, 1=gray, 2=black
        let mut post = Vec::with_capacity(n);
        let mut stack: Vec<(Vertex, usize)> = Vec::new();
        for root in 0..n {
            if color[root] != 0 {
                continue;
            }
            color[root] = 1;
            stack.push((root, 0));
            while let Some(last) = stack.last_mut() {
                let (v, ref mut next) = *last;
                if *next < self.out_edges[v].len() {
                    let u = self.out_edges[v][*next];
                    *next += 1;
                    if color[u] == 0 {
                        color[u] = 1;
                        stack.push((u, 0));
                    }
                } else {
                    color[v] = 2;
                    post.push(v);
                    stack.pop();
                }
            }
        }
        post.reverse();
        post
    }

    /// 对应 `sequenceAlongPath`
    pub fn sequence_along_path(&self, path: &[Vertex]) -> Vec<u8> {
        path.iter().map(|&v| self.nodes[v].base).collect()
    }

    // ---------------- 首条 read 穿线 ----------------

    /// 对应 `threadFirstRead`
    fn thread_first_read(&mut self, sequence: &[u8], mut out_path: Option<&mut Vec<Vertex>>) {
        let mut u = NULL_VERTEX;
        let mut read_pos = 0usize;
        let mut start_span_vertex = NULL_VERTEX;

        if let Some(p) = out_path.as_deref_mut() {
            p.clear();
        }

        // 与 C++ threadFirstRead 的 readPos 递增结构对应
        #[allow(clippy::explicit_counter_loop)]
        for &base in sequence {
            let v = self.add_vertex(base, 1, 0);
            if let Some(p) = out_path.as_deref_mut() {
                p.push(v);
            }
            if read_pos == 0 {
                self.add_edge(self.enter_vertex, v);
                start_span_vertex = v;
            } else {
                self.add_edge(u, v);
            }
            u = v;
            read_pos += 1;
        }
        debug_assert!(start_span_vertex != NULL_VERTEX && u != NULL_VERTEX);
        let end_span_vertex = u;
        self.add_edge(u, self.exit_vertex);
        self.tag_span(start_span_vertex, end_span_vertex);
    }

    // ---------------- SpanningDFS / tagSpan ----------------

    /// 对应 `SpanningDFS`:start 可达 ∩ 可达 end
    fn spanning_dfs(&self, start: Vertex, end: Vertex) -> HashSet<Vertex> {
        let mut stack = vec![start];
        let mut fwd: HashSet<Vertex> = HashSet::new();
        let mut rev: HashSet<Vertex> = HashSet::new();
        // start 可达的所有顶点
        while let Some(v) = stack.pop() {
            if !fwd.insert(v) {
                continue;
            }
            for &t in &self.out_edges[v] {
                stack.push(t);
            }
        }
        // 能到达 end(经由 start 可达集)的顶点
        stack.push(end);
        while let Some(v) = stack.pop() {
            if !fwd.contains(&v) || !rev.insert(v) {
                continue;
            }
            for &s in &self.in_edges[v] {
                stack.push(s);
            }
        }
        rev
    }

    fn tag_span(&mut self, start: Vertex, end: Vertex) {
        let vertices = self.spanning_dfs(start, end);
        for v in vertices {
            self.nodes[v].spanning_reads += 1;
        }
    }

    // ---------------- DP 列 ----------------

    fn predecessor_columns<'a>(
        &self,
        v: Vertex,
        col_map: &'a AlignmentColumnMap,
    ) -> Vec<&'a AlignmentColumn> {
        self.sorted_in_edges(v)
            .iter()
            .map(|&u| col_map.get(u))
            .collect()
    }

    /// 对应 `makeAlignmentColumn`(Sparse DP 列;begin_row > end_row 表示
    /// 无锚点,退化为单行 START 列)
    fn make_alignment_column(
        &self,
        v: Vertex,
        col_map: &AlignmentColumnMap,
        sequence: &[u8],
        config: &AlignConfig,
        begin_row: i32,
        end_row: i32,
    ) -> AlignmentColumn {
        if begin_row > end_row {
            // 无锚点:以 START move 穿成单例路径(仅 LOCAL 模式可达)
            let mut col = AlignmentColumn::new(v, 0, 1);
            col.reaching_move.set(0, MoveType::Start);
            col.previous_vertex.set(0, self.enter_vertex);
            col.score.set(0, 0.0);
            return col;
        }

        let mut col = AlignmentColumn::new(v, begin_row, end_row);
        let vertex_base = self.nodes[v].base;
        let predecessor_columns = self.predecessor_columns(v, col_map);

        for i in begin_row..end_row {
            let i_us = i as usize;

            let (mut best_score, mut prev_vertex, mut reaching_move) =
                if config.mode == AlignMode::Local {
                    (0.0f32, self.enter_vertex, MoveType::Start)
                } else {
                    (f32::NEG_INFINITY, NULL_VERTEX, MoveType::Invalid)
                };

            if i == 0 {
                if predecessor_columns.is_empty() {
                    // enter 顶点:无 reaching move
                    debug_assert_eq!(v, self.enter_vertex);
                    col.score.set(0, 0.0);
                    col.reaching_move.set(0, MoveType::Invalid);
                    col.previous_vertex.set(0, NULL_VERTEX);
                    continue;
                } else if config.mode == AlignMode::SemiGlobal || config.mode == AlignMode::Local {
                    col.score.set(0, 0.0);
                    col.reaching_move.set(0, MoveType::Start);
                    col.previous_vertex.set(0, self.enter_vertex);
                } else {
                    for prev_col in &predecessor_columns {
                        let candidate = prev_col.score.get(0) + config.delete as f32;
                        if candidate > best_score {
                            best_score = candidate;
                            prev_vertex = prev_col.current_vertex;
                            reaching_move = MoveType::Delete;
                        }
                    }
                    debug_assert!(reaching_move != MoveType::Invalid);
                    col.score.set(0, best_score);
                    col.reaching_move.set(0, reaching_move);
                    col.previous_vertex.set(0, prev_vertex);
                }
            } else {
                for prev_col in &predecessor_columns {
                    // Match / Mismatch
                    if prev_col.has_row(i_us - 1) {
                        let is_match = sequence[i_us - 1] == vertex_base;
                        let candidate = prev_col.score.get(i_us - 1)
                            + (if is_match {
                                config.match_score
                            } else {
                                config.mismatch_score
                            }) as f32;
                        if candidate > best_score {
                            best_score = candidate;
                            prev_vertex = prev_col.current_vertex;
                            reaching_move = if is_match {
                                MoveType::Match
                            } else {
                                MoveType::Mismatch
                            };
                        }
                    }
                    // Delete
                    if prev_col.has_row(i_us) {
                        let candidate = prev_col.score.get(i_us) + config.delete as f32;
                        if candidate > best_score {
                            best_score = candidate;
                            prev_vertex = prev_col.current_vertex;
                            reaching_move = MoveType::Delete;
                        }
                    }
                }
                // Extra(本列自引用)
                if col.has_row(i_us - 1) {
                    let candidate = col.score.get(i_us - 1) + config.insert as f32;
                    if candidate > best_score {
                        best_score = candidate;
                        prev_vertex = v;
                        reaching_move = MoveType::Extra;
                    }
                }
                debug_assert!(reaching_move != MoveType::Invalid);
                col.score.set(i_us, best_score);
                col.reaching_move.set(i_us, reaching_move);
                col.previous_vertex.set(i_us, prev_vertex);
            }
        }

        col
    }

    /// 对应 `makeAlignmentColumnForExit`
    fn make_alignment_column_for_exit(
        &self,
        col_map: &AlignmentColumnMap,
        sequence: &[u8],
        config: &AlignConfig,
    ) -> AlignmentColumn {
        let v = self.exit_vertex;
        let i_seq = sequence.len();
        let mut col = AlignmentColumn::new(v, 0, i_seq as i32 + 1);

        let mut best_score = f32::NEG_INFINITY;
        let mut prev_vertex = NULL_VERTEX;

        if config.mode == AlignMode::SemiGlobal || config.mode == AlignMode::Local {
            // 顶点按**插入序**遍历(= BGL vertices(g_) 的 listS 序),
            // 严格大于 first-wins —— 顺序敏感,不可改
            for u in 0..self.nodes.len() {
                if u != self.exit_vertex {
                    let pred_col = col_map.get(u);
                    let prev_row = if config.mode == AlignMode::Local {
                        pred_col.score.arg_max()
                    } else {
                        i_seq
                    };
                    if pred_col.has_row(prev_row) && pred_col.score.get(prev_row) > best_score {
                        best_score = pred_col.score.get(prev_row);
                        prev_vertex = pred_col.current_vertex;
                    }
                }
            }
        } else {
            let predecessor_columns = self.predecessor_columns(v, col_map);
            for pred_col in predecessor_columns {
                if pred_col.has_row(i_seq) && pred_col.score.get(i_seq) > best_score {
                    best_score = pred_col.score.get(i_seq);
                    prev_vertex = pred_col.current_vertex;
                }
            }
        }
        debug_assert!(prev_vertex != NULL_VERTEX);
        col.score.set(i_seq, best_score);
        col.previous_vertex.set(i_seq, prev_vertex);
        col.reaching_move.set(i_seq, MoveType::End);
        col
    }

    // ---------------- read 增量 ----------------

    /// 对应 `AddFirstRead`
    pub fn add_first_read(&mut self, sequence: &[u8], out_path: Option<&mut Vec<Vertex>>) {
        debug_assert!(self.num_reads == 0);
        self.thread_first_read(sequence, out_path);
        self.num_reads += 1;
    }

    /// 对应 `TryAddRead`(range finder 先行初始化,再按拓扑序铺列)
    pub fn try_add_read(
        &self,
        read_seq: &[u8],
        config: &AlignConfig,
        range_finder: &mut dyn SdpRangeFinder,
    ) -> PoaAlignmentMatrix {
        debug_assert!(self.num_reads > 0);

        // "中间共识"路径(无 minCoverage;首尾可能带额外序列)
        let css_path = self.consensus_path(config.mode, -INT_MAX_HALF);
        let css_seq = self.sequence_along_path(&css_path);
        range_finder.init(self, &css_path, &css_seq, read_seq);

        let mut mat = PoaAlignmentMatrix {
            read_sequence: read_seq.to_vec(),
            mode: config.mode,
            columns: AlignmentColumnMap::new(self.nodes.len()),
            score: 0.0,
        };

        let topo = self.topological_order();
        for &v in &topo {
            let cur_col = if v != self.exit_vertex {
                let (start_row, end_row) = range_finder.find_alignable_range(v);
                let end_row = if end_row == -INT_MAX_HALF {
                    end_row
                } else {
                    end_row + 1
                };
                self.make_alignment_column(v, &mat.columns, read_seq, config, start_row, end_row)
            } else {
                self.make_alignment_column_for_exit(&mat.columns, read_seq, config)
            };
            mat.columns.insert(v, cur_col);
        }

        mat.score = mat.columns.get(self.exit_vertex).score.get(read_seq.len());
        mat
    }

    /// 对应 `CommitAdd`
    pub fn commit_add(&mut self, mat: PoaAlignmentMatrix, out_path: Option<&mut Vec<Vertex>>) {
        self.traceback_and_thread(&mat.read_sequence, &mat.columns, mat.mode, out_path);
        self.num_reads += 1;
    }

    // ---------------- traversals ----------------

    /// 对应 `tracebackAndThread` 状态机
    fn traceback_and_thread(
        &mut self,
        sequence: &[u8],
        columns: &AlignmentColumnMap,
        align_mode: AlignMode,
        mut out_path: Option<&mut Vec<Vertex>>,
    ) {
        let i_seq = sequence.len() as i32;
        let mut i = i_seq;
        let mut fork_vertex = NULL_VERTEX;
        let mut u = self.exit_vertex;
        let mut v = NULL_VERTEX;
        let mut span = 0i32;
        let mut start_span_vertex;
        let end_span_vertex = columns
            .get(self.exit_vertex)
            .previous_vertex
            .get(i_seq as usize);

        if let Some(p) = out_path.as_deref_mut() {
            p.resize(i_seq as usize, NULL_VERTEX);
        }

        while !(u == self.enter_vertex && i == 0) {
            let cur_col = columns.get(u);
            let cur_node_spanning = self.nodes[u].spanning_reads;
            let prev_vertex = cur_col.previous_vertex.get(i as usize);
            let reaching_move = cur_col.reaching_move.get(i as usize);

            match reaching_move {
                MoveType::Start => {
                    debug_assert!(v != NULL_VERTEX);
                    if fork_vertex == NULL_VERTEX {
                        fork_vertex = v;
                    }
                    while i > 0 {
                        debug_assert!(align_mode == AlignMode::Local);
                        let read_pos = (i - 1) as usize;
                        let new_fork = self.add_vertex(sequence[read_pos], 1, span);
                        self.add_edge(new_fork, fork_vertex);
                        if let Some(pp) = out_path.as_deref_mut() {
                            pp[read_pos] = new_fork;
                        }
                        fork_vertex = new_fork;
                        i -= 1;
                    }
                }
                MoveType::End => {
                    debug_assert!(
                        fork_vertex == NULL_VERTEX && u == self.exit_vertex && v == NULL_VERTEX
                    );
                    fork_vertex = self.exit_vertex;

                    if align_mode == AlignMode::Local {
                        let prev_col = columns.get(prev_vertex);
                        let prev_row = prev_col.score.arg_max();
                        while i > prev_row as i32 {
                            let read_pos = (i - 1) as usize;
                            let new_fork = self.add_vertex(sequence[read_pos], 1, span);
                            self.add_edge(new_fork, fork_vertex);
                            if let Some(pp) = out_path.as_deref_mut() {
                                pp[read_pos] = new_fork;
                            }
                            fork_vertex = new_fork;
                            i -= 1;
                        }
                    }
                }
                MoveType::Match => {
                    if let Some(pp) = out_path.as_deref_mut() {
                        pp[(i - 1) as usize] = u;
                    }
                    if fork_vertex != NULL_VERTEX {
                        self.add_edge(u, fork_vertex);
                        fork_vertex = NULL_VERTEX;
                    }
                    self.nodes[u].reads += 1;
                    i -= 1;
                }
                MoveType::Delete => {
                    if fork_vertex == NULL_VERTEX {
                        fork_vertex = v;
                    }
                }
                MoveType::Extra | MoveType::Mismatch => {
                    let read_pos = (i - 1) as usize;
                    let new_fork = self.add_vertex(sequence[read_pos], 1, span);
                    if fork_vertex == NULL_VERTEX {
                        fork_vertex = v;
                    }
                    self.add_edge(new_fork, fork_vertex);
                    if let Some(pp) = out_path.as_deref_mut() {
                        pp[read_pos] = new_fork;
                    }
                    fork_vertex = new_fork;
                    i -= 1;
                }
                MoveType::Invalid => panic!("unreachable: Invalid move in traceback"),
            }

            v = u;
            u = prev_vertex;
            span = cur_node_spanning;
        }
        start_span_vertex = v;

        if fork_vertex != NULL_VERTEX {
            self.add_edge(self.enter_vertex, fork_vertex);
            start_span_vertex = fork_vertex;
            fork_vertex = NULL_VERTEX;
        }

        if start_span_vertex != self.exit_vertex {
            self.tag_span(start_span_vertex, end_span_vertex);
        }
        let _ = fork_vertex; // 与 C++ 结构对应(结束前置空)

        // 全部填充校验(C++ 的 assert)
        debug_assert!(
            out_path.is_none() || out_path.as_ref().unwrap().iter().all(|&x| x != NULL_VERTEX)
        );
    }

    /// 对应 `consensusPath`。结果与拓扑序实现无关(依赖显式 tie-break),
    /// 但 bestReachingScore 的"严格大于"更新依赖遍历序 —— 这里用与 BGL
    /// 一致的 DFS 拓扑序 + 排序 in-edges。
    pub fn consensus_path(&self, mode: AlignMode, min_coverage: i32) -> Vec<Vertex> {
        let total_reads = self.num_reads as i32;

        let mut topo = self.topological_order();
        // 排除首尾(^ 与 $;拓扑序中 enter 必为第一个、exit 必为最后)
        debug_assert!(topo[0] == self.enter_vertex && topo[topo.len() - 1] == self.exit_vertex);
        topo.remove(topo.len() - 1);
        topo.remove(0);

        let n = self.nodes.len();
        let mut best_prev_vertex: Vec<Vertex> = vec![NULL_VERTEX; n];
        // C++ 把 Score/ReachingScore 存在节点上(副作用仅供 GraphViz 展示);
        // 这里用局部数组承载,避免 &self 下的内部可变性
        let mut reaching: Vec<f32> = vec![0.0; n];

        let mut best_vertex = NULL_VERTEX;
        let mut best_reaching_score = f32::NEG_INFINITY;
        for &v in &topo {
            let containing_reads = self.nodes[v].reads;
            let spanning_reads = self.nodes[v].spanning_reads;
            let score = if mode != AlignMode::Global {
                2.0 * containing_reads as f32
                    - 1.0 * (spanning_reads.max(min_coverage)) as f32
                    - 0.0001f32
            } else {
                2.0 * containing_reads as f32 - 1.0 * total_reads as f32 - 0.0001f32
            };
            reaching[v] = score;
            best_prev_vertex[v] = NULL_VERTEX;
            for &source_vertex in &self.sorted_in_edges(v) {
                let rsc = score + reaching[source_vertex];
                if rsc > reaching[v] {
                    reaching[v] = rsc;
                    best_prev_vertex[v] = source_vertex;
                }
                if rsc > best_reaching_score {
                    best_vertex = v;
                    best_reaching_score = rsc;
                } else if rsc == best_reaching_score {
                    // tie-break:顶点 index 小者胜(与遍历顺序无关)
                    if v < best_vertex {
                        best_vertex = v;
                    }
                }
            }
        }
        debug_assert!(best_vertex != NULL_VERTEX);

        // 自 best 回溯
        let mut path = Vec::new();
        let mut v = best_vertex;
        while v != NULL_VERTEX {
            path.push(v);
            v = best_prev_vertex[v];
        }
        path.reverse();
        path
    }

    /// 对应 `PoaGraphImpl::FindConsensus`
    pub fn find_consensus(
        &self,
        config: &AlignConfig,
        min_coverage: i32,
    ) -> (Vec<u8>, Vec<Vertex>) {
        let best_path = self.consensus_path(config.mode, min_coverage);
        let consensus = self.sequence_along_path(&best_path);
        (consensus, best_path)
    }
}

impl Default for PoaGraph {
    fn default() -> Self {
        Self::new()
    }
}
