//! Sparse POA 编排(对应 cpp_code/SparsePoa.h / SparsePoa.cpp)。
//!
//! 单一职责:read 的增量编排(定向、评分阈值)与共识摘要;
//! 图算法在 `poa_graph`,区间定位在 `range_finder`。

use std::collections::HashMap;

use crate::native::poa_graph::{AlignConfig, AlignMode, PoaGraph, Vertex};
use crate::native::range_finder::{SdpRangeFinder, SdpRangeFinderV1, SdpRangeFinderV2};
use crate::native::sequence::reverse_complement;
use crate::setting::PoaSetting;

/// 对应 `SparsePoa::ReadKey`(-1 = 被拒)
pub type ReadKey = i32;

/// 对应 `PoaAlignmentSummary`(仅保留 FFI 路径消费的字段)
#[derive(Debug, Clone)]
pub struct PoaAlignmentSummary {
    pub reverse_complemented_read: bool,
    pub extent_on_read: (usize, usize),
    pub extent_on_consensus: (usize, usize),
    pub alignment_identity: f32,
}

/// settings → LOCAL 模式 AlignConfig(对应 DefaultPoaConfig(LOCAL) + 覆盖)
fn local_config(setting: &PoaSetting) -> AlignConfig {
    AlignConfig {
        match_score: setting.match_score,
        mismatch_score: setting.mismatch_score,
        insert: setting.insertion_score,
        delete: setting.deletion_score,
        mode: AlignMode::Local,
    }
}

/// 对应 `ConsensusSpan`:计算 read 路径与共识的交集与错误数
fn consensus_span(
    read_path: &[Vertex],
    css_position: &HashMap<Vertex, usize>,
) -> (usize, usize, usize, usize, usize) {
    let (mut css_s, mut css_e, mut read_s, mut read_e) = (0usize, 0usize, 0usize, 0usize);
    let mut n_err = 0usize;
    let mut found_start = false;
    for (read_pos, &v) in read_path.iter().enumerate() {
        if let Some(&css_pos) = css_position.get(&v) {
            if !found_start {
                css_s = css_pos;
                read_s = read_pos;
                found_start = true;
            }
            css_e = css_pos + 1;
            read_e = read_pos + 1;
        } else {
            n_err += 1;
        }
    }
    (css_s, css_e, read_s, read_e, n_err)
}

/// Sparse partial-order aligner(对应 `PacBio::Poa::SparsePoa`)
pub struct SparsePoa {
    graph: PoaGraph,
    read_paths: Vec<Vec<Vertex>>,
    reverse_complemented: Vec<bool>,
    range_finder: Box<dyn SdpRangeFinder>,
}

impl SparsePoa {
    /// 对应默认构造 `SparsePoa()`:固定 V1 —— **FwdStrand 入口的怪癖来源**
    pub fn new_v1() -> Self {
        Self {
            graph: PoaGraph::new(),
            read_paths: Vec::new(),
            reverse_complemented: Vec::new(),
            range_finder: Box::new(SdpRangeFinderV1::new()),
        }
    }

    /// 对应 `SparsePoa(int version)`;非 1/2 返回 `Err`(C++ 抛异常)
    pub fn new(version: i32) -> Result<Self, crate::error::XpoaError> {
        if version != 1 && version != 2 {
            return Err(crate::error::XpoaError::InvalidVersion { version });
        }
        let rf: Box<dyn SdpRangeFinder> = if version == 1 {
            Box::new(SdpRangeFinderV1::new())
        } else {
            Box::new(SdpRangeFinderV2::new())
        };
        Ok(Self {
            graph: PoaGraph::new(),
            read_paths: Vec::new(),
            reverse_complemented: Vec::new(),
            range_finder: rf,
        })
    }

    pub fn num_reads(&self) -> usize {
        self.graph.num_reads()
    }

    /// 只读图访问(调试/验证用)
    pub fn graph(&self) -> &PoaGraph {
        &self.graph
    }

    /// 对应 `AddRead`(方向已确定的 read;`rc` 仅记录,不参与定向)
    pub fn add_read(&mut self, read_sequence: &[u8], setting: &PoaSetting, rc: bool) -> ReadKey {
        let config = local_config(setting);
        // C++ minScoreToAdd 默认 0
        const MIN_SCORE_TO_ADD: f32 = 0.0;

        if self.graph.num_reads() == 0 {
            let mut path = Vec::new();
            self.graph.add_first_read(read_sequence, Some(&mut path));
            self.read_paths.push(path);
            self.reverse_complemented.push(false);
            (self.graph.num_reads() - 1) as i32
        } else {
            let c1 = self
                .graph
                .try_add_read(read_sequence, &config, self.range_finder.as_mut());
            if c1.score >= MIN_SCORE_TO_ADD {
                let mut path = Vec::new();
                self.graph.commit_add(c1, Some(&mut path));
                self.read_paths.push(path);
                self.reverse_complemented.push(rc);
                (self.graph.num_reads() - 1) as i32
            } else {
                -1
            }
        }
    }

    /// 对应 `OrientAndAddRead`:正反向各试一次,取分高者
    pub fn orient_and_add_read(
        &mut self,
        read_sequence: &[u8],
        setting: &PoaSetting,
    ) -> Result<ReadKey, crate::error::XpoaError> {
        let config = local_config(setting);
        const MIN_SCORE_TO_ADD: f32 = 0.0;

        if self.graph.num_reads() == 0 {
            let mut path = Vec::new();
            self.graph.add_first_read(read_sequence, Some(&mut path));
            self.read_paths.push(path);
            self.reverse_complemented.push(false);
            Ok((self.graph.num_reads() - 1) as i32)
        } else {
            let c1 = self
                .graph
                .try_add_read(read_sequence, &config, self.range_finder.as_mut());
            let rc_seq = reverse_complement(read_sequence)?;
            let c2 = self
                .graph
                .try_add_read(&rc_seq, &config, self.range_finder.as_mut());

            if c1.score >= c2.score && c1.score >= MIN_SCORE_TO_ADD {
                let mut path = Vec::new();
                self.graph.commit_add(c1, Some(&mut path));
                self.read_paths.push(path);
                self.reverse_complemented.push(false);
                Ok((self.graph.num_reads() - 1) as i32)
            } else if c2.score >= c1.score && c2.score >= MIN_SCORE_TO_ADD {
                let mut path = Vec::new();
                self.graph.commit_add(c2, Some(&mut path));
                self.read_paths.push(path);
                self.reverse_complemented.push(true);
                Ok((self.graph.num_reads() - 1) as i32)
            } else {
                Ok(-1)
            }
        }
    }

    /// 对应 `FindConsensus`:返回 (共识序列, 各 read 摘要)
    pub fn find_consensus(
        &self,
        min_coverage: i32,
        _setting: &PoaSetting,
    ) -> (Vec<u8>, Vec<PoaAlignmentSummary>) {
        let config = local_config(_setting);
        let (css, css_path) = self.graph.find_consensus(&config, min_coverage);

        let mut summaries = Vec::with_capacity(self.read_paths.len());
        let css_position: HashMap<Vertex, usize> =
            css_path.iter().enumerate().map(|(i, &v)| (v, i)).collect();

        for (read_id, read_path) in self.read_paths.iter().enumerate() {
            let (css_s, css_e, read_s, read_e, n_err) = consensus_span(read_path, &css_position);
            let identity = (1.0 - 1.0 * n_err as f32 / css_position.len() as f32).max(0.0);
            summaries.push(PoaAlignmentSummary {
                reverse_complemented_read: self.reverse_complemented[read_id],
                extent_on_read: (read_s, read_e),
                extent_on_consensus: (css_s, css_e),
                alignment_identity: identity,
            });
        }
        (css, summaries)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_reads_chain() {
        let setting = crate::xpoa::get_default_poa_setting();
        let mut poa = SparsePoa::new(1).unwrap();
        for _ in 0..4 {
            let key = poa.add_read(b"AACGGATCGGA", &setting, false);
            assert!(key >= 0);
        }
        let (css, summaries) = poa.find_consensus(1, &setting);
        assert_eq!(css, b"AACGGATCGGA".to_vec());
        assert_eq!(summaries.len(), 4);
        assert!(summaries.iter().all(|s| s.alignment_identity > 0.99));
    }

    #[test]
    fn consensus_span_counts_errors() {
        let mut css_position = HashMap::new();
        css_position.insert(10usize, 0usize);
        css_position.insert(11usize, 1usize);
        // read 路径:两个共识顶点 + 1 个插入(错误)
        let (css_s, css_e, read_s, read_e, n_err) = consensus_span(&[10, 99, 11], &css_position);
        assert_eq!((css_s, css_e), (0, 2));
        assert_eq!((read_s, read_e), (0, 3));
        assert_eq!(n_err, 1);
    }
}
