//! 种子类型(对应 cpp_code/Seed.h / Seed-inl.h)。
//!
//! H = query(read),V = reference(consensus)。

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Seed {
    begin_pos_h: u64,
    begin_pos_v: u64,
    end_pos_h: u64,
    end_pos_v: u64,
    #[allow(dead_code)]
    lower_diagonal: i64,
    #[allow(dead_code)]
    upper_diagonal: i64,
    score: i32,
}

impl Seed {
    /// 3 参构造(begin + length),`lower/upper_diagonal = beginH - beginV`
    pub fn new(begin_pos_h: u64, begin_pos_v: u64, seed_length: u64) -> Self {
        let diag = (begin_pos_h as i64) - (begin_pos_v as i64);
        Self {
            begin_pos_h,
            begin_pos_v,
            end_pos_h: begin_pos_h + seed_length,
            end_pos_v: begin_pos_v + seed_length,
            lower_diagonal: diag,
            upper_diagonal: diag,
            score: 0,
        }
    }

    #[inline]
    pub fn begin_pos_h(&self) -> u64 {
        self.begin_pos_h
    }
    #[inline]
    pub fn begin_pos_v(&self) -> u64 {
        self.begin_pos_v
    }
    #[inline]
    pub fn end_pos_h(&self) -> u64 {
        self.end_pos_h
    }
    #[inline]
    pub fn end_pos_v(&self) -> u64 {
        self.end_pos_v
    }
    #[inline]
    pub fn score(&self) -> i32 {
        self.score
    }

    /// `BeginDiagonal = beginH - beginV`
    #[inline]
    pub fn begin_diagonal(&self) -> i64 {
        self.begin_pos_h as i64 - self.begin_pos_v as i64
    }

    /// `Diagonal`:带符号的对角距离(C++ 的 Diagonal() 特殊语义)
    #[inline]
    pub fn diagonal(&self) -> i64 {
        if self.begin_pos_v > self.begin_pos_h {
            -((self.begin_pos_v - self.begin_pos_h) as i64)
        } else {
            (self.begin_pos_h - self.begin_pos_v) as i64
        }
    }

    /// `Size = max(endH-beginH, endV-beginV)`
    #[inline]
    pub fn size(&self) -> u64 {
        (self.end_pos_h - self.begin_pos_h).max(self.end_pos_v - self.begin_pos_v)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seed_geometry() {
        let s = Seed::new(10, 4, 6);
        assert_eq!(s.end_pos_h(), 16);
        assert_eq!(s.end_pos_v(), 10);
        assert_eq!(s.begin_diagonal(), 6);
        assert_eq!(s.size(), 6);
        let t = Seed::new(3, 9, 5);
        assert_eq!(t.begin_diagonal(), -6);
        assert_eq!(t.diagonal(), -6);
        assert_eq!(Seed::new(5, 5, 4).diagonal(), 0);
    }
}
