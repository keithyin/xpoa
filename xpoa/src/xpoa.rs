use crate::setting::PoaSetting;

pub trait TSubread {
    fn get_seq(&self) -> &str;
    fn get_cx(&self) -> u8;
    fn is_fwd(&self) -> Option<bool> {
        None
    }
}

impl TSubread for &String {
    fn get_seq(&self) -> &str {
        self
    }

    fn get_cx(&self) -> u8 {
        3
    }
}

impl TSubread for String {
    fn get_seq(&self) -> &str {
        self
    }

    fn get_cx(&self) -> u8 {
        3
    }
}

impl TSubread for &str {
    fn get_seq(&self) -> &str {
        self
    }

    fn get_cx(&self) -> u8 {
        3
    }
}

pub fn get_default_poa_setting() -> PoaSetting {
    PoaSetting {
        min_identity: 0.82,
        match_score: 3,
        mismatch_score: -5,
        insertion_score: -2,
        deletion_score: -2,
        ed_unify_strand: 1,
        version: 1, // version 2 意味着使用 range_finder2
    }
}
