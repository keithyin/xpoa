pub mod xpoa_sys;

pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

#[cfg(test)]
mod tests {
    use crate::xpoa_sys::PoaDraftGen;

    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }

    fn test_poa_consensus() {
        let sequences = vec!["AACGGATCGGA", "AACGGATCGGA", "AACGGATCGGA", "AACGGATCGGA"];
        
        // let res = PoaDraftGen(&sequences, &PoaSetting::default());
        println!("res:{:?}", res);
    }
}
