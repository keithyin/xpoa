use std::mem::MaybeUninit;

pub static COMPLEMENT_TABLE: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'A' as usize] = b'T';
    table[b'T' as usize] = b'A';
    table[b'C' as usize] = b'G';
    table[b'G' as usize] = b'C';
    table[b'-' as usize] = b'-';
    table[b'*' as usize] = b'*';
    table[b'N' as usize] = b'N';

    table[b'a' as usize] = b't';
    table[b't' as usize] = b'a';
    table[b'c' as usize] = b'g';
    table[b'g' as usize] = b'c';
    table[b'-' as usize] = b'-';
    table[b'*' as usize] = b'*';
    table[b'n' as usize] = b'n';

    table
};

/// b"ACGTAA" -> b"TTACGT"
///
/// 表外字节(非 ACGTN/-/*)无法互补,返回 `Err`(对应 C++ 的 throw)。
pub fn reverse_complement(dna_seq: &[u8]) -> Result<Vec<u8>, crate::error::XpoaError> {
    let mut result: Vec<MaybeUninit<u8>> = vec![MaybeUninit::uninit(); dna_seq.len()];
    for (&base, slot) in dna_seq.iter().rev().zip(result.iter_mut()) {
        let complement = COMPLEMENT_TABLE[base as usize];
        if complement == 0 {
            return Err(crate::error::XpoaError::InvalidBase { base });
        }
        slot.write(complement);
    }
    Ok(unsafe { std::mem::transmute::<Vec<MaybeUninit<u8>>, Vec<u8>>(result) })
}
