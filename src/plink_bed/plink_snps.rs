use analytic::traits::ToIterator;

use crate::plink_bed::{
    geno_to_lowest_two_bits, get_num_people_last_byte, lowest_two_bits_to_geno,
    usize_div_ceil,
};

pub struct PlinkSnps {
    bytes: Vec<u8>,
    num_snps: usize,
}

impl PlinkSnps {
    pub fn new(bytes: Vec<u8>, num_snps: usize) -> PlinkSnps {
        assert!(num_snps <= bytes.len() * 4,
                format!("num_snps ({}) > bytes.len() * 4 = {}", num_snps, bytes.len() * 4));
        PlinkSnps {
            bytes,
            num_snps,
        }
    }

    pub fn from_geno(geno: Vec<u8>) -> PlinkSnps {
        let num_bytes = usize_div_ceil(geno.len(), 4);
        let mut bytes: Vec<u8> = Vec::with_capacity(num_bytes);
        let mut snp_index = 0;
        if let Some(num_people_last_byte) = get_num_people_last_byte(geno.len()) {
            for _ in 0..num_bytes - 1 {
                bytes.push(
                    geno_to_lowest_two_bits(geno[snp_index])
                        | (geno_to_lowest_two_bits(geno[snp_index + 1]) << 2)
                        | (geno_to_lowest_two_bits(geno[snp_index + 2]) << 4)
                        | (geno_to_lowest_two_bits(geno[snp_index + 3]) << 6)
                );
                snp_index += 4;
            }
            // last byte
            let mut last_byte = 0u8;
            for j in 0..num_people_last_byte {
                last_byte |= geno_to_lowest_two_bits(geno[snp_index + j]) << (j * 2);
            }
            bytes.push(last_byte);
        }
        PlinkSnps::new(bytes, geno.len())
    }

    #[inline]
    pub fn get_num_snps(&self) -> usize {
        self.num_snps
    }

    #[inline]
    pub fn to_bytes(&self) -> &Vec<u8> {
        &self.bytes
    }

    #[inline]
    pub fn into_bytes(self) -> Vec<u8> {
        self.bytes
    }
}

impl ToIterator<'_, PlinkSnpsIter, u8> for PlinkSnps {
    fn to_iter(&self) -> PlinkSnpsIter {
        PlinkSnpsIter {
            bytes: self.bytes.clone(),
            num_snps: self.num_snps,
            byte_cursor: 0,
            cursor: 0,
        }
    }
}

impl IntoIterator for PlinkSnps {
    type Item = <PlinkSnpsIter as Iterator>::Item;
    type IntoIter = PlinkSnpsIter;
    fn into_iter(self) -> Self::IntoIter {
        PlinkSnpsIter {
            bytes: self.bytes,
            num_snps: self.num_snps,
            byte_cursor: 0,
            cursor: 0,
        }
    }
}

pub struct PlinkSnpsIter {
    bytes: Vec<u8>,
    num_snps: usize,
    byte_cursor: usize,
    cursor: usize,
}

impl Iterator for PlinkSnpsIter {
    type Item = u8;
    fn next(&mut self) -> Option<Self::Item> {
        if self.cursor >= self.num_snps {
            None
        } else {
            let snp = lowest_two_bits_to_geno(self.bytes[self.byte_cursor] >> (2 * (self.cursor % 4) as u8));
            self.cursor += 1;
            if self.cursor % 4 == 0 {
                self.byte_cursor += 1;
            }
            Some(snp)
        }
    }
}

#[cfg(test)]
mod tests {
    use analytic::traits::ToIterator;

    use super::PlinkSnps;

    #[test]
    fn test_plink_snps() {
        let expected_num_snps = 12;
        let snps = PlinkSnps::new(
            vec![0b10_00_11_00, 0b00_00_11_11, 0b11_10_10_10, 0b11001100],
            expected_num_snps,
        );
        let expected: Vec<u8> = vec![2, 0, 2, 1, 0, 0, 2, 2, 1, 1, 1, 0];
        let mut num_snps = 0;
        for (s, e) in snps.to_iter().zip(expected.iter()) {
            assert_eq!(s, *e);
            num_snps += 1;
        }
        assert_eq!(num_snps, expected_num_snps);
    }

    #[test]
    fn test_plink_snps_from_geno() {
        fn test(geno: Vec<u8>) {
            let plink_snps = PlinkSnps::from_geno(geno.clone());
            let num_snps = geno.len();
            let mut iter = plink_snps.into_iter();
            for i in 0..num_snps {
                assert_eq!(Some(geno[i]), iter.next());
            }
            assert_eq!(None, iter.next());
        }
        test(vec![]);
        test(vec![0]);
        test(vec![1]);
        test(vec![2]);
        test(vec![0, 1, 1]);
        test(vec![2, 2, 0, 1]);
        test(vec![1, 2, 0, 1, 1, 0]);
        test(vec![1, 2, 0, 1, 1, 0, 2, 1]);
        test(vec![1, 2, 2, 0, 0, 1, 0, 0, 2, 0]);
    }
}
