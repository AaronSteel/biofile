//! # General traits

pub mod trait_impl;

/// Data type of the start and end coordinates
pub type Coordinate = i64;

/// Data type of the chromosomes
pub type Chrom = String;

/// Returns a tuple of the form
/// (chromosome, start_inclusive, end_exclusive, value)
pub trait ToChromStartEndVal<V> {
    fn to_chrom_start_end_val(
        &self,
    ) -> (Chrom, Coordinate, Coordinate, Option<V>);
}
