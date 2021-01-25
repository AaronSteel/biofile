use crate::error::Error;
use math::{
    partition::integer_interval_map::IntegerIntervalMap,
    set::{
        contiguous_integer_set::ContiguousIntegerSet,
        ordered_integer_set::OrderedIntegerSet, traits::Intersect,
    },
};
use num::{Integer, Num};
use std::collections::HashMap;

pub trait ChromIntervalValue<T, V>
where
    T: Copy + Integer, {
    fn chrom_interval_value(&self) -> (String, ContiguousIntegerSet<T>, V);
}

pub trait ToChromIntervalValueIter<I, C, T, V>
where
    I: Iterator<Item = C>,
    C: ChromIntervalValue<T, V>,
    T: Copy + Integer, {
    fn to_chrom_interval_value_iter(&self) -> I;
}

pub type Chrom = String;

type Coordinate = i64;

impl<I, C, V> dyn ToChromIntervalValueIter<I, C, Coordinate, V>
where
    I: Iterator<Item = C>,
    C: ChromIntervalValue<Coordinate, V>,
    V: Copy + Num,
{
    /// Will discard the lines in the bed file if the corresponding range has a
    /// non-empty intersection with any of the intervals in `exclude`.
    pub fn get_chrom_to_interval_to_val(
        &self,
        exclude: Option<&HashMap<Chrom, OrderedIntegerSet<Coordinate>>>,
    ) -> Result<HashMap<String, IntegerIntervalMap<V>>, Error> {
        let mut chrom_to_interval_map = HashMap::new();
        for item in self.to_chrom_interval_value_iter() {
            let (chrom, interval, value) = item.chrom_interval_value();
            if let Some(chrom_to_excluded_intervals) = exclude {
                if let Some(excluded_intervals) =
                    chrom_to_excluded_intervals.get(&chrom)
                {
                    if interval
                        .has_non_empty_intersection_with(excluded_intervals)
                    {
                        continue;
                    }
                }
            }

            let interval_map = chrom_to_interval_map
                .entry(chrom)
                .or_insert_with(IntegerIntervalMap::new);

            interval_map.aggregate(interval, value);
        }
        Ok(chrom_to_interval_map)
    }
}
