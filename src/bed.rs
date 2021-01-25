//! An interface to the BED track format file as specified in
//! https://genome.ucsc.edu/FAQ/FAQformat.html#format1

use crate::util::{get_buf, Strand};
use math::{
    partition::integer_interval_map::IntegerIntervalMap,
    set::{
        contiguous_integer_set::ContiguousIntegerSet,
        ordered_integer_set::OrderedIntegerSet,
    },
    traits::ToIterator,
};
use num::Float;
use std::{
    collections::HashMap,
    fmt::Debug,
    fs::File,
    io::{BufRead, BufReader},
    marker::PhantomData,
    str::FromStr,
};

pub mod bed_writer;
pub mod paired_end_collator;

use crate::iter::{ChromIntervalValue, ToChromIntervalValueIter};
pub use bed_writer::BedWriter;

pub struct Bed {
    filepath: String,

    /// Every line in the BED file will contribute a unit score for the
    /// corresponding interval.
    binarize_score: bool,
}

impl Bed {
    pub fn new(filepath: &str, use_binary_score: bool) -> Bed {
        Bed {
            filepath: filepath.to_string(),
            binarize_score: use_binary_score,
        }
    }

    #[inline]
    pub fn get_filepath(&self) -> &str {
        &self.filepath
    }

    // /// Will discard the lines in the bed file if the corresponding range has
    // a /// non-empty intersection with any of the intervals in `exclude`.
    // pub fn get_chrom_to_interval_to_val<D, E>(
    //     &self,
    //     exclude: Option<&HashMap<Chrom, OrderedIntegerSet<Coordinate>>>,
    // ) -> Result<HashMap<String, IntegerIntervalMap<D>>, Error>
    // where
    //     D: Float + FromStr<Err = E>,
    //     E: Debug, {
    //     let mut chrom_to_interval_map = HashMap::new();
    //     for BedDataLine {
    //         chrom,
    //         start,
    //         end,
    //         name: _,
    //         score,
    //         strand: _,
    //     } in self.to_iter(): BedDataLineIter<D>
    //     {
    //         let score = if let Some(score) = score {
    //             score
    //         } else {
    //             return Err(Error::Generic(
    //                 "the BED file does not have a score field".into(),
    //             ));
    //         };
    //
    //         let interval = ContiguousIntegerSet::new(start, end - 1);
    //         if let Some(chrom_to_excluded_intervals) = exclude {
    //             if let Some(excluded_intervals) =
    //                 chrom_to_excluded_intervals.get(&chrom)
    //             {
    //                 if interval
    //                     .has_non_empty_intersection_with(excluded_intervals)
    //                 {
    //                     continue;
    //                 }
    //             }
    //         }
    //
    //         let interval_map = chrom_to_interval_map
    //             .entry(chrom)
    //             .or_insert_with(IntegerIntervalMap::new);
    //
    //         interval_map.aggregate(interval, score);
    //     }
    //     Ok(chrom_to_interval_map)
    // }

    pub fn get_chrom_to_intervals(
        &self,
    ) -> HashMap<Chrom, OrderedIntegerSet<Coordinate>> {
        let mut chrom_to_interval_map = HashMap::new();
        for (chrom, start, end) in self.to_coord_iter() {
            let interval_map = chrom_to_interval_map
                .entry(chrom)
                .or_insert_with(IntegerIntervalMap::new);
            interval_map
                .aggregate(ContiguousIntegerSet::new(start, end - 1), 1);
        }
        chrom_to_interval_map
            .into_iter()
            .map(|(chrom, interval_map)| {
                let intervals: Vec<ContiguousIntegerSet<Coordinate>> =
                    interval_map
                        .into_map()
                        .into_iter()
                        .map(|(k, _)| k)
                        .collect();
                (chrom, OrderedIntegerSet::from(intervals))
            })
            .collect()
    }

    pub fn to_coord_iter(&self) -> BedCoordinateIter {
        BedCoordinateIter {
            buf: get_buf(&self.filepath).unwrap(),
            filename: self.filepath.clone(),
        }
    }
}

impl<D, E>
    ToIterator<'_, BedDataLineIter<D>, <BedDataLineIter<D> as Iterator>::Item>
    for Bed
where
    D: Float + FromStr<Err = E>,
    E: Debug,
{
    fn to_iter(&self) -> BedDataLineIter<D> {
        BedDataLineIter {
            buf: get_buf(&self.filepath).unwrap(),
            filename: self.filepath.clone(),
            binarize_score: self.binarize_score,
            phantom: PhantomData,
        }
    }
}

impl<V: Float + FromStr<Err = E>, E: Debug>
    ToChromIntervalValueIter<BedDataLineIter<V>, BedDataLine<V>, Coordinate, V>
    for Bed
{
    fn to_chrom_interval_value_iter(&self) -> BedDataLineIter<V> {
        self.to_iter()
    }
}

impl<V: Copy + Float + FromStr<Err = E>, E: Debug>
    ChromIntervalValue<Coordinate, V> for BedDataLine<V>
{
    fn chrom_interval_value(
        &self,
    ) -> (Chrom, ContiguousIntegerSet<Coordinate>, V) {
        let score = self.score.expect("BED file is missing the score field");
        (
            self.chrom.clone(),
            ContiguousIntegerSet::new(self.start, self.end - 1),
            score,
        )
    }
}

/// Data type of the Bed coordinates
pub type Coordinate = i64;

/// Data type of the chromosome names
pub type Chrom = String;

/// `BedDataLine` corresponds to a line of data in the Bed
/// file, where each line is of the form
/// `chrom start end name score strand ...`,
/// where the first three fields are required, and the remaining 9 fields are
/// optional.
///
/// The [start, end) is a zero-based left-closed right-open coordinate range.
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct BedDataLine<D> {
    pub chrom: Chrom,
    pub start: Coordinate,
    pub end: Coordinate,
    pub name: Option<String>,
    pub score: Option<D>,
    pub strand: Option<Strand>,
}

pub struct BedDataLineIter<D> {
    buf: BufReader<File>,
    filename: String,

    /// Every line in the BED file will contribute a unit score for the
    /// corresponding interval.
    binarize_score: bool,
    phantom: PhantomData<D>,
}

impl<D> BedDataLineIter<D> {
    pub fn get_filename(&self) -> &str {
        &self.filename
    }
}

impl<D: Float + FromStr<Err = E>, E: Debug> Iterator for BedDataLineIter<D> {
    type Item = BedDataLine<D>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let mut line = String::new();
            return if self.buf.read_line(&mut line).unwrap() == 0 {
                None
            } else {
                let mut toks = line.split_whitespace();
                let chrom = {
                    let chrom = toks.next().unwrap();
                    if chrom.starts_with('#') || chrom == "track" {
                        continue;
                    }
                    chrom.to_string()
                };
                let start = toks.next().unwrap().parse::<Coordinate>().unwrap();
                let end = toks.next().unwrap().parse::<Coordinate>().unwrap();

                // optional fields
                let name = toks.next().map(|name| name.to_string());
                let score = if self.binarize_score {
                    toks.next();
                    Some(D::one())
                } else {
                    toks.next().map(|score| score.parse::<D>().unwrap())
                };
                let strand = toks.next().and_then(|strand| {
                    Strand::new(strand)
                        .expect("failed to parse the strand symbol")
                });
                Some(BedDataLine {
                    chrom,
                    start,
                    end,
                    name,
                    score,
                    strand,
                })
            };
        }
    }
}

pub struct BedCoordinateIter {
    buf: BufReader<File>,
    filename: String,
}

impl BedCoordinateIter {
    pub fn get_filename(&self) -> &str {
        &self.filename
    }
}

impl Iterator for BedCoordinateIter {
    type Item = (Chrom, Coordinate, Coordinate);

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let mut line = String::new();
            return if self.buf.read_line(&mut line).unwrap() == 0 {
                None
            } else {
                let mut toks = line.split_whitespace();
                let chrom = {
                    let chrom = toks.next().unwrap();
                    if chrom.starts_with('#') || chrom == "track" {
                        continue;
                    }
                    chrom.to_string()
                };
                let start = toks.next().unwrap().parse::<Coordinate>().unwrap();
                let end = toks.next().unwrap().parse::<Coordinate>().unwrap();
                Some((chrom, start, end))
            };
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        bed::{Bed, Chrom, Coordinate},
        iter::{ChromIntervalValue, ToChromIntervalValueIter},
    };
    use math::{
        partition::integer_interval_map::IntegerIntervalMap,
        set::{
            contiguous_integer_set::ContiguousIntegerSet,
            ordered_integer_set::OrderedIntegerSet,
        },
    };
    use std::{
        collections::HashMap,
        io::{BufWriter, Write},
    };
    use tempfile::NamedTempFile;

    #[test]
    fn test_get_chrom_to_interval_to_val() {
        let file = NamedTempFile::new().unwrap();
        {
            let mut writer = BufWriter::new(&file);
            writer
                .write_fmt(format_args!(
                    "chr1 100 200 name_1 3.5\n\
                    chr1 150 250 name_2 2\n\
                    chr1 200 350 name_3 4.0\n\
                    chr3 1000 3000 name_4 -0.3\n\
                    chr1 400 450 name_5 -0.9\n\
                    chr3 2500 3000 name_6 0.3\n"
                ))
                .unwrap();
        }
        let bed = Bed::new(file.path().to_str().unwrap(), false);
        {
            let chrom_to_interval_to_val =
                ToChromIntervalValueIter::get_chrom_to_interval_to_val(
                    &bed, None,
                )
                .unwrap();
            {
                let mut expected_chr1 = IntegerIntervalMap::<f64>::new();
                expected_chr1
                    .aggregate(ContiguousIntegerSet::new(100, 149), 3.5);
                expected_chr1
                    .aggregate(ContiguousIntegerSet::new(150, 199), 5.5);
                expected_chr1
                    .aggregate(ContiguousIntegerSet::new(200, 249), 6.);
                expected_chr1
                    .aggregate(ContiguousIntegerSet::new(250, 349), 4.);
                expected_chr1
                    .aggregate(ContiguousIntegerSet::new(400, 449), -0.9);
                assert_eq!(chrom_to_interval_to_val["chr1"], expected_chr1);
            }
            {
                let mut expected_chr3 = IntegerIntervalMap::<f64>::new();
                expected_chr3
                    .aggregate(ContiguousIntegerSet::new(1000, 2499), -0.3);
                expected_chr3
                    .aggregate(ContiguousIntegerSet::new(2500, 2999), 0.);
                assert_eq!(chrom_to_interval_to_val["chr3"], expected_chr3);
            }
        }

        {
            let exclude: HashMap<Chrom, OrderedIntegerSet<Coordinate>> = [
                (
                    "chr1".into(),
                    OrderedIntegerSet::from_slice(&[[80, 100], [190, 220]]),
                ),
                (
                    "chr3".into(),
                    OrderedIntegerSet::from_slice(&[[1010, 1020]]),
                ),
            ]
            .iter()
            .cloned()
            .collect();

            let chrom_to_interval_to_val =
                ToChromIntervalValueIter::get_chrom_to_interval_to_val(
                    &bed,
                    Some(exclude).as_ref(),
                )
                .unwrap();

            {
                let mut expected_chr1 = IntegerIntervalMap::<f64>::new();
                expected_chr1
                    .aggregate(ContiguousIntegerSet::new(400, 449), -0.9);
                assert_eq!(chrom_to_interval_to_val["chr1"], expected_chr1);
            }

            {
                let mut expected_chr3 = IntegerIntervalMap::<f64>::new();
                expected_chr3
                    .aggregate(ContiguousIntegerSet::new(2500, 2999), 0.3);
                assert_eq!(chrom_to_interval_to_val["chr3"], expected_chr3);
            }
        }
    }

    #[test]
    fn test_get_chrom_to_intervals() {
        let file = NamedTempFile::new().unwrap();
        {
            let mut writer = BufWriter::new(&file);
            writer
                .write_fmt(format_args!(
                    "chr1 100 200 name_1 3.5\n\
                    chr1 150 250 name_2 2\n\
                    chr1 200 350 name_3 4.0\n\
                    chr3 1000 3000 name_4 -0.3\n\
                    chr1 400 450 name_5 -0.9\n\
                    chr3 2500 3000 name_6 0.3\n"
                ))
                .unwrap();
        }
        let bed = Bed::new(file.path().to_str().unwrap(), false);
        let chrom_to_intervals = bed.get_chrom_to_intervals();

        let expected: HashMap<Chrom, OrderedIntegerSet<Coordinate>> = [
            (
                "chr1".into(),
                OrderedIntegerSet::from_slice(&[[100, 349], [400, 449]]),
            ),
            (
                "chr3".into(),
                OrderedIntegerSet::from_slice(&[[1000, 2999]]),
            ),
        ]
        .iter()
        .cloned()
        .collect();

        assert_eq!(chrom_to_intervals, expected);
    }

    #[test]
    fn test_bed_to_chrom_interval_value_iter() {
        let file = NamedTempFile::new().unwrap();
        {
            let mut writer = BufWriter::new(&file);
            writer
                .write_fmt(format_args!(
                    "chr1 100 200 name_1 3.5\n\
                    chr1 150 250 name_2 2\n\
                    chr1 200 350 name_3 4.0\n\
                    chr3 1000 3000 name_4 -0.3\n\
                    chr1 400 450 name_5 -0.9\n\
                    chr3 2500 3000 name_6 0.3\n"
                ))
                .unwrap();
        }
        let bed = Bed::new(file.path().to_str().unwrap(), false);
        let mut iter = bed.to_chrom_interval_value_iter();
        assert_eq!(
            iter.next().unwrap().chrom_interval_value(),
            ("chr1".to_string(), ContiguousIntegerSet::new(100, 199), 3.5)
        );
        assert_eq!(
            iter.next().unwrap().chrom_interval_value(),
            ("chr1".to_string(), ContiguousIntegerSet::new(150, 249), 2.)
        );
        assert_eq!(
            iter.next().unwrap().chrom_interval_value(),
            ("chr1".to_string(), ContiguousIntegerSet::new(200, 349), 4.)
        );
        assert_eq!(
            iter.next().unwrap().chrom_interval_value(),
            (
                "chr3".to_string(),
                ContiguousIntegerSet::new(1000, 2999),
                -0.3
            )
        );
        assert_eq!(
            iter.next().unwrap().chrom_interval_value(),
            (
                "chr1".to_string(),
                ContiguousIntegerSet::new(400, 449),
                -0.9
            )
        );
        assert_eq!(
            iter.next().unwrap().chrom_interval_value(),
            (
                "chr3".to_string(),
                ContiguousIntegerSet::new(2500, 2999),
                0.3
            )
        );
        assert_eq!(iter.next(), None);
    }
}
