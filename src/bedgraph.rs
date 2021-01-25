//! An interface to the BedGraph track format file as specified in
//! https://genome.ucsc.edu/goldenPath/help/bedgraph.html

use math::traits::ToIterator;
use num::Float;
use std::{
    fmt::Debug,
    fs::File,
    io::{BufRead, BufReader},
    marker::PhantomData,
    str::FromStr,
};

use crate::{
    iter::{ChromIntervalValue, ToChromIntervalValueIter},
    util::get_buf,
};
use math::set::contiguous_integer_set::ContiguousIntegerSet;

pub struct BedGraph {
    filepath: String,

    /// Every line in the bedgraph file will contribute a unit score for the
    /// corresponding interval.
    binarize_score: bool,
}

impl BedGraph {
    pub fn new(filepath: &str, binarize_score: bool) -> BedGraph {
        BedGraph {
            filepath: filepath.to_string(),
            binarize_score,
        }
    }

    #[inline]
    pub fn get_filepath(&self) -> &str {
        &self.filepath
    }
}

impl<D, E>
    ToIterator<
        '_,
        BedGraphDataLineIter<D>,
        <BedGraphDataLineIter<D> as Iterator>::Item,
    > for BedGraph
where
    D: Float + FromStr<Err = E>,
    E: Debug,
{
    fn to_iter(&self) -> BedGraphDataLineIter<D> {
        let buf = get_buf(&self.filepath).unwrap();
        BedGraphDataLineIter {
            buf,
            filename: self.filepath.clone(),
            binarize_score: self.binarize_score,
            phantom: PhantomData,
        }
    }
}

impl<V: Float + FromStr<Err = E>, E: Debug>
    ToChromIntervalValueIter<
        BedGraphDataLineIter<V>,
        BedGraphDataLine<V>,
        Coordinate,
        V,
    > for BedGraph
{
    fn to_chrom_interval_value_iter(&self) -> BedGraphDataLineIter<V> {
        self.to_iter()
    }
}

impl<V: Copy + Float + FromStr<Err = E>, E: Debug>
    ChromIntervalValue<Coordinate, V> for BedGraphDataLine<V>
{
    fn chrom_interval_value(
        &self,
    ) -> (Chrom, ContiguousIntegerSet<Coordinate>, V) {
        (
            self.chrom.clone(),
            ContiguousIntegerSet::new(self.start, self.end_exclusive - 1),
            self.value,
        )
    }
}

/// Data type of the Bedgraph coordinates
pub type Coordinate = i64;

/// Data type of the Bedgraph chromosomes
pub type Chrom = String;

/// The four-element tuple in the `BedGraphDataLine` corresponds to a line of
/// data in the BedGraph file, where each line is of the form
/// `chrom start end value`
///
/// The [start, end) is a zero-based left-closed right-open coordinate range
#[derive(Clone, PartialEq, Debug)]
pub struct BedGraphDataLine<D> {
    pub chrom: Chrom,
    pub start: Coordinate,
    pub end_exclusive: Coordinate,
    pub value: D,
}

pub struct BedGraphDataLineIter<D> {
    buf: BufReader<File>,
    filename: String,

    /// Every line in the bedgraph file will contribute a unit score for the
    /// corresponding interval.
    binarize_score: bool,
    phantom: PhantomData<D>,
}

impl<D> BedGraphDataLineIter<D> {
    pub fn get_filename(&self) -> &str {
        &self.filename
    }
}

impl<D: Float + FromStr<Err = E>, E: Debug> Iterator
    for BedGraphDataLineIter<D>
{
    type Item = BedGraphDataLine<D>;

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
                let end_exclusive =
                    toks.next().unwrap().parse::<Coordinate>().unwrap();

                let value = if self.binarize_score {
                    D::one()
                } else {
                    toks.next().unwrap().parse::<D>().unwrap()
                };

                Some(BedGraphDataLine {
                    chrom,
                    start,
                    end_exclusive,
                    value,
                })
            };
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        bedgraph::{BedGraph, BedGraphDataLineIter},
        iter::{ChromIntervalValue, ToChromIntervalValueIter},
    };
    use math::{
        partition::integer_interval_map::IntegerIntervalMap,
        set::contiguous_integer_set::ContiguousIntegerSet,
    };
    use std::io::{BufWriter, Write};
    use tempfile::NamedTempFile;

    #[test]
    fn test_get_chrom_to_interval_to_val() {
        let file = NamedTempFile::new().unwrap();
        {
            let mut writer = BufWriter::new(&file);
            writer
                .write_fmt(format_args!(
                    "chr1 100 200 3.5\n\
                    chr1 150 250 2\n\
                    chr1 200 350 4.0\n\
                    chr3 1000 3000 -0.3\n\
                    chr1 400 450 -0.9\n\
                    chr3 2500 3000 0.3\n"
                ))
                .unwrap();
        }
        let bedgraph = BedGraph::new(file.path().to_str().unwrap(), false);
        let chrom_to_interval_to_val =
            ToChromIntervalValueIter::get_chrom_to_interval_to_val(
                &bedgraph, None,
            )
            .unwrap();

        let mut expected_chr1 = IntegerIntervalMap::<f64>::new();
        expected_chr1.aggregate(ContiguousIntegerSet::new(100, 149), 3.5);
        expected_chr1.aggregate(ContiguousIntegerSet::new(150, 199), 5.5);
        expected_chr1.aggregate(ContiguousIntegerSet::new(200, 249), 6.);
        expected_chr1.aggregate(ContiguousIntegerSet::new(250, 349), 4.);
        expected_chr1.aggregate(ContiguousIntegerSet::new(400, 449), -0.9);

        let mut expected_chr3 = IntegerIntervalMap::<f64>::new();
        expected_chr3.aggregate(ContiguousIntegerSet::new(1000, 2499), -0.3);
        expected_chr3.aggregate(ContiguousIntegerSet::new(2500, 2999), 0.);
        assert_eq!(chrom_to_interval_to_val["chr1"], expected_chr1);
        assert_eq!(chrom_to_interval_to_val["chr3"], expected_chr3);
    }

    #[test]
    fn test_bedgraph_to_chrom_interval_value_iter() {
        let file = NamedTempFile::new().unwrap();
        {
            let mut writer = BufWriter::new(&file);
            writer
                .write_fmt(format_args!(
                    "chr1 100 200 3.5\n\
                    chr1 150 250 2\n\
                    chr1 200 350 4.0\n\
                    chr3 1000 3000 -0.3\n\
                    chr1 400 450 -0.9\n\
                    chr3 2500 3000 0.3\n"
                ))
                .unwrap();
        }
        let bedgraph = BedGraph::new(file.path().to_str().unwrap(), false);
        let mut iter: BedGraphDataLineIter<f64> =
            bedgraph.to_chrom_interval_value_iter();

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

    // TODO: test binarize_score
}
