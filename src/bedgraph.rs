//! An interface to the BedGraph track format file as specified in
//! https://genome.ucsc.edu/goldenPath/help/bedgraph.html

use analytic::{
    partition::integer_interval_map::IntegerIntervalMap,
    set::contiguous_integer_set::ContiguousIntegerSet, traits::ToIterator,
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

use crate::{error::Error, util::get_buf};

pub struct BedGraph {
    filepath: String,
}

impl BedGraph {
    pub fn new(filepath: &str) -> Result<BedGraph, Error> {
        Ok(BedGraph {
            filepath: filepath.to_string(),
        })
    }

    #[inline]
    pub fn get_filepath(&self) -> &str {
        &self.filepath
    }

    pub fn get_chrom_to_interval_to_val<D, E>(&self) -> HashMap<String, IntegerIntervalMap<D>>
    where
        D: Float + FromStr<Err = E>,
        E: Debug, {
        let mut chrom_to_interval_to_val = HashMap::new();
        for (chrom, start, end, val) in self.to_iter(): BedGraphDataLineIter<D> {
            let interval_to_val = chrom_to_interval_to_val
                .entry(chrom)
                .or_insert_with(IntegerIntervalMap::new);
            interval_to_val.aggregate(ContiguousIntegerSet::new(start, end - 1), val);
        }
        chrom_to_interval_to_val
    }
}

impl<D, E> ToIterator<'_, BedGraphDataLineIter<D>, <BedGraphDataLineIter<D> as Iterator>::Item>
    for BedGraph
where
    D: Float + FromStr<Err = E>,
    E: Debug,
{
    fn to_iter(&self) -> BedGraphDataLineIter<D> {
        let buf = get_buf(&self.filepath).unwrap();
        BedGraphDataLineIter {
            buf,
            filename: self.filepath.clone(),
            phantom: PhantomData,
        }
    }
}

/// Data type of the Bedgraph coordinates
pub type Coordinate = i64;

/// The four-element tuple in the `BedGraphDataLine` corresponds to a line of data in the BedGraph
/// file, where each line is of the form
/// `chrom start end value`
///
/// The [start, end) is a zero-based left-closed right-open coordinate range
pub type BedGraphDataLine<D> = (String, Coordinate, Coordinate, D);

pub struct BedGraphDataLineIter<D> {
    buf: BufReader<File>,
    filename: String,
    phantom: PhantomData<D>,
}

impl<D> BedGraphDataLineIter<D> {
    pub fn get_filename(&self) -> &str {
        &self.filename
    }
}

impl<D: Float + FromStr<Err = E>, E: Debug> Iterator for BedGraphDataLineIter<D> {
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
                let end = toks.next().unwrap().parse::<Coordinate>().unwrap();
                let value = toks.next().unwrap().parse::<D>().unwrap();
                Some((chrom, start, end, value))
            };
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::bedgraph::BedGraph;
    use analytic::{
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
        let bedgraph = BedGraph::new(file.path().to_str().unwrap()).unwrap();
        let chrom_to_interval_to_val = bedgraph.get_chrom_to_interval_to_val();

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
}
