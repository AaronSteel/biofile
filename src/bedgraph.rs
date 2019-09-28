//! An interface to the BedGraph track format file as specified in
//! https://genome.ucsc.edu/goldenPath/help/bedgraph.html

use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::marker::PhantomData;
use std::str::FromStr;

use analytic::set::ordered_integer_set::ContiguousIntegerSet;
use analytic::traits::ToIterator;
use num::Float;

use crate::error::Error;
use crate::util::get_buf;

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

    pub fn get_chrom_to_interval_to_val<D, E>(&self)
        -> HashMap<String, HashMap<ContiguousIntegerSet<usize>, D>>
        where D: Float + FromStr<Err=E>, E: fmt::Debug {
        let mut chrom_to_interval_to_val = HashMap::new();
        for (chrom, start, end, val) in self.to_iter(): BedGraphDataLineIter<D> {
            let interval_to_val = chrom_to_interval_to_val.entry(chrom).or_insert(HashMap::new());
            interval_to_val.insert(ContiguousIntegerSet::new(start, end - 1), val);
        }
        chrom_to_interval_to_val
    }
}

impl<D, E> ToIterator<'_, BedGraphDataLineIter<D>, <BedGraphDataLineIter<D> as Iterator>::Item> for BedGraph
    where D: Float + FromStr<Err=E>, E: fmt::Debug {
    fn to_iter(&self) -> BedGraphDataLineIter<D> {
        let buf = get_buf(&self.filepath).unwrap();
        BedGraphDataLineIter {
            buf,
            filename: self.filepath.clone(),
            phantom: PhantomData,
        }
    }
}

/// The four-element tuple in the `BedGraphDataLine` corresponds to a line of data in the BedGraph file,
/// where each line is of the form
/// `chrom start end value`
///
/// The [start, end) is a zero-based left-closed right-open coordinate range
pub type BedGraphDataLine<D> = (String, usize, usize, D);

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

impl<D: Float + FromStr<Err=E>, E: fmt::Debug> Iterator for BedGraphDataLineIter<D> {
    type Item = BedGraphDataLine<D>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let mut line = String::new();
            if self.buf.read_line(&mut line).unwrap() == 0 {
                return None;
            } else {
                let mut toks = line.split_whitespace();
                let chrom = {
                    let chrom = toks.next().unwrap();
                    if chrom.starts_with("#") || chrom == "track" {
                        continue;
                    }
                    chrom.to_string()
                };
                let start = toks.next().unwrap().parse::<usize>().unwrap();
                let end = toks.next().unwrap().parse::<usize>().unwrap();
                let value = toks.next().unwrap().parse::<D>().unwrap();
                return Some((chrom, start, end, value));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::io::{BufWriter, Write};

    use analytic::set::ordered_integer_set::ContiguousIntegerSet;
    use tempfile::NamedTempFile;

    use crate::bedgraph::BedGraph;

    #[test]
    fn test_get_chrom_to_interval_to_val() {
        let file = NamedTempFile::new().unwrap();
        {
            let mut writer = BufWriter::new(&file);
            writer.write_fmt(
                format_args!(
                    "chr1 100 200 3.5\n\
                    chr1 200 350 4.0\n\
                    chr3 1000 3000 -0.3\n\
                    chr1 400 450 -0.9"
                )
            ).unwrap();
        }
        let bedgraph = BedGraph::new(file.path().to_str().unwrap()).unwrap();
        let chrom_to_interval_to_val = bedgraph.get_chrom_to_interval_to_val();
        let expected_chr1: HashMap<ContiguousIntegerSet<usize>, f64> = [
            (ContiguousIntegerSet::new(100, 199usize), 3.5),
            (ContiguousIntegerSet::new(200, 349), 4.),
            (ContiguousIntegerSet::new(400, 449), -0.9),
        ].into_iter().map(|(a, b)| (*a, *b)).collect();
        let expected_chr3: HashMap<ContiguousIntegerSet<usize>, f64> = [
            (ContiguousIntegerSet::new(1000, 2999usize), -0.3)
        ].into_iter().map(|(a, b)| (*a, *b)).collect();
        assert_eq!(
            chrom_to_interval_to_val["chr1"],
            expected_chr1
        );
        assert_eq!(
            chrom_to_interval_to_val["chr3"],
            expected_chr3
        );
    }
}
