//! An interface to the BedGraph track format file as specified in
//! https://genome.ucsc.edu/goldenPath/help/bedgraph.html

use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::marker::PhantomData;

use num::Float;

use crate::error::Error;
use crate::util::get_buf;
use std::str::FromStr;
use analytic::traits::ToIterator;

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
/// chrom start end value
///
/// The [start, end) is a zero-based left-closed right-open coordinate range
#[derive(Clone, Debug)]
pub struct BedGraphDataLine<D: Float>(pub String, pub usize, pub usize, pub D);

pub struct BedGraphDataLineIter<D> {
    buf: BufReader<File>,
    filename: String,
    phantom: PhantomData<D>,
}

impl<D: Float + FromStr<Err=E>, E: fmt::Debug> Iterator for BedGraphDataLineIter<D> {
    type Item = BedGraphDataLine<D>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut line = String::new();
        loop {
            line.clear();
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
                return Some(BedGraphDataLine(chrom, start, end, value));
            }
        }
    }
}
