use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use analytic::partition::ordered_interval_partitions::OrderedIntervalPartitions;
use analytic::set::ordered_integer_set::{ContiguousIntegerSet, OrderedIntegerSet};

use crate::error::Error;
use crate::util::get_buf;

pub struct PeakFile {
    filepath: String
}

impl PeakFile {
    pub fn new(filepath: String) -> PeakFile {
        PeakFile {
            filepath
        }
    }

    pub fn iter(&self) -> Result<PeakFileIter, Error> {
        Ok(PeakFileIter::new(get_buf(&self.filepath)?))
    }

    pub fn get_chrom_to_peaks(&self)
        -> Result<HashMap<String, OrderedIntervalPartitions<usize>>, Error> {
        let mut chrom_to_intervals = HashMap::new();
        for line in self.iter()? {
            let chrom = line.chrom;
            let start = line.start;
            let end = line.end;
            if end > start {
                let interval_to_val = chrom_to_intervals.entry(chrom).or_insert(Vec::new());
                interval_to_val.push(ContiguousIntegerSet::new(start, end - 1));
            }
        }
        let mut chrom_to_partitions = HashMap::new();
        for (chrom, intervals) in chrom_to_intervals.into_iter() {
            chrom_to_partitions.insert(
                chrom,
                OrderedIntervalPartitions::from_vec_with_trusted_order(
                    OrderedIntegerSet::from_contiguous_integer_sets(intervals).into_intervals()
                ),
            );
        }
        Ok(chrom_to_partitions)
    }

    pub fn get_filepath(&self) -> &str {
        &self.filepath
    }
}

#[derive(PartialEq, Copy, Clone, Debug)]
pub enum Strand {
    Positive,
    Negative,
}

/// The [start, end) is a zero-based left-closed right-open coordinate range
#[derive(PartialEq, Clone, Debug)]
pub struct PeakFileDataLine {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    pub name: String,
    pub score: f64,
    pub strand: Option<Strand>,
    pub signal_value: f64,
    pub p_value: f64,
    pub q_value: f64,
    pub peak: Option<usize>,
}

pub struct PeakFileIter {
    buf: BufReader<File>
}

impl PeakFileIter {
    pub fn new(buf: BufReader<File>) -> PeakFileIter {
        PeakFileIter {
            buf
        }
    }
}

impl Iterator for PeakFileIter {
    type Item = PeakFileDataLine;
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

                let name = toks.next().unwrap().to_string();
                let score = toks.next().unwrap().parse::<f64>().unwrap();
                let strand = match toks.next().unwrap() {
                    "+" => Some(Strand::Positive),
                    "-" => Some(Strand::Negative),
                    "." => None,
                    s => panic!(format!("unrecognized strand symbol {}", s)),
                };

                let signal_value = toks.next().unwrap().parse::<f64>().unwrap();
                let p_value = toks.next().unwrap().parse::<f64>().unwrap();
                let q_value = toks.next().unwrap().parse::<f64>().unwrap();
                let peak = match toks.next() {
                    Some(p) => Some(p.parse::<usize>().unwrap()),
                    None => None,
                };
                return Some(PeakFileDataLine {
                    chrom,
                    start,
                    end,
                    name,
                    score,
                    strand,
                    signal_value,
                    p_value,
                    q_value,
                    peak,
                });
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::io::{BufWriter, Write};

    use tempfile::NamedTempFile;
    use std::collections::HashMap;
    use crate::peak_file::{PeakFile, PeakFileDataLine};
    use analytic::partition::ordered_interval_partitions::OrderedIntervalPartitions;
    use analytic::set::ordered_integer_set::ContiguousIntegerSet;

    #[test]
    fn test_get_chrom_to_interval_to_val() {
        let file = NamedTempFile::new().unwrap();
        {
            let mut writer = BufWriter::new(&file);
            writer.write_fmt(
                format_args!(
                    "track type=narrowPeak name=\"H9_H3K4ME3_H3K79ME2.phased.1\" description=\"H9_H3K4ME3_H3K79ME2.phased.1\" nextItemButton=on\n\
                    chr1 10050 10500 H9_H3K4ME3_H3K79ME2.phased.1_narrowPeak1 153 . 0.00000 0.00000 0.00000 125\n\
                    chr1 28650 28900 H9_H3K4ME3_H3K79ME2.phased.1_narrowPeak2 96 . 0.00000 0.00000 0.00000 175\n\
                    chr1 29000 29250 H9_H3K4ME3_H3K79ME2.phased.1_narrowPeak3 153 . 0.00000 0.00000 0.00000 125\n"
                )
            ).unwrap();
        }
        let peak_file = PeakFile::new(file.path().to_str().unwrap().to_string());
        let mut iter = peak_file.iter().unwrap();
        assert_eq!(
            Some(PeakFileDataLine {
                chrom: "chr1".to_string(),
                start: 10050,
                end: 10500,
                name: "H9_H3K4ME3_H3K79ME2.phased.1_narrowPeak1".to_string(),
                score: 153.,
                strand: None,
                signal_value: 0.,
                p_value: 0.,
                q_value: 0.,
                peak: Some(125),
            }),
            iter.next()
        );

        assert_eq!(
            Some(PeakFileDataLine {
                chrom: "chr1".to_string(),
                start: 28650,
                end: 28900,
                name: "H9_H3K4ME3_H3K79ME2.phased.1_narrowPeak2".to_string(),
                score: 96.,
                strand: None,
                signal_value: 0.,
                p_value: 0.,
                q_value: 0.,
                peak: Some(175),
            }),
            iter.next()
        );

        assert_eq!(
            Some(PeakFileDataLine {
                chrom: "chr1".to_string(),
                start: 29000,
                end: 29250,
                name: "H9_H3K4ME3_H3K79ME2.phased.1_narrowPeak3".to_string(),
                score: 153.,
                strand: None,
                signal_value: 0.,
                p_value: 0.,
                q_value: 0.,
                peak: Some(125),
            }),
            iter.next()
        );
    }

    #[test]
    fn test_get_chrom_to_peaks() {
        let file = NamedTempFile::new().unwrap();
        {
            let mut writer = BufWriter::new(&file);
            writer.write_fmt(
                format_args!(
                    "track type=narrowPeak name=\"H9_H3K4ME3_H3K79ME2.phased.1\" description=\"H9_H3K4ME3_H3K79ME2.phased.1\" nextItemButton=on\n\
                    chr1 10050 10500 H9_H3K4ME3_H3K79ME2.phased.1_narrowPeak1 153 . 0.00000 0.00000 0.00000 125\n\
                    chr1 28650 28900 H9_H3K4ME3_H3K79ME2.phased.1_narrowPeak2 96 . 0.00000 0.00000 0.00000 175\n\
                    chr1 1000 2000 H9_H3K4ME3_H3K79ME2.phased.1_narrowPeak3 153 . 0.00000 0.00000 0.00000 125\n\
                    chr2 0 2000 H9_H3K4ME3_H3K79ME2.phased.1_narrowPeak4 153 . 0.00000 0.00000 0.00000 125\n\
                    chr2 100 2500 H9_H3K4ME3_H3K79ME2.phased.1_narrowPeak5 153 . 0.00000 0.00000 0.00000 125\n\
                    chr2 3000 10000 H9_H3K4ME3_H3K79ME2.phased.1_narrowPeak6 153 . 0.00000 0.00000 0.00000 125\n\
                    chr3 0 1000 H9_H3K4ME3_H3K79ME2.phased.1_narrowPeak7 153 . 0.00000 0.00000 0.00000 125\n\
                    chr4 10050 10500 H9_H3K4ME3_H3K79ME2.phased.1_narrowPeak8 153 . 0.00000 0.00000 0.00000 125\n\
                    chr4 28650 28900 H9_H3K4ME3_H3K79ME2.phased.1_narrowPeak9 153 . 0.00000 0.00000 0.00000 125\n"
                )
            ).unwrap();
        }
        let peak_file = PeakFile::new(file.path().to_str().unwrap().to_string());
        let expected: HashMap<String, OrderedIntervalPartitions<usize>> = vec![
            (
                "chr1",
                OrderedIntervalPartitions::from_vec_with_trusted_order(
                    vec![
                        ContiguousIntegerSet::new(1000usize, 1999),
                        ContiguousIntegerSet::new(10050, 10499),
                        ContiguousIntegerSet::new(28650, 28899)
                    ]
                )
            ),
            (
                "chr2",
                OrderedIntervalPartitions::from_vec_with_trusted_order(
                    vec![
                        ContiguousIntegerSet::new(0usize, 2499),
                        ContiguousIntegerSet::new(3000, 9999)
                    ]
                )
            ),
            (
                "chr3",
                OrderedIntervalPartitions::from_vec_with_trusted_order(
                    vec![ContiguousIntegerSet::new(0usize, 999)]
                )
            ),
            (
                "chr4",
                OrderedIntervalPartitions::from_vec_with_trusted_order(
                    vec![
                        ContiguousIntegerSet::new(10050usize, 10499),
                        ContiguousIntegerSet::new(28650, 28899)
                    ]
                )
            ),
        ].into_iter()
         .map(|(c, p)| (c.to_string(), p))
         .collect();
        assert_eq!(expected, peak_file.get_chrom_to_peaks().unwrap());
    }
}
