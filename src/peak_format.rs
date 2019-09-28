use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::error::Error;
use crate::util::get_buf;

pub struct PeakFormat {
    filepath: String
}

impl PeakFormat {
    pub fn new(filepath: String) -> PeakFormat {
        PeakFormat {
            filepath
        }
    }

    pub fn iter(&self) -> Result<PeakFormatIter, Error> {
        Ok(PeakFormatIter::new(get_buf(&self.filepath)?))
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

/// Browser Extensible Data (BED) with 6 fields
/// The [start, end) is a zero-based left-closed right-open coordinate range
#[derive(PartialEq, Clone, Debug)]
pub struct Bed6Line {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    pub name: String,
    pub score: f64,
    pub strand: Option<Strand>,
}

#[derive(PartialEq, Clone, Debug)]
pub struct PeakFormatDataLine {
    pub bed_6: Bed6Line,
    pub signal_value: f64,
    pub p_value: f64,
    pub q_value: f64,
    pub peak: Option<usize>,
}

pub struct PeakFormatIter {
    buf: BufReader<File>
}

impl PeakFormatIter {
    pub fn new(buf: BufReader<File>) -> PeakFormatIter {
        PeakFormatIter {
            buf
        }
    }
}

impl Iterator for PeakFormatIter {
    type Item = PeakFormatDataLine;
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
                return Some(PeakFormatDataLine {
                    bed_6: Bed6Line {
                        chrom,
                        start,
                        end,
                        name,
                        score,
                        strand,
                    },
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

    use crate::peak_format::{Bed6Line, PeakFormat, PeakFormatDataLine};

    #[test]
    fn test_get_chrom_to_interval_to_val() {
        let file = NamedTempFile::new().unwrap();
        {
            let mut writer = BufWriter::new(&file);
            writer.write_fmt(
                format_args!(
                    "track type=narrowPeak name=\"H9_H3K4ME3_H3K79ME2.phased.1\" description=\"H9_H3K4ME3_H3K79ME2.phased.1\" nextItemButton=on\n\
                    chr1	10050	10500	H9_H3K4ME3_H3K79ME2.phased.1_narrowPeak1	153	.	0.00000	0.00000	0.00000	125\n\
                    chr1	28650	28900	H9_H3K4ME3_H3K79ME2.phased.1_narrowPeak2	96	.	0.00000	0.00000	0.00000	175\n\
                    chr1	29000	29250	H9_H3K4ME3_H3K79ME2.phased.1_narrowPeak3	153	.	0.00000	0.00000	0.00000	125\n"
                )
            ).unwrap();
        }
        let peak_format = PeakFormat::new(file.path().to_str().unwrap().to_string());
        let mut iter = peak_format.iter().unwrap();
        assert_eq!(
            Some(PeakFormatDataLine {
                bed_6: Bed6Line {
                    chrom: "chr1".to_string(),
                    start: 10050,
                    end: 10500,
                    name: "H9_H3K4ME3_H3K79ME2.phased.1_narrowPeak1".to_string(),
                    score: 153.,
                    strand: None,
                },
                signal_value: 0.,
                p_value: 0.,
                q_value: 0.,
                peak: Some(125),
            }),
            iter.next()
        );

        assert_eq!(
            Some(PeakFormatDataLine {
                bed_6: Bed6Line {
                    chrom: "chr1".to_string(),
                    start: 28650,
                    end: 28900,
                    name: "H9_H3K4ME3_H3K79ME2.phased.1_narrowPeak2".to_string(),
                    score: 96.,
                    strand: None,
                },
                signal_value: 0.,
                p_value: 0.,
                q_value: 0.,
                peak: Some(175),
            }),
            iter.next()
        );

        assert_eq!(
            Some(PeakFormatDataLine {
                bed_6: Bed6Line {
                    chrom: "chr1".to_string(),
                    start: 29000,
                    end: 29250,
                    name: "H9_H3K4ME3_H3K79ME2.phased.1_narrowPeak3".to_string(),
                    score: 153.,
                    strand: None,
                },
                signal_value: 0.,
                p_value: 0.,
                q_value: 0.,
                peak: Some(125),
            }),
            iter.next()
        );
    }
}
