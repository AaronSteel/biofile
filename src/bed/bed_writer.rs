use crate::{bed::BedDataLine, bedgraph::BedGraphDataLine, util::Strand};
use num::Float;
use std::{
    fs::{File, OpenOptions},
    io::{BufWriter, Write},
};

pub struct BedWriter {
    path: String,
    writer: BufWriter<File>,
    num_lines_written: i64,
}

impl BedWriter {
    pub fn new(path: &str) -> Result<Self, std::io::Error> {
        let file = OpenOptions::new()
            .create(true)
            .truncate(true)
            .write(true)
            .open(path)?;
        let path_buf = std::fs::canonicalize(path)?;
        Ok(BedWriter {
            path: path_buf.to_str().unwrap().to_string(),
            writer: BufWriter::new(file),
            num_lines_written: 0,
        })
    }

    pub fn canonical_path(&self) -> &str {
        &self.path
    }

    pub fn write_bed_line<D: Float + std::fmt::Display>(
        &mut self,
        data: &BedDataLine<D>,
    ) -> std::io::Result<usize> {
        let mut line = format!("{}\t{}\t{}", data.chrom, data.start, data.end);
        match &data.name {
            None => line.push_str(&format!("\tid_{}", self.num_lines_written)),
            Some(name) => line.push_str(&format!("\t{}", name)),
        };
        match data.score {
            None => line.push_str("\t0"),
            Some(score) => line.push_str(&format!("\t{}", score)),
        }
        match data.strand {
            None => line.push_str("\t."),
            Some(strand) => match strand {
                Strand::Positive => {
                    line.push_str("\t+");
                }
                Strand::Negative => {
                    line.push_str("\t-");
                }
            },
        }
        line.push_str("\n");
        self.num_lines_written += 1;
        self.writer.write(line.as_bytes())
    }

    pub fn write_bedgraph_line<D: Float + std::fmt::Display>(
        &mut self,
        data: &BedGraphDataLine<D>,
    ) -> std::io::Result<usize> {
        let line = format!(
            "{}\t{}\t{}\t{}\n",
            data.chrom, data.start, data.end_exclusive, data.value
        );
        self.num_lines_written += 1;
        self.writer.write(line.as_bytes())
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        bed::{bed_writer::BedWriter, Bed, BedDataLine},
        bedgraph::{BedGraph, BedGraphDataLine, BedGraphDataLineIter},
        util::Strand,
    };
    use math::traits::ToIterator;
    use num::Num;
    use tempfile::NamedTempFile;

    #[test]
    fn test_bed_writer() {
        let data_lines = vec![
            BedDataLine {
                chrom: "chr1".into(),
                start: 0,
                end: 20,
                name: None,
                score: Some(9.),
                strand: None,
            },
            BedDataLine {
                chrom: "chr1".into(),
                start: 10,
                end: 50,
                name: Some("some_id".into()),
                score: None::<f64>,
                strand: Some(Strand::Positive),
            },
            BedDataLine {
                chrom: "chr2".into(),
                start: 0,
                end: 100,
                name: None,
                score: Some(10.),
                strand: Some(Strand::Negative),
            },
        ];
        let file = NamedTempFile::new().unwrap();
        let path = file.into_temp_path();
        {
            let mut writer = BedWriter::new(path.to_str().unwrap()).unwrap();
            for data in data_lines.iter() {
                writer.write_bed_line(data).unwrap();
            }
        }
        {
            let bed = Bed::new(path.to_str().unwrap(), false);
            for (line, expected) in bed
                .to_iter()
                .zip(get_expected_bed_data_lines(&data_lines).iter())
            {
                assert_eq!(&line, expected);
            }
        }
        path.close().unwrap();
    }

    #[test]
    fn test_write_bedgraph() {
        let data_lines = vec![
            BedGraphDataLine::<f32> {
                chrom: "chr1".into(),
                start: 3,
                end_exclusive: 10,
                value: 5.,
            },
            BedGraphDataLine {
                chrom: "chr1".into(),
                start: 0,
                end_exclusive: 23,
                value: 6.,
            },
            BedGraphDataLine {
                chrom: "chr3".into(),
                start: 100,
                end_exclusive: 115,
                value: 11.,
            },
            BedGraphDataLine {
                chrom: "chr1".into(),
                start: 100,
                end_exclusive: 150,
                value: 18.,
            },
        ];
        let file = NamedTempFile::new().unwrap();
        let path = file.into_temp_path();
        {
            let mut writer = BedWriter::new(path.to_str().unwrap()).unwrap();
            for data in data_lines.iter() {
                writer.write_bedgraph_line(data).unwrap();
            }
        }
        {
            let bedgraph = BedGraph::new(path.to_str().unwrap());
            for (expected, line) in data_lines
                .iter()
                .zip(bedgraph.to_iter(): BedGraphDataLineIter<f32>)
            {
                assert_eq!(expected, &line);
            }
        }
    }

    fn get_expected_bed_data_lines<D: Copy + Num>(
        original_data_lines: &[BedDataLine<D>],
    ) -> Vec<BedDataLine<D>> {
        original_data_lines
            .iter()
            .enumerate()
            .map(|(i, l)| BedDataLine {
                chrom: l.chrom.clone(),
                start: l.start,
                end: l.end,
                name: Some(
                    l.name.clone().unwrap_or_else(|| format!("id_{}", i)),
                ),
                score: Some(l.score.unwrap_or_else(D::zero)),
                strand: l.strand,
            })
            .collect()
    }
}
