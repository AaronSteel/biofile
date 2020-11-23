use crate::{
    bed::{
        bed_writer::BedWriter, Bed, BedDataLine, BedDataLineIter, Coordinate,
    },
    error::Error,
    util::Strand,
};
use math::{
    interval::traits::Interval,
    set::{contiguous_integer_set::ContiguousIntegerSet, traits::Set},
    traits::ToIterator,
};
use std::{
    collections::{HashMap, HashSet},
    fmt::Debug,
};

pub struct PairedEndCollator {
    path: String,
}

impl PairedEndCollator {
    pub fn new(path: &str) -> PairedEndCollator {
        PairedEndCollator {
            path: path.to_string(),
        }
    }

    pub fn write(
        &self,
        out_path: &str,
        left_offset: Coordinate,
        right_offset: Coordinate,
    ) -> Result<PairedEndCollatorResult, Error>
where {
        let mut visited = HashSet::new();
        let mut chrom_to_id_to_pair = HashMap::new();
        let mut num_pcr_duplicates = 0i64;

        let bed = Bed::new(&self.path);
        for BedDataLine {
            chrom,
            start,
            end,
            name,
            score: _,
            strand,
        } in bed.to_iter(): BedDataLineIter<f64>
        {
            // only process reads with either a + or a - strand
            let strand = match strand {
                None => continue,
                Some(s) => s,
            };
            if !visited.insert((chrom.clone(), start, end, strand)) {
                // duplicate PCR reads
                num_pcr_duplicates += 1;
                continue;
            }
            match name {
                None => {
                    // the paired-end info is assumed to be contained in the
                    // name field
                    continue;
                }
                Some(name) => match PairedEndCollator::get_pair_id(&name) {
                    None => continue,
                    Some((pair_id, _index)) => {
                        let id_to_pair = chrom_to_id_to_pair
                            .entry(chrom.to_string())
                            .or_insert_with(|| HashMap::new());

                        let pair = id_to_pair
                            .entry(pair_id.clone())
                            .or_insert_with(|| PairedEndReadPair::new());

                        match strand {
                            Strand::Positive => match pair.postive_strand {
                                None => {
                                    pair.postive_strand = Some(
                                        ContiguousIntegerSet::new(start, end),
                                    );
                                }
                                Some(old_strand) => {
                                    return Err(Error::Generic(format!(
                                        "conflicting read with id {} on the \
                                        positive strand with coordinates \
                                        ({}, {}) and ({}, {})",
                                        pair_id,
                                        old_strand.get_start(),
                                        old_strand.get_end(),
                                        start,
                                        end
                                    )));
                                }
                            },
                            Strand::Negative => match pair.negative_strand {
                                None => {
                                    pair.negative_strand = Some(
                                        ContiguousIntegerSet::new(start, end),
                                    );
                                }
                                Some(old_strand) => {
                                    return Err(Error::Generic(format!(
                                        "conflicting read with id {} on the \
                                        negative strand with coordinates \
                                        ({}, {}) and ({}, {})",
                                        pair_id,
                                        old_strand.get_start(),
                                        old_strand.get_end(),
                                        start,
                                        end
                                    )));
                                }
                            },
                        };
                    }
                },
            }
        }

        let chrom_to_sorted_id_pairs: HashMap<
            String,
            Vec<(String, PairedEndReadPair)>,
        > = chrom_to_id_to_pair
            .into_iter()
            .map(|(chrom, id_to_pair)| {
                let mut id_pair_list: Vec<(String, PairedEndReadPair)> =
                    id_to_pair
                        .into_iter()
                        .map(|(id, pair)| (id, pair))
                        .collect();
                id_pair_list.sort_by(|(id_1, _), (id_2, _)| id_1.cmp(id_2));
                (chrom, id_pair_list)
            })
            .collect();

        let ordered_chroms: Vec<String> = {
            let mut ordered_chroms: Vec<String> = chrom_to_sorted_id_pairs
                .keys()
                .map(|k| k.to_string())
                .collect();
            ordered_chroms.sort();
            ordered_chroms
        };

        let mut writer = BedWriter::new(out_path)?;
        for chrom in ordered_chroms.into_iter() {
            let sorted_id_pairs = chrom_to_sorted_id_pairs.get(&chrom).unwrap();
            for (id, pair) in sorted_id_pairs.into_iter() {
                let midpoint = match pair.mid_point() {
                    None => continue,
                    Some(mid) => mid,
                };
                let start = midpoint - left_offset;
                let end = midpoint + right_offset;
                let data_line = BedDataLine {
                    chrom: chrom.to_string(),
                    start,
                    end,
                    name: Some(id.to_string()),
                    score: None::<f64>,
                    strand: None,
                };
                writer.write_bed_line(&data_line)?;
            }
        }
        Ok(PairedEndCollatorResult {
            num_pcr_duplicates,
        })
    }

    fn get_pair_id(name: &str) -> Option<(String, String)> {
        let toks: Vec<String> =
            name.split('/').map(|tok| tok.to_string()).collect();

        match toks.split_last() {
            None => None,
            Some((last, first)) => Some((first.join("/"), last.to_string())),
        }
    }
}

/// Each pair of coordinates represent the start, and end field read from
/// the bed file, where the start is incluisve, and the end is exclusive.
#[derive(Copy, Clone, PartialOrd, PartialEq, Hash, Debug)]
pub struct PairedEndReadPair {
    pub postive_strand: Option<ContiguousIntegerSet<Coordinate>>,
    pub negative_strand: Option<ContiguousIntegerSet<Coordinate>>,
}

impl PairedEndReadPair {
    pub fn new() -> PairedEndReadPair {
        PairedEndReadPair {
            postive_strand: None,
            negative_strand: None,
        }
    }

    fn is_strand_valid(
        strand: &Option<ContiguousIntegerSet<Coordinate>>,
    ) -> bool {
        match strand {
            None => false,
            Some(strand) => !strand.is_empty(),
        }
    }

    pub fn has_valid_positive_strand(&self) -> bool {
        PairedEndReadPair::is_strand_valid(&self.postive_strand)
    }

    pub fn has_valid_negative_strand(&self) -> bool {
        PairedEndReadPair::is_strand_valid(&self.negative_strand)
    }

    pub fn both_strands_valid(&self) -> bool {
        self.has_valid_positive_strand() && self.has_valid_negative_strand()
    }

    pub fn min_coordiante(&self) -> Option<Coordinate> {
        match self.postive_strand {
            None => self.negative_strand.map(|interval| interval.get_start()),
            Some(positive_strand) => match self.negative_strand {
                None => Some(positive_strand.get_start()),
                Some(negative_strand) => Some(std::cmp::min(
                    positive_strand.get_start(),
                    negative_strand.get_start(),
                )),
            },
        }
    }

    pub fn max_coordinate(&self) -> Option<Coordinate> {
        match self.postive_strand {
            None => self.negative_strand.map(|interval| interval.get_end()),
            Some(positive_strand) => match self.negative_strand {
                None => Some(positive_strand.get_end()),
                Some(negative_strand) => Some(std::cmp::max(
                    positive_strand.get_end(),
                    negative_strand.get_end(),
                )),
            },
        }
    }

    pub fn mid_point(&self) -> Option<Coordinate> {
        let min = match self.min_coordiante() {
            None => return None,
            Some(min) => min,
        };
        let max = match self.max_coordinate() {
            None => return None,
            Some(max) => max,
        };
        Some(min + (max - min) / 2)
    }
}

#[derive(Copy, Clone, Debug)]
pub struct PairedEndCollatorResult {
    pub num_pcr_duplicates: i64,
}

#[cfg(test)]
mod tests {
    use crate::bed::{
        paired_end_collator::PairedEndCollator, Bed, BedDataLine,
        BedDataLineIter,
    };
    use math::traits::ToIterator;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_collator() {
        let mut file = NamedTempFile::new().unwrap();
        write!(
            &mut file,
            "chr1 100 158 id_1/1 4 +\n\
            chr1 200 250 id_2/1 5 +\n\
            chr1 200 250 id_2/1 5 +\n\
            chr1 350 400 id_1/2 5 -\n\
            chr1 400 450 id_2/2 5 -\n\
            chr1 1000 1040 id_4/1 9 +\n\
            chr2 200 250 id_3/1 2 +\n\
            chr1 1200 1243 id_4/2 8 -\n\
            chr2 400 450 id_3/2 2 -\n\
            chr1 10000 10050 id_3/2 2 -\n\
            chr3 203 230 id_5/2 4 +\n"
        )
        .unwrap();
        let file_path = file.into_temp_path();
        let path_str = file_path.to_str().unwrap().to_string();
        let collator = PairedEndCollator::new(&path_str);

        let out_path = NamedTempFile::new().unwrap().into_temp_path();
        let out_path_str = out_path.to_str().unwrap().to_string();
        collator.write(&out_path_str, 75, 75).unwrap();

        let expected_data_lines = vec![
            BedDataLine {
                chrom: "chr1".into(),
                start: 175,
                end: 325,
                name: Some("id_1".into()),
                score: Some(0.),
                strand: None,
            },
            BedDataLine {
                chrom: "chr1".into(),
                start: 250,
                end: 400,
                name: Some("id_2".into()),
                score: Some(0.),
                strand: None,
            },
            BedDataLine {
                chrom: "chr1".into(),
                start: 9950,
                end: 10100,
                name: Some("id_3".into()),
                score: Some(0.),
                strand: None,
            },
            BedDataLine {
                chrom: "chr1".into(),
                start: 1046,
                end: 1196,
                name: Some("id_4".into()),
                score: Some(0.),
                strand: None,
            },
            BedDataLine {
                chrom: "chr2".into(),
                start: 250,
                end: 400,
                name: Some("id_3".into()),
                score: Some(0.),
                strand: None,
            },
            BedDataLine {
                chrom: "chr3".into(),
                start: 141,
                end: 291,
                name: Some("id_5".into()),
                score: Some(0.),
                strand: None,
            },
        ];

        let bed = Bed::new(&out_path_str);
        for (line, expected) in (bed.to_iter(): BedDataLineIter<f64>)
            .zip(expected_data_lines.into_iter())
        {
            assert_eq!(line, expected);
        }
    }
}
