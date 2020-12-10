use crate::{
    bed::{
        bed_writer::BedWriter, Bed, BedDataLine, BedDataLineIter, Coordinate,
    },
    error::Error,
    util::Strand,
};
use math::{
    histogram::Histogram,
    interval::traits::Interval,
    set::{
        contiguous_integer_set::ContiguousIntegerSet,
        traits::{Intersect, Set},
    },
    traits::{Collecting, ToIterator},
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
        min_distance: Option<Coordinate>,
        max_distance: Option<Coordinate>,
        debug: bool,
    ) -> Result<PairedEndCollatorResult, Error>
where {
        let mut visited = HashSet::new();
        let mut chrom_to_id_to_pair = HashMap::new();
        let mut num_conflicting_reads = 0i64;
        let mut num_pcr_duplicates = 0i64;
        let mut num_overlapping_pairs = 0i64;

        let bed = Bed::new(&self.path, true);
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
                if debug {
                    eprintln!(
                        "PCR duplciate (chrom, start, end, strand): \
                        ({}, {}, {}, {:?})",
                        chrom, start, end, strand
                    )
                }
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
                            .or_insert_with(HashMap::new);

                        let pair = id_to_pair
                            .entry(pair_id.clone())
                            .or_insert_with(PairedEndReadPair::new);

                        match strand {
                            Strand::Positive => match pair.positive_strand {
                                None => {
                                    pair.positive_strand =
                                        Some(ContiguousIntegerSet::new(
                                            start,
                                            end - 1,
                                        ));
                                }
                                Some(old_strand) => {
                                    num_conflicting_reads += 1;
                                    eprintln!(
                                        "{}",
                                        format!(
                                            "conflicting read with id {} on \
                                            the positive strand with \
                                            coordinates ({}, {}) and ({}, {})",
                                            pair_id,
                                            old_strand.get_start(),
                                            old_strand.get_end(),
                                            start,
                                            end
                                        )
                                    );
                                }
                            },
                            Strand::Negative => match pair.negative_strand {
                                None => {
                                    pair.negative_strand =
                                        Some(ContiguousIntegerSet::new(
                                            start,
                                            end - 1,
                                        ));
                                }
                                Some(old_strand) => {
                                    num_conflicting_reads += 1;
                                    eprintln!(
                                        "{}",
                                        format!(
                                            "conflicting read with id {} on \
                                            the negative strand with \
                                            coordinates ({}, {}) and ({}, {})",
                                            pair_id,
                                            old_strand.get_start(),
                                            old_strand.get_end(),
                                            start,
                                            end
                                        )
                                    );
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

        let mut distance_histogram = Histogram::new(None, 20, 0, 1000)?;
        let mut writer = BedWriter::new(out_path)?;
        for chrom in ordered_chroms.into_iter() {
            let sorted_id_pairs = chrom_to_sorted_id_pairs.get(&chrom).unwrap();
            for (id, pair) in sorted_id_pairs.iter() {
                // the pair must cover a valid distance
                if let Some(distance) = pair.distance_covered() {
                    distance_histogram.collect(distance);

                    if let Some(min_distance) = min_distance {
                        if distance < min_distance {
                            continue;
                        }
                    }

                    if let Some(max_distance) = max_distance {
                        if distance > max_distance {
                            continue;
                        }
                    }

                    let midpoint = match pair.paired_midpoint() {
                        None => continue,
                        Some(mid) => mid,
                    };

                    if let Some(is_overlapping) = pair.is_overlapping() {
                        num_overlapping_pairs += is_overlapping as i64;
                    }

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
        }
        Ok(PairedEndCollatorResult {
            num_conflicting_reads,
            num_overlapping_pairs,
            num_pcr_duplicates,
            distance_histogram,
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
    // note that ContiguousIntegerSet is right inclusive
    pub positive_strand: Option<ContiguousIntegerSet<Coordinate>>,
    pub negative_strand: Option<ContiguousIntegerSet<Coordinate>>,
}

impl PairedEndReadPair {
    pub fn new() -> PairedEndReadPair {
        PairedEndReadPair {
            positive_strand: None,
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
        PairedEndReadPair::is_strand_valid(&self.positive_strand)
    }

    pub fn has_valid_negative_strand(&self) -> bool {
        PairedEndReadPair::is_strand_valid(&self.negative_strand)
    }

    pub fn both_strands_valid(&self) -> bool {
        self.has_valid_positive_strand() && self.has_valid_negative_strand()
    }

    pub fn get_valid_positive_strand(
        &self,
    ) -> Option<ContiguousIntegerSet<Coordinate>> {
        match self.positive_strand {
            None => None,
            Some(strand) => {
                if strand.is_empty() {
                    None
                } else {
                    Some(strand)
                }
            }
        }
    }

    pub fn get_valid_negative_strand(
        &self,
    ) -> Option<ContiguousIntegerSet<Coordinate>> {
        match self.negative_strand {
            None => None,
            Some(strand) => {
                if strand.is_empty() {
                    None
                } else {
                    Some(strand)
                }
            }
        }
    }

    /// Returns the minimum coordinate (inclusive) across the positive and
    /// the negative strands, i.e., if the positive strand is [a, b) and the
    /// negative strand is [c, d), the result is min(a, c).
    pub fn min_coordinate(&self) -> Option<Coordinate> {
        match self.positive_strand {
            None => match self.negative_strand {
                None => None,
                Some(negative_strand) => {
                    negative_strand.get_start_if_nonempty()
                }
            },
            Some(positive_strand) => match self.negative_strand {
                None => positive_strand.get_start_if_nonempty(),
                Some(negative_strand) => {
                    let positive_start =
                        positive_strand.get_start_if_nonempty();
                    let negative_start =
                        negative_strand.get_start_if_nonempty();

                    match positive_start {
                        None => negative_start,
                        Some(p) => match negative_start {
                            None => Some(p),
                            Some(n) => Some(std::cmp::min(p, n)),
                        },
                    }
                }
            },
        }
    }

    /// Returns the maximum coordinate (inclusive) across the positive and
    /// the negative strands, i.e., if the positive strand is [a, b) and the
    /// negative strand is [c, d), the result is max(b, d) - 1.
    pub fn max_coordinate(&self) -> Option<Coordinate> {
        match self.positive_strand {
            None => match self.negative_strand {
                None => None,
                Some(negative_strand) => negative_strand.get_end_if_nonempty(),
            },
            Some(positive_strand) => match self.negative_strand {
                None => positive_strand.get_end_if_nonempty(),
                Some(negative_strand) => {
                    let positive_end = positive_strand.get_end_if_nonempty();
                    let negative_end = negative_strand.get_end_if_nonempty();

                    match positive_end {
                        None => negative_end,
                        Some(p) => match negative_end {
                            None => Some(p),
                            Some(n) => Some(std::cmp::max(p, n)),
                        },
                    }
                }
            },
        }
    }

    /// If the positive strand is [a, b) and the negative strand is [c, d), the
    /// mid-point will be (max(a, c) + max(b, d) - 1) / 2
    pub fn paired_midpoint(&self) -> Option<Coordinate> {
        if self.both_strands_valid() {
            let min = match self.min_coordinate() {
                None => return None,
                Some(min) => min,
            };
            let max = match self.max_coordinate() {
                None => return None,
                Some(max) => max,
            };
            Some(min + (max - min) / 2)
        } else {
            None
        }
    }

    /// Let the positive strand be [a, b), let the negative strand be [c, d).
    /// Returns `Some(d - a)`, if both the positive and negative strands
    /// are present, have nonempty interval ranges, and a <= c. Otherwise
    /// returns `None`.
    pub fn distance_covered(&self) -> Option<Coordinate> {
        match self.positive_strand {
            None => None,
            Some(positive) => match self.negative_strand {
                None => None,
                Some(negative) => match positive.get_start_if_nonempty() {
                    None => None,
                    Some(p) => match negative.get_end_if_nonempty() {
                        None => None,
                        Some(n) => {
                            if positive.get_start() <= negative.get_start()
                                && positive.get_end() <= negative.get_end()
                            {
                                Some(n - p + 1)
                            } else {
                                None
                            }
                        }
                    },
                },
            },
        }
    }

    pub fn is_overlapping(&self) -> Option<bool> {
        if let Some(positive) = self.get_valid_positive_strand() {
            if let Some(negative) = self.get_valid_negative_strand() {
                return Some(
                    positive.has_non_empty_intersection_with(&negative),
                );
            }
        }
        None
    }
}

#[derive(Clone, Debug)]
pub struct PairedEndCollatorResult {
    // Two reads are considered conflicting if they have the same name field
    // but different (start, end, strand) tuples.
    pub num_conflicting_reads: i64,

    // The number of paired-end reads where the pair overlap each other.
    pub num_overlapping_pairs: i64,

    // Two reads are considered as PCR duplicates if they have the same
    // (start, end, strand) tuple but a different name field.
    pub num_pcr_duplicates: i64,

    // A histogram for the distance between the maximum and the minimum
    // basepair coordinates covered by each pair of paired reads.
    pub distance_histogram: Histogram<Coordinate>,
}

#[cfg(test)]
mod tests {
    use crate::bed::{
        paired_end_collator::{PairedEndCollator, PairedEndReadPair},
        Bed, BedDataLine, BedDataLineIter,
    };
    use math::{
        set::contiguous_integer_set::ContiguousIntegerSet, traits::ToIterator,
    };
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_valid_strands() {
        macro_rules! expected_bools {
            ($pair:expr, $valid_pos:expr, $valid_neg:expr) => {
                assert_eq!($pair.has_valid_positive_strand(), $valid_pos);
                assert_eq!($pair.has_valid_negative_strand(), $valid_neg);
                assert_eq!(
                    $pair.both_strands_valid(),
                    $valid_pos && $valid_neg
                );
            };
        }

        let mut pair = PairedEndReadPair::new();
        expected_bools!(pair, false, false);

        pair.positive_strand = Some(ContiguousIntegerSet::new(100, 109));
        expected_bools!(pair, true, false);

        pair.negative_strand = Some(ContiguousIntegerSet::new(200, 249));
        expected_bools!(pair, true, true);

        pair.positive_strand = None;
        expected_bools!(pair, false, true);

        pair.positive_strand = Some(ContiguousIntegerSet::new(200, 109));
        expected_bools!(pair, false, true);

        pair.negative_strand = Some(ContiguousIntegerSet::new(250, 249));
        expected_bools!(pair, false, false);
    }

    #[test]
    fn test_get_valid_strands() {
        macro_rules! expect_strands {
            ($pair:expr, $pos_strand:expr, $neg_strand:expr) => {
                assert_eq!($pair.get_valid_positive_strand(), $pos_strand);
                assert_eq!($pair.get_valid_negative_strand(), $neg_strand);
            };
        }
        let mut pair = PairedEndReadPair::new();
        expect_strands!(pair, None, None);

        let invalid_1 = Some(ContiguousIntegerSet::new(200, 109));
        let invalid_2 = Some(ContiguousIntegerSet::new(400, 249));
        pair.positive_strand = invalid_1;
        pair.negative_strand = invalid_2;
        expect_strands!(pair, None, None);

        let s1 = Some(ContiguousIntegerSet::new(100, 109));
        pair.positive_strand = s1;
        expect_strands!(pair, s1, None);

        let s2 = Some(ContiguousIntegerSet::new(200, 249));
        pair.negative_strand = s2;
        expect_strands!(pair, s1, s2);

        pair.positive_strand = invalid_1;
        expect_strands!(pair, None, s2);
    }

    #[test]
    fn test_min_max_coordinates() {
        macro_rules! expect_min_max {
            ($pair:expr, $min:expr, $max:expr) => {
                assert_eq!($pair.min_coordinate(), $min);
                assert_eq!($pair.max_coordinate(), $max);
            };
        }
        let mut pair = PairedEndReadPair::new();
        expect_min_max!(pair, None, None);

        // invalid, None
        pair.positive_strand = Some(ContiguousIntegerSet::new(200, 109));
        expect_min_max!(pair, None, None);

        // invalid, invalid
        pair.negative_strand = Some(ContiguousIntegerSet::new(100, 29));
        expect_min_max!(pair, None, None);

        // None, invalid
        pair.positive_strand = None;
        expect_min_max!(pair, None, None);

        // valid, invalid
        pair.positive_strand = Some(ContiguousIntegerSet::new(100, 109));
        expect_min_max!(pair, Some(100), Some(109));

        // valid, valid
        pair.negative_strand = Some(ContiguousIntegerSet::new(300, 350));
        expect_min_max!(pair, Some(100), Some(350));

        // valid, valid
        pair.negative_strand = Some(ContiguousIntegerSet::new(10, 20));
        expect_min_max!(pair, Some(10), Some(109));

        // valid, valid
        pair.negative_strand = Some(ContiguousIntegerSet::new(105, 105));
        expect_min_max!(pair, Some(100), Some(109));

        // valid, valid
        pair.negative_strand = Some(ContiguousIntegerSet::new(10, 400));
        expect_min_max!(pair, Some(10), Some(400));

        // valid, valid
        pair.negative_strand = Some(ContiguousIntegerSet::new(0, 105));
        expect_min_max!(pair, Some(0), Some(109));

        // valid, valid
        pair.negative_strand = Some(ContiguousIntegerSet::new(105, 1000));
        expect_min_max!(pair, Some(100), Some(1000));

        // None, valid
        pair.positive_strand = None;
        expect_min_max!(pair, Some(105), Some(1000));

        // valid, None
        pair.positive_strand = Some(ContiguousIntegerSet::new(200, 240));
        pair.negative_strand = None;
        expect_min_max!(pair, Some(200), Some(240));
    }

    #[test]
    fn test_paired_midpoint() {
        let mut pair = PairedEndReadPair::new();
        assert_eq!(pair.paired_midpoint(), None);

        pair.positive_strand = Some(ContiguousIntegerSet::new(100, 109));
        assert_eq!(pair.paired_midpoint(), None);

        pair.negative_strand = Some(ContiguousIntegerSet::new(100, 109));
        assert_eq!(pair.paired_midpoint(), Some(104));

        pair.negative_strand = Some(ContiguousIntegerSet::new(200, 249));
        assert_eq!(pair.paired_midpoint(), Some(174));

        pair.negative_strand = Some(ContiguousIntegerSet::new(200, 250));
        assert_eq!(pair.paired_midpoint(), Some(175));

        pair.positive_strand = Some(ContiguousIntegerSet::new(300, 300));
        assert_eq!(pair.paired_midpoint(), Some(250));

        pair.positive_strand = Some(ContiguousIntegerSet::new(100, 300));
        assert_eq!(pair.paired_midpoint(), Some(200));

        pair.positive_strand = Some(ContiguousIntegerSet::new(100, 301));
        assert_eq!(pair.paired_midpoint(), Some(200));

        pair.positive_strand = Some(ContiguousIntegerSet::new(10, 10));
        pair.negative_strand = Some(ContiguousIntegerSet::new(20, 20));
        assert_eq!(pair.paired_midpoint(), Some(15));

        pair.negative_strand = Some(ContiguousIntegerSet::new(20, 21));
        assert_eq!(pair.paired_midpoint(), Some(15));

        pair.positive_strand = None;
        assert_eq!(pair.paired_midpoint(), None);
    }

    #[test]
    fn test_pair_distance_covered() {
        let mut pair = PairedEndReadPair::new();
        assert_eq!(pair.distance_covered(), None);

        // only positive strand present
        pair.positive_strand = Some(ContiguousIntegerSet::new(0, 49));
        assert_eq!(pair.distance_covered(), None);

        // empty negative strand ragne
        pair.negative_strand = Some(ContiguousIntegerSet::new(100, 90));
        assert_eq!(pair.distance_covered(), None);

        pair.negative_strand = Some(ContiguousIntegerSet::new(100, 199));
        assert_eq!(pair.distance_covered(), Some(200));

        // same start coordinates for the positive and the negative strands
        pair.positive_strand = Some(ContiguousIntegerSet::new(100, 149));
        assert_eq!(pair.distance_covered(), Some(100));

        // positive strand start is after the negative starnd start
        pair.positive_strand = Some(ContiguousIntegerSet::new(101, 149));
        assert_eq!(pair.distance_covered(), None);

        // empty positive strand ragne
        pair.positive_strand = Some(ContiguousIntegerSet::new(10, 5));
        assert_eq!(pair.distance_covered(), None);

        // positive strand contains the negative strand
        pair.positive_strand = Some(ContiguousIntegerSet::new(0, 1000));
        assert_eq!(pair.distance_covered(), None);
    }

    #[test]
    fn test_is_overlapping() {
        let mut pair = PairedEndReadPair::new();
        assert_eq!(pair.is_overlapping(), None);

        pair.positive_strand = Some(ContiguousIntegerSet::new(0, 49));
        assert_eq!(pair.is_overlapping(), None);

        // empty negative strand
        pair.negative_strand = Some(ContiguousIntegerSet::new(40, 30));
        assert_eq!(pair.is_overlapping(), None);

        // valid overlapping negative strand
        pair.negative_strand = Some(ContiguousIntegerSet::new(49, 50));
        assert_eq!(pair.is_overlapping(), Some(true));

        pair.negative_strand = Some(ContiguousIntegerSet::new(40, 100));
        assert_eq!(pair.is_overlapping(), Some(true));

        pair.positive_strand = None;
        assert_eq!(pair.is_overlapping(), None);

        // is_overlapping will return true is one strand's range is a subset of
        // the other's
        pair.positive_strand = Some(ContiguousIntegerSet::new(0, 1000));
        assert_eq!(pair.is_overlapping(), Some(true));

        pair.positive_strand = Some(ContiguousIntegerSet::new(0, 29));
        assert_eq!(pair.is_overlapping(), Some(false));

        pair.positive_strand = Some(ContiguousIntegerSet::new(0, 39));
        assert_eq!(pair.is_overlapping(), Some(false));

        pair.positive_strand = Some(ContiguousIntegerSet::new(101, 109));
        assert_eq!(pair.is_overlapping(), Some(false));

        // is_overlapping will return true even if the positive strand's
        // start is larger than the negative strand's start
        pair.positive_strand = Some(ContiguousIntegerSet::new(100, 109));
        assert_eq!(pair.is_overlapping(), Some(true));

        pair.negative_strand = Some(ContiguousIntegerSet::new(0, 1000));
        assert_eq!(pair.is_overlapping(), Some(true));
    }

    #[test]
    fn test_collator() {
        macro_rules! data_line {
            ($chrom:expr, $start:expr, $end:expr, $name:expr) => {
                BedDataLine {
                    chrom: $chrom.into(),
                    start: $start,
                    end: $end,
                    name: Some($name.into()),
                    score: Some(1.),
                    strand: None,
                }
            };
        }

        let mut file = NamedTempFile::new().unwrap();
        write!(
            &mut file,
            "chr1 100 158 id_1/1 4 +\n\
            chr1 200 250 id_2/1 5 +\n\
            chr1 200 250 id_duplicate_1/1 5 +\n\
            chr1 200 250 id_duplicate_2/1 3 +\n\
            chr1 350 400 id_1/2 5 -\n\
            chr1 400 450 id_2/2 5 -\n\
            chr1 1000 1040 id_4/1 9 +\n\
            chr2 200 250 id_3/1 2 +\n\
            chr1 1200 1243 id_4/2 8 -\n\
            chr2 400 450 id_3/2 2 -\n\
            chr1 10000 10050 id_6/2 2 -\n\
            chr2 10000 10010 id_conflict_1/1 2 +\n\
            chr2 10000 10012 id_conflict_1/1 2 +\n\
            chr2 900 1000 id_conflict_2/1 2 -\n\
            chr2 920 1000 id_conflict_2/1 2 -\n\
            chr2 1000 10012 id_conflict_2/1 2 -\n\
            chr3 203 230 id_5/2 4 +\n"
        )
        .unwrap();
        let file_path = file.into_temp_path();
        let path_str = file_path.to_str().unwrap().to_string();
        let collator = PairedEndCollator::new(&path_str);

        let out_path = NamedTempFile::new().unwrap().into_temp_path();
        let out_path_str = out_path.to_str().unwrap().to_string();
        let result = collator
            .write(&out_path_str, 75, 75, None, None, false)
            .unwrap();

        let expected_data_lines = vec![
            data_line!("chr1", 174, 324, "id_1"),
            data_line!("chr1", 249, 399, "id_2"),
            data_line!("chr1", 1046, 1196, "id_4"),
            data_line!("chr2", 249, 399, "id_3"),
        ];

        let bed = Bed::new(&out_path_str, true);
        let actual_data_lines: Vec<BedDataLine<f64>> =
            (bed.to_iter(): BedDataLineIter<f64>).collect();

        for (line, expected) in
            actual_data_lines.iter().zip(expected_data_lines.iter())
        {
            assert_eq!(line, expected);
        }
        assert_eq!(actual_data_lines.len(), expected_data_lines.len());

        assert_eq!(result.num_pcr_duplicates, 2);
        assert_eq!(result.num_conflicting_reads, 3);
    }
}
