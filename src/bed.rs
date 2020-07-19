//! An interface to the BED track format file as specified in
//! https://genome.ucsc.edu/FAQ/FAQformat.html#format1

use std::{
    collections::HashMap,
    fmt::Debug,
    fs::File,
    io::{BufRead, BufReader},
    marker::PhantomData,
    str::FromStr,
};

use analytic::{set::ordered_integer_set::ContiguousIntegerSet, traits::ToIterator};
use num::Float;

use crate::{
    error::Error,
    util::{get_buf, Strand},
};

pub struct Bed {
    filepath: String,
}

impl Bed {
    pub fn new(filepath: &str) -> Result<Bed, Error> {
        Ok(Bed {
            filepath: filepath.to_string(),
        })
    }

    #[inline]
    pub fn get_filepath(&self) -> &str {
        &self.filepath
    }

    pub fn get_chrom_to_interval_to_val<D, E>(
        &self,
    ) -> Result<HashMap<String, HashMap<ContiguousIntegerSet<Coordinate>, D>>, Error>
    where
        D: Float + FromStr<Err = E>,
        E: Debug, {
        let mut chrom_to_interval_to_val = HashMap::new();
        for BedDataLine {
            chrom,
            start,
            end,
            name: _,
            score,
            strand: _,
        } in self.to_iter(): BedDataLineIter<D>
        {
            let score = if let Some(score) = score {
                score
            } else {
                return Err(Error::Generic(
                    "the BED file does not have a score field".into(),
                ));
            };
            let interval_to_val = chrom_to_interval_to_val
                .entry(chrom)
                .or_insert_with(HashMap::new);
            interval_to_val.insert(ContiguousIntegerSet::new(start, end - 1), score);
        }
        Ok(chrom_to_interval_to_val)
    }
}

impl<D, E> ToIterator<'_, BedDataLineIter<D>, <BedDataLineIter<D> as Iterator>::Item> for Bed
where
    D: Float + FromStr<Err = E>,
    E: Debug,
{
    fn to_iter(&self) -> BedDataLineIter<D> {
        let buf = get_buf(&self.filepath).unwrap();
        BedDataLineIter {
            buf,
            filename: self.filepath.clone(),
            phantom: PhantomData,
        }
    }
}

/// Data type of the Bedgraph coordinates
pub type Coordinate = i64;

/// `BedDataLine` corresponds to a line of data in the Bed
/// file, where each line is of the form
/// `chrom start end name score strand ...`,
/// where the first three fields are required, and the remaining 9 fields are optional.
///
/// The [start, end) is a zero-based left-closed right-open coordinate range.
pub struct BedDataLine<D> {
    pub chrom: String,
    pub start: Coordinate,
    pub end: Coordinate,
    pub name: Option<String>,
    pub score: Option<D>,
    pub strand: Option<Strand>,
}

pub struct BedDataLineIter<D> {
    buf: BufReader<File>,
    filename: String,
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
                let score = toks.next().map(|score| score.parse::<D>().unwrap());
                let strand = toks.next().and_then(|strand| {
                    Strand::new(strand).expect("failed to parse the strand symbol")
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

#[cfg(test)]
mod tests {
    use std::{
        collections::HashMap,
        io::{BufWriter, Write},
    };

    use analytic::set::ordered_integer_set::ContiguousIntegerSet;
    use tempfile::NamedTempFile;

    use crate::bed::{Bed, Coordinate};

    #[test]
    fn test_get_chrom_to_interval_to_val() {
        let file = NamedTempFile::new().unwrap();
        {
            let mut writer = BufWriter::new(&file);
            writer
                .write_fmt(format_args!(
                    "chr1 100 200 name_1 3.5\n\
                    chr1 200 350 name_2 4.0\n\
                    chr3 1000 3000 name_3 -0.3\n\
                    chr1 400 450 name_4 -0.9\n"
                ))
                .unwrap();
        }
        let bed = Bed::new(file.path().to_str().unwrap()).unwrap();
        let chrom_to_interval_to_val = bed.get_chrom_to_interval_to_val().unwrap();
        let expected_chr1: HashMap<ContiguousIntegerSet<Coordinate>, f64> = [
            (ContiguousIntegerSet::new(100, 199), 3.5),
            (ContiguousIntegerSet::new(200, 349), 4.),
            (ContiguousIntegerSet::new(400, 449), -0.9),
        ]
        .iter()
        .map(|(a, b)| (*a, *b))
        .collect();
        let expected_chr3: HashMap<ContiguousIntegerSet<Coordinate>, f64> =
            [(ContiguousIntegerSet::new(1000, 2999), -0.3)]
                .iter()
                .map(|(a, b)| (*a, *b))
                .collect();
        assert_eq!(chrom_to_interval_to_val["chr1"], expected_chr1);
        assert_eq!(chrom_to_interval_to_val["chr3"], expected_chr3);
    }
}
