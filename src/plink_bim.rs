use std::collections::{HashMap, HashSet};
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, Lines, Read, Seek, SeekFrom};

use analytic::set::ordered_integer_set::OrderedIntegerSet;
use analytic::traits::Collecting;

use crate::error::Error;

pub const CHROM_FIELD_INDEX: usize = 0;
pub const VARIANT_ID_FIELD_INDEX: usize = 1;
pub const COORDINATE_FIELD_INDEX: usize = 3;
pub const FIRST_ALLELE_FIELD_INDEX: usize = 4;
pub const SECOND_ALLELE_FIELD_INDEX: usize = 5;

pub type PartitionKeyType = String;

pub struct PlinkBim {
    filepath: String,
    buf: BufReader<File>,
    // maps partition_id to the file line indices
    fileline_partitions: Option<HashMap<PartitionKeyType, OrderedIntegerSet<usize>>>,
}

impl PlinkBim {
    pub fn new(filepath: &str) -> Result<PlinkBim, Error> {
        let buf = BufReader::new(OpenOptions::new().read(true).open(filepath)?);
        Ok(PlinkBim {
            filepath: filepath.to_string(),
            buf,
            fileline_partitions: None,
        })
    }

    pub fn new_with_partitions(filepath: &str, partitions: HashMap<PartitionKeyType, OrderedIntegerSet<usize>>) -> Result<PlinkBim, Error> {
        let mut bim = PlinkBim::new(filepath)?;
        bim.set_fileline_partitions(Some(partitions));
        Ok(bim)
    }

    pub fn new_with_partition_file(bim_filepath: &str, partition_filepath: &str) -> Result<PlinkBim, Error> {
        let bim = PlinkBim::new(bim_filepath)?;
        bim.into_partitioned_by_file(partition_filepath)
    }

    fn get_id_and_partition_from_partition_fileline_iter<'a, T: Iterator<Item=&'a str>>(iter: &'a mut T) -> Option<(String, String)> {
        let id = match iter.next() {
            None => return None,
            Some(id) => id.to_string(),
        };
        match iter.next() {
            None => None,
            Some(partition) => Some((id, partition.to_string()))
        }
    }

    pub fn get_fileline_partitions(&self) -> Option<HashMap<String, OrderedIntegerSet<usize>>> {
        match &self.fileline_partitions {
            None => None,
            Some(p) => Some(p.clone())
        }
    }

    pub fn get_fileline_partitions_by_ref(&self) -> Option<&HashMap<String, OrderedIntegerSet<usize>>> {
        match &self.fileline_partitions {
            None => None,
            Some(p) => Some(p)
        }
    }

    /// each line in the partition file has two fields separated by space:
    /// variant_id assigned_partition
    pub fn get_fileline_partitions_from_partition_file(&mut self, partition_file_path: &str) -> Result<HashMap<PartitionKeyType, OrderedIntegerSet<usize>>, Error> {
        let mut id_to_partition: HashMap<String, PartitionKeyType> = HashMap::new();
        for line in BufReader::new(OpenOptions::new().read(true).open(partition_file_path)?).lines() {
            match PlinkBim::get_id_and_partition_from_partition_fileline_iter(&mut line?.split_whitespace()) {
                None => return Err(Error::BadFormat(String::from(
                    "each line in the partition file should have two fields: id partition_assignment"))),
                Some((id, partition)) => id_to_partition.insert(id.to_string(), partition.to_string())
            };
        }
        let mut visited_ids = HashSet::new();
        self.reset_buf()?;
        let mut partitions: HashMap<PartitionKeyType, OrderedIntegerSet<usize>> = HashMap::new();
        for (i, line) in self.buf.by_ref().lines().enumerate() {
            let id = match line?.split_whitespace().nth(VARIANT_ID_FIELD_INDEX) {
                None => return Err(Error::BadFormat(format!(
                    "failed to read the variant id in line {} in bim file: {}", i + 1, self.filepath))),
                Some(id) => id.to_string()
            };
            visited_ids.insert(id.to_string());
            match id_to_partition.get(&id) {
                None => {}
                Some(p) => partitions.entry(p.to_owned()).or_insert(OrderedIntegerSet::new()).collect(i)
            }
        }
        let num_ids_not_in_bim = id_to_partition.keys().filter(|&k| !visited_ids.contains(k)).count();
        if num_ids_not_in_bim > 0 {
            Err(Error::Generic(format!("{} ID(s) from the partition file {} are not in the bim file {}",
                                       num_ids_not_in_bim, partition_file_path, self.filepath)))
        } else {
            Ok(partitions)
        }
    }

    #[inline]
    pub fn get_filepath(&self) -> &str {
        &self.filepath
    }

    #[inline]
    pub fn reset_buf(&mut self) -> Result<(), Error> {
        self.buf.seek(SeekFrom::Start(0))?;
        Ok(())
    }

    pub fn lines(&mut self) -> Lines<&mut BufReader<File>> {
        self.buf.by_ref().lines()
    }

    pub fn get_all_chroms(&mut self) -> Result<HashSet<String>, Error> {
        self.reset_buf()?;
        return Ok(self.buf.by_ref().lines().filter_map(|l|
            l.unwrap().split_whitespace().nth(CHROM_FIELD_INDEX).and_then(|s| Some(s.to_string()))
        ).collect());
    }

    pub fn get_chrom_fileline_positions(&mut self, chrom: &str) -> Result<OrderedIntegerSet<usize>, Error> {
        let mut set = OrderedIntegerSet::new();
        self.reset_buf()?;
        for (i, l) in self.buf.by_ref().lines().enumerate() {
            if l.unwrap()
                .split_whitespace()
                .nth(CHROM_FIELD_INDEX).unwrap() == chrom {
                set.collect(i);
            }
        }
        Ok(set)
    }

    pub fn get_chrom_to_fileline_positions(&mut self) -> Result<HashMap<String, OrderedIntegerSet<usize>>, Error> {
        let mut chrom_to_positions = HashMap::new();
        for chrom in self.get_all_chroms()? {
            let positions = self.get_chrom_fileline_positions(&chrom)?;
            chrom_to_positions.insert(chrom, positions);
        }
        Ok(chrom_to_positions)
    }

    #[inline]
    pub fn set_fileline_partitions(&mut self, partitions: Option<HashMap<PartitionKeyType, OrderedIntegerSet<usize>>>) {
        self.fileline_partitions = partitions;
    }

    pub fn into_partitioned_by_chrom(mut self) -> Result<PlinkBim, Error> {
        let partitions = self.get_chrom_to_fileline_positions()?;
        self.set_fileline_partitions(Some(partitions));
        Ok(self)
    }

    pub fn into_partitioned_by_file(mut self, partition_file: &str) -> Result<PlinkBim, Error> {
        let partitions = self.get_fileline_partitions_from_partition_file(partition_file)?;
        self.set_fileline_partitions(Some(partitions));
        Ok(self)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::{HashMap, HashSet};
    use std::fs::OpenOptions;
    use std::io::{BufWriter, Write};

    use tempfile::NamedTempFile;

    use analytic::set::ordered_integer_set::{ContiguousIntegerSet, OrderedIntegerSet};

    use super::PlinkBim;

    #[inline]
    fn write_bim_line<W: Write>(buf_writer: &mut BufWriter<W>, chrom: &str, id: &str, coordinate: u64,
                                first_allele: char, second_allele: char) {
        buf_writer.write_fmt(format_args!("{} {} {} {} {} {}\n", chrom, id, 0, coordinate, first_allele, second_allele)).unwrap();
    }

    #[test]
    fn test_get_chrom_to_positions() {
        let file = NamedTempFile::new().unwrap();
        {
            let mut writer = BufWriter::new(&file);
            for (chrom, coord, id) in &[
                (1, 2, "chr1:2"), (1, 3, "chr1:3"), (1, 4, "chr1:4"), (1, 7, "chr1:7"), (1, 8, "chr1:8"), (1, 20, "chr1:20"),
                (3, 4, "chr3:4"), (3, 5, "chr3:5"), (3, 6, "chr3:6"), (3, 10, "chr3:10"),
                (4, 100, "chr4:100"),
                (5, 1, "chr5:1"), (5, 10, "chr5:10"),
                (3, 32, "chr3:32"), (3, 2, "chr3:2"), (3, 1, "chr3:1"),
                (5, 4, "chr5:4"), (5, 8, "chr5:8"),
            ] {
                write_bim_line(&mut writer, &chrom.to_string(), id, *coord, 'A', 'C');
            }
        }
        let mut bim = PlinkBim::new(file.into_temp_path().to_str().unwrap()).unwrap();
        let positions = bim.get_chrom_to_fileline_positions().unwrap();
        let expected: HashMap<String, OrderedIntegerSet<usize>> = vec![
            ("1".to_string(), OrderedIntegerSet::from(vec![ContiguousIntegerSet::new(0, 5)])),
            ("3".to_string(), OrderedIntegerSet::from(vec![
                ContiguousIntegerSet::new(6, 9),
                ContiguousIntegerSet::new(13, 15)
            ])),
            ("4".to_string(), OrderedIntegerSet::from(vec![ContiguousIntegerSet::new(10, 10)])),
            ("5".to_string(), OrderedIntegerSet::from(vec![
                ContiguousIntegerSet::new(11, 12),
                ContiguousIntegerSet::new(16, 17),
            ]))
        ].into_iter().collect();
        assert_eq!(positions, expected);
    }

    #[test]
    fn test_get_all_chroms() {
        let file = NamedTempFile::new().unwrap();
        {
            let mut writer = BufWriter::new(&file);
            for chrom in &[1, 2, 1, 3, 2, 5] {
                write_bim_line(&mut writer, &chrom.to_string(), "id", 0, 'A', 'C');
            }
        }
        let mut bim = PlinkBim::new(file.into_temp_path().to_str().unwrap()).unwrap();
        let chrom_set = bim.get_all_chroms().unwrap();
        let expected: HashSet<String> = vec!["1", "2", "3", "5"].into_iter().map(|x| x.to_string()).collect();
        assert_eq!(chrom_set, expected);
    }

    #[test]
    fn test_get_fileline_partitions() {
        let file = NamedTempFile::new().unwrap();
        {
            let mut writer = BufWriter::new(&file);
            for (chrom, coord, id) in &[
                (1, 2, "chr1:2"), (1, 3, "chr1:3"), (1, 4, "chr1:4"), (1, 7, "chr1:7"), (1, 8, "chr1:8"), (1, 20, "chr1:20"),
                (3, 4, "chr3:4"), (3, 5, "chr3:5"), (3, 6, "chr3:6"), (3, 10, "chr3:10"),
                (4, 100, "chr4:100"),
                (5, 1, "chr5:1"), (5, 10, "chr5:10"),
                (3, 32, "chr3:32"), (3, 2, "chr3:2"), (3, 1, "chr3:1"),
                (5, 4, "chr5:4"), (5, 8, "chr5:8"),
            ] {
                write_bim_line(&mut writer, &chrom.to_string(), id, *coord, 'A', 'C');
            }
        }
        let mut bim = PlinkBim::new(file.into_temp_path().to_str().unwrap()).unwrap();

        let partition_file = NamedTempFile::new().unwrap();
        {
            let mut writer = BufWriter::new(&partition_file);
            for (id, partition) in &[
                ("chr1:4", "p1"), ("chr1:7", "p2"), ("chr1:8", "p2"), ("chr1:20", "p1"),
                ("chr5:1", "p3"), ("chr5:10", "p1"),
                ("chr3:4", "p2"), ("chr3:5", "p3"), ("chr3:6", "p2"), ("chr3:10", "p4"),
                ("chr3:32", "p1"), ("chr3:2", "p1"), ("chr3:1", "p1"),
                ("chr5:4", "p3"), ("chr5:8", "p2"),
                ("chr4:100", "p1"),
                ("chr1:2", "p3"), ("chr1:3", "p1"),
            ] {
                writer.write_fmt(format_args!("{} {}\n", id, partition)).unwrap();
            }
        }
        let partition_file_path = partition_file.into_temp_path();
        let partitions = bim.get_fileline_partitions_from_partition_file(partition_file_path.to_str().unwrap()).unwrap();
        assert_eq!(partitions.get("p1").unwrap(), &OrderedIntegerSet::from_slice(&[[1, 2], [5, 5], [10, 10], [12, 15]]));
        assert_eq!(partitions.get("p2").unwrap(), &OrderedIntegerSet::from_slice(&[[3, 4], [6, 6], [8, 8], [17, 17]]));
        assert_eq!(partitions.get("p3").unwrap(), &OrderedIntegerSet::from_slice(&[[0, 0], [7, 7], [11, 11], [16, 16]]));
        assert_eq!(partitions.get("p4").unwrap(), &OrderedIntegerSet::from_slice(&[[9, 9]]));

        let mut new_bim = bim.into_partitioned_by_file(partition_file_path.to_str().unwrap()).unwrap();
        assert_eq!(new_bim.get_fileline_partitions_by_ref().unwrap(), &partitions);

        {
            let mut writer = BufWriter::new(
                OpenOptions::new().write(true).append(true).open(partition_file_path.to_str().unwrap()).unwrap());
            writer.write_fmt(format_args!("{} {}\n", "extra_id", "p2")).unwrap();
        }
        assert!(new_bim.get_fileline_partitions_from_partition_file(partition_file_path.to_str().unwrap()).is_err());
    }
}
