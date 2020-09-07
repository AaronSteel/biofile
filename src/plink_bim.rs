use std::{
    collections::{HashMap, HashSet},
    fs::{File, OpenOptions},
    io::{BufRead, BufReader},
    iter::FromIterator,
    slice::Iter,
};

use math::{
    partition::integer_partitions::Partition,
    set::ordered_integer_set::OrderedIntegerSet, traits::Collecting,
};

use crate::error::Error;

pub const CHROM_FIELD_INDEX: usize = 0;
pub const VARIANT_ID_FIELD_INDEX: usize = 1;
pub const COORDINATE_FIELD_INDEX: usize = 3;
pub const FIRST_ALLELE_FIELD_INDEX: usize = 4;
pub const SECOND_ALLELE_FIELD_INDEX: usize = 5;

pub type PartitionKey = String;
pub type Coordinate = i64;

pub struct PlinkBim {
    bim_path_list: Vec<String>,
    // maps partition_id to the file line indices
    fileline_partitions: Option<FilelinePartitions>,
}

impl PlinkBim {
    pub fn new(bim_path_list: Vec<String>) -> Result<PlinkBim, Error> {
        Ok(PlinkBim {
            bim_path_list,
            fileline_partitions: None,
        })
    }

    pub fn new_with_partitions(
        bim_path_list: Vec<String>,
        partitions: HashMap<PartitionKey, OrderedIntegerSet<Coordinate>>,
    ) -> Result<PlinkBim, Error> {
        let mut bim = PlinkBim::new(bim_path_list)?;
        bim.set_fileline_partitions(Some(FilelinePartitions::new(partitions)));
        Ok(bim)
    }

    pub fn new_with_partition_file(
        bim_path_list: Vec<String>,
        partition_filepath: &str,
    ) -> Result<PlinkBim, Error> {
        let bim = PlinkBim::new(bim_path_list)?;
        bim.into_partitioned_by_file(partition_filepath)
    }

    fn get_buf_list(&self) -> Result<Vec<BufReader<File>>, Error> {
        self.bim_path_list
            .iter()
            .map(|p| Ok(BufReader::new(OpenOptions::new().read(true).open(p)?)))
            .collect::<Result<Vec<BufReader<File>>, Error>>()
    }

    pub fn into_partitioned_by_file(
        mut self,
        partition_file: &str,
    ) -> Result<PlinkBim, Error> {
        let partitions =
            self.get_fileline_partitions_from_partition_file(partition_file)?;
        self.set_fileline_partitions(Some(partitions));
        Ok(self)
    }

    pub fn into_partitioned_by_chrom(mut self) -> Result<PlinkBim, Error> {
        let partitions = self.get_chrom_to_fileline_positions()?;
        self.set_fileline_partitions(Some(FilelinePartitions::new(partitions)));
        Ok(self)
    }

    #[inline]
    pub fn set_fileline_partitions(
        &mut self,
        partitions: Option<FilelinePartitions>,
    ) {
        self.fileline_partitions = partitions;
    }

    pub fn get_chrom_to_fileline_positions(
        &mut self,
    ) -> Result<HashMap<PartitionKey, OrderedIntegerSet<Coordinate>>, Error>
    {
        Ok(self
            .get_all_chroms()?
            .into_iter()
            .map(|chrom| {
                let partition = self.get_chrom_fileline_positions(&chrom)?;
                Ok((chrom, partition))
            })
            .collect::<Result<
                HashMap<PartitionKey, OrderedIntegerSet<Coordinate>>,
                Error,
            >>()?)
    }

    pub fn get_fileline_partitions(&self) -> Option<FilelinePartitions> {
        match &self.fileline_partitions {
            None => None,
            Some(p) => Some(p.clone()),
        }
    }

    pub fn get_fileline_partitions_by_ref(
        &self,
    ) -> Option<&FilelinePartitions> {
        match &self.fileline_partitions {
            None => None,
            Some(p) => Some(p),
        }
    }

    pub fn get_fileline_partitions_or(
        &self,
        default_partition_name: &str,
        default_partition_range: OrderedIntegerSet<Coordinate>,
    ) -> FilelinePartitions {
        match &self.fileline_partitions {
            None => FilelinePartitions::new(HashMap::from_iter(
                vec![(
                    default_partition_name.to_string(),
                    default_partition_range,
                )]
                .into_iter(),
            )),
            Some(p) => p.clone(),
        }
    }

    /// each line in the partition file has two fields separated by space:
    /// variant_id assigned_partition
    pub fn get_fileline_partitions_from_partition_file(
        &mut self,
        partition_file_path: &str,
    ) -> Result<FilelinePartitions, Error> {
        let id_to_partition_key: HashMap<String, PartitionKey> =
            BufReader::new(
                OpenOptions::new().read(true).open(partition_file_path)?,
            )
            .lines()
            .map(|line| PlinkBim::get_id_and_partition_from_line(&line?))
            .collect::<Result<HashMap<String, PartitionKey>, Error>>()?;

        let mut partitions: HashMap<
            PartitionKey,
            OrderedIntegerSet<Coordinate>,
        > = HashMap::new();
        let mut visited_ids = HashSet::new();
        let mut num_ignored_ids = 0;

        let mut file_line_offset = 0;
        for (b, buf) in self.get_buf_list()?.into_iter().enumerate() {
            let mut last_line_index = 0;
            for (i, line) in buf.lines().enumerate() {
                last_line_index = i as Coordinate;
                let id = match line?.split_whitespace().nth(VARIANT_ID_FIELD_INDEX) {
                    None => {
                        return Err(Error::BadFormat(format!(
                            "failed to read the variant id in line {} in bim file: {}",
                            i + 1,
                            self.bim_path_list[b]
                        )))
                    }
                    Some(id) => id.to_string(),
                };
                visited_ids.insert(id.to_string());
                match id_to_partition_key.get(&id) {
                    None => {
                        if !visited_ids.contains(&id) {
                            num_ignored_ids += 1;
                        }
                    }
                    Some(key) => partitions
                        .entry(key.to_owned())
                        .or_insert(OrderedIntegerSet::new())
                        .collect(file_line_offset + i as Coordinate),
                }
            }
            file_line_offset += last_line_index + 1;
        }

        if num_ignored_ids > 0 {
            println!("PlinkBim: num_ignored_ids: {}", num_ignored_ids);
        }
        let num_ids_not_in_bim = id_to_partition_key
            .keys()
            .filter(|&k| !visited_ids.contains(k))
            .count();
        if num_ids_not_in_bim > 0 {
            Err(Error::Generic(format!(
                "{} ID(s) from the partition file {} are not in the bim files {:?}",
                num_ids_not_in_bim, partition_file_path, self.bim_path_list
            )))
        } else {
            Ok(FilelinePartitions::new(partitions))
        }
    }

    fn get_id_and_partition_from_line(
        line: &str,
    ) -> Result<(String, PartitionKey), Error> {
        let mut iter = line.split_whitespace();
        let id = match iter.next() {
            None => return Err(Error::BadFormat(String::from(
                "each line in the partition file should have two fields: id partition_assignment",
            ))),
            Some(id) => id.to_string(),
        };
        match iter.next() {
            None => Err(Error::BadFormat(String::from(
                "each line in the partition file should have two fields: id partition_assignment",
            ))),
            Some(partition) => Ok((id, partition.to_string())),
        }
    }

    #[inline]
    pub fn get_bim_path_list(&self) -> &Vec<String> {
        &self.bim_path_list
    }

    pub fn get_all_chroms(&mut self) -> Result<HashSet<String>, Error> {
        Ok(self
            .get_buf_list()?
            .into_iter()
            .enumerate()
            .map(|(b, buf)| {
                buf.lines()
                    .enumerate()
                    .map(|(i, l)| {
                        Ok(l?
                            .split_whitespace()
                            .nth(CHROM_FIELD_INDEX)
                            .ok_or_else(|| {
                                Error::Generic(format!(
                                    "failed to extract chromosome id from line {} in file {}",
                                    i + 1,
                                    self.bim_path_list[b]
                                ))
                            })?
                            .to_string())
                    })
                    .collect::<Result<HashSet<String>, Error>>()
            })
            .collect::<Result<Vec<HashSet<String>>, Error>>()?
            .into_iter()
            .flatten()
            .collect())
    }

    pub fn get_chrom_fileline_positions(
        &mut self,
        chrom: &str,
    ) -> Result<OrderedIntegerSet<Coordinate>, Error> {
        let mut set = OrderedIntegerSet::new();
        let mut file_line_offset = 0;
        for (b, buf) in self.get_buf_list()?.into_iter().enumerate() {
            let mut last_line_index = 0;
            for (i, line) in buf.lines().enumerate() {
                last_line_index = i as Coordinate;
                if line?
                    .split_whitespace()
                    .nth(CHROM_FIELD_INDEX)
                    .ok_or_else(|| {
                        Error::Generic(format!(
                            "failed to extract chromosome id from line {} in file {}",
                            i + 1,
                            self.bim_path_list[b]
                        ))
                    })?
                    == chrom
                {
                    set.collect(file_line_offset + i as Coordinate);
                }
            }
            file_line_offset += last_line_index + 1;
        }
        Ok(set)
    }
}

#[derive(Clone, Debug)]
pub struct FilelinePartitions {
    partitions: HashMap<PartitionKey, Partition>,
    ordered_partition_keys: Vec<String>,
}

impl FilelinePartitions {
    pub fn new(
        partitions: HashMap<PartitionKey, Partition>,
    ) -> FilelinePartitions {
        let ordered_partition_keys =
            FilelinePartitions::get_ordered_keys(&partitions);
        FilelinePartitions {
            partitions,
            ordered_partition_keys,
        }
    }

    fn get_ordered_keys(
        partitions: &HashMap<PartitionKey, Partition>,
    ) -> Vec<PartitionKey> {
        let mut keys: Vec<PartitionKey> =
            partitions.keys().map(|s| s.to_string()).collect();
        if keys.iter().filter(|&k| k.parse::<i32>().is_err()).count() > 0 {
            keys.sort();
        } else {
            keys.sort_by_key(|k| k.parse::<i32>().unwrap());
        }
        keys
    }

    #[inline]
    pub fn ordered_partition_keys(&self) -> &Vec<PartitionKey> {
        &self.ordered_partition_keys
    }

    pub fn ordered_partition_array(&self) -> Vec<Partition> {
        self.iter().map(|(_, p)| p.clone()).collect()
    }

    #[inline]
    pub fn to_hash_map(&self) -> HashMap<PartitionKey, Partition> {
        self.partitions.clone()
    }

    #[inline]
    pub fn into_hash_map(self) -> HashMap<PartitionKey, Partition> {
        self.partitions
    }

    pub fn iter(&self) -> FilelinePartitionsIter {
        FilelinePartitionsIter {
            partitions: &self.partitions,
            key_iter: self.ordered_partition_keys.iter(),
        }
    }
}

pub struct FilelinePartitionsIter<'a> {
    partitions: &'a HashMap<String, OrderedIntegerSet<Coordinate>>,
    key_iter: Iter<'a, String>,
}

impl<'a> Iterator for FilelinePartitionsIter<'a> {
    type Item = (&'a str, &'a OrderedIntegerSet<Coordinate>);

    fn next(&mut self) -> Option<Self::Item> {
        match self.key_iter.next() {
            None => None,
            Some(key) => Some((key, &self.partitions[key])),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::plink_bim::{Coordinate, PlinkBim};
    use math::set::{
        contiguous_integer_set::ContiguousIntegerSet,
        ordered_integer_set::OrderedIntegerSet,
    };
    use std::{
        collections::{HashMap, HashSet},
        fs::OpenOptions,
        io::{BufWriter, Write},
    };
    use tempfile::NamedTempFile;

    #[inline]
    fn write_bim_line<W: Write>(
        buf_writer: &mut BufWriter<W>,
        chrom: &str,
        id: &str,
        coordinate: u64,
        first_allele: char,
        second_allele: char,
    ) {
        buf_writer
            .write_fmt(format_args!(
                "{} {} {} {} {} {}\n",
                chrom, id, 0, coordinate, first_allele, second_allele
            ))
            .unwrap();
    }

    #[test]
    fn test_get_chrom_to_positions() {
        let bim_file_1 = NamedTempFile::new().unwrap();
        {
            let mut writer = BufWriter::new(&bim_file_1);
            for (chrom, coord, id) in &[
                (1, 2, "chr1:2"),
                (1, 3, "chr1:3"),
                (1, 4, "chr1:4"),
                (1, 7, "chr1:7"),
                (1, 8, "chr1:8"),
                (1, 20, "chr1:20"),
                (3, 4, "chr3:4"),
                (3, 5, "chr3:5"),
                (3, 6, "chr3:6"),
                (3, 10, "chr3:10"),
                (4, 100, "chr4:100"),
                (5, 1, "chr5:1"),
                (5, 10, "chr5:10"),
                (3, 32, "chr3:32"),
                (3, 2, "chr3:2"),
                (3, 1, "chr3:1"),
                (5, 4, "chr5:4"),
                (5, 8, "chr5:8"),
            ] {
                write_bim_line(
                    &mut writer,
                    &chrom.to_string(),
                    id,
                    *coord,
                    'A',
                    'C',
                );
            }
        }
        let bim_file_2 = NamedTempFile::new().unwrap();
        {
            let mut writer = BufWriter::new(&bim_file_2);
            for (chrom, coord, id) in &[
                (1, 0, "chr1:0"),
                (1, 200, "chr1:200"),
                (5, 1000, "chr5:1000"),
                (12, 100, "chr12:100"),
            ] {
                write_bim_line(
                    &mut writer,
                    &chrom.to_string(),
                    id,
                    *coord,
                    'A',
                    'C',
                );
            }
        }
        let bim_temp_path_1 = bim_file_1.into_temp_path();
        let bim_temp_path_2 = bim_file_2.into_temp_path();
        let mut bim = PlinkBim::new(vec![
            bim_temp_path_1.to_str().unwrap().to_string(),
            bim_temp_path_2.to_str().unwrap().to_string(),
        ])
        .unwrap();
        let positions = bim.get_chrom_to_fileline_positions().unwrap();
        let expected: HashMap<String, OrderedIntegerSet<Coordinate>> = vec![
            (
                "1".to_string(),
                OrderedIntegerSet::from(vec![
                    ContiguousIntegerSet::new(0, 5),
                    ContiguousIntegerSet::new(18, 19),
                ]),
            ),
            (
                "3".to_string(),
                OrderedIntegerSet::from(vec![
                    ContiguousIntegerSet::new(6, 9),
                    ContiguousIntegerSet::new(13, 15),
                ]),
            ),
            (
                "4".to_string(),
                OrderedIntegerSet::from(vec![ContiguousIntegerSet::new(
                    10, 10,
                )]),
            ),
            (
                "5".to_string(),
                OrderedIntegerSet::from(vec![
                    ContiguousIntegerSet::new(11, 12),
                    ContiguousIntegerSet::new(16, 17),
                    ContiguousIntegerSet::new(20, 20),
                ]),
            ),
            (
                "12".to_string(),
                OrderedIntegerSet::from(vec![ContiguousIntegerSet::new(
                    21, 21,
                )]),
            ),
        ]
        .into_iter()
        .collect();
        assert_eq!(positions, expected);
    }

    #[test]
    fn test_get_all_chroms() {
        let file = NamedTempFile::new().unwrap();
        {
            let mut writer = BufWriter::new(&file);
            for chrom in &[1, 2, 1, 3, 2, 5] {
                write_bim_line(
                    &mut writer,
                    &chrom.to_string(),
                    "id",
                    0,
                    'A',
                    'C',
                );
            }
        }
        let bim_temp_path = file.into_temp_path();
        let mut bim =
            PlinkBim::new(vec![bim_temp_path.to_str().unwrap().to_string()])
                .unwrap();
        let chrom_set = bim.get_all_chroms().unwrap();
        let expected: HashSet<String> = vec!["1", "2", "3", "5"]
            .into_iter()
            .map(|x| x.to_string())
            .collect();
        assert_eq!(chrom_set, expected);
    }

    fn create_dummy_bim() -> (NamedTempFile, NamedTempFile) {
        let bim_temp_file = NamedTempFile::new().unwrap();
        {
            let mut writer = BufWriter::new(&bim_temp_file);
            for (chrom, coord, id) in &[
                (1, 2, "chr1:2"),
                (1, 3, "chr1:3"),
                (1, 4, "chr1:4"),
                (1, 7, "chr1:7"),
                (1, 8, "chr1:8"),
                (1, 20, "chr1:20"),
                (3, 4, "chr3:4"),
                (3, 5, "chr3:5"),
                (3, 6, "chr3:6"),
                (3, 10, "chr3:10"),
                (4, 100, "chr4:100"),
                (5, 1, "chr5:1"),
                (5, 10, "chr5:10"),
                (3, 32, "chr3:32"),
                (3, 2, "chr3:2"),
                (3, 1, "chr3:1"),
                (5, 4, "chr5:4"),
                (5, 8, "chr5:8"),
            ] {
                write_bim_line(
                    &mut writer,
                    &chrom.to_string(),
                    id,
                    *coord,
                    'A',
                    'C',
                );
            }
        }

        let partition_file = NamedTempFile::new().unwrap();
        {
            let mut writer = BufWriter::new(&partition_file);
            for (id, partition) in &[
                ("chr1:4", "p1"),
                ("chr1:7", "p2"),
                ("chr1:8", "p2"),
                ("chr1:20", "p1"),
                ("chr5:1", "p3"),
                ("chr5:10", "p1"),
                ("chr3:4", "p2"),
                ("chr3:5", "p3"),
                ("chr3:6", "p2"),
                ("chr3:10", "p4"),
                ("chr3:32", "p1"),
                ("chr3:2", "p1"),
                ("chr3:1", "p1"),
                ("chr5:4", "p3"),
                ("chr5:8", "p2"),
                ("chr4:100", "p1"),
                ("chr1:2", "p3"),
                ("chr1:3", "p1"),
            ] {
                writer
                    .write_fmt(format_args!("{} {}\n", id, partition))
                    .unwrap();
            }
        }
        (partition_file, bim_temp_file)
    }

    #[test]
    fn test_get_fileline_partitions() {
        let (partition_file, bim_temp_file) = create_dummy_bim();
        let bim_temp_path = bim_temp_file.into_temp_path();
        let mut bim =
            PlinkBim::new(vec![bim_temp_path.to_str().unwrap().to_string()])
                .unwrap();

        let partition_file_path = partition_file.into_temp_path();
        let partitions = bim
            .get_fileline_partitions_from_partition_file(
                partition_file_path.to_str().unwrap(),
            )
            .unwrap()
            .partitions;
        assert_eq!(
            partitions.get("p1").unwrap(),
            &OrderedIntegerSet::from_slice(&[[1, 2], [5, 5], [10, 10], [
                12, 15
            ]])
        );
        assert_eq!(
            partitions.get("p2").unwrap(),
            &OrderedIntegerSet::from_slice(&[[3, 4], [6, 6], [8, 8], [17, 17]])
        );
        assert_eq!(
            partitions.get("p3").unwrap(),
            &OrderedIntegerSet::from_slice(&[[0, 0], [7, 7], [11, 11], [
                16, 16
            ]])
        );
        assert_eq!(
            partitions.get("p4").unwrap(),
            &OrderedIntegerSet::from_slice(&[[9, 9]])
        );

        let mut new_bim = bim
            .into_partitioned_by_file(partition_file_path.to_str().unwrap())
            .unwrap();
        assert_eq!(
            new_bim.get_fileline_partitions_by_ref().unwrap().partitions,
            partitions
        );

        {
            let mut writer = BufWriter::new(
                OpenOptions::new()
                    .write(true)
                    .append(true)
                    .open(partition_file_path.to_str().unwrap())
                    .unwrap(),
            );
            writer
                .write_fmt(format_args!("{} {}\n", "extra_id", "p2"))
                .unwrap();
        }
        assert!(new_bim
            .get_fileline_partitions_from_partition_file(
                partition_file_path.to_str().unwrap()
            )
            .is_err());
    }

    #[test]
    fn test_fileline_partitions_iter() {
        let (partition_file, bim_temp_file) = create_dummy_bim();
        let bim_temp_path = bim_temp_file.into_temp_path();
        let mut bim =
            PlinkBim::new(vec![bim_temp_path.to_str().unwrap().to_string()])
                .unwrap();

        let partition_file_path = partition_file.into_temp_path();
        let partitions = bim
            .get_fileline_partitions_from_partition_file(
                partition_file_path.to_str().unwrap(),
            )
            .unwrap();
        let mut iter = partitions.iter();
        assert_eq!(
            iter.next(),
            Some((
                "p1",
                &OrderedIntegerSet::from_slice(&[[1, 2], [5, 5], [10, 10], [
                    12, 15
                ]])
            ))
        );
        assert_eq!(
            iter.next(),
            Some((
                "p2",
                &OrderedIntegerSet::from_slice(&[[3, 4], [6, 6], [8, 8], [
                    17, 17
                ]])
            ))
        );
        assert_eq!(
            iter.next(),
            Some((
                "p3",
                &OrderedIntegerSet::from_slice(&[[0, 0], [7, 7], [11, 11], [
                    16, 16
                ]])
            ))
        );
        assert_eq!(
            iter.next(),
            Some(("p4", &OrderedIntegerSet::from_slice(&[[9, 9]])))
        );
        assert_eq!(iter.next(), None);
    }
}
