use std::cmp::min;
use std::fs::{File, OpenOptions};
use std::io;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};

use analytic::set::ordered_integer_set::OrderedIntegerSet;
use analytic::set::traits::{Finite, Set};
use analytic::stats::sum_f32;
use analytic::traits::ToIterator;
use ndarray::{Array, Ix2, ShapeBuilder};
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use rayon::iter::plumbing::{bridge, Consumer, Producer, ProducerCallback, UnindexedConsumer};

use crate::byte_chunk_iter::ByteChunkIter;
use crate::error::Error;
use crate::util::get_buf;

const NUM_PEOPLE_PER_BYTE: usize = 4;

#[inline]
fn get_num_people_last_byte(total_num_people: usize) -> Option<usize> {
    if total_num_people == 0 {
        None
    } else {
        match total_num_people % NUM_PEOPLE_PER_BYTE {
            0 => Some(NUM_PEOPLE_PER_BYTE),
            x => Some(x)
        }
    }
}

pub struct PlinkBed {
    bed_buf: BufReader<File>,
    pub num_people: usize,
    pub num_snps: usize,
    pub num_bytes_per_snp: usize,
    pub filepath: String,
}

impl PlinkBed {
    #[inline]
    pub fn get_magic_bytes() -> [u8; 3] {
        [0x6c_u8, 0x1b_u8, 0x01_u8]
    }

    #[inline]
    pub fn get_num_magic_bytes() -> usize {
        PlinkBed::get_magic_bytes().len()
    }

    fn get_line_count(filename: &str) -> Result<usize, Error> {
        let fam_buf = get_buf(filename)?;
        Ok(fam_buf.lines().count())
    }

    #[inline]
    fn usize_div_ceil(a: usize, divisor: usize) -> usize {
        a / divisor + (a % divisor != 0) as usize
    }

    pub fn new(bed_filename: &str, bim_filename: &str, fam_filename: &str) -> Result<PlinkBed, Error> {
        let mut bed_buf = get_buf(bed_filename)?;

        // check if PLINK bed file has the correct file signature
        let expected_bytes = PlinkBed::get_magic_bytes();
        let mut magic_bytes = [0u8; 3];
        if let Err(io_error) = bed_buf.read_exact(&mut magic_bytes) {
            return Err(
                Error::IO {
                    why: format!("Failed to read the first three bytes of {}: {}", bed_filename, io_error),
                    io_error,
                }
            );
        }
        if magic_bytes != expected_bytes {
            return Err(Error::BadFormat(
                format!("The first three bytes of the PLINK bed file are supposed to be 0x{:x?}, but found 0x{:x?}",
                        expected_bytes, magic_bytes)
            ));
        }

        let num_people = PlinkBed::get_line_count(fam_filename)?;
        let num_snps = PlinkBed::get_line_count(bim_filename)?;
        let num_bytes_per_snp = PlinkBed::usize_div_ceil(num_people, 4);
        println!("{} stats:\nnum_snps: {}\nnum_people: {}\nnum_bytes_per_block: {}\n----------",
                 bed_filename, num_snps, num_people, num_bytes_per_snp);

        if num_snps == 0 {
            return Err(Error::Generic("cannot create PlinkBed with 0 SNPs".to_string()));
        }
        if num_people == 0 {
            return Err(Error::Generic("cannot create PlinkBed with 0 people".to_string()));
        }
        Ok(PlinkBed { bed_buf, num_people, num_snps, num_bytes_per_snp, filepath: bed_filename.to_string() })
    }

    // the first person is the lowest two bits
    // 00 -> 2 homozygous for the first allele in the .bim file (usually the minor allele)
    // 01 -> 0 missing genotype
    // 10 -> 1 heterozygous
    // 11 -> 0 homozygous for the second allele in the .bim file (usually the major allele)
    pub fn create_bed(arr: &Array<u8, Ix2>, out_path: &str) -> Result<(), Error> {
        let (num_people, _num_snps) = arr.dim();
        let mut buf_writer = BufWriter::new(OpenOptions::new().create(true).truncate(true).write(true).open(out_path)?);
        buf_writer.write(&[0x6c, 0x1b, 0x1])?;
        for col in arr.gencolumns() {
            let mut i = 0;
            for _ in 0..num_people / 4 {
                buf_writer.write(&[
                    PlinkBed::geno_to_lowest_two_bits(col[i])
                        | (PlinkBed::geno_to_lowest_two_bits(col[i + 1]) << 2)
                        | (PlinkBed::geno_to_lowest_two_bits(col[i + 2]) << 4)
                        | (PlinkBed::geno_to_lowest_two_bits(col[i + 3]) << 6)
                ])?;
                i += 4;
            }
            let remainder = num_people % 4;
            if remainder > 0 {
                let mut byte = 0u8;
                for j in 0..remainder {
                    byte |= PlinkBed::geno_to_lowest_two_bits(col[i + j]) << (j * 2);
                }
                buf_writer.write(&[byte])?;
            }
        }
        Ok(())
    }

    pub fn reset_bed_buf(&mut self) -> Result<(), io::Error> {
        // the first three bytes are the file signature
        self.bed_buf.seek(SeekFrom::Start(3))?;
        Ok(())
    }

    #[inline]
    fn set_buf_to_snp_start(buf: &mut BufReader<File>) -> io::Result<u64> {
        buf.seek(SeekFrom::Start(3))
    }

    #[inline]
    fn get_buf_pointed_to_first_snp(&self) -> Result<BufReader<File>, Error> {
        let mut buf = get_buf(&self.filepath)?;
        PlinkBed::set_buf_to_snp_start(&mut buf)?;
        Ok(buf)
    }

    /// makes the BufReader point to the start of the byte containing the SNP i individual j
    /// 0-indexing
    fn seek_to_byte_containing_snp_i_person_j(&mut self, snp_i: usize, person_j: usize) -> Result<(), io::Error> {
        // the first three bytes are the file signature
        self.bed_buf.seek(SeekFrom::Start((3 + self.num_bytes_per_snp * snp_i + person_j / 4) as u64))?;
        Ok(())
    }

    pub fn geno_to_lowest_two_bits(geno: u8) -> u8 {
        // 00 -> 2 homozygous for the first allele in the .bim file (usually the minor allele)
        // 01 -> 0 missing genotype
        // 10 -> 1 heterozygous
        // 11 -> 0 homozygous for the second allele in the .bim file (usually the major allele)
        let not_a = ((geno & 0b10) >> 1) ^ 1;
        let not_b = (geno & 1) ^ 1;
        (not_a << 1) | (not_b & not_a)
    }

    pub fn lowest_two_bits_to_geno(byte: u8) -> u8 {
        // 00 -> 2 homozygous for the first allele in the .bim file (usually the minor allele)
        // 01 -> 0 missing genotype
        // 10 -> 1 heterozygous
        // 11 -> 0 homozygous for the second allele in the .bim file (usually the major allele)
        let a = (byte & 0b10) >> 1;
        let b = byte & 1;
        (((a | b) ^ 1) << 1) | (a & (!b))
    }

    pub fn get_genotype_matrix(
        &self,
        snps_range: Option<OrderedIntegerSet<usize>>,
    ) -> Result<Array<f32, Ix2>, Error> {
        let mut bed_buf = self.get_buf_pointed_to_first_snp()?;

        let num_snps = match &snps_range {
            None => self.num_snps,
            Some(range) => range.size()
        };
        let mut v = Vec::with_capacity(self.num_people * num_snps);
        let mut vi = 0usize;

        if let Some(num_people_last_byte) = get_num_people_last_byte(self.num_people) {
            let last_byte_index = self.num_bytes_per_snp - 1;
            let mut snp_bytes = vec![0u8; self.num_bytes_per_snp];
            match snps_range {
                None => {
                    unsafe {
                        v.set_len(self.num_people * num_snps);
                    }
                    for _ in 0..self.num_snps {
                        bed_buf.read_exact(&mut snp_bytes)?;
                        for i in 0..last_byte_index {
                            v[vi] = PlinkBed::lowest_two_bits_to_geno(snp_bytes[i]) as f32;
                            v[vi + 1] = PlinkBed::lowest_two_bits_to_geno(snp_bytes[i] >> 2) as f32;
                            v[vi + 2] = PlinkBed::lowest_two_bits_to_geno(snp_bytes[i] >> 4) as f32;
                            v[vi + 3] = PlinkBed::lowest_two_bits_to_geno(snp_bytes[i] >> 6) as f32;
                            vi += 4;
                        }
                        // last byte
                        for k in 0..num_people_last_byte {
                            v[vi] = PlinkBed::lowest_two_bits_to_geno(snp_bytes[last_byte_index] >> (k << 1)) as f32;
                            vi += 1;
                        }
                    }
                }
                Some(range) => {
                    for snp_chunk in self.col_chunk_iter(1000, Some(range)) {
                        v.extend_from_slice(snp_chunk.t().to_owned().as_slice().unwrap());
                    }
                }
            };
        }
        let geno_arr = Array::from_shape_vec(
            (self.num_people, num_snps).strides((1, self.num_people)),
            v,
        ).unwrap();
        Ok(geno_arr)
    }

    pub fn col_chunk_iter(
        &self,
        num_snps_per_iter: usize,
        range: Option<OrderedIntegerSet<usize>>,
    ) -> PlinkColChunkIter {
        let buf = get_buf(&self.filepath).unwrap();
        match range {
            Some(range) => PlinkColChunkIter::new(
                buf,
                range,
                num_snps_per_iter,
                self.num_people,
                &self.filepath,
            ),
            None => {
                PlinkColChunkIter::new(
                    buf,
                    OrderedIntegerSet::from_slice(&[[0, self.num_snps - 1]]),
                    num_snps_per_iter,
                    self.num_people,
                    &self.filepath,
                )
            }
        }
    }

    pub fn get_minor_allele_frequencies(&self, chunk_size: Option<usize>) -> Vec<f32> {
        let num_alleles = (self.num_people * 2) as f32;
        self.col_chunk_iter(chunk_size.unwrap_or(50), None)
            .into_par_iter()
            .flat_map(|snps| {
                snps.gencolumns()
                    .into_iter()
                    .map(|col| sum_f32(col.iter()) / num_alleles)
                    .collect::<Vec<f32>>()
            })
            .collect()
    }

    pub fn byte_chunk_iter(
        &self,
        start_byte_index: usize,
        end_byte_index_exclusive: usize,
        chunk_size: usize,
    ) -> Result<ByteChunkIter<File>, Error> {
        let buf = match OpenOptions::new().read(true).open(&self.filepath) {
            Ok(file) => BufReader::new(file),
            Err(io_error) => return Err(
                Error::IO {
                    why: format!("failed to open {}: {}", self.filepath, io_error),
                    io_error,
                }
            ),
        };
        Ok(ByteChunkIter::new(buf, start_byte_index, end_byte_index_exclusive, chunk_size))
    }

    /// save the transpose of the BED file into `out_path`, which should have an extension of .bedt
    /// wherein the n-th sequence of bytes corresponds to the SNPs for the n-th person
    /// larger values of `snp_byte_chunk_size` lead to faster performance, at the cost of higher memory requirement
    pub fn create_bed_t(&mut self, out_path: &str, snp_byte_chunk_size: usize) -> Result<(), io::Error> {
        let mut buf_writer = BufWriter::new(
            OpenOptions::new().create(true).truncate(true).write(true).open(out_path)?
        );
        let num_bytes_per_person = PlinkBed::usize_div_ceil(self.num_snps, 4);

        let people_stride = snp_byte_chunk_size * 4;
        let mut snp_bytes = vec![0u8; snp_byte_chunk_size];

        // write people_stride people at a time
        for j in (0..self.num_people).step_by(people_stride) {
            let mut people_buf = vec![vec![0u8; num_bytes_per_person]; people_stride];
            if self.num_people - j < people_stride {
                let remaining_people = self.num_people % people_stride;
                snp_bytes = vec![0u8; PlinkBed::usize_div_ceil(remaining_people, 4)];
            }
            let relative_seek_offset = (self.num_bytes_per_snp - snp_bytes.len()) as i64;
            // read 4 SNPs to the buffers at a time
            self.seek_to_byte_containing_snp_i_person_j(0, j)?;
            for (snp_byte_index, k) in (0..self.num_snps).step_by(4).enumerate() {
                for (snp_offset, _) in (k..min(k + 4, self.num_snps)).enumerate() {
                    self.bed_buf.read_exact(&mut snp_bytes)?;
                    for w in 0..snp_bytes.len() {
                        people_buf[w + 0][snp_byte_index] |= (snp_bytes[w] & 0b11) << (snp_offset << 1);
                        people_buf[w + 1][snp_byte_index] |= ((snp_bytes[w] >> 2) & 0b11) << (snp_offset << 1);
                        people_buf[w + 2][snp_byte_index] |= ((snp_bytes[w] >> 4) & 0b11) << (snp_offset << 1);
                        people_buf[w + 3][snp_byte_index] |= ((snp_bytes[w] >> 6) & 0b11) << (snp_offset << 1);
                    }
                    self.bed_buf.seek_relative(relative_seek_offset)?;
                }
            }
            for (p, buf) in people_buf.iter().enumerate() {
                if j + p < self.num_people {
                    buf_writer.write(buf.as_slice())?;
                }
            }
        }
        Ok(())
    }

    pub fn create_dominance_geno_bed(&self, out_path: &str) -> Result<(), Error> {
        let mut writer = BufWriter::new(
            OpenOptions::new().create(true).truncate(true).write(true).open(out_path)?
        );
        writer.write(&PlinkBed::get_magic_bytes())?;
        let num_magic_bytes = PlinkBed::get_num_magic_bytes();
        for bytes in self.byte_chunk_iter(
            num_magic_bytes,
            num_magic_bytes + self.num_snps * self.num_bytes_per_snp,
            self.num_bytes_per_snp,
        )? {
            let out_bytes = PlinkSnps::from_geno(
                PlinkSnps::new(bytes, self.num_people)
                    .into_iter()
                    .map(|s| match s {
                        2 => 1,
                        s => s,
                    })
                    .collect()
            ).into_bytes();
            writer.write(&out_bytes)?;
        }
        Ok(())
    }
}

pub struct PlinkSnps {
    bytes: Vec<u8>,
    num_snps: usize,
}

impl PlinkSnps {
    pub fn new(bytes: Vec<u8>, num_snps: usize) -> PlinkSnps {
        assert!(num_snps <= bytes.len() * 4,
                format!("num_snps ({}) > bytes.len() * 4 = {}", num_snps, bytes.len() * 4));
        PlinkSnps {
            bytes,
            num_snps,
        }
    }

    pub fn from_geno(geno: Vec<u8>) -> PlinkSnps {
        let num_bytes = PlinkBed::usize_div_ceil(geno.len(), 4);
        let mut bytes: Vec<u8> = Vec::with_capacity(num_bytes);
        let mut snp_index = 0;
        if let Some(num_people_last_byte) = get_num_people_last_byte(geno.len()) {
            for _ in 0..num_bytes - 1 {
                bytes.push(
                    PlinkBed::geno_to_lowest_two_bits(geno[snp_index])
                        | (PlinkBed::geno_to_lowest_two_bits(geno[snp_index + 1]) << 2)
                        | (PlinkBed::geno_to_lowest_two_bits(geno[snp_index + 2]) << 4)
                        | (PlinkBed::geno_to_lowest_two_bits(geno[snp_index + 3]) << 6)
                );
                snp_index += 4;
            }
            // last byte
            let mut last_byte = 0u8;
            for j in 0..num_people_last_byte {
                last_byte |= PlinkBed::geno_to_lowest_two_bits(geno[snp_index + j]) << (j * 2);
            }
            bytes.push(last_byte);
        }
        PlinkSnps::new(bytes, geno.len())
    }

    #[inline]
    pub fn get_num_snps(&self) -> usize {
        self.num_snps
    }

    #[inline]
    pub fn to_bytes(&self) -> &Vec<u8> {
        &self.bytes
    }

    #[inline]
    pub fn into_bytes(self) -> Vec<u8> {
        self.bytes
    }
}

impl ToIterator<'_, PlinkSnpsIter, u8> for PlinkSnps {
    fn to_iter(&self) -> PlinkSnpsIter {
        PlinkSnpsIter {
            bytes: self.bytes.clone(),
            num_snps: self.num_snps,
            byte_cursor: 0,
            cursor: 0,
        }
    }
}

impl IntoIterator for PlinkSnps {
    type Item = <PlinkSnpsIter as Iterator>::Item;
    type IntoIter = PlinkSnpsIter;
    fn into_iter(self) -> Self::IntoIter {
        PlinkSnpsIter {
            bytes: self.bytes,
            num_snps: self.num_snps,
            byte_cursor: 0,
            cursor: 0,
        }
    }
}

pub struct PlinkSnpsIter {
    bytes: Vec<u8>,
    num_snps: usize,
    byte_cursor: usize,
    cursor: usize,
}

impl Iterator for PlinkSnpsIter {
    type Item = u8;
    fn next(&mut self) -> Option<Self::Item> {
        if self.cursor >= self.num_snps {
            None
        } else {
            let snp = PlinkBed::lowest_two_bits_to_geno(self.bytes[self.byte_cursor] >> (2 * (self.cursor % 4) as u8));
            self.cursor += 1;
            if self.cursor % 4 == 0 {
                self.byte_cursor += 1;
            }
            Some(snp)
        }
    }
}

pub struct PlinkColChunkIter {
    buf: BufReader<File>,
    range: OrderedIntegerSet<usize>,
    num_snps_per_iter: usize,
    num_people: usize,
    num_snps_in_range: usize,
    range_cursor: usize,
    last_read_snp_index: Option<usize>,
    bed_filename: String,
}

impl PlinkColChunkIter {
    pub fn new(buf: BufReader<File>,
               range: OrderedIntegerSet<usize>,
               num_snps_per_iter: usize,
               num_people: usize,
               bed_filename: &str,
    ) -> PlinkColChunkIter {
        let num_snps_in_range = range.size();
        let first = range.first();
        let mut iter = PlinkColChunkIter {
            buf,
            range,
            num_snps_per_iter,
            num_people,
            num_snps_in_range,
            range_cursor: 0,
            last_read_snp_index: None,
            bed_filename: bed_filename.to_string(),
        };
        if let Some(start) = first {
            iter.seek_to_snp(start).unwrap();
        }
        iter
    }

    #[inline]
    fn num_bytes_per_snp(&self) -> usize {
        PlinkBed::usize_div_ceil(self.num_people, NUM_PEOPLE_PER_BYTE)
    }

    fn seek_to_snp(&mut self, snp_index: usize) -> Result<(), Error> {
        if !self.range.contains(snp_index) {
            return Err(Error::Generic(format!("SNP index {} is not in the interator range", snp_index)));
        }
        // skip the first 3 magic bytes
        self.buf.seek(SeekFrom::Start(3 + (self.num_bytes_per_snp() * snp_index) as u64)).unwrap();
        Ok(())
    }

    /// indices are 0 based
    #[inline]
    fn clone_with_range(&self, range: OrderedIntegerSet<usize>) -> PlinkColChunkIter {
        PlinkColChunkIter::new(
            get_buf(&self.bed_filename).unwrap(),
            range,
            self.num_snps_per_iter,
            self.num_people,
            &self.bed_filename,
        )
    }

    fn read_chunk(&mut self, chunk_size: usize) -> Array<f32, Ix2> {
        let num_bytes_per_snp = self.num_bytes_per_snp();
        let num_people_last_byte = get_num_people_last_byte(self.num_people).unwrap_or(0);

        let snp_indices = self.range.slice(self.range_cursor..self.range_cursor + chunk_size);
        self.range_cursor += chunk_size;

        let mut v = Vec::with_capacity(self.num_people * chunk_size);
        unsafe {
            v.set_len(self.num_people * chunk_size);
        }
        let mut acc_i = 0usize;

        let mut snp_bytes = vec![0u8; num_bytes_per_snp];
        for index in snp_indices.to_iter() {
            if let Some(last_read_snp_index) = self.last_read_snp_index {
                let snp_index_gap = index - last_read_snp_index;
                if snp_index_gap > 1 {
                    self.buf.seek_relative(((snp_index_gap - 1) * self.num_bytes_per_snp()) as i64).unwrap();
                }
            }
            self.buf.read_exact(&mut snp_bytes).unwrap();
            for i in 0..num_bytes_per_snp - 1 {
                v[acc_i] = PlinkBed::lowest_two_bits_to_geno(snp_bytes[i]) as f32;
                v[acc_i + 1] = PlinkBed::lowest_two_bits_to_geno(snp_bytes[i] >> 2) as f32;
                v[acc_i + 2] = PlinkBed::lowest_two_bits_to_geno(snp_bytes[i] >> 4) as f32;
                v[acc_i + 3] = PlinkBed::lowest_two_bits_to_geno(snp_bytes[i] >> 6) as f32;
                acc_i += 4;
            }
            // last byte
            for k in 0..num_people_last_byte {
                v[acc_i] = PlinkBed::lowest_two_bits_to_geno(snp_bytes[num_bytes_per_snp - 1] >> (k << 1)) as f32;
                acc_i += 1;
            }
            self.last_read_snp_index = Some(index);
        }
        Array::from_shape_vec((self.num_people, chunk_size)
                                  .strides((1, self.num_people)), v).unwrap()
    }
}

impl IntoParallelIterator for PlinkColChunkIter {
    type Iter = PlinkColChunkParallelIter;
    type Item = <PlinkColChunkParallelIter as ParallelIterator>::Item;

    fn into_par_iter(self) -> Self::Iter {
        PlinkColChunkParallelIter { iter: self }
    }
}

impl Iterator for PlinkColChunkIter {
    type Item = Array<f32, Ix2>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.range_cursor >= self.num_snps_in_range {
            return None;
        }
        let chunk_size = min(self.num_snps_per_iter, self.num_snps_in_range - self.range_cursor);
        Some(self.read_chunk(chunk_size))
    }
}

impl ExactSizeIterator for PlinkColChunkIter {
    fn len(&self) -> usize {
        PlinkBed::usize_div_ceil(self.num_snps_in_range - self.range_cursor, self.num_snps_per_iter)
    }
}

impl DoubleEndedIterator for PlinkColChunkIter {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.range_cursor >= self.num_snps_in_range {
            return None;
        }
        let chunk_size = min(self.num_snps_per_iter, self.num_snps_in_range - self.range_cursor);
        // reading from the back is equivalent to reducing the number of SNPs in range
        self.num_snps_in_range -= chunk_size;

        // save and restore self.last_read_snp_index after the call to self.read_chunk
        // we set the self.last_read_snp_index to None to prevent self.read_chunk from performing seek_relative on the buffer
        let last_read_snp_index = self.last_read_snp_index;
        self.last_read_snp_index = None;

        let snp = self.range.slice(self.num_snps_in_range..self.num_snps_in_range + 1).first().unwrap();
        self.seek_to_snp(snp).unwrap();
        let chunk = self.read_chunk(chunk_size);
        self.seek_to_snp(last_read_snp_index.unwrap_or(0)).unwrap();

        self.last_read_snp_index = last_read_snp_index;
        Some(chunk)
    }
}

struct ColChunkIterProducer {
    iter: PlinkColChunkIter,
}

impl Producer for ColChunkIterProducer {
    type Item = <PlinkColChunkIter as Iterator>::Item;
    type IntoIter = PlinkColChunkIter;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.iter
    }

    fn split_at(self, index: usize) -> (Self, Self) {
        let mid_range_index = min(self.iter.num_snps_per_iter * index, self.iter.range.size());
        (
            ColChunkIterProducer {
                iter: self.iter.clone_with_range(self.iter.range.slice(0..mid_range_index))
            },
            ColChunkIterProducer {
                iter: self.iter.clone_with_range(self.iter.range.slice(mid_range_index..self.iter.range.size()))
            }
        )
    }
}

impl IntoIterator for ColChunkIterProducer {
    type Item = <PlinkColChunkIter as Iterator>::Item;
    type IntoIter = PlinkColChunkIter;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.iter
    }
}

pub struct PlinkColChunkParallelIter {
    iter: PlinkColChunkIter
}

impl ParallelIterator for PlinkColChunkParallelIter {
    type Item = <PlinkColChunkIter as Iterator>::Item;

    fn drive_unindexed<C>(self, consumer: C) -> C::Result
        where C: UnindexedConsumer<Self::Item>
    {
        bridge(self, consumer)
    }

    fn opt_len(&self) -> Option<usize> {
        Some(self.iter.len())
    }
}

impl IndexedParallelIterator for PlinkColChunkParallelIter {
    fn len(&self) -> usize {
        self.iter.len()
    }

    fn drive<C>(self, consumer: C) -> C::Result
        where C: Consumer<Self::Item>,
    {
        bridge(self, consumer)
    }

    fn with_producer<CB>(self, callback: CB) -> CB::Output
        where CB: ProducerCallback<Self::Item>,
    {
        callback.callback(ColChunkIterProducer { iter: self.iter })
    }
}

#[cfg(test)]
mod tests {
    use std::cmp::min;
    use std::io;
    use std::io::Write;

    use analytic::set::ordered_integer_set::OrderedIntegerSet;
    use analytic::traits::ToIterator;
    use ndarray::{array, Array, Ix2, s};
    use ndarray_rand::RandomExt;
    use rand::distributions::Uniform;
    use tempfile::NamedTempFile;

    use super::PlinkBed;
    use crate::plink_bed::PlinkSnps;

    fn create_dummy_bim_fam(
        bim: &mut NamedTempFile,
        fam: &mut NamedTempFile,
        num_people: usize,
        num_snps: usize,
    ) -> Result<(), io::Error> {
        for i in 1..=num_people {
            fam.write_fmt(format_args!("{}\n", i))?;
        }
        for i in 1..=num_snps {
            bim.write_fmt(format_args!("{}\n", i))?;
        }
        Ok(())
    }

    #[test]
    fn test_create_bed() {
        fn test(geno: &Array<u8, Ix2>) {
            let mut bim = NamedTempFile::new().unwrap();
            let mut fam = NamedTempFile::new().unwrap();
            create_dummy_bim_fam(&mut bim, &mut fam, geno.dim().0, geno.dim().1).unwrap();
            let path = NamedTempFile::new().unwrap().into_temp_path().to_str().unwrap().to_string();
            PlinkBed::create_bed(&geno, &path).unwrap();
            let geno_bed = PlinkBed::new(&path,
                                         bim.into_temp_path().to_str().unwrap(),
                                         fam.into_temp_path().to_str().unwrap()).unwrap();
            assert_eq!(geno.mapv(|x| x as f32), geno_bed.get_genotype_matrix(None).unwrap());
        }
        test(&array![
            [0],
        ]);
        test(&array![
            [1],
        ]);
        test(&array![
            [2],
        ]);
        test(&array![
            [0, 1, 2],
        ]);
        test(&array![
            [0],
            [1],
            [2],
        ]);
        test(&array![
            [0, 0, 1],
            [1, 1, 2],
            [0, 2, 1],
        ]);
        test(&array![
            [0, 0, 1, 2],
            [1, 1, 2, 1],
            [2, 0, 0, 0],
            [1, 0, 0, 2],
            [0, 2, 1, 0],
        ]);
        test(&array![
            [0, 0, 1, 2, 1],
            [1, 1, 2, 1, 2],
            [2, 0, 0, 0, 0],
            [1, 0, 0, 2, 2],
            [0, 2, 1, 0, 1],
        ]);
        test(&array![
            [0, 0, 1, 2, 1],
            [1, 0, 0, 2, 1],
            [2, 0, 2, 0, 0],
            [1, 1, 0, 2, 2],
            [0, 2, 2, 1, 1],
            [2, 1, 2, 0, 0],
            [1, 2, 0, 1, 2],
            [2, 0, 1, 0, 1],
        ]);
        test(&array![
            [0, 0, 1, 2, 1, 2, 2, 0],
            [1, 0, 0, 2, 1, 2, 1, 1],
            [2, 0, 2, 0, 0, 0, 2, 1],
            [1, 1, 0, 2, 2, 1, 1, 1],
            [0, 2, 2, 1, 1, 2, 0, 2],
            [2, 1, 2, 0, 0, 0, 2, 2],
            [1, 2, 0, 1, 2, 1, 1, 0],
            [2, 0, 1, 0, 1, 0, 0, 2],
        ]);
    }

    #[test]
    fn test_chunk_iter() {
        let (num_people, num_snps) = (137usize, 71usize);
        let geno = Array::random((num_people, num_snps), Uniform::from(0..3));

        let mut bim = NamedTempFile::new().unwrap();
        let mut fam = NamedTempFile::new().unwrap();
        create_dummy_bim_fam(&mut bim, &mut fam, num_people, num_snps).unwrap();
        let bed_file = NamedTempFile::new().unwrap();
        let bed_path = bed_file.into_temp_path().to_str().unwrap().to_string();
        PlinkBed::create_bed(&geno, &bed_path).unwrap();

        let bed = PlinkBed::new(
            &bed_path,
            bim.into_temp_path().to_str().unwrap(),
            fam.into_temp_path().to_str().unwrap(),
        ).unwrap();
        let true_geno_arr = bed.get_genotype_matrix(None).unwrap();
        let chunk_size = 5;
        for (i, snps) in bed.col_chunk_iter(chunk_size, None).enumerate() {
            let end_index = min((i + 1) * chunk_size, true_geno_arr.dim().1);
            assert!(true_geno_arr.slice(s![..,i * chunk_size..end_index]) == snps);
        }

        let snp_index_slices = OrderedIntegerSet::from_slice(&[[2, 4], [6, 9], [20, 46], [70, 70]]);
        for (i, snps) in bed.col_chunk_iter(chunk_size, Some(snp_index_slices.clone())).enumerate() {
            let end_index = min((i + 1) * chunk_size, true_geno_arr.dim().1);
            let snp_indices = snp_index_slices.slice(i * chunk_size..end_index);
            for (k, j) in snp_indices.to_iter().enumerate() {
                assert_eq!(true_geno_arr.slice(s![.., j]), snps.slice(s![.., k]));
            }
        }
        let geno = bed.get_genotype_matrix(Some(snp_index_slices.clone())).unwrap();
        let mut arr = Array::zeros((num_people, 35));
        let mut jj = 0;
        for j in snp_index_slices.to_iter() {
            for i in 0..num_people {
                arr[[i, jj]] = true_geno_arr[[i, j]];
            }
            jj += 1;
        }
        assert_eq!(arr, geno);
    }

    #[test]
    fn test_plink_snps() {
        let expected_num_snps = 12;
        let snps = PlinkSnps::new(
            vec![0b10_00_11_00, 0b00_00_11_11, 0b11_10_10_10, 0b11001100],
            expected_num_snps,
        );
        let expected: Vec<u8> = vec![2, 0, 2, 1, 0, 0, 2, 2, 1, 1, 1, 0];
        let mut num_snps = 0;
        for (s, e) in snps.to_iter().zip(expected.iter()) {
            assert_eq!(s, *e);
            num_snps += 1;
        }
        assert_eq!(num_snps, expected_num_snps);
    }

    #[test]
    fn test_plink_snps_from_geno() {
        fn test(geno: Vec<u8>) {
            let plink_snps = PlinkSnps::from_geno(geno.clone());
            let num_snps = geno.len();
            let mut iter = plink_snps.into_iter();
            for i in 0..num_snps {
                assert_eq!(Some(geno[i]), iter.next());
            }
            assert_eq!(None, iter.next());
        }
        test(vec![]);
        test(vec![0]);
        test(vec![1]);
        test(vec![2]);
        test(vec![0, 1, 1]);
        test(vec![2, 2, 0, 1]);
        test(vec![1, 2, 0, 1, 1, 0]);
        test(vec![1, 2, 0, 1, 1, 0, 2, 1]);
        test(vec![1, 2, 2, 0, 0, 1, 0, 0, 2, 0]);
    }

    #[test]
    fn test_create_dominance_geno_bed() {
        fn test(geno: &Array<u8, Ix2>) {
            let mut bim = NamedTempFile::new().unwrap();
            let mut fam = NamedTempFile::new().unwrap();
            create_dummy_bim_fam(&mut bim, &mut fam, geno.dim().0, geno.dim().1).unwrap();
            let bed_path = NamedTempFile::new().unwrap().into_temp_path().to_str().unwrap().to_string();
            let bim_path = bim.into_temp_path();
            let fam_path = fam.into_temp_path();
            PlinkBed::create_bed(&geno, &bed_path).unwrap();
            let geno_bed = PlinkBed::new(
                &bed_path,
                bim_path.to_str().unwrap(),
                fam_path.to_str().unwrap(),
            ).unwrap();
            let dominance_path = NamedTempFile::new().unwrap().into_temp_path();
            geno_bed.create_dominance_geno_bed(dominance_path.to_str().unwrap()).unwrap();
            let dominance_geno = PlinkBed::new(
                dominance_path.to_str().unwrap(),
                bim_path.to_str().unwrap(),
                fam_path.to_str().unwrap(),
            ).unwrap();
            assert_eq!(
                geno_bed.get_genotype_matrix(None)
                        .unwrap()
                        .mapv(|s| match s as u8 {
                            2 => 1u8,
                            s => s,
                        }),
                dominance_geno.get_genotype_matrix(None)
                              .unwrap()
                              .mapv(|s| s as u8)
            );
        }
        test(
            &array![
            [0, 0, 1, 2],
            [1, 1, 2, 1],
            [2, 0, 0, 0],
            [1, 0, 0, 2],
            [0, 2, 1, 0],
        ]);

        test(
            &array![
            [0, 0, 1, 2, 2],
            [1, 1, 2, 1, 0],
            [2, 0, 0, 0, 2],
            [1, 0, 0, 2, 1],
            [0, 2, 1, 0, 1],
        ]);

        test(
            &array![
            [0, 0, 1, 2, 2],
            [1, 1, 2, 1, 0],
            [2, 0, 0, 0, 2],
            [1, 0, 0, 2, 1],
            [0, 1, 2, 1, 2],
            [2, 1, 2, 0, 1],
            [1, 0, 1, 1, 0],
            [2, 1, 0, 2, 0],
        ]);

        test(
            &array![
            [0, 0, 1, 2, 2, 1, 1, 0],
            [1, 1, 2, 1, 0, 0, 0, 0],
            [2, 0, 0, 0, 2, 1, 0, 1],
            [1, 0, 0, 2, 1, 1, 2, 0],
            [0, 1, 2, 1, 2, 1, 1, 2],
            [2, 1, 2, 0, 1, 0, 2, 0],
            [1, 0, 1, 1, 0, 0, 0, 2],
            [2, 1, 0, 2, 0, 0, 1, 1],
        ]);

        test(
            &array![
            [0, 0, 1, 2, 2, 1, 1, 0, 2],
            [1, 1, 2, 1, 0, 0, 0, 0, 1],
            [2, 0, 0, 0, 2, 1, 0, 1, 1],
            [1, 0, 0, 2, 1, 1, 2, 0, 2],
            [0, 1, 2, 1, 2, 1, 1, 2, 2],
            [2, 1, 2, 0, 1, 0, 2, 0, 0],
            [1, 0, 1, 1, 0, 0, 0, 2, 0],
            [2, 1, 0, 2, 0, 0, 1, 1, 2],
        ]);
    }
}
