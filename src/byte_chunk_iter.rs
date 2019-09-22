use std::io;
use std::cmp::min;
use std::io::{BufReader, Read, Seek, SeekFrom};

pub struct ByteChunkIter<R> {
    pub start_byte_index: usize,
    end_byte_index_exclusive: usize,
    current_byte_index: usize,
    pub chunk_size: usize,
    buf: BufReader<R>,
}

impl<R: Seek> ByteChunkIter<R> {
    pub fn new(
        mut buf: BufReader<R>,
        start_byte_index: usize,
        end_byte_index_exclusive: usize,
        chunk_size: usize
    ) -> ByteChunkIter<R> {
        let offset = buf.seek(SeekFrom::Start(start_byte_index as u64)).unwrap() as usize;
        assert_eq!(offset, start_byte_index);
        ByteChunkIter {
            start_byte_index,
            end_byte_index_exclusive,
            current_byte_index: start_byte_index,
            chunk_size,
            buf,
        }
    }
}

impl<R: Seek> Seek for ByteChunkIter<R> {
    fn seek(&mut self, seek_from: SeekFrom) -> io::Result<u64> {
        let offset_from_start = self.buf.seek(seek_from)?;
        self.current_byte_index = offset_from_start as usize;
        Ok(offset_from_start)
    }
}

impl<R: Read> Iterator for ByteChunkIter<R> {
    type Item = Vec<u8>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_byte_index >= self.end_byte_index_exclusive {
            None
        } else {
            let len = min(self.end_byte_index_exclusive - self.current_byte_index, self.chunk_size);
            let mut bytes = vec![0u8; len];
            self.buf.read_exact(bytes.as_mut_slice()).unwrap();
            self.current_byte_index += len;
            Some(bytes)
        }
    }
}
