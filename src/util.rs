use std::fs::{File, OpenOptions};
use std::io::BufReader;

use crate::error::Error;

pub fn get_buf(filename: &str) -> Result<BufReader<File>, Error> {
    match OpenOptions::new().read(true).open(filename) {
        Err(io_error) => Err(Error::IO { why: format!("failed to open {}: {}", filename, io_error), io_error }),
        Ok(f) => Ok(BufReader::new(f))
    }
}
