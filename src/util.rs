use std::{
    fs::{File, OpenOptions},
    io::BufReader,
};

use crate::error::Error;

pub fn get_buf(filename: &str) -> Result<BufReader<File>, Error> {
    match OpenOptions::new().read(true).open(filename) {
        Err(io_error) => Err(Error::IO {
            why: format!("failed to open {}: {}", filename, io_error),
            io_error,
        }),
        Ok(f) => Ok(BufReader::new(f)),
    }
}

#[derive(PartialEq, Copy, Clone, Debug)]
pub enum Strand {
    Positive,
    Negative,
}

impl Strand {
    pub fn new(strand: &str) -> Result<Option<Strand>, Error> {
        match strand {
            "+" => Ok(Some(Strand::Positive)),
            "-" => Ok(Some(Strand::Negative)),
            "." => Ok(None),
            s => Err(Error::BadFormat(format!(
                "unrecognized strand symbol: {}",
                s
            ))),
        }
    }
}
