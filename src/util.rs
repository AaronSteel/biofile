use std::{
    fs::{File, OpenOptions},
    io::BufReader,
};

use crate::{bed::Bed, bedgraph::BedGraph, error::Error};

pub enum TrackVariant {
    Bed(Bed),
    BedGraph(BedGraph),
}

pub fn get_buf(filename: &str) -> Result<BufReader<File>, Error> {
    match OpenOptions::new().read(true).open(filename) {
        Err(io_error) => Err(Error::IO {
            why: format!("failed to open {}: {}", filename, io_error),
            io_error,
        }),
        Ok(f) => Ok(BufReader::new(f)),
    }
}

#[derive(Eq, PartialEq, Copy, Clone, Debug, Hash)]
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
