use std::{fmt, io};

pub enum Error {
    IO { why: String, io_error: io::Error },
    BadFormat(String),
    Generic(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Error::BadFormat(why) => write!(f, "Bad format: {}", why),
            Error::IO {
                why, ..
            } => write!(f, "IO error: {}", why),
            Error::Generic(why) => write!(f, "Generic error: {}", why),
        }
    }
}

impl fmt::Debug for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Error::BadFormat(why) => write!(f, "Bad format: {}", why),
            Error::IO {
                why, ..
            } => write!(f, "IO error: {}", why),
            Error::Generic(why) => write!(f, "Generic error: {}", why),
        }
    }
}

impl From<io::Error> for Error {
    fn from(io_error: io::Error) -> Error {
        Error::IO {
            why: io_error.to_string(),
            io_error,
        }
    }
}

impl From<String> for Error {
    fn from(err: String) -> Error {
        Error::Generic(err)
    }
}
