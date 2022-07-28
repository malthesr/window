use std::io;

use crate::Record;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Status {
    Done,
    NotDone,
}

impl Status {
    pub fn is_done(&self) -> bool {
        matches!(self, Self::Done)
    }
}

#[derive(Debug)]
pub struct Reader<R> {
    inner: R,
    buf: String,
}

impl<R> Reader<R>
where
    R: io::BufRead,
{
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            buf: String::new(),
        }
    }

    pub fn read_record(&mut self, record: &mut Record) -> io::Result<Status> {
        self.buf.clear();

        if self.inner.read_line(&mut self.buf)? == 0 {
            return Ok(Status::Done);
        }

        if self.buf.ends_with('\n') {
            self.buf.pop();
            if self.buf.ends_with('\r') {
                self.buf.pop();
            }
        }

        record
            .parse_into(&self.buf)
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;

        Ok(Status::NotDone)
    }
}
