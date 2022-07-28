use std::{cmp::Ordering, io};

use angsd_io::glf;

use rand::{rngs::StdRng, SeedableRng};

use crate::{Cli, Reader, Record, Simulator};

pub type DynRunner = Runner<Box<(dyn io::BufRead)>, Box<(dyn io::Write)>>;

pub struct Runner<R, W>
where
    R: io::BufRead,
    W: io::Write,
{
    reader: Reader<R>,
    writer: glf::BgzfWriter<W>,
    in_buf: Record,
    out_buf: Vec<glf::Record>,
    simulator: Simulator,
    interleaver: Option<Interleaver>,
    rng: StdRng,
}

impl<R, W> Runner<R, W>
where
    R: io::BufRead,
    W: io::Write,
{
    pub fn run(&mut self) -> io::Result<()> {
        // Use Option::take here only to appease the borrow checker
        if let Some(mut interleaver) = self.interleaver.take() {
            self.run_with_interleaver(&mut interleaver)?;

            self.interleaver = Some(interleaver);
        } else {
            self.run_without_interleaver()?;
        }

        Ok(())
    }

    pub fn setup(reader: R, writer: W, args: &Cli) -> io::Result<Self> {
        let mut reader = Reader::new(reader);

        let writer = glf::BgzfWriter::from_bgzf(writer);

        let mut in_buf = Record::default();
        reader.read_record(&mut in_buf)?;

        let n = in_buf.genotypes().len();
        let out_buf = vec![glf::Record::default(); n];

        let simulator = Simulator::setup(n, &args.mean_depths, &args.error_rates);

        let interleaver = args.interleave.then(|| Interleaver::setup(&in_buf));

        let rng = match args.seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };

        Ok(Self {
            reader,
            writer,
            in_buf,
            out_buf,
            simulator,
            interleaver,
            rng,
        })
    }

    fn run_with_interleaver(&mut self, interleaver: &mut Interleaver) -> io::Result<()> {
        loop {
            match interleaver.next_record_type() {
                Some(RecordType::Interleaved) => {
                    self.simulator
                        .simulate_monomorphic(&mut self.out_buf, &mut self.rng);

                    self.writer.write_records(self.out_buf.as_slice())?;
                }
                Some(RecordType::External) => {
                    self.simulator.simulate(
                        &mut self.out_buf,
                        self.in_buf.genotypes(),
                        &mut self.rng,
                    );

                    self.writer.write_records(self.out_buf.as_slice())?;

                    if self.reader.read_record(&mut self.in_buf)?.is_done() {
                        break;
                    } else {
                        interleaver.update(&self.in_buf);
                    }
                }
                None => return Err(io::Error::new(io::ErrorKind::InvalidData, "unsorted input")),
            }
        }

        Ok(())
    }

    fn run_without_interleaver(&mut self) -> io::Result<()> {
        loop {
            self.simulator
                .simulate(&mut self.out_buf, self.in_buf.genotypes(), &mut self.rng);

            self.writer.write_records(self.out_buf.as_slice())?;

            if self.reader.read_record(&mut self.in_buf)?.is_done() {
                break;
            }
        }

        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum RecordType {
    Interleaved,
    External,
}

#[derive(Clone, Debug)]
struct Interleaver {
    contig: String,
    current_position: u64,
    next_position: u64,
}

impl Interleaver {
    pub fn next_record_type(&mut self) -> Option<RecordType> {
        self.current_position += 1;

        match self.current_position.cmp(&self.next_position) {
            Ordering::Less => Some(RecordType::Interleaved),
            Ordering::Equal => Some(RecordType::External),
            Ordering::Greater => None,
        }
    }

    pub fn setup(record: &Record) -> Self {
        Self {
            contig: record.contig().to_string(),
            current_position: 0,
            next_position: record.position(),
        }
    }

    pub fn update(&mut self, record: &Record) {
        if record.contig() != self.contig {
            self.contig = record.contig().to_string();
            self.current_position = 0;
        }

        self.next_position = record.position();
    }
}
