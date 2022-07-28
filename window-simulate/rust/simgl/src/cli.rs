use std::{
    fs, io,
    path::{Path, PathBuf},
};

use clap::{CommandFactory, Parser};

use crate::{ErrorRate, MeanDepth};

const AUTHOR: &str = env!("CARGO_PKG_AUTHORS");
const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Create ANGSD GLF file from simulate genotyped.
///
/// Program takes a simple TSV format containing simulated diploid, diallelic genotypes, and writes
/// a GLF based on simulating NGS data.
#[derive(Debug, Parser)]
#[clap(author = AUTHOR, version = VERSION, about)]
pub struct Cli {
    #[clap(long, hide = true)]
    pub debug: bool,

    /// Path to input file. Defaults to reading from stdin if no argument is provided.
    ///
    /// Format is plain text lines, with each line consisting of '{chrom}\t{pos}\t{genotypes}', and
    /// genotypes in turn are ':'-delimited VCF style genotypes. For example, a valid input line
    /// would be 'chr1\t1\t0/1:1|1:0/0'.
    pub input: Option<PathBuf>,

    /// Sequencing error rates used for simulating sequencing data.
    ///
    /// Multiple error rates can be provided. If the number of error rates do not match the number of
    ///  samples, values will be cycled for each site. For example, with three samples, and
    /// '-e 0.1 0.01 0.1', samples one and three will have error rate 0.1, while the second
    /// individual will have error rate 0.01. Each value must be in the closed unit interval.
    #[clap(short = 'e', long, multiple_values = true, default_value = "0.01")]
    pub error_rates: Vec<ErrorRate>,

    /// Interleave monomorphic sites.
    ///
    /// By default, sites not present in the input will not be included in the output. By enabling
    /// this flag, these sites will be assumed to be monomorphic and interleaved between the input
    /// records.
    #[clap(short = 'i', long)]
    pub interleave: bool,

    /// Mean depths used for simulating sequencing data.
    ///
    /// Depths are assumed to be Poisson distributed. Mean depth values will cycle in the same
    /// manner as error rates, see '--error-rates' for details. Each value must be non-negative.
    #[clap(short = 'd', long, multiple_values = true, default_value = "2")]
    pub mean_depths: Vec<MeanDepth>,

    /// Path to output GLF file. Defaults to writing to stdout if no argument is provided.
    #[clap(short = 'o', long)]
    pub output: Option<PathBuf>,

    /// Seed used for simulating sequencing data.
    ///
    /// Random seed used by default, set for reproducibility.
    #[clap(short = 's', long)]
    pub seed: Option<u64>,
}

pub fn input_file_or_stdin<P>(input: Option<P>) -> clap::Result<Box<dyn io::BufRead + 'static>>
where
    P: AsRef<Path>,
{
    let stdin_atty = atty::is(atty::Stream::Stdin);

    match (input, stdin_atty) {
        (Some(path), true) => {
            let reader = fs::File::open(path).map(io::BufReader::new)?;

            Ok(Box::new(reader))
        }
        (None, false) => {
            let reader = Box::leak(Box::new(io::stdin())).lock();

            Ok(Box::new(reader))
        }
        (Some(_path), false) => Err(Cli::command().error(
            clap::ErrorKind::ArgumentConflict,
            "Provided input from stdin and path",
        )),
        (None, true) => {
            Err(Cli::command().error(clap::ErrorKind::MissingRequiredArgument, "Missing input"))
        }
    }
}

pub fn output_file_or_stdout<P>(output: Option<P>) -> clap::Result<Box<dyn io::Write + 'static>>
where
    P: AsRef<Path>,
{
    match output {
        Some(path) => {
            let writer = fs::File::create(path).map(io::BufWriter::new)?;

            Ok(Box::new(writer))
        }
        None => {
            let writer = Box::leak(Box::new(io::stdout())).lock();

            Ok(Box::new(writer))
        }
    }
}

// Don't panic on closed output pipe errors on UNIX (courtesy of BurntSushi, see
// stackoverflow.com/a/65760807/11890789)
#[cfg(unix)]
pub fn reset_sigpipe() {
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }
}

// No-op on non-UNIX platforms
#[cfg(not(unix))]
pub fn reset_sigpipe() {}

#[cfg(test)]
mod test {
    use super::*;

    fn parse_args(cmd: &str) -> Cli {
        Parser::parse_from(cmd.split_whitespace())
    }

    #[test]
    fn test_input() {
        let args = parse_args("gladd /path/to/input");
        assert_eq!(args.input, Some(PathBuf::from("/path/to/input")));
    }

    #[test]
    fn test_error_rates() {
        let args = parse_args("gladd --error-rates 0.1");
        assert_eq!(args.error_rates, vec![ErrorRate::try_from(0.1).unwrap()]);

        let args = parse_args("gladd -e 0.1 0.01 0.2");
        assert_eq!(
            args.error_rates,
            vec![
                ErrorRate::try_from(0.1).unwrap(),
                ErrorRate::try_from(0.01).unwrap(),
                ErrorRate::try_from(0.2).unwrap()
            ]
        );

        let args = parse_args("gladd --error-rates 0 1");
        assert_eq!(
            args.error_rates,
            vec![
                ErrorRate::try_from(0.0).unwrap(),
                ErrorRate::try_from(1.0).unwrap(),
            ]
        );
    }

    #[test]
    fn test_interleave() {
        let args = parse_args("gladd");
        assert_eq!(args.interleave, false);

        let args = parse_args("gladd --interleave");
        assert_eq!(args.interleave, true);

        let args = parse_args("gladd -i");
        assert_eq!(args.interleave, true);
    }

    #[test]
    fn test_mean_depths() {
        let args = parse_args("gladd");
        assert_eq!(args.mean_depths, vec![MeanDepth::try_from(2.0).unwrap()]);

        let args = parse_args("gladd --mean-depths 1");
        assert_eq!(args.mean_depths, vec![MeanDepth::try_from(1.0).unwrap()]);

        let args = parse_args("gladd -d 2 5.5 22.0");
        assert_eq!(
            args.mean_depths,
            vec![
                MeanDepth::try_from(2.0).unwrap(),
                MeanDepth::try_from(5.5).unwrap(),
                MeanDepth::try_from(22.0).unwrap()
            ]
        );
    }

    #[test]
    fn test_seed() {
        let args = parse_args("gladd");
        assert_eq!(args.seed, None);

        let args = parse_args("gladd -s 1");
        assert_eq!(args.seed, Some(1));

        let args = parse_args("gladd --seed 112233");
        assert_eq!(args.seed, Some(112233));
    }
}
