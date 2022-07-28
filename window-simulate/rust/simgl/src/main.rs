use std::process::exit;

mod cli;
use cli::Cli;

pub mod reader;
pub use reader::Reader;

pub mod record;
pub use record::Record;

mod runner;
use runner::{DynRunner, Runner};

mod simulate;
pub use simulate::{ErrorRate, MeanDepth, Simulator};

fn try_setup(args: &Cli) -> clap::Result<DynRunner> {
    let input = cli::input_file_or_stdin(args.input.as_deref())?;
    let output = cli::output_file_or_stdout(args.output.as_deref())?;

    Ok(Runner::setup(input, output, args)?)
}

fn main() {
    let args: Cli = clap::Parser::parse();

    if args.debug {
        eprintln!("{args:#?}");
        exit(0);
    }

    cli::reset_sigpipe();

    let mut runner = try_setup(&args).unwrap_or_else(|error| {
        eprintln!("{error}");
        exit(64);
    });

    match runner.run() {
        Ok(()) => (),
        Err(error) => {
            eprintln!("{error}");
            exit(1);
        }
    }
}
