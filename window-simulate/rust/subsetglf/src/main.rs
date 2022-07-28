use std::{env::args, io, str::FromStr};

use angsd_io::glf;

fn main() -> io::Result<()> {
    let args: Vec<String> = args().skip(1).collect();
    assert_eq!(args.len(), 3);

    let subset = args[0]
        .split(',')
        .map(usize::from_str)
        .collect::<Result<Vec<_>, _>>()
        .expect("failed to parse subset");
    let n = usize::from_str(&args[1]).expect("failed to parse total number of individuals");

    let glf_path = &args[2];
    let mut reader = glf::BgzfReader::from_bgzf_path(glf_path)?;

    let mut records_buf = vec![None; n];
    for i in subset {
        if i >= n {
            panic!("invalid subset index")
        }
        records_buf[i] = Some(glf::Record::default());
    }

    let stdout = io::stdout();
    let mut writer = glf::BgzfWriter::from_bgzf(stdout.lock());

    while reader.read_some_records(&mut records_buf)?.is_not_done() {
        for record in records_buf.iter().filter_map(|x| x.as_ref()) {
            writer.write_record(record)?;
        }
    }

    Ok(())
}
