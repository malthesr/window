use std::{
    error, fmt,
    num::ParseIntError,
    ops::{Deref, DerefMut},
    str::FromStr,
};

const RECORD_SEPARATOR: char = '\t';
const GENOTYPE_SEPARATOR: char = ':';
const ALLELE_SEPARATOR: &[char] = &['|', '/'];

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Record {
    contig: String,
    position: u64,
    genotypes: Genotypes,
}

impl Record {
    pub fn new(contig: String, position: u64, genotypes: Genotypes) -> Self {
        Self {
            contig,
            position,
            genotypes,
        }
    }

    pub fn contig(&self) -> &str {
        &self.contig
    }

    pub fn genotypes(&self) -> &Genotypes {
        &self.genotypes
    }

    pub fn parse_into(&mut self, s: &str) -> Result<(), ParseRecordError> {
        let create_error = |error_type| ParseRecordError::new(s.to_string(), error_type);

        let (raw_contig, (raw_position, raw_genotypes)) = s
            .split_once(RECORD_SEPARATOR)
            .and_then(|(l, mr)| mr.split_once(RECORD_SEPARATOR).map(|mr| (l, mr)))
            .ok_or_else(|| create_error(ParseRecordErrorType::MissingField))?;

        self.contig.clear();
        self.contig.push_str(raw_contig);

        self.position = match raw_position.parse() {
            Ok(v) => v,
            Err(error) => {
                return Err(create_error(ParseRecordErrorType::ParsePositionError(
                    error,
                )))
            }
        };

        self.genotypes
            .parse_into(raw_genotypes)
            .map_err(|error| create_error(ParseRecordErrorType::ParseGenotypeError(error)))?;

        Ok(())
    }

    pub fn position(&self) -> u64 {
        self.position
    }
}

impl FromStr for Record {
    type Err = ParseRecordError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut new = Record::default();

        new.parse_into(s)?;

        Ok(new)
    }
}

#[derive(Debug)]
pub struct ParseRecordError {
    raw_record: String,
    error_type: ParseRecordErrorType,
}

impl ParseRecordError {
    fn new(raw_record: String, error_type: ParseRecordErrorType) -> Self {
        Self {
            raw_record,
            error_type,
        }
    }
}

impl fmt::Display for ParseRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "failed to parse record '{}' with the following reason: {}",
            self.raw_record, self.error_type
        )
    }
}

impl error::Error for ParseRecordError {}

#[derive(Clone, Debug, Eq, PartialEq)]
enum ParseRecordErrorType {
    MissingField,
    ParseGenotypeError(ParseGenotypeError),
    ParsePositionError(ParseIntError),
}

impl fmt::Display for ParseRecordErrorType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingField => f.write_str("missing record fields"),
            Self::ParseGenotypeError(error) => write!(f, "{error}"),
            Self::ParsePositionError(error) => write!(f, "{error}"),
        }
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Genotypes(Vec<Genotype>);

impl Genotypes {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn parse_into(&mut self, s: &str) -> Result<(), ParseGenotypeError> {
        self.0.clear();

        for raw_genotype in s.split(GENOTYPE_SEPARATOR) {
            self.0.push(raw_genotype.parse()?)
        }

        Ok(())
    }
}

impl Deref for Genotypes {
    type Target = [Genotype];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Genotypes {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl From<Vec<Genotype>> for Genotypes {
    fn from(genotypes: Vec<Genotype>) -> Self {
        Self(genotypes)
    }
}

impl FromStr for Genotypes {
    type Err = ParseGenotypeError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut new = Self::new();

        new.parse_into(s)?;

        Ok(new)
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
#[repr(u8)]
pub enum Genotype {
    Zero = 0,
    One = 1,
    Two = 2,
}

impl FromStr for Genotype {
    type Err = ParseGenotypeError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let allele_sum = s.split_once(ALLELE_SEPARATOR).map(|(fst, snd)| {
            fst.parse::<u8>()
                .and_then(|fst| snd.parse::<u8>().map(|snd| fst + snd))
        });

        match allele_sum {
            Some(Ok(0)) => Ok(Self::Zero),
            Some(Ok(1)) => Ok(Self::One),
            Some(Ok(2)) => Ok(Self::Two),
            Some(Ok(_)) | Some(Err(_)) | None => Err(ParseGenotypeError::new(s.to_string())),
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ParseGenotypeError {
    raw_genotype: String,
}

impl ParseGenotypeError {
    fn new(raw_genotype: String) -> Self {
        Self { raw_genotype }
    }
}

impl fmt::Display for ParseGenotypeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "failed to parse genotype '{}'", self.raw_genotype)
    }
}

impl error::Error for ParseGenotypeError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genotype_from_str() {
        assert_eq!(Genotype::from_str("0/0").unwrap(), Genotype::Zero);
        assert_eq!(Genotype::from_str("0|1").unwrap(), Genotype::One);
        assert_eq!(Genotype::from_str("1/1").unwrap(), Genotype::Two);
    }

    #[test]
    fn test_genotypes_from_str() -> Result<(), ParseGenotypeError> {
        assert_eq!(
            "0/1:1|1".parse::<Genotypes>()?,
            Genotypes::from(vec![Genotype::One, Genotype::Two]),
        );
        assert!(matches!(
            "".parse::<Genotypes>(),
            Err(ParseGenotypeError { .. })
        ));

        Ok(())
    }

    #[test]
    fn test_record_from_str() {
        assert_eq!(
            Record::from_str("chr1\t1\t0/0:1|1").unwrap(),
            Record::new(
                String::from("chr1"),
                1,
                Genotypes::from(vec![Genotype::Zero, Genotype::Two])
            )
        );
    }
}
