use std::{error, fmt, num::ParseFloatError, str::FromStr};

use angsd_io::glf::{self, Genotype::*};

use rand::{
    distributions::{Distribution, Standard},
    Rng,
};

use rand_distr::Poisson;

use crate::record;

const DEFAULT_GENOTYPE_BASES: [[Base; 2]; 3] =
    [[Base::A, Base::A], [Base::A, Base::C], [Base::C, Base::C]];

const A_HOMOZYGOUS: glf::Genotype = AA;
const C_HOMOZYGOUS: glf::Genotype = CC;
const G_HOMOZYGOUS: glf::Genotype = GG;
const T_HOMOZYGOUS: glf::Genotype = TT;

const A_HETEROZYGOUS: [glf::Genotype; 3] = [AC, AG, AT];
const C_HETEROZYGOUS: [glf::Genotype; 3] = [AC, CG, CT];
const G_HETEROZYGOUS: [glf::Genotype; 3] = [AG, CG, GT];
const T_HETEROZYGOUS: [glf::Genotype; 3] = [AT, CT, GT];

const A_REMAINING: [glf::Genotype; 6] = [CC, CG, CT, GG, GT, TT];
const C_REMAINING: [glf::Genotype; 6] = [AA, AG, AT, GG, GT, TT];
const G_REMAINING: [glf::Genotype; 6] = [AA, AC, AT, CC, CT, TT];
const T_REMAINING: [glf::Genotype; 6] = [AA, AC, AG, CC, CG, GG];

#[derive(Debug)]
pub struct Simulator {
    simulators: Vec<IndividualSimulator>,
}

impl Simulator {
    pub fn setup(n: usize, mean_depths: &[MeanDepth], error_rates: &[ErrorRate]) -> Self {
        let simulators = mean_depths
            .iter()
            .cycle()
            .zip(error_rates.iter().cycle())
            .take(n)
            .map(|(mean_depth, error_rate)| IndividualSimulator::setup(*mean_depth, *error_rate))
            .collect();

        Self { simulators }
    }

    pub fn simulate<R>(
        &self,
        records: &mut [glf::Record],
        genotypes: &[record::Genotype],
        rng: &mut R,
    ) where
        R: Rng,
    {
        assert_eq!(records.len(), genotypes.len());

        records
            .iter_mut()
            .zip(genotypes.iter())
            .zip(self.simulators.iter().cycle())
            .for_each(|((record, genotype), simulator)| simulator.simulate(record, *genotype, rng))
    }

    pub fn simulate_monomorphic<R>(&self, records: &mut [glf::Record], rng: &mut R)
    where
        R: Rng,
    {
        records
            .iter_mut()
            .zip(self.simulators.iter().cycle())
            .for_each(|(record, simulator)| simulator.simulate(record, record::Genotype::Zero, rng))
    }
}

#[derive(Debug)]
struct IndividualSimulator {
    depth_distribution: DepthDistribution,
    log_genotype_probabilities: LogGenotypeProbabilities,
    error_rate: ErrorRate,
}

impl IndividualSimulator {
    pub fn setup(mean_depth: MeanDepth, error_rate: ErrorRate) -> Self {
        Self {
            depth_distribution: DepthDistribution::from_mean_depth(mean_depth),
            log_genotype_probabilities: LogGenotypeProbabilities::from_error_rate(error_rate),
            error_rate,
        }
    }

    pub fn simulate<R>(&self, record: &mut glf::Record, genotype: record::Genotype, rng: &mut R)
    where
        R: Rng,
    {
        let depth = self.depth_distribution.sample(rng);

        let pileup = Pileup::simulate(genotype, depth, self.error_rate, rng);

        pileup.calculate_likelihood(record, self.log_genotype_probabilities);

        let max = record
            .to_array()
            .into_iter()
            .fold(f64::NEG_INFINITY, f64::max);
        record.as_mut_slice().iter_mut().for_each(|x| *x -= max);
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
#[repr(u8)]
enum Base {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
}

impl Base {
    pub fn apply_error<R>(&self, error_rate: ErrorRate, rng: &mut R) -> Self
    where
        R: Rng + ?Sized,
    {
        let mut base = *self;

        if rng.gen::<f64>() < f64::from(error_rate) {
            while base == *self {
                base = rng.gen::<Self>();
            }
        }

        base
    }

    #[inline]
    pub fn homozygous_index(&self) -> glf::Genotype {
        match self {
            Self::A => A_HOMOZYGOUS,
            Self::C => C_HOMOZYGOUS,
            Self::G => G_HOMOZYGOUS,
            Self::T => T_HOMOZYGOUS,
        }
    }

    #[inline]
    pub fn heterozygous_indexes(&self) -> [glf::Genotype; 3] {
        match self {
            Self::A => A_HETEROZYGOUS,
            Self::C => C_HETEROZYGOUS,
            Self::G => G_HETEROZYGOUS,
            Self::T => T_HETEROZYGOUS,
        }
    }

    #[inline]
    pub fn remaining_indexes(&self) -> [glf::Genotype; 6] {
        match self {
            Self::A => A_REMAINING,
            Self::C => C_REMAINING,
            Self::G => G_REMAINING,
            Self::T => T_REMAINING,
        }
    }
}

impl Distribution<Base> for Standard {
    fn sample<R>(&self, rng: &mut R) -> Base
    where
        R: Rng + ?Sized,
    {
        match rng.gen_range(0u8..4) {
            0 => Base::A,
            1 => Base::C,
            2 => Base::G,
            _ => Base::T,
        }
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
struct Pileup(Vec<Base>);

impl Pileup {
    pub fn calculate_likelihood(
        &self,
        record: &mut glf::Record,
        log_genotype_probabilities: LogGenotypeProbabilities,
    ) {
        record.as_mut().fill(0.0);

        for base in self.0.iter() {
            record[base.homozygous_index()] += log_genotype_probabilities.homozygous();

            for idx in base.heterozygous_indexes() {
                record[idx] += log_genotype_probabilities.heterozygous();
            }

            for idx in base.remaining_indexes() {
                record[idx] += log_genotype_probabilities.remaining();
            }
        }
    }

    pub fn simulate<R>(
        genotype: record::Genotype,
        depth: u64,
        error_rate: ErrorRate,
        rng: &mut R,
    ) -> Self
    where
        R: Rng,
    {
        let bases = (0..depth)
            .map(|_| {
                let genotype_bases = DEFAULT_GENOTYPE_BASES[genotype as usize];

                genotype_bases[rng.gen_range::<usize, _>(0..2)].apply_error(error_rate, rng)
            })
            .collect();

        Self(bases)
    }
}

#[derive(Debug)]
struct DepthDistribution(Poisson<f64>);

impl DepthDistribution {
    pub fn from_mean_depth(mean_depth: MeanDepth) -> Self {
        Self(Poisson::new(f64::from(mean_depth)).unwrap())
    }
}

impl Distribution<u64> for DepthDistribution {
    fn sample<R>(&self, rng: &mut R) -> u64
    where
        R: Rng + ?Sized,
    {
        // See rust-random.github.io/book/guide-dist.html#integers for notes on casting
        self.0.sample(rng) as u64
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct LogGenotypeProbabilities([f64; 3]);

impl LogGenotypeProbabilities {
    pub fn from_error_rate(error_rate: ErrorRate) -> Self {
        let e = f64::from(error_rate);

        Self([
            (1.0 - e).ln(),
            ((1.0 - e) / 2.0 + (e / 6.0)).ln(),
            (e / 3.0).ln(),
        ])
    }

    #[inline]
    pub fn homozygous(&self) -> f64 {
        self.0[0]
    }

    #[inline]
    pub fn heterozygous(&self) -> f64 {
        self.0[1]
    }

    #[inline]
    pub fn remaining(&self) -> f64 {
        self.0[2]
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, PartialOrd)]
pub struct MeanDepth(f64);

impl From<MeanDepth> for f64 {
    fn from(error_rate: MeanDepth) -> Self {
        error_rate.0
    }
}

impl TryFrom<f64> for MeanDepth {
    type Error = MeanDepthError;

    fn try_from(value: f64) -> Result<Self, Self::Error> {
        if value.is_sign_positive() {
            Ok(Self(value))
        } else {
            Err(MeanDepthError::Negative(value))
        }
    }
}

impl FromStr for MeanDepth {
    type Err = MeanDepthError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let v = s.parse::<f64>()?;

        Self::try_from(v)
    }
}

#[derive(Debug)]
pub enum MeanDepthError {
    ParseFailed(ParseFloatError),
    Negative(f64),
}

impl fmt::Display for MeanDepthError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ParseFailed(error) => write!(f, "{error}"),
            Self::Negative(v) => write!(f, "expected a positive value, found {}", v),
        }
    }
}

impl error::Error for MeanDepthError {}

impl From<ParseFloatError> for MeanDepthError {
    fn from(e: ParseFloatError) -> Self {
        Self::ParseFailed(e)
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, PartialOrd)]
pub struct ErrorRate(f64);

impl From<ErrorRate> for f64 {
    fn from(error_rate: ErrorRate) -> Self {
        error_rate.0
    }
}

impl TryFrom<f64> for ErrorRate {
    type Error = ErrorRateError;

    fn try_from(value: f64) -> Result<Self, Self::Error> {
        if (0.0..=1.0).contains(&value) {
            Ok(Self(value))
        } else {
            Err(ErrorRateError::NotProbability(value))
        }
    }
}

impl FromStr for ErrorRate {
    type Err = ErrorRateError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let v = s.parse::<f64>()?;

        Self::try_from(v)
    }
}

#[derive(Debug)]
pub enum ErrorRateError {
    ParseFailed(ParseFloatError),
    NotProbability(f64),
}

impl fmt::Display for ErrorRateError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ParseFailed(error) => write!(f, "{error}"),
            Self::NotProbability(v) => write!(f, "expected a probability, found {}", v),
        }
    }
}

impl error::Error for ErrorRateError {}

impl From<ParseFloatError> for ErrorRateError {
    fn from(e: ParseFloatError) -> Self {
        Self::ParseFailed(e)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_setup() {
        let mean_depths = &[
            MeanDepth::try_from(1.0).unwrap(),
            MeanDepth::try_from(2.0).unwrap(),
        ];
        let error_rates = &[
            ErrorRate::try_from(0.1).unwrap(),
            ErrorRate::try_from(0.2).unwrap(),
            ErrorRate::try_from(0.3).unwrap(),
        ];
        let simulator = Simulator::setup(10, mean_depths, error_rates);
        assert_eq!(simulator.simulators.len(), 10);
    }
}
