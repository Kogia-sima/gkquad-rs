use super::algorithm::*;
use super::common::{Integrand2, Integration2Config, Points2, Range2};

use crate::{IntegrationResult, Tolerance};

/// 2-dimentional integration Executor
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Integrator2<F: Integrand2, A: Algorithm<F>> {
    integrand: F,
    algorithm: A,
    config: Integration2Config,
}

impl<F: Integrand2, A: Algorithm<F>> Integrator2<F, A> {
    #[inline]
    pub fn new(integrand: F, algorithm: A) -> Integrator2<F, A> {
        Self {
            integrand,
            algorithm,
            config: Integration2Config::default(),
        }
    }

    /// Set tolerance
    #[inline]
    pub fn tolerance(mut self, t: Tolerance) -> Self {
        assert!(!t.contains_nan(), "Tolerance must not contain NAN.");

        self.config.tolerance = t;
        self
    }

    /// Set maximum number of subranges
    #[inline]
    pub fn limit(mut self, limit: usize) -> Self {
        self.config.limit = limit;
        self
    }

    /// Set singular points
    pub fn points(mut self, pts: &[(f64, f64)]) -> Self {
        assert!(
            pts.iter().all(|(x, y)| !x.is_nan() && !y.is_nan()),
            "cannot include NAN value in singular points: {:?}",
            pts
        );

        self.config.points = Points2::from(pts);

        self
    }

    #[inline]
    pub fn get_algorithm(&self) -> &A {
        &self.algorithm
    }

    #[inline]
    pub fn run<T: Into<Range2>>(&mut self, range: T) -> IntegrationResult {
        self.algorithm
            .integrate(&mut self.integrand, &range.into(), &self.config)
    }
}
