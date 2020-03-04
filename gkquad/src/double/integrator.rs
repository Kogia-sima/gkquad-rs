use super::algorithm::*;
use super::common::{Integrand2, IntegrationConfig2, Points2};
use super::range::IntoRange2;

use crate::common::{IntegrationResult, Tolerance};

/// 2-dimentional integration Executor
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Integrator2<F: Integrand2, A> {
    integrand: F,
    algorithm: A,
    config: IntegrationConfig2,
}

impl<F: Integrand2> Integrator2<F, AUTO2> {
    pub fn new(integrand: F) -> Integrator2<F, AUTO2> {
        Self {
            integrand,
            algorithm: AUTO2::new(),
            config: IntegrationConfig2::default(),
        }
    }
}

impl<F: Integrand2, A> Integrator2<F, A> {
    #[inline]
    pub fn with_algorithm(integrand: F, algorithm: A) -> Integrator2<F, A> {
        Self {
            integrand,
            algorithm,
            config: IntegrationConfig2::default(),
        }
    }

    /// Set algorithm
    #[inline]
    pub fn algorithm<B>(self, algorithm: B) -> Integrator2<F, B> {
        Integrator2 {
            integrand: self.integrand,
            algorithm,
            config: self.config,
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
    pub fn max_evals(mut self, max_evals: usize) -> Self {
        self.config.max_evals = max_evals;
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
    pub fn run<'a, T>(&mut self, range: T) -> IntegrationResult
    where
        T: IntoRange2,
        A: Algorithm2<F, T::IntoRange>,
    {
        self.algorithm
            .integrate(&mut self.integrand, &range.into_range(), &self.config)
    }
}
