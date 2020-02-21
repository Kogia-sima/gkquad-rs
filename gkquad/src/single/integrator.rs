use super::algorithm::*;
use super::common::{Integrand, IntegrationConfig, Points, Range};

use crate::{IntegrationResult, Tolerance};

/// Integration Executor
///
/// This object is useful when you want to re-use the configuration.
///
/// ```
/// use std::cell::Cell;
///
/// use gkquad::Tolerance;
/// use gkquad::single::Integrator;
///
/// let m = Cell::new(1);
/// let mut result = 1.0;
///
/// let mut integrator = Integrator::new(|x: f64| x.powi(m.get()))
///     .tolerance(Tolerance::Relative(1e-7))
///     .max_iters(100);
///
/// while m.get() <= 10 {
///     result *= integrator.run(0.0..1.0).estimate().unwrap();
///     
///     // increment the exponent
///     m.set(m.get() + 1);
/// }
/// ```
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Integrator<F: Integrand, A: Algorithm<F>> {
    integrand: F,
    algorithm: A,
    config: IntegrationConfig,
}

impl<F: Integrand> Integrator<F, AUTO> {
    #[inline]
    pub fn new(integrand: F) -> Integrator<F, AUTO> {
        Self {
            integrand,
            algorithm: AUTO::new(),
            config: IntegrationConfig::default(),
        }
    }
}

impl<F: Integrand, A: Algorithm<F>> Integrator<F, A> {
    #[inline]
    pub fn with_algorithm(integrand: F, algorithm: A) -> Integrator<F, A> {
        Self {
            integrand,
            algorithm,
            config: IntegrationConfig::default(),
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
    pub fn max_iters(mut self, max_iters: usize) -> Self {
        self.config.max_iters = max_iters;
        self
    }

    /// Set singular points
    pub fn points(mut self, pts: &[f64]) -> Self {
        assert!(
            pts.iter().all(|&x| !x.is_nan()),
            "cannot include NAN value in singular points: {:?}",
            pts
        );

        self.config.points = Points::from(pts);

        self
    }

    #[inline]
    pub fn get_algorithm(&self) -> &A {
        &self.algorithm
    }

    #[inline]
    pub fn run<T: Into<Range>>(&mut self, range: T) -> IntegrationResult {
        self.algorithm
            .integrate(&mut self.integrand, &range.into(), &self.config)
    }
}

// Integrator can safely implement Eq because Nan value of tolerance is always checked.
impl<F: Integrand + Eq, A: Algorithm<F> + Eq> Eq for Integrator<F, A> {}

impl<F: Integrand> From<F> for Integrator<F, AUTO> {
    #[inline]
    fn from(integrand: F) -> Self {
        Self::new(integrand)
    }
}

impl<F: Integrand + Clone> From<&F> for Integrator<F, AUTO> {
    #[inline]
    fn from(integrand: &F) -> Self {
        Self::new(integrand.clone())
    }
}
