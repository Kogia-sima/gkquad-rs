use alloc::boxed::Box;

use std::fmt::{self, Debug};

use super::algorithm::*;
use super::common::{Integrand, IntegrationConfig, Interval, Points};

use crate::{IntegrationResult, Tolerance};

/// Integration Executor
///
/// This object is useful when you want to re-use the configuration.
///
/// ```
/// let m = Cell::new(1);
/// let mut result = 1.0;
///
/// let mut integrator = Integrator::new(|x| x.powi(m.get()))
///     .tolerance(Tolerance::Relative(1e-7))
///     .limit(100);
///
/// while m <= 10 {
///     result *= integrator.run(0.0..1.0).estimate().unwrap();
///     
///     // increment the exponent
///     result.set(result.get() + 1);
/// }
/// ```
pub struct Integrator<F: Integrand> {
    integrand: F,
    algorithm: Box<dyn Algorithm<F>>,
    config: IntegrationConfig,
}

impl<F: Integrand> Integrator<F> {
    #[inline]
    pub fn new(f: F) -> Integrator<F> {
        Self {
            integrand: f,
            algorithm: Box::new(AUTO::new()),
            config: IntegrationConfig::default(),
        }
    }

    #[inline]
    pub(crate) fn with_config(f: F, config: IntegrationConfig) -> Integrator<F> {
        Self {
            integrand: f,
            algorithm: Box::new(AUTO::new()),
            config,
        }
    }

    /// Set integration algorithm
    #[inline]
    pub fn algorithm<A: Algorithm<F> + 'static>(mut self, algorithm: A) -> Self {
        self.algorithm = Box::new(algorithm);
        self
    }

    /// Set integration algorithm
    #[inline]
    pub fn boxed_algorithm(mut self, algorithm: Box<dyn Algorithm<F>>) -> Self {
        self.algorithm = algorithm;
        self
    }

    /// Set tolerance
    #[inline]
    pub fn tolerance(mut self, t: Tolerance) -> Self {
        assert!(!t.contains_nan(), "Tolerance must not contain NAN.");

        self.config.tolerance = t;
        self
    }

    /// Set maximum number of subintervals
    #[inline]
    pub fn limit(mut self, limit: usize) -> Self {
        self.config.limit = limit;
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
    pub fn get_algorithm(&self) -> &dyn Algorithm<F> {
        self.algorithm.as_ref()
    }

    pub fn run<T: Into<Interval>>(&mut self, interval: T) -> IntegrationResult {
        self.algorithm
            .integrate(&mut self.integrand, &interval.into(), &self.config)
    }
}

impl<F: Integrand + Debug> Debug for Integrator<F> {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        fmt.debug_struct("Integrator")
            .field("integrand", &self.integrand)
            .field("config", &self.config)
            .finish()
    }
}

unsafe impl<F: Integrand + Send> Send for Integrator<F> {}
