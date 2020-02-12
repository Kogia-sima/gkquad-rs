use std::cell::UnsafeCell;
use std::fmt::{self, Debug};

use super::algorithm::*;
use super::common::{Integrand, IntegrationConfig, Interval, Points};

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
/// let integrator = Integrator::new(|x: f64| x.powi(m.get()))
///     .tolerance(Tolerance::Relative(1e-7))
///     .limit(100);
///
/// while m.get() <= 10 {
///     result *= integrator.run(0.0..1.0).estimate().unwrap();
///     
///     // increment the exponent
///     m.set(m.get() + 1);
/// }
/// ```
pub struct Integrator<F: Integrand, A: Algorithm<F>> {
    integrand: UnsafeCell<F>,
    algorithm: A,
    config: IntegrationConfig,
}

impl<F: Integrand, A: Algorithm<F>> Integrator<F, A> {
    #[inline]
    pub fn new(integrand: F, algorithm: A) -> Integrator<F, A> {
        Self {
            integrand: UnsafeCell::new(integrand),
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
    pub fn get_algorithm(&self) -> &A {
        &self.algorithm
    }

    #[inline]
    pub fn run<T: Into<Interval>>(&self, interval: T) -> IntegrationResult {
        self.algorithm.integrate(
            unsafe { &mut *self.integrand.get() },
            &interval.into(),
            &self.config,
        )
    }
}

impl<F, A> Debug for Integrator<F, A>
where
    F: Integrand + Debug,
    A: Algorithm<F> + Debug,
{
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        fmt.debug_struct("Integrator")
            .field("integrand", &self.integrand)
            .field("algorithm", &self.algorithm)
            .field("config", &self.config)
            .finish()
    }
}

impl<F, A> Default for Integrator<F, A>
where
    F: Integrand + Default,
    A: Algorithm<F> + Default,
{
    fn default() -> Self {
        Self {
            integrand: UnsafeCell::new(F::default()),
            algorithm: A::default(),
            config: IntegrationConfig::default(),
        }
    }
}
