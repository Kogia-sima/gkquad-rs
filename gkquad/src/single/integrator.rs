use smallbox::{SmallBox, smallbox};
use smallbox::space::S4;
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
pub struct Integrator<F: Integrand> {
    integrand: UnsafeCell<F>,
    algorithm: SmallBox<dyn Algorithm<F>, S4>,
    config: IntegrationConfig,
}

impl<F: Integrand> Integrator<F> {
    #[inline]
    pub fn new(f: F) -> Integrator<F> {
        Self {
            integrand: UnsafeCell::new(f),
            algorithm: smallbox!(AUTO::new()),
            config: IntegrationConfig::default(),
        }
    }

    #[inline]
    pub fn with_config(f: F, config: IntegrationConfig) -> Integrator<F> {
        Self {
            integrand: UnsafeCell::new(f),
            algorithm: smallbox!(AUTO::new()),
            config,
        }
    }

    /// Set integration algorithm
    #[inline]
    pub fn algorithm<A: Algorithm<F> + 'static>(mut self, algorithm: A) -> Self {
        self.algorithm = smallbox!(algorithm);
        self
    }

    /// Set integration algorithm
    #[inline]
    pub fn boxed_algorithm(mut self, algorithm: SmallBox<dyn Algorithm<F>, S4>) -> Self {
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
        &*self.algorithm
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

impl<F: Integrand + Debug> Debug for Integrator<F> {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        fmt.debug_struct("Integrator")
            .field("integrand", &self.integrand)
            .field("config", &self.config)
            .finish()
    }
}

impl<F: Integrand + Default> Default for Integrator<F> {
    fn default() -> Self {
        Self {
            integrand: UnsafeCell::new(F::default()),
            algorithm: smallbox!(AUTO::new()),
            config: IntegrationConfig::default()
        }
    }
}

unsafe impl<F: Integrand + Send> Send for Integrator<F> {}
