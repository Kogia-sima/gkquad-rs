use alloc::boxed::Box;

use std::fmt::{self, Debug};

use super::algorithm::*;
use super::common::{Integrand, IntegrationConfig, IntegrationWrapper, Interval, Points};
use super::util::transform_param;

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
    integrand: IntegrationWrapper<F>,
    algorithm: Box<dyn Algorithm<F>>,
    config: IntegrationConfig,
}

impl<F: Integrand> Integrator<F> {
    #[inline]
    pub fn new(f: F) -> Integrator<F> {
        Self {
            integrand: IntegrationWrapper {
                inner: f,
                transform: false,
            },
            algorithm: Box::new(AUTO::new()),
            config: IntegrationConfig::default(),
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
        let mut interval = interval.into();

        // transform interval
        if !interval.begin.is_finite() || !interval.end.is_finite() {
            self.integrand.transform = true;
            interval.begin = transform_param(interval.begin);
            interval.end = transform_param(interval.end);

            // transform singular points
            let old_points = Some(self.config.points.clone());
            self.config
                .points
                .iter_mut()
                .for_each(|x| *x = transform_param(*x));

            let result =
                self.algorithm
                    .integrate(&mut self.integrand, &interval.into(), &self.config);

            // reset state
            self.integrand.transform = false;
            self.config.points = old_points.unwrap();

            result
        } else {
            // compiler hint
            self.integrand.transform = false;

            self.algorithm
                .integrate(&mut self.integrand, &interval.into(), &self.config)
        }
    }
}

impl<F: Integrand + Debug> Debug for Integrator<F> {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        fmt.debug_struct("Integrator")
            .field("integrand", &self.integrand.inner)
            .field("config", &self.config)
            .finish()
    }
}

unsafe impl<F: Integrand + Send> Send for Integrator<F> {}
