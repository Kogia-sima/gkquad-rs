use super::super::common::{Integrand2, IntegrationConfig2, Range2};
use super::Algorithm2;
use crate::single::algorithm::{Algorithm, QAG};
use crate::single::IntegrationConfig;
use crate::single::WorkSpace;
use crate::IntegrationResult;

pub struct QAG2;

impl QAG2 {
    pub fn new() -> Self {
        Self
    }
}

impl<F: Integrand2> Algorithm2<F> for QAG2 {
    fn integrate(
        &self,
        f: &mut F,
        range: &Range2,
        config: &IntegrationConfig2,
    ) -> IntegrationResult {
        let config1 = IntegrationConfig {
            tolerance: config.tolerance.clone(),
            max_iters: config.max_iters,
            ..Default::default()
        };

        let mut inner_ws = WorkSpace::new();
        let mut inner = QAG::with_workspace(&mut inner_ws);
        let mut error = None;

        match range {
            &Range2::Square {
                ref xrange,
                ref yrange,
            } => {
                let mut integrand = |x: f64| -> f64 {
                    let mut integrand2 = |y: f64| f.apply((x, y));
                    let result = inner.integrate(&mut integrand2, yrange, &config1);
                    if result.has_err() {
                        error = result.err();
                    }
                    result.estimate().unwrap_or(std::f64::NAN)
                };

                let mut result = QAG::new().integrate(&mut integrand, xrange, &config1);
                if error.is_some() {
                    result.error = error;
                }

                result
            }
            &Range2::Custom {
                ref xrange,
                ref yrange,
            } => {
                let mut integrand = |x: f64| -> f64 {
                    let mut integrand2 = |y: f64| f.apply((x, y));
                    let result = inner.integrate(&mut integrand2, &yrange(x), &config1);
                    if result.has_err() {
                        error = result.err();
                    }
                    result.estimate().unwrap_or(std::f64::NAN)
                };

                let mut result = QAG::new().integrate(&mut integrand, xrange, &config1);
                if error.is_some() {
                    result.error = error;
                }

                result
            }
        }
    }
}

extra_traits!(QAG2);
