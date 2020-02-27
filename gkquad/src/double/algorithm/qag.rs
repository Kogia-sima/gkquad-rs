use super::super::common::{Integrand2, IntegrationConfig2};
use super::super::range::{DynamicY, Rectangle};
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

impl<F: Integrand2> Algorithm2<F, Rectangle> for QAG2 {
    fn integrate(
        &mut self,
        f: &mut F,
        square: &Rectangle,
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

        let mut integrand = |x: f64| -> f64 {
            let mut integrand2 = |y: f64| f.apply((x, y));
            let result = inner.integrate(&mut integrand2, &square.yrange, &config1);
            if result.has_err() {
                error = result.err();
            }
            result.estimate().unwrap_or(core::f64::NAN)
        };

        let mut result = QAG::new().integrate(&mut integrand, &square.xrange, &config1);
        if error.is_some() {
            result.error = error;
        }

        result
    }
}

impl<'a, F: Integrand2> Algorithm2<F, DynamicY<'a>> for QAG2 {
    fn integrate(
        &mut self,
        f: &mut F,
        range: &DynamicY<'a>,
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

        let mut integrand = |x: f64| -> f64 {
            let mut integrand2 = |y: f64| f.apply((x, y));
            let result = inner.integrate(&mut integrand2, &(range.yrange)(x), &config1);
            if result.has_err() {
                error = result.err();
            }
            result.estimate().unwrap_or(core::f64::NAN)
        };

        let mut result = QAG::new().integrate(&mut integrand, &range.xrange, &config1);
        if error.is_some() {
            result.error = error;
        }

        result
    }
}

extra_traits!(QAG2);
