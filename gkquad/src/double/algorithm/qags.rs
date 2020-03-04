use super::super::common::{Integrand2, IntegrationConfig2};
use super::super::range::{DynamicX, DynamicY, Rectangle};
use super::Algorithm2;
use crate::common::IntegrationResult;
use crate::single::algorithm::{Algorithm as Algorithm1, QAGS};
use crate::single::IntegrationConfig;
use crate::single::WorkSpace;

pub struct QAGS2;

impl QAGS2 {
    pub fn new() -> Self {
        Self
    }
}

impl<F: Integrand2> Algorithm2<F, Rectangle> for QAGS2 {
    fn integrate(
        &mut self,
        f: &mut F,
        range: &Rectangle,
        config: &IntegrationConfig2,
    ) -> IntegrationResult {
        self.integrate(f, &DynamicY::from(range.clone()), config)
    }
}

impl<'a, F: Integrand2> Algorithm2<F, DynamicX<'a>> for QAGS2 {
    fn integrate(
        &mut self,
        f: &mut F,
        range: &DynamicX<'a>,
        config: &IntegrationConfig2,
    ) -> IntegrationResult {
        // swap x and y
        let mut g = |x: f64, y: f64| f.apply((y, x));
        let range = DynamicY {
            xrange: range.yrange.clone(),
            yrange: range.xrange.clone(),
        };
        self.integrate(&mut g, &range, config)
    }
}

impl<'a, F: Integrand2> Algorithm2<F, DynamicY<'a>> for QAGS2 {
    fn integrate(
        &mut self,
        f: &mut F,
        range: &DynamicY<'a>,
        config: &IntegrationConfig2,
    ) -> IntegrationResult {
        let mut config1 = IntegrationConfig {
            tolerance: config.tolerance.clone(),
            max_evals: config.max_evals / 17,
            ..Default::default()
        };
        let config2 = config1.clone();

        let mut inner_ws = WorkSpace::new();
        let mut inner = QAGS::with_workspace(&mut inner_ws);
        let mut error = None;
        let mut nevals = 0usize;

        let mut integrand = |x: f64| -> f64 {
            let mut integrand2 = |y: f64| f.apply((x, y));
            config1.max_evals = if error.is_some() {
                0
            } else {
                config.max_evals - nevals
            };
            let result = inner.integrate(&mut integrand2, &(range.yrange)(x), &config1);

            unsafe {
                if result.has_err() {
                    if error.is_none() {
                        error = result.error;
                    }
                    nevals += result.unwrap_unchecked().nevals;
                    core::f64::NAN
                } else {
                    let result = result.unwrap_unchecked();
                    nevals += result.nevals;
                    result.estimate
                }
            }
        };

        let mut result = QAGS::new().integrate(&mut integrand, &range.xrange, &config2);
        result.value.nevals = nevals;
        if error.is_some() {
            result.error = error;
        }

        result
    }
}

extra_traits!(QAGS2);
