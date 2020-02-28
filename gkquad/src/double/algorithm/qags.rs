use super::super::common::{Integrand2, IntegrationConfig2};
use super::super::range::{DynamicX, DynamicY, Rectangle};
use super::Algorithm2;
use crate::single::algorithm::{Algorithm as Algorithm1, QAGS};
use crate::single::IntegrationConfig;
use crate::single::WorkSpace;
use crate::IntegrationResult;

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
        square: &Rectangle,
        config: &IntegrationConfig2,
    ) -> IntegrationResult {
        self.integrate(f, &DynamicY::from(square.clone()), config)
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
        let config1 = IntegrationConfig {
            tolerance: config.tolerance.clone(),
            max_iters: config.max_iters,
            ..Default::default()
        };

        let mut inner_ws = WorkSpace::new();
        let mut inner = QAGS::with_workspace(&mut inner_ws);
        let mut error = None;

        let mut integrand = |x: f64| -> f64 {
            let mut integrand2 = |y: f64| f.apply((x, y));
            let result = inner.integrate(&mut integrand2, &(range.yrange)(x), &config1);
            if result.has_err() {
                error = result.err();
            }
            result.estimate().unwrap_or(core::f64::NAN)
        };

        let mut result = QAGS::new().integrate(&mut integrand, &range.xrange, &config1);
        if error.is_some() {
            result.error = error;
        }

        result
    }
}

extra_traits!(QAGS2);
