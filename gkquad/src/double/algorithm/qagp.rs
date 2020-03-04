use alloc::borrow::Cow;

use super::super::common::{Integrand2, IntegrationConfig2};
use super::super::range::{DynamicX, DynamicY, Rectangle};
use super::Algorithm2;
use crate::common::IntegrationResult;
use crate::single::algorithm::{Algorithm, QAGP};
use crate::single::{IntegrationConfig, Points, Range, WorkSpace};

#[cfg(not(feature = "std"))]
use crate::float::Float;

pub struct QAGP2;

impl QAGP2 {
    pub fn new() -> Self {
        Self
    }
}

impl<F: Integrand2> Algorithm2<F, Rectangle> for QAGP2 {
    fn integrate(
        &mut self,
        f: &mut F,
        range: &Rectangle,
        config: &IntegrationConfig2,
    ) -> IntegrationResult {
        let yrange = |_: f64| Cow::Borrowed(&range.yrange);
        integrate_impl(f, &range.xrange, yrange, config)
    }
}

impl<'a, F: Integrand2> Algorithm2<F, DynamicX<'a>> for QAGP2 {
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
        let config = IntegrationConfig2 {
            tolerance: config.tolerance.clone(),
            max_calls: config.max_calls,
            points: config.points.iter().map(|&(x, y)| (y, x)).collect(),
        };
        self.integrate(&mut g, &range, &config)
    }
}

impl<'a, F: Integrand2> Algorithm2<F, DynamicY<'a>> for QAGP2 {
    fn integrate(
        &mut self,
        f: &mut F,
        range: &DynamicY<'a>,
        config: &IntegrationConfig2,
    ) -> IntegrationResult {
        let yrange = |x: f64| Cow::Owned((range.yrange)(x));
        integrate_impl(f, &range.xrange, yrange, config)
    }
}

extra_traits!(QAGP2);

fn integrate_impl<'a, F, G>(
    f: &mut F,
    xrange: &Range,
    yrange: G,
    config: &IntegrationConfig2,
) -> IntegrationResult
where
    F: Integrand2,
    G: Fn(f64) -> Cow<'a, Range>,
{
    let mut inner_config = IntegrationConfig {
        tolerance: config.tolerance.clone(),
        max_calls: config.max_calls,
        points: Points::with_capacity(config.points.len()),
    };

    let mut outer_config = IntegrationConfig {
        tolerance: config.tolerance.clone(),
        max_calls: config.max_calls,
        points: Points::with_capacity(config.points.len()),
    };

    let mut inner_ws = WorkSpace::new();
    let mut inner = QAGP::with_workspace(&mut inner_ws);
    let mut error = None;

    let xtransform = !xrange.begin.is_finite() || !xrange.end.is_finite();
    config.points.iter().for_each(|&(x, y)| {
        if xtransform {
            outer_config.points.push(transform_point(x));
        } else {
            outer_config.points.push(x)
        }
        inner_config.points.push(y);
    });

    let mut integrand = |x: f64| -> f64 {
        let mut integrand2 = |y: f64| f.apply((x, y));
        let result = inner.integrate(&mut integrand2, &*yrange(x), &inner_config);

        unsafe {
            if result.has_err() {
                error = result.err();
                std::f64::NAN
            } else {
                result.unwrap_unchecked().estimate
            }
        }
    };

    let mut result = QAGP::new().integrate(&mut integrand, xrange, &outer_config);
    if error.is_some() {
        result.error = error;
    }

    result
}

#[inline]
fn transform_point(x: f64) -> f64 {
    if x == core::f64::NEG_INFINITY {
        -1.0
    } else if x == core::f64::INFINITY {
        1.0
    } else {
        x / (1.0 + x.abs())
    }
}
