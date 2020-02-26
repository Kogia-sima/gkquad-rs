use alloc::borrow::Cow;

use super::super::common::{Integrand2, IntegrationConfig2, Range2};
use super::Algorithm2;
use crate::single::algorithm::{Algorithm, QAGP};
use crate::single::{IntegrationConfig, Points, WorkSpace};
use crate::IntegrationResult;

pub struct QAGP2;

impl QAGP2 {
    pub fn new() -> Self {
        Self
    }
}

impl<F: Integrand2> Algorithm2<F> for QAGP2 {
    fn integrate(
        &self,
        f: &mut F,
        range: &Range2,
        config: &IntegrationConfig2,
    ) -> IntegrationResult {
        let mut inner_config = IntegrationConfig {
            tolerance: config.tolerance.clone(),
            max_iters: config.max_iters,
            points: Points::with_capacity(config.points.len()),
        };

        let mut outer_config = IntegrationConfig {
            tolerance: config.tolerance.clone(),
            max_iters: config.max_iters,
            points: Points::with_capacity(config.points.len()),
        };

        let mut inner_ws = WorkSpace::new();
        let mut inner = QAGP::with_workspace(&mut inner_ws);
        let mut error = None;

        let xrange = match range {
            &Range2::Square { ref xrange, .. } => xrange,
            &Range2::Custom { ref xrange, .. } => xrange,
        };

        let yrange = |x: f64| match range {
            &Range2::Square { ref yrange, .. } => Cow::Borrowed(yrange),
            &Range2::Custom { ref yrange, .. } => Cow::Owned(yrange(x)),
        };

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
            if result.has_err() {
                error = result.err();
            }

            result.estimate().unwrap_or(std::f64::NAN)
        };

        let mut result = QAGP::new().integrate(&mut integrand, xrange, &outer_config);
        if error.is_some() {
            result.error = error;
        }

        result
    }
}

extra_traits!(QAGP2);

#[inline]
fn transform_point(x: f64) -> f64 {
    if x == std::f64::NEG_INFINITY {
        -1.0
    } else if x == std::f64::INFINITY {
        1.0
    } else {
        x / (1.0 + x.abs())
    }
}
