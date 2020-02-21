use alloc::borrow::Cow;

use super::super::common::{Integrand2, IntegrationConfig2, Range2};
use super::Algorithm2;
use crate::single::algorithm::{Algorithm, QAGP};
use crate::single::IntegrationConfig;
use crate::single::WorkSpaceId;
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
        let mut config1 = IntegrationConfig {
            tolerance: config.tolerance.clone(),
            limit: config.limit,
            ..Default::default()
        };

        let config2 = config1.clone();

        let inner = QAGP::with_id(WorkSpaceId::Single);
        let mut error = None;

        let xrange = match range {
            &Range2::Square { ref xrange, .. } => xrange,
            &Range2::Custom { ref xrange, .. } => xrange,
        };

        let yrange = |x: f64| match range {
            &Range2::Square { ref yrange, .. } => Cow::Borrowed(yrange),
            &Range2::Custom { ref yrange, .. } => Cow::Owned(yrange(x)),
        };

        let transform = !xrange.begin.is_finite() || !xrange.end.is_finite();

        // look up table
        let point_lut = if transform {
            Cow::Owned(
                config
                    .points
                    .iter()
                    .map(|(x, y)| (transform_point(*x), *y))
                    .collect(),
            )
        } else {
            Cow::Borrowed(&config.points)
        };

        let mut integrand = |x: f64| -> f64 {
            let ypoints = point_lut
                .iter()
                .filter_map(|(x2, y2)| if *x2 == x { Some(*y2) } else { None });

            if ypoints.clone().count() > 0 {
                config1.points.extend(ypoints);
            }

            let mut integrand2 = |y: f64| f.apply((x, y));
            let result = inner.integrate(&mut integrand2, &yrange(x), &config1);
            if result.has_err() {
                error = result.err();
            }

            if !config1.points.is_empty() {
                config1.points.clear();
            }

            result.estimate().unwrap_or(std::f64::NAN)
        };

        let mut result =
            QAGP::with_id(WorkSpaceId::Double).integrate(&mut integrand, xrange, &config2);
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
