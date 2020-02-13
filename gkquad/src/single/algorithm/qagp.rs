use crate::error::IntegrationResult;
use crate::single::algorithm::{qagp_finite::QAGP_FINITE, Algorithm};
use crate::single::common::{Integrand, IntegrationConfig, Points, Range};
use crate::single::util::{transform_point, transform_range, IntegrandWrapper};
use crate::single::workspace::WorkSpaceId;

#[derive(Clone)]
pub struct QAGP {
    inner: QAGP_FINITE,
}

impl QAGP {
    #[inline]
    pub fn new() -> Self {
        Self {
            inner: QAGP_FINITE::new(),
        }
    }

    #[inline]
    pub(crate) fn with_id(id: WorkSpaceId) -> Self {
        Self {
            inner: QAGP_FINITE::with_id(id),
        }
    }
}

impl<F: Integrand> Algorithm<F> for QAGP {
    #[inline]
    fn integrate(&self, f: &mut F, range: &Range, config: &IntegrationConfig) -> IntegrationResult {
        let transform = !range.begin.is_finite() || !range.end.is_finite();
        let mut wrapper = IntegrandWrapper {
            inner: f,
            transform,
        };

        if transform {
            let range = transform_range(range);

            // transform singular points
            let mut points = Points::with_capacity(config.points.len());
            config
                .points
                .iter()
                .zip(points.iter_mut())
                .for_each(|(x, y)| {
                    *y = transform_point(*x);
                });

            let new_config = IntegrationConfig {
                tolerance: config.tolerance.clone(),
                limit: config.limit,
                points,
            };

            self.inner.integrate(&mut wrapper, &range, &new_config)
        } else {
            self.inner.integrate(&mut wrapper, range, config)
        }
    }
}

extra_traits!(QAGP);
