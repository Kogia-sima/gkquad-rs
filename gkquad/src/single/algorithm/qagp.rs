use crate::error::IntegrationResult;
use crate::single::algorithm::{qagp_finite::QAGP_FINITE, Algorithm};
use crate::single::common::{ITransform, Integrand, IntegrationConfig, Points, Range};
use crate::single::util::{transform_point, transform_range};

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
}

impl<F: Integrand> Algorithm<F> for QAGP {
    #[inline]
    fn integrate(&self, f: &mut F, range: &Range, config: &IntegrationConfig) -> IntegrationResult {
        if !range.begin.is_finite() || !range.end.is_finite() {
            let mut f = ITransform(f);
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

            self.inner.integrate(&mut f, &range, &new_config)
        } else {
            self.inner.integrate(f, range, config)
        }
    }
}

extra_traits!(QAGP);
