use crate::error::IntegrationResult;
use crate::single::algorithm::{Algorithm, qagp_finite::QAGP_FINITE};
use crate::single::common::{ITransform, Integrand, IntegrationConfig, Interval, Points};
use crate::single::util::{transform_interval, transform_point};

#[derive(Clone)]
pub struct QAGP {
    inner: QAGP_FINITE
}

impl QAGP {
    #[inline]
    pub fn new() -> Self {
        Self {
            inner: QAGP_FINITE::new()
        }
    }
}

impl<F: Integrand> Algorithm<F> for QAGP {
    #[inline]
    fn integrate(
        &self,
        f: &mut F,
        interval: &Interval,
        config: &IntegrationConfig,
    ) -> IntegrationResult {
        if !interval.begin.is_finite() || !interval.end.is_finite() {
            let mut f = ITransform(f);
            let interval = transform_interval(interval);

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

            self.inner.integrate(&mut f, &interval, &new_config)
        } else {
            self.inner.integrate(f, interval, config)
        }
    }
}

extra_traits!(QAGP);
