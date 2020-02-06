use crate::error::IntegrationResult;
use crate::single::algorithm::{Algorithm, QAGP_FINITE};
use crate::single::common::{Integrand, IntegrationConfig, Interval, ITransform, Points};
use crate::single::util::{transform_interval, transform_point};

#[derive(Clone)]
pub struct QAGP;

impl QAGP {
    #[inline]
    pub fn new() -> Self {
        Self
    }
}

impl<F: Integrand> Algorithm<F> for QAGP {
    fn integrate(
        &self,
        f: &mut F,
        interval: &Interval,
        config: &IntegrationConfig,
    ) -> IntegrationResult {
        let qagp_finite = QAGP_FINITE::new();
        if !interval.begin.is_finite() || !interval.end.is_finite() {
            let mut f = ITransform(f);
            let interval = transform_interval(interval);

            // transform singular points
            let mut points = Points::with_capacity(config.points.len());
            config.points.iter().zip(points.iter_mut()).for_each(|(x, y)| {
                *y = transform_point(*x);
            });

            let new_config = IntegrationConfig {
                tolerance: config.tolerance.clone(),
                limit: config.limit,
                points,
            };

            qagp_finite.integrate(&mut f, &interval, &new_config)
        } else {
            qagp_finite.integrate(f, interval, config)
        }
    }
}
