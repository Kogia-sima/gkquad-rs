use crate::error::IntegrationResult;
use crate::single::algorithm::{Algorithm, QAGS_FINITE};
use crate::single::common::{Integrand, IntegrationConfig, Interval, ITransform};
use crate::single::util::transform_interval;

#[derive(Clone)]
pub struct QAGS;

impl QAGS {
    #[inline]
    pub fn new() -> Self {
        Self
    }
}

impl<F: Integrand> Algorithm<F> for QAGS {
    fn integrate(
        &self,
        f: &mut F,
        interval: &Interval,
        config: &IntegrationConfig,
    ) -> IntegrationResult {
        let qags_finite = QAGS_FINITE::new();

        if !interval.begin.is_finite() || !interval.end.is_finite() {
            let mut f = ITransform(f);
            let interval = transform_interval(interval);
            qags_finite.integrate(&mut f, &interval, config)
        } else {
            qags_finite.integrate(f, interval, config)
        }
    }
}
