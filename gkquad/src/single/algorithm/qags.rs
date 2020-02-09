use crate::error::IntegrationResult;
use crate::single::algorithm::{Algorithm, qags_finite::QAGS_FINITE};
use crate::single::common::{ITransform, Integrand, IntegrationConfig, Interval};
use crate::single::util::transform_interval;

#[derive(Clone)]
pub struct QAGS {
    inner: QAGS_FINITE
}

impl QAGS {
    #[inline]
    pub fn new() -> Self {
        Self {
            inner: QAGS_FINITE::new()
        }
    }
}

impl<F: Integrand> Algorithm<F> for QAGS {
    fn integrate(
        &self,
        f: &mut F,
        interval: &Interval,
        config: &IntegrationConfig,
    ) -> IntegrationResult {
        if !interval.begin.is_finite() || !interval.end.is_finite() {
            let mut f = ITransform(f);
            let interval = transform_interval(interval);
            self.inner.integrate(&mut f, &interval, config)
        } else {
            self.inner.integrate(f, interval, config)
        }
    }
}

extra_traits!(QAGS);
