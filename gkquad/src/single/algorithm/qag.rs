#[cfg(not(feature = "std"))]
use crate::float::Float;

use crate::error::IntegrationResult;
use crate::single::algorithm::{Algorithm, QAG_FINITE};
use crate::single::common::{ITransform, Integrand, IntegrationConfig, Interval};
use crate::single::util::transform_interval;

#[derive(Clone)]
pub struct QAG {
    #[doc(hidden)]
    inner: QAG_FINITE
}

impl QAG {
    #[inline]
    pub fn new() -> Self {
        Self {
            inner: QAG_FINITE::new()
        }
    }
}

impl<F: Integrand> Algorithm<F> for QAG {
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
