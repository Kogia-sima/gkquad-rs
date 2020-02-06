use crate::error::IntegrationResult;
use crate::single::algorithm::{Algorithm, QAG_FINITE};
use crate::single::common::{Integrand, IntegrationConfig, Interval, ITransform};
use crate::single::util::transform_interval;

#[derive(Clone)]
pub struct QAG;

impl QAG {
    #[inline]
    pub fn new() -> Self {
        Self
    }
}

impl<F: Integrand> Algorithm<F> for QAG {
    fn integrate(
        &self,
        f: &mut F,
        interval: &Interval,
        config: &IntegrationConfig,
    ) -> IntegrationResult {
        let qag_finite = QAG_FINITE::new();

        if !interval.begin.is_finite() || !interval.end.is_finite() {
            let mut f = ITransform(f);
            let interval = transform_interval(interval);
            qag_finite.integrate(&mut f, &interval, config)
        } else {
            qag_finite.integrate(f, interval, config)
        }
    }
}
