use crate::error::IntegrationResult;
use crate::single::algorithm::*;
use crate::single::common::{Integrand, IntegrationConfig, Interval, IntegrationWrapper};

#[derive(Clone)]
pub struct AUTO;

impl AUTO {
    #[inline]
    pub fn new() -> Self {
        Self
    }
}

impl<F: Integrand> Algorithm<F> for AUTO {
    #[inline]
    fn integrate(
        &self,
        f: &mut IntegrationWrapper<F>,
        interval: &Interval,
        config: &IntegrationConfig,
    ) -> IntegrationResult {
        if config.points.is_empty() {
            QAGS::new().integrate(f, interval, config)
        } else {
            QAGP::new().integrate(f, interval, config)
        }
    }
}
