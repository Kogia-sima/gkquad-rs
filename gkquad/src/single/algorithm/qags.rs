use crate::error::IntegrationResult;
use crate::single::algorithm::{qags_finite::QAGS_FINITE, Algorithm};
use crate::single::common::{ITransform, Integrand, IntegrationConfig, Range};
use crate::single::util::transform_range;

#[derive(Clone)]
pub struct QAGS {
    inner: QAGS_FINITE,
}

impl QAGS {
    #[inline]
    pub fn new() -> Self {
        Self {
            inner: QAGS_FINITE::new(),
        }
    }
}

impl<F: Integrand> Algorithm<F> for QAGS {
    fn integrate(&self, f: &mut F, range: &Range, config: &IntegrationConfig) -> IntegrationResult {
        if !range.begin.is_finite() || !range.end.is_finite() {
            let mut f = ITransform(f);
            let range = transform_range(range);
            self.inner.integrate(&mut f, &range, config)
        } else {
            self.inner.integrate(f, range, config)
        }
    }
}

extra_traits!(QAGS);
