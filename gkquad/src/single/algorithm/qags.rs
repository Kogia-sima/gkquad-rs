use crate::error::IntegrationResult;
use crate::single::algorithm::{qags_finite::QAGS_FINITE, Algorithm};
use crate::single::common::{Integrand, IntegrationConfig, Range};
use crate::single::util::{transform_range, IntegrandWrapper};
use crate::single::workspace::WorkSpaceId;

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

    #[inline]
    pub(crate) fn with_id(id: WorkSpaceId) -> Self {
        Self {
            inner: QAGS_FINITE::with_id(id),
        }
    }
}

impl<F: Integrand> Algorithm<F> for QAGS {
    fn integrate(&self, f: &mut F, range: &Range, config: &IntegrationConfig) -> IntegrationResult {
        let transform = !range.begin.is_finite() || !range.end.is_finite();
        let mut wrapper = IntegrandWrapper {
            inner: f,
            transform,
        };
        if transform {
            let range = transform_range(range);
            self.inner.integrate(&mut wrapper, &range, config)
        } else {
            self.inner.integrate(&mut wrapper, range, config)
        }
    }
}

extra_traits!(QAGS);
