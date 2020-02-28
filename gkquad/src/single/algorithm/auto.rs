use crate::error::IntegrationResult;
use crate::single::algorithm::*;
use crate::single::common::{Integrand, IntegrationConfig, Range};

/// Automatically select algorithm based on configuration
#[derive(Clone)]
pub struct AUTO {
    qags: QAGS<'static>,
    qagp: QAGP<'static>,
}

impl AUTO {
    #[inline(always)]
    pub fn new() -> Self {
        Self {
            qags: QAGS::new(),
            qagp: QAGP::new(),
        }
    }
}

impl<F: Integrand> Algorithm<F> for AUTO {
    fn integrate(
        &mut self,
        f: &mut F,
        range: &Range,
        config: &IntegrationConfig,
    ) -> IntegrationResult {
        if config.points.is_empty() {
            self.qags.integrate(f, range, config)
        } else {
            self.qagp.integrate(f, range, config)
        }
    }
}

extra_traits!(AUTO);
