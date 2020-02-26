use crate::error::IntegrationResult;
use crate::single::algorithm::*;
use crate::single::common::{Integrand, IntegrationConfig, Range};

/// Automatically select algorithm based on configuration
#[derive(Clone)]
#[deprecated(since = "0.0.3", note = "Use QAGS or QAGP method instead.")]
pub struct AUTO;

impl AUTO {
    #[inline(always)]
    pub const fn new() -> Self {
        Self
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
            let mut qags = QAGS::new();
            qags.integrate(f, range, config)
        } else {
            let mut qagp = QAGP::new();
            qagp.integrate(f, range, config)
        }
    }
}

extra_traits!(AUTO);
