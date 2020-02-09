use std::cell::RefCell;

use crate::error::IntegrationResult;
use crate::single::algorithm::*;
use crate::single::common::{Integrand, IntegrationConfig, Interval};

/// Automatically select algorithm based on configuration
#[derive(Clone)]
pub struct AUTO {
    // initialize lazily
    qags: RefCell<Option<QAGS>>,
    qagp: RefCell<Option<QAGP>>,
}

impl AUTO {
    #[inline]
    pub fn new() -> Self {
        Self {
            qags: RefCell::new(None),
            qagp: RefCell::new(None)
        }
    }
}

impl<F: Integrand> Algorithm<F> for AUTO {
    fn integrate(
        &self,
        f: &mut F,
        interval: &Interval,
        config: &IntegrationConfig,
    ) -> IntegrationResult {
        if config.points.is_empty() {
            let mut qags = self.qags.borrow_mut();
            if qags.is_none() {
                *qags = Some(QAGS::new())
            }
            qags.as_ref().unwrap().integrate(f, interval, config)
        } else {
            let mut qagp = self.qagp.borrow_mut();
            if qagp.is_none() {
                *qagp = Some(QAGP::new())
            }
            qagp.as_ref().unwrap().integrate(f, interval, config)
        }
    }
}

extra_traits!(AUTO);
