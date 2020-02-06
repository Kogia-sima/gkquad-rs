#[cfg(not(feature = "std"))]
use crate::float::Float;

use std::cell::RefCell;
use std::rc::Rc;

use crate::error::IntegrationResult;
use crate::single::algorithm::*;
use crate::single::common::{Integrand, IntegrationConfig, Interval};
use crate::single::workspace::WorkSpace;

#[derive(Clone)]
pub struct AUTO {
    #[doc(hidden)]
    pub workspace: Rc<RefCell<WorkSpace>>,
}

impl AUTO {
    #[inline]
    pub fn new() -> Self {
        Self {
            workspace: Rc::new(RefCell::new(WorkSpace::new())),
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
            let qags = QAGS {
                workspace: self.workspace.clone(),
            };
            qags.integrate(f, interval, config)
        } else {
            let qagp = QAGP {
                workspace: self.workspace.clone(),
            };
            qagp.integrate(f, interval, config)
        }
    }
}
