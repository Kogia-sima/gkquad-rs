use std::cell::RefCell;
use std::rc::Rc;

#[cfg(not(feature = "std"))]
use crate::float::Float;

use crate::error::IntegrationResult;
use crate::single::algorithm::{Algorithm, QAG_FINITE};
use crate::single::common::{ITransform, Integrand, IntegrationConfig, Interval};
use crate::single::util::transform_interval;
use crate::single::workspace::WorkSpace;

#[derive(Clone)]
pub struct QAG {
    #[doc(hidden)]
    pub workspace: Rc<RefCell<WorkSpace>>,
}

impl QAG {
    #[inline]
    pub fn new() -> Self {
        Self {
            workspace: Rc::new(RefCell::new(WorkSpace::new())),
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
        let qag_finite = QAG_FINITE {
            workspace: self.workspace.clone(),
        };

        if !interval.begin.is_finite() || !interval.end.is_finite() {
            let mut f = ITransform(f);
            let interval = transform_interval(interval);
            qag_finite.integrate(&mut f, &interval, config)
        } else {
            qag_finite.integrate(f, interval, config)
        }
    }
}
