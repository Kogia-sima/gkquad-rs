//! Algorithms for 1-dimentional numerical integration
//!
//! Reference:
//! * [Numerical Integration — GNU GSL documentation](https://www.gnu.org/software/gsl/doc/html/integration.html)
//! * [Netlib quadpack library](http://www.netlib.org/quadpack/)

use std::any::Any;

use super::common::{Integrand, IntegrationConfig, Interval};
use crate::error::IntegrationResult;

/// 1-dimentional integration algorithm API
///
/// # Notes
///
/// This API is still unstable, and may changes dramatically in the future.
pub trait Algorithm<F: Integrand>: Downcast {
    fn integrate(
        &self,
        f: &mut F,
        interval: &Interval,
        config: &IntegrationConfig,
    ) -> IntegrationResult;

    #[doc(hidden)]
    fn get_workspace(&self) -> Option<std::cell::Ref<super::workspace::WorkSpace>> {
        None
    }
}

#[doc(hidden)]
pub trait Downcast: Any {
    fn as_any(&self) -> &dyn Any;
    fn as_any_mut(&mut self) -> &mut dyn Any;
}

impl<T: Any> Downcast for T {
    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

impl<F> dyn Algorithm<F>
where F: Integrand + Any + 'static {
    pub fn downcast_ref<T: Algorithm<F>>(&self) -> Option<&T> {
        Downcast::as_any(self).downcast_ref::<T>()
    }

    pub fn downcast_mut<T: Algorithm<F>>(&mut self) -> Option<&mut T> {
        Downcast::as_any_mut(self).downcast_mut::<T>()
    }
}

mod qag;
pub use qag::*;

mod qags;
pub use qags::*;

mod qagp;
pub use qagp::*;

mod auto;
pub use auto::*;
