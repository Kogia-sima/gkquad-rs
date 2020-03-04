//! Algorithms for 1-dimentional numerical integration
//!
//! # References
//!
//! * [Numerical Integration — GNU GSL documentation](https://www.gnu.org/software/gsl/doc/html/integration.html)
//! * [Netlib quadpack library](http://www.netlib.org/quadpack/)

use super::common::{Integrand, IntegrationConfig, Range};
use super::workspace::WorkSpace;
use crate::common::IntegrationResult;

/// 1-dimentional integration algorithm API
///
/// # Notes
///
/// This API is still unstable, and may changes dramatically in the future.
pub trait Algorithm<F: Integrand> {
    fn integrate(
        &mut self,
        f: &mut F,
        range: &Range,
        config: &IntegrationConfig,
    ) -> IntegrationResult;
}

#[doc(hidden)]
pub trait AlgorithmWithWorkSpace {
    fn workspace(&self) -> &WorkSpace;
}

macro_rules! extra_traits {
    ($name:ident) => {
        extra_traits!(@INNER $name [<>]);
    };
    ($name:ident<$($lifetimes:tt>),*) => {
        extra_traits!(@INNER $name [<$($lifetimes),*>]);
    };
    (@INNER $name:ident [$($lifetimes:tt)*]) => {
        impl $($lifetimes)* Default for $name $($lifetimes)* {
            #[inline]
            fn default() -> Self {
                Self::new()
            }
        }

        impl $($lifetimes)* core::fmt::Debug for $name $($lifetimes)* {
            fn fmt(&self, f: &mut core::fmt::Formatter) -> core::fmt::Result {
                f.write_str(stringify!($name))
            }
        }

        impl $($lifetimes)* PartialEq<$name $($lifetimes)*> for $name $($lifetimes)* {
            #[inline]
            fn eq(&self, _: &Self) -> bool {
                true
            }
        }

        impl $($lifetimes)* Eq for $name $($lifetimes)* {}

        impl $($lifetimes)* PartialOrd for $name $($lifetimes)* {
            #[inline]
            fn partial_cmp(&self, _: &Self) -> Option<core::cmp::Ordering> {
                Some(core::cmp::Ordering::Equal)
            }
        }

        impl $($lifetimes)* Ord for $name $($lifetimes)* {
            #[inline]
            fn cmp(&self, _: &Self) -> core::cmp::Ordering {
                core::cmp::Ordering::Equal
            }
        }

        impl $($lifetimes)* core::hash::Hash for $name $($lifetimes)* {
            fn hash<H: core::hash::Hasher>(&self, _: &mut H) {}
        }

        #[cfg(feature = "std")]
        impl $($lifetimes)* std::panic::UnwindSafe for $name $($lifetimes)* {}
    };
}

mod qag;
pub use qag::*;

mod qags;
pub use qags::*;

mod qagp;
pub use qagp::*;

mod auto;
pub use auto::*;
