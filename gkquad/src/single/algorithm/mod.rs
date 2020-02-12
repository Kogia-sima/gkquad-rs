//! Algorithms for 1-dimentional numerical integration
//!
//! # Finite Algorithms
//!
//! Default algorithms (AUTO, QAG, QAGS, QAGP) will be applicable to both finite
//! and infinite range, but will produce multiple assembly for a single
//! integrand, which results in large binary size.
//!
//! If you know the range is always finite, then you should use *_FINITE
//! algorithms.
//!
//! # References
//!
//! * [Numerical Integration — GNU GSL documentation](https://www.gnu.org/software/gsl/doc/html/integration.html)
//! * [Netlib quadpack library](http://www.netlib.org/quadpack/)

use super::common::{Integrand, IntegrationConfig, Range};
use crate::error::IntegrationResult;

/// 1-dimentional integration algorithm API
///
/// # Notes
///
/// This API is still unstable, and may changes dramatically in the future.
pub trait Algorithm<F: Integrand> {
    fn integrate(&self, f: &mut F, range: &Range, config: &IntegrationConfig) -> IntegrationResult;
}

macro_rules! extra_traits {
    ($name:ident) => {
        impl Default for $name {
            #[inline]
            fn default() -> Self {
                Self::new()
            }
        }

        impl std::fmt::Debug for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                f.write_str(stringify!($name))
            }
        }

        impl PartialEq<$name> for $name {
            #[inline]
            fn eq(&self, _: &Self) -> bool {
                true
            }
        }

        impl Eq for $name {}

        impl PartialOrd for $name {
            #[inline]
            fn partial_cmp(&self, _: &Self) -> Option<std::cmp::Ordering> {
                Some(std::cmp::Ordering::Equal)
            }
        }

        impl Ord for $name {
            #[inline]
            fn cmp(&self, _: &Self) -> std::cmp::Ordering {
                std::cmp::Ordering::Equal
            }
        }

        impl std::hash::Hash for $name {
            fn hash<H: std::hash::Hasher>(&self, _: &mut H) {}
        }
    };
}

mod qag;
pub use qag::*;

mod qag_finite;
pub use qag_finite::*;

mod qags;
pub use qags::*;

mod qags_finite;
pub use qags_finite::*;

mod qagp;
pub use qagp::*;

mod qagp_finite;
pub use qagp_finite::*;

mod auto;
pub use auto::*;
