//! Algorithms for 2-dimentional numerical integration

use super::common::{Integrand2, IntegrationConfig2, Range2};
use crate::IntegrationResult;

pub trait Algorithm2<F: Integrand2> {
    fn integrate(
        &self,
        f: &mut F,
        range: &Range2,
        config: &IntegrationConfig2,
    ) -> IntegrationResult;
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

mod qags;
pub use qags::*;

mod qagp;
pub use qagp::*;
