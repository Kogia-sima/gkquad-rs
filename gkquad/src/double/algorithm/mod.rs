//! Algorithms for 2-dimentional numerical integration

use super::common::{Integrand2, IntegrationConfig2};
use super::range::Range2;
use crate::common::IntegrationResult;

pub trait Algorithm2<F: Integrand2 + ?Sized, R: Range2> {
    fn integrate(&mut self, f: &mut F, range: &R, config: &IntegrationConfig2)
        -> IntegrationResult;
}

macro_rules! extra_traits {
    ($name:ident) => {
        impl Default for $name {
            #[inline]
            fn default() -> Self {
                Self::new()
            }
        }

        impl core::fmt::Debug for $name {
            fn fmt(&self, f: &mut core::fmt::Formatter) -> core::fmt::Result {
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
            fn partial_cmp(&self, _: &Self) -> Option<core::cmp::Ordering> {
                Some(core::cmp::Ordering::Equal)
            }
        }

        impl Ord for $name {
            #[inline]
            fn cmp(&self, _: &Self) -> core::cmp::Ordering {
                core::cmp::Ordering::Equal
            }
        }

        impl core::hash::Hash for $name {
            fn hash<H: core::hash::Hasher>(&self, _: &mut H) {}
        }
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
