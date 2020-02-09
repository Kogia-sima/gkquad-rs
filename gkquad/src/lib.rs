#![cfg_attr(docsrs, feature(doc_cfg, optin_builtin_traits))]

//! Perform general numeric integration

#![cfg_attr(not(feature = "std"), no_std)]

#[cfg(not(feature = "std"))]
extern crate core as std;

extern crate alloc;

mod common;
mod error;

pub use common::*;
pub use error::*;

#[cfg(feature = "single")]
#[cfg_attr(docsrs, doc(cfg(feature = "single")))]
pub mod single;

pub mod prelude;
