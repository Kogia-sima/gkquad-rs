#![cfg_attr(docsrs, feature(doc_cfg))]

//! Perform general numeric integration

#![cfg_attr(not(feature = "std"), no_std)]

#[cfg(not(feature = "std"))]
extern crate core as std;

extern crate alloc;

mod common;
mod error;

pub use common::*;
pub use error::*;

pub mod single;
pub mod prelude;
