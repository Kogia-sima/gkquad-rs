//! Performs 1-dimentional numerical integration

pub mod algorithm;
mod common;
mod integrate;
mod integrator;
mod workspace;
mod util;
mod qk_impl;
mod qk;
mod qelg;

pub use common::*;
pub use integrate::*;
pub use integrator::*;
pub use qk::*;
#[doc(hidden)]
pub use workspace::*;
