//! Performs 2-dimentional numerical integration

pub mod algorithm;
mod common;
mod integral;
mod integrator;
pub mod range;

pub use common::*;
pub use integral::*;
pub use integrator::*;
