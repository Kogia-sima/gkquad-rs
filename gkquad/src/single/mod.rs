//! Performs 1-dimentional numerical integration
//!
//! ## Examples
//!
//! In order to calculate the simple integral, you can use `integral` function.
//!
//! ```
//! use gkquad::single::integral;
//!
//! // calculate the integral over range (0.0, 1.0)
//! let result = integral(|x: f64| x.sin() * (- x * x).exp(), 0.0..1.0)
//!         .estimate()
//!         .unwrap();
//! ```
//!
//! If you want to calculate more complicated integral, you can use `Integrator` object.
//!
//! ```
//! use std::f64::{INFINITY, NEG_INFINITY};
//!
//! use gkquad::single::Integrator;
//! use gkquad::single::algorithm::QAGP;
//!
//! // calculate the integral over range (-∞, ∞), with QAGP algorithm,
//! // maximum iteration limit being 100, singular point on the origin of coordinate.
//! let result = Integrator::new(|x: f64| 1. - (-(x.abs() / 1.6).powf(-2.3)).exp())
//!         .max_iters(100)
//!         .points(&[0.])
//!         .run(NEG_INFINITY..INFINITY)
//!         .estimate()
//!         .unwrap();
//! ```

pub mod algorithm;
mod common;
mod integral;
mod integrator;
mod qelg;
mod qk;
mod qk_impl;
mod util;
mod workspace;

pub use common::*;
pub use integral::*;
pub use integrator::*;
pub use qk::*;
#[doc(hidden)]
pub use workspace::*;
