//! Numerical Integradion Library for Rust
//!
//! This crate offers utilities for calculating 1-dimentional integral by default.
//! Also, multi-dimentional integral toolbox are available (but experimental).
//!
//! Experimental utilities are not tested well, and then quite buggy at the moment.
//! If you found any bugs, feel free to create an issue or send a pull request.
//! This crate is looking for contributors.
//!
//! # Features
//!
//! `gkquad` aims to offer **easy**, **extensible**, and **fast** utilities.
//! Not only users can calculate the integral for various functions, but
//! developpers are able to write a new integration algorithms.

#![cfg_attr(docsrs, feature(doc_cfg, optin_builtin_traits))]
#![cfg_attr(not(feature = "std"), no_std)]
#![allow(deprecated)]

#[cfg(not(feature = "std"))]
extern crate core as std;

extern crate alloc;

#[cfg(all(
    not(target_arch = "x86"),
    not(target_arch = "x86_64"),
    feature = "simd"
))]
compile_error!("`simd` feature flag requires x86 or x86_64 architecture.");

mod common;
mod error;
mod utils;

pub use common::*;
pub use error::*;

pub mod single;

#[cfg(feature = "double")]
#[cfg_attr(docsrs, doc(cfg(feature = "double")))]
pub mod double;

pub mod prelude;

/// Deallocate the cache memory
///
/// `gkquad` keeps the allocated memory on the thread-local storage until program ends.
/// However, if you use this crate in the multi-thread application or set `limit`
/// option too high, then that storage leaks and occupies large amount of memory.
///
/// `free_memory` function will deallocate this storage of the `current thread`. If you
/// want to deallocate the storages on the all threads, you have to call this function
/// from all threads.
pub fn free_memory() {
    use single::borrow_workspace;

    let mut ws = borrow_workspace();
    ws.release();
}
