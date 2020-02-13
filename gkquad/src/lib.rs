//! Perform general numeric integration

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

pub mod prelude;

/// Deallocate the cache memory
pub fn free_memory() {
    use single::WorkSpaceProvider;

    let provider = WorkSpaceProvider::new();
    let mut ws = provider.get_mut();
    ws.release();
}
