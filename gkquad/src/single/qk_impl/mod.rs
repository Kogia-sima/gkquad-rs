#[cfg(target_feature = "avx")]
#[path = "avx.rs"]
mod imp;

#[cfg(not(target_feature = "avx"))]
#[path = "naive.rs"]
mod imp;

pub use imp::qk;
