#[cfg(all(feature = "simd", target_feature = "avx"))]
#[path = "avx.rs"]
mod imp;

#[cfg(not(all(feature = "simd", target_feature = "avx")))]
#[path = "naive.rs"]
mod imp;

pub use imp::qk;
