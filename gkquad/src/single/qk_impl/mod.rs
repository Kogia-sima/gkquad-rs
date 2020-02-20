#[cfg(all(feature = "simd", target_feature = "avx"))]
mod avx;
#[cfg(all(feature = "simd", target_feature = "avx"))]
pub use avx::qk;

#[cfg(not(all(feature = "simd", target_feature = "avx")))]
mod naive;
#[cfg(not(all(feature = "simd", target_feature = "avx")))]
pub use naive::qk;
