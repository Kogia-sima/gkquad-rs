#[cfg(target_feature = "avx")]
mod imp {
    use super::super::common::{Integrand, Interval};
    use super::super::qk::QKResult;
    use super::super::util::rescale_error;
    use super::Array;

    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    #[inline]
    unsafe fn abs256(x: __m256d) -> __m256d {
        let bitmask = _mm256_set1_pd(f64::from_bits(0x7FFFFFFFFFFFFFFF));
        _mm256_and_pd(x, bitmask)
    }

    #[inline]
    unsafe fn sum256(x: __m256d) -> f64 {
        let mut vlow = _mm256_castpd256_pd128(x);
        let vhigh = _mm256_extractf128_pd(x, 1);
        vlow = _mm_add_pd(vlow, vhigh);
        let high64 = _mm_unpackhi_pd(vlow, vlow);
        _mm_cvtsd_f64(_mm_add_sd(vlow, high64))
    }

    #[cfg(target_feature = "fma")]
    #[inline]
    unsafe fn fmadd256(a: __m256d, b: __m256d, c: __m256d) -> __m256d {
        _mm256_fmadd_pd(a, b, c)
    }

    #[cfg(target_feature = "fma")]
    #[inline]
    unsafe fn fmadd128(a: __m128d, b: __m128d, c: __m128d) -> __m128d {
        _mm_fmadd_pd(a, b, c)
    }

    #[cfg(not(target_feature = "fma"))]
    #[inline]
    unsafe fn fmadd256(a: __m256d, b: __m256d, c: __m256d) -> __m256d {
        _mm256_add_pd(_mm256_mul_pd(a, b), c)
    }

    #[cfg(not(target_feature = "fma"))]
    #[inline]
    unsafe fn fmadd128(a: __m128d, b: __m128d, c: __m128d) -> __m128d {
        _mm_add_pd(_mm_mul_pd(a, b), c)
    }

    #[cfg(target_feature = "avx2")]
    #[inline]
    unsafe fn compact(x: __m256d) -> __m128d {
        _mm256_castpd256_pd128(_mm256_permute4x64_pd(x, 0x88))
    }

    #[cfg(not(target_feature = "avx2"))]
    #[inline]
    unsafe fn compact(x: __m256d) -> __m128d {
        let vlow = _mm256_castpd256_pd128(x);
        let vhigh = _mm256_extractf128_pd(x, 1);
        _mm_unpacklo_pd(vlow, vhigh)
    }

    pub unsafe fn qk<F, K, G>(
        f: &mut F,
        range: &Interval,
        xgk: &K,
        wg: &G,
        wgk: &K,
        buf: &mut [f64],
    ) -> QKResult
    where
        F: Integrand,
        K: Array<Item = f64>,
        G: Array<Item = f64>,
    {
        let xgk = xgk.as_slice();
        let wg = wg.as_slice();
        let wgk = wgk.as_slice();

        debug_assert!(!xgk.is_empty());
        debug_assert!(xgk.len() / 2 == wg.len());
        debug_assert!(buf.len() >= xgk.len() * 2 - 1);

        let n = K::CAPACITY;
        let center = 0.5 * (range.begin + range.end);
        let half_length = 0.5 * (range.end - range.begin);
        let abs_half_length = half_length.abs();

        let xgkp = xgk.as_ptr();
        let wgp = wg.as_ptr();
        let wgkp = wgk.as_ptr();
        let bufp = buf.as_mut_ptr();
        let bufp2 = bufp.add(n - 1);

        *bufp = center;

        // The following lines may prevent optimization by compiler

        // let center_simd = _mm256_set1_pd(center);
        // let half_length_simd = _mm256_set1_pd(half_length);

        // for j in (1..n).step_by(4) {
        //     let abscissa = _mm256_mul_pd(half_length_simd, _mm256_loadu_pd(xgkp.add(j)));
        //     _mm256_storeu_pd(bufp.add(j), _mm256_sub_pd(center_simd, abscissa));
        //     _mm256_storeu_pd(bufp2.add(j), _mm256_add_pd(center_simd, abscissa));
        // }

        for j in 1..n {
            let abscissa = half_length * *xgkp.add(j);
            *bufp.add(j) = center - abscissa;
            *bufp2.add(j) = center + abscissa;
        }

        f.apply_to_slice(buf);

        let f_center = *bufp;
        let tmp = f_center * *wgkp;
        let mut result_gauss = _mm_setzero_pd();
        let mut result_kronrod = _mm256_set_pd(tmp, 0., 0., 0.);
        let mut result_abs = _mm256_set_pd(tmp.abs(), 0., 0., 0.);

        for j in (1..n).step_by(4) {
            let fval1 = _mm256_loadu_pd(bufp.add(j));
            let fval2 = _mm256_loadu_pd(bufp2.add(j));
            let fsum = _mm256_add_pd(fval1, fval2);
            let abssum = _mm256_add_pd(abs256(fval1), abs256(fval2));
            let wgk = _mm256_loadu_pd(wgkp.add(j));
            result_kronrod = fmadd256(wgk, fsum, result_kronrod);
            result_abs = fmadd256(wgk, abssum, result_abs);

            let fsum_compact = compact(fsum);
            let wg = _mm_loadu_pd(wgp.add(j >> 1));
            result_gauss = fmadd128(wg, fsum_compact, result_gauss);
        }

        let mut result_kronrod = sum256(result_kronrod);
        let mut result_abs = sum256(result_abs);
        let result_gauss = _mm_cvtsd_f64(_mm_hadd_pd(result_gauss, result_gauss));

        let mean = result_kronrod * 0.5;
        let mean_simd = _mm256_set1_pd(mean);
        let mut result_asc = _mm256_set_pd(*wgkp * (f_center - mean).abs(), 0., 0., 0.);

        for j in (1..n).step_by(4) {
            let fval1 = _mm256_loadu_pd(bufp.add(j));
            let fval2 = _mm256_loadu_pd(bufp2.add(j));
            let diff1 = abs256(_mm256_sub_pd(fval1, mean_simd));
            let diff2 = abs256(_mm256_sub_pd(fval2, mean_simd));
            let diff_sum = _mm256_add_pd(diff1, diff2);
            let wgk = _mm256_loadu_pd(wgkp.add(j));
            result_asc = fmadd256(wgk, diff_sum, result_asc);
        }

        let mut result_asc = sum256(result_asc);
        let err = (result_kronrod - result_gauss) * half_length;

        result_kronrod *= half_length;
        result_abs *= abs_half_length;
        result_asc *= abs_half_length;

        QKResult {
            estimate: result_kronrod,
            delta: rescale_error(err, result_abs, result_asc),
            absvalue: result_abs,
            asc: result_asc,
        }
    }
}

#[cfg(not(target_feature = "avx"))]
mod imp {
    use super::super::common::{Integrand, Interval};
    use super::super::qk::QKResult;
    use super::super::util::rescale_error;
    use super::Array;

    /// perform Gauss-Kronrod integration with custom points
    ///
    /// # Parameters
    ///
    /// * xgk: abscissae of the Gauss-Kronrod rule
    /// * wg: weight of the Gauss rule
    /// * wgk: weight of the Kronrod rule
    /// * buf: buffer to hold return value of integrand, may be not initialized.
    ///
    /// # Precondition
    ///
    /// The arguments must satisfy the following conditions, otherwise this
    /// function causes an undefined behaviour.
    ///
    /// * xgk.len() > 0
    /// * xgk.len() / 2 == wg.len()
    /// * buf.len() >= xgk.len() * 2 - 1
    ///
    /// # TODO
    ///
    /// - SIMD implementation
    pub unsafe fn qk<F, K, G>(
        f: &mut F,
        range: &Interval,
        xgk: &K,
        wg: &G,
        wgk: &K,
        buf: &mut [f64],
    ) -> QKResult
    where
        F: Integrand,
        K: Array<Item = f64>,
        G: Array<Item = f64>,
    {
        let xgk = xgk.as_slice();
        let wg = wg.as_slice();
        let wgk = wgk.as_slice();

        debug_assert!(!xgk.is_empty());
        debug_assert!(xgk.len() / 2 == wg.len());
        debug_assert!(buf.len() >= xgk.len() * 2 - 1);

        let n = K::CAPACITY;
        let center = 0.5 * (range.begin + range.end);
        let half_length = 0.5 * (range.end - range.begin);
        let abs_half_length = half_length.abs();

        *buf.get_unchecked_mut(0) = center;

        for j in 1..n {
            let abscissa = half_length * xgk.get_unchecked(j);
            *buf.get_unchecked_mut(j) = center - abscissa;
            *buf.get_unchecked_mut(j + n - 1) = center + abscissa;
        }

        f.apply_to_slice(buf);

        let f_center = buf.get_unchecked(0);
        let mut result_gauss = 0.;
        let mut result_kronrod = f_center * wgk.get_unchecked(0);
        let mut result_abs = result_kronrod.abs();

        for j in 1..n {
            let fval1 = *buf.get_unchecked(j);
            let fval2 = *buf.get_unchecked(j + n - 1);
            let fsum = fval1 + fval2;
            result_kronrod += wgk.get_unchecked(j) * fsum;
            result_abs += wgk.get_unchecked(j) * (fval1.abs() + fval2.abs());

            if j % 2 == 1 {
                result_gauss += *wg.get_unchecked(j / 2) * fsum;
            }
        }

        let mean = result_kronrod * 0.5;
        let mut result_asc = wgk.get_unchecked(0) * (f_center - mean).abs();

        for j in 1..n {
            result_asc += wgk.get_unchecked(j)
                * ((buf.get_unchecked(j) - mean).abs()
                    + (buf.get_unchecked(j + n - 1) - mean).abs());
        }

        let err = (result_kronrod - result_gauss) * half_length;
        result_kronrod *= half_length;
        result_abs *= abs_half_length;
        result_asc *= abs_half_length;

        QKResult {
            estimate: result_kronrod,
            delta: rescale_error(err, result_abs, result_asc),
            absvalue: result_abs,
            asc: result_asc,
        }
    }
}

pub use imp::qk;

pub trait Array {
    type Item;
    const CAPACITY: usize;

    fn as_slice(&self) -> &[Self::Item];
    fn as_mut_slice(&mut self) -> &mut [Self::Item];
}

macro_rules! impl_array {
    ($($N:tt)*) => {
        $(
            impl Array for [f64; $N] {
                type Item = f64;
                const CAPACITY: usize = $N;

                #[inline]
                fn as_slice(&self) -> &[f64] {
                    self.as_ref()
                }

                #[inline]
                fn as_mut_slice(&mut self) -> &mut [f64] {
                    self.as_mut()
                }
            }
        )*
    };
}

impl_array!(4 6 8 10 12 14 9 13 17 21 25 29);
