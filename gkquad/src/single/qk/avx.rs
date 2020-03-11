use super::super::common::{Integrand, Range};
use super::super::qk::QKResult;
use super::super::util::{rescale_error, Array};

#[cfg(target_arch = "x86")]
use core::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use core::arch::x86_64::*;

#[inline(always)]
unsafe fn abs256(x: __m256d) -> __m256d {
    let bitmask = _mm256_set1_pd(f64::from_bits(0x7FFFFFFFFFFFFFFF));
    _mm256_and_pd(x, bitmask)
}

#[inline(always)]
unsafe fn sum256(x: __m256d) -> f64 {
    let mut vlow = _mm256_castpd256_pd128(x);
    let vhigh = _mm256_extractf128_pd(x, 1);
    vlow = _mm_add_pd(vlow, vhigh);
    let high64 = _mm_unpackhi_pd(vlow, vlow);
    _mm_cvtsd_f64(_mm_add_sd(vlow, high64))
}

#[cfg(target_feature = "fma")]
#[inline(always)]
unsafe fn fmadd256(a: __m256d, b: __m256d, c: __m256d) -> __m256d {
    _mm256_fmadd_pd(a, b, c)
}

#[cfg(target_feature = "fma")]
#[inline(always)]
unsafe fn fmadd128(a: __m128d, b: __m128d, c: __m128d) -> __m128d {
    _mm_fmadd_pd(a, b, c)
}

#[cfg(not(target_feature = "fma"))]
#[inline(always)]
unsafe fn fmadd256(a: __m256d, b: __m256d, c: __m256d) -> __m256d {
    _mm256_add_pd(_mm256_mul_pd(a, b), c)
}

#[cfg(not(target_feature = "fma"))]
#[inline(always)]
unsafe fn fmadd128(a: __m128d, b: __m128d, c: __m128d) -> __m128d {
    _mm_add_pd(_mm_mul_pd(a, b), c)
}

#[inline(always)]
unsafe fn compact(x: __m256d) -> __m128d {
    let vlow = _mm256_castpd256_pd128(x);
    let vhigh = _mm256_extractf128_pd(x, 1);
    _mm_unpacklo_pd(vlow, vhigh)
}

#[inline(always)]
pub unsafe fn qk<F, K, G>(
    f: &mut F,
    range: &Range,
    xgk: &K,
    wg: &G,
    wgk: &K,
    wck: f64,
    buf: &mut [f64],
) -> QKResult
where
    F: Integrand + ?Sized,
    K: Array<Item = f64>,
    G: Array<Item = f64>,
{
    let xgk = xgk.as_slice();
    let wg = wg.as_slice();
    let wgk = wgk.as_slice();

    debug_assert!(!xgk.is_empty());
    debug_assert!(xgk.len() == wg.len() * 2);
    debug_assert!(buf.len() >= xgk.len() * 2 + 1);
    debug_assert!(range.begin.is_finite() && range.end.is_finite());

    let n = K::CAPACITY;
    let center = 0.5 * (range.begin + range.end);
    let half_length = 0.5 * (range.end - range.begin);
    let abs_half_length = half_length.abs();

    let xgkp = xgk.as_ptr();
    let wgp = wg.as_ptr();
    let wgkp = wgk.as_ptr();
    let bufp = buf.as_mut_ptr();

    let center_simd = _mm256_set1_pd(center);
    let half_length_simd = _mm256_set1_pd(half_length);

    for j in (0..n).step_by(4) {
        let abscissa = _mm256_mul_pd(half_length_simd, _mm256_load_pd(xgkp.add(j)));
        _mm256_store_pd(bufp.add(j), _mm256_sub_pd(center_simd, abscissa));
        _mm256_store_pd(bufp.add(j + n), _mm256_add_pd(center_simd, abscissa));
    }

    *bufp.add(n << 1) = center;

    f.apply_to_slice(buf);

    let mut result_gauss = _mm_setzero_pd();
    let mut result_kronrod = _mm256_setzero_pd();
    let mut result_abs = _mm256_setzero_pd();

    for j in (0..n).step_by(4) {
        let fval1 = _mm256_load_pd(bufp.add(j));
        let fval2 = _mm256_load_pd(bufp.add(j + n));
        let fsum = _mm256_add_pd(fval1, fval2);
        let abssum = _mm256_add_pd(abs256(fval1), abs256(fval2));
        let wgk = _mm256_load_pd(wgkp.add(j));
        result_kronrod = fmadd256(wgk, fsum, result_kronrod);
        result_abs = fmadd256(wgk, abssum, result_abs);

        let fsum_compact = compact(fsum);
        let wg = _mm_load_pd(wgp.add(j >> 1));
        result_gauss = fmadd128(wg, fsum_compact, result_gauss);
    }

    let f_center = *bufp.add(n << 1);
    let tmp = f_center * wck;
    let mut result_kronrod = sum256(result_kronrod) + tmp;
    let mut result_abs = sum256(result_abs) + tmp.abs();
    let result_gauss = _mm_cvtsd_f64(_mm_hadd_pd(result_gauss, result_gauss));

    let mean = result_kronrod * 0.5;
    let mean_simd = _mm256_set1_pd(mean);
    let mut result_asc = _mm256_set_pd(wck * (f_center - mean).abs(), 0., 0., 0.);

    for j in (0..n).step_by(4) {
        let fval1 = _mm256_load_pd(bufp.add(j));
        let fval2 = _mm256_load_pd(bufp.add(j + n));
        let diff1 = abs256(_mm256_sub_pd(fval1, mean_simd));
        let diff2 = abs256(_mm256_sub_pd(fval2, mean_simd));
        let diff_sum = _mm256_add_pd(diff1, diff2);
        let wgk = _mm256_load_pd(wgkp.add(j));
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

#[inline(always)]
pub unsafe fn qk17<F>(
    f: &mut F,
    range: &Range,
    xgk: &[f64; 8],
    wg: &[f64; 4],
    wgk: &[f64; 8],
    wck: f64,
    buf: &mut [f64; 17],
) -> QKResult
where
    F: Integrand + ?Sized,
{
    debug_assert!(range.begin.is_finite() && range.end.is_finite());

    let center = 0.5 * (range.begin + range.end);
    let half_length = 0.5 * (range.end - range.begin);
    let abs_half_length = half_length.abs();

    let center_simd = _mm256_set1_pd(center);
    let half_length_simd = _mm256_set1_pd(half_length);

    {
        let abscissa = _mm256_mul_pd(half_length_simd, _mm256_load_pd(&xgk[0]));
        _mm256_store_pd(&mut buf[0], _mm256_sub_pd(center_simd, abscissa));
        _mm256_store_pd(&mut buf[8], _mm256_add_pd(center_simd, abscissa));
    }

    {
        let abscissa = _mm256_mul_pd(half_length_simd, _mm256_load_pd(&xgk[4]));
        _mm256_store_pd(&mut buf[4], _mm256_sub_pd(center_simd, abscissa));
        _mm256_store_pd(&mut buf[12], _mm256_add_pd(center_simd, abscissa));
    }

    buf[16] = center;

    f.apply_to_slice(buf);
    // buf.iter_mut().for_each(|v| *v = f.apply(*v));

    let mut result_kronrod: __m256d;
    let mut result_abs: __m256d;
    let mut result_gauss: __m128d;

    {
        let fval1 = _mm256_load_pd(&buf[0]);
        let fval2 = _mm256_load_pd(&buf[8]);
        let fsum = _mm256_add_pd(fval1, fval2);
        let abssum = _mm256_add_pd(abs256(fval1), abs256(fval2));
        let wgk_simd = _mm256_load_pd(&wgk[0]);
        result_kronrod = _mm256_mul_pd(wgk_simd, fsum);
        result_abs = _mm256_mul_pd(wgk_simd, abssum);

        let fsum_compact = compact(fsum);
        let wg_simd = _mm_load_pd(&wg[0]);
        result_gauss = _mm_mul_pd(wg_simd, fsum_compact);
    }

    {
        let fval1 = _mm256_load_pd(&buf[4]);
        let fval2 = _mm256_load_pd(&buf[12]);
        let fsum = _mm256_add_pd(fval1, fval2);
        let abssum = _mm256_add_pd(abs256(fval1), abs256(fval2));
        let wgk_simd = _mm256_load_pd(&wgk[4]);
        result_kronrod = fmadd256(wgk_simd, fsum, result_kronrod);
        result_abs = fmadd256(wgk_simd, abssum, result_abs);

        let fsum_compact = compact(fsum);
        let wg_simd = _mm_load_pd(&wg[2]);
        result_gauss = fmadd128(wg_simd, fsum_compact, result_gauss);
    }

    let f_center = buf[16];
    let tmp = f_center * wck;
    let mut result_kronrod = sum256(result_kronrod) + tmp;
    let mut result_abs = sum256(result_abs) + tmp.abs();
    let result_gauss = _mm_cvtsd_f64(_mm_hadd_pd(result_gauss, result_gauss));

    let mean = result_kronrod * 0.5;
    let mean_simd = _mm256_set1_pd(mean);
    let mut result_asc = _mm256_set_pd(wck * (f_center - mean).abs(), 0., 0., 0.);

    {
        let fval1 = _mm256_load_pd(&buf[0]);
        let fval2 = _mm256_load_pd(&buf[8]);
        let diff1 = abs256(_mm256_sub_pd(fval1, mean_simd));
        let diff2 = abs256(_mm256_sub_pd(fval2, mean_simd));
        let diff_sum = _mm256_add_pd(diff1, diff2);
        let wgk = _mm256_load_pd(&wgk[0]);
        result_asc = fmadd256(wgk, diff_sum, result_asc);
    }

    {
        let fval1 = _mm256_load_pd(&buf[4]);
        let fval2 = _mm256_load_pd(&buf[12]);
        let diff1 = abs256(_mm256_sub_pd(fval1, mean_simd));
        let diff2 = abs256(_mm256_sub_pd(fval2, mean_simd));
        let diff_sum = _mm256_add_pd(diff1, diff2);
        let wgk = _mm256_load_pd(&wgk[4]);
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
