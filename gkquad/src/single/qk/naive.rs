use super::super::common::{Integrand, Range};
use super::super::qk::QKResult;
use super::super::util::{rescale_error, Array};

#[cfg(not(feature = "std"))]
use crate::float::Float;

/// perform Gauss-Kronrod integration with custom points
///
/// # Parameters
///
/// * xgk: abscissae of the Gauss-Kronrod rule
/// * wg: weight of the Gauss rule
/// * wgk: weight of the Kronrod rule
/// * wck: Kronrod weight for center point
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
/// * range is finite
///
/// # TODO
///
/// - SIMD implementation
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
    debug_assert!(buf.len() > xgk.len() * 2);
    debug_assert!(range.begin.is_finite() && range.end.is_finite());

    let n = K::CAPACITY;
    let center = 0.5 * (range.begin + range.end);
    let half_length = 0.5 * (range.end - range.begin);
    let abs_half_length = half_length.abs();

    *buf.get_unchecked_mut(n << 1) = center;

    for j in 0..n {
        let abscissa = half_length * xgk.get_unchecked(j);
        *buf.get_unchecked_mut(j) = center - abscissa;
        *buf.get_unchecked_mut(j + n) = center + abscissa;
    }

    f.apply_to_slice(buf);

    let f_center = buf.get_unchecked(n << 1);
    let mut result_gauss = 0.;
    let mut result_kronrod = f_center * wck;
    let mut result_abs = result_kronrod.abs();

    for j in 0..n {
        let fval1 = *buf.get_unchecked(j);
        let fval2 = *buf.get_unchecked(j + n);
        let fsum = fval1 + fval2;
        result_kronrod += wgk.get_unchecked(j) * fsum;
        result_abs += wgk.get_unchecked(j) * (fval1.abs() + fval2.abs());

        if j % 2 == 0 {
            result_gauss += *wg.get_unchecked(j / 2) * fsum;
        }
    }

    let mean = result_kronrod * 0.5;
    let mut result_asc = wck * (f_center - mean).abs();

    for j in 0..n {
        result_asc += wgk.get_unchecked(j)
            * ((buf.get_unchecked(j) - mean).abs() + (buf.get_unchecked(j + n) - mean).abs());
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

#[inline(always)]
pub unsafe fn qk17<F: Integrand + ?Sized>(
    f: &mut F,
    range: &Range,
    xgk: &[f64; 8],
    wg: &[f64; 4],
    wgk: &[f64; 8],
    wck: f64,
    buf: &mut [f64; 17],
) -> QKResult {
    qk(f, range, xgk, wg, wgk, wck, buf)
}
