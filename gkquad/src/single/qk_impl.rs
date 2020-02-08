use super::common::{Integrand, Interval};
use super::qk::QKResult;
use super::util::rescale_error;

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

    if n % 2 == 0 {
        result_gauss = f_center * wg.get_unchecked(0);
    }

    for j in 1..n {
        let fval1 = *buf.get_unchecked(j);
        let fval2 = *buf.get_unchecked(j + n - 1);
        let fsum = fval1 + fval2;
        result_kronrod += wgk.get_unchecked(j) * fsum;
        result_abs += wgk.get_unchecked(j) * (fval1.abs() + fval2.abs());

        if j % 2 == n % 2 {
            result_gauss += *wg.get_unchecked(j / 2) * fsum;
        }
    }

    let mean = result_kronrod * 0.5;
    let mut result_asc = wgk.get_unchecked(0) * (f_center - mean).abs();

    for j in 1..n {
        result_asc += wgk.get_unchecked(j)
            * ((buf.get_unchecked(j) - mean).abs() + (buf.get_unchecked(j + n - 1) - mean).abs());
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

impl_array!(4 8 5 11 16 10 21 13 26 15 31);
