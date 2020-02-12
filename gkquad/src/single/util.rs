use std::mem::MaybeUninit;
use std::ops::{Deref, DerefMut};

use crate::single::common::Interval;

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

impl_array!(4 6 8 10 12 14 16 20 24 28);

#[repr(align(32))]
pub struct Aligned<T: ?Sized> {
    value: T,
}

impl<T> Aligned<T> {
    pub const fn new(value: T) -> Aligned<T> {
        Aligned { value }
    }
}

impl<T: Sized> Aligned<T> {
    #[inline]
    pub unsafe fn uninit() -> Self {
        MaybeUninit::uninit().assume_init()
    }
}

impl<T: ?Sized> Deref for Aligned<T> {
    type Target = T;

    #[inline]
    fn deref(&self) -> &T {
        &self.value
    }
}

impl<T: ?Sized> DerefMut for Aligned<T> {
    #[inline]
    fn deref_mut(&mut self) -> &mut T {
        &mut self.value
    }
}

/// 区間幅が中央値の値に対して狭すぎる場合trueを返す
///
/// 例えば、区間[1e20, 1e20 + 1]は浮動小数点の桁落ちにより台形公式による分割を
/// 行ったときのxの値の誤差が大きくなるため、これ以上分割できない
#[inline]
pub fn subinterval_too_small(a1: f64, a2: f64, b2: f64) -> bool {
    let tmp = (1. + 100. * core::f64::EPSILON) * (a2.abs() + 1000. * core::f64::MIN_POSITIVE);

    // a1, b2は昇順とは限らないため両方チェックする:
    a1.abs() <= tmp && b2.abs() <= tmp
}

#[inline]
pub fn rescale_error(mut err: f64, result_abs: f64, result_asc: f64) -> f64 {
    err = err.abs();

    if result_asc != 0.0 && err != 0.0 {
        let mut scale = 200.0 * err;

        if scale < result_asc {
            scale /= result_asc;
            err = result_asc * scale * scale.sqrt();
        } else {
            err = result_asc;
        }
    }

    if result_abs > std::f64::MIN_POSITIVE / (50.0 * std::f64::EPSILON) {
        let min_err = 50.0 * std::f64::EPSILON * result_abs;

        if min_err > err {
            err = min_err;
        }
    };

    err
}

/// Compare the integral of f(x) with the integral of |f(x)| to determine if
/// f(x) covers both positive and negative values
#[inline]
pub fn test_positivity(result: f64, resabs: f64) -> bool {
    result.abs() >= (1.0 - 50.0 * core::f64::EPSILON) * resabs
}

#[inline]
pub fn bisect(interval: &Interval) -> (Interval, Interval) {
    let center = (interval.begin + interval.end) * 0.5;
    unsafe {
        (
            Interval::new_unchecked(interval.begin, center),
            Interval::new_unchecked(center, interval.end),
        )
    }
}

#[inline]
pub fn transform_point(x: f64) -> f64 {
    if x == std::f64::NEG_INFINITY {
        -1.0
    } else if x == std::f64::INFINITY {
        1.0
    } else {
        x / (1.0 + x.abs())
    }
}

// transform infinite interval to finite
#[inline]
pub fn transform_interval(interval: &Interval) -> Interval {
    unsafe {
        Interval::new_unchecked(
            transform_point(interval.begin),
            transform_point(interval.end),
        )
    }
}
