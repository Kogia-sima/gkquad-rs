use smallvec::SmallVec;
use std::marker::PhantomData;
use std::ops::{Bound, RangeBounds};

use crate::Tolerance;

pub type Points = SmallVec<[f64; 8]>;

/// Represent the interval for which the integral is estimated.
///
/// Both `begin` and `end` are not NaN values (but may be infinite).
#[derive(Debug, Clone, PartialEq)]
pub struct Interval {
    /// beginning of the interval
    pub begin: f64,
    /// end of the interval
    pub end: f64,
    _private: PhantomData<()>,
}

impl Interval {
    /// Create a new `Interval` object
    ///
    /// Return `None` if either begin or end is NaN.
    #[inline]
    pub fn new(begin: f64, end: f64) -> Option<Interval> {
        if begin.is_nan() || end.is_nan() {
            None
        } else {
            unsafe { Some(Self::new_unchecked(begin, end)) }
        }
    }

    /// Create a new `Interval` object without NaN check
    #[inline]
    pub unsafe fn new_unchecked(begin: f64, end: f64) -> Interval {
        Interval {
            begin,
            end,
            _private: PhantomData,
        }
    }
}

impl Eq for Interval {}

impl<R: RangeBounds<f64>> From<R> for Interval {
    fn from(r: R) -> Interval {
        let a = match r.start_bound() {
            Bound::Excluded(&x) | Bound::Included(&x) => x.into(),
            Bound::Unbounded => std::f64::NEG_INFINITY,
        };

        let b = match r.end_bound() {
            Bound::Excluded(&x) | Bound::Included(&x) => x.into(),
            Bound::Unbounded => std::f64::INFINITY,
        };

        Interval::new(a, b)
            .expect("cannot create Interval object from Range which contains NaN value.")
    }
}

/// Integration configuration
#[non_exhaustive]
#[derive(Clone, Debug)]
pub struct IntegrationConfig {
    /// the tolerance to be satisfied
    pub tolerance: Tolerance,
    /// maximum number of subdivisions
    pub limit: usize,
    /// specify singular points
    pub points: Points,
}

impl Default for IntegrationConfig {
    #[inline]
    fn default() -> Self {
        Self {
            tolerance: Tolerance::default(),
            limit: 50,
            points: Points::new(),
        }
    }
}

/// The function that is to be integrated
pub trait Integrand {
    /// apply function to explanatory variable `x`
    fn apply(&mut self, x: f64) -> f64;

    /// apply function to each of the elements of `s`.
    fn apply_to_slice(&mut self, s: &mut [f64]) {
        s.iter_mut().for_each(|x| *x = self.apply(*x));
    }
}

impl<F: FnMut(f64) -> f64> Integrand for F {
    fn apply(&mut self, x: f64) -> f64 {
        (*self)(x)
    }
}

/// This struct is used in order to allow integration over infinite interval
pub(crate) struct ITransform<'a, F: Integrand>(pub &'a mut F);

impl<'a, F: Integrand> Integrand for ITransform<'a, F> {
    #[inline]
    fn apply(&mut self, x: f64) -> f64 {
        let coef = 1. / (1. - x.abs());
        let x2 = x * coef;
        self.0.apply(x2) * coef * coef
    }

    #[inline]
    fn apply_to_slice(&mut self, s: &mut [f64]) {
        s.iter_mut().for_each(|x| *x = self.apply(*x))
    }
}
