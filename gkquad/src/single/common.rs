#![allow(clippy::derive_hash_xor_eq)]

use smallvec::SmallVec;
use std::fmt::{self, Debug, Display};
use std::hash::{Hash, Hasher};
use std::marker::PhantomData;
use std::ops::{Bound, RangeBounds};

use crate::Tolerance;

/// Singular points
pub type Points = SmallVec<[f64; 8]>;

/// Represent the range for which the integral is estimated.
///
/// Both `begin` and `end` are not NaN values (but may be infinite).
#[derive(Clone, PartialEq)]
pub struct Range {
    /// beginning of the range
    pub begin: f64,
    /// end of the range
    pub end: f64,
    _private: PhantomData<()>,
}

impl Range {
    /// Create a new `Range` object
    ///
    /// Return `None` if either begin or end is NaN.
    #[inline]
    pub fn new(begin: f64, end: f64) -> Option<Range> {
        if begin.is_nan() || end.is_nan() {
            None
        } else {
            unsafe { Some(Self::new_unchecked(begin, end)) }
        }
    }

    /// Create a new `Range` object without NaN check
    ///
    /// # Safety
    ///
    /// Arguments must not be a NaN value, otherwise causes an undefined behaviour
    #[inline]
    pub unsafe fn new_unchecked(begin: f64, end: f64) -> Range {
        Range {
            begin,
            end,
            _private: PhantomData,
        }
    }
}

impl Eq for Range {}

impl Display for Range {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[{:?}, {:?}]", self.begin, self.end)
    }
}

impl Debug for Range {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        <Self as Display>::fmt(self, f)
    }
}

impl Hash for Range {
    #[inline]
    fn hash<H: Hasher>(&self, h: &mut H) {
        for &x in &[self.begin, self.end] {
            let bits = if x == 0.0 {
                0 // this accounts for +0.0 and -0.0
            } else {
                x.to_bits()
            };
            bits.hash(h);
        }
    }
}

impl<R: RangeBounds<f64>> From<R> for Range {
    fn from(r: R) -> Range {
        let a = match r.start_bound() {
            Bound::Excluded(&x) | Bound::Included(&x) => x,
            Bound::Unbounded => std::f64::NEG_INFINITY,
        };

        let b = match r.end_bound() {
            Bound::Excluded(&x) | Bound::Included(&x) => x,
            Bound::Unbounded => std::f64::INFINITY,
        };

        Range::new(a, b).expect("cannot create Range object from Range which contains NaN value.")
    }
}

impl<'a> From<&'a Range> for Range {
    fn from(other: &'a Range) -> Self {
        other.clone()
    }
}

/// Integration configuration
#[non_exhaustive]
#[derive(Clone, Debug, PartialEq)]
pub struct IntegrationConfig {
    /// the tolerance to be satisfied
    pub tolerance: Tolerance,
    /// maximum number of subdivisions
    pub max_iters: usize,
    /// specify singular points
    pub points: Points,
}

impl Default for IntegrationConfig {
    #[inline]
    fn default() -> Self {
        Self {
            tolerance: Tolerance::default(),
            max_iters: 50,
            points: Points::new(),
        }
    }
}

/// The function that is to be integrated
pub trait Integrand {
    /// apply function to explanatory variable `x`
    fn apply(&mut self, x: f64) -> f64;

    /// apply function to each of the elements of `s`.
    #[inline]
    fn apply_to_slice(&mut self, s: &mut [f64]) {
        s.iter_mut().for_each(|x| *x = self.apply(*x));
    }
}

impl<F: FnMut(f64) -> f64> Integrand for F {
    #[inline]
    fn apply(&mut self, x: f64) -> f64 {
        (*self)(x)
    }
}
