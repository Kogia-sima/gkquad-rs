use super::error::RuntimeError;
use std::fmt::Debug;

/// Specify the tolerance which must be satisfied after calculation
///
/// It is recommended to Relative tolerance higher than 1e-14, since roundoff
/// error prevents the calculation error from converging into 0.
#[derive(Debug, Clone, PartialEq)]
pub enum Tolerance {
    /// `delta < absolute`
    Absolute(f64),
    /// `delta < relative * result`
    Relative(f64),
    /// `delta < max(absolute, relative * result)`
    AbsOrRel(f64, f64),
    /// `delta < min(absolute, relative * result)`
    AbsAndRel(f64, f64),
}

impl Tolerance {
    #[inline]
    pub(crate) fn contains_nan(&self) -> bool {
        match *self {
            Tolerance::Absolute(x) | Tolerance::Relative(x) if x.is_nan() => true,
            Tolerance::AbsOrRel(x, y) | Tolerance::AbsAndRel(x, y) if x.is_nan() || y.is_nan() => {
                true
            }
            _ => false,
        }
    }

    /// calculate the absolute tolerance from estimation
    #[inline]
    pub fn to_abs(&self, value: f64) -> f64 {
        match *self {
            Tolerance::Absolute(x) => x,
            Tolerance::Relative(y) => value * y,
            Tolerance::AbsOrRel(x, y) => f64::max(x, value * y),
            Tolerance::AbsAndRel(x, y) => f64::min(x, value * y),
        }
    }
}

impl Default for Tolerance {
    #[inline]
    fn default() -> Self {
        Tolerance::AbsOrRel(1.49e-8, 1.49e-8)
    }
}

/// `ValueWithError` is a type that holds both partial result (value) and failure
/// (error).
///
/// It is quite similar to `std::result::Result` type, but you can always
/// extract the partial result by calling `unwrap_unchecked()` method.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd, Eq, Ord, Hash)]
pub struct ValueWithError<T, E> {
    pub(crate) value: T,
    pub(crate) error: Option<E>,
}

impl<T, E> ValueWithError<T, E> {
    #[inline]
    pub fn new(value: T) -> Self {
        Self { value, error: None }
    }

    #[inline]
    pub fn with_error(value: T, error: E) -> Self {
        Self {
            value,
            error: Some(error),
        }
    }

    #[inline]
    pub fn has_err(&self) -> bool {
        self.error.is_some()
    }

    #[inline]
    pub fn ok(self) -> Option<T> {
        match self.error {
            Some(_) => None,
            None => Some(self.value),
        }
    }

    #[inline]
    pub fn err(self) -> Option<E> {
        self.error
    }

    #[inline]
    pub fn as_ref(&self) -> ValueWithError<&T, &E> {
        ValueWithError {
            value: &self.value,
            error: self.error.as_ref(),
        }
    }

    #[inline]
    pub fn as_mut(&mut self) -> ValueWithError<&mut T, &mut E> {
        ValueWithError {
            value: &mut self.value,
            error: self.error.as_mut(),
        }
    }

    #[inline]
    pub unsafe fn unwrap_unchecked(self) -> T {
        self.value
    }

    #[inline]
    pub fn unwrap_or(self, optb: T) -> T {
        match self.error {
            Some(_) => optb,
            None => self.value,
        }
    }

    #[inline]
    pub fn unwrap_or_else<F: FnOnce(E) -> T>(self, op: F) -> T {
        match self.error {
            Some(e) => op(e),
            None => self.value,
        }
    }

    #[inline]
    pub fn unwrap_err(self) -> E {
        self.expect_err("called `ValueWithError::unwrap_err() but it does not contains an error")
    }

    #[inline]
    pub fn expect_err(self, msg: &str) -> E {
        self.error.expect(msg)
    }
}

impl<T, E: Debug> ValueWithError<T, E> {
    #[inline]
    pub fn unwrap(self) -> T {
        self.expect("called `ValueWithError::unwrap()` but it contains an error")
    }

    pub fn expect(self, msg: &str) -> T {
        if let Some(e) = self.error {
            panic!("{}: {:?}", msg, &e);
        }

        self.value
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Solution {
    pub estimate: f64,
    pub delta: f64,
    pub nevals: usize,
}

impl Default for Solution {
    fn default() -> Self {
        Self {
            estimate: 0.0,
            delta: core::f64::MAX,
            nevals: 0,
        }
    }
}

pub type IntegrationResult = ValueWithError<Solution, RuntimeError>;
