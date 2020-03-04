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
///
/// This type is useful when the calculation must satisfy a some post-condition
/// which is checked at runtime.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd, Eq, Ord, Hash)]
pub struct ValueWithError<T, E> {
    pub(crate) value: T,
    pub(crate) error: Option<E>,
}

impl<T, E> ValueWithError<T, E> {
    /// Create a new instance with no error
    #[inline]
    pub fn new(value: T) -> Self {
        Self { value, error: None }
    }

    /// Create a new instance with given error
    #[inline]
    pub fn with_error(value: T, error: E) -> Self {
        Self {
            value,
            error: Some(error),
        }
    }

    /// Return true if the instance contains an error
    #[inline]
    pub fn has_err(&self) -> bool {
        self.error.is_some()
    }

    /// Convert into the inner value.
    ///
    /// Return `None` if the instance has an error.
    #[inline]
    pub fn ok(self) -> Option<T> {
        match self.error {
            Some(_) => None,
            None => Some(self.value),
        }
    }

    /// Convert into the inner error
    ///
    /// Return `None` if the instance does not contains an error
    #[inline]
    pub fn err(self) -> Option<E> {
        self.error
    }

    /// Converts from `&ValueWithError<T, E>` to `ValueWithError<&T, &E>`.
    #[inline]
    pub fn as_ref(&self) -> ValueWithError<&T, &E> {
        ValueWithError {
            value: &self.value,
            error: self.error.as_ref(),
        }
    }

    /// Converts from `&mut ValueWithError<T, S>` to
    /// `ValueWithError<&mut T, &mut E>`.
    #[inline]
    pub fn as_mut(&mut self) -> ValueWithError<&mut T, &mut E> {
        ValueWithError {
            value: &mut self.value,
            error: self.error.as_mut(),
        }
    }

    /// Extract the inner value regardless of the error.
    ///
    /// # Safety
    ///
    /// If the instance has an error, the returned value might not satisfy the
    /// post-condition.
    #[inline]
    pub unsafe fn unwrap_unchecked(self) -> T {
        self.value
    }

    /// Extract the inner value if the instance does not have an error.
    /// Else, it returns optb.
    #[inline]
    pub fn unwrap_or(self, optb: T) -> T {
        match self.error {
            Some(_) => optb,
            None => self.value,
        }
    }

    /// Extract the inner value if the instance does not have an error.
    /// Else, it calls `op` with the error value.
    #[inline]
    pub fn unwrap_or_else<F: FnOnce(E) -> T>(self, op: F) -> T {
        match self.error {
            Some(e) => op(e),
            None => self.value,
        }
    }

    /// Extract the content of error.
    ///
    /// # Panics
    ///
    /// Panics if the instance does not have an error.
    #[inline]
    pub fn unwrap_err(self) -> E {
        self.expect_err("called `ValueWithError::unwrap_err() but it does not contains an error")
    }

    /// Extract the content of error.
    ///
    /// # Panics
    ///
    /// Panics if the instance does not have an error, with given panic message.
    #[inline]
    pub fn expect_err(self, msg: &str) -> E {
        self.error.expect(msg)
    }

    /// Converts into `std::result::Result<T, E>`.
    #[inline]
    pub fn into_result(self) -> Result<T, E> {
        match self.error {
            Some(e) => Err(e),
            None => Ok(self.value),
        }
    }
}

impl<T, E: Debug> ValueWithError<T, E> {
    /// Extract the inner value.
    ///
    /// # Panics
    ///
    /// Panics if the instance has an error.
    #[inline]
    pub fn unwrap(self) -> T {
        self.expect("called `ValueWithError::unwrap()` but it contains an error")
    }

    /// Extract the inner value
    ///
    /// # Panics
    ///
    /// Panics if the instance has an error, with given panic message.
    pub fn expect(self, msg: &str) -> T {
        if let Some(e) = self.error {
            panic!("{}: {:?}", msg, &e);
        }

        self.value
    }
}

/// Estimation result for integral.
#[derive(Clone, Debug, PartialEq)]
pub struct Solution {
    /// Estimation for the integral
    pub estimate: f64,
    /// Estimated maximum absolute error for the estimation
    pub delta: f64,
    /// What times the integrand was evaluated
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

/// Result of numerical integration
pub type IntegrationResult = ValueWithError<Solution, RuntimeError>;
