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
