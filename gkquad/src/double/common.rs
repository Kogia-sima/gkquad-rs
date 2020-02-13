use alloc::sync::Arc;

use smallvec::SmallVec;

use crate::single::Range;
use crate::Tolerance;

/// Point in 2-dimension
pub type Point2 = (f64, f64);

/// Singular points
pub type Points2 = SmallVec<[(f64, f64); 8]>;

/// Represent the range over which the integral is estimated.
#[non_exhaustive]
#[derive(Clone)]
pub enum Range2 {
    /// Square (determined by ranges for each coordinate)
    Square { xrange: Range, yrange: Range },
    /// custom range
    Custom {
        xrange: Range,
        yrange: Arc<dyn Fn(f64) -> Range>,
    },
}

impl Range2 {
    pub fn square(x1: f64, x2: f64, y1: f64, y2: f64) -> Option<Range2> {
        let xrange = Range::new(x1, x2)?;
        let yrange = Range::new(y1, y2)?;
        if !(x1.is_finite() && x2.is_finite() && y1.is_finite() && y2.is_finite()) {
            panic!("Infinite interval in 2-dimension is not already supported.");
        }
        Some(Range2::Square { xrange, yrange })
    }

    #[inline]
    pub fn custom<F: Fn(f64) -> Range + Clone + 'static>(
        x1: f64,
        x2: f64,
        yrange: F,
    ) -> Option<Range2> {
        let xrange = Range::new(x1, x2)?;
        Some(Range2::Custom {
            xrange,
            yrange: Arc::new(yrange),
        })
    }

    #[inline]
    pub fn from_ranges(xrange: Range, yrange: Range) -> Range2 {
        Range2::Square { xrange, yrange }
    }
}

/// Integration configuration
#[non_exhaustive]
#[derive(Clone, Debug, PartialEq)]
pub struct Integration2Config {
    /// the tolerance to be satisfied
    pub tolerance: Tolerance,
    /// maximum number of subdivisions
    pub limit: usize,
    /// specify singular points
    pub points: Points2,
}

impl Default for Integration2Config {
    #[inline]
    fn default() -> Self {
        Self {
            tolerance: Tolerance::default(),
            limit: 50,
            points: Points2::new(),
        }
    }
}

/// The function that is to be integrated
pub trait Integrand2 {
    /// apply function to explanatory variable `x`
    fn apply(&mut self, x: Point2) -> f64;
}

impl<F: FnMut(f64, f64) -> f64> Integrand2 for F {
    #[inline]
    fn apply(&mut self, x: Point2) -> f64 {
        (*self)(x.0, x.1)
    }
}
