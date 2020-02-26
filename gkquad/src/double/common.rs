use alloc::sync::Arc;
use smallvec::SmallVec;
use std::fmt::{self, Debug};
use std::ops::RangeBounds;

use crate::single::Range;
use crate::Tolerance;

/// Point in 2-dimension
pub type Point2 = (f64, f64);

/// Singular points
pub type Points2 = SmallVec<[(f64, f64); 8]>;

/// Represent the range over which the integral is estimated.
#[non_exhaustive]
#[derive(Clone)]
pub enum Range2<'a> {
    /// Square (determined by ranges for each coordinate)
    Square { xrange: Range, yrange: Range },
    /// custom range
    Custom {
        xrange: Range,
        yrange: Arc<dyn Fn(f64) -> Range + 'a>,
    },
}

impl<'a> Range2<'a> {
    pub fn square(x1: f64, x2: f64, y1: f64, y2: f64) -> Option<Range2<'static>> {
        let xrange = Range::new(x1, x2)?;
        let yrange = Range::new(y1, y2)?;
        if !(x1.is_finite() && x2.is_finite() && y1.is_finite() && y2.is_finite()) {
            panic!("Infinite interval in 2-dimension is not already supported.");
        }
        Some(Range2::Square { xrange, yrange })
    }

    #[inline]
    pub fn custom<F>(x1: f64, x2: f64, yrange: F) -> Option<Range2<'a>>
    where
        F: Fn(f64) -> Range + 'a,
    {
        let xrange = Range::new(x1, x2)?;
        Some(Range2::Custom {
            xrange,
            yrange: Arc::new(yrange),
        })
    }
}

impl<R1: RangeBounds<f64>, R2: RangeBounds<f64>> From<(R1, R2)> for Range2<'static> {
    fn from(r: (R1, R2)) -> Range2<'static> {
        Range2::Square {
            xrange: r.0.into(),
            yrange: r.1.into(),
        }
    }
}

impl<'a> Debug for Range2<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            &Range2::Square {
                ref xrange,
                ref yrange,
            } => f
                .debug_struct("Square")
                .field("xrange", xrange)
                .field("yrange", yrange)
                .finish(),
            &Range2::Custom { ref xrange, .. } => f
                .debug_struct("Custom")
                .field("xrange", xrange)
                .field("yrange", &"<Function>")
                .finish(),
        }
    }
}

impl<'a, 'b> From<&'b Range2<'a>> for Range2<'a> {
    fn from(other: &'b Range2<'a>) -> Self {
        other.clone()
    }
}

/// Integration configuration
#[non_exhaustive]
#[derive(Clone, Debug, PartialEq)]
pub struct IntegrationConfig2 {
    /// the tolerance to be satisfied
    pub tolerance: Tolerance,
    /// maximum number of subdivisions
    pub max_iters: usize,
    /// specify singular points
    pub points: Points2,
}

impl Default for IntegrationConfig2 {
    #[inline]
    fn default() -> Self {
        Self {
            tolerance: Tolerance::default(),
            max_iters: 50,
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
