use smallvec::SmallVec;

use crate::single::Interval;
use crate::Tolerance;

/// Point in 2-dimension
pub type Point2 = (f64, f64);

/// Singular points
pub type Points2 = SmallVec<[(f64, f64); 8]>;

pub(crate) type PhantomFn = fn(f64) -> Interval;

pub(crate) struct Square {
    xrange: Interval,
    yrange: Interval,
}

/// Represent the range over which the integral is estimated.
#[non_exhaustive]
pub enum Range2 {
    /// Square (determined by ranges for each coordinate)
    Square { xrange: Interval, yrange: Interval },
    /// custom range
    Custom {
        xrange: Interval,
        yrange: Box<dyn Fn(f64) -> f64>,
    },
}

impl Range2 {
    pub fn square(x1: f64, x2: f64, y1: f64, y2: f64) -> Option<Range2> {
        if !(x1.is_finite() && x2.is_finite() && y1.is_finite() && y2.is_finite()) {
            panic!("Infinite interval in 2-dimension is not already supported.");
        }
        let xrange = Interval::new(x1, x2)?;
        let yrange = Interval::new(y1, y2)?;
        Some(Range2::Square { xrange, yrange })
    }

    #[inline]
    pub fn custom<F: Fn(f64) -> f64 + 'static>(x1: f64, x2: f64, yrange: F) -> Option<Range2> {
        let xrange = Interval::new(x1, x2)?;
        Some(Range2::Custom {
            xrange,
            yrange: Box::new(yrange),
        })
    }

    #[inline]
    pub fn from_ranges(xrange: Interval, yrange: Interval) -> Range2 {
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

    /// apply function to each of the elements of `s`.
    #[inline]
    fn apply_to_slice(&mut self, s: &mut [Point2]) {
        s.iter_mut().for_each(|x| x.0 = self.apply(*x));
    }
}

impl<F: FnMut(Point2) -> f64> Integrand2 for F {
    #[inline]
    fn apply(&mut self, x: Point2) -> f64 {
        (*self)(x)
    }
}
