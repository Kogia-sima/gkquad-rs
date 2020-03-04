use smallvec::SmallVec;

use crate::Tolerance;

/// Point in 2-dimension
pub type Point2 = (f64, f64);

/// Singular points
pub type Points2 = SmallVec<[(f64, f64); 8]>;

/// Integration configuration
#[non_exhaustive]
#[derive(Clone, Debug, PartialEq)]
pub struct IntegrationConfig2 {
    /// the tolerance to be satisfied
    pub tolerance: Tolerance,
    /// maximum number of subdivisions
    pub max_evals: usize,
    /// specify singular points
    pub points: Points2,
}

impl Default for IntegrationConfig2 {
    #[inline]
    fn default() -> Self {
        Self {
            tolerance: Tolerance::default(),
            max_evals: 2000,
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
