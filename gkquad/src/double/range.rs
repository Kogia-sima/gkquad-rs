//! 2-dimentional range types

use alloc::sync::Arc;
use core::fmt::{self, Debug};
use core::ops::RangeBounds;

use crate::single::Range;

/// Rectangle range
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Rectangle {
    pub xrange: Range,
    pub yrange: Range,
}

impl Rectangle {
    pub fn new(x1: f64, x2: f64, y1: f64, y2: f64) -> Option<Rectangle> {
        let xrange = Range::new(x1, x2)?;
        let yrange = Range::new(y1, y2)?;
        if !(x1.is_finite() && x2.is_finite() && y1.is_finite() && y2.is_finite()) {
            panic!("Infinite interval in 2-dimension is not already supported.");
        }
        Some(Rectangle { xrange, yrange })
    }
}

impl<R1: RangeBounds<f64>, R2: RangeBounds<f64>> From<(R1, R2)> for Rectangle {
    fn from(r: (R1, R2)) -> Rectangle {
        Rectangle {
            xrange: r.0.into(),
            yrange: r.1.into(),
        }
    }
}

/// Range of `x` varies based on `y` value
#[derive(Clone)]
pub struct DynamicX<'a> {
    pub xrange: Arc<dyn Fn(f64) -> Range + Send + Sync + 'a>,
    pub yrange: Range,
}

impl<'a> DynamicX<'a> {
    pub fn new<F>(xrange: F, y1: f64, y2: f64) -> Option<DynamicX<'a>>
    where
        F: Fn(f64) -> Range + Send + Sync + 'a,
    {
        let yrange = Range::new(y1, y2)?;
        Some(DynamicX {
            xrange: Arc::new(xrange),
            yrange,
        })
    }
}

impl<'a> Debug for DynamicX<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.debug_struct("DynamicX")
            .field("xrange", &"<Function>")
            .field("yrange", &self.yrange)
            .finish()
    }
}

impl From<Rectangle> for DynamicX<'static> {
    fn from(square: Rectangle) -> DynamicX<'static> {
        let Rectangle { xrange, yrange } = square;
        Self {
            xrange: Arc::new(move |_| xrange.clone()),
            yrange,
        }
    }
}

/// Range of `y` varies based on `x` value
#[derive(Clone)]
pub struct DynamicY<'a> {
    pub xrange: Range,
    pub yrange: Arc<dyn Fn(f64) -> Range + Send + Sync + 'a>,
}

impl<'a> DynamicY<'a> {
    pub fn new<F>(x1: f64, x2: f64, yrange: F) -> Option<DynamicY<'a>>
    where
        F: Fn(f64) -> Range + Send + Sync + 'a,
    {
        let xrange = Range::new(x1, x2)?;
        Some(DynamicY {
            xrange,
            yrange: Arc::new(yrange),
        })
    }
}

impl<'a> Debug for DynamicY<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.debug_struct("DynamicY")
            .field("xrange", &self.xrange)
            .field("yrange", &"<Function>")
            .finish()
    }
}

impl From<Rectangle> for DynamicY<'static> {
    fn from(square: Rectangle) -> DynamicY<'static> {
        let Rectangle { xrange, yrange } = square;
        Self {
            xrange,
            yrange: Arc::new(move |_| yrange.clone()),
        }
    }
}

/// 2-dimentional range type API
///
/// This is a marker trait, and does not implement anything. If you generalize
/// function with `Range2`, use [IntoRange2](./trait.IntoRange2.html) trait instead.
pub trait Range2 {}

impl Range2 for Rectangle {}
impl<'a> Range2 for DynamicX<'a> {}
impl<'a> Range2 for DynamicY<'a> {}

/// Conversion into `Range2`
pub trait IntoRange2 {
    type IntoRange: Range2;

    fn into_range(self) -> Self::IntoRange;
}

impl<T: Range2> IntoRange2 for T {
    type IntoRange = T;

    #[inline]
    fn into_range(self) -> T {
        self
    }
}

impl<'a> IntoRange2 for &'a Rectangle {
    type IntoRange = Rectangle;

    #[inline]
    fn into_range(self) -> Rectangle {
        self.clone()
    }
}

impl<'a, 'b> IntoRange2 for &'a DynamicX<'b> {
    type IntoRange = DynamicX<'b>;

    #[inline]
    fn into_range(self) -> DynamicX<'b> {
        self.clone()
    }
}

impl<'a, 'b> IntoRange2 for &'a DynamicY<'b> {
    type IntoRange = DynamicY<'b>;

    #[inline]
    fn into_range(self) -> DynamicY<'b> {
        self.clone()
    }
}

impl<R1: RangeBounds<f64>, R2: RangeBounds<f64>> IntoRange2 for (R1, R2) {
    type IntoRange = Rectangle;

    fn into_range(self) -> Rectangle {
        Rectangle::from(self)
    }
}
