use super::super::common::{Integrand2, IntegrationConfig2};
use super::super::range::{DynamicX, DynamicY, Rectangle};
use super::{Algorithm2, QAGP2, QAGS2};
use crate::common::IntegrationResult;

#[derive(Clone)]
pub struct AUTO2;

impl AUTO2 {
    pub fn new() -> Self {
        Self
    }
}

macro_rules! impl_algorithm2 {
    ($range:ty) => {
        impl<'a, F: Integrand2 + ?Sized> Algorithm2<F, $range> for AUTO2 {
            fn integrate(
                &mut self,
                f: &mut F,
                range: &$range,
                config: &IntegrationConfig2,
            ) -> IntegrationResult {
                if config.points.is_empty() {
                    QAGS2::new().integrate(f, range, config)
                } else {
                    QAGP2::new().integrate(f, range, config)
                }
            }
        }
    };
}

impl_algorithm2!(Rectangle);
impl_algorithm2!(DynamicX<'a>);
impl_algorithm2!(DynamicY<'a>);
