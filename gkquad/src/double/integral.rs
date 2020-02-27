use super::algorithm::*;
use super::common::{Integrand2, IntegrationConfig2};
use super::range::IntoRange2;
use crate::error::IntegrationResult;

/// Performs integration using `QAGS2` algorithm, which achieves greate performance for
/// many kinds of functions.
#[inline]
pub fn integral2<F, R>(mut f: F, r: R) -> IntegrationResult
where
    F: Integrand2,
    R: IntoRange2,
    QAGS2: Algorithm2<F, R::IntoRange>,
{
    QAGS2::new().integrate(&mut f, &r.into_range(), &IntegrationConfig2::default())
}

/// Performs the integration with custom configuration
#[inline]
pub fn integral2_with_config<'a, F, R>(
    mut f: F,
    r: R,
    config: &IntegrationConfig2,
) -> IntegrationResult
where
    F: Integrand2,
    R: IntoRange2,
    QAGS2: Algorithm2<F, R::IntoRange>,
{
    QAGS2::new().integrate(&mut f, &r.into_range(), config)
}
