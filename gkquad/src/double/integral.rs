use super::algorithm::*;
use super::common::{Integrand2, IntegrationConfig2, Range2};
use crate::error::IntegrationResult;

#[inline]
pub fn integral2<'a, F: Integrand2, R: AsRef<Range2<'a>>>(mut f: F, r: R) -> IntegrationResult {
    QAGS2::new().integrate(&mut f, r.as_ref(), &IntegrationConfig2::default())
}

#[inline]
pub fn integral2_with_config<'a, F: Integrand2, R: AsRef<Range2<'a>>>(
    mut f: F,
    r: R,
    config: &IntegrationConfig2,
) -> IntegrationResult {
    QAGS2::new().integrate(&mut f, r.as_ref(), config)
}
