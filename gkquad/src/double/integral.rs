use super::algorithm::*;
use super::common::{Integrand2, Integration2Config, Range2};
use crate::error::IntegrationResult;

#[inline]
pub fn integral2<F: Integrand2, R: AsRef<Range2>>(mut f: F, r: R) -> IntegrationResult {
    QAGS::new().integrate(&mut f, r.as_ref(), &Integration2Config::default())
}

#[inline]
pub fn integral2_with_config<F: Integrand2, R: AsRef<Range2>>(
    mut f: F,
    r: R,
    config: &Integration2Config,
) -> IntegrationResult {
    QAGS::new().integrate(&mut f, r.as_ref(), config)
}
