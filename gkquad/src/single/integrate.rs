use super::algorithm::*;
use super::common::{Integrand, IntegrationConfig, Interval};

use crate::IntegrationResult;

/// Performs integration using `QAGS` algorithm,
/// which achieves great performance for many kinds of functions.
///
/// # Examples
///
/// ```
/// let result = integrate(|x| x.sqrt(), 1.0..2.0).estimate();
/// ```
pub fn integrate<F: Integrand, I: Into<Interval>>(mut f: F, interval: I) -> IntegrationResult {
    QAGS::new().integrate(&mut f, &interval.into(), &IntegrationConfig::default())
}

/// Performs the integration with custom configuration.
///
/// The algorithm will be automatically selected to achieve the greatest performance.
pub fn integrate_with_config<F: Integrand, I: Into<Interval>>(
    mut f: F,
    interval: I,
    config: IntegrationConfig,
) -> IntegrationResult {
    AUTO::new().integrate(&mut f, &interval.into(), &config)
}
