pub use crate::{IntegrationResult, RuntimeError, Tolerance};

#[cfg(feature = "single")]
#[cfg_attr(docsrs, doc(cfg(feature = "single")))]
pub use crate::single::{integral, Integrand, Integrator};
