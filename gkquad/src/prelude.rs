pub use crate::{IntegrationResult, RuntimeError, Tolerance};

pub use crate::single::{algorithm::*, integral, Integrand, Integrator};

#[cfg(feature = "double")]
#[cfg_attr(docsrs, doc(cfg(feature = "double")))]
pub use crate::double::{algorithm::*, integral2, Integrand2, Integrator2};
