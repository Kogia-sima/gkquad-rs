use super::super::common::{Integrand2, Integration2Config, Range2};
use super::Algorithm;
use crate::single::algorithm::{Algorithm as Algorithm1, QAG as QAG1};
use crate::single::IntegrationConfig;
use crate::single::WorkSpaceId;
use crate::IntegrationResult;

pub struct QAG;

impl QAG {
    pub fn new() -> Self {
        Self
    }
}

impl<F: Integrand2> Algorithm<F> for QAG {
    fn integrate(
        &self,
        f: &mut F,
        range: &Range2,
        config: &Integration2Config,
    ) -> IntegrationResult {
        let config1 = IntegrationConfig {
            tolerance: config.tolerance.clone(),
            limit: config.limit,
            ..Default::default()
        };

        let inner = QAG1::with_id(WorkSpaceId::Single);

        match range {
            &Range2::Square {
                ref xrange,
                ref yrange,
            } => {
                let mut integrand = |x: f64| -> f64 {
                    let mut integrand2 = |y: f64| f.apply((x, y));
                    inner
                        .integrate(&mut integrand2, yrange, &config1)
                        .estimate()
                        .unwrap()
                };

                QAG1::with_id(WorkSpaceId::Double).integrate(&mut integrand, xrange, &config1)
            }
            &Range2::Custom {
                ref xrange,
                ref yrange,
            } => {
                let mut integrand = |x: f64| -> f64 {
                    let mut integrand2 = |y: f64| f.apply((x, y));
                    inner
                        .integrate(&mut integrand2, &yrange(x), &config1)
                        .estimate()
                        .unwrap()
                };

                QAG1::with_id(WorkSpaceId::Double).integrate(&mut integrand, xrange, &config1)
            }
        }
    }
}

extra_traits!(QAG);
