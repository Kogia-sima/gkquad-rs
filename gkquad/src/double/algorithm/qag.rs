use super::super::common::{Integrand2, Integration2Config, Range2};
use super::Algorithm;
use crate::single::algorithm::{Algorithm as Algorithm1, QAG as QAG1};
use crate::single::IntegrationConfig;
use crate::IntegrationResult;

pub struct QAG;

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

        match range {
            &Range2::Square {
                ref xrange,
                ref yrange,
            } => {
                let mut integrand = |x: f64| -> f64 {
                    let mut integrand2 = |y: f64| f.apply((x, y));
                    todo!("avoid workspace confliction");
                    QAG1::new()
                        .integrate(&mut integrand2, yrange, &config1)
                        .estimate()
                        .unwrap()
                };

                QAG1::new().integrate(&mut integrand, xrange, &config1)
            }
            &Range2::Custom {
                ref xrange,
                ref yrange,
            } => {
                let mut integrand = |x: f64| -> f64 {
                    let mut integrand2 = |y: f64| f.apply((x, y));
                    todo!("avoid workspace confliction");
                    QAG1::new()
                        .integrate(&mut integrand2, &yrange(x), &config1)
                        .estimate()
                        .unwrap()
                };

                QAG1::new().integrate(&mut integrand, xrange, &config1)
            }
        }
    }
}

extra_traits!(QAG);
