#![allow(non_camel_case_types)]

use std::cell::RefCell;
use std::fmt::{self, Debug};
use std::rc::Rc;

#[cfg(not(feature = "std"))]
use crate::float::Float;

use crate::error::{IntegrationResult, RuntimeError::*};
use crate::single::algorithm::Algorithm;
use crate::single::common::{Integrand, IntegrationConfig, Interval};
use crate::single::qk::qk21;
use crate::single::util::{bisect, subinterval_too_small};
use crate::single::workspace::{SubIntervalInfo, WorkSpace};

#[derive(Clone)]
pub struct QAG_FINITE {
    #[doc(hidden)]
    pub workspace: Rc<RefCell<WorkSpace>>
}

impl QAG_FINITE {
    #[inline]
    pub fn new() -> Self {
        Self {
            workspace: Rc::new(RefCell::new(WorkSpace::new()))
        }
    }
}

impl Debug for QAG_FINITE {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("QAG_FINITE")
    }
}

impl<F: Integrand> Algorithm<F> for QAG_FINITE {
    fn integrate(
        &self,
        f: &mut F,
        interval: &Interval,
        config: &IntegrationConfig,
    ) -> IntegrationResult {
        let mut ws = self.workspace.borrow_mut();
        ws.clear();
        ws.reserve(config.limit);

        let (mut roundoff_type1, mut roundoff_type2) = (0_i32, 0_i32);
        let mut error = None;

        let result0 = qk21(f, &interval);

        if result0.estimate.is_nan() {
            return IntegrationResult::new(
                result0.estimate,
                result0.delta,
                Some(NanValueEncountered)
            );
        }

        ws.push(SubIntervalInfo::new(
            interval.clone(),
            result0.estimate,
            result0.delta,
            0,
        ));

        let mut tolerance = config.tolerance.to_abs(result0.estimate.abs());
        let round_off = 50. * std::f64::EPSILON * result0.absvalue;

        if result0.delta <= round_off && result0.delta > tolerance {
            // 精度の限界によりこれ以上誤差を減らすことは不可能
            return IntegrationResult::new(result0.estimate, result0.delta, Some(RoundoffError));
        } else if result0.delta <= tolerance && result0.delta != result0.asc || result0.delta == 0.0
        {
            return IntegrationResult::new(result0.estimate, result0.delta, None);
        } else if config.limit == 1 {
            return IntegrationResult::new(
                result0.estimate,
                result0.delta,
                Some(InsufficientIteration),
            );
        }

        // sum of the integral estimates for each interval
        let mut area = result0.estimate;

        // sum of the errors for each interval
        let mut deltasum = result0.delta;

        for _ in 1..config.limit {
            // 最も誤差が大きい部分区間を取り出す
            let info = ws.get();
            let current_level = info.level + 1;

            // 区間を半分に分割
            let (il1, il2) = bisect(&info.interval);

            // 各部分区間でGauss-Kronrod積分
            let result1 = qk21(f, &il1);
            let result2 = qk21(f, &il2);

            if result1.estimate.is_nan() || result2.estimate.is_nan() {
                error = Some(NanValueEncountered);
                break;
            }

            // もとの区間での推定値を部分区間ごとの推定値の和で置き換え
            let area12 = result1.estimate + result2.estimate;
            let delta12 = result1.delta + result2.delta;
            deltasum += delta12 - info.delta;
            area += area12 - info.estimate;

            // resascの値とerrorの値は理論上一致するはず
            // => しかし丸め誤差により異なる値になる場合がある
            if result1.asc != result1.delta && result2.asc != result2.delta {
                if (info.estimate - area12).abs() <= 1e-5 * area12.abs()
                    && delta12 >= 0.99 * info.delta
                {
                    roundoff_type1 += 1;
                } else {
                    roundoff_type2 += 1;
                }
            }

            tolerance = config.tolerance.to_abs(area.abs());

            // 丸め誤差が多数発生してなおかつ収束しない場合、即座にエラー終了する
            if deltasum > tolerance {
                if roundoff_type1 >= 6 || roundoff_type2 >= 20 {
                    error = Some(RoundoffError);
                }

                // very small interval cannot achive further precision
                if subinterval_too_small(il1.begin, il1.end, il2.end) {
                    error = Some(SubintervalTooSmall);
                }
            }

            // 部分区間における積分結果を保存する
            ws.update(
                SubIntervalInfo::new(il1, result1.estimate, result1.delta, current_level),
                SubIntervalInfo::new(il2, result2.estimate, result2.delta, current_level),
            );

            if error.is_some() || deltasum <= tolerance {
                break;
            }
        }

        // 再度結果を足し合わせて正確な推定値を得る
        IntegrationResult::new(ws.sum_results(), deltasum, error)
    }

    #[doc(hidden)]
    fn get_workspace(&self) -> Option<std::cell::Ref<WorkSpace>> {
        Some(self.workspace.borrow())
    }
}
