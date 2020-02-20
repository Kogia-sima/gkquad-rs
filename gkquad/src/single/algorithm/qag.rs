use std::cell::UnsafeCell;

use crate::error::{IntegrationResult, RuntimeError::*};
use crate::single::algorithm::Algorithm;
use crate::single::common::{Integrand, IntegrationConfig, Range};
use crate::single::util::{bisect, subrange_too_small, transform_range, IntegrandWrapper};
use crate::single::workspace::{SubRangeInfo, WorkSpaceId, WorkSpaceProvider};
use crate::single::{qk17, qk25, QKResult};

#[derive(Clone)]
#[deprecated(since = "0.0.3", note = "QAG algorithm is always worse than QAGS.")]
pub struct QAG {
    id: WorkSpaceId,
}

impl QAG {
    #[inline]
    pub fn new() -> Self {
        Self::with_id(WorkSpaceId::Single)
    }

    #[inline]
    pub(crate) fn with_id(id: WorkSpaceId) -> Self {
        Self { id }
    }
}

impl<F: Integrand> Algorithm<F> for QAG {
    fn integrate(&self, f: &mut F, range: &Range, config: &IntegrationConfig) -> IntegrationResult {
        let transform = !range.begin.is_finite() || !range.end.is_finite();
        let wrapper = UnsafeCell::new(IntegrandWrapper {
            inner: f,
            transform,
        });

        let qk17 = |r: &Range| unsafe { qk17(&mut *wrapper.get(), r) };
        let qk25 = |r: &Range| unsafe { qk25(&mut *wrapper.get(), r) };

        if transform {
            let range = transform_range(range);
            integrate_impl(&qk17, &qk25, &range, config, self.id)
        } else {
            integrate_impl(&qk17, &qk25, range, config, self.id)
        }
    }
}

extra_traits!(QAG);

fn integrate_impl(
    qk17: &dyn Fn(&Range) -> QKResult,
    qk25: &dyn Fn(&Range) -> QKResult,
    range: &Range,
    config: &IntegrationConfig,
    id: WorkSpaceId,
) -> IntegrationResult {
    let (mut roundoff_type1, mut roundoff_type2) = (0_i32, 0_i32);
    let mut error = None;

    // initial integral
    let (result0, finished) = initial_integral(qk17, qk25, range, config);
    if finished {
        return result0;
    }

    // sum of the integral estimates for each range
    let mut area = result0.estimate;

    // sum of the errors for each range
    let mut deltasum = result0.delta;

    let provider = WorkSpaceProvider::new(id);
    let mut ws = provider.get_mut();
    ws.clear();
    ws.reserve(config.limit);

    ws.push(SubRangeInfo::new(
        range.clone(),
        result0.estimate,
        result0.delta,
        0,
    ));

    for _ in 1..config.limit {
        // 最も誤差が大きい部分区間を取り出す
        let info = ws.get();
        let current_level = info.level + 1;

        // 区間を半分に分割
        let (r1, r2) = bisect(&info.range);

        // 各部分区間でGauss-Kronrod積分
        let result1 = qk25(&r1);
        let result2 = qk25(&r2);

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
            if (info.estimate - area12).abs() <= 1e-5 * area12.abs() && delta12 >= 0.99 * info.delta
            {
                roundoff_type1 += 1;
            } else {
                roundoff_type2 += 1;
            }
        }

        let tolerance = config.tolerance.to_abs(area.abs());

        // 丸め誤差が多数発生してなおかつ収束しない場合、即座にエラー終了する
        if deltasum > tolerance {
            if roundoff_type1 >= 6 || roundoff_type2 >= 20 {
                error = Some(RoundoffError);
            } else if subrange_too_small(r1.begin, r1.end, r2.end) {
                error = Some(SubrangeTooSmall);
            }

            if error.is_some() {
                return IntegrationResult::new(
                    ws.sum_results() - info.estimate + result1.estimate + result2.estimate,
                    deltasum,
                    error,
                );
            }
        }

        // 部分区間における積分結果を保存する
        ws.update(
            SubRangeInfo::new(r1, result1.estimate, result1.delta, current_level),
            SubRangeInfo::new(r2, result2.estimate, result2.delta, current_level),
        );

        if deltasum <= tolerance {
            break;
        }
    }

    // 再度結果を足し合わせて正確な推定値を得る
    IntegrationResult::new(ws.sum_results(), deltasum, error)
}

// initial integral
fn initial_integral(
    qk17: &dyn Fn(&Range) -> QKResult,
    qk25: &dyn Fn(&Range) -> QKResult,
    range: &Range,
    config: &IntegrationConfig,
) -> (IntegrationResult, bool) {
    for i in 0..2 {
        let result0 = if i == 0 {
            qk17(&range)
        } else if i == 1 {
            qk25(&range)
        } else {
            unreachable!();
        };

        if result0.estimate.is_nan() {
            return (
                IntegrationResult::new(result0.estimate, result0.delta, Some(NanValueEncountered)),
                true,
            );
        }

        let tolerance = config.tolerance.to_abs(result0.estimate.abs());
        if result0.delta <= tolerance && result0.delta != result0.asc || result0.delta == 0.0 {
            return (
                IntegrationResult::new(result0.estimate, result0.delta, None),
                true,
            );
        } else if config.limit == 1 {
            return (
                IntegrationResult::new(
                    result0.estimate,
                    result0.delta,
                    Some(InsufficientIteration),
                ),
                true,
            );
        }

        let round_off = 50. * std::f64::EPSILON * result0.absvalue;
        if result0.delta <= round_off && result0.delta > tolerance {
            // 精度の限界によりこれ以上誤差を減らすことは不可能
            return (
                IntegrationResult::new(result0.estimate, result0.delta, Some(RoundoffError)),
                true,
            );
        }

        if i == 1 || result0.delta > tolerance * 1024. {
            return (
                IntegrationResult::new(result0.estimate, result0.delta, None),
                false,
            );
        }
    }

    unreachable!();
}
