use alloc::borrow::Cow;
use core::cell::UnsafeCell;

use crate::common::{IntegrationResult, Solution};
use crate::error::RuntimeError::{self, *};
use crate::single::algorithm::Algorithm;
use crate::single::common::{Integrand, IntegrationConfig, Range};
use crate::single::qelg::ExtrapolationTable;
use crate::single::qk::{qk17, qk25, QKResult};
use crate::single::util::{
    bisect, subrange_too_small, test_positivity, transform_range, IntegrandWrapper,
};
use crate::single::workspace::{SubRangeInfo, WorkSpace};
use crate::utils::CowMut;

#[cfg(not(feature = "std"))]
use crate::float::Float;

#[derive(Clone)]
pub struct QAGS<'a> {
    workspace: CowMut<'a, WorkSpace>,
}

impl<'a> QAGS<'a> {
    #[inline]
    pub fn new() -> Self {
        Self {
            workspace: CowMut::Owned(WorkSpace::new()),
        }
    }

    #[inline]
    #[doc(hidden)]
    pub fn with_workspace(ws: &'a mut WorkSpace) -> Self {
        Self {
            workspace: CowMut::Borrowed(ws),
        }
    }
}

#[doc(hidden)]
impl<'a> super::AlgorithmWithWorkSpace for QAGS<'a> {
    fn workspace(&self) -> &WorkSpace {
        &*self.workspace
    }
}

impl<'a, F: Integrand> Algorithm<F> for QAGS<'a> {
    fn integrate(
        &mut self,
        f: &mut F,
        range: &Range,
        config: &IntegrationConfig,
    ) -> IntegrationResult {
        let transform = !range.begin.is_finite() || !range.end.is_finite();
        let wrapper = UnsafeCell::new(IntegrandWrapper {
            inner: f,
            transform,
        });
        let range = if transform {
            Cow::Owned(transform_range(range))
        } else {
            Cow::Borrowed(range)
        };

        let qk17 = |r: &Range| unsafe { qk17(&mut *wrapper.get(), r) };
        let qk25 = |r: &Range| unsafe { qk25(&mut *wrapper.get(), r) };
        integrate_impl(&qk17, &qk25, &range, config, &mut *self.workspace)
    }
}

extra_traits!(QAGS<'a>);

fn integrate_impl(
    qk17: &dyn Fn(&Range) -> QKResult,
    qk25: &dyn Fn(&Range) -> QKResult,
    range: &Range,
    config: &IntegrationConfig,
    ws: &mut WorkSpace,
) -> IntegrationResult {
    let mut ertest = 0f64;
    let mut error_over_large_ranges = 0f64;
    let mut correc = 0.;

    let mut ktmin = 0usize;
    let (mut roundoff_type1, mut roundoff_type2, mut roundoff_type3) = (0i32, 0i32, 0i32);
    let mut error = None;
    let mut error2 = 0i32;

    let mut extrapolate = false;
    let mut disallow_extrapolation = false;

    if config.max_evals < 17 {
        return IntegrationResult::with_error(Solution::default(), InsufficientIteration);
    }

    let (result0, absvalue, finished) = initial_integral(qk17, qk25, range, config);
    if finished {
        return result0;
    }

    let result0 = result0.unwrap();
    let mut nevals = result0.nevals;

    ws.clear();
    ws.reserve((config.max_evals - nevals) / 50 + 1);

    ws.push(SubRangeInfo::new(
        range.clone(),
        result0.estimate,
        result0.delta,
        0,
    ));

    // 計算結果を補外用の配列に加える
    let mut table = ExtrapolationTable::default();
    table.append(result0.estimate);

    let mut area = result0.estimate;
    let mut errsum = result0.delta;

    // 現在の計算結果を保存
    let mut res_ext = result0.estimate;
    let mut err_ext = core::f64::MAX;
    let max_iters = (config.max_evals - nevals) / 50;

    for iteration in 1..=max_iters {
        // Bisect the subrange with the largest error estimate
        let info = ws.get();
        let current_level = info.level + 1;

        let (r1, r2) = bisect(&info.range);

        // 各部分区間でGauss-Kronrod積分
        let result1 = qk25(&r1);
        let result2 = qk25(&r2);
        nevals += 50;

        if result1.estimate.is_nan() || result2.estimate.is_nan() {
            error = Some(NanValueEncountered);
            break;
        }

        let area12 = result1.estimate + result2.estimate;
        let error12 = result1.delta + result2.delta;
        let last_e_i = info.delta;

        // もとの区間での推定値を部分区間ごとの推定値の和で置き換え
        errsum += error12 - info.delta;
        area += area12 - info.estimate;

        let tolerance = config.tolerance.to_abs(area.abs());

        // resascの値とerrorの値は理論上一致するはず
        // => しかし丸め誤差により異なる値になる場合がある
        if result1.asc != result1.delta && result2.asc != result2.delta {
            if (info.estimate - area12).abs() <= 1e-5 * area12.abs() && error12 >= 0.99 * info.delta
            {
                if !extrapolate {
                    roundoff_type1 += 1;
                } else {
                    roundoff_type2 += 1;
                }
            }

            if iteration > 10 && error12 > info.delta {
                roundoff_type3 += 1;
            }
        }

        // Test for roundoff and eventually set error flag
        if roundoff_type1 + roundoff_type2 >= 10 || roundoff_type3 >= 20 {
            error = Some(RoundoffError);
        }

        if roundoff_type2 >= 5 {
            error2 = 1;
        }

        // set error flag in the case of bad integrand behaviour at a point of
        // the integration range
        if subrange_too_small(r1.begin, r1.end, r2.end) {
            error = Some(SubrangeTooSmall);
        }

        // 要求精度を下回った場合即座にreturnする
        if errsum <= tolerance {
            return finish(
                ws.sum_results() - info.estimate + result1.estimate + result2.estimate,
                errsum,
                nevals,
                error,
            );
        }

        // append the newly-created ranges to the list
        ws.update(
            SubRangeInfo::new(r1, result1.estimate, result1.delta, current_level),
            SubRangeInfo::new(r2, result2.estimate, result2.delta, current_level),
        );

        if error.is_some() {
            break;
        }

        // 最終ループでは補外を行う必要がないため即座にreturnする
        if iteration >= config.max_evals - 1 {
            error = Some(InsufficientIteration);
            break;
        }

        // 最初のループで補外用のパラメータを初期化する
        if iteration == 2 {
            error_over_large_ranges = errsum;
            ertest = tolerance;
            table.append(area);
            continue;
        };

        if disallow_extrapolation {
            continue;
        }

        error_over_large_ranges -= last_e_i;

        if current_level < ws.maximum_level() {
            error_over_large_ranges += error12;
        }
        if !extrapolate {
            // 次に分割する区間が最小区間である場合のみ、補外を行う
            if ws.get().level < ws.maximum_level() {
                continue;
            }

            extrapolate = true;
            ws.nrmax = 1;
        }

        // 大区間のみの誤差がまだ要求値を上回っている場合、大区間の分割を優先する
        if error_over_large_ranges > ertest && ws.increase_nrmax() {
            continue;
        }

        // 今までの計算結果から収束値を推定する
        let (mut reseps, mut abseps) = (0., 0.);
        table.append(area);
        table.qelg(&mut reseps, &mut abseps);

        ktmin += 1;
        if ktmin > 5 && err_ext < 0.001 * errsum {
            error = Some(RoundoffError);
        }

        // 補外の精度が前回（の補外）を上回った場合、結果を置き換える
        if abseps < err_ext {
            ktmin = 0;
            err_ext = abseps;
            res_ext = reseps;
            correc = error_over_large_ranges;
            ertest = config.tolerance.to_abs(reseps.abs());
            if err_ext <= ertest {
                break;
            }
        }

        // Prepare bisection of the smallest range.
        if table.n == 1 {
            disallow_extrapolation = true;
        }

        if error.is_some() {
            break;
        }

        // work on range with largest error
        ws.reset_nrmax();
        extrapolate = false;
        error_over_large_ranges = errsum;
    }

    if err_ext == core::f64::MAX {
        return finish(ws.sum_results(), errsum, nevals, error);
    }

    if error.is_some() || error2 > 0 {
        if error2 > 0 {
            err_ext += correc;
        }

        if error.is_none() {
            error = Some(RoundoffError);
        }

        if res_ext != 0.0 && area != 0.0 {
            if err_ext / res_ext.abs() > errsum / area.abs() {
                return finish(ws.sum_results(), errsum, nevals, error);
            }
        } else if err_ext > errsum {
            return finish(ws.sum_results(), errsum, nevals, error);
        } else if area == 0.0 {
            return finish(res_ext, err_ext, nevals, error);
        }
    }

    //  Test on divergence.
    let positive_integrand = test_positivity(result0.estimate, absvalue);

    if !positive_integrand && f64::max(res_ext.abs(), area.abs()) < 0.01 * absvalue {
        return finish(res_ext, err_ext, nevals, error);
    }

    let ratio = res_ext / area;
    if (ratio < 0.01 || ratio > 100.0 || errsum > area.abs()) && error.is_none() {
        error = Some(Divergent);
    }

    finish(res_ext, err_ext, nevals, error)
}

// initial integral
#[inline]
fn initial_integral(
    qk17: &dyn Fn(&Range) -> QKResult,
    qk25: &dyn Fn(&Range) -> QKResult,
    range: &Range,
    config: &IntegrationConfig,
) -> (IntegrationResult, f64, bool) {
    let mut solution = Solution::default();
    let mut absvalue = 0.0;

    for i in 0..2 {
        let result0 = if i == 0 {
            solution.nevals += 17;
            qk17(&range)
        } else if i == 1 {
            solution.nevals += 25;
            qk25(&range)
        } else {
            unreachable!();
        };

        solution.estimate = result0.estimate;
        solution.delta = result0.delta;
        absvalue = result0.absvalue;

        if result0.estimate.is_nan() {
            return (
                IntegrationResult::with_error(solution, NanValueEncountered),
                0.0,
                true,
            );
        }

        let tolerance = config.tolerance.to_abs(result0.estimate.abs());
        if result0.delta <= tolerance && result0.delta != result0.asc || result0.delta == 0.0 {
            return (IntegrationResult::new(solution), 0.0, true);
        }

        let round_off = 100. * core::f64::EPSILON * result0.absvalue;
        if result0.delta <= round_off && result0.delta > tolerance {
            // 精度の限界によりこれ以上誤差を減らすことは不可能
            return (
                IntegrationResult::with_error(solution, RoundoffError),
                0.0,
                true,
            );
        }

        if config.max_evals < 42 + i * 25 {
            return (
                IntegrationResult::with_error(solution, InsufficientIteration),
                0.0,
                true,
            );
        }

        if i == 0 && result0.delta > tolerance * 1024. {
            break;
        }
    }

    return (IntegrationResult::new(solution), absvalue, false);
}

#[inline]
#[must_use]
fn finish(
    estimate: f64,
    delta: f64,
    nevals: usize,
    error: Option<RuntimeError>,
) -> IntegrationResult {
    IntegrationResult {
        value: Solution {
            estimate,
            delta,
            nevals,
        },
        error,
    }
}
