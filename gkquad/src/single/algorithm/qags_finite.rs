#![allow(non_camel_case_types, clippy::float_cmp)]

use crate::error::{IntegrationResult, RuntimeError::*};
use crate::single::algorithm::Algorithm;
use crate::single::common::{Integrand, IntegrationConfig, Interval};
use crate::single::qelg::ExtrapolationTable;
use crate::single::qk::{qk15, qk21};
use crate::single::util::{bisect, subinterval_too_small, test_positivity};
use crate::single::workspace::{SubIntervalInfo, WorkSpaceProvider};

/// QAGS algorithm over finite interval
#[derive(Clone)]
pub struct QAGS_FINITE {
    provider: WorkSpaceProvider,
}

impl QAGS_FINITE {
    #[inline]
    pub fn new() -> Self {
        Self {
            provider: WorkSpaceProvider::new(),
        }
    }

    // initial integral
    fn initial_integral<F: Integrand>(
        &self,
        f: &mut F,
        interval: &Interval,
        config: &IntegrationConfig,
    ) -> (IntegrationResult, f64, bool) {
        for i in 0..2 {
            let result0 = if i == 0 {
                qk15(f, &interval)
            } else if i == 1 {
                qk21(f, &interval)
            } else {
                unreachable!();
            };

            if result0.estimate.is_nan() {
                return (
                    IntegrationResult::new(
                        result0.estimate,
                        result0.delta,
                        Some(NanValueEncountered),
                    ),
                    0.,
                    true,
                );
            }

            let tolerance = config.tolerance.to_abs(result0.estimate.abs());
            if result0.delta <= tolerance && result0.delta != result0.asc || result0.delta == 0.0 {
                return (
                    IntegrationResult::new(result0.estimate, result0.delta, None),
                    0.,
                    true,
                );
            } else if config.limit == 1 {
                return (
                    IntegrationResult::new(
                        result0.estimate,
                        result0.delta,
                        Some(InsufficientIteration),
                    ),
                    0.,
                    true,
                );
            }

            let round_off = 100. * std::f64::EPSILON * result0.absvalue;
            if result0.delta <= round_off && result0.delta > tolerance {
                // 精度の限界によりこれ以上誤差を減らすことは不可能
                return (
                    IntegrationResult::new(result0.estimate, result0.delta, Some(RoundoffError)),
                    0.,
                    true,
                );
            }

            if i == 1 || result0.delta > tolerance * 1024. {
                return (
                    IntegrationResult::new(result0.estimate, result0.delta, None),
                    result0.absvalue,
                    false,
                );
            }
        }

        unreachable!();
    }
}

impl<F: Integrand> Algorithm<F> for QAGS_FINITE {
    fn integrate(
        &self,
        f: &mut F,
        interval: &Interval,
        config: &IntegrationConfig,
    ) -> IntegrationResult {
        let mut ertest = 0f64;
        let mut error_over_large_intervals = 0f64;
        let mut correc = 0.;

        let mut ktmin = 0usize;
        let (mut roundoff_type1, mut roundoff_type2, mut roundoff_type3) = (0i32, 0i32, 0i32);
        let mut error = None;
        let mut error2 = 0i32;

        let mut extrapolate = false;
        let mut disallow_extrapolation = false;

        let (result0, absvalue, finished) = self.initial_integral(f, interval, config);
        if finished {
            return result0;
        }

        let mut ws = unsafe { self.provider.get_mut() };
        ws.clear();
        ws.reserve(config.limit);

        ws.push(SubIntervalInfo::new(
            interval.clone(),
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

        for iteration in 2..=config.limit {
            // Bisect the subinterval with the largest error estimate
            let info = ws.get();
            let current_level = info.level + 1;

            let (il1, il2) = bisect(&info.interval);

            // 各部分区間でGauss-Kronrod積分
            let result1 = qk21(f, &il1);
            let result2 = qk21(f, &il2);

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
                if (info.estimate - area12).abs() <= 1e-5 * area12.abs()
                    && error12 >= 0.99 * info.delta
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
            if subinterval_too_small(il1.begin, il1.end, il2.end) {
                error = Some(SubintervalTooSmall);
            }

            // append the newly-created intervals to the list
            ws.update(
                SubIntervalInfo::new(il1, result1.estimate, result1.delta, current_level),
                SubIntervalInfo::new(il2, result2.estimate, result2.delta, current_level),
            );

            // 要求精度を下回った場合即座にreturnする
            if errsum <= tolerance {
                return IntegrationResult::new(ws.sum_results(), errsum, error);
            }

            if error.is_some() {
                break;
            }

            // 最終ループでは補外を行う必要がないため即座にreturnする
            if iteration >= config.limit - 1 {
                error = Some(InsufficientIteration);
                break;
            }

            // 最初のループで補外用のパラメータを初期化する
            if iteration == 2 {
                error_over_large_intervals = errsum;
                ertest = tolerance;
                table.append(area);
                continue;
            };

            if disallow_extrapolation {
                continue;
            }

            error_over_large_intervals -= last_e_i;

            if current_level < ws.maximum_level() {
                error_over_large_intervals += error12;
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
            if error_over_large_intervals > ertest && ws.increase_nrmax() {
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
                correc = error_over_large_intervals;
                ertest = config.tolerance.to_abs(reseps.abs());
                if err_ext <= ertest {
                    break;
                }
            }

            // Prepare bisection of the smallest interval.
            if table.n == 1 {
                disallow_extrapolation = true;
            }

            if error.is_some() {
                break;
            }

            // work on interval with largest error
            ws.reset_nrmax();
            extrapolate = false;
            error_over_large_intervals = errsum;
        }

        if err_ext == core::f64::MAX {
            return IntegrationResult::new(ws.sum_results(), errsum, error);
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
                    return IntegrationResult::new(ws.sum_results(), errsum, error);
                }
            } else if err_ext > errsum {
                return IntegrationResult::new(ws.sum_results(), errsum, error);
            } else if area == 0.0 {
                return IntegrationResult::new(res_ext, err_ext, error);
            }
        }

        //  Test on divergence.
        let positive_integrand = test_positivity(result0.estimate, absvalue);

        if !positive_integrand && f64::max(res_ext.abs(), area.abs()) < 0.01 * absvalue {
            return IntegrationResult::new(res_ext, err_ext, error);
        }

        let ratio = res_ext / area;
        if (ratio < 0.01 || ratio > 100.0 || errsum > area.abs()) && error.is_none() {
            error = Some(Divergent);
        }

        IntegrationResult::new(res_ext, err_ext, error)
    }
}

extra_traits!(QAGS_FINITE);
