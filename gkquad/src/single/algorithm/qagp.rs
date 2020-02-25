use std::cell::UnsafeCell;

use crate::error::{IntegrationResult, RuntimeError::*};
use crate::single::algorithm::Algorithm;
use crate::single::common::{Integrand, IntegrationConfig, Points, Range};
use crate::single::qelg::ExtrapolationTable;
use crate::single::qk::{qk25, QKResult};
use crate::single::util::{
    bisect, insert_sort, subrange_too_small, test_positivity, transform_point, transform_range,
    IntegrandWrapper,
};
use crate::single::workspace::{SubRangeInfo, WorkSpaceId, WorkSpaceProvider};

#[derive(Clone)]
pub struct QAGP {
    id: WorkSpaceId,
}

impl QAGP {
    #[inline]
    pub fn new() -> Self {
        Self::with_id(WorkSpaceId::Single)
    }

    #[inline]
    pub(crate) fn with_id(id: WorkSpaceId) -> Self {
        Self { id }
    }
}

impl<F: Integrand> Algorithm<F> for QAGP {
    #[inline]
    fn integrate(&self, f: &mut F, range: &Range, config: &IntegrationConfig) -> IntegrationResult {
        let transform = !range.begin.is_finite() || !range.end.is_finite();
        let wrapper = UnsafeCell::new(IntegrandWrapper {
            inner: f,
            transform,
        });

        let qk25 = |r: &Range| unsafe { qk25(&mut *wrapper.get(), r) };

        if transform {
            let range = transform_range(range);
            integrate_impl(&qk25, &range, config, self.id, transform)
        } else {
            integrate_impl(&qk25, range, config, self.id, transform)
        }
    }
}

extra_traits!(QAGP);

fn integrate_impl(
    qk25: &dyn Fn(&Range) -> QKResult,
    range: &Range,
    config: &IntegrationConfig,
    id: WorkSpaceId,
    transform: bool,
) -> IntegrationResult {
    let pts = make_sorted_points(range, &config.points, transform);
    let nint = pts.len() - 1; // number of ranges

    let provider = WorkSpaceProvider::new(id);
    let mut ws = provider.get_mut();
    ws.clear();
    ws.reserve(usize::max(config.max_iters, pts.len()));

    let (mut reseps, mut abseps, mut correc) = (0.0, 0.0, 0.0);
    let mut ktmin = 0;
    let (mut roundoff_type1, mut roundoff_type2, mut roundoff_type3) = (0, 0, 0);
    let mut error = None;
    let mut error2 = false;

    let mut extrapolate = false;
    let mut disallow_extrapolation = false;

    let mut result0 = QKResult {
        estimate: 0.,
        delta: 0.,
        absvalue: 0.,
        asc: 0.,
    };

    for w in pts.windows(2) {
        // ignore small range
        if (w[1] - w[0]).abs() < 100. * std::f64::MIN_POSITIVE {
            continue;
        }

        let range = unsafe { Range::new_unchecked(w[0], w[1]) };
        let result1 = qk25(&range);

        if result1.estimate.is_nan() {
            return IntegrationResult::new(
                result0.estimate,
                result0.delta,
                Some(NanValueEncountered),
            );
        }

        let current_level = (result1.delta == result1.asc && result1.delta != 0.0) as usize;
        add_qkresult(&mut result0, &result1);

        ws.push(SubRangeInfo::new(
            range,
            result1.estimate,
            result1.delta,
            current_level,
        ));
    }

    //# Compute the initial error estimate
    let mut deltasum = 0.;

    for si in ws.subranges.iter_mut() {
        if si.level > 0 {
            si.delta = result0.delta;
        }

        deltasum += si.delta;
    }

    ws.subranges.iter_mut().for_each(|si| si.level = 0);

    // Sort results into order of decreasing error via the indirection
    // array order[]

    ws.sort_results();

    // Test on accuracy

    let tolerance = config.tolerance.to_abs(result0.estimate.abs());

    let round_off = 100. * std::f64::EPSILON * result0.absvalue;

    if result0.delta <= round_off && result0.delta > tolerance {
        return IntegrationResult::new(result0.estimate, result0.delta, Some(RoundoffError));
    } else if result0.delta <= tolerance && result0.delta != result0.asc || result0.delta == 0.0 {
        return IntegrationResult::new(result0.estimate, result0.delta, None);
    } else if config.max_iters == 1 {
        return IntegrationResult::new(
            result0.estimate,
            result0.delta,
            Some(InsufficientIteration),
        );
    }

    // Initialization

    let mut table = ExtrapolationTable::default();
    table.append(result0.estimate);

    let mut area = result0.estimate;
    let mut res_ext = result0.estimate;
    let mut err_ext = std::f64::MAX;
    let mut error_over_large_ranges = deltasum;
    let mut ertest = tolerance;

    for iteration in nint..=config.max_iters {
        let info = ws.get();

        let current_level = info.level + 1;
        let (r1, r2) = bisect(&info.range);

        let result1 = qk25(&r1);
        let result2 = qk25(&r2);

        if result1.estimate.is_nan() || result2.estimate.is_nan() {
            error = Some(NanValueEncountered);
            break;
        }

        let area12 = result1.estimate + result2.estimate;
        let error12 = result1.delta + result2.delta;
        let last_e_i = info.delta;

        // Improve previous approximations to the integral and test for
        // accuracy.

        deltasum += error12 - info.delta;
        area += area12 - info.estimate;
        let tolerance = config.tolerance.to_abs(area.abs());

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
            error2 = true;
        }

        // set error flag in the case of bad integrand behaviour at
        // a point of the integration range

        if subrange_too_small(r1.begin, r1.end, r2.end) {
            error = Some(SubrangeTooSmall);
        }

        if deltasum <= tolerance {
            return IntegrationResult::new(
                ws.sum_results() - info.estimate + result1.estimate + result2.estimate,
                deltasum,
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

        if iteration >= config.max_iters - 1 {
            error = Some(InsufficientIteration);
            break;
        }

        if disallow_extrapolation {
            continue;
        }

        error_over_large_ranges += -last_e_i;

        if current_level < ws.maximum_level {
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

        // The smallest range has the largest error.  Before
        // bisecting decrease the sum of the errors over the larger
        // ranges (error_over_large_ranges) and perform
        // extrapolation.
        if !error2 && error_over_large_ranges > ertest && ws.increase_nrmax() {
            continue;
        }

        // Perform extrapolation

        table.append(area);
        if table.n < 3 {
            ws.reset_nrmax();
            extrapolate = false;
            error_over_large_ranges = deltasum;
            continue;
        }

        table.qelg(&mut reseps, &mut abseps);
        ktmin += 1;

        if ktmin > 5 && err_ext < 0.001 * deltasum {
            error = Some(RoundoffError);
        }

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

        ws.reset_nrmax();
        extrapolate = false;
        error_over_large_ranges = deltasum;
    }

    if err_ext == std::f64::MAX {
        return IntegrationResult::new(ws.sum_results(), deltasum, error);
    }
    if error.is_some() || error2 {
        if error2 {
            err_ext += correc;
        }

        if error.is_none() {
            error = Some(RoundoffError);
        }

        if res_ext != 0. && area != 0. {
            if err_ext / res_ext.abs() > deltasum / area.abs() {
                return IntegrationResult::new(ws.sum_results(), deltasum, error);
            }
        } else if err_ext > deltasum {
            return IntegrationResult::new(ws.sum_results(), deltasum, error);
        } else if area == 0.0 {
            return IntegrationResult::new(res_ext, err_ext, error);
        }
    }

    let positive_integrand = test_positivity(result0.estimate, result0.absvalue);
    if !positive_integrand && f64::max(res_ext.abs(), area.abs()) < 0.01 * result0.absvalue {
        return IntegrationResult::new(res_ext, err_ext, error);
    }

    let ratio = res_ext / area;
    if (ratio < 0.01 || ratio > 100. || deltasum > area.abs()) && error.is_none() {
        error = Some(Divergent);
    }

    IntegrationResult::new(res_ext, err_ext, error)
}

fn make_sorted_points(range: &Range, pts: &[f64], transform: bool) -> Points {
    let (min, max) = if range.begin < range.end {
        (range.begin, range.end)
    } else {
        (range.end, range.begin)
    };

    let mut pts2 = Points::with_capacity(pts.len() + 2);
    pts2.push(range.begin);

    if transform {
        pts2.extend(pts.iter().map(|v| transform_point(*v)));
    } else {
        pts2.extend_from_slice(pts);
    }

    if range.begin < range.end {
        insert_sort(&mut pts2[1..], &mut |a, b| a < b);
    } else {
        insert_sort(&mut pts2[1..], &mut |a, b| a > b);
    };

    pts2.retain(|&mut x| min <= x && x <= max);
    pts2.push(range.end);
    pts2
}

fn add_qkresult(result1: &mut QKResult, result2: &QKResult) {
    result1.estimate += result2.estimate;
    result1.delta += result2.delta;
    result1.absvalue += result2.absvalue;
    result1.asc += result2.asc;
}
