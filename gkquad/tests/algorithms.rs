#[macro_use]
mod common;
use common::functions::*;

use gkquad::single::algorithm::*;
use gkquad::single::Integrator;
use gkquad::single::WorkSpaceProvider;
use gkquad::RuntimeError;
use gkquad::Tolerance::{self, *};

struct Expect<'a> {
    value: f64,
    delta: f64,
    order: &'a [usize],
    error: Option<RuntimeError>,
}

fn test_algorithm<A: Algorithm<fn(f64) -> f64>>(
    f: fn(f64) -> f64,
    a: f64,
    b: f64,
    pts: &[f64],
    algorithm: A,
    tol: Tolerance,
    expect: Expect,
) {
    let integrator = Integrator::new(f)
        .algorithm(algorithm)
        .tolerance(tol)
        .points(pts);
    let provider = WorkSpaceProvider::new();

    let result = integrator.run(a..b);
    unsafe {
        assert_eq!(result.err(), expect.error);
        assert_rel!(result.estimate_unchecked(), expect.value, 1e-15);
        assert_rel!(result.delta_unchecked(), expect.delta, 1e-7);

        if cfg!(feature = "std") && !expect.order.is_empty() {
            let ws = provider.get_mut();
            assert_eq!(ws.size(), expect.order.len());
            assert_eq!(&ws.order, &expect.order);
        }
    }

    let result = integrator.run(b..a);

    unsafe {
        assert_eq!(result.err(), expect.error);
        assert_rel!(result.estimate_unchecked(), -expect.value, 1e-15);
        assert_rel!(result.delta_unchecked(), expect.delta, 1e-7);

        if cfg!(feature = "std") && !expect.order.is_empty() && pts.is_empty() {
            let ws = provider.get_mut();
            assert_eq!(ws.size(), expect.order.len());
            assert_eq!(&ws.order, &expect.order);
        }
    }
}

#[test]
#[ignore]
fn qag_f1_15pt() {
    let expect = Expect {
        value: 7.716049382715854665E-02,
        delta: 6.679384885865053037E-12,
        order: &[0, 1, 2, 3, 4, 5],
        error: None,
    };
    test_algorithm(f1, 0.0, 1.0, &[], QAG::new(), Relative(1e-10), expect);
}

#[test]
fn qag_f1_21pt() {
    let expect = Expect {
        value: 7.716049382716050342E-02,
        delta: 2.2279695259216852E-15,
        order: &[0, 1, 2, 3, 4, 5, 6, 7],
        error: None,
    };
    test_algorithm(f1, 0.0, 1.0, &[], QAG::new(), Absolute(1e-14), expect);
}

#[test]
#[ignore]
#[rustfmt::skip]
fn qag_f2_15pt() {
    let expect = Expect {
        value: 9.955884613859257E+01,
        delta: 3.986587136565759E-1,
        order: &[
            0, 7, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90, 89, 88, 87, 86, 85, 84, 83, 82, 81, 80, 79, 78, 77, 76, 75, 74, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 50, 51, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
        error: Some(RuntimeError::InsufficientIteration),
    };
    test_algorithm(f2, 0.0, 1.0, &[], QAG::new(), Absolute(1e-2), expect);
}

#[test]
fn qag_f3_21pt() {
    let expect = Expect {
        value: -9.449347361165283E-01,
        delta: 1.2850150257194227E-14,
        order: &[1, 4, 3, 2, 5, 0, 6],
        error: Some(RuntimeError::RoundoffError),
    };
    test_algorithm(f3, 0.3, 2.71, &[], QAG::new(), Absolute(1e-14), expect);
}

#[test]
fn qags_f1() {
    let expect = Expect {
        value: 7.716049382715789440E-02,
        delta: 2.2163949776216234E-12,
        order: &[0, 1, 2, 3, 4],
        error: None,
    };
    test_algorithm(f1, 0.0, 1.0, &[], QAGS::new(), Relative(1e-10), expect);
}

#[test]
fn qags_f2() {
    let expect = Expect {
        value: 9.99999999999887E+01,
        delta: 8.43698444441543E-11,
        order: &[0, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1],
        error: None,
    };
    test_algorithm(f2, 0.0, 1.0, &[], QAGS::new(), Absolute(1e-10), expect);
}

#[test]
fn qags_f3() {
    let expect = Expect {
        value: -7.238969575482959717E-01,
        delta: 1.2854641464232937E-14,
        order: &[],
        error: None,
    };
    test_algorithm(f3, 0.3, 2.71, &[], QAGS::new(), Absolute(1e-10), expect);
}

#[test]
fn qags_f4() {
    let expect = Expect {
        value: -5.908755278982136588E+03,
        delta: 1.3007539610690693E-10,
        order: &[0, 1, 2, 3, 4, 5, 6, 7, 8],
        error: None,
    };
    test_algorithm(f4, 1.0, 1000.0, &[], QAGS::new(), Absolute(1e-7), expect);
}

#[test]
fn qagp_f5() {
    let expect = Expect {
        value: 2.63346535755973E+02,
        delta: 2.1000791359861604E-01,
        order: &[0, 5, 1, 3, 2, 10, 7, 9, 8, 6, 4],
        error: None,
    };
    test_algorithm(f5, 0., 4., &[1., 2.], QAGP::new(), Relative(1e-3), expect);
}

#[test]
fn qagp_f6() {
    let expect = Expect {
        value: -9.559338370055698E-01,
        delta: 3.530509218307998E-13,
        order: &[0, 1, 3, 4, 6, 8, 10, 11, 9, 7, 5, 2],
        error: Some(RuntimeError::Divergent),
    };
    test_algorithm(f6, 0., 2., &[1.], QAGP::new(), Absolute(1e-10), expect);
}

#[test]
fn qagp_f7() {
    let expect = Expect {
        value: 0E+00,
        delta: 0E+00,
        order: &[0, 1, 5, 4, 3, 2],
        error: None,
    };
    test_algorithm(f7, -1., 1., &[0.], QAGP::new(), Absolute(1e-10), expect);
}
