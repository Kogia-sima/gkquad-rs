#[macro_use]
mod common;
use common::functions::*;

use gkquad::single::algorithm::*;
use gkquad::single::Integrator;
use gkquad::RuntimeError;
use gkquad::Tolerance::{self, *};

struct Expect {
    value: f64,
    delta: f64,
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

    let result = integrator.run(a..b);
    unsafe {
        assert_eq!(result.err(), expect.error);
        assert_rel!(result.estimate_unchecked(), expect.value, 1e-15);
        assert_rel!(result.delta_unchecked(), expect.delta, 1e-7);
    }

    let result = integrator.run(b..a);

    unsafe {
        assert_eq!(result.err(), expect.error);
        assert_rel!(result.estimate_unchecked(), -expect.value, 1e-15);
        assert_rel!(result.delta_unchecked(), expect.delta, 1e-7);
    }
}

#[test]
#[ignore]
fn qag_f1_15pt() {
    let expect = Expect {
        value: 7.716049382715854665E-02,
        delta: 6.679384885865053037E-12,
        error: None,
    };
    test_algorithm(f1, 0.0, 1.0, &[], QAG::new(), Relative(1e-10), expect);
}

#[test]
fn qag_f1_21pt() {
    let expect = Expect {
        value: 7.716049382716050342E-02,
        delta: 2.2279695259216852E-15,
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
        error: Some(RuntimeError::InsufficientIteration),
    };
    test_algorithm(f2, 0.0, 1.0, &[], QAG::new(), Absolute(1e-2), expect);
}

#[test]
fn qag_f3_21pt() {
    let expect = Expect {
        value: -7.238969575482959717E-01,
        delta: 1.2850150257194227E-14,
        error: Some(RuntimeError::RoundoffError),
    };
    test_algorithm(f3, 0.3, 2.71, &[], QAG::new(), Absolute(1e-14), expect);
}

#[test]
fn qags_f1() {
    let expect = Expect {
        value: 7.716049382715789440E-02,
        delta: 2.2163949776216234E-12,
        error: None,
    };
    test_algorithm(f1, 0.0, 1.0, &[], QAGS::new(), Relative(1e-10), expect);
}

#[test]
fn qags_f2() {
    let expect = Expect {
        value: 9.999999999999434E+01,
        delta: 7.87281351222191E-11,
        error: None,
    };
    test_algorithm(f2, 0.0, 1.0, &[], QAGS::new(), Absolute(1e-10), expect);
}

#[test]
fn qags_f3() {
    let expect = Expect {
        value: -7.238969575482959717E-01,
        delta: 7.999214913217825E-11,
        error: None,
    };
    test_algorithm(f3, 0.3, 2.71, &[], QAGS::new(), Absolute(1e-10), expect);
}

#[test]
fn qags_f4() {
    let expect = Expect {
        value: -5.908755278982136588E+03,
        delta: 1.3007539610690693E-10,
        error: None,
    };
    test_algorithm(f4, 1.0, 1000.0, &[], QAGS::new(), Absolute(1e-7), expect);
}

#[test]
fn qagp_f5() {
    let expect = Expect {
        value: 2.635960527469052E+02,
        delta: 2.1000791359861604E-01,
        error: None,
    };
    test_algorithm(f5, 0., 4., &[1., 2.], QAGP::new(), Relative(1e-3), expect);
}

#[test]
fn qagp_f6() {
    let expect = Expect {
        value: -9.559338370055698E-01,
        delta: 3.530509218307998E-13,
        error: Some(RuntimeError::Divergent),
    };
    test_algorithm(f6, 0., 2., &[1.], QAGP::new(), Absolute(1e-10), expect);
}

#[test]
fn qagp_f7() {
    let expect = Expect {
        value: 0E+00,
        delta: 0E+00,
        error: None,
    };
    test_algorithm(f7, -1., 1., &[0.], QAGP::new(), Absolute(1e-10), expect);
}
