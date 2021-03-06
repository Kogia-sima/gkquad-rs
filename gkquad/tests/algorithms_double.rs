#[macro_use]
mod common;
use common::functions::*;

use gkquad::double::algorithm::*;
use gkquad::double::range::*;
use gkquad::double::Integrator2;
use gkquad::single::Range;
use gkquad::RuntimeError;
use gkquad::Tolerance::{self, *};

struct Expect {
    value: f64,
    delta: f64,
    nevals: usize,
    error: Option<RuntimeError>,
}

fn test_algorithm<A, R>(
    f: fn(f64, f64) -> f64,
    r: R,
    pts: &[(f64, f64)],
    algorithm: A,
    tol: Tolerance,
    expect: Expect,
) where
    R: Range2 + Clone,
    A: Algorithm2<fn(f64, f64) -> f64, R>,
{
    let mut integrator = Integrator2::with_algorithm(f, algorithm)
        .tolerance(tol)
        .max_evals(100000)
        .points(pts);

    let result = integrator.run(r);
    assert_eq!(result.as_ref().err(), expect.error.as_ref());

    let result = unsafe { result.unwrap_unchecked() };
    assert_rel!(result.estimate, expect.value, 1e-15);
    assert_rel!(result.delta, expect.delta, 1e-7);
    assert_eq!(result.nevals, expect.nevals);
}

#[test]
fn qag_g1() {
    let expect = Expect {
        value: 7.500000000000000000e-01,
        delta: 8.326672684688676E-15,
        nevals: 289,
        error: None,
    };
    let range = Rectangle::new(0., 1., 0., 1.).unwrap();
    test_algorithm(g1, range, &[], QAG2::new(), Relative(1e-10), expect);
}

#[test]
fn qags_g1() {
    let expect = Expect {
        value: 7.500000000000000000e-01,
        delta: 8.326672684688676E-15,
        nevals: 289,
        error: None,
    };
    let range = Rectangle::new(0., 1., 0., 1.).unwrap();
    test_algorithm(g1, range, &[], QAGS2::new(), Relative(1e-10), expect);
}

#[test]
fn qag_g2() {
    let expect = Expect {
        value: 8.591409142295225e-01,
        delta: 9.538380243751227E-15,
        nevals: 289,
        error: None,
    };
    let range = DynamicY::new(0., 1., |x| (0.0..x).into()).unwrap();
    test_algorithm(g2, range, &[], QAG2::new(), Relative(1e-10), expect);
}

#[test]
fn qags_g2() {
    let expect = Expect {
        value: 8.591409142295225e-01,
        delta: 9.538380243751227E-15,
        nevals: 289,
        error: None,
    };
    let range = DynamicY::new(0., 1., |x| (0.0..x).into()).unwrap();
    test_algorithm(g2, range, &[], QAGS2::new(), Relative(1e-10), expect);
}

#[test]
fn qags_g3() {
    let expect = Expect {
        value: 2.617993877991505e-01,
        delta: 2.906557081673861E-15,
        nevals: 34239,
        error: None,
    };
    let yrange = |x: f64| {
        let ymax = (0.25 - x * x).sqrt();
        Range::new(-ymax, ymax).unwrap()
    };
    let range = DynamicY::new(-0.5, 0.5, yrange).unwrap();
    test_algorithm(g3, range, &[], QAGS2::new(), Absolute(1e-12), expect);
}

#[test]
fn qag_g4() {
    let expect = Expect {
        value: 4.9999999999953537e-01,
        delta: 9.058237841799824E-13,
        nevals: 4889,
        error: None,
    };
    let range = Rectangle::from((0.0.., 0.0..));
    test_algorithm(g4, range, &[], QAG2::new(), Absolute(1e-8), expect);
}

#[test]
fn qags_g4() {
    let expect = Expect {
        value: 4.9999999999953537e-01,
        delta: 9.058231687644665E-13,
        nevals: 4889,
        error: None,
    };
    let range = Rectangle::from((0.0.., 0.0..));
    test_algorithm(g4, range, &[], QAGS2::new(), Absolute(1e-8), expect);
}

#[test]
fn qagp_gp1() {
    let expect = Expect {
        value: 4.418277998646525e+00,
        delta: 5.178080186851730e-13,
        nevals: 84700,
        error: None,
    };
    let range = Rectangle::new(-1.0, 1.0, -1.0, 1.0).unwrap();
    test_algorithm(
        gp1,
        range,
        &[(0.0, 0.0)],
        QAGP2::new(),
        Relative(1e-6),
        expect,
    )
}

#[test]
fn qagp_gp2() {
    let expect = Expect {
        value: 0.0,
        delta: 0.0,
        nevals: 22700,
        error: None,
    };
    let yrange = |x: f64| {
        let ymax = (1.0 - x * x).sqrt();
        Range::new(-ymax, ymax).unwrap()
    };
    let range = DynamicY::new(-1.0, 1.0, yrange).unwrap();
    test_algorithm(
        gp2,
        range,
        &[(0.0, 0.0)],
        QAGP2::new(),
        Absolute(1e-5),
        expect,
    )
}

#[test]
fn qagp_gp3() {
    let expect = Expect {
        value: 3.1415174675291286e+00,
        delta: 8.343456282745318e-02,
        nevals: 100000,
        error: Some(RuntimeError::InsufficientIteration),
    };
    let yrange = |x: f64| {
        let ymax = (1.0 - x * x).sqrt();
        Range::new(-ymax, ymax).unwrap()
    };
    let range = DynamicY::new(-1.0, 1.0, yrange).unwrap();
    test_algorithm(
        gp3,
        range,
        &[(0.0, 0.0)],
        QAGP2::new(),
        Absolute(1e-5),
        expect,
    )
}
