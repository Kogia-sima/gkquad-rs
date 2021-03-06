#![allow(deprecated)]

#[macro_use]
mod common;
use common::functions::*;

use gkquad::single::algorithm::*;
use gkquad::single::{Integrator, WorkSpace};
use gkquad::RuntimeError;
use gkquad::Tolerance::{self, *};

trait AlgorithmWithWorkSpace {
    fn with_workspace(ws: &mut WorkSpace) -> Self;
}

impl AlgorithmWithWorkSpace for QAG<'_> {
    fn with_workspace(ws: &mut WorkSpace) -> Self {
        let ws = unsafe { std::mem::transmute(ws) };
        Self::with_workspace(ws)
    }
}

impl AlgorithmWithWorkSpace for QAGS<'_> {
    fn with_workspace(ws: &mut WorkSpace) -> Self {
        let ws = unsafe { std::mem::transmute(ws) };
        Self::with_workspace(ws)
    }
}

impl AlgorithmWithWorkSpace for QAGP<'_> {
    fn with_workspace(ws: &mut WorkSpace) -> Self {
        let ws = unsafe { std::mem::transmute(ws) };
        Self::with_workspace(ws)
    }
}

struct Expect<'a> {
    value: f64,
    delta: f64,
    order: &'a [usize],
    nevals: usize,
    error: Option<RuntimeError>,
}

fn test_algorithm<A: Algorithm<fn(f64) -> f64> + AlgorithmWithWorkSpace>(
    f: fn(f64) -> f64,
    a: f64,
    b: f64,
    pts: &[f64],
    tol: Tolerance,
    expect: Expect,
) {
    let mut ws = WorkSpace::with_capacity(50);
    let algorithm = A::with_workspace(&mut ws);
    let mut integrator = Integrator::with_algorithm(f, algorithm)
        .tolerance(tol)
        .points(pts);
    let result = integrator.run(a..b);
    assert_eq!(result.as_ref().err(), expect.error.as_ref());

    let result = unsafe { result.unwrap_unchecked() };
    assert_rel!(result.estimate, expect.value, 1e-15);
    assert_rel!(result.delta, expect.delta, 1e-7);
    assert_eq!(result.nevals, expect.nevals);

    if cfg!(feature = "std") && !expect.order.is_empty() {
        assert_eq!(&ws.order, &expect.order);
    }

    let result = integrator.run(b..a);
    assert_eq!(result.as_ref().err(), expect.error.as_ref());

    let result = unsafe { result.unwrap_unchecked() };
    assert_rel!(result.estimate, -expect.value, 1e-15);
    assert_rel!(result.delta, expect.delta, 1e-7);
    assert_eq!(result.nevals, expect.nevals);

    if cfg!(feature = "std") && !expect.order.is_empty() && pts.is_empty() {
        assert_eq!(&ws.order, &expect.order);
    }
}

#[test]
#[ignore]
fn qag_f1_15pt() {
    let expect = Expect {
        value: 7.716049382715854665E-02,
        delta: 6.679384885865053037E-12,
        order: &[0, 1, 2, 3, 4, 5],
        nevals: 317,
        error: None,
    };
    test_algorithm::<QAG>(f1, 0.0, 1.0, &[], Relative(1e-10), expect);
}

#[test]
fn qag_f1() {
    let expect = Expect {
        value: 7.716049382716048954E-2,
        delta: 3.385161349035466596E-15,
        order: &[0, 1, 2, 3, 4, 5, 6],
        nevals: 317,
        error: None,
    };
    test_algorithm::<QAG>(f1, 0.0, 1.0, &[], Absolute(1e-14), expect);
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
        nevals: 1967,
        error: Some(RuntimeError::InsufficientIteration),
    };
    test_algorithm::<QAG>(f2, 0.0, 1.0, &[], Absolute(1e-2), expect);
}

#[test]
fn qag_f3() {
    let expect = Expect {
        value: -7.238969575482963E-1,
        delta: 1.285829033513453162E-14,
        order: &[1, 4, 3, 2, 5, 0, 6],
        nevals: 367,
        error: Some(RuntimeError::RoundoffError),
    };
    test_algorithm::<QAG>(f3, 0.3, 2.71, &[], Absolute(1e-14), expect);
}

#[test]
fn qags_f1() {
    let expect = Expect {
        value: 7.716049382715210736E-2,
        delta: 4.281425050711728165E-12,
        order: &[0, 1, 2],
        nevals: 167,
        error: None,
    };
    test_algorithm::<QAGS>(f1, 0.0, 1.0, &[], Relative(1e-10), expect);
}

#[test]
fn qags_f2() {
    let expect = Expect {
        value: 1.000000000000036806E2,
        delta: 7.300116067199269310E-11,
        order: &[0, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1],
        nevals: 767,
        error: None,
    };
    test_algorithm::<QAGS>(f2, 0.0, 1.0, &[], Absolute(1e-10), expect);
}

#[test]
fn qags_f3() {
    let expect = Expect {
        value: -7.238969575482958607E-1,
        delta: 1.842466275093522320E-14,
        order: &[],
        nevals: 42,
        error: None,
    };
    test_algorithm::<QAGS>(f3, 0.3, 2.71, &[], Absolute(1e-10), expect);
}

#[test]
fn qags_f4() {
    let expect = Expect {
        value: -5.908755278982136588E3,
        delta: 3.845076707914035937E-9,
        order: &[0, 1, 2, 3, 4, 5, 6],
        nevals: 367,
        error: None,
    };
    test_algorithm::<QAGS>(f4, 1.0, 1000.0, &[], Absolute(1e-7), expect);
}

#[test]
fn qagp_f5() {
    let expect = Expect {
        value: 2.635888729963342E2,
        delta: 2.439296220664418646E-1,
        order: &[1, 5, 0, 3, 2, 7, 8, 9, 6, 4],
        nevals: 475,
        error: None,
    };
    test_algorithm::<QAGP>(f5, 0., 4., &[1., 2.], Relative(1e-3), expect);
}

#[test]
fn qagp_f6() {
    let expect = Expect {
        value: -9.559338370056563727E-1,
        delta: 2.436939539052218606E-13,
        order: &[0, 1, 3, 4, 6, 8, 10, 11, 9, 7, 5, 2],
        nevals: 550,
        error: Some(RuntimeError::Divergent),
    };
    test_algorithm::<QAGP>(f6, 0., 2., &[1.], Absolute(1e-10), expect);
}

#[test]
fn qagp_f7() {
    let expect = Expect {
        value: 0E+00,
        delta: 0E+00,
        order: &[0, 1, 5, 4, 3, 2],
        nevals: 250,
        error: None,
    };
    test_algorithm::<QAGP>(f7, -1., 1., &[0.], Absolute(1e-10), expect);
}
