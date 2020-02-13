#[macro_use]
mod common;
use common::functions::*;

use gkquad::double::algorithm::*;
use gkquad::double::{Integrator2, Range2};
use gkquad::RuntimeError;
use gkquad::Tolerance::{self, *};

struct Expect {
    value: f64,
    delta: f64,
    error: Option<RuntimeError>,
}

fn test_algorithm<A: Algorithm<fn(f64, f64) -> f64>>(
    f: fn(f64, f64) -> f64,
    r: Range2,
    _: &[(f64, f64)],
    algorithm: A,
    tol: Tolerance,
    expect: Expect,
) {
    let mut integrator = Integrator2::new(f, algorithm).tolerance(tol);

    let result = integrator.run(r.clone());
    unsafe {
        assert_eq!(result.err(), expect.error);
        assert_rel!(result.estimate_unchecked(), expect.value, 1e-15);
        assert_rel!(result.delta_unchecked(), expect.delta, 1e-7);
    }
}

#[test]
fn qag_g1() {
    let expect = Expect {
        value: 7.500000000000000000e-01,
        delta: 8.326672684688676E-15,
        error: None,
    };
    let range = Range2::square(0., 1., 0., 1.).unwrap();
    test_algorithm(g1, range, &[], QAG::new(), Relative(1e-10), expect);
}

#[test]
fn qags_g1() {
    let expect = Expect {
        value: 7.500000000000000000e-01,
        delta: 8.326672684688676E-15,
        error: None,
    };
    let range = Range2::square(0., 1., 0., 1.).unwrap();
    test_algorithm(g1, range, &[], QAGS::new(), Relative(1e-10), expect);
}
