#[macro_use]
extern crate bencher;

extern crate gkquad;

use bencher::{black_box, Bencher};
use std::f64::{INFINITY, NEG_INFINITY};
use std::f64::consts::PI;

use gkquad::single::algorithm::*;
use gkquad::single::Integrator;
use gkquad::Tolerance;

fn simple(b: &mut Bencher) {
    let integrator = Integrator::new(|x: f64| x * x).algorithm(QAG::new());
    b.iter(|| {
        let interval = black_box(0.0..1.0);
        let result = integrator.run(interval).estimate().unwrap();
        assert!((result - 1.0 / 3.0).abs() <= 1.49e-8);
    });
}

fn singular_points(b: &mut Bencher) {
    let integrator = Integrator::new(|x: f64| 1.0 / x)
        .algorithm(QAGP::new())
        .points(&[0.])
        .tolerance(Tolerance::Absolute(1.49e-8));

    b.iter(|| {
        let interval = black_box(-1.0..1.0);
        let result = integrator.run(interval).estimate().unwrap();
        assert!(result.abs() <= 1.49e-8);
    });
}

fn infinite_interval(b: &mut Bencher) {
    let integrator = Integrator::new(|x: f64| (-0.5 * x * x).exp()).algorithm(QAGS::new());

    b.iter(|| {
        let interval = black_box(NEG_INFINITY..INFINITY);
        let result = integrator.run(interval).estimate().unwrap();
        assert!((result - (2. * PI).sqrt()).abs() <= 1.49e-8);
    });
}

benchmark_group!(benches, simple, singular_points, infinite_interval);
benchmark_main!(benches);
