#[macro_use]
extern crate bencher;

extern crate gkquad;

use bencher::{black_box, Bencher};
use std::f64::consts::PI;
use std::f64::{INFINITY, NEG_INFINITY};

use gkquad::single::algorithm::*;
use gkquad::single::{integral, Integrator};
use gkquad::Tolerance;

fn simple(b: &mut Bencher) {
    b.iter(|| {
        let interval = black_box(0.0..1.0);
        let result = integral(|x: f64| x * x, interval).estimate().unwrap();
        assert!((result - 1.0 / 3.0).abs() <= 1.49e-8);
    });
}

fn singular_points(b: &mut Bencher) {
    let integrator = Integrator::new(|x: f64| x.recip())
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
    let integrand = |x: f64| (1.0 + x * x).recip();

    b.iter(|| {
        let interval = black_box(NEG_INFINITY..INFINITY);
        let result = integral(integrand, interval).estimate().unwrap();
        assert!((result - PI).abs() <= 1.49e-8);
    });
}

benchmark_group!(benches, simple, singular_points, infinite_interval);
benchmark_main!(benches);
