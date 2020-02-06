#[macro_use]
extern crate bencher;

extern crate gkquad;

use bencher::{black_box, Bencher};
use std::f64::{INFINITY, NEG_INFINITY};
use std::f64::consts::PI;

use gkquad::single::algorithm::*;
use gkquad::single::Integrator;

fn square_qag(b: &mut Bencher) {
    let mut integrator = Integrator::new(|x: f64| x * x).algorithm(QAG::new());
    b.iter(|| {
        let interval = black_box(0.0..1.0);
        let result = integrator.run(interval).estimate().unwrap();
        assert!((result - 1.0 / 3.0).abs() <= 1.49e-8);
    });
}

fn square_qags(b: &mut Bencher) {
    let mut integrator = Integrator::new(|x: f64| x * x).algorithm(QAGS::new());
    b.iter(|| {
        let interval = black_box(0.0..1.0);
        let result = integrator.run(interval).estimate().unwrap();
        assert!((result - 1.0 / 3.0).abs() <= 1.49e-8);
    });
}

fn normal_distribution_qags(b: &mut Bencher) {
    let mut integrator = Integrator::new(|x: f64| (-0.5 * x * x).exp()).algorithm(QAGS::new());

    b.iter(|| {
        let interval = black_box(NEG_INFINITY..INFINITY);
        let result = integrator.run(interval).estimate().unwrap() / (2. * PI).sqrt();
        assert!((result - 1.0).abs() <= 1.49e-8);
    });
}

benchmark_group!(benches, square_qag, square_qags, normal_distribution_qags);
benchmark_main!(benches);
