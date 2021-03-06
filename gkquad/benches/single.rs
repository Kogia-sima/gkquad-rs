extern crate gkquad;

use smbench::*;
use std::f64::consts::PI;
use std::f64::{INFINITY, NEG_INFINITY};

use gkquad::single::Integrator;
use gkquad::Tolerance;

pub fn simple(b: &mut Bencher) {
    let mut integrator = Integrator::new(|x: f64| x * x);
    b.iter(|| {
        let range = black_box(0.0..1.0);
        let result = integrator.run(range).unwrap().estimate;
        assert!((result - 1.0 / 3.0).abs() <= 1.49e-8);
    });
}

pub fn singular_points(b: &mut Bencher) {
    let mut integrator = Integrator::new(|x: f64| x.recip())
        .points(&[0.])
        .tolerance(Tolerance::Absolute(1.49e-8));

    b.iter(|| {
        let range = black_box(-1.0..1.0);
        let result = integrator.run(range).unwrap().estimate;
        assert!(result.abs() <= 1.49e-8);
    });
}

pub fn infinite_range(b: &mut Bencher) {
    let mut integrator = Integrator::new(|x: f64| (1.0 + x * x).recip());

    b.iter(|| {
        let range = black_box(NEG_INFINITY..INFINITY);
        let result = integrator.run(range).unwrap().estimate;
        assert!((result - PI).abs() <= 1.49e-8);
    });
}

smbench_group!(single, simple, singular_points, infinite_range);
