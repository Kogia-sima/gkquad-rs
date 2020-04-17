extern crate gkquad;

use smbench::*;
use std::f64::consts::PI;

use gkquad::double::Integrator2;
use gkquad::Tolerance;

fn simple(b: &mut Bencher) {
    let mut integrator = Integrator2::new(|x: f64, y: f64| x * y);
    b.iter(|| {
        let range = black_box((0.0..1.0, 0.0..1.0));
        let result = integrator.run(range).unwrap().estimate;
        assert!((result - 1.0 / 4.0).abs() <= 1.49e-8);
    });
}

fn singular_points(b: &mut Bencher) {
    let mut integrator = Integrator2::new(|x: f64, y: f64| (x * y).recip())
        .points(&[(0., 0.)])
        .tolerance(Tolerance::Absolute(1.49e-8));

    b.iter(|| {
        let range = black_box((-1.0..1.0, -1.0..1.0));
        let result = integrator.run(range).unwrap().estimate;
        assert!(result.abs() <= 1.49e-8);
    });
}

fn infinite_range(b: &mut Bencher) {
    let mut integrator = Integrator2::new(|x: f64, y: f64| {
        let denom = 1.0 + x * x + y * y;
        (denom * denom.sqrt()).recip()
    });

    b.iter(|| {
        let range = black_box((.., ..));
        let result = integrator.run(range).unwrap().estimate;
        assert!((result - 2.0 * PI).abs() <= 1.49e-8);
    });
}

smbench_trace_memory!();
smbench_group!(benches, simple, singular_points, infinite_range);
smbench_main!(benches);
