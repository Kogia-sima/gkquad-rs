extern crate gkquad;

use gkquad::single::Integrator;
use gkquad::Tolerance;
use std::cell::Cell;
use std::f64::consts::PI;
use std::time::Instant;

fn main() {
    let timer = Instant::now();
    let m = Cell::new(1);
    let mut integrator =
        Integrator::new(|x: f64| (x * PI).sin().powi(m.get())).tolerance(Tolerance::Absolute(1.));
    let expected = 36722.362174233774391;

    for _ in 0..100000 {
        m.set(1);
        let mut product = 1.0;

        for _ in 0..10 {
            let estimate = integrator.run(0.0..1.0).unwrap().estimate;
            product *= estimate;
            m.set(m.get() + 1);
        }

        assert!((1. / product - expected).abs() < 1.);
    }

    println!("Elapsed time: {:?}", timer.elapsed() / 100000);
}
