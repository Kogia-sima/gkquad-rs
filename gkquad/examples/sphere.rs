extern crate gkquad;

use gkquad::prelude::{integral, integral2, DynamicY};

fn main() {
    println!("1D: {}", 2.0 * 1.0);

    println!(
        "2D: {}",
        2.0 * integral(|x: f64| (1.0 - x * x).sqrt(), -1.0..1.0)
            .unwrap()
            .estimate
    );

    let range = DynamicY::new(-1.0, 1.0, |x| {
        let ymax = (1.0 - x * x).sqrt();
        (-ymax..ymax).into()
    })
    .unwrap();

    println!(
        "3D: {}",
        2.0 * integral2(|x: f64, y: f64| (1.0 - x * x - y * y).sqrt(), range)
            .unwrap()
            .estimate
    );
}
