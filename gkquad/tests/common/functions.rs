#![allow(dead_code)]

pub fn f1(x: f64) -> f64 {
    x.powf(2.6) * (1. / x).ln()
}

// singular points on x = 0
pub fn f2(x: f64) -> f64 {
    x.powf(-0.9) * (1. / x).ln()
}

// oscillation
pub fn f3(x: f64) -> f64 {
    (2.0f64.powf(1.3) * x.sin()).cos()
}

pub fn f4(x: f64) -> f64 {
    (1. / x).ln()
}

// singular points in [1.0, 2.0]
pub fn f5(x: f64) -> f64 {
    let x2 = x * x;
    x2 * x * (((x2 - 1.) * (x2 - 2.)).abs()).ln()
}

pub fn f6(x: f64) -> f64 {
    let x2 = x * x;
    let x3 = x2 * x;
    1.0 / (x3 - x2 + x - 1.)
}

pub fn f7(x: f64) -> f64 {
    1. / x
}

pub fn g1(x: f64, y: f64) -> f64 {
    1.0 - x * y
}

pub fn g2(x: f64, y: f64) -> f64 {
    (y / x).exp()
}

pub fn g3(x: f64, y: f64) -> f64 {
    (0.25 - x * x - y * y).sqrt()
}

pub fn g4(x: f64, y: f64) -> f64 {
    1.0 / (1.0 + x + y).powi(3)
}

// singular points in (0.0, 0.0)
pub fn gp1(x: f64, y: f64) -> f64 {
    1.0 / f64::sqrt(f64::abs(x) + f64::abs(y))
}

pub fn gp2(x: f64, y: f64) -> f64 {
    x / (x * x + y * y)
}

pub fn gp3(x: f64, y: f64) -> f64 {
    0.5 / f64::sqrt(x * x + y * y)
}
