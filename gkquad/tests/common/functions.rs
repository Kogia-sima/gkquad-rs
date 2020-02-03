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
