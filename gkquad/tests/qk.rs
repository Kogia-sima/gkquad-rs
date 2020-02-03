#[macro_use]
mod common;
use common::functions::*;

use gkquad::single::{qk15, qk21, qk31, qk41, qk51, qk61, QKResult, Interval};

struct Expect {
    estimate: f64,
    delta: f64,
    absvalue: f64,
    asc: f64
}

fn test_qk(
    qk: fn(f: &mut fn(f64) -> f64, r: &Interval) -> QKResult,
    mut f: fn(f64) -> f64,
    a: f64,
    b: f64,
    expect: Expect
) {
    let mut interval = Interval::new(a, b).unwrap();
    let mut result = qk(&mut f, &interval);

    assert_rel!(result.estimate, expect.estimate, 1e-15);
    assert_rel!(result.delta, expect.delta, 1e-7);
    assert_rel!(result.absvalue, expect.absvalue, 1e-15);
    assert_rel!(result.asc, expect.asc, 1e-15);

    interval = Interval::new(b, a).unwrap();
    result = qk(&mut f, &interval);

    assert_rel!(result.estimate, -expect.estimate, 1e-15);
    assert_rel!(result.delta, expect.delta, 1e-7);
    assert_rel!(result.absvalue, expect.absvalue, 1e-15);
    assert_rel!(result.asc, expect.asc, 1e-15);

}

#[test]
fn qk15_f1() {
    test_qk(qk15, f1, 0.0, 1.0, Expect {
        estimate: 7.716049357767090777E-02,
        delta: 2.990224871000550874E-06,
        absvalue: 7.716049357767090777E-02,
        asc: 4.434273814139995384E-02,
    });
}

#[test]
fn qk15_f2() {
    test_qk(qk15, f2, 0.0, 1.0, Expect {
        estimate: 1.555688196612745777E+01,
        delta: 2.350164577239293706E+01,
        absvalue: 1.555688196612745777E+01,
        asc: 2.350164577239293706E+01,
    });
}

#[test]
fn qk15_f3() {
    test_qk(qk15, f3, 0.3, 2.71, Expect {
        estimate: -7.238969575483799046E-01,
        delta: 8.760080200939757174E-06,
        absvalue: 1.165564172429140788E+00,
        asc: 9.334560307787327371E-01,
    });
}

#[test]
fn qk21_f1() {
    test_qk(qk21, f1, 0.0, 1.0, Expect {
        estimate: 7.716049379303084599E-02,
        delta: 9.424302194248481445E-08,
        absvalue: 7.716049379303084599E-02,
        asc: 4.434311425038358484E-02,
    });
}

#[test]
fn qk21_f2() {
    test_qk(qk21, f2, 0.0, 1.0, Expect {
        estimate: 1.799045317938126232E+01,
        delta: 2.782360287710622515E+01,
        absvalue: 1.799045317938126232E+01,
        asc: 2.782360287710622515E+01,
    });
}

#[test]
fn qk21_f3() {
    test_qk(qk21, f3, 0.3, 2.71, Expect {
        estimate: -7.238969575482959717E-01,
        delta: 7.999214913217825E-11,
        absvalue: 1.150829032708484023E+00,
        asc: 9.297591249133687619E-01,
    });
}

#[test]
fn qk31_f1() {
    test_qk(qk31, f1, 0.0, 1.0, Expect {
        estimate: 7.716049382494900855E-02,
        delta: 1.713503193600029893E-09,
        absvalue: 7.716049382494900855E-02,
        asc: 4.427995051868838933E-02,
    });
}

#[test]
fn qk31_f2() {
    test_qk(qk31, f2, 0.0, 1.0, Expect {
        estimate: 2.081873305159121657E+01,
        delta: 3.296500137482590276E+01,
        absvalue: 2.081873305159121301E+01,
        asc: 3.296500137482590276E+01,
    });
}

#[test]
fn qk31_f3() {
    test_qk(qk31, f3, 0.3, 2.71, Expect {
        estimate: -7.238969575482959717E-01,
        delta: 1.285805464427459261E-14,
        absvalue: 1.158150602093290571E+00,
        asc: 9.277828092501518853E-01,
    });
}

#[test]
fn qk41_f1() {
    test_qk(qk41, f1, 0.0, 1.0, Expect {
        estimate: 7.716049382681375302E-02,
        delta: 9.576386660975511224E-11,
        absvalue: 7.716049382681375302E-02,
        asc: 4.421521169637691873E-02,
    });
}

#[test]
fn qk41_f2() {
    test_qk(qk41, f2, 0.0, 1.0, Expect {
        estimate: 2.288677623903126701E+01,
        delta: 3.671538820274916048E+01,
        absvalue: 2.288677623903126701E+01,
        asc: 3.671538820274916048E+01,
    });
}

#[test]
fn qk41_f3() {
    test_qk(qk41, f3, 0.3, 2.71, Expect {
        estimate: -7.238969575482959717E-01,
        delta: 1.286535726271015626E-14,
        absvalue: 1.158808363486595328E+00,
        asc: 9.264382258645686985E-01,
    });
}

#[test]
fn qk51_f1() {
    test_qk(qk51, f1, 0.0, 1.0, Expect {
        estimate: 7.716049382708510540E-02,
        delta: 1.0020804881809666E-11,
        absvalue: 7.716049382708510540E-02,
        asc: 4.416474291216854892E-02,
    });
}

#[test]
fn qk51_f2() {
    test_qk(qk51, f2, 0.0, 1.0, Expect {
        estimate: 2.449953612016972215E+01,
        delta: 3.967771249391228849E+01,
        absvalue: 2.449953612016972215E+01,
        asc: 3.967771249391228849E+01,
    });
}

#[test]
fn qk51_f3() {
    test_qk(qk51, f3, 0.3, 2.71, Expect {
        estimate: -7.238969575482961938E-01,
        delta: 1.285290995039385778E-14,
        absvalue: 1.157687209264406381E+00,
        asc: 9.264666884071264263E-01,
    });
}

#[test]
fn qk61_f1() {
    test_qk(qk61, f1, 0.0, 1.0, Expect {
        estimate: 7.716049382713800753E-02,
        delta: 1.5660589947977072E-12,
        absvalue: 7.716049382713800753E-02,
        asc: 4.419287685934316506E-02,
    });
}

#[test]
fn qk61_f2() {
    test_qk(qk61, f2, 0.0, 1.0, Expect {
        estimate: 2.583030240976628988E+01,
        delta: 4.213750493076978643E+01,
        absvalue: 2.583030240976628988E+01,
        asc: 4.213750493076978643E+01,
    });
}

#[test]
fn qk61_f3() {
    test_qk(qk61, f3, 0.3, 2.71, Expect {
        estimate: -7.238969575482959717E-01,
        delta: 1.286438572027470736E-14,
        absvalue: 1.158720854723590099E+00,
        asc: 9.270469641771273972E-01,
    });
}
