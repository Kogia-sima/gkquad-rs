#[macro_use]
mod common;
use common::functions::*;

use gkquad::single::{qk17, qk25, qk33, qk41, qk49, qk57, QKResult, Range};

struct Expect {
    estimate: f64,
    delta: f64,
    absvalue: f64,
    asc: f64,
}

fn test_qk(
    qk: fn(f: &mut fn(f64) -> f64, r: &Range) -> QKResult,
    mut f: fn(f64) -> f64,
    a: f64,
    b: f64,
    expect: Expect,
) {
    let mut range = Range::new(a, b).unwrap();
    let mut result = qk(&mut f, &range);

    assert_rel!(result.estimate, expect.estimate, 1e-15);
    assert_rel!(result.delta, expect.delta, 1e-7);
    assert_rel!(result.absvalue, expect.absvalue, 1e-15);
    assert_rel!(result.asc, expect.asc, 1e-15);

    range = Range::new(b, a).unwrap();
    result = qk(&mut f, &range);

    assert_rel!(result.estimate, -expect.estimate, 1e-15);
    assert_rel!(result.delta, expect.delta, 1e-7);
    assert_rel!(result.absvalue, expect.absvalue, 1e-15);
    assert_rel!(result.asc, expect.asc, 1e-15);
}

#[test]
fn qk17_f1() {
    test_qk(
        qk17,
        f1,
        0.0,
        1.0,
        Expect {
            estimate: 7.716049368182086032E-02,
            delta: 8.277426086317404754E-07,
            absvalue: 7.716049368182086032E-02,
            asc: 4.421466357285612492E-02,
        },
    );
}

#[test]
fn qk17_f2() {
    test_qk(
        qk17,
        f2,
        0.0,
        1.0,
        Expect {
            estimate: 1.647598925451970331E+01,
            delta: 2.506857790293670263E+01,
            absvalue: 1.647598925451970331E+01,
            asc: 2.506857790293670263E+01,
        },
    );
}

#[test]
fn qk17_f3() {
    test_qk(
        qk17,
        f3,
        0.3,
        2.71,
        Expect {
            estimate: -7.238969575482938623E-01,
            delta: 6.992691324590760438E-08,
            absvalue: 1.148366070353963986E+00,
            asc: 9.201990573819788244E-01,
        },
    );
}

#[test]
fn qk25_f1() {
    test_qk(
        qk25,
        f1,
        0.0,
        1.0,
        Expect {
            estimate: 7.716049381682653363E-02,
            delta: 1.573492378686401831E-08,
            absvalue: 7.716049381682653363E-02,
            asc: 4.408419558833566454E-02,
        },
    );
}

#[test]
fn qk25_f2() {
    test_qk(
        qk25,
        f2,
        0.0,
        1.0,
        Expect {
            estimate: 1.925410329897458439E+01,
            delta: 3.013694371714961306E+01,
            absvalue: 1.925410329897458439E+01,
            asc: 3.013694371714961306E+01,
        },
    );
}

#[test]
fn qk25_f3() {
    test_qk(
        qk25,
        f3,
        0.3,
        2.71,
        Expect {
            estimate: -7.238969575482958607E-01,
            delta: 1.842466275093522320E-14,
            absvalue: 1.156037579315941199E+00,
            asc: 9.306835267690729552E-01,
        },
    );
}

#[test]
fn qk33_f1() {
    test_qk(
        qk33,
        f1,
        0.0,
        1.0,
        Expect {
            estimate: 7.716049382562058245E-02,
            delta: 9.003490364414140615E-10,
            absvalue: 7.716049382562058245E-02,
            asc: 4.419742537125531667E-02,
        },
    );
}

#[test]
fn qk33_f2() {
    test_qk(
        qk33,
        f2,
        0.0,
        1.0,
        Expect {
            estimate: 2.128610145597360059E+01,
            delta: 3.378402800601331535E+01,
            absvalue: 2.128610145597360059E+01,
            asc: 3.378402800601331535E+01,
        },
    );
}

#[test]
fn qk33_f3() {
    test_qk(
        qk33,
        f3,
        0.3,
        2.71,
        Expect {
            estimate: -7.238969575482959717E-01,
            delta: 1.286818478735271366E-14,
            absvalue: 1.159063044265127296E+00,
            asc: 9.265867066767684568E-01,
        },
    );
}

#[test]
fn qk41_f1() {
    test_qk(
        qk41,
        f1,
        0.0,
        1.0,
        Expect {
            estimate: 7.716049382681375302E-02,
            delta: 9.576386660975511224E-11,
            absvalue: 7.716049382681375302E-02,
            asc: 4.421521169637691873E-02,
        },
    );
}

#[test]
fn qk41_f2() {
    test_qk(
        qk41,
        f2,
        0.0,
        1.0,
        Expect {
            estimate: 2.288677623903126701E+01,
            delta: 3.671538820274916048E+01,
            absvalue: 2.288677623903126701E+01,
            asc: 3.671538820274916048E+01,
        },
    );
}

#[test]
fn qk41_f3() {
    test_qk(
        qk41,
        f3,
        0.3,
        2.71,
        Expect {
            estimate: -7.238969575482959717E-01,
            delta: 1.286535726271015626E-14,
            absvalue: 1.158808363486595328E+00,
            asc: 9.264382258645686985E-01,
        },
    );
}

#[test]
fn qk49_f1() {
    test_qk(
        qk49,
        f1,
        0.0,
        1.0,
        Expect {
            estimate: 7.716049382705876536E-02,
            delta: 1.516099996697578150E-11,
            absvalue: 7.716049382705876536E-02,
            asc: 4.418030645992132577E-02,
        },
    );
}

#[test]
fn qk49_f2() {
    test_qk(
        qk49,
        f2,
        0.0,
        1.0,
        Expect {
            estimate: 2.420606744674078925E+01,
            delta: 3.914885058892737391E+01,
            absvalue: 2.420606744674078925E+01,
            asc: 3.914885058892737391E+01,
        },
    );
}

#[test]
fn qk49_f3() {
    test_qk(
        qk49,
        f3,
        0.3,
        2.71,
        Expect {
            estimate: -7.238969575482961938E-01,
            delta: 1.285880207548217569E-14,
            absvalue: 1.158217924711449687E+00,
            asc: 9.270623808448203995E-01,
        },
    );
}

#[test]
fn qk57_f1() {
    test_qk(
        qk57,
        f1,
        0.0,
        1.0,
        Expect {
            estimate: 7.716049382712462934E-02,
            delta: 3.163998417528912935E-12,
            absvalue: 7.716049382712462934E-02,
            asc: 4.419993837904739181E-02,
        },
    );
}

#[test]
fn qk57_f2() {
    test_qk(
        qk57,
        f2,
        0.0,
        1.0,
        Expect {
            estimate: 2.532725054513154817E+01,
            delta: 4.120456621814102505E+01,
            absvalue: 2.532725054513154817E+01,
            asc: 4.120456621814102505E+01,
        },
    );
}

#[test]
fn qk57_f3() {
    test_qk(
        qk57,
        f3,
        0.3,
        2.71,
        Expect {
            estimate: -7.238969575482961938E-01,
            delta: 1.286489289344262938E-14,
            absvalue: 1.158766536821391302E+00,
            asc: 9.270091578774536378E-01,
        },
    );
}
