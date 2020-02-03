macro_rules! assert_rel {
    ($left:expr, $right:expr, $tolerance:expr) => ({
        match (&$left, &$right) {
            (left_val, right_val) => {
                if (*right_val == 0.0 && *left_val >= $tolerance) || (*right_val != 0.0 && (*left_val - *right_val).abs() / right_val.abs() > $tolerance) {
                    // The reborrows below are intentional. Without them, the stack slot for the
                    // borrow is initialized even before the values are compared, leading to a
                    // noticeable slow down.
                    panic!(r#"assertion failed: `(left == right)`
  left: `{:?}`,
 right: `{:?}`"#, &*left_val, &*right_val)
                }
            }
        }
    });
}

