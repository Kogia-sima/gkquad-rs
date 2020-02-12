# gkquad
[![Version](https://img.shields.io/crates/v/gkquad)](https://crates.io/crates/gkquad)
[![docs](https://docs.rs/gkquad/badge.svg)](https://docs.rs/gkquad)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/Kogia-sima/gkquad-rs/blob/master/LICENSE)
[![Twitter: Kogia_sima](https://img.shields.io/twitter/follow/Kogia\_sima.svg?style=social)](https://twitter.com/Kogia\_sima)

Numerical Integration Library for Rust

## Features

- Compatible with latest stable/beta/nightly Rust compiler
- Compatible with `no_std`
- Super fast and simple API
- Semi verified computation (You can specify the maximum calculation tolerance)
- Lightweight (small dependencies)
- Highly extensible (you can implement a new algorithm)

#### Note: no\_std compatibility

`gkquad` depends [alloc](https://doc.rust-lang.org/alloc/) crate, so you have to specify the global allocator in order to use `gkquad`.

If you want to use this crate in no\_std environment, you must disable the `std` feature flag.

```toml
[dependencies.gkquad]
version = "0.0.1"
default-features = false
features = ["single"]
```

## Performance

Here is the benchmark measured on Intel Core(TM) i5 @ 1.60GHz

```console
$ cargo bench

     Running /tmp/gkquad-rs/target/release/deps/single-3b52efd7f739cf4b

running 3 tests
test infinite_interval ... bench:         211 ns/iter (+/- 4)
test simple            ... bench:          38 ns/iter (+/- 1)
test singular_points   ... bench:         621 ns/iter (+/- 42)

test result: ok. 0 passed; 0 failed; 0 ignored; 3 measured
```

Source code can be found [here](https://github.com/Kogia-sima/gkquad-rs/blob/master/gkquad/benches/single.rs).

## Author

👤 **Kogia-sima**

* Twitter: [@Kogia\_sima](https://twitter.com/Kogia\_sima)
* Github: [@Kogia-sima](https://github.com/Kogia-sima)

## 🤝 Contributing

Contributions, issues and feature requests are welcome!

Feel free to check [issues page](https://github.com/Kogia-sima/gkquad-rs/issues). 

## Show your support

Give a ⭐️ if this project helped you!


## 📝 License

Copyright © 2020 [Kogia-sima](https://github.com/Kogia-sima).

This project is [MIT](https://github.com/Kogia-sima/gkquad-rs/blob/master/LICENSE) licensed.

***
_This README was generated with ❤️ by [readme-md-generator](https://github.com/kefranabg/readme-md-generator)_
