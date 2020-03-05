# gkquad
[![Build Status](https://travis-ci.org/Kogia-sima/gkquad-rs.svg?branch=master)](https://travis-ci.org/Kogia-sima/gkquad-rs)
[![Build status](https://ci.appveyor.com/api/projects/status/lw1w7sgf5fnrg9fg/branch/master?svg=true)](https://ci.appveyor.com/project/Kogiasima/gkquad-rs/branch/master)
[![Version](https://img.shields.io/crates/v/gkquad)](https://crates.io/crates/gkquad)
[![docs](https://docs.rs/gkquad/badge.svg)](https://docs.rs/gkquad)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/Kogia-sima/gkquad-rs/blob/master/LICENSE)
[![Twitter: Kogia_sima](https://img.shields.io/twitter/follow/Kogia\_sima.svg?style=social)](https://twitter.com/Kogia\_sima)

Numerical Integration Library for Rust

## Features

- Compatible with latest stable/beta/nightly Rust compiler
- Compatible with `no_std`
- Extremely fast and simple API
- Semi-verified computation (You can specify the maximum calculation tolerance)
- Lightweight (small dependencies)
- Highly extensible (you can implement a new algorithm)

#### Note: no\_std compatibility

`gkquad` depends on [alloc](https://doc.rust-lang.org/alloc/) crate, so you have to specify the global allocator in order to use `gkquad`.

If you want to use this crate in no\_std environment, you must disable the `std` feature flag.

```toml
[dependencies.gkquad]
version = "0.0.4"
default-features = false
features = ["single"]
```

## Performance

Here is the benchmark measured on the Intel Core(TM) i5 @ 1.60GHz (without turbo boost)

#### 1-dimentional integration

```
test infinite_range  ... bench:         467 ns/iter (+/- 13)
test simple          ... bench:          88 ns/iter (+/- 1)
test singular_points ... bench:       1,416 ns/iter (+/- 12)
```

#### 2-dimentional integration

```
test infinite_range  ... bench:     165,486 ns/iter (+/- 1,535)
test simple          ... bench:       1,761 ns/iter (+/- 11)
test singular_points ... bench:      88,945 ns/iter (+/- 3,548)
```

Source code can be found [here](https://github.com/Kogia-sima/gkquad-rs/blob/master/gkquad/benches).

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
