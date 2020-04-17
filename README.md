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
OS Type: linux
CPU Architecture: x86_64
CPU Model Name: Intel(R) Core(TM) i5-8265U CPU @ 1.60GHz
SMBench Version: 0.1.0

simple:          88.420 ns [88.362 ns, 88.479 ns]    0 B (0 allocs)
singular_points: 1.3792 us [1.3790 us, 1.3794 us]   1 KB (2 allocs)
infinite_range:  500.75 ns [500.65 ns, 500.85 ns]   1 KB (2 allocs)
```

#### 2-dimentional integration

```
OS Type: linux
CPU Architecture: x86_64
CPU Model Name: Intel(R) Core(TM) i5-8265U CPU @ 1.60GHz
SMBench Version: 0.1.0

simple:          1.7229 us [1.7227 us, 1.7232 us]   32 B (1 allocs)
singular_points: 88.261 us [88.248 us, 88.274 us] 101 KB (4 allocs)
infinite_range:  166.71 us [166.68 us, 166.75 us] 101 KB (5 allocs)
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
