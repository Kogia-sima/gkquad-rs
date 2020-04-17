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

```
OS Type: linux
CPU Architecture: x86_64
CPU Model Name: Intel(R) Core(TM) i5-8265U CPU @ 1.60GHz
Intel Turbo Boost Technology: disabled
SMBench Version: 0.1.0

# single (gkquad/benches/single.rs)
Benchmark              Time                  95% CI         Allocation
----------------------------------------------------------------------
simple            88.351 ns  [88.223 ns, 88.479 ns]     0 B (0 allocs)
singular_points   1.3804 us  [1.3799 us, 1.3810 us]    1 KB (2 allocs)
infinite_range    500.17 ns  [499.87 ns, 500.47 ns]    1 KB (2 allocs)

# double (gkquad/benches/double.rs)
Benchmark              Time                  95% CI         Allocation
----------------------------------------------------------------------
simple            1.7205 us  [1.7197 us, 1.7213 us]    32 B (1 allocs)
singular_points   88.415 us  [88.354 us, 88.476 us]  101 KB (4 allocs)
infinite_range    167.16 us  [166.76 us, 167.55 us]  101 KB (5 allocs)
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
