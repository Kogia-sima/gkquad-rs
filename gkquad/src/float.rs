#![allow(clippy::many_single_char_names)]

pub trait Float: Copy {
    fn abs(self) -> Self;
    fn sqrt(self) -> Self;
}

impl Float for f64 {
    #[inline]
    fn abs(self: f64) -> f64 {
        f64::from_bits(self.to_bits() & (std::u64::MAX / 2))
    }

    fn sqrt(self: f64) -> f64 {
        #[cfg(target_feature = "sse2")]
        {
            #[cfg(target_arch = "x86")]
            use core::arch::x86::*;
            #[cfg(target_arch = "x86_64")]
            use core::arch::x86_64::*;
            unsafe {
                let m = _mm_set_sd(self);
                let m_sqrt = _mm_sqrt_pd(m);
                _mm_cvtsd_f64(m_sqrt)
            }
        }
        #[cfg(not(target_feature = "sse2"))]
        {
            use core::num::Wrapping;

            const TINY: f64 = 1.0e-300;

            let mut z: f64;
            let sign: Wrapping<u32> = Wrapping(0x80000000);
            let mut ix0: i32;
            let mut s0: i32;
            let mut q: i32;
            let mut m: i32;
            let mut t: i32;
            let mut i: i32;
            let mut r: Wrapping<u32>;
            let mut t1: Wrapping<u32>;
            let mut s1: Wrapping<u32>;
            let mut ix1: Wrapping<u32>;
            let mut q1: Wrapping<u32>;

            ix0 = (self.to_bits() >> 32) as i32;
            ix1 = Wrapping(self.to_bits() as u32);

            /* take care of Inf and NaN */
            if (ix0 & 0x7ff00000) == 0x7ff00000 {
                return self * self + self; /* sqrt(NaN)=NaN, sqrt(+inf)=+inf, sqrt(-inf)=sNaN */
            }
            /* take care of zero */
            if ix0 <= 0 {
                if ((ix0 & !(sign.0 as i32)) | ix1.0 as i32) == 0 {
                    return self; /* sqrt(+-0) = +-0 */
                }
                if ix0 < 0 {
                    return (self - self) / (self - self); /* sqrt(-ve) = sNaN */
                }
            }
            /* normalize x */
            m = ix0 >> 20;
            if m == 0 {
                /* subnormal x */
                while ix0 == 0 {
                    m -= 21;
                    ix0 |= (ix1 >> 11).0 as i32;
                    ix1 <<= 21;
                }
                i = 0;
                while (ix0 & 0x00100000) == 0 {
                    i += 1;
                    ix0 <<= 1;
                }
                m -= i - 1;
                ix0 |= (ix1 >> (32 - i) as usize).0 as i32;
                ix1 = ix1 << i as usize;
            }
            m -= 1023; /* unbias exponent */
            ix0 = (ix0 & 0x000fffff) | 0x00100000;
            if (m & 1) == 1 {
                /* odd m, double x to make it even */
                ix0 += ix0 + ((ix1 & sign) >> 31).0 as i32;
                ix1 += ix1;
            }
            m >>= 1; /* m = [m/2] */

            /* generate sqrt(x) bit by bit */
            ix0 += ix0 + ((ix1 & sign) >> 31).0 as i32;
            ix1 += ix1;
            q = 0; /* [q,q1] = sqrt(x) */
            q1 = Wrapping(0);
            s0 = 0;
            s1 = Wrapping(0);
            r = Wrapping(0x00200000); /* r = moving bit from right to left */

            while r != Wrapping(0) {
                t = s0 + r.0 as i32;
                if t <= ix0 {
                    s0 = t + r.0 as i32;
                    ix0 -= t;
                    q += r.0 as i32;
                }
                ix0 += ix0 + ((ix1 & sign) >> 31).0 as i32;
                ix1 += ix1;
                r >>= 1;
            }

            r = sign;
            while r != Wrapping(0) {
                t1 = s1 + r;
                t = s0;
                if t < ix0 || (t == ix0 && t1 <= ix1) {
                    s1 = t1 + r;
                    if (t1 & sign) == sign && (s1 & sign) == Wrapping(0) {
                        s0 += 1;
                    }
                    ix0 -= t;
                    if ix1 < t1 {
                        ix0 -= 1;
                    }
                    ix1 -= t1;
                    q1 += r;
                }
                ix0 += ix0 + ((ix1 & sign) >> 31).0 as i32;
                ix1 += ix1;
                r >>= 1;
            }

            /* use floating add to find out rounding direction */
            if (ix0 as u32 | ix1.0) != 0 {
                z = 1.0 - TINY; /* raise inexact flag */
                if z >= 1.0 {
                    z = 1.0 + TINY;
                    if q1.0 == 0xffffffff {
                        q1 = Wrapping(0);
                        q += 1;
                    } else if z > 1.0 {
                        if q1.0 == 0xfffffffe {
                            q += 1;
                        }
                        q1 += Wrapping(2);
                    } else {
                        q1 += q1 & Wrapping(1);
                    }
                }
            }
            ix0 = (q >> 1) + 0x3fe00000;
            ix1 = q1 >> 1;
            if (q & 1) == 1 {
                ix1 |= sign;
            }
            ix0 += m << 20;
            f64::from_bits((ix0 as u64) << 32 | ix1.0 as u64)
        }
    }
}
