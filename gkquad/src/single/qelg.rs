//# struct extrapolation_table
//#   {
//#     size_t n;
//#     double rlist2[52];
//#     size_t nres;
//#     double res3la[3];
//#   };

pub struct ExtrapolationTable {
    /// rlist2\[n\] contains the new element in the first
    /// column of the epsilon table
    pub n: usize,
    /// the vector containing the elements of the two lower
    /// diagonals of the triangular epsilon table. the elements
    /// are numbered starting at the right^hand corner of the
    /// triangle
    pub rlist2: [f64; 52],
    /// number of calls to the routine
    pub nres: usize,
    /// the vector containing the last 3 results
    pub res3la: [f64; 3],
}

impl ExtrapolationTable {
    /// append calculation result.
    ///
    /// size of rlist2 is limited, but when it get full, the contents are
    /// shifted and the oldest result are removed from array
    #[inline]
    pub fn append(&mut self, y: f64) {
        let n = self.n;
        self.rlist2[n] = y;
        self.n += 1;
    }

    /// the routine determines the limit of a given sequence of
    /// approximations, by means of the epsilon algorithm of
    /// p. wynn. an estimate of the absolute error is also given.
    /// the condensed epsilon table is computed. only those
    /// elements needed for the computation of the next diagonal
    /// are preserved.
    ///
    /// * `result` - resulting approximation to the integral
    /// * `abserr` - estimate of the absolute error computed from
    ///              result and the 3 previous results
    ///
    pub fn qelg(&mut self, result: &mut f64, abserr: &mut f64) {
        //# double *epstab = table->rlist2;
        //# double *res3la = table->res3la;
        //# const size_t n = table->n - 1;
        //#
        //# const double current = epstab[n];
        //#
        //# double absolute = GSL_DBL_MAX;
        //# double relative = 5 * GSL_DBL_EPSILON * fabs (current);
        //#
        //# const size_t newelm = n / 2;
        //# const size_t n_orig = n;
        //# size_t n_final = n;
        //# size_t i;
        //#
        //# const size_t nres_orig = table->nres;
        //#
        //# *result = current;
        //# *abserr = GSL_DBL_MAX;

        let epstab = &mut self.rlist2;
        let res3la = &mut self.res3la;
        let n = self.n - 1;

        let current = unsafe { *epstab.get_unchecked(n) };

        let newelm = n / 2;
        let n_orig = n;
        let mut n_final = n;

        let nres_orig = self.nres;

        *result = current;
        *abserr = core::f64::MAX;

        //# if (n < 2)
        //#   {
        //#     *result = current;
        //#     *abserr = GSL_MAX_DBL (absolute, relative);
        //#     return;
        //#   }

        if n < 2 {
            *result = current;
            *abserr = core::f64::MAX;
            return;
        }

        unsafe {
            //# epstab[n + 2] = epstab[n];
            //# epstab[n] = GSL_DBL_MAX;

            let mut ep = epstab.as_mut_ptr().add(n);
            *ep.add(2) = *ep;
            *ep = core::f64::MAX;

            for i in 0..newelm {
                //# double res = epstab[n - 2 * i + 2];
                //# double e0 = epstab[n - 2 * i - 2];
                //# double e1 = epstab[n - 2 * i - 1];
                //# double e2 = res;

                let mut res = *ep.add(2);
                let e0 = *ep.sub(2);
                let e1 = *ep.sub(1);
                let e2 = res;

                //# double e1abs = fabs (e1);
                //# double delta2 = e2 - e1;
                //# double err2 = fabs (delta2);
                //# double tol2 = GSL_MAX_DBL (fabs (e2), e1abs) * GSL_DBL_EPSILON;
                //# double delta3 = e1 - e0;
                //# double err3 = fabs (delta3);
                //# double tol3 = GSL_MAX_DBL (e1abs, fabs (e0)) * GSL_DBL_EPSILON;

                let elabs = e1.abs();
                let delta2 = e2 - e1;
                let err2 = delta2.abs();
                let tol2 = f64::max(e2.abs(), elabs) * core::f64::EPSILON;
                let delta3 = e1 - e0;
                let err3 = delta3.abs();
                let tol3 = f64::max(elabs, e0.abs()) * core::f64::EPSILON;

                //# if (err2 <= tol2 && err3 <= tol3)
                //#   {
                //#     /* If e0, e1 and e2 are equal to within machine accuracy,
                //#        convergence is assumed. */
                //#     *result = res;
                //#     absolute = err2 + err3;
                //#     relative = 5 * GSL_DBL_EPSILON * fabs (res);
                //#     *abserr = GSL_MAX_DBL (absolute, relative);
                //#     return;
                //#   }

                if err2 <= tol2 && err3 <= tol3 {
                    *result = res;
                    let absolute = err2 + err3;
                    let relative = 5. * core::f64::EPSILON * res.abs();
                    *abserr = f64::max(absolute, relative);
                    return;
                }

                //# e3 = epstab[n - 2 * i];
                //# epstab[n - 2 * i] = e1;
                //# delta1 = e1 - e3;
                //# err1 = fabs (delta1);
                //# tol1 = GSL_MAX_DBL (e1abs, fabs (e3)) * GSL_DBL_EPSILON;

                let e3 = *ep;
                *ep = e1;
                let delta1 = e1 - e3;
                let err1 = delta1.abs();
                let tol1 = f64::max(elabs, e3.abs()) * core::f64::EPSILON;

                //# if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3)
                //#   {
                //#     n_final = 2 * i;
                //#     break;
                //#   }
                //# ss = (1 / delta1 + 1 / delta2) - 1 / delta3;

                if err1 <= tol1 || err2 <= tol2 || err3 <= tol3 {
                    n_final = 2 * i;
                    break;
                }

                let ss = (1. / delta1 + 1. / delta2) - 1. / delta3;

                // Test to detect irregular behaviour in the table, and
                // eventually omit a part of the table by adjusting the value of
                // n.

                //# if (fabs (ss * e1) <= 0.0001)
                //#   {
                //#     n_final = 2 * i;
                //#     break;
                //#   }

                if (ss * e1).abs() <= 0.0001 {
                    n_final = 2 * i;
                    break;
                }

                // Compute a new element and eventually adjust the value of
                // result.

                //# res = e1 + 1 / ss;
                //# epstab[n - 2 * i] = res;

                res = e1 + 1. / ss;
                *ep = res;

                //# {
                //#   const double error = err2 + fabs (res - e2) + err3;
                //#   if (error <= *abserr)
                //#     {
                //#       *abserr = error;
                //#       *result = res;
                //#     }
                //# }

                let error = err2 + (res - e2).abs() + err3;

                if error <= *abserr {
                    *abserr = error;
                    *result = res;
                }

                ep = ep.sub(2);
            }

            // shift table

            //# const size_t limexp = 50 - 1;

            //# if (n_final == limexp)
            //#   {
            //#     n_final = 2 * (limexp / 2);
            //#   }

            const LIMEXP: usize = 50 - 1;

            if n_final == LIMEXP {
                n_final = 2 * (LIMEXP / 2);
            }

            //# if (n_orig % 2 == 1)
            //#   {
            //#     for (i = 0; i <= newelm; i++)
            //#       {
            //#         epstab[1 + i * 2] = epstab[i * 2 + 3];
            //#       }
            //#   }
            //# else
            //#   {
            //#     for (i = 0; i <= newelm; i++)
            //#       {
            //#         epstab[i * 2] = epstab[i * 2 + 2];
            //#       }
            //#   }

            let mut ep = if n_orig & 1 == 1 {
                epstab.as_mut_ptr().add(1)
            } else {
                epstab.as_mut_ptr()
            };
            let ep_end = ep.add(newelm * 2);

            while ep <= ep_end {
                *ep = *ep.add(2);
                ep = ep.add(2);
            }

            //# if (n_orig != n_final)
            //#   {
            //#     for (i = 0; i <= n_final; i++)
            //#       {
            //#         epstab[i] = epstab[n_orig - n_final + i];
            //#       }
            //#   }
            //# table->n = n_final + 1;

            if n_orig != n_final {
                std::ptr::copy(
                    epstab.as_ptr().add(n_orig - n_final),
                    epstab.as_mut_ptr(),
                    n_final + 1,
                );
            }
            self.n = n_final + 1;

            //# if (nres_orig < 3)
            //#   {
            //#     res3la[nres_orig] = *result;
            //#     *abserr = GSL_DBL_MAX;
            //#   }
            //# else
            //#   {                           /* Compute error estimate */
            //#     *abserr = (fabs (*result - res3la[2]) + fabs (*result - res3la[1])
            //#                + fabs (*result - res3la[0]));

            //#     res3la[0] = res3la[1];
            //#     res3la[1] = res3la[2];
            //#     res3la[2] = *result;
            //#   }

            if let Some(res) = res3la.get_mut(nres_orig) {
                *res = *result;
                *abserr = core::f64::MAX;
            } else {
                /* Compute error estimate */
                *abserr = (*result - res3la[2]).abs()
                    + (*result - res3la[1]).abs()
                    + (*result - res3la[0]).abs();

                /* append result at the tail of res3la */
                res3la[0] = res3la[1];
                res3la[1] = res3la[2];
                res3la[2] = *result;
            }
        }

        // In QUADPACK the variable table->nres is incremented at the top of
        // qelg, so it increases on every call. This leads to the array
        // res3la being accessed when its elements are still undefined, so I
        // have moved the update to this point so that its value more
        // useful.

        //# table->nres = nres_orig + 1;
        //# *abserr = GSL_MAX_DBL (*abserr, 5 * GSL_DBL_EPSILON * fabs (*result));

        self.nres = nres_orig + 1;
        *abserr = f64::max(*abserr, 5. * core::f64::EPSILON * result.abs());
    }
}

impl Default for ExtrapolationTable {
    #[inline]
    fn default() -> Self {
        Self {
            n: 0,
            rlist2: [0f64; 52],
            nres: 0,
            res3la: [0f64; 3],
        }
    }
}
