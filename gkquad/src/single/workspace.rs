use alloc::vec::Vec;
use std::cell::UnsafeCell;

use super::Interval;

/// Representing the subinterval and the integral estimates
#[derive(Clone, Debug)]
pub struct SubIntervalInfo {
    /// Subinterval range
    pub interval: Interval,
    /// Result of Gauss-Kronrod integration for subinterval
    pub estimate: f64,
    /// Absolute estimation error
    pub delta: f64,
    /// Recursion depth of subinterval
    pub level: usize
}

impl SubIntervalInfo {
    #[inline]
    pub fn new(interval: Interval, estimate: f64, delta: f64, level: usize) -> Self {
        Self {interval, estimate, delta, level}
    }
}

/// handles the memory for the subinterval ranges, results, and error estimates
#[derive(Clone, Debug)]
pub struct WorkSpace {
    /// maxerr = `subintervals[order[nrmax]].delta`. nrmax is normally 0 but will be
    /// positive if subdivision increased error estimate
    pub nrmax: usize,
    /// the partition index to be devided into sub partitions next
    pub i: usize,
    /// current maximum recursion depth
    pub maximum_level: usize,
    /// vector of dimension at least limit, the elements of which are the
    /// subintervals
    pub subintervals: Vec<SubIntervalInfo>,
    /// vector of dimension at least limit, the first k elements of which are
    /// indices to the error estimates over the subintervals, such that
    /// `subintervals[order[0]].delta, ..., subintervals[order[n - 1]].delta`
    /// form a decreasing sequence
    pub order: Vec<usize>,
}

impl WorkSpace {
    #[inline]
    pub fn new() -> WorkSpace {
        WorkSpace::with_capacity(0)
    }

    #[inline]
    pub fn with_capacity(n: usize) -> WorkSpace {
        WorkSpace {
            nrmax: 0,
            i: 0,
            maximum_level: 0,
            subintervals: Vec::with_capacity(n),
            order: Vec::with_capacity(n)
        }
    }

    /// return the number of subintervals
    #[inline]
    pub fn size(&self) -> usize {
        debug_assert_eq!(self.subintervals.len(), self.order.len());
        self.subintervals.len()
    }

    #[inline]
    pub fn capacity(&self) -> usize {
        debug_assert_eq!(self.subintervals.capacity(), self.order.capacity());
        self.subintervals.capacity()
    }

    #[inline]
    pub fn reserve(&mut self, n: usize) {
        self.subintervals.reserve(n);
        self.order.reserve(n);
    }

    #[inline]
    pub fn clear(&mut self) {
        self.subintervals.clear();
        self.order.clear();
        self.i = 0;
        self.nrmax = 0;
        self.maximum_level = 0;
    }

    #[inline]
    pub fn push(&mut self, subinterval: SubIntervalInfo) {
        self.subintervals.push(subinterval);
        self.order.push(self.order.len());
    }

    pub fn sort_results(&mut self) {
        debug_assert_eq!(self.subintervals.len(), self.order.len());
        let nint = self.size();
        if nint == 0 {
            return;
        }

        let subintervals = &mut self.subintervals;
        let order = &mut self.order;

        for i in 0..nint {
            let i1 = order[i];
            let mut e1 = subintervals[i1].delta;
            let mut i_max = i1;

            for j in i + 1..nint {
                let i2 = order[j];
                let e2 = subintervals[i2].delta;

                if e2 >= e1 {
                    i_max = i2;
                    e1 = e2;
                }
            }

            if i_max != i1 {
                order[i] = order[i_max];
                order[i_max] = i1;
            }
        }

        self.i = order[0];
    }

    /// append the newly-created subintervals to the list
    pub fn update(self: &mut WorkSpace, s1: SubIntervalInfo, s2: SubIntervalInfo) {
        debug_assert_eq!(self.subintervals.len(), self.order.len());
        let new_level = self.subintervals[self.i].level + 1;

        if s2.delta > s1.delta {
            self.subintervals[self.i] = s2;
            self.subintervals.push(s1);
        } else {
            self.subintervals[self.i] = s1;
            self.subintervals.push(s2);
        }
        self.order.push(self.order.len());

        if new_level > self.maximum_level {
            self.maximum_level = new_level;
        }

        self.qpsrt();
        debug_assert_eq!(self.subintervals.len(), self.order.len());
    }

    fn qpsrt(&mut self) {
        let last = self.size() - 1;
        let limit = self.capacity();
        let subintervals = &self.subintervals;
        let order = &mut self.order;
        let mut i_nrmax = self.nrmax;
        let mut i_maxdelta = order[i_nrmax];

        // Check whether the list contains more than two error estimates

        if last < 2 {
            order[0] = 0;
            order[1] = 1;
            self.i = i_maxdelta;
            return;
        }

        // search the position for inserting the `order[i_nrmax]`
        let deltamax = subintervals[i_maxdelta].delta;
        while i_nrmax > 0 && deltamax > subintervals[order[i_nrmax - 1]].delta {
            order[i_nrmax] = order[i_nrmax - 1];
            i_nrmax -= 1;
        }

        // If last < (limit / 2 + 2), then the remaining subintervals will not
        // divided any more, since the algorithm will split the subinterval with
        // largest delta.

        let top = if last < (limit / 2 + 2) {
            last
        } else {
            limit - last + 1
        };

        // search the position for inserting the `order[i_nrmax]`
        let mut i = i_nrmax + 1;
        while i < top && deltamax < subintervals[order[i]].delta {
            order[i - 1] = order[i];
            i += 1;
        }
        order[i - 1] = i_maxdelta;

        // search the position for inserting the `order[last]`

        let errmin = subintervals[last].delta;
        let mut k = top as i64 - 1;

        while k >= i as i64 - 1 && errmin >= subintervals[order[k as usize]].delta {
            order[k as usize + 1] = order[k as usize];

            k -= 1;
        }

        order[(k + 1) as usize] = last;

        // Set i_max and e_max

        i_maxdelta = order[i_nrmax];
        self.i = i_maxdelta;
        self.nrmax = i_nrmax;
    }

    #[inline]
    pub fn maximum_level(&self) -> usize {
        self.maximum_level
    }

    /// The smallest interval has the largest error. Before bisecting decrease the
    /// sum of the errors over the larger intervals (error_over_large_intervals)
    /// and perform extrapolation.
    pub(crate) fn increase_nrmax(&mut self) -> bool {
        let id = self.nrmax;
        let order = &self.order;
        let limit = self.capacity();
        let last = self.size() - 1;

        let jupbnd = if last > (1 + limit / 2) {
            limit + 1 - last
        } else {
            last
        };

        // 最小でない部分区間のうち、最も誤差が大きい部分を次に分割する
        for _ in id..=jupbnd {
            let i_max = order[self.nrmax];
            self.i = i_max;

            if self.subintervals[i_max].level < self.maximum_level {
                return true;
            }

            self.nrmax += 1;
        }

        // large interval not found
        false
    }

    #[inline]
    pub(crate) fn reset_nrmax(&mut self) {
        self.nrmax = 0;
        self.i = self.order[0];
    }

    /// retrieve the next subinterval
    #[inline]
    pub fn get(&self) -> &SubIntervalInfo {
        &self.subintervals[self.i]
    }

    /// calculate the sum of integral estimates for all subintervals
    #[inline]
    pub fn sum_results(self: &WorkSpace) -> f64 {
        self.subintervals.iter().map(|s| s.estimate).sum()
    }
}

impl Default for WorkSpace {
    fn default() -> WorkSpace {
        WorkSpace::new()
    }
}

#[cfg(feature = "std")]
thread_local! {
    static WORKSPACE: UnsafeCell<WorkSpace> = UnsafeCell::new(WorkSpace::new());
}

#[cfg(feature = "std")]
#[derive(Clone)]
pub struct WorkSpaceProvider;

#[cfg(feature = "std")]
impl WorkSpaceProvider {
    pub fn new() -> Self {
        Self
    }

    pub unsafe fn get_mut(&self) -> &mut WorkSpace {
        let ptr = WORKSPACE.with(|v| v.get());
        &mut *ptr
    }
}

#[cfg(not(feature = "std"))]
pub struct WorkSpaceProvider {
    workspace: UnsafeCell<WorkSpace>
}

#[cfg(not(feature = "std"))]
impl WorkSpaceProvider {
    pub fn new() -> Self {
        Self {
            workspace: UnsafeCell::new(WorkSpace::new())
        }
    }

    pub unsafe fn get_mut(&self) -> &mut WorkSpace {
        let ptr = self.workspace.get();
        &mut *ptr
    }
}

#[cfg(not(feature = "std"))]
impl Clone for WorkSpaceProvider {
    fn clone(&self) -> Self {
        let ptr = self.workspace.get();
        unsafe {
            Self {
                workspace: UnsafeCell::new((*ptr).clone())
            }
        }
    }
}

#[cfg(docsrs)]
impl !Sync for WorkSpaceProvider {}
