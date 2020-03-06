use alloc::vec::Vec;

use super::Range;

/// Representing the subrange and the integral estimates
#[derive(Clone, Debug)]
pub struct SubRangeInfo {
    /// Subrange
    pub range: Range,
    /// Result of Gauss-Kronrod integration for subrange
    pub estimate: f64,
    /// Absolute estimation error
    pub delta: f64,
    /// Recursion depth of subrange
    pub level: usize,
}

impl SubRangeInfo {
    #[inline]
    pub fn new(range: Range, estimate: f64, delta: f64, level: usize) -> Self {
        Self {
            range,
            estimate,
            delta,
            level,
        }
    }
}

/// handles the memory for the subrange ranges, results, and error estimates
#[derive(Clone, Debug)]
pub struct WorkSpace {
    /// maxerr = `subranges[order[nrmax]].delta`. nrmax is normally 0 but will be
    /// positive if subdivision increased error estimate
    pub nrmax: usize,
    /// the partition index to be devided into sub partitions next
    pub i: usize,
    /// current maximum recursion depth
    pub maximum_level: usize,
    /// vector of dimension at least limit, the elements of which are the
    /// subranges
    pub subranges: Vec<SubRangeInfo>,
    /// vector of dimension at least limit, the first k elements of which are
    /// indices to the error estimates over the subranges, such that
    /// `subranges[order[0]].delta, ..., subranges[order[n - 1]].delta`
    /// form a decreasing sequence
    pub order: Vec<usize>,
}

impl WorkSpace {
    #[inline]
    pub const fn new() -> WorkSpace {
        WorkSpace {
            nrmax: 0,
            i: 0,
            maximum_level: 0,
            subranges: Vec::new(),
            order: Vec::new(),
        }
    }

    #[inline]
    pub fn with_capacity(n: usize) -> WorkSpace {
        WorkSpace {
            nrmax: 0,
            i: 0,
            maximum_level: 0,
            subranges: Vec::with_capacity(n),
            order: Vec::with_capacity(n),
        }
    }

    /// return the number of subranges
    #[inline]
    pub fn size(&self) -> usize {
        debug_assert_eq!(self.subranges.len(), self.order.len());
        self.subranges.len()
    }

    #[inline]
    pub fn capacity(&self) -> usize {
        debug_assert_eq!(self.subranges.capacity(), self.order.capacity());
        self.subranges.capacity()
    }

    #[inline]
    pub fn reserve(&mut self, n: usize) {
        self.subranges.reserve(n);
        self.order.reserve(n);
    }

    #[inline]
    pub fn clear(&mut self) {
        self.subranges.clear();
        self.order.clear();
        self.i = 0;
        self.nrmax = 0;
        self.maximum_level = 0;
    }

    #[inline]
    pub fn push(&mut self, subrange: SubRangeInfo) {
        self.subranges.push(subrange);
        self.order.push(self.order.len());
    }

    pub fn sort_results(&mut self) {
        debug_assert_eq!(self.subranges.len(), self.order.len());
        let nint = self.size();
        if nint == 0 {
            return;
        }

        let subranges = &mut self.subranges;
        let order = &mut self.order;

        for i in 0..nint {
            let i1 = order[i];
            let mut e1 = subranges[i1].delta;
            let mut i_max = i1;

            for &i2 in order.iter().skip(i + 1) {
                let e2 = subranges[i2].delta;

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

    /// append the newly-created subranges to the list
    pub fn update(self: &mut WorkSpace, s1: SubRangeInfo, s2: SubRangeInfo) {
        debug_assert_eq!(self.subranges.len(), self.order.len());
        let new_level = self.subranges[self.i].level + 1;

        if s2.delta > s1.delta {
            self.subranges[self.i] = s2;
            self.subranges.push(s1);
        } else {
            self.subranges[self.i] = s1;
            self.subranges.push(s2);
        }
        self.order.push(self.order.len());

        if new_level > self.maximum_level {
            self.maximum_level = new_level;
        }

        self.qpsrt();
        debug_assert_eq!(self.subranges.len(), self.order.len());
    }

    #[inline(always)]
    fn qpsrt(&mut self) {
        let last = self.size() - 1;
        let limit = self.capacity();
        let subranges = &self.subranges;
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
        let deltamax = subranges[i_maxdelta].delta;
        while i_nrmax > 0 && deltamax > subranges[order[i_nrmax - 1]].delta {
            order[i_nrmax] = order[i_nrmax - 1];
            i_nrmax -= 1;
        }

        // If last < (limit / 2 + 2), then the remaining subranges will not
        // divided any more, since the algorithm will split the subrange with
        // largest delta.

        let top = if last < (limit / 2 + 2) {
            last
        } else {
            limit - last + 1
        };

        // search the position for inserting the `order[i_nrmax]`
        let mut i = i_nrmax + 1;
        while i < top && deltamax < subranges[order[i]].delta {
            order[i - 1] = order[i];
            i += 1;
        }
        order[i - 1] = i_maxdelta;

        // search the position for inserting the `order[last]`

        let errmin = subranges[last].delta;
        let mut k = top as i64 - 1;

        while k >= i as i64 - 1 && errmin >= subranges[order[k as usize]].delta {
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

    /// The smallest range has the largest error. Before bisecting decrease the
    /// sum of the errors over the larger ranges (error_over_large_ranges)
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

            if self.subranges[i_max].level < self.maximum_level {
                return true;
            }

            self.nrmax += 1;
        }

        // large range not found
        false
    }

    #[inline]
    pub(crate) fn reset_nrmax(&mut self) {
        self.nrmax = 0;
        self.i = self.order[0];
    }

    /// retrieve the next subrange
    #[inline]
    pub fn get(&self) -> &SubRangeInfo {
        &self.subranges[self.i]
    }

    /// calculate the sum of integral estimates for all subranges
    #[inline]
    pub fn sum_results(self: &WorkSpace) -> f64 {
        self.subranges.iter().map(|s| s.estimate).sum()
    }
}

impl Default for WorkSpace {
    #[inline]
    fn default() -> WorkSpace {
        WorkSpace::new()
    }
}
