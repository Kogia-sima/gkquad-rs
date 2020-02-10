#[cfg(all(not(feature = "std"), feature = "single"))]
mod mutex {
    use core::cell::UnsafeCell;
    use core::fmt;
    use core::ops::{Deref, DerefMut, Drop};
    use core::sync::atomic::{spin_loop_hint as cpu_relax, AtomicBool, Ordering};

    pub struct Mutex<T: ?Sized> {
        lock: AtomicBool,
        data: UnsafeCell<T>,
    }

    #[derive(Debug)]
    pub struct MutexGuard<'a, T: ?Sized + 'a> {
        lock: &'a AtomicBool,
        data: &'a mut T,
    }

    // Same unsafe impls as `std::sync::Mutex`
    unsafe impl<T: ?Sized + Send> Sync for Mutex<T> {}
    unsafe impl<T: ?Sized + Send> Send for Mutex<T> {}

    impl<T> Mutex<T> {
        pub const fn new(user_data: T) -> Mutex<T> {
            Mutex {
                lock: AtomicBool::new(false),
                data: UnsafeCell::new(user_data),
            }
        }
    }

    impl<T: ?Sized> Mutex<T> {
        fn obtain_lock(&self) {
            while self.lock.compare_and_swap(false, true, Ordering::Acquire) != false {
                // Wait until the lock looks unlocked before retrying
                while self.lock.load(Ordering::Relaxed) {
                    cpu_relax();
                }
            }
        }

        pub fn lock(&self) -> MutexGuard<T> {
            self.obtain_lock();
            MutexGuard {
                lock: &self.lock,
                data: unsafe { &mut *self.data.get() },
            }
        }

        pub fn try_lock(&self) -> Option<MutexGuard<T>> {
            if self.lock.compare_and_swap(false, true, Ordering::Acquire) == false {
                Some(MutexGuard {
                    lock: &self.lock,
                    data: unsafe { &mut *self.data.get() },
                })
            } else {
                None
            }
        }
    }

    impl<T: ?Sized + fmt::Debug> fmt::Debug for Mutex<T> {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            match self.try_lock() {
                Some(guard) => write!(f, "Mutex {{ data: ")
                    .and_then(|()| (&*guard).fmt(f))
                    .and_then(|()| write!(f, "}}")),
                None => write!(f, "Mutex {{ <locked> }}"),
            }
        }
    }

    impl<'a, T: ?Sized> Deref for MutexGuard<'a, T> {
        type Target = T;
        fn deref<'b>(&'b self) -> &'b T {
            &*self.data
        }
    }

    impl<'a, T: ?Sized> DerefMut for MutexGuard<'a, T> {
        fn deref_mut<'b>(&'b mut self) -> &'b mut T {
            &mut *self.data
        }
    }

    impl<'a, T: ?Sized> Drop for MutexGuard<'a, T> {
        fn drop(&mut self) {
            self.lock.store(false, Ordering::Release);
        }
    }
}

#[cfg(all(not(feature = "std"), feature = "single"))]
pub use mutex::*;
