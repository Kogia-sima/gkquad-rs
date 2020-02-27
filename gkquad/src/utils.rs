mod cowmut {
    use core::fmt::{self, Debug, Display};
    pub use core::ops::{Deref, DerefMut};

    pub enum CowMut<'a, T: 'a> {
        Owned(T),
        Borrowed(&'a mut T),
    }

    impl<'a, T: 'a> Deref for CowMut<'a, T> {
        type Target = T;

        #[inline]
        fn deref(&self) -> &T {
            self.as_ref()
        }
    }

    impl<'a, T: 'a> DerefMut for CowMut<'a, T> {
        #[inline]
        fn deref_mut(&mut self) -> &mut T {
            self.as_mut()
        }
    }

    impl<'a, T: 'a> AsRef<T> for CowMut<'a, T> {
        #[inline]
        fn as_ref(&self) -> &T {
            match *self {
                CowMut::Owned(ref t) => t,
                CowMut::Borrowed(&mut ref t) => t,
            }
        }
    }

    impl<'a, T: 'a> AsMut<T> for CowMut<'a, T> {
        #[inline]
        fn as_mut(&mut self) -> &mut T {
            match *self {
                CowMut::Owned(ref mut t) => t,
                CowMut::Borrowed(&mut ref mut t) => t,
            }
        }
    }

    impl<'a, T: Clone + 'a> Clone for CowMut<'a, T> {
        #[inline]
        fn clone(&self) -> Self {
            match *self {
                CowMut::Owned(ref t) => CowMut::Owned(t.clone()),
                CowMut::Borrowed(&mut ref t) => CowMut::Owned(t.clone()),
            }
        }
    }

    impl<'a, T: Display + 'a> Display for CowMut<'a, T> {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            Display::fmt(self.as_ref(), f)
        }
    }

    impl<'a, T: Debug + 'a> Debug for CowMut<'a, T> {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            Debug::fmt(self.as_ref(), f)
        }
    }

    impl<'a, T: PartialEq + 'a> PartialEq for CowMut<'a, T> {
        #[inline]
        fn eq(&self, other: &Self) -> bool {
            self.as_ref().eq(other.as_ref())
        }
    }

    impl<'a, T: Eq + 'a> Eq for CowMut<'a, T> {}

    impl<'a, T: Default + 'a> Default for CowMut<'a, T> {
        fn default() -> Self {
            CowMut::Owned(T::default())
        }
    }
}

pub use cowmut::*;
