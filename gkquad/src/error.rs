use std::fmt::Debug;

macro_rules! impl_error {
    (
        $(#[$outer:meta])*
        pub enum $name:ident {
            $($chunk:tt)*
        }
    ) => {
        impl_error!(
            @PARSE_INNER
            meta $(#[$outer])*
            name $name
            queue [ $($chunk)* ]
            items []
        );
    };
    (
        @PARSE_INNER
        meta $(#[$outer:meta])*
        name $name:ident
        queue []
        items [ $($item:tt)* ]
    ) => {
        impl_error!(
            @WRITE
            meta $(#[$outer])*
            name $name
            items [ $($item)* ]
        );
    };
    (
        @PARSE_INNER
        meta $(#[$outer:meta])*
        name $name:ident
        queue [ $var:ident { $desc:expr, $msg:expr, }, $($chunk:tt)* ]
        items [ $($item:tt)* ]
    ) => {
        impl_error!(
            @PARSE_INNER
            meta $(#[$outer])*
            name $name
            queue [ $($chunk)* ]
            items [
                $($item)*
                {
                    var $var
                    args []
                    args2 []
                    args3 []
                    desc ($desc)
                    msg ($msg)
                }
            ]
        );
    };
    (
        @PARSE_INNER
        meta $(#[$outer:meta])*
        name $name:ident
        queue [ $var:ident ($($arg:ident: $argtype:ty),*) { $desc:expr, $msg:expr, }, $($chunk:tt)* ]
        items [ $($item:tt)* ]
    ) => {
        impl_error!(
            @PARSE_INNER
            meta $(#[$outer])*
            name $name
            queue [ $($chunk)* ]
            items [
                $($item)*
                {
                    var $var
                    args [
                        $({
                            arg $arg
                            argtype $argtype
                        })*
                    ]
                    args2 [ ($($arg),*) ]
                    args3 [ ($($argtype),*) ]
                    desc ($desc)
                    msg ($msg)
                }
            ]
        );
    };
    (
        @WRITE
        meta $(#[$outer:meta])*
        name $name:ident
        items [
            $({
                var $var:ident
                args [
                    $({
                        arg $arg:ident
                        argtype $argtype:ty
                    })*
                ]
                args2 [ $($args2:tt)* ]
                args3 [ $($args3:tt)* ]
                desc ($desc:expr)
                msg ($msg:expr)
            })*
        ]
    ) => {
        $(#[$outer])*
        pub enum $name {
            $(
                #[doc=$msg]
                $var $($args3)*,
            )*
        }

        impl ::std::fmt::Display for $name {
            #[allow(unused_variables)]
            fn fmt(&self, f: &mut ::std::fmt::Formatter<'_>) -> ::std::fmt::Result {
                match self {
                    $(
                        Self::$var $($args2)* => write!(f, $desc),
                    )*
                }
            }
        }

        impl ::std::fmt::Debug for $name {
            fn fmt(&self, f: &mut ::std::fmt::Formatter<'_>) -> ::std::fmt::Result {
                match self {
                    $(
                        Self::$var $($args2)* => write!(f, $msg $(,$arg)*),
                    )*
                }
            }
        }

        #[cfg(feature = "std")]
        #[cfg_attr(docsrs, doc(cfg(feature = "single")))]
        impl ::std::error::Error for $name {}
    };
}

impl_error!(
    #[doc = "calculation error information occured during integration."]
    #[derive(Clone, Copy, PartialEq, Eq, Hash)]
    #[non_exhaustive]
    pub enum RuntimeError {
        InsufficientIteration {
            "number of iteration was insufficient",
            "The maximum number of subdivisions has been achieved.\n\
            If increasing the limit results in no further improvement, \
            check the integrand in order to determine the difficulties. \
            If the function contains the singular points in the interval, \
            you should specify the singular points by hand, or transform \
            the function to eliminate the singular points.",
        },
        RoundoffError {
            "cannot reach tolerance because of roundoff error",
            "Cannot reach tolerance because of roundoff error, \
            which prevents the given tolerance from being achieved.\n\
            It is assumed that the requested tolerance cannot be achieved, and \
            that the returned result is the bst which can be obtained.",
        },
        SubintervalTooSmall {
            "subinterval is too small to calculate the integral",
            "Subinterval was too small to calculate the integral.\n\
            Maybe you should specify the singular points, or transform the \
            function to eliminate the singular points.",
        },
        Divergent {
            "integral is divergent, or slowly convergent",
            "Integral is divergent, or slowly convergent.\n\
            Delta (estimation of absolute error) may be underestimated.",
        },
        NanValueEncountered {
            "integrand has returned a NAN value",
            "Integrand has returned a NAN value, so the algorithm cannot \
            continue the calculation.",
        },
    }
);

/// Store the estimation of integral and estimated absolute error (delta) of
/// calculation.
#[derive(Debug, Clone)]
pub struct IntegrationResult {
    pub(crate) estimate: f64,
    pub(crate) delta: f64,
    pub(crate) error: Option<RuntimeError>,
}

impl IntegrationResult {
    /// Create `IntegrationResult` object
    #[inline]
    pub fn new(estimate: f64, delta: f64, error: Option<RuntimeError>) -> Self {
        Self {
            estimate,
            delta,
            error,
        }
    }

    /// Get the estimation of integral.
    ///
    /// Returns `Err` if integration error has occured.
    #[inline]
    pub fn estimate(&self) -> Result<f64, RuntimeError> {
        match self.error {
            Some(error) => Err(error),
            None => Ok(self.estimate),
        }
    }

    /// Get the estimated absolute error of integral.
    ///
    /// Returns `Err` if integration error has occured.
    #[inline]
    pub fn delta(&self) -> Result<f64, RuntimeError> {
        match self.error {
            Some(error) => Err(error),
            None => Ok(self.delta),
        }
    }

    /// Get the estimation of integral and absolute error.
    ///
    /// Returns `Err` if integration error has occured.
    ///
    /// # Examples
    ///
    /// ```
    /// use gkquad::single::integral;
    /// use std::f64::NEG_INFINITY;
    ///
    /// let result = integral(|x: f64| x.exp(), NEG_INFINITY..0.0);
    /// let (estimate, delta) = result.estimate_delta().unwrap();
    /// ```
    #[inline]
    pub fn estimate_delta(&self) -> Result<(f64, f64), RuntimeError> {
        match self.error {
            Some(error) => Err(error),
            None => Ok((self.estimate, self.delta)),
        }
    }

    /// Get the estimation of integral even if integration error occured.
    ///
    /// # Notes
    ///
    /// If integration error occured during computation, integration function
    /// will return the current estimation of integral anyway.
    /// In this case the estimation does not achieve the given tolerance.
    ///
    /// You should use this function **only if** you don't care about the
    /// estimation error.
    #[inline]
    pub unsafe fn estimate_unchecked(&self) -> f64 {
        self.estimate
    }

    /// Get the estimated absolute error even if integration error has occured.
    ///
    /// # Notes
    ///
    /// If integration error occured during computation, integration function
    /// will return the current estimation of integral anyway.
    /// In this case the estimation does not achieve the given tolerance, and
    /// the error estimation (delta) may not be reliable.
    #[inline]
    pub unsafe fn delta_unchecked(&self) -> f64 {
        self.delta
    }

    /// Get the estimation of integral and absolute error even if integration
    /// error has occured.
    ///
    /// # Notes
    ///
    /// See notes of `estimate_unchecked` method and `delta_unchecked` method.
    #[inline]
    pub unsafe fn estimate_delta_unchecked(&self) -> (f64, f64) {
        (self.estimate, self.delta)
    }

    /// Return true if integration error has occured.
    #[inline]
    pub fn has_err(&self) -> bool {
        self.error.is_some()
    }

    /// Get integration error.
    ///
    /// Return `None` if any integration error has not occured.
    #[inline]
    pub fn err(&self) -> Option<RuntimeError> {
        self.error
    }
}
