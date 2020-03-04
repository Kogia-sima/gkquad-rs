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

        impl ::core::fmt::Display for $name {
            #[allow(unused_variables)]
            fn fmt(&self, f: &mut ::core::fmt::Formatter<'_>) -> ::core::fmt::Result {
                match self {
                    $(
                        Self::$var $($args2)* => write!(f, $desc),
                    )*
                }
            }
        }

        impl ::core::fmt::Debug for $name {
            fn fmt(&self, f: &mut ::core::fmt::Formatter<'_>) -> ::core::fmt::Result {
                match self {
                    $(
                        Self::$var $($args2)* => write!(f, $msg $(,$arg)*),
                    )*
                }
            }
        }

        #[cfg(feature = "std")]
        #[cfg_attr(docsrs, doc(cfg(feature = "std")))]
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
            If the function contains the singular points in the range, \
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
        SubrangeTooSmall {
            "subrange is too small to calculate the integral",
            "Subrange was too small to calculate the integral.\n\
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
