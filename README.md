akima_spline
============

[`AkimaSplines`]: https://docs.rs/var_quantity/0.1.0/akima_spline/struct.AkimaSpline.html


Implementation of an Akima spline according to [Aki1970]
In addition to the original method, this implementation also allows for customized extrapolation at both spline sides.
The extrapolation behaviour is defined by the coefficients of a polynomial

Case 1: Extrapolation is set to "nothing"
================================================================================
In this case, the default quadratic extrapolation according to section 2.3 of [1]
is used to construct the spline and the spline throws an error when evaluated
outside its boundaries.

Case 2: Extrapolation is set to a vector
================================================================================
The container is given as a n-element collection of values:
    `extrapl = [p1, p2, p3, ..., pn]`
which are used to evaluate the following function:
    `y = p0 + p1*(x-x3) + p2*(x-x3)^2 + p3*(x-x3)^3 + ... + pn*(x-x3)^n`
The constant p0 is equal to y3 by definition (continuous spline). x3 and y3
are the first/last datapoint values respectively. If `extrapl` is an empty vector `Vec::new()`,
the extrapolation is therefore the constant `p0`.

[Aki1970]:  Akima, Hiroshi (1970). "A new method of interpolation and smooth curve fitting
based on local procedures". Journal of the ACM. 17: 589–602.

A spline models a x-y relationship via multiple polynomials.
The first step of evaluating it at a position `x` is to find the correct polynomial.
For this purpose, this implementation uses a binary search since it is pretty much always
the optimal algorithm, see for example the following two blog posts:
https://dirtyhandscoding.github.io/posts/performance-comparison-linear-search-vs-binary-search.html
https://pvk.ca/Blog/2012/07/03/binary-search-star-eliminates-star-branch-mispredictions/


Akima, Hiroshi: A new method of interpolation and smooth curve fitting based on
local procedures. Journal of the Association for Computing Machinery, Vol. 17,
No. 4, October 1970, pp.589-602.

# Serialization and deserialization

The [`AkimaSplines`] struct can be serialized / deserialized if the `serde`
feature is enabled.

TODO

# Documentation

The full API documentation is available at
[https://docs.rs/akima_spline/0.1.0/akima_spline/](https://docs.rs/akima_spline/0.1.0/akima_spline/).