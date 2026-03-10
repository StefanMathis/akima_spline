This crate implements a 1d Akima spline as described in:
> Akima, Hiroshi: A new method of interpolation and smooth curve fitting based on
local procedures. Journal of the Association for Computing Machinery, Vol. 17,
No. 4, October 1970, pp.589-602.

A spline can be constructed by providing x/y data points in the form of two
vectors `xs` and `ys` holding the respective coordinates:

```rust
use akima_spline::AkimaSpline;
use approx;

// Coordinates used in the plot above
let xs = vec![-2.0, -1.0, -0.8, 0.25, 1.0, 2.0];
let ys = vec![0.0, 1.0, 2.0, -1.5, 1.0, 3.0];
let spline = AkimaSpline::new(xs, ys, None, None).expect("valid input data");

// Data inside the spline
approx::assert_abs_diff_eq!(spline.eval(0.2).expect("0.5 is within bounds"), -1.582, epsilon=1e-3);

// Outside the defined points:
assert!(spline.eval(2.2).is_none());

// But it is also possible to evaluate the spline "infallible" via a simple 
// flat interpolation using the last known datapoint:
assert_eq!(spline.eval_infallible(2.2), 3.0);
```

[`AkimaSpline`] also offers the option of providing left-side (`x < xs[0]`)
and right-side (`x > xs[xs.len() - 1]`) extrapolation polynoms. To make sure the
transition between extrapolation and spline is "smooth" (first derivatives equal
at the transition point), the polynoms within the spline are adjusted. This can
be clearly seen when comparing splines made from the same datapoints with and
without extrapolation: