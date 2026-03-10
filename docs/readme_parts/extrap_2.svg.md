The extrapolation polynoms can be of any degree `n` defined by the length of the
given vector with the first value being the coefficient of the `x^n` term. The
constant part of the polynom is omitted, since it is inferred from the given
datapoints.

```rust
use akima_spline::AkimaSpline;

let xs = vec![-2.0, -1.0, -0.8, 0.25, 1.0, 2.0];
let ys = vec![0.0, 1.0, 2.0, -1.5, 1.0, 3.0];

// Polynom 2(x-x0) + k, where k is 0 (first ys-value) and x0 is -2 (first xs-value)
let extrapl = Some(vec![2.0]);

// Polynom 3(x-x0)² - 2(x-x0) + k, where k is 3 (last ys-value) and x0 is 2 (last xs-value)
let extrapr = Some(vec![3.0, -2.0]); 

let spline = AkimaSpline::new(xs, ys, extrapl, extrapr).expect("valid input data");

// Left-side extrapolation
assert_eq!(spline.eval(-2.5).expect("covered by extrapolation polynom"), -1.0);

// Right-side extrapolation
approx::assert_abs_diff_eq!(spline.eval(2.2).expect("covered by extrapolation polynom"), 2.72, epsilon=1e-3);
```

# Derivatives

The spline can be differentiated by an arbitrary degree at any position `x`
using the [`derivative`] method.