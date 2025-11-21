#![cfg_attr(docsrs, doc = include_str!("../README.md"))]
#![cfg_attr(not(docsrs), doc = include_str!("../README_local.md"))]
#![deny(missing_docs)]

use horner::eval_polynomial;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/**
An Akima spline with optional extrapolation.

See the [module-level documentation](crate) for more.
 */
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Deserialize, Serialize))]
#[cfg_attr(feature = "serde", serde(try_from = "serde_impl::Alias"))]
#[cfg_attr(feature = "serde", serde(into = "serde_impl::Alias"))]
pub struct AkimaSpline {
    xs: Vec<f64>,
    ys: Vec<f64>,
    ps: Vec<[f64; 4]>,
    extrapl: Option<Vec<f64>>,
    extrapr: Option<Vec<f64>>,
}

impl AkimaSpline {
    /**
    Constructs an [`AkimaSpline`] out of valid input data.

    The input data is valid if the following conditions are fulfilled.
    - The lengths of `xs` and `ys` must be equal.
    - `xs` (and `ys`) must have at least five elements.
    - `xs` needs to be strictly monotonic increasing.
    See [`BuildError`] documentation.

    If coefficients for a left- or right-side extrapolation polynom are given
    (arguments `extrapl` and `extrapr`), they are interpreted as polynom
    coefficients of descending degree, where the maximum degree is the length
    of the coefficient vector. For example, the vector `[a, b, c]` is converted
    into the following polynom:

    `y = a(x-x0)³ + b(x-x0)² + c(x-x0) + d`,

    where x0 is the first / last value of `xs`. The coefficient `d` is set to
    the y-value of the first / last value of `ys` to make the transition between
    inter- and extrapolation "smooth" (values of inter- and extrapolation
    polynoms and their first derivatives are equal at the transition point).
    `d` is then pushed to the end of `extrapl` / `extrapr` and is therefore the
    last value returned from [`AkimaSpline::extrapl`] / [`AkimaSpline::extrapr`]
    (see examples).

    # Examples

    ## Valid / invalid input data
    ```
    use akima_spline::AkimaSpline;

    // Valid input data
    let xdata = vec![-2.0, -1.0, 0.0, 1.0, 2.0];
    let ydata = vec![0.0, 1.0, 0.0, 1.0, 0.0];
    assert!(AkimaSpline::new(xdata, ydata, None, None).is_ok());

    // Lengths of xs and ys unequal
    let xdata = vec![-2.0, -1.0, 0.0, 1.0, 2.0, 3.0];
    let ydata = vec![0.0, 1.0, 0.0, 1.0, 0.0];
    assert!(AkimaSpline::new(xdata, ydata, None, None).is_err());

    // Less than five datapoints
    let xdata = vec![-2.0, -1.0, 0.0, 1.0];
    let ydata = vec![0.0, 1.0, 0.0, 1.0];
    assert!(AkimaSpline::new(xdata, ydata, None, None).is_err());

    // xs not strictly monotonic increasing
    let xdata = vec![-2.0, -1.0, 0.0, 2.0, 1.0];
    let ydata = vec![0.0, 1.0, 0.0, 1.0, 0.0];
    assert!(AkimaSpline::new(xdata, ydata, None, None).is_err());
    ```

    ## Modification of extrapolation coefficient vector

    ```
    use akima_spline::AkimaSpline;

    let xdata = vec![-2.0, -1.0, 0.0, 1.0, 2.0];
    let ydata = vec![0.0, 1.0, 0.0, 1.0, 3.0];
    let extrapl = vec![1.0, 2.0];
    let extrapr = vec![0.5];

    let spline = AkimaSpline::new(xdata, ydata, Some(extrapl.clone()), Some(extrapr.clone())).expect("valid data");

    let extrapl_mod = spline.extrapl().expect("has left-side extrapolation polynom");
    let extrapr_mod = spline.extrapr().expect("has right-side extrapolation polynom");

    assert_eq!(extrapl_mod[0], extrapl[0]);
    assert_eq!(extrapl_mod[1], extrapl[1]);
    assert_eq!(extrapl_mod[2], 0.0); // First y-value

    assert_eq!(extrapr_mod[0], extrapr[0]);
    assert_eq!(extrapr_mod[1], 3.0); // Last y-value
    ```
     */
    pub fn new(
        xs: Vec<f64>,
        ys: Vec<f64>,
        mut extrapl: Option<Vec<f64>>,
        mut extrapr: Option<Vec<f64>>,
    ) -> Result<Self, BuildError> {
        // Sanity checks for the input data
        if xs.len() < 5 {
            return Err(BuildError::MinFivePointsNeeded);
        }
        if xs.len() != ys.len() {
            return Err(BuildError::InequalSliceLength {
                xs_len: xs.len(),
                ys_len: ys.len(),
            });
        }

        // Make sure the x-vector is strictly monotonic increasing
        if !(xs.windows(2).all(|w| w[0] < w[1])) {
            return Err(BuildError::XsNotStrictlyMonotonicIncreasing);
        }

        // Section 2.3 from Akima's paper: Use end point and two adjacent points
        // to calculate slopes "beyond" the datapoints. If no extrapolation
        // polynom is given, a quadratic extrapolation is used. If one is given,
        // the constant component k of the polynom is added to the polynom.

        // Left side
        let xs_l = [xs[2], xs[1], xs[0]];
        let ys_l = [ys[2], ys[1], ys[0]];
        let (xb, xa, yb, ya) = extrapolate(xs_l, ys_l, &mut extrapl);

        // Right side
        let xs_r = [xs[xs.len() - 3], xs[xs.len() - 2], xs[xs.len() - 1]];
        let ys_r = [ys[xs.len() - 3], ys[xs.len() - 2], ys[xs.len() - 1]];
        let (xy, xz, yy, yz) = extrapolate(xs_r, ys_r, &mut extrapr);

        // Extend the datapoint vectors xs and ys with the extrapolated points
        let mut xs_extrap: Vec<f64> = vec![0.0; xs.len() + 4];
        xs_extrap[0] = xa;
        xs_extrap[1] = xb;
        xs_extrap[xs.len() + 2] = xy;
        xs_extrap[xs.len() + 3] = xz;

        let mut ys_extrap: Vec<f64> = vec![0.0; xs.len() + 4];
        ys_extrap[0] = ya;
        ys_extrap[1] = yb;
        ys_extrap[xs.len() + 2] = yy;
        ys_extrap[xs.len() + 3] = yz;

        for ii in 0..xs.len() {
            xs_extrap[ii + 2] = xs[ii];
            ys_extrap[ii + 2] = ys[ii];
        }

        // Local segment slopes dy / dx
        let mut ms: Vec<f64> = vec![0.0; xs.len() + 3];

        let x_iter_0 = [xa, xb]
            .into_iter()
            .chain(xs.iter().cloned())
            .chain([xy, xz].into_iter());
        let x_iter_1 = x_iter_0.clone().skip(1);
        let dx_iter = x_iter_1.zip(x_iter_0).map(|(x1, x0)| x1 - x0);

        let y_iter_0 = [ya, yb]
            .into_iter()
            .chain(ys.iter().cloned())
            .chain([yy, yz].into_iter());
        let y_iter_1 = y_iter_0.clone().skip(1);
        let dy_iter = y_iter_1.zip(y_iter_0).map(|(y1, y0)| y1 - y0);

        for ((m, dx), dy) in ms.iter_mut().zip(dx_iter).zip(dy_iter) {
            *m = dy / dx;
        }

        // Weighted segment slope calculation
        let mut ts: Vec<f64> = vec![0.0; xs.len()];
        for (elem, ms_slice) in ts.iter_mut().zip(ms.windows(4)) {
            let m1 = ms_slice[0];
            let m2 = ms_slice[1];
            let m3 = ms_slice[2];
            let m4 = ms_slice[3];

            // As described in Akima's paper, p.591 (parentheses block), this is an arbitrary
            // convention to guarantee uniqueness of the solution
            if m1 == m2 && m3 == m4 {
                *elem = (m2 + m3) / 2.0;
            } else {
                let numer = (m4 - m3).abs() * m2 + (m2 - m1).abs() * m3;
                let denom = (m4 - m3).abs() + (m2 - m1).abs();
                *elem = numer / denom;
            }
        }

        /*
        If an extrapolation polynom is given, the weighted segment slope at
        the respective end point is equal to the first derivative of that
        polynom (which in turn is simply the first-degree term of the
        polynom). If the polynom does not have a first-degree term, use the
        previously calculated weighted slope
        */
        if let Some(ext) = extrapl.as_ref() {
            if ext.len() > 1 {
                ts[0] = ext[ext.len() - 2];
            }
        }
        if let Some(ext) = extrapr.as_ref() {
            if ext.len() > 1 {
                // No overflow risk here, since ts has xs.len() elements and
                // xs.len() is guaranteed to be at least 5.
                let l = ts.len() - 1;
                ts[l] = ext[ext.len() - 2];
            }
        }

        // Calculate the spline coefficients
        let mut ps: Vec<[f64; 4]> = Vec::with_capacity(xs.len() - 1);
        for ii in 0..(xs.len() - 1) {
            let x1 = xs[ii];
            let x2 = xs[ii + 1];

            let y1 = ys[ii];
            let y2 = ys[ii + 1];

            let t1 = ts[ii];
            let t2 = ts[ii + 1];

            let p0 = y1;
            let p1 = t1;
            let p2 = (3.0 * (y2 - y1) / (x2 - x1) - 2.0 * t1 - t2) / (x2 - x1);
            let p3 = (t1 + t2 - 2.0 * (y2 - y1) / (x2 - x1)) / (x2 - x1).powi(2);
            ps.push([p3, p2, p1, p0]);
        }

        // Now the spline can be constructed
        return Ok(AkimaSpline {
            xs,
            ys,
            ps,
            extrapl,
            extrapr,
        });
    }

    /**
    Evaluate the spline at position `x` and return either the corresponding
    `y`-value or `None`, if `x` is out of bounds.

    There are three possible scenarios:
    - `x` is between [`xmin`](AkimaSpline::xmin) and [`xmax`](AkimaSpline::xmax):
    The polynom of the corresponding interval is used to calculate the return
    value.
    - `x` is smaller than [`xmin`](AkimaSpline::xmin): If left-side
    extrapolation polynom coefficients have been provided (see
    [`AkimaSpline::extrapl`]), those are used to calculate the return value,
    otherwise [`None`] is returned.
    - `x` is larger than [`xmax`](AkimaSpline::xmax): If right-side
    extrapolation polynom coefficients have been provided (see
    [`AkimaSpline::extrapr`]), those are used to calculate the return value,
    otherwise [`None`] is returned.

    # Examples

    ```
    use akima_spline::AkimaSpline;
    use approx;

    let xdata = vec![-2.0, -1.0, 0.0, 1.0, 2.0];
    let ydata = vec![0.0, 1.0, 0.0, 1.0, 0.0];

    // No extrapolation for either left or right side
    let spline = AkimaSpline::new(xdata, ydata, None, None).unwrap();

    // Within bounds
    assert_eq!(spline.eval(0.5).unwrap(), 0.5);
    approx::assert_abs_diff_eq!(spline.eval(1.9).unwrap(), 0.19, epsilon = 1e-10);
    assert_eq!(spline.eval(2.0).unwrap(), 0.0);

    // Out of bounds
    assert!(spline.eval(2.5).is_none());
    assert!(spline.eval(-3.6).is_none());
    ```
     */
    pub fn eval(&self, x: f64) -> Option<f64> {
        match self.eval_priv(x) {
            EvalResult::Value(v) => Some(v),
            EvalResult::OutOfBoundsLeft => None,
            EvalResult::OutOfBoundsRight => None,
        }
    }

    /**
    Like [`AkimaSpline::eval`], but returns the `y`-value of the leftmost
    datapoint if `x` is smaller than [`xmin`](AkimaSpline::xmin) and no
    left-side extrapolation coefficients have been provided. Likewise, the
    function returns the `y`-value of the rightmost datapoint if `x` is larger
    than [`xmax`](AkimaSpline::xmax) and no right-side extrapolation
    coefficients have been provided.

    # Examples

    Same example as used in [`AkimaSpline::eval`]:
    ```
    use akima_spline::AkimaSpline;
    use approx;

    let xdata = vec![-2.0, -1.0, 0.0, 1.0, 2.0];
    let ydata = vec![0.0, 1.0, 0.0, 1.0, 0.0];

    // No extrapolation for either left or right side
    let spline = AkimaSpline::new(xdata, ydata, None, None).unwrap();

    // Within bounds
    assert_eq!(spline.eval_infallible(0.5), 0.5);
    approx::assert_abs_diff_eq!(spline.eval_infallible(1.9), 0.19, epsilon = 1e-10);
    assert_eq!(spline.eval_infallible(2.0), 0.0);

    // Out of bounds
    assert_eq!(spline.eval_infallible(2.5), 0.0); // y-value of rightmost datapoint
    assert_eq!(spline.eval_infallible(-3.6), 0.0); // y-value of leftmost datapoint
    ```
     */
    pub fn eval_infallible(&self, x: f64) -> f64 {
        match self.eval_priv(x) {
            EvalResult::Value(v) => v,
            EvalResult::OutOfBoundsLeft => *unsafe { self.ys.get_unchecked(0) },
            EvalResult::OutOfBoundsRight => *unsafe { self.ys.get_unchecked(self.ys.len() - 1) },
        }
    }

    fn eval_priv(&self, x: f64) -> EvalResult {
        let idx = self.xs.partition_point(|&val| val <= x);

        // Extrapolation to the "left"
        if idx == 0 {
            match &self.extrapl {
                None => return EvalResult::OutOfBoundsLeft,
                // Calculate the extrapolation values with the extrapolation vector
                Some(extrap_vec) => {
                    // SAFETY: xs is guaranteed to contain at least one value by the constructor
                    return EvalResult::Value(unsafe {
                        eval_polynomial(x - self.xs.get_unchecked(0), &extrap_vec)
                            .unwrap_unchecked()
                    });
                }
            }

        // Extrapolation to the "right"
        } else if idx == self.xs.len() {
            let last_x_val = *self.xs.last().unwrap();
            if x > last_x_val {
                match &self.extrapr {
                    // Quadratic extrapolation of x4 and x5
                    None => return EvalResult::OutOfBoundsRight,
                    // Calculate the extrapolation values with the extrapolation vector
                    Some(extrap_vec) => {
                        // SAFETY: xs is guaranteed to contain at least one value by the constructor
                        return EvalResult::Value(unsafe {
                            eval_polynomial(x - last_x_val, &extrap_vec).unwrap_unchecked()
                        });
                    }
                }
            // x is the last datapoint -> simply return the corresponding y-value
            } else {
                // SAFETY: ys is guaranteed to contain at least one value by the constructor
                return EvalResult::Value(unsafe { *self.ys.get_unchecked(self.ys.len() - 1) });
            }
        }

        // If no extrapolation occurs, evaluate the found spline segment
        return EvalResult::Value(
            // Safety: xs and ps are guaranteed to contain the value idx-1
            unsafe {
                eval_polynomial(
                    x - self.xs.get_unchecked(idx - 1),
                    self.ps.get_unchecked(idx - 1),
                )
                .unwrap_unchecked()
            },
        );
    }

    /**
    Returns the `diff_degree`th derivative of the spline at position `x` if `x`
    is within bounds, otherwise returns [`None`]. See [`AkimaSpline::eval`] for
    the definition of "within bounds".

    # Examples

    ```
    use akima_spline::AkimaSpline;

    let xdata = vec![-2.0, -1.0, 0.0, 1.0, 2.0];
    let ydata = vec![0.0, 1.0, 0.0, 1.0, 0.0];
    let spline = AkimaSpline::new(xdata, ydata, None, None).unwrap();

    assert_eq!(spline.derivative(0.5, 1).unwrap(), 1.5);
    assert_eq!(spline.derivative(0.5, 2).unwrap(), 0.0);

    assert_eq!(spline.derivative(1.0, 1).unwrap(), 0.0);
    assert_eq!(spline.derivative(1.0, 2).unwrap(), -2.0);

    assert_eq!(spline.derivative(1.5, 1).unwrap(), -1.0);
    assert_eq!(spline.derivative(1.5, 2).unwrap(), -2.0);

    assert_eq!(spline.derivative(2.0, 1).unwrap(), -2.0);
    assert_eq!(spline.derivative(2.0, 2).unwrap(), -2.0);
    ```
     */
    pub fn derivative(&self, x: f64, diff_degree: usize) -> Option<f64> {
        let x_ref: f64;
        let coeffs: &[f64];
        let idx = self.xs.partition_point(|&val| val <= x);

        // Extrapolation to the "left"
        if idx == 0 {
            coeffs = self.extrapl.as_ref()?.as_slice();
            x_ref = self.xs[0];
        }
        // Extrapolation to the "right"
        else if idx == self.xs.len() {
            // SAFETY: xs is guaranteed to have at least one value
            let x_check = *unsafe { self.xs.get_unchecked(self.xs.len() - 1) };
            if x > x_check {
                coeffs = self.extrapr.as_ref()?;
                x_ref = x_check;
            } else {
                // SAFETY: ps has at least one element (this is assured in the constructor)
                coeffs = unsafe { self.ps.get_unchecked(self.ps.len() - 1) };
                x_ref = *unsafe { self.xs.get_unchecked(idx - 2) };
            }
        } else {
            // SAFETY: idx is in bounds
            coeffs = unsafe { self.ps.get_unchecked(idx - 1) }.as_slice();
            x_ref = *unsafe { self.xs.get_unchecked(idx - 1) };
        }

        return Some(derivative(coeffs, x - x_ref, diff_degree));
    }

    /**
    Returns the underlying `xs` vector.
     */
    pub fn xs(&self) -> &[f64] {
        return self.xs.as_slice();
    }

    /**
    Returns the underlying `ys` vector.
     */
    pub fn ys(&self) -> &[f64] {
        return self.ys.as_slice();
    }

    /**
    Returns the calculated polynomial coefficients `ps`.

    The `n` points defined by `xs` and `ys` (see [`AkimaSpline`] docstring)
    result in `n - 1` interpolation polynoms which cover the intervals between
    them. These polynoms are of third degree (`ax³ + bx² + cx + d`).
    The coefficients correspond to the arrays `[a, b, c, d]` contained within
    the returned slice. See [`AkimaSpline`] docstring for more.

    # Examples

    ```
    use akima_spline::AkimaSpline;

    let xdata = vec![-2.0, -1.0, 0.0, 1.0, 2.0];
    let ydata = vec![0.0, 1.0, 0.0, 1.0, 0.0];
    let num_intervals = xdata.len() - 1;

    let spline = AkimaSpline::new(xdata, ydata, None, None).expect("valid data");

    // 5 points -> 4 intervals
    assert_eq!(spline.ps().len(), num_intervals);

    // Manually evaluate one of the polynoms and compare its value to that
    // returned by "spline".
    let i = spline.ps()[2]; // Interval from x = 0 to x = 1

    // Evaluate manually for x = 0.5
    let x: f64 = 0.5;
    let pol_eval = i[0] * x.powi(3) + i[1] * x.powi(2) + i[2] * x.powi(1) + i[3];
    assert_eq!(spline.eval_infallible(x), pol_eval);
    ```
     */
    pub fn ps(&self) -> &[[f64; 4]] {
        return self.ps.as_slice();
    }

    /**
    Returns the modified left-side (`x < xs.max()`) extrapolation coefficients,
    if given during construction. See [`AkimaSpline::new`] documentation.
     */
    pub fn extrapl(&self) -> Option<&[f64]> {
        return self.extrapl.as_ref().map(|val| val.as_slice());
    }

    /**
    Returns the modified right-side (`x > xs.max()`) extrapolation coefficients,
    if given during construction. See [`AkimaSpline::new`] documentation.
     */
    pub fn extrapr(&self) -> Option<&[f64]> {
        return self.extrapr.as_ref().map(|val| val.as_slice());
    }

    /**
    Returns the `xmin` value of the first data point / smallest value of `xs`.
    When trying the evaluate the spline for a smaller `x`-value, the left-side
    interpolation [`AkimaSpline::extrapl`] is used if available.
     */
    pub fn xmin(&self) -> f64 {
        // SAFETY: xs has at least one value
        return unsafe { *self.xs.get_unchecked(0) };
    }

    /**
    Returns the `xmax` value of the last data point / largest value of `xs`.
    When trying the evaluate the spline for a larger `x`-value, the right-side
    interpolation [`AkimaSpline::extrapr`] is used if available.
     */
    pub fn xmax(&self) -> f64 {
        // SAFETY: xs has at least one value
        return unsafe { *self.xs.get_unchecked(self.xs.len() - 1) };
    }
}

/// Helper enum used in the internal eval_priv implementation
enum EvalResult {
    Value(f64),
    OutOfBoundsLeft,
    OutOfBoundsRight,
}

fn extrapolate(xs: [f64; 3], ys: [f64; 3], extrap: &mut Option<Vec<f64>>) -> (f64, f64, f64, f64) {
    let [x1, x2, x3] = xs;
    let [y1, y2, y3] = ys;

    // Extrapolate the x-coordinates of two additional datapoints with (8) from Akima's paper
    let x4 = x3 - x1 + x2;
    let x5 = x3 - x1 + x3;

    let y4: f64;
    let y5: f64;

    match extrap {
        // Quadratic extrapolation of x4 and x5
        None => {
            y4 = y3 + 2.0 * (x4 - x3) * (y3 - y2) / (x3 - x2) - (x4 - x3) * (y2 - y1) / (x2 - x1);
            y5 = y4 + 2.0 * (x5 - x4) * (y4 - y3) / (x4 - x3) - (x5 - x4) * (y3 - y2) / (x3 - x2);
        }
        // Calculate the extrapolation values with the extrapolation polynom
        Some(extrap_vec) => {
            extrap_vec.push(y3);

            y4 = eval_polynomial(x4 - x3, extrap_vec.as_slice())
                .expect("guaranteed to have at least one value");
            y5 = eval_polynomial(x5 - x3, extrap_vec.as_slice())
                .expect("guaranteed to have at least one value");
        }
    }
    return (x4, x5, y4, y5);
}

/**
Calculates the `diff_degree` derivative of a polynom with the given `coeffs`
for `x`.

Example: Let coeff be [30, 5, 1].
The corresponding polynom would be 30x^2 + 5x^1 + 1x^0
Differentiating this polynom by x once yields:
30*2*x^1 + 5*1*x + 1*0*x
The corresponding coefficient vector is therefore coeff_diff = [60, 5, 0]
The length of the coefficient vector is num_elems = 3, therefore the elements are calculated as:
coeff_diff[0] = coeff[0]*(3-1)
coeff_diff[1] = coeff[1]*(3-2)
coeff_diff[2] = coeff[2]*(3-3)
*/
fn derivative(coeffs: &[f64], x: f64, diff_degree: usize) -> f64 {
    let num_coeffs = coeffs.len();
    if diff_degree > num_coeffs {
        return 0.0;
    }

    // All coefficients beyond `num_coeffs - diff_degree` will be zero and can
    // therefore be ignored
    let coeffs = &coeffs[..(num_coeffs - diff_degree)];

    return coeffs
        .iter()
        .enumerate()
        .map(|(idx, coeff)| {
            let exponent = num_coeffs - idx - 1;

            // Differentiate n times
            let mut diff_coeff = coeff.clone();
            for k in 0..diff_degree {
                if exponent > k {
                    let factor: f64 = (exponent - k) as f64;
                    diff_coeff *= factor;
                } else {
                    diff_coeff = 0.0;
                    break;
                }
            }
            diff_coeff * x.powi((exponent - diff_degree) as i32)
        })
        .sum();
}

/**
Errors which can happen during [`AkimaSpline`] construction due to bad input
data.
 */
#[derive(Debug)]
pub enum BuildError {
    /**
    Less than five points where given, making construction of an [`AkimaSpline`]
    impossible.
     */
    MinFivePointsNeeded,
    /**
    The given vectors `xs` and `ys` do not have the same length
     */
    InequalSliceLength {
        /// Length of the xs-slice.
        xs_len: usize,
        /// Length of the ys-slice.
        ys_len: usize,
    },
    /**
    The `xs`-values are not strictly monotonic increasing.
     */
    XsNotStrictlyMonotonicIncreasing,
}

impl std::fmt::Display for BuildError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BuildError::MinFivePointsNeeded => write!(
                f,
                "at least five datapoints are needed to create an Akima spline"
            ),
            BuildError::InequalSliceLength { xs_len, ys_len } => write!(
                f,
                "xs has {} elements, while ys has {} elements (must be equal)",
                xs_len, ys_len
            ),
            BuildError::XsNotStrictlyMonotonicIncreasing => write!(
                f,
                "values in xs-slice must be strictly monotonic increasing"
            ),
        }
    }
}

impl std::error::Error for BuildError {}

#[cfg(feature = "serde")]
mod serde_impl {
    use super::*;

    #[derive(Deserialize, Serialize)]
    pub(super) struct Alias {
        pub(super) xs: Vec<f64>,
        pub(super) ys: Vec<f64>,
        #[serde(default)]
        pub(super) extrapl: Option<Vec<f64>>,
        #[serde(default)]
        pub(super) extrapr: Option<Vec<f64>>,
    }

    impl From<AkimaSpline> for Alias {
        fn from(value: AkimaSpline) -> Self {
            let extrapl = value.extrapl.map(|mut vec| {
                let _ = vec.pop(); // Remove the last element of the extrapolation vector
                vec
            });
            let extrapr = value.extrapr.map(|mut vec| {
                let _ = vec.pop(); // Remove the last element of the extrapolation vector
                vec
            });
            return Self {
                xs: value.xs,
                ys: value.ys,
                extrapl,
                extrapr,
            };
        }
    }

    impl TryFrom<Alias> for AkimaSpline {
        type Error = BuildError;

        fn try_from(value: Alias) -> Result<Self, Self::Error> {
            return Self::new(value.xs, value.ys, value.extrapl, value.extrapr);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_derivative() {
        {
            let coeffs = [30.0, 5.0, 1.0];
            // n = 0: 30x² + 5x + 1
            // n = 1: 60x + 5
            // n = 2: 60
            // n = 3: 0
            // n = 3: 0

            let x = 1.0;
            assert_eq!(derivative(coeffs.as_slice(), x, 0), 36.0);
            assert_eq!(derivative(coeffs.as_slice(), x, 1), 65.0);
            assert_eq!(derivative(coeffs.as_slice(), x, 2), 60.0);
            assert_eq!(derivative(coeffs.as_slice(), x, 3), 0.0);
            assert_eq!(derivative(coeffs.as_slice(), x, 4), 0.0);

            let x = 2.0;
            assert_eq!(derivative(coeffs.as_slice(), x, 0), 131.0);
            assert_eq!(derivative(coeffs.as_slice(), x, 1), 125.0);
            assert_eq!(derivative(coeffs.as_slice(), x, 2), 60.0);
            assert_eq!(derivative(coeffs.as_slice(), x, 3), 0.0);
            assert_eq!(derivative(coeffs.as_slice(), x, 4), 0.0);
        }

        {
            let coeffs = [2.0, 30.0, 5.0, 1.0];
            // n = 0: 2x³ + 30x² + 5x + 1
            // n = 1: 6x² + 60x + 5
            // n = 2: 12x + 60
            // n = 3: 12
            // n = 3: 0

            let x = 1.0;
            assert_eq!(derivative(coeffs.as_slice(), x, 0), 38.0);
            assert_eq!(derivative(coeffs.as_slice(), x, 1), 71.0);
            assert_eq!(derivative(coeffs.as_slice(), x, 2), 72.0);
            assert_eq!(derivative(coeffs.as_slice(), x, 3), 12.0);
            assert_eq!(derivative(coeffs.as_slice(), x, 4), 0.0);

            let x = 2.0;
            assert_eq!(derivative(coeffs.as_slice(), x, 0), 147.0);
            assert_eq!(derivative(coeffs.as_slice(), x, 1), 149.0);
            assert_eq!(derivative(coeffs.as_slice(), x, 2), 84.0);
            assert_eq!(derivative(coeffs.as_slice(), x, 3), 12.0);
            assert_eq!(derivative(coeffs.as_slice(), x, 4), 0.0);
        }
    }
}
