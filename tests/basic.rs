use akima_spline::*;

#[test]
fn test_interpolation() {
    // Create data that is symmetric about zero
    let xdata = vec![-2.0, -1.0, 0.0, 1.0, 2.0];
    let ydata = vec![0.0, 1.0, 0.0, 1.0, 0.0];

    // Test that the spline is symmetric about zero
    let spline = AkimaSpline::new(xdata, ydata, None, None).unwrap();

    // Request some values and assert that the spline is symmetrical
    let result = 0.5;
    approx::assert_abs_diff_eq!(spline.eval(0.5).unwrap(), result, epsilon = 0.0001);
    approx::assert_abs_diff_eq!(spline.eval(-0.5).unwrap(), result, epsilon = 0.0001);

    let result = 0.896;
    approx::assert_abs_diff_eq!(spline.eval(0.8).unwrap(), result, epsilon = 0.0001);
    approx::assert_abs_diff_eq!(spline.eval(-0.8).unwrap(), result, epsilon = 0.0001);

    // Interpolate a sine curve with binary search
    let xdata = vec![-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0];
    let ydata: Vec<f64> = xdata.clone().into_iter().map(|x: f64| x.sin()).collect();

    let spline = AkimaSpline::new(xdata, ydata, None, None).unwrap();
    approx::assert_abs_diff_eq!(spline.eval(0.5).unwrap(), (0.5f64).sin(), epsilon = 0.0001);
    approx::assert_abs_diff_eq!(spline.eval(0.25).unwrap(), 0.239712769, epsilon = 0.0001);

    // Previously found bug: Request value directly on the data limits
    approx::assert_abs_diff_eq!(spline.eval(-2.0).unwrap(), -0.909297, epsilon = 0.0001);
    approx::assert_abs_diff_eq!(spline.eval(2.0).unwrap(), 0.909297, epsilon = 0.0001);

    // Check that no result is returned when a value outside of the xdata-range is requested, since we specified no extrapolation rule on either side
    assert!(spline.eval(-2.1).is_none()); // Left side
    assert!(spline.eval(2.1).is_none()); // Right side
}

#[test]
fn test_extrapolation_test() {
    {
        // Check left-hand side linear extrapolation (right-hand side is not extrapolated)
        let xdata = vec![0.0, 0.5, 1.0, 1.5, 2.0];
        let ydata: Vec<f64> = xdata.clone().into_iter().map(|x: f64| x.sin()).collect();
        let extrapl = Some(vec![1.0]);

        let spline = AkimaSpline::new(xdata, ydata, extrapl, None).unwrap();
        assert!(spline.eval(2.1).is_none()); // Right-hand side still not extrapolated
        approx::assert_abs_diff_eq!(spline.eval(-1.0).unwrap(), -1.0, epsilon = 0.0001);
        // Slope of 1, starting at (0,0)
    }
    {
        // Right-hand side quadratic extrapolation, left-hand side linear interpolation
        let xdata = vec![-0.5, 0.0, 0.5, 1.0, 1.5];
        let ydata: Vec<f64> = xdata.clone().into_iter().map(|x: f64| x.sin()).collect();
        let extrapl = Some(vec![1.0]); // Lin polynom
        let extrapr = Some(vec![1.0, 0.0]); // Quadratic polynom: f(x) = 1.0*x^2 + 0.0*x^0
        let spline = AkimaSpline::new(xdata, ydata, extrapl, extrapr).unwrap();
        approx::assert_abs_diff_eq!(spline.eval(-1.0).unwrap(), -0.9794, epsilon = 0.0001);
        approx::assert_abs_diff_eq!(
            spline.eval(2.0).unwrap(),
            0.99749 + 0.5 * 0.5,
            epsilon = 0.0001
        );
        approx::assert_abs_diff_eq!(
            spline.eval(4.0).unwrap(),
            0.99749 + 2.5 * 2.5,
            epsilon = 0.0001
        );
    }
    {
        // Extrapolate with a constant (slope = 0)
        let xdata = vec![-0.5, 0.0, 0.5, 1.0, 1.5];
        let ydata: Vec<f64> = xdata.clone(); // Linear curve
        let extrapl = Some(Vec::new()); // Constant
        let extrapr = Some(Vec::new()); // Constant
        let spline = AkimaSpline::new(xdata, ydata, extrapl, extrapr).unwrap();
        assert_eq!(spline.eval(-1.0).unwrap(), -0.5);
        assert_eq!(spline.eval(2.0).unwrap(), 1.5);
        assert_eq!(spline.eval(10.0).unwrap(), 1.5);
    }
}

#[test]
fn test_derivative_test() {
    let xdata = vec![-0.5, 0.0, 0.5, 1.0, 1.5];
    let ydata: Vec<f64> = xdata.clone().into_iter().map(|x: f64| x.sin()).collect();
    let extrapl = Some(vec![1.0]); // Lin polynom
    let extrapr = Some(vec![1.0, 0.0]); // Quadratic polynom: f(x) = 1.0*x^2 + 0.0*x^0
    let spline = AkimaSpline::new(xdata, ydata, extrapl, extrapr).unwrap();

    // Derivative from the spline itself
    approx::assert_abs_diff_eq!(
        spline.derivative(-0.5, 1).unwrap(),
        0.9794255,
        epsilon = 0.0001
    );

    // Show continuity of the first derivative between two spline pieces
    approx::assert_abs_diff_eq!(
        spline.derivative(-1e-8, 1).unwrap(),
        0.95885,
        epsilon = 0.0001
    );
    approx::assert_abs_diff_eq!(
        spline.derivative(0.0, 1).unwrap(),
        0.95885,
        epsilon = 0.0001
    );
    approx::assert_abs_diff_eq!(
        spline.derivative(1e-8, 1).unwrap(),
        0.95885,
        epsilon = 0.0001
    );

    // First derivative on the left extrapolation side (which is known to be 1.0 from extrapl for all x <= -0.5)
    assert_eq!(spline.derivative(-0.51, 1).unwrap(), 1.0);
    assert_eq!(spline.derivative(-1.0, 1).unwrap(), 1.0);

    // Second derivative
    approx::assert_abs_diff_eq!(
        spline.derivative(-0.5, 2).unwrap(),
        -0.1645956,
        epsilon = 0.0001
    );

    // Third derivative
    approx::assert_abs_diff_eq!(
        spline.derivative(-0.5, 3).unwrap(),
        0.493787,
        epsilon = 0.0001
    );

    // First derivative on the right extrapolation side: f(x) = 1.0*x^2 + 0.0*x^0 => f'(x) = 2.0*(2.0 - 1.5) (the 1.5 is the last x-datapoint of the spline)
    assert_eq!(spline.derivative(2.0, 1).unwrap(), 1.0);
    assert_eq!(spline.derivative(3.0, 1).unwrap(), 3.0); // f'(x) = 2.0*(3.0 - 1.5)

    // Second derivative
    assert_eq!(spline.derivative(2.0, 2).unwrap(), 2.0);

    // Almost the same spline, but the slope of the extrapolation should be 2 this time. On the right side, no interpolation is used, hence the derivative is None.
    let xdata = vec![-0.5, 0.0, 0.5, 1.0, 1.5];
    let ydata: Vec<f64> = xdata.clone().into_iter().map(|x: f64| x.sin()).collect();
    let extrapl = Some(vec![2.0]); // Lin polynom
    let extrapr = None; // Quadratic polynom: f(x) = 1.0*x^2 + 0.0*x^0
    let spline = AkimaSpline::new(xdata, ydata, extrapl, extrapr).unwrap();

    // First derivative on the left extrapolation side (which is known to be 2.0 from extrapl for all x < -0.5)
    assert_eq!(spline.derivative(-0.6, 1).unwrap(), 2.0);

    // First derivative on the right extrapolation side is known to be None
    assert!(spline.derivative(2.0, 1).is_none());
}

// Now check that the spline constructor panics if given wrong inputs
#[test]
#[should_panic]
fn too_few_values() {
    let xdata = vec![0.0, 0.5, 1.0, 1.5]; // Too few input values, 5 are needed at minimum
    let ydata: Vec<f64> = xdata.clone().into_iter().map(|x: f64| x.sin()).collect();
    AkimaSpline::new(xdata, ydata, None, None).unwrap();
}

#[test]
#[should_panic]
fn ineq_length() {
    let xdata = vec![0.0, 0.5, 1.0, 1.5, 2.0]; // 5 values
    let ydata: Vec<f64> = vec![0.0, 0.5, 1.0, 1.5, 2.0, 2.5]; // 6 values
    AkimaSpline::new(xdata, ydata, None, None).unwrap();
}
