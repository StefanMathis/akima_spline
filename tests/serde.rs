use akima_spline::*;
use serde_yaml;

#[test]
fn test_serialize_and_deserialize_no_extrap() {
    {
        // Create data that is symmetric about zero
        let xdata = vec![-2.0, -1.0, 0.0, 1.0, 2.0];
        let ydata = vec![0.0, 1.0, 0.0, 1.0, 0.0];

        let spline = AkimaSpline::new(xdata, ydata, None, None).unwrap();

        let serialized = serde_yaml::to_string(&spline).unwrap();
        let deserialized: AkimaSpline = serde_yaml::from_str(&serialized).unwrap();
        assert_eq!(deserialized.eval(1.0), spline.eval(1.0));
        assert_eq!(deserialized.eval(-1.0), spline.eval(-1.0));
    }
    {
        // Sinusoidal curve
        let xdata = vec![0.0, 0.5, 1.0, 1.5, 2.0];
        let ydata: Vec<f64> = xdata.clone().into_iter().map(|x: f64| x.sin()).collect();
        let spline = AkimaSpline::new(xdata, ydata, None, None).unwrap();

        let serialized = serde_yaml::to_string(&spline).unwrap();
        let deserialized: AkimaSpline = serde_yaml::from_str(&serialized).unwrap();
        assert_eq!(deserialized.eval(0.25), spline.eval(0.25));
        assert_eq!(deserialized.eval(1.25), spline.eval(1.25));
    }
    {
        // Example from magnetic permeability
        let xdata = vec![
            0.6066441589356562,
            0.7683640964032328,
            0.8340098434493579,
            0.8876048324612676,
            0.9309554970819128,
            0.9677761672163732,
            0.9998276617824616,
            1.0286465517935117,
            1.0571552831541977,
            1.0821887085641997,
            1.0997841549055147,
        ];
        let ydata = vec![
            8045.868049369589,
            7643.059002307005,
            7374.266065125728,
            7063.334829923217,
            6734.825870870863,
            6417.765034125603,
            6120.2890235650175,
            5846.935121870383,
            5608.382964323428,
            5382.365072694591,
            5148.120134922217,
        ];
        let spline = AkimaSpline::new(xdata, ydata, None, None).unwrap();

        let serialized = serde_yaml::to_string(&spline).unwrap();
        let deserialized: AkimaSpline = serde_yaml::from_str(&serialized).unwrap();

        assert_eq!(deserialized.eval(0.8), spline.eval(0.8));
        assert_eq!(deserialized.eval(1.0), spline.eval(1.0));
    }
}

#[test]
fn test_serialize_and_deserialize_extrap() {
    {
        // Check left-hand side linear extrapolation (right-hand side is not extrapolated)
        let xdata = vec![0.0, 0.5, 1.0, 1.5, 2.0];
        let ydata: Vec<f64> = xdata.clone().into_iter().map(|x: f64| x.sin()).collect();
        let extrapl = Some(vec![1.0]);

        let spline = AkimaSpline::new(xdata, ydata, extrapl, None).unwrap();

        let serialized = serde_yaml::to_string(&spline).unwrap();
        let deserialized: AkimaSpline = serde_yaml::from_str(&serialized).unwrap();

        assert_eq!(deserialized.eval(0.25), spline.eval(0.25));
        assert_eq!(deserialized.eval(-1.0), spline.eval(-1.0));
    }
    {
        // Check left-hand side linear extrapolation (right-hand side is not extrapolated)
        let xdata = vec![0.0, 0.5, 1.0, 1.5, 2.0];
        let ydata: Vec<f64> = xdata.clone().into_iter().map(|x: f64| x.sin()).collect();
        let extrapr = Some(vec![1.0]);

        let spline = AkimaSpline::new(xdata, ydata, None, extrapr).unwrap();

        let serialized = serde_yaml::to_string(&spline).unwrap();
        let deserialized: AkimaSpline = serde_yaml::from_str(&serialized).unwrap();

        assert_eq!(deserialized.eval(0.25), spline.eval(0.25));
        assert_eq!(deserialized.eval(5.0), spline.eval(5.0));
    }
}
