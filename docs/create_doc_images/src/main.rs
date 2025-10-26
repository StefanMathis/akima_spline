use akima_spline::AkimaSpline;
use plotters::{prelude::*, style::full_palette::GREEN_800};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Raw dataponts
    let xs = vec![-2.0, -1.0, -0.8, 0.25, 1.0, 2.0];
    let ys = vec![0.0, 1.0, 2.0, -1.5, 1.0, 3.0];

    // Extrapolation

    let spline_no_extrap =
        AkimaSpline::new(xs.clone(), ys.clone(), None, None).expect("valid input data");

    draw_no_extrap(&spline_no_extrap, "no_extrap")?;

    let extrapl = vec![2.0];
    let extrapr = vec![3.0, -2.0];
    let spline_extrapl =
        AkimaSpline::new(xs.clone(), ys.clone(), Some(extrapl), None).expect("valid input data");
    let spline_extrapr =
        AkimaSpline::new(xs.clone(), ys.clone(), None, Some(extrapr)).expect("valid input data");

    draw_extrap(
        &spline_no_extrap,
        &spline_extrapl,
        &spline_extrapr,
        "extrap_1",
        "left-side extrapolation (y = 2x + k)",
        "right-side extrapolation (y = 3x² - 2x + k)",
    )?;

    let extrapl = vec![];
    let extrapr = vec![-2.0, 0.0, 3.0, 0.0, 0.0 - 6.0];
    let spline_extrapl =
        AkimaSpline::new(xs.clone(), ys.clone(), Some(extrapl), None).expect("valid input data");
    let spline_extrapr =
        AkimaSpline::new(xs.clone(), ys.clone(), None, Some(extrapr)).expect("valid input data");

    draw_extrap(
        &spline_no_extrap,
        &spline_extrapl,
        &spline_extrapr,
        "extrap_2",
        "left-side extrapolation (y = k)",
        "right-side extrapolation (y = -2x⁶ + 3x⁴ -6x + k)",
    )?;

    return Ok(());
}

fn draw_no_extrap(
    spline_no_extrap: &AkimaSpline,
    filename: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    // General config
    // =========================================================================
    let size = (600, 400);
    let font_size_labels = 18;
    let font_size_ticks = 16;
    let font_size_legend = 18;

    // Calculate data
    // =========================================================================

    let resolution = 1e-3;
    let xmin = -3.0;
    let xmax = 3.0;
    let capacity = ((xmax - xmin) / resolution) as usize;

    let mut x_no_extrap = Vec::with_capacity(capacity);
    let mut y_no_extrap = Vec::with_capacity(capacity);

    let mut x = xmin;
    while x <= xmax {
        if let Some(y) = spline_no_extrap.eval(x) {
            x_no_extrap.push(x);
            y_no_extrap.push(y);
        }

        x += resolution;
    }

    let file_path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(&format!("../{filename}.svg"));
    let root = SVGBackend::new(&file_path, size).into_drawing_area();

    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(xmin..xmax, -2.0..4.0)?;

    chart
        .configure_mesh()
        .axis_desc_style(("sans-serif", font_size_labels))
        .label_style(("sans-serif", font_size_ticks))
        .draw()?;

    // Interpolated data
    chart
        .draw_series(LineSeries::new(
            x_no_extrap.iter().cloned().zip(y_no_extrap.iter().cloned()),
            &BLUE,
        ))?
        .label("spline interpolation")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    // Original data points
    chart
        .draw_series(
            spline_no_extrap
                .xs()
                .iter()
                .cloned()
                .zip(spline_no_extrap.ys().iter().cloned())
                .map(|(x, y)| Circle::new((x, y), 5.0, BLACK.filled())),
        )?
        .label("original data")
        .legend(|(x, y)| Circle::new((x + 10, y), 5.0, BLACK.filled()));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8)) // semi-transparent background
        .border_style(&BLACK)
        .label_font(("sans-serif", font_size_legend))
        .position(SeriesLabelPosition::UpperLeft) // position on the chart
        .draw()?;

    root.present().expect(&format!(
        "Unable to write result to file, please make sure you have write permissions for {}",
        file_path.display()
    ));

    return Ok(());
}

fn draw_extrap(
    spline_no_extrap: &AkimaSpline,
    spline_extrapl: &AkimaSpline,
    spline_extrapr: &AkimaSpline,
    filename: &str,
    legend_left: &str,
    legend_right: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    // General config
    // =========================================================================
    let size = (600, 400);
    let font_size_labels = 18;
    let font_size_ticks = 16;
    let font_size_legend = 18;

    // Calculate data
    // =========================================================================

    let resolution = 1e-3;
    let xmin = -3.0;
    let xmax = 3.0;
    let capacity = ((xmax - xmin) / resolution) as usize;

    let mut x_no_extrap = Vec::with_capacity(capacity);
    let mut y_no_extrap = Vec::with_capacity(capacity);

    let mut x_extrapl = Vec::with_capacity(capacity);
    let mut y_extrapl = Vec::with_capacity(capacity);

    let mut x_extrapr = Vec::with_capacity(capacity);
    let mut y_extrapr = Vec::with_capacity(capacity);

    let mut x = xmin;
    while x <= xmax {
        if let Some(y) = spline_no_extrap.eval(x) {
            x_no_extrap.push(x);
            y_no_extrap.push(y);
        }

        if let Some(y) = spline_extrapl.eval(x) {
            x_extrapl.push(x);
            y_extrapl.push(y);
        }

        if let Some(y) = spline_extrapr.eval(x) {
            x_extrapr.push(x);
            y_extrapr.push(y);
        }

        x += resolution;
    }

    let file_path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(&format!("../{filename}.svg"));
    let root = SVGBackend::new(&file_path, size).into_drawing_area();

    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(xmin..xmax, -2.0..4.0)?;

    chart
        .configure_mesh()
        .axis_desc_style(("sans-serif", font_size_labels))
        .label_style(("sans-serif", font_size_ticks))
        .draw()?;

    // Left-side extrapolation
    chart
        .draw_series(LineSeries::new(
            x_extrapl.iter().cloned().zip(y_extrapl.iter().cloned()),
            &RED,
        ))?
        .label(legend_left)
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    // Right-side extrapolation
    chart
        .draw_series(LineSeries::new(
            x_extrapr.iter().cloned().zip(y_extrapr.iter().cloned()),
            &GREEN_800,
        ))?
        .label(legend_right)
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN_800));

    // No extrapolation
    chart
        .draw_series(LineSeries::new(
            x_no_extrap.iter().cloned().zip(y_no_extrap.iter().cloned()),
            &BLUE,
        ))?
        .label("no extrapolation")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    // Original data points
    chart
        .draw_series(
            spline_no_extrap
                .xs()
                .iter()
                .cloned()
                .zip(spline_no_extrap.ys().iter().cloned())
                .map(|(x, y)| Circle::new((x, y), 5.0, BLACK.filled())),
        )?
        .label("original data")
        .legend(|(x, y)| Circle::new((x + 10, y), 5.0, BLACK.filled()));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8)) // semi-transparent background
        .border_style(&BLACK)
        .label_font(("sans-serif", font_size_legend))
        .position(SeriesLabelPosition::UpperLeft) // position on the chart
        .draw()?;

    root.present().expect(&format!(
        "Unable to write result to file, please make sure you have write permissions for {}",
        file_path.display()
    ));

    return Ok(());
}
