use akima_spline::AkimaSpline;
use plotters::{prelude::*, style::full_palette::GREEN_800};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // General config
    // =========================================================================
    let size = (600, 400);
    let font_size_labels = 18;
    let font_size_ticks = 16;
    let font_size_legend = 18;

    // Spline definition
    // =========================================================================

    let xs = vec![-2.0, -1.0, -0.8, 0.25, 1.0, 2.0];
    let ys = vec![0.0, 1.0, 2.0, -1.5, 1.0, 3.0];
    let extrapl = vec![2.0];
    let extrapr = vec![3.0, -2.0];

    let spline_no_extrap =
        AkimaSpline::new(xs.clone(), ys.clone(), None, None).expect("valid input data");
    let spline_extrapl =
        AkimaSpline::new(xs.clone(), ys.clone(), Some(extrapl), None).expect("valid input data");
    let spline_extrapr =
        AkimaSpline::new(xs.clone(), ys.clone(), None, Some(extrapr)).expect("valid input data");

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

    // readme_example_no_extrap
    // =========================================================================

    {
        let file_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../readme_example_no_extrap.svg");
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
                xs.iter()
                    .cloned()
                    .zip(ys.iter().cloned())
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
    }

    // readme_example_extrap
    // =========================================================================

    {
        let file_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../readme_example_extrap.svg");
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
            .label("left-side extrapolation (y = 2x + k)")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

        // Right-side extrapolation
        chart
            .draw_series(LineSeries::new(
                x_extrapr.iter().cloned().zip(y_extrapr.iter().cloned()),
                &GREEN_800,
            ))?
            .label("right-side extrapolation (y = 3x² - 2x + k)")
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
                xs.iter()
                    .cloned()
                    .zip(ys.iter().cloned())
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
    }

    return Ok(());
}
