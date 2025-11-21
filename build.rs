fn main() {
    // If building for docs.rs, DO NOT create the README files from the template
    if let Ok(env) = std::env::var("DOCS_RS") {
        if &env == "1" {
            return ();
        }
    }

    let mut readme = std::fs::read_to_string("README.template.md").unwrap();
    readme = readme.replace(
        "{{VERSION}}",
        std::env::var("CARGO_PKG_VERSION")
            .expect("version is available in build.rs")
            .as_str(),
    );

    // Generate README_local.md using local images
    let mut local = readme.replace("{{no_extrap.svg}}", "docs/no_extrap.svg");
    local = local.replace("{{extrap_1.svg}}", "docs/extrap_1.svg");
    local = local.replace("{{extrap_2.svg}}", "docs/extrap_2.svg");
    std::fs::write("README_local.md", local).unwrap();

    // Generate README,md using online hosted images
    let mut docsrs = readme.replace(
        "{{example.svg}}",
        "https://raw.githubusercontent.com/StefanMathis/akima_spline/refs/heads/main/docs/no_extrap.svg",
    );
    docsrs = docsrs.replace("{{extrap_1.svg}}", "https://raw.githubusercontent.com/StefanMathis/akima_spline/refs/heads/main/docs/extrap_1.svg");
    docsrs = docsrs.replace("{{extrap_2.svg}}", "https://raw.githubusercontent.com/StefanMathis/akima_spline/refs/heads/main/docs/extrap_2.svg");
    std::fs::write("README.md", docsrs).unwrap();
}
