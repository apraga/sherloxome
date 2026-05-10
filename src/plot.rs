//! Visualise benchmarking results as interactive Vega-Lite boxplots.
//!
//! Reads the `merged.csv` produced by [`crate::analyze`] and opens a browser window
//! showing F1-score distributions broken down by patient, sequencer, depth, and capture kit.

// Column indices in merged.csv:
// 0=Type, 1=Filter, 13=METRIC.F1_Score, 18=patient, 19=capture, 20=sequencer, 21=depth

use std::path::PathBuf;
use vega_lite_4::*;

/// Build one horizontal row of the faceted boxplot for a given factor column index.
fn make_row(
    factor_idx: &str,
    show_column_labels: bool,
) -> Result<NormalizedSpec, Box<dyn std::error::Error>> {
    // let sort = if factor_idx == "21" {
    //     serde_json::from_value(serde_json::json!(["50x", "75x", "100x"]))?
    // } else {
    //     serde_json::from_value(serde_json::json!("ascending"))?
    // };

    Ok(NormalizedSpecBuilder::default()
        .mark(Mark::Boxplot)
        .encoding(
            EdEncodingBuilder::default()
                .y(YClassBuilder::default()
                    .field(factor_idx)
                    .position_def_type(Type::Nominal)
                    // .sort(sort)
                    .build()?)
                .x(XClassBuilder::default()
                    .field("13")
                    .position_def_type(Type::Quantitative)
                    .scale(ScaleBuilder::default().zero(false).build()?)
                    .title("F1 score")
                    .build()?)
                .color(
                    ColorClassBuilder::default()
                        .field(factor_idx)
                        .mark_prop_def_gradient_string_null_type(Type::Nominal)
                        .legend(RemovableValue::Remove)
                        .build()?,
                )
                .column(
                    RowColumnEncodingFieldDefBuilder::default()
                        .field("0")
                        .title(RemovableValue::Remove)
                        .header(
                            HeaderBuilder::default()
                                .labels(show_column_labels)
                                .build()?,
                        )
                        .build()?,
                )
                .build()?,
        )
        .resolve(
            ResolveBuilder::default()
                .scale(
                    ScaleResolveMapBuilder::default()
                        .x(ResolveMode::Independent)
                        .build()?,
                )
                .build()?,
        )
        .build()?)
}

/// Read `merged.csv` and open a browser with F1-score boxplots faceted by variant type.
pub fn plot(input: PathBuf, _output: PathBuf) -> Result<(), Box<dyn std::error::Error>> {
    let chart = VegaliteBuilder::default()
        .data(csv::Reader::from_path(&input)?)
        .transform(vec![
            TransformBuilder::default()
                .filter("datum[1]==='PASS'")
                .build()?,
        ])
        .vconcat(vec![
            make_row("18", true)?,  // patient
            make_row("20", false)?, // sequencer
            make_row("21", false)?, // depth
            make_row("19", false)?, // capture
        ])
        .spacing(Spacing::Double(30.0))
        .config(
            ConfigBuilder::default()
                .facet(CompositionConfigBuilder::default().spacing(5.0).build()?)
                .view(
                    ViewConfigBuilder::default()
                        .stroke(RemovableValue::Remove)
                        .build()?,
                )
                .build()?,
        )
        .build()?;

    chart.show()?;
    Ok(())
}
