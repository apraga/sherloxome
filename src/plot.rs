//! Visualise benchmarking results as interactive Vega-Lite boxplots.
//!
//! Reads the `merged.csv` produced by [`crate::benchmark`] and opens a browser window
//! showing F1-score distributions broken down by patient, sequencer, depth, and capture kit.

// Column indices are harcoded to match merged.csv:
// 0=Type, 1=Filter, 13=METRIC.F1_Score, 19=patient, 20=capture, 21=sequencer, 22=depth

use std::path::PathBuf;
use vega_lite_4::*;

/// Build one horizontal row of the faceted boxplot for a given factor column index.
fn make_row(
    factor_idx: &str,
    show_column_labels: bool,
) -> Result<NormalizedSpec, Box<dyn std::error::Error>> {
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
            make_row("19", true)?,  // patient
            make_row("20", false)?, // capture
            make_row("21", false)?, // sequencer
            make_row("22", false)?, // depth
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
