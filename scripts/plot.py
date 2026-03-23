import altair as alt
import polars as pl

# Remove duplicates
df = pl.read_csv("../test/merged.csv").filter(pl.col('Filter') == 'PASS')


def make_row(factor):
    sort = ['50x', '75x', '100x'] if factor == 'depth' else 'ascending'
    return alt.Chart(df).mark_boxplot().encode(
        alt.Y(f'{factor}:N',sort=sort),
        alt.X('METRIC\\.F1_Score:Q', scale=alt.Scale(zero=False), title='F1 score'),
        color=alt.Color(f'{factor}:N', legend=None),
        column=alt.Column('Type:N', title=None, header=alt.Header(
            labels=True if factor == 'patient' else False  # only show SNP/INDEL on top row
        ))
    ).resolve_scale(x='independent')

chart = alt.vconcat(
    *[make_row(f) for f in ['patient', 'sequencer', 'depth', 'capture']],
    spacing=30
).configure_facet(spacing=5).configure_view(stroke=None)
chart.save('plot.html')
