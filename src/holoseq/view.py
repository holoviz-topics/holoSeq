from bisect import bisect_left
from pathlib import Path

import holoviews as hv
import numpy as np
import panel as pn
from holoviews.operation.datashader import dynspread, rasterize
from holoviews.operation.resample import ResampleOperation2D

from holo_seq.data import load_2d

hv.extension("bokeh")
dynspread.max_px = 8
dynspread.threshold = 0.5
ResampleOperation2D.width = 1500
ResampleOperation2D.height = 1500


def show_H1(names, starts, x, y, width):
    if np.isnan(x) or np.isnan(y):
        s = "Mouse click on image for current location"
    else:
        index_x = bisect_left(starts, x)
        chromosome_x = names[index_x - 1]
        offset_x = x - starts[index_x - 1]
        index_y = bisect_left(starts, y)
        chromosome_y = names[index_y - 1]
        offset_y = y - starts[index_y - 1]
        s = f"X axis {chromosome_x}:{offset_x} Y axis {chromosome_y}:{offset_y}"
    str_pane = pn.pane.Str(
        s,
        styles={"font-size": "10pt", "color": "darkblue", "text-align": "center"},
        width=width,
    )
    return str_pane


def maker(path: Path, width: int):
    haploids, x_coords, y_coords, annotations, title = load_2d(path=path)
    points = hv.Points(np.column_stack((x_coords, y_coords)))
    tick_names = list(zip(haploids["H1"]["positions"], haploids["H1"]["contig_names"]))
    stream = hv.streams.Tap(source=points, x=np.nan, y=np.nan)
    show_location = pn.bind(
        show_H1,
        names=haploids["H1"]["contig_names"],
        starts=haploids["H1"]["positions"],
        x=stream.param.x,
        y=stream.param.y,
        width=width,
    )

    holo_seq_heatmap = pn.Column(
        show_location,
        pn.pane.HoloViews(
            dynspread(
                rasterize(points)
                .relabel("%s" % title)
                .opts(
                    cmap="inferno",
                    cnorm="log",
                    colorbar=True,
                    width=width,
                    height=width,
                    xticks=tick_names,
                    yticks=tick_names,
                    xrotation=45,
                    fontsize={"xticks": 5, "yticks": 5},
                    tools=["tap"],
                    shared_axes=False,
                    scalebar=True,
                    scalebar_range="x",
                    scalebar_location="top_left",
                    scalebar_unit=("bp"),
                    show_grid=True,
                )
            )
        ),
    )
    return holo_seq_heatmap, title
