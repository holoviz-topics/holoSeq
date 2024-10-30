# see https://github.com/fubar2/holoSeq
# random data 2d example
from bisect import bisect_left
import itertools
import numpy as np

import holoviews as hv
import panel as pn

from holoviews.operation.datashader import (
    rasterize,
    dynspread,
)

hv.extension("bokeh")
pn.extension()


def showTap(x, y):
    # populate a string widget with tap coordinates translated into contig and offset
    if np.isnan(x) or np.isnan(y):
        s = "Click on plot for coordinates as contig:offset"
    else:
        i = bisect_left(hstarts, x)
        chrx = hnames[i - 1]
        offsx = x - hstarts[i - 1]
        i = bisect_left(hstarts, y)
        chry = hnames[i - 1]
        offsy = y - hstarts[i - 1]
        s = "X axis %s:%d Y axis %s:%d" % (chrx, offsx, chry, offsy)
    str_pane = pn.pane.Str(
        s,
        styles={"font-size": "10pt", "color": "darkblue", "text-align": "center"},
        width=width,
    )
    return str_pane


xmax = 10000
width = 1200
rng = np.random.default_rng(1)  # all plots will be identical !
hlen = sorted(rng.exponential(10000000,10), reverse=True)
hstarts = [int(x) for x in list(itertools.accumulate(hlen))]
xwidth = hstarts[-1] + hlen[-1]
hstarts.insert(0, 0)
hnames = ["chr%d" % i for i in range(1, 11)]
xcoords = np.array([rng.uniform(0, xwidth) for i in range(xmax)])
yoffset = rng.normal(0, 1000000, xmax)
ycoords = np.array(xcoords + yoffset)
title = "%d random points" % xmax
ticks = [(hstarts[i], hnames[i]) for i in range(len(hnames))]
print(ticks)
points = hv.Points((xcoords, ycoords))
stream = hv.streams.Tap(source=points, x=np.nan, y=np.nan)
showloc = pn.bind(showTap, x=stream.param.x, y=stream.param.y)
pnc = pn.Column(
    showloc,
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
                xticks=ticks,
                yticks=ticks,
                xrotation=45,
                fontsize={"xticks": 7, "yticks": 7},
                tools=["tap"],
                scalebar=True,
                scalebar_range="x",
                scalebar_location="top_left",
                scalebar_unit=("bp"),
                show_grid=True
            )
        )
    ),
)


pnc.servable(title=title)
