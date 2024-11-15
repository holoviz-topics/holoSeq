# see https://github.com/fubar2/holoSeq
# illustrates some of the basic ideas in converting
# a set of features that have been mapped to a genome into a 
# linear or in this case 2D display.
# a tap will reveal the current x and y positions backcalculated from the
# cumulated start positions of the contigs.
# run with
# panel serve holoSeq_random.py --show
# should pop up 10000 points.
# try 10 million - still works smoothly
# Needs depenencies
# pip install datashader dask[dataframe] holoviews[recommended] pandas matplotlib bokeh
#

from bisect import bisect_left
from collections import OrderedDict
import gzip
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

xwidth = 3000000
xmax = 10000
width = 1000
rng = np.random.default_rng(1) # all plots will be identical !
hlen = [xwidth/2, xwidth/4, xwidth/8, xwidth/8]
hstarts = list(itertools.accumulate(hlen))
hstarts.insert(0,0)
hnames = ['chr%d' % i for i in range(1,5)]
xcoords = np.array([rng.uniform(0,xwidth) for i in range(xmax)])
yoffset = rng.normal(0,100000,xmax)
ycoords = np.array(xcoords + yoffset)
title = "%d random points" % xmax
ticks = [(hstarts[i], hnames[i]) for i in range(4)]
points = hv.Points((xcoords,ycoords))
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
            )
        )
    ),
)


pnc.servable(title=title, )
