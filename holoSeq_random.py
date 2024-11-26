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
import math
import numpy as np

import holoviews as hv
import pandas as pd
import panel as pn

from holoviews.operation.datashader import (
    rasterize,
    dynspread,
)

hv.extension("bokeh")
pn.extension()



def rotatecoords(xin, yin, radians=0.7853981633974483, origin=(0, 0), xwidth=3000000, ywidth=3000000):
    # make a rotated heatmap where the diagonal becomes the horizontal x axis
    # https://gist.github.com/LyleScott/d17e9d314fbe6fc29767d8c5c029c362
    # this inflates the xaxis by sqrt(2) but can rescale and it seems to work (TM)
    # xcis1r, ycis1r = rotatecoords(xcis1, ycis1, radians=0.7853981633974483, origin=(max(xcis1),max(ycis1)))
    # pafxycis1 = pd.DataFrame(np.vstack([xcis1r,ycis1r]).T, columns = ['x', 'y'])
    # y height is complicated - h^2 + (1/2*sqrt2*xwdith)^2 = y^2 
    sqrt2 = math.sqrt(2)
    offset_x, offset_y = origin
    adjusted_x = xin - offset_x
    adjusted_y = yin - offset_y
    cos_rad = math.cos(radians)
    sin_rad = math.sin(radians)
    qx = offset_x + cos_rad * adjusted_x + sin_rad * adjusted_y
    qy = offset_y + -sin_rad * adjusted_x + cos_rad * adjusted_y
    print('max/min qy=', min(qy), max(qy))
    # must rescale x after inflation from rotation
    xdelta = xwidth*sqrt2
    xmin = xwidth - xdelta
    ydelta = (ywidth**2 - (xdelta/2)**2)**0.5
    ymin = (1.0 - math.sin(radians)) * ywidth
    xr = [(x - xmin)/xdelta*xwidth for x in qx]
    yr = [(x - ymin)/ydelta*ywidth for x in qy]
    xyr = pd.DataFrame(np.vstack([xr,yr]).T, columns = ['x', 'y'])
    return xyr


def showTap(x, y):
    # populate a string widget with tap coordinates translated into contig and offset
    if np.isnan(x) or np.isnan(y):
        s = "Click on plot for coordinates as contig:offset"
    else:
        chrx = "Out of range"
        offsx = 0
        chry= "Out of range"
        offsy = 0
        i = bisect_left(hstarts, x)
        if i > 0 and i <= len(hnames):
            chrx = hnames[i - 1]
            offsx = x - hstarts[i - 1]
        i = bisect_left(hstarts, y)
        if i > 0 and i <= len(hnames):
            chry = hnames[i - 1]
            offsy = y - hstarts[i - 1]
        s = "X axis %s:%d Y axis %s:%d" % (chrx, offsx, chry, offsy)
    str_pane = pn.pane.Str(
        s,
        styles={"font-size": "10pt", "color": "darkblue", "text-align": "center"},
        width=width,
    )
    return str_pane


def showRotTap(x, y):
    # populate a string widget with tap coordinates translated into contig and offset
    if np.isnan(x) or np.isnan(y):
        s = "Click on plot for coordinates as contig:offset"
    else:
        chrx = "Out of range"
        offsx = 0
        chry= "Out of range"
        offsy = 0
        i = bisect_left(hstarts, x)
        if i > 0 and i <= len(hnames):
            chrx = hnames[i - 1]
            offsx = x - hstarts[i - 1]
        i = bisect_left(hstarts, y)
        if i > 0 and i <= len(hnames):
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
height = 400
rng = np.random.default_rng(1) # all plots will be identical !
hlen = [xwidth/2, xwidth/4, xwidth/8, xwidth/8]
hstarts = list(itertools.accumulate(hlen))
hstarts.insert(0,0)
hnames = ['chr%d' % i for i in range(1,5)]
title = "%d random points" % xmax
ticks = [(hstarts[i], hnames[i]) for i in range(4)]
xcoords = np.array([rng.uniform(0,xwidth) for i in range(xmax)])
ycoords = np.array([xcoords[i] - abs(rng.uniform(0,xcoords[i])) for i in range(xmax)])
points = hv.Points((xcoords,ycoords))
stream = hv.streams.Tap(source=points, x=np.nan, y=np.nan)
showloc = pn.bind(showTap, x=stream.param.x, y=stream.param.y)
xyr = rotatecoords(xcoords,ycoords, radians=0.7853981633974483, origin=(max(xcoords),max(ycoords)), xwidth=xwidth, ywidth=xwidth)
rotpoints = hv.Points(xyr, kdims=["x","y"])
streamr = hv.streams.Tap(source=rotpoints, x=np.nan, y=np.nan)
showlocr = pn.bind(showTap, x=streamr.param.x, y=streamr.param.y)
rot = pn.pane.HoloViews(dynspread(rasterize(rotpoints).relabel("%s" % title)
            .opts(
                #ylim=(0,xwidth),
                framewise=True,
                autorange=None,
                cmap="inferno",
                cnorm="log",
                colorbar=True,
                width=width,
                height=height,
                xticks=ticks,
                yticks=ticks,
                xrotation=45,
                fontsize={"xticks": 7, "yticks": 7},
                tools=["tap"],
                scalebar=True,
                scalebar_range="x",
                scalebar_location="top_left",
                scalebar_unit=("bp"),
                shared_axes=False
            )))

unrot = pn.pane.HoloViews(
        dynspread(
            rasterize(points)
            .relabel("%s" % title)
            .opts(
                cmap="inferno",
                cnorm="log",
                colorbar=True,
                width=width,
                height=height,
                xticks=ticks,
                yticks=ticks,
                xrotation=45,
                fontsize={"xticks": 7, "yticks": 7},
                tools=["tap"],
                scalebar=True,
                scalebar_range="x",
                scalebar_location="top_left",
                scalebar_unit=("bp"),
                shared_axes=False
            )
        )
    )
pnc = pn.Column(
    showlocr,
    rot,
    showloc,
    unrot
)


pnc.servable(title=title, )
