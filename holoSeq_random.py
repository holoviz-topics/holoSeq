    

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

SQRT2 = math.sqrt(2)

class rotater():
    """
    moved into a class so the constants are only calculated once
    if single pairs are being rotated
    x 0 y 0 xr -2121320 yr 878679
    x 3000000 y 0 xr 0 yr -1242640
    x 1500000 y 1500000 xr 0 yr 878679
    x 3000000 y 3000000 xr 2121320 yr 878679
    """
    def __init__(self, xwidth, ywidth):
        self.rotate = True
        self.onepointRot = True
        self.radians=0.7853981633974483
        self.origin=(0, ywidth)
        self.cos_rad = math.cos(self.radians)
        self.sin_rad = math.sin(self.radians)
        (self.xmin, self.ymax) = self.rotatecoords(0,0, adjust=False)
        self.ymin = self.rotatecoords(xwidth,0, adjust=False)[1]
        self.xmax = self.rotatecoords(xwidth,ywidth, adjust=False)[0]
        self.xnew = SQRT2*xwidth
        self.ynew = ywidth/SQRT2
        self.xwidth = xwidth
        self.ywidth = ywidth
        self.xscalefact = self.xwidth/self.xnew
        self.yscalefact = self.ywidth/self.ynew


    def rotatecoords(self,
            xin=0, yin=0, adjust=True,
        ):
            """
            This version optimised to operate on pair at a time so precomputed constants
            # make a rotated heatmap where the diagonal becomes the horizontal x axis
            # https://gist.github.com/LyleScott/d17e9d314fbe6fc29767d8c5c029c362
            # this inflates the xaxis by sqrt(2) but can rescale and it seems to work (TM)
            # xcis1r, ycis1r = rotatecoords(xcis1, ycis1, radians=0.7853981633974483, origin=(max(xcis1),max(ycis1)))
            # pafxycis1 = pd.DataFrame(np.vstack([xcis1r,ycis1r]).T, columns = ['x', 'y'])
            # y height is complicated - h^2 + (1/2*sqrt2*xwdith)^2 = y^2
            """
            self.offset_x, self.offset_y = self.origin
            adjusted_x = xin - self.offset_x
            adjusted_y = yin - self.offset_y
            qx = self.offset_x + self.cos_rad * adjusted_x + self.sin_rad * adjusted_y
            qy = self.offset_y + -self.sin_rad * adjusted_x + self.cos_rad * adjusted_y
            if adjust:
                xrs = (qx-self.xmin)*self.xscalefact
                yrs = (qy-self.ymin)*self.yscalefact
            else:
                xrs = qx
                yrs = qy
            if self.onepointRot:
                return (xrs, yrs)
            else:
                xyr = pd.DataFrame(np.vstack([xrs,yrs]).T, columns=["x", "y"])
                return xyr

    def unrotatecoords(self, xr, yr):
        """
        rotated coords back calculate - not needed - rescaled rotated coords are correct
        this is simplified because cos_rad == sin_rad for this 45 degree rotation
        """
        qx = xr/self.xscalefact + self.xmin
        qy = yr/self.yscalefact + self.ymin
        adjx = (self.offset_y - qy + qx - self.offset_x)/2
        adjy = (qx - self.offset_x - adjx)
        x = adjx/self.cos_rad + self.offset_x
        y = adjy/self.cos_rad + self.offset_y
        return x,y



def showTap(x, y):
    # populate a string widget with tap coordinates translated into contig and offset
    if np.isnan(x) or np.isnan(y):
        s = "Click on plot for coordinates as contig:offset"
    else:
        chrx = "Out of range"
        offsx = 0
        chry = "Out of range"
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


def showRotTap(x, y, rot):
    # populate a string widget with tap coordinates translated into contig and offset
    if np.isnan(x) or np.isnan(y):
        s = "Click on plot for coordinates as contig:offset"
    else:
        chrx = "Out of range"
        offsx = 0
        chry = "Out of range"
        offsy = 0
        xur, yur = rot.unrotatecoords(x,y)
        if yur <= xur:
            i = bisect_left(hstarts, xur)
            if i > 0 and i <= len(hnames):
                chrx = hnames[i - 1]
                offsx = xur - hstarts[i - 1]
            i = bisect_left(hstarts, yur)
            if i > 0 and i <= len(hnames):
                chry = hnames[i - 1]
                offsy = yur - hstarts[i - 1]
        s = "X axis %s:%d Y axis %s:%d x %d y %d xur %d yur %d" % (chrx, offsx, chry, offsy, x, y, xur, yur)
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
rot = rotater(xwidth, xwidth)
rot.onepointRot = True
rng = np.random.default_rng(1)  # all plots will be identical !
hlen = [xwidth / 2, xwidth / 4, xwidth / 8, xwidth / 8]
hstarts = list(itertools.accumulate(hlen))
hstarts.insert(0, 0)
hnames = ["chr%d" % i for i in range(1, 5)]
title = "%d random points" % xmax
ticks = [(hstarts[i], hnames[i]) for i in range(4)]
xcoords = np.array([rng.uniform(0, xwidth) for i in range(xmax)])
ycoords = np.array([xcoords[i] - abs(rng.uniform(0, xcoords[i])) for i in range(xmax)])
points = hv.Points((xcoords, ycoords))
stream = hv.streams.Tap(source=points, x=np.nan, y=np.nan)
showloc = pn.bind(showTap, x=stream.param.x, y=stream.param.y)
xr = []
yr = []
for i, x in enumerate(xcoords):
    xrr, yrr = rot.rotatecoords(x,ycoords[i])
    xr.append(xrr)
    yr.append(yrr)
print('maxxr %f minxxr %f maxyr %f minyr %f' % (max(xr),min(xr), max(yr),min(yr)))
xyr = pd.DataFrame(np.vstack([xr, yr]).T, columns=["x", "y"])
rotpoints = hv.Points(xyr, kdims=["x", "y"])
streamr = hv.streams.Tap(source=rotpoints, x=np.nan, y=np.nan)
showlocr = pn.bind(showRotTap, x=streamr.param.x, y=streamr.param.y, rot=rot)
rot = pn.pane.HoloViews(
    dynspread(
        rasterize(rotpoints)
        .relabel("%s" % title)
        .opts(
            # ylim=(0,xwidth),
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
            shared_axes=False,
        )
    )
)

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
            shared_axes=False,
        )
    )
)
pnc = pn.Column(showlocr, rot, showloc, unrot)


pnc.servable(
    title=title,
)
