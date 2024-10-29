# see https://github.com/fubar2/holoSeq
# RUN THIS to create an interactive displays in an IPython notebook
# after loading dependencies with
# ! pip install datashader dask[dataframe] holoviews[recommended] pandas matplotlib bokeh
#
# The data are Arima HiC reads from the Arctic Ground Squirrel mUroPar1 VGP genomeArk repository
# processed with Dephine's Pretext workflow using Bellerophon to remove chimeric reads
# The paired bam output is converted to PAF with an awk script (!) and that's what is read
# in this code.
# The pairs are parsed to extract the haplotype designator from the contig name
# typically a suffix like H1 extracted in getHap - rewrite that to suit your names.
# Sorting by contig name is based on VGP conventions - SUPER_ first, then scaffolds
# One problem to watch out for is that any differences in ordering of the X and Y contigs can make all sorts of
# artifacts appear such as the kaleidoscopic patterns seen in Pretextviewer.
# Ross Lazarus October 2024
# This holoviews application is mostly monolithic because it cannot easily be
# split up without passing lots of parameters AFAIK.
# it works. be happy.
from bisect import bisect_left
from collections import OrderedDict
import gzip
import numpy as np
import os
import sys

import holoviews as hv
import panel as pn
import pandas as pd

from holoviews.operation.datashader import (
    rasterize,
    dynspread,
)
from holoviews.operation.resample import ResampleOperation2D
from holoviews.operation import decimate

# inFile = "galaxy_inputs/paf/bothmap.paf.tab.tabular"
inFile = "/home/ross/rossgit/holoviews-examples/holoSeqtest.gz"
holoSeqHeaders = ["@v1HoloSeq1D", "@v1HoloSeq2D"]
hv.extension("bokeh")
pn.extension()


def getHap(contig):
    """
    function to return suffix H1 from chrH1 - adjust to suit.
    """
    return contig[-2:]


def xportHtml(fname, hObj):
    "save a holoview object to an interactive but not adaptive scaling HTML page"
    hv.save(filename=fname, obj=hObj)

def showH1(x, y):
    if np.isnan(x) or np.isnan(y):
        s = "Mouse click on image for location"
    else:
        i = bisect_left(h1starts, x)
        chrx = h1names[i - 1]
        offsx = x - h1starts[i - 1]
        i = bisect_left(h1starts, y)
        chry = h1names[i - 1]
        offsy = y - h1starts[i - 1]
        s = "X axis %s:%d Y axis %s:%d" % (chrx, offsx, chry, offsy)
    str_pane = pn.pane.Str(
        s,
        styles={"font-size": "10pt", "color": "darkblue", "text-align": "center"},
        width=pcwidth,
    )
    return str_pane


def showH2(x, y):
    if np.isnan(x) or np.isnan(y):
        s = "Mouse click on image for location"
    else:
        i = bisect_left(h2starts, x)
        chrx = h2names[i - 1]
        offsx = x - h2starts[i - 1]
        i = bisect_left(h2starts, y)
        chry = h2names[i - 1]
        offsy = y - h2starts[i - 1]
        s = "X axis %s:%d Y axis %s:%d" % (chrx, offsx, chry, offsy)
    str_pane = pn.pane.Str(
        s,
        styles={"font-size": "10pt", "color": "darkblue", "text-align": "center"},
        width=pcwidth,
    )
    return str_pane


def showTrans(x, y):
    if np.isnan(x) or np.isnan(y):
        s = "Mouse click on image for location"
    else:
        i = bisect_left(h1starts, x)
        chrx = h1names[i - 1]
        offsx = x - h1starts[i - 1]
        i = bisect_left(h2starts, y)
        chry = h2names[i - 1]
        offsy = y - h2starts[i - 1]
        s = "X axis %s:%d Y axis %s:%d" % (chrx, offsx, chry, offsy)
    str_pane = pn.pane.Str(
        s,
        styles={"font-size": "10pt", "color": "darkblue", "text-align": "center"},
        width=pcwidth,
    )
    return str_pane

def import_holoSeq_data(inFile):
    """
    reverse process of dumping the data in holoSeq format from a converter 
    see https://github.com/fubar2/holoSeq/blob/main/HoloSeqOverview.md 
    """
    haps = {}
    xcoords = []
    ycoords = []
    annos = []
    hsDims = None
    plotType = None
    title = "Plot"
    with gzip.open(inFile, 'rt') as f:

        for i,row in enumerate(f.readlines()):
            if i == 0:
                hseqformat = row.split()
                if hseqformat[0] not in holoSeqHeaders:
                    print("Supplied input", inFile, "has a first row =", row,"so is not a valid holoSeq input file")
                    print("First row must start with one of these:", holoSeqHeaders)
                    return None
                hsDims = holoSeqHeaders.index(hseqformat[0]) + 1
                print('hsDims=', hsDims)
                if hsDims == "1":
                    plotType = 'bar'
                    if len(hseqformat) > 1:
                        plotType = hseqformat[1].strip()
            else:
                if row.startswith('@'):
                    if row.startswith('@title'):
                        title = row.replace('@title','').strip()
                    else:
                        srow = row.split()
                        if len(srow) >= 3: 
                            hap, cname, clen = row.split()[:3]
                            if not haps.get(hap, None):
                                haps[hap] = {'cn': [], 'cl': []}
                            haps[hap]['cn'].append(cname.strip())
                            haps[hap]['cl'].append(int(clen.strip()))
                        else:
                            print("Supplied input", inFile, "at row", i, "=", row, "lacking the required reference name, contig name and contig length. Not a valid holoSeq input file")
                            return None
                else:
                    srow = [x.strip() for x in row.split()]
                    lrow = len(srow)
                    if hsDims == 2:
                        if lrow < 2:
                            print("At row", i, "Supplied 2D input", inFile, "has", row,"which does not parse into x and ycoordinates and optional annotation")
                            return None
                        else:
                            if srow[0].isdigit() and srow[1].isdigit():
                                xcoords.append(int(srow[0]))
                                ycoords.append(int(srow[1]))
                                if lrow > 2:
                                    annos.append(srow[2:])
                            else:
                                print("At row", i, "Supplied 2D input", inFile, "has", row,"which does not parse into x and ycoordinates and optional annotation")
                                return None
                    else:
                        if srow[0].isdigit():
                                xcoords.append(int(srow[0]))
                                if lrow > 1:
                                    annos.append(srow[1:])
                        else:
                            print("At row", i, "Supplied 1D input", inFile, "has", row,"which does not parse into x coordinate and optional annotation")
                            return None
    return((hsDims, haps, xcoords, ycoords, annos, plotType, title))

pcwidth = 800
# width settings for plots and location bars
# Default values suitable for this notebook
decimate.max_samples = 10000
dynspread.max_px = 10
dynspread.threshold = 0.6
ResampleOperation2D.width = 1500
ResampleOperation2D.height = 1500
# need to convert the categorical contigs into a sequence for holoviews to munch
# use the contig length from the paf to figure out the cumulative start for each contig
# contigs are length ordered - that does not always work well when the haplotypes differ widely
inFile = 'mUroPar1_cis1.hseq.gz'
print('Infile = ', inFile)
hqstarts = OrderedDict()
hlens = {}
haps = []
(hsDims, hapsread, xcoords, ycoords, annos, plotType, plotTitle) = import_holoSeq_data(inFile)
print(len(xcoords), len(ycoords))
for hap in hapsread.keys():
    haps.append(hap)
    cum = 1
    hqstarts[hap] = OrderedDict()
    for i, contig in enumerate(hapsread[hap]['cn']):
        clen = hapsread[hap]['cl'][i]
        hqstarts[hap][contig] = cum
        cum += clen
hap = haps[0]
h1starts = [hqstarts[hap][x] for x in hqstarts[haps[0]].keys()]
h1names = list(hqstarts[hap].keys())
qtic1 = [(hqstarts[hap][x], x) for x in hqstarts[hap].keys()]
isTrans = False
if hsDims == "2": # may be one or two distinct haplotype identifiers - H1 +/- H2
    if len(haps) > 1:
        hap = haps[1]
        isTrans = True
h2starts = [hqstarts[hap][x] for x in hqstarts[hap].keys()]
h2names = list(hqstarts[hap].keys())
qtic2 = [(hqstarts[hap][x], x) for x in hqstarts[hap].keys()]
# can take the np.tril or filter the upper triangle while processing pairs
# and rotate so the diagonal becomes the x axis but need some kind of
# sideways scroller to work right
# xcis1r, ycis1r = rotatecoords(xcis1, ycis1, radians=0.7853981633974483, origin=(max(xcis1),max(ycis1)))
# pafxycis1 = pd.DataFrame(np.vstack([xcis1r,ycis1r]).T, columns = ['x', 'y'])
# --------------------cut here---------------------------------
# once the pairs have been read and mapped into a grid, the code
# below does the plotting.
# it can be copied, edited to suit your needs and
# run repeatedly without waiting for the data to be mapped.
pafxy = pd.DataFrame.from_dict({'x': xcoords, 'y': ycoords})
pafp = hv.Points(pafxy)
stream = hv.streams.Tap(source=pafp, x=np.nan, y=np.nan)
if isTrans:
    showloc = pn.bind(showTrans, x=stream.param.x, y=stream.param.y)
else:
    showloc = pn.bind(showH1, x=stream.param.x, y=stream.param.y)
# to rotate so the diagonal becomes the x axis
# xcis1r, ycis1r = rotatecoords(xcis1, ycis1, radians=0.7853981633974483, origin=(max(xcis1),max(ycis1)))
# pafxycis1 = pd.DataFrame(np.vstack([xcis1r,ycis1r]).T, columns = ['x', 'y'])
# bisect.bisect_left(a, x, lo=0, hi=len(a), *, key=None)
# prepare and show the 3 plots
p1 = pn.Column(
    showloc,
    pn.pane.HoloViews(
        dynspread(
            rasterize(pafp)
            .relabel("%s" % plotTitle)
            .opts(
                cmap="inferno",
                cnorm="log",
                colorbar=True,
                width=pcwidth,
                height=pcwidth,
                xticks=qtic1,
                yticks=qtic2,
                xrotation=45,
                fontsize={"xticks": 5, "yticks": 5},
                tools=["tap"],
                shared_axes=False,
                scalebar=True,
                scalebar_range="x",
                scalebar_location="top_left",
                scalebar_unit=("bp"),
            )
        )
    ),
)
pn.panel(p1).servable(title=plotTitle)
