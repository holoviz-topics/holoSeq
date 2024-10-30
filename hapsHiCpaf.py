# RUN THIS to create three interactive HiC displays in a Jypyter notebook
# after loading dependencies with
# ! pip install datashader dask[dataframe] holoviews[recommended] pandas matplotlib bokeh
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
from functools import cmp_to_key
import math
import numpy as np

import holoviews as hv
import panel as pn
import pandas as pd

from holoviews.operation.datashader import (
    rasterize,
    dynspread,
)
from holoviews.operation.resample import ResampleOperation2D
from holoviews.operation import decimate

useDecimate = False  # will rasterize instead
# inFile = "galaxy_inputs/paf/bothmap.paf.tab.tabular"
inFile = "/home/ross/rossgit/holoviews-examples/huge.paf"
ptwidth = 1000
pcwidth = 800
# width settings for plots and location bars
hv.extension("bokeh")
pn.extension()


def rotatecoords(x, y, radians=0.7853981633974483, origin=(0, 0)):
    # https://gist.github.com/LyleScott/d17e9d314fbe6fc29767d8c5c029c362
    offset_x, offset_y = origin
    adjusted_x = x - offset_x
    adjusted_y = y - offset_y
    cos_rad = math.cos(radians)
    sin_rad = math.sin(radians)
    qx = offset_x + cos_rad * adjusted_x + sin_rad * adjusted_y
    qy = offset_y + -sin_rad * adjusted_x + cos_rad * adjusted_y
    return qx, qy


def getHap(contig):
    """
    function to return suffix H1 from chrH1 - adjust to suit.
    """
    return contig[-2:]


def xportHtml():
    "save a holoview object to an interactive but not adaptive scaling HTML page"
    p = (
        rasterize(pafcis1)
        .relabel("%s Cis HiC interactions" % haps[0])
        .opts(
            cmap="inferno",
            cnorm="log",
            colorbar=True,
            width=1000,
            height=600,
            xticks=qtic1,
            yticks=qtic1,
            xrotation=45,
            fontsize={"xticks": 6, "yticks": 6},
        )
    )
    hv_pane = pn.pane.HoloViews(p, height=2000, width=2000)
    hv.save(filename="H1cis.html", obj=hv_pane)


def export_mapping():

    for ofn, d, hstarts in [
        ("cis1.tab", cis1, hqstarts[0]),
        ("cis2.tab", cis2, hqstarts[1]),
        ("trans1.tab", trans1, None),
    ]:
        with open(ofn, "w") as f:
            for i in range(len(d["x"])):
                row = "%d\t%d\n" % (d["x"][i], d["y"][i])
                f.write(row)
        if hstarts is not None:
            with open("%s.len" % ofn, "w") as f2:
                clens = ["%s\t%d\n" % (x, hstarts[x]) for x in inhstarts.keys()]
                f2.write("".join(clens))


def sorthapqname(s1, s2):
    """
    Ultimately futile unless there are consistent contig  naming rules and there are not - there is super_2 and super13_unloc_5 FFS.
    Best if the assemblers and curators named things the way they want them to be sorted.
    fugly hack to sort super contigs before anything else
    then by contig number or if they are the same offset
    ('SUPER_2H1', 226668729) , ('SUPER_1H1', 284260672), ('SUPER13_unloc_5H1',..), (Scaffold_1116H2, ...)
    """
    if s1[0] == s2[0]:  # simplest case - same contig, sort on offset
        return s1[1] - s2[1]  # neg if left sorts before
    s11, s12 = s1[0].split("_", 1)
    s1n = s12.split("_")[-1][:-2]
    s1_super = (s11.upper() == "SUPER") or (s11.upper().startswith("CHR"))
    s21, s22 = s2[0].split("_", 1)
    s2n = s22.split("_")[-1][:-2]
    s2_super = (s21.upper() == "SUPER") or (s21.upper().startswith("CHR"))
    if s1n.isdigit():
        s1nn = int(s1n)
    else:
        s1nn = ord(s1n[0]) * 1000
    if s2n.isdigit():
        s2nn = int(s2n)
    else:
        s2nn = ord(s2n[0]) * 1000
    if s1_super == s2_super:
        nunder1 = len(s1[0].split("_"))
        nunder2 = len(s2[0].split("_"))  # _unloc or whatever
        if nunder1 == nunder2:
            return s1nn - s2nn
        else:
            return nunder1 - nunder2
    elif s1_super:
        return -1
    elif s2_super:
        return 1
    else:
        return s1nn - s2nn


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
        width=ptwidth,
    )
    return str_pane


# Default values suitable for this notebook
decimate.max_samples = 10000
dynspread.max_px = 10
dynspread.threshold = 0.6
ResampleOperation2D.width = 1500
ResampleOperation2D.height = 1500
# need to convert the categorical contigs into a sequence for holoviews to munch
# use the contig length from the paf to figure out the cumulative start for each contig
# contigs are length ordered - that does not always work well when the haplotypes differ widely

hpstarts = OrderedDict()
hqstarts = OrderedDict()
hlens = {}
haps = []
hpsorts = {}
hqsorts = {}
with open(inFile, "r") as f:
    for i, rows in enumerate(f):
        row = rows.split()
        c1 = row[0]
        c2 = row[5]
        hap = getHap(c1)  # define your own function - this caters to a suffix H1 or H2
        if hap not in haps:
            haps.append(hap)
            hlens[hap] = {}
            hpsorts[hap] = []
            hqsorts[hap] = []
        if not hlens[hap].get(c1, None):
            hlens[hap][c1] = int(row[1])
            hpsorts[hap].append((int(row[1]), c1))  # length sort gives horrible results
            hqsorts[hap].append(
                (c1, int(row[1]))
            )  # query name sort with a special sort function is better
for hap in haps:
    cum = 1
    hpsorts[hap].sort(reverse=True)
    hqsorts[hap].sort(key=cmp_to_key(sorthapqname))
    hpstarts[hap] = OrderedDict()
    hqstarts[hap] = OrderedDict()
    for clen, contig in hpsorts[hap]:
        hpstarts[hap][contig] = cum
        cum += clen
    cum = 1
    for contig, clen in hqsorts[hap]:
        hqstarts[hap][contig] = cum
        cum += clen
h1starts = [hqstarts[haps[0]][x] for x in hqstarts[haps[0]].keys()]
h1names = list(hqstarts[haps[0]].keys())
h2starts = [hqstarts[haps[1]][x] for x in hqstarts[haps[1]].keys()]
h2names = list(hqstarts[haps[1]].keys())
# have the axes set up so prepare the three plot x/y vectors
# for a second pass to calculate all the coordinates.
cis1 = {"x": [], "y": []}
cis2 = {"x": [], "y": []}
trans1 = {"x": [], "y": []}
with open(inFile, "r") as f:
    for rows in f:
        row = rows.split()
        c1 = row[0]
        c2 = row[5]
        H1 = getHap(c1)
        H2 = getHap(c2)
        if H1 != H2:  # trans
            if H1 == haps[0]:  # x is h1 for trans - otherwise ignore
                trans1["x"].append(hqstarts[H1][c1] + int(row[2]))
                trans1["y"].append(hqstarts[H2][c2] + int(row[7]))
            else:
                trans1["y"].append(hqstarts[H1][c1] + int(row[2]))
                trans1["x"].append(hqstarts[H2][c2] + int(row[7]))
        else:  # cis
            if H1 == haps[0]:
                cis1["x"].append(hqstarts[H1][c1] + int(row[2]))
                cis1["y"].append(hqstarts[H2][c2] + int(row[7]))
            else:
                cis2["x"].append(hqstarts[H1][c1] + int(row[2]))
                cis2["y"].append(hqstarts[H2][c2] + int(row[7]))
hap = haps[0]
qtic1 = [(hqstarts[hap][x], x) for x in hqstarts[hap].keys()]
hap = haps[1]
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
pafxycis1 = pd.DataFrame.from_dict(cis1)
pafxycis2 = pd.DataFrame.from_dict(cis2)
pafxytrans = pd.DataFrame.from_dict(trans1)
pafcis1 = hv.Points(pafxycis1)
pafcis2 = hv.Points(pafxycis2)
paftrans = hv.Points(pafxytrans)
streamcis1 = hv.streams.Tap(source=pafcis1, x=np.nan, y=np.nan)
streamcis2 = hv.streams.Tap(source=pafcis2, x=np.nan, y=np.nan)
streamtrans = hv.streams.Tap(source=paftrans, x=np.nan, y=np.nan)
showloc1 = pn.bind(showH1, x=streamcis1.param.x, y=streamcis1.param.y)
showloc2 = pn.bind(showH2, x=streamcis2.param.x, y=streamcis2.param.y)
showloctrans = pn.bind(showTrans, x=streamtrans.param.x, y=streamtrans.param.y)
hap = haps[0]
qtic1 = [(hqstarts[hap][x], x) for x in hqstarts[hap].keys()]
hap = haps[1]
qtic2 = [(hqstarts[hap][x], x) for x in hqstarts[hap].keys()]
# to rotate so the diagonal becomes the x axis
# xcis1r, ycis1r = rotatecoords(xcis1, ycis1, radians=0.7853981633974483, origin=(max(xcis1),max(ycis1)))
# pafxycis1 = pd.DataFrame(np.vstack([xcis1r,ycis1r]).T, columns = ['x', 'y'])
# bisect.bisect_left(a, x, lo=0, hi=len(a), *, key=None)
pafxycis1 = pd.DataFrame.from_dict(cis1)
pafxycis2 = pd.DataFrame.from_dict(cis2)
pafxytrans = pd.DataFrame.from_dict(trans1)
pafcis1 = hv.Points(pafxycis1)
pafcis2 = hv.Points(pafxycis2)
paftrans = hv.Points(pafxytrans)
streamcis1 = hv.streams.Tap(source=pafcis1, x=np.nan, y=np.nan)
streamcis2 = hv.streams.Tap(source=pafcis2, x=np.nan, y=np.nan)
streamtrans = hv.streams.Tap(source=paftrans, x=np.nan, y=np.nan)
showloc1 = pn.bind(showH1, x=streamcis1.param.x, y=streamcis1.param.y)
showloc2 = pn.bind(showH2, x=streamcis2.param.x, y=streamcis2.param.y)
showloctrans = pn.bind(showTrans, x=streamtrans.param.x, y=streamtrans.param.y)
# prepare and show the 3 plots
p1 = pn.Column(
    showloc1,
    pn.pane.HoloViews(
        dynspread(
            rasterize(pafcis1)
            .relabel("%s Cis HiC interactions" % haps[0])
            .opts(
                cmap="inferno",
                cnorm="log",
                colorbar=True,
                width=pcwidth,
                height=pcwidth,
                xticks=qtic1,
                yticks=qtic1,
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

p2 = pn.Column(
    showloc2,
    pn.pane.HoloViews(
        dynspread(
            rasterize(pafcis2)
            .relabel("%s Cis HiC interactions" % haps[1])
            .opts(
                cmap="inferno",
                cnorm="log",
                colorbar=True,
                width=pcwidth,
                height=pcwidth,
                xticks=qtic2,
                yticks=qtic2,
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

p12 = pn.Row(p1, p2)

p3 = pn.Row(
    pn.Column(
        showloctrans,
        pn.pane.HoloViews(
            dynspread(
                rasterize(paftrans)
                .relabel("%s/%s Trans HiC interactions" % (haps[0], haps[1]))
                .opts(
                    cmap="inferno",
                    cnorm="log",
                    colorbar=True,
                    width=ptwidth,
                    height=ptwidth,
                    xticks=qtic1,
                    yticks=qtic2,
                    xrotation=45,
                    shared_axes=False,
                    fontsize={"xticks": 5, "yticks": 5},
                    tools=["tap"],
                    scalebar=True,
                    scalebar_range="x",
                    scalebar_location="top_left",
                    scalebar_unit=("bp"),
                 show_grid=True,
                )
            )
        ),
    )
)

pn.Column(p12, p3).servable(title=inFile)
