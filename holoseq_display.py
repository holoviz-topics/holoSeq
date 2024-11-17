# see https://github.com/fubar2/holoSeq
# pip install datashader dask[dataframe] holoviews[recommended] pandas matplotlib bokeh
#
# panel serve --address 0.0.0.0 --port 8080 --show --session-token-expiration 9999999 --args --inFile ../hg002_bothHiC.paf_cisH1_hseq.gz 
# This is a generic plot re-creator for holoSeq hseq compressed coordinate data.
# It presents interactive plots at scale using panel.
# Plot coordinates and axis metadata are prepared from genomic data such as PAF or bigwig files,
# using companion converters. They output all the information needed to re-create an interactive plot, as coordinates and axis metadata.
# The intended use case is making interactive specialised browsers for VGP and other genome assemblies and their annotation
# easily browsed and explored using a web browser.
# Annotation files with 10s of millions of data points can be zoomed smoothly down from the whole assembly scale, down to individual points,
# on a laptop.
# Ross Lazarus October 2024

import argparse
from bisect import bisect_left
from collections import OrderedDict
import gzip
import logging
import numpy as np
import os

import holoviews as hv
import panel as pn
import pandas as pd

from holoviews.operation.datashader import (
    rasterize,
    dynspread,
    datashade,
)
from holoviews.operation.element import apply_when
from holoviews.operation.resample import ResampleOperation2D
from holoviews.operation import decimate


logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger("holoseq_display")

from holoviews import opts
from holoviews.operation.datashader import datashade, rasterize, shade, dynspread, spread


hv.extension('bokeh','matplotlib', width=100)

# Default values suitable for this notebook
decimate.max_samples=1000
dynspread.max_px=8
dynspread.threshold=0.75
ResampleOperation2D.width=250
ResampleOperation2D.height=250


# inFile = "galaxy_inputs/paf/bothmap.paf.tab.tabular"
inFile = "/home/ross/rossgit/holoviews-examples/holoSeqtest.gz"
holoSeqHeaders = ["@v1HoloSeq1D", "@v1HoloSeq2D"]
hv.extension("bokeh")
pn.extension()

dynspread.max_px = 7
dynspread.threshold = 0.6


def xportHtml(fname, hObj):
    "save a holoview object to an interactive but not adaptive scaling HTML page"
    hv.save(filename=fname, obj=hObj)


class holoSeq_maker:
    """
    returns panels reconstructed from hseq data - coordinates with axis and other metadata
    """

    def __init__(self, width):
        """ """
        self.pwidth = width

    def xportHtml(self, fname, hObj):
        "save a holoview object to an interactive but not adaptive scaling HTML page"
        hv.save(filename=fname, obj=hObj)

    def import_holoSeq_data(self, inFile):
        """
        reverse process of dumping the data in holoSeq format from a converter
        see https://github.com/fubar2/holoSeq/blob/main/HoloSeqOverview.md
        """
        haps = {}
        hh = []
        xcoords = []
        ycoords = []
        annos = []
        metadata = {}
        gffdata = []
        hsDims = None
        plotType = None
        with gzip.open(inFile, "rt") as f:
            for i, trow in enumerate(f):
                if i == 0:
                    hseqformat = trow.split()[0]
                    if hseqformat not in holoSeqHeaders:
                        log.warn(
                            f"Supplied input {inFile} has first row {trow} so is not a valid holoSeq input file"
                        )
                        log.warn("First row must start with one of these:%s" % holoSeqHeaders)
                        return
                    hsDims = holoSeqHeaders.index(hseqformat) + 1
                    if hsDims == 1:
                        plotType = "bar"
                        if len(hseqformat) > 1:
                            plotType = hseqformat[1].strip()
                elif trow.startswith("@"):
                    row = trow[1:]
                    if row.startswith("@"):
                        srow = row[1:].split()
                        metadata[srow[0]] = srow[1:]
                        if srow[0] == "GFF":
                            isGFF = True
                    else:
                        srow = row.split()
                        if hsDims == 2:
                            if len(srow) >= 3:
                                hap, cname, cstart = srow[:3]
                                if not haps.get(hap, None):
                                    print("adding hap", hap)
                                    hh.append(hap)
                                    haps[hap] = {"cn": [], "startpos": []}
                                haps[hap]["cn"].append(cname.strip())
                                haps[hap]["startpos"].append(int(cstart.strip()))
                            else:
                                log.warn(f"Supplied input {inFile} at row {i} = {row} lacks reference name, contig name and contig length. Not a valid holoSeq input file")
                                return
                        else:
                            if len(srow) >= 2:
                                hap, cname, cstart = srow[:3]
                                if not haps.get(hap, None):
                                    log.debug("adding hap %s" % hap)
                                    hh.append(hap)
                                    haps[hap] = {"cn": [], "startpos": []}
                                haps[hap]["cn"].append(cname)
                                haps[hap]["startpos"].append(int(cstart))
                            else:
                                log.warn(
                                    f"Supplied input {inFile} at row {i} = {row} lacks reference name, contig name and contig length. Not a valid holoSeq input file")
                                return
                else:  # not header row
                    srow = trow.split()
                    lrow = len(srow)
                    if hsDims == 2:
                        if lrow < 2:
                            log.warn(
                                    f"Supplied 2D input {inFile} at row {i} = {trow} - needs at least two coordinates to be a valid holoSeq input file")
                            return
                        else:
                            if srow[0].isdigit() and srow[1].isdigit():
                                xcoords.append(int(srow[0]))
                                ycoords.append(int(srow[1]))
                                if lrow > 2:
                                    annos.append(srow[2:])
                            else:
                                log.warn(
                                    f"Supplied 2D input {inFile} at row {i} = {trow} - needs at least two integer coordinates to be a valid holoSeq input file")
                                return
                    else:
                        if isGFF:
                            gffdata.append(
                                srow
                            )  # mrna XP_026254554.1 SUPER_5 1006698964 1006703995 100 100 + 3226 1006703993

                        else:
                            if srow[0].isdigit():
                                xcoords.append(int(srow[0]))
                                if lrow > 1:
                                    ycoords.append(int(srow[1]))
                                if lrow > 2:
                                    annos.append(srow[2:])
                            else:
                                log.warn(
                                    f"Supplied 1D input {inFile} at row {i} = {trow} - needs at least one integer coordinate to be a valid holoSeq input file")
                                return
        if len(hh) < 2:
            log.debug('extending haps %s' % hh)
            hh.append(hh[0])
        hh.sort()
        return (hsDims, haps, xcoords, ycoords, annos, plotType, metadata, gffdata, hh)

    def makePafPanel(self, inFile, pwidth):
        """
        prepare a complete panel for the final display
        """

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
                styles={
                    "font-size": "10pt",
                    "color": "darkblue",
                    "text-align": "center",
                },
                width=pwidth,
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
                styles={
                    "font-size": "10pt",
                    "color": "darkblue",
                    "text-align": "center",
                },
                width=pwidth,
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
                styles={
                    "font-size": "10pt",
                    "color": "darkblue",
                    "text-align": "center",
                },
                width=pwidth,
            )
            return str_pane

        (hsDims, hapsread, xcoords, ycoords, annos, plotType, metadata, gffdata, hh) = (
            self.import_holoSeq_data(inFile)
        )
        title = " ".join(metadata["title"])
        hqstarts = OrderedDict()
        haps = []
        print("Read nx=", len(xcoords), "ny=", len(ycoords))
        h1starts = []
        h1names = []
        h2starts = []
        h2names = []
        for i, hap in enumerate(hapsread.keys()):
            haps.append(hap)
            hqstarts[hap] = OrderedDict()
            for j, contig in enumerate(hapsread[hap]["cn"]):
                cstart = hapsread[hap]["startpos"][j]
                hqstarts[hap][contig] = cstart
                if i == 0:
                    h1starts.append(cstart)
                    h1names.append(contig)
                else:
                    h2starts.append(cstart)
                    h2names.append(contig)
        hap = hh[0]
        if (len(h2starts) == 0):
            h2starts = h1starts
            h2names = h1names
            log.warn('only one haplotype read for %s' % title)
        qtic1 = [(h1starts[i], h1names[i]) for i in range(len(h1starts))]
        hap = hh[1]
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
        xcf = os.path.splitext(metadata["xclenfile"][0])[0]
        ycf = 'Y:' + os.path.splitext(metadata["yclenfile"][0])[0]
        print("xcf", xcf, "ycf", ycf)
        pafxy = pd.DataFrame.from_dict({xcf: xcoords, ycf: ycoords})
        pafp = hv.Points(pafxy, kdims=[xcf, ycf])

        # apply_when(pafp, operation=rasterize, predicate=lambda x: len(x) > 5000)
        stream = hv.streams.Tap(x=0, y=0)
        ax = metadata.get('axes', [None])[0]
        log.debug('axes = %s' % ax)
        if ax == "BOTH":
            showloc = pn.bind(showTrans, x=stream.param.x, y=stream.param.y)
        elif ax == haps[0]:
            showloc = pn.bind(showH1, x=stream.param.x, y=stream.param.y)
        elif ax == haps[1]:
            showloc = pn.bind(showH2, x=stream.param.x, y=stream.param.y)
        else:
            log.warn('ax = %s for title = %s - cannot assign axes' % (ax, title))
            showloc = pn.bind(showTrans, x=stream.param.x, y=stream.param.y)
        p1 = pn.Column(
            showloc,
            pn.pane.HoloViews(
                dynspread(rasterize(pafp), streams=[stream])
                .relabel("%s" % title)
                .opts(
                    shared_axes=False,
                    cmap="inferno",
                    cnorm="log",
                    colorbar=True,
                    width=self.pwidth,
                    height=self.pwidth,
                    xticks=qtic1,
                    yticks=qtic2,
                    xrotation=45,
                    fontsize={"xticks": 5, "yticks": 5},
                    tools=["tap"],
                    scalebar=True,
                    scalebar_range="x",
                    scalebar_location="top_left",
                    scalebar_unit=("bp"),
                    show_grid=True,
                )
            ),
        )
        return p1, title

    def makeBWPanel(self, inFile, pwidth):
        """
        prepare a complete panel for the final display
        """

        def showX(x, y):
            if np.isnan(x):
                s = "Mouse click on image for location"
            else:
                i = bisect_left(h1starts, x)
                chrx = h1names[i - 1]
                offsx = x - h1starts[i - 1]
                s = "%s:%d" % (chrx, offsx)
            str_pane = pn.pane.Str(
                s,
                styles={
                    "font-size": "10pt",
                    "color": "darkblue",
                    "text-align": "center",
                },
                width=pwidth,
            )
            return str_pane

        (hsDims, hapsread, xcoords, ycoords, annos, plotType, metadata, gffdata, hh) = (
            self.import_holoSeq_data(inFile)
        )
        title = " ".join(metadata["title"])
        haps = []
        print("Read nx=", len(xcoords), "ny=", len(ycoords))
        h1starts = []
        h1names = []
        for i, hap in enumerate(hapsread.keys()):
            haps.append(hap)
            for j, contig in enumerate(hapsread[hap]["cn"]):
                cstart = hapsread[hap]["startpos"][j]
                h1starts.append(cstart)
                h1names.append(contig)
        hap = haps[0]
        log.debug("h1names=%s" % h1names[:20])
        # qtic1 = [(hqstarts[hap][x], x) for x in hqstarts[hap].keys()]
        qtic1 = [(h1starts[i], h1names[i]) for i in range(len(h1starts))]
        log.debug("qtic1=%s" % qtic1[:20])
        xax = metadata["xclenfile"][0]
        yax = metadata["yclenfile"][0] + "Bigwig value"
        pafxy = pd.DataFrame.from_dict({xax: xcoords, yax: ycoords})
        taps = hv.streams.Tap(x=0, y=0)
        showloc = pn.bind(showX, x=taps.param.x, y=taps.param.y)
        bigw = pn.pane.HoloViews(
            decimate(hv.Curve(pafxy), streams=[taps])
            .opts(interpolation="steps-pre", color="darkblue")
            .relabel("%s" % title)
            .opts(
                width=pwidth,
                height=300,
                xticks=qtic1,
                xrotation=45,
                fontsize={"xticks": 8, "yticks": 10},
                scalebar=True,
                scalebar_range="x",
                scalebar_location="top_left",
                scalebar_unit=("bp"),
                show_grid=True,
                ylim=(-0.1, 200),
                tools=[
                    "xwheel_zoom",
                    "tap",
                    "xpan",
                    "reset",
                ],
                default_tools=[],
                active_tools=["xwheel_zoom", "tap", "pan"],
            )
        )

        p1 = pn.Column(showloc, bigw)

        return p1, title

    def makeGFFPanel(self, inFile, pwidth):
        """
        prepare a complete panel for the final display
https://www.ncbi.nlm.nih.gov/gene/?term=XP_026235740.1

import urllib.request
xpuri = 'https://www.ncbi.nlm.nih.gov/gene/?term=XP_026235740.1'
req = urllib.request.Request(xpuri)
with urllib.request.urlopen(req) as response:
   apage = response.read()
escaped_html = html.escape(apage)

# Create iframe embedding the escaped HTML and display it
iframe_html = f'<iframe srcdoc="{escaped_html}" style="height:100%; width:100%" frameborder="0"></iframe>'

# Display iframe in a Panel HTML pane
pn.pane.HTML(iframe_html, height=350, sizing_mode="stretch_width")
        """

        def showX(x, y):
            if np.isnan(x):
                s = "Mouse click on image for location"
            else:
                i = bisect_left(h1starts, x)
                chrx = h1names[i - 1]
                offsx = x - h1starts[i - 1]
                s = "%s:%d" % (chrx, offsx)
                xi = bisect_left(segs[xcf], x)
                xtarget = segs['target'][xi]
                s += ' x %s' % (xtarget )

            str_pane = pn.pane.Str(
                s,
                styles={
                    "font-size": "10pt",
                    "color": "darkblue",
                    "text-align": "center",
                },
                width=pwidth,
                
            )
            return str_pane

        (hsDims, hapsread, xcoords, ycoords, annos, plotType, metadata, gffdata, hh) = self.import_holoSeq_data(inFile)
        xcf = os.path.splitext(metadata["xclenfile"][0])[0]
        segs = {
            xcf: [],
            "x2": [],
            "wy1": [],
            "y2": [],
            "target": [],
            "colour": [],
            "thickness": [],
            "alpha": []
        }
        """
cds XP_026238700.1 1401967516 1401967635 100 100 - 204
mrna XP_026248570.1 SUPER_3H1 531341254 531595863 100 100 + 1102 -1
cds XP_026248570.1 531341254 531341334 100 100 + 134
        """
        mthick = 3
        cdthick = 50
        for i, rows in enumerate(gffdata):
            if rows[0].lower() == "mrna":
                (kind, targ,  contig, startp, endp, y1, y2, strand, score) = rows[:10]
                startp = int(startp)
                endp = int(endp)
                y1 = int(y1)
                y2 = int(y2)
                colr = "blue"
                if strand == "-":
                    colr = "maroon"        
                segs["target"].append(targ)
                segs[xcf].append(startp)
                segs["x2"].append(endp)
                segs["wy1"].append(y1)
                segs["y2"].append(y2)
                segs["colour"].append(colr)
                segs["thickness"].append(mthick)
                segs["alpha"].append(1.0)
                #segs["stopc"].append(stopc)
            elif rows[0].lower() == "cds": #  f"cds {targ} {con} {startp} {endp} {y} {y} {strand} {score}\n"
                (kind, targ, contig, startp, endp , y1, y2, strand, score) = rows[:10]
                startp = int(startp)
                endp = int(endp)
                y = int(y1)
                colr = "blue"
                if strand == "-":
                    colr = "maroon"        
                segs["target"].append(targ)
                segs[xcf].append(startp)
                segs["x2"].append(endp)
                segs["wy1"].append(y)
                segs["y2"].append(y)
                segs["colour"].append(colr)
                segs["thickness"].append(cdthick)
                segs["alpha"].append(1.0)
        xmin = min(segs[xcf])
        xmax = max(segs["x2"])
        ymin = min(segs["wy1"])
        ymax = max(segs["wy1"])
        title = " ".join(metadata["title"])
        haps = []
        print("GFF rows read =", len(gffdata))
        h1starts = []
        h1names = []
        qtic1 = []
        for i, hap in enumerate(hapsread.keys()):
            haps.append(hap)
            for j, contig in enumerate(hapsread[hap]["cn"]):
                cstart = hapsread[hap]["startpos"][j]
                h1starts.append(cstart)
                h1names.append(contig)
                qtic1.append((cstart, contig))
        hap = haps[0]
        # print("h1names=", h1names[:20])
        # qtic1 = [(hqstarts[hap][x], x) for x in hqstarts[hap].keys()]
        # print("qtic1=", qtic1[:20])
        gffp = hv.Segments(
            segs, [xcf, 'wy1', 'x2', 'y2'], vdims=["target", "colour", "thickness", "alpha"]
            )

        gffp.opts(title="title",
            color="colour",
            line_width="thickness",
            alpha="alpha",
            width=pwidth,
            height=300,
            xticks=qtic1,
            xrotation=45,
            scalebar=True,
            scalebar_range="x",
            scalebar_location="bottom_left",
            scalebar_unit=("bp"),
            fontsize={"xticks": 8, "yticks": 10},
            show_grid=True,
            autorange='y',
            tools=[
                "xwheel_zoom",
                "box_zoom",
                "tap",
                "xpan",
                "reset",
            ],
            default_tools=[],
            active_tools=["xwheel_zoom", "tap", "xpan"],
            shared_axes=True
            )
        apply_when(gffp, operation=rasterize, predicate=lambda x: len(x) > 5000)
        taps = hv.streams.Tap(source=gffp, x=0, y=0)
        showloc = pn.bind(showX, x=taps.param.x, y=taps.param.y)
        gp = pn.pane.HoloViews(gffp)
        p = pn.Column(showloc, gp)

        return p, title


parser = argparse.ArgumentParser(description="", epilog="")
parser.add_argument(
    "--inFile",
    help="gzipped hseq coordinates and contigs",
    default="mUroPar1_cis1.hseq.gz",
    nargs="+",
)
parser.add_argument(
    "--size", help="Display size in pixels. Default is 800", default=1000
)
parser.add_argument("--version", "-V", action="version", version="0.1")
args = parser.parse_args()
pwidth = int(args.size)
hsm = holoSeq_maker(pwidth)
for i, infile in enumerate(args.inFile):
    print("Infile = ", infile)
    if "bw.hseq.gz" in infile:
        p1, title = hsm.makeBWPanel(infile, pwidth)
    elif "gff.hseq.gz" in infile:
        p1, title = hsm.makeGFFPanel(infile, pwidth)
    else:
        p1, title = hsm.makePafPanel(infile, pwidth)
    if i == 0:
        outp = p1
    else:
        outp = outp + p1
pn.Row(outp).servable(title=title)
