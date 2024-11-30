# for Mashmap paf, python holoSeq_prepare_paf.py --inFile  hg002_2k99.paf --title "hg002 Mashmap" --hap_indicator None --contig_sort length
# for HiC pairs
# python holoSeq_prepare_paf.py --inFile mUroPar1H1H2.paf --xclenfile mUroPar1H1suffix.len --yclenfile mUroPar1H2suffix.len --contig_sort VGPname --hap_indicator Suffix --title "VGP mUroPar1 HiC data"
# panel serve holoseq_display.py --show --args --inFile mUroPar1H1H2.paf_cisH1_hseq.gz mUroPar1H1H2.paf_cisH2_hseq.gz mUroPar1H1H2.paf_trans_hseq.gz  --size 1000
# or
# panel serve holoseq_display.py --show --args --inFile mUroPar1H1H2.paf_cisH1_hseq.gz --size 1000
#
# python holoSeq_prepare_paf.py --inFile mUroPar1_protH1.gff --xclenfile mUroPar1H1suffix.len --contig_sort VGPname --title "mUroPar1 NCBI protein GFF"
#
# python holoSeq_prepare_paf.py --inFile ../hg002_bothHiC.paf --xclenfile hg002H1_suffixed.len --yclenfile hg002H2_suffixed.len --contig_sort VGPname --hap_indicator Suffix --title "T2T HG002 HiC data"
#
# Proof of cocept data are Arima HiC reads from the Arctic Ground Squirrel mUroPar1 VGP genomeArk repository
# processed with Dephine's Pretext workflow using Bellerophon to remove chimeric reads
# The paired bam output is converted to PAF with an awk script (!) and that's what is read
# in this code.
# The pairs are parsed to extract the haplotype designator from the contig name
# typically a suffix like H1 extracted in getHap - rewrite that to suit your names.
# Contig ordering really matters for the plots to make any sense.
# Ideally, curators name them so they sort alphanumerically without effort.
# There's a sorthapqname function that is used here. Designed for the VGP data
# Will need to be replaced for other naming conventions.
# It's a mess.
# Sorting by contig name is based on VGP conventions - SUPER_ first, then scaffolds
# One problem to watch out for is that any differences in ordering of the X and Y contigs can make all sorts of
# artifacts appear such as the kaleidoscopic patterns seen in Pretextviewer.

# Ross Lazarus October 2024

import argparse

from collections import OrderedDict
from functools import cmp_to_key
from pathlib import Path

import gzip
import io
import itertools
import logging
import math
import numpy as np

import re
import os

import pandas as pd
import pybigtools

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger("holoseq_prepare")

# inFile = "galaxy_inputs/paf/bothmap.paf.tab.tabular"
inFile = "/home/ross/rossgit/holoviews-examples/huge.paf"

holoSeqHeaders = ["@v1HoloSeq1D", "@v1HoloSeq2D"]


def getHap(contig):
    """
    function to return suffix H1 from "chrH1" - adjust to suit the specific conventions used by
    whoever assembles your genomes.
     help="None, Suffix (H[1,2]) Dashsuffix (_H...)"
    """
    hi = args.hap_indicator
    if hi == "None":
        return "H1"
    elif hi == "Suffix":
        return contig[-2:]
    elif hi == "Dashsuffix":
        return contig.split("_")[-1]


def getContigs(lenFile):
    # samtools faidx will make one of these from a genome fasta
    # whitespace delimited contig names and lengths.
    contigs = []
    seen = {}
    haps = []
    with open(lenFile) as lf:
        for i, row in enumerate(lf):
            row = [x.strip() for x in row.strip().split()]
            if len(row) > 1:
                c, clen = row[:2]
                if seen.get(c, None):
                    log.debug("Contig %s seen again at row %d of %s" % (c, i, lenFile))
                else:
                    seen[c] = c
                    contigs.append((c, int(clen)))
                h = getHap(c)
                if h not in haps:
                    haps.append(h)
    return contigs, haps


def VGPsortfunc(s1, s2):
    """
    # big fugly hack to sort super contigs before anything else
    # then by contig number or if they are the same offset
    # ('SUPER_2H1', 226668729) , ('SUPER_1H1', 284260672), ('SUPER13_unloc_5H1',..), (Scaffold_aa16H2, ...)
    # or chr10H1 etc
    # always work with uppercase
    # ^([^_]+)_(\S+)H[12]{1}$ gives group1 SUPER or Scaffold and group2 X or 22
    #     ^([^_]+)_([^_]+)_unloc_(\S+)H[12]{1}$ gives super aa x for super_aa_unloc_XH2
    """
    ssc = re.compile(r"^(CHR)_*(\d+|[a-zA-Z0-9]+)_*(\S*)$")
    # should match chrY or chr_y_paternal or chr333
    ss1 = re.compile(
        r"^([^_]+)_(\S+)H[12]{1}$"
    )  # should match VGP non-unloc - super/chr/scaffold. Always uppercase at entry.
    ss2 = re.compile(r"^([^_]+)_([^_]+)_UNLOC_(\S+)H[12]{1}$")  # for the unloc

    def intOrd(n):
        "x, y, z, w typically sex chromosomes - there be some weird beasts"
        if n:
            n = n.replace("_", "").replace("chr", "")  # in case chr_aa or something
            if n.isdigit():
                return int(n)
            else:
                return ord(n[0])
        else:
            return None

    def matchme(s):
        c1 = n1 = n2 = None
        found = ssc.search(s)
        if found:
            c1, n1, n2 = found.groups()[:3]
        else:
            found = ss2.search(s)
            if found:
                c1, n1, n2 = found.groups()[:3]
            else:
                found = ss1.search(s)
                if found:
                    c1, n1 = found.groups()[:2]
                    n2 = None
        if n1.isdigit():
            n1 = int(n1)
        else:
            n1 = ord(n1[0])
        if n2:
            if n2.isdigit():
                n2 = int(n2)
            else:
                n2 = ord(n2[0])
        return c1, n1, n2

    s1 = [s1[0].upper(), s1[1]]
    s2 = [s2[0].upper(), s2[1]]
    if s1[0] == s2[0]:  # simplest case - same contig, sort on offset
        return s1[1] - s2[1]  # neg if left sorts before
    u1 = "UNLOC" in s1[0]
    u2 = "UNLOC" in s2[0]
    if u1 and not u2:
        return 1  # u1 unloc goes after
    elif u2 and not u1:
        return -1  # u1 goes before unloc
    isSuper1 = (not u1) and (("SUPER" in s1[0]) or ("CHR" in s1[0]))
    isSuper2 = (not u2) and (("SUPER" in s2[0]) or ("CHR" in s2[0]))
    isScaff1 = (not u1) and (("SCAFFOLD" in s1[0]))
    isScaff2 = (not u2) and (("SCAFFOLD" in s2[0]))
    if isSuper1 and not isSuper2:
        return -1
    elif isSuper2 and not isSuper1:
        return 1
    # Must parse

    c1, naa, nab = matchme(s1[0].upper())
    c2, nba, nbb = matchme(s2[0].upper())
    if not c1 or not c2:
        log.debug("no match for ss1 %s and/or ss2 %s" % (ss1, ss2))
        return 0
    else:
        if isSuper1 or (
            isScaff1 and isScaff2
        ):  # if a super must both be supers or if both are scaffolds
            return naa - nba
        else:  # must both be unlocs
            if naa == nba:
                return nab - nbb
            else:
                return naa - nba


def Lengthsortfunc(s1, s2):
    """ """
    return s1[1] - s2[1]  # neg if left sorts before


def contsort(contigs, args):
    # sort and return offsets to starts of each contig
    # hstarts = list(itertools.accumulate(hlens))
    if args.contig_sort.lower() == "vgpname":
        contigs.sort(key=cmp_to_key(VGPsortfunc))
    elif args.contig_sort.lower() == "name":
        contigs.sort()
    elif args.contig_sort.lower() == "length":
        contigs.sort(key=cmp_to_key(Lengthsortfunc), reverse=True)
    clens = [x[1] for x in contigs]
    cnames = [x[0] for x in contigs]
    cstarts = list(itertools.accumulate(clens))
    cstarts.insert(0, 0)  # first one starts at 0
    scont = OrderedDict(zip(cnames, cstarts))
    maxpos = cstarts[-1] + clens[-1]
    return scont, maxpos

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
        log.debug("xmin %d xnew %d ymin %d ynew %d, cos %f, sin %f" % (self.xmin, self.xnew, self.ymin, self.ynew, self.cos_rad, self.sin_rad))
        


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


class gffConvert:
    """
            Only care about mRNA cds and stop codons initally. Turn into segments. Filter so only data in contigs is retained from input.
    SUPER_1 miniprot        mRNA    139006290       139072696       22660   -       .       ID=MP000006;Rank=1;Identity=0.9979;Positive=0.9984;Target=XP_026244093.1 1 4350
        fix positions as we go with a lookup contig -> cumulated offset
        can pop open https://www.ncbi.nlm.nih.gov/protein/XP_026244093.1
    """

    def __init__(self, gff, outFname, contigs, args):
        mrnaseen = {}
        segs = {}
        comment = "#"
        self.hsId = "@v1HoloSeq1D"
        self.inFname = gff
        log.debug("contigs=%s" % str(contigs)[:1000])
        with open(gff) as g:
            for i, row in enumerate(g):
                if not row.startswith(comment):
                    (id, name, kind, startp, endp, score, strand, phase, text) = [
                        x.strip() for x in row.split()[:9]
                    ]
                    if not segs.get(id, None):
                        segs[id] = []
                    if kind.lower() in ["cds", "mrna"]:
                        anno = text.split(";")
                        tanno = [
                            x.strip()[7:]
                            for x in anno
                            if x.lower().startswith("target=")
                        ]
                        target = tanno[0]
                    startp = int(startp)
                    endp = int(endp)
                    offset = contigs.get(id, -1)
                    if args.addH1 and offset < 0:
                        offset = contigs.get(id + "H1", -1)
                        id = id + "H1"
                    if offset < 0:
                        log.warn(
                            "Ignored gff3 id %s missing from supplied xcontigs, in row %d %s of %s with addH1=%s"
                            % (id, i, row, gff, args.addH1)
                        )
                    else:
                        if kind.lower() == "mrna":
                            if target:
                                if mrnaseen.get(target, None):
                                    log.debug(
                                        "Seeing mrna target %s again at row %d"
                                        % (target, i)
                                    )
                                else:
                                    mrnaseen[target] = target
                                segs[id].append(
                                    (
                                        startp + offset,
                                        endp + offset,
                                        strand,
                                        score,
                                        target,
                                        "mrna",
                                    )
                                )
                            else:
                                log.warn("no target found in %s at row %d" % (text, i))
                        elif kind.lower() == "stop_codon":
                            segs[id].append((startp + offset, target, "stopc"))
                        elif kind.lower() == "cds":
                            segs[id].append(
                                (
                                    startp + offset,
                                    endp + offset,
                                    strand,
                                    score,
                                    target,
                                    "cds",
                                )
                            )

        self.export_mapping(outFname, contigs, segs, args)

    def export_mapping(self, outFname, contigs, segs, args):
        """
        for GFF
        @v1HoloSeq2D for example
        A default  Y value of 100 is set for each mRNA's extent, but often there are dozens of different named sequences in the databases that will overlap.
        Assuming the GFF is sorted, overlapping mRNA is identified as ending or starting in the previous extent and the y value is decremented to avoid overlap
        Y is reset if no overlap
        """

        def prepHeader(contigs, args):
            """
            holoSeq output format
            """
            h = ["@%s %s %d" % (getHap(k), k, contigs[k]) for k in contigs.keys()]
            metah = [
                self.hsId,
                "@@GFF 1",
                "@@title %s" % args.title,
                "@@datasource GFF",
                "@@datafile %s" % self.inFname,
                "@@refURI %s" % args.refURI,
                "@@xclenfile %s" % args.xclenfile,
            ]

            return metah + h

        def ranges_overlap(x1, x2, y1, y2):
            # buried in https://stackoverflow.com/questions/6821156/how-to-find-range-overlap-in-python
            if x1 == x2 or y1 == y2:
                return False
            return x1 <= y2 and y1 <= x2

        hdr = prepHeader(contigs, args)

        with gzip.open(outFname, mode="wb") as ofn:
            ofn.write(str.encode("\n".join(hdr) + "\n"))
            y = 100
            for con in contigs.keys():
                lastseg = ("", 0, 0)
                subs = segs.get(con, [])
                if len(subs) > 0:
                    subs.sort(key=lambda x: x[0])
                    for i, m in enumerate(subs):
                        kind = m[-1]
                        if kind == "mrna":
                            (startp, endp, strand, score, targ, _) = m
                            if ranges_overlap(lastseg[1], lastseg[2], startp, endp):
                                y -= 1
                            else:
                                y = 100
                                lastseg = (targ, startp, endp)
                            row = str.encode(
                                f"mrna {targ} {con} {startp} {endp} {y} {y} {strand} {score}\n"
                            )
                            ofn.write(row)
                        elif kind == "cds":
                            (startp, endp, strand, score, targ, _) = m
                            row = str.encode(
                                f"cds {targ} {con} {startp} {endp} {y} {y} {strand} {score}\n"
                            )
                            ofn.write(row)
                        elif kind == "stopc":
                            (startp, targ, _) = m
                            row = str.encode(f"stopc {targ} {con} {startp}\n")
                            ofn.write(row)


class bwConvert:
    """
    points = [(50*i, 100+random.random()) for i in range(10000)]
    hv.Curve(points).opts(interpolation='steps-post').opts(width=1000)
    bigwig to barchart - will autoscale to line?
    """

    def __init__(self, inFname, outFname, args, contigs):
        self.inFname = inFname
        self.hsId = "@v1HoloSeq2D"
        fakepath = "in.bw"
        if os.path.isfile(fakepath):
            os.remove(fakepath)
        p = Path(fakepath)
        p.symlink_to(inFname)  # required by pybigtools (!)
        bwf = pybigtools.open(fakepath)
        bchrlist = bwf.chroms()
        bwchrs = list(bchrlist.keys())
        data = {}
        for i, bchr in enumerate(bwchrs):
            cchr = bchr
            if (not contigs.get(bchr, None)) and args.addH1:
                cchr = cchr + "H1"
            if contigs.get(cchr, None):
                cstart = contigs[cchr]
                data[cchr] = {}
                bw = bwf.records(bchr)
                # Return the records of a given range on a chromosome. The result is an iterator of tuples. For BigWigs, these tuples are in the format (start: int, end: int, value: float).
                data[cchr]["xstart"] = [x[0] + cstart for x in bw]
                bw = bwf.records(bchr)
                data[cchr]["xend"] = [x[1] + cstart for x in bw]
                bw = bwf.records(bchr)
                data[cchr]["xval"] = [x[2] for x in bw]
            else:
                log.warn(
                    "Bigwig contig %s not found in supplied X axis lengths file" % cchr
                )
        self.export_mapping(outFname, contigs, data, args)

    def export_mapping(self, outFname, contigs, data, args):
        """
        for bigwig
        @v1HoloSeq2D for example
        """

        def prepHeader(contigs, args):
            """
            holoSeq output format
            """
            h = ["@%s %s %d" % (getHap(k), k, contigs[k]) for k in contigs.keys()]
            metah = [
                self.hsId,
                "@@bigwig 1",
                "@@title %s" % args.title,
                "@@datasource %s" % "bigwig",
                "@@datafile %s" % self.inFname,
                "@@refURI %s" % args.refURI,
                "@@xclenfile %s" % args.xclenfile,
            ]

            return metah + h

        hdr = prepHeader(contigs, args)

        with gzip.open(outFname, mode="wb") as ofn:
            ofn.write(str.encode("\n".join(hdr) + "\n"))
            for chr in data.keys():
                for i in range(len(data[chr]["xstart"])):
                    row = str.encode(
                        "%d %d\n" % (data[chr]["xstart"][i], data[chr]["xval"][i])
                    )
                    ofn.write(row)


class pafConvert:
    """
    updated to stream row at a time -
    slower but no room for the entire output if the input is 60GB+

    paf to xy and axis metadata
    Assumes pairs of points representing HiC contact pairs, or Mashmap sequence similarity hits. HiC data typically comes from pairs of haplotypes and is used to help assemble all
    the contigs into chromosomes

    These coordinate pairs are of 3 types - both ends on one or the other reference sequence, or one end on each for HiC pairs.
    Let's call these cis if the same reference and trans if between the two references. The ones in cis are likely to be consistent-ish between the two haplotypes, but the trans ones turn
    out to be very interesting...

    For HiC data, there should really be little difference - whether a haplotype contacts another haplotype or itself depends on the 3D folding and since the haplotypes are wound together into
    a helix, then the whole total contig length - about 2 meters for mammals are folded up into a 10μ^3 ball.

    python holoSeq_prepare_paf.py --inFile mUroPar1H1H2.paf --xclenfile mUroPar1H1suffix.len --yclenfile mUroPar1H2suffix.len --contig_sort VGPname --hap_indicator Suffix
    """

    def __init__(self, inFname, args, xcontigs, ycontigs, haps, xwidth, ywidth):
        """
        if rotating paf line at a time, better to pre-calculate the constants once rather than once for each point
        """
        self.rot = rotater(xwidth, ywidth)
        self.rotate = args.rotate
        self.rot.onepointRot = True
        self.rot.origin = (0, ywidth)
        self.outFprefix = inFname
        self.xcontigs = xcontigs
        self.ycontigs = ycontigs
        # have the axes set up so prepare the three plot x/y vectors
        # for a second pass to calculate all the coordinates.
        # adding tooltips just does not scale so abandoned - see the tooltip old version
        ncis1 = ncis2 = ntrans = 0
        hsId = holoSeqHeaders[1]
        self.inFname = inFname
        self.prepPafGZ(hsId, haps, xcontigs, ycontigs, args)
        if self.isGzip(inFname):
            with gzip.open(inFname, "rt") as f:
                self.readPAF(f)
        else:
            with open(inFname) as f:
                self.readPAF(f)
        self.cis1f.close()
        self.cis2f.close()
        self.transf.close()



    def readPAF(self, f):
        # might be a gz
        ncis1 = ncis2 = ntrans = 0
        for rowi, rows in enumerate(f):
            row = rows.strip().split()
            if len(row) > 7:
                c1 = row[0]
                c2 = row[5]
                n1 = int(row[2])
                n2 = int(row[7])
                H1 = getHap(c1)
                H2 = getHap(c2)
                if H1 != H2:  # trans
                    if H1 == haps[0]:  # x is h1 for trans - otherwise ignore
                        x = self.xcontigs[c1] + n1
                        y = self.ycontigs[c2] + n2
                        if self.rotate:
                            if y <= x: # lower triangle
                                x, y = self.rot.rotatecoords(x, y, True)
                            else:
                                x = None
                        if x is not None:
                            row = str.encode("%d %d\n" % (x, y))
                            self.transf.write(row)
                            ntrans += 1
                    else:
                        x = self.xcontigs[c2] + n2
                        y = self.ycontigs[c1] + n1
                        if self.rotate:
                            if y <= x: # lower triangle
                                x, y = self.rot.rotatecoords(xin=x, yin=y, adjust=True)
                            else:
                                x = None
                        if x is not None:
                            row = str.encode("%d %d\n" % (x, y))
                            self.transf.write(row)
                            ntrans += 1
                else:  # cis
                    if H1 == haps[0]:
                        x = self.xcontigs[c1] + n1
                        y = self.xcontigs[c2] + n2
                        if self.rotate:
                            if y <= x: # lower triangle
                                x, y = self.rot.rotatecoords(xin=x, yin=y, adjust=True)
                            else:
                                x = None
                        if x is not None:
                            row = str.encode("%d %d\n" % (x, y))
                            self.cis1f.write(row)
                            ncis1 += 1
                    else:
                        x = self.ycontigs[c1] + n1
                        y = self.ycontigs[c2] + n2
                        if self.rotate:
                            if y <= x: # lower triangle
                                x, y = self.rot.rotatecoords(xin=x, yin=y, adjust=True)
                            else:
                                x = None
                        if x is not None:
                            row = str.encode("%d %d\n" % (x, y))
                            self.cis2f.write(row)
                            ncis2 += 1
        log.debug("ncis1=%d, ncis2=%d, ntrans=%d" % (ncis1, ncis2, ntrans))

    def isGzip(self, inFname):
        with gzip.open(inFname, "r") as fh:
            try:
                fh.read(1)
                return True
            except gzip.BadGzipFile:
                log.info(
                    "inFname %s is not a gzip so will read as text" % inFname
                )
                return False

    def prepPafGZ(self, hsId, haps, xcontigs, ycontigs, args):
        """
        @v1HoloSeq2D for example
        """

        def prepHeader(
            haps,
            hsId,
            xcontigs,
            ycontigs,
            args,
            outf,
            subtitle,
            xclenfile,
            yclenfile,
            ax,
        ):
            """
            holoSeq output format - prepare gzip output channels
            """
            h = ["@%s %s %d" % (getHap(k), k, xcontigs[k]) for k in xcontigs.keys()]
            if len(haps) > 1:
                h += [
                    "@%s %s %d" % (getHap(k), k, ycontigs[k]) for k in ycontigs.keys()
                ]
            metah = [
                hsId,
                "@@heatmap",
                "@@title %s" % args.title + subtitle,
                "@@datasource %s" % "paf",
                "@@datafile %s" % self.inFname,
                "@@refURI %s" % args.refURI,
                "@@xclenfile %s" % xclenfile,
                "@@yclenfile %s" % yclenfile,
                "@@axes %s" % ax,
                "@@rotated %s" % self.rotate,
            ]

            outs = "\n".join(metah + h) + "\n"
            outf.write(str.encode(outs))
        
        fn1 = "%s_cis%s_hseq.gz" % (self.inFname, haps[0])
        if self.rotate:
            fn1 = "%s_rotated_cis%s_hseq.gz" % (self.inFname, haps[0])
        f1 = gzip.open(fn1, mode="wb")
        self.cis1f = io.BufferedWriter(f1, buffer_size=1024 * 1024)
        prepHeader(
            haps[0],
            hsId,
            xcontigs,
            xcontigs,
            args,
            self.cis1f,
            " Pairs on %s" % haps[0],
            xclenfile=args.xclenfile,
            yclenfile=args.xclenfile,
            ax=haps[0],
        )
        fn2 = "%s_cis%s_hseq.gz" % (self.inFname, haps[1])
        if self.rotate:
            fn2 = "%s_rotated_cis%s_hseq.gz" % (self.inFname, haps[1])
        f2 = gzip.open(fn2, mode="wb")
        self.cis2f = io.BufferedWriter(f2, buffer_size=1024 * 1024)
        prepHeader(
            haps[1],
            hsId,
            ycontigs,
            ycontigs,
            args,
            self.cis2f,
            " Pairs on %s" % haps[1],
            xclenfile=args.yclenfile,
            yclenfile=args.yclenfile,
            ax=haps[1],
        )
        fn3 = "%s_trans_hseq.gz" % (self.inFname)
        if self.rotate:
            fn3 = "%s_rotated_trans_hseq.gz" % (self.inFname)
        f3 = gzip.open(fn3, mode="wb")
        self.transf = io.BufferedWriter(f3, buffer_size=1024 * 1024)
        prepHeader(
            haps,
            hsId,
            xcontigs,
            ycontigs,
            args,
            self.transf,
            " Pairs on different haplotypes",
            xclenfile=args.xclenfile,
            yclenfile=args.yclenfile,
            ax="BOTH",
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="", epilog="")
    parser.add_argument(
        "--inFile",
        help="PAF with paired alignments, bigwig or gff3. Separate multiple with space. Must have matching .paf, .bw/.bigwig, or .gff/.gff3 extension to parse",
        nargs="+",
        default=[],
    )
    parser.add_argument(
        "--xclenfile",
        help="X axis contig names and lengths, whitespace delimited",
        required=True,
    )
    parser.add_argument(
        "--addH1",
        help="Bigwig and gff contigs can have H1 added if that matches the supplied xclenfile contig names. Not recommended - best to map against the right fasta",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--yclenfile",
        help="Optional Y axis contig names and lengths, whitespace delimited for different reference sequences",
        required=False,
    )
    
    parser.add_argument(
        "--title", help="Title for the plot", default="Plot title goes here"
    )
    parser.add_argument(
        "--contig_sort", help="VGPname, name, length, none", default="length"
    )
    parser.add_argument(
        "--refURI",
        help="URI for the genome reference sequence used for the coordinates for metadata",
        default="Unknown",
    )
    parser.add_argument(
        "--hap_indicator",
        help="None, Suffix (H[1,2]) Dashsuffix (_H...)",
        default="None",
    )  
    parser.add_argument(
        "--rotate",
        help="Rotate the 2D plot so the diagonal becomes the x axis",
        action="store_true",
        default=False,
    )
    parser.add_argument("--version", "-V", action="version", version="0.1")
    args = parser.parse_args()
    haps = []
    yhaps = []
    xcontigs, xhaps = getContigs(args.xclenfile)
    sxcontigs, xwidth = contsort(xcontigs, args)
    if args.yclenfile:
        ycontigs, yhaps = getContigs(args.yclenfile)
        sycontigs, ywidth = contsort(ycontigs, args)
    else:
        sycontigs = sxcontigs
        ywidth = xwidth
    for h in xhaps + yhaps:
        if h not in haps:
            haps.append(h)
    if len(haps) == 1:
        log.debug("extending haps %s" % haps)
        haps.append(haps[0])
    haps.sort()
    for f in args.inFile:
        ps = Path(f).suffix.lower()
        log.debug("inFile=%s, ftype = %s" % (f, ps))

        if ps in [".paf", ".paf.gz"]:
            p = pafConvert(f, args, sxcontigs, sycontigs, haps, xwidth, ywidth)
        elif ps in [".bw", ".bigwig"]:
            outf = "%s.hseq.gz" % f
            p = bwConvert(f, outf, args, sxcontigs)
        elif ps in [".gff3", ".gff"]:
            outf = "%s.hseq.gz" % f
            p = gffConvert(f, outf, sxcontigs, args)
        else:
            log.warn("%s unknown type - cannot process" % ps)
    logging.shutdown()
